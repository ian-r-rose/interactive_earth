#include "stokes.h"
#include "color.h"
#include <iostream>

GLuint program;
GLuint vbo_triangle;
GLuint vbo_triangle_colors;
GLint attribute_coord2d;
GLint attribute_v_color;

void print_log(GLuint object)
{
  GLint log_length = 0;
  if (glIsShader(object))
    glGetShaderiv(object, GL_INFO_LOG_LENGTH, &log_length);
  else if (glIsProgram(object))
    glGetProgramiv(object, GL_INFO_LOG_LENGTH, &log_length);
  else {
    fprintf(stderr, "printlog: Not a shader or a program\n");
    return;
  }

  char* log = (char*)malloc(log_length);

  if (glIsShader(object))
    glGetShaderInfoLog(object, log_length, NULL, log);
  else if (glIsProgram(object))
    glGetProgramInfoLog(object, log_length, NULL, log);

  fprintf(stderr, "%s", log);
  free(log);
}

int StokesSolver::setup_opengl()
{
  //Setup the vertices, indices, and colors
  {
    GLfloat DX = 2.0/grid.nx;
    GLfloat DY = 2.0/grid.ny;
    const short triangles_per_quad = 2;
    const short vertices_per_triangle = 3;
    const short coordinates_per_vertex = 2;
    const short colors_per_vertex = 3;
    const unsigned long n_triangles = grid.nx * (grid.ny-1) * triangles_per_quad;
    

    triangle_vertices = new GLfloat[ n_triangles * vertices_per_triangle * coordinates_per_vertex ];
    vertex_colors = new GLfloat[ n_triangles * vertices_per_triangle * colors_per_vertex ];
//    vertex_indices = new GLint[ n_triangles * vertices_per_triangle ];
    
    unsigned long i=0;
    for( StaggeredGrid::iterator cell = grid.begin(); !cell->at_top_boundary(); ++cell)
    {
      //Triangle 1, vertex 1
      triangle_vertices[i + 0] = cell->xindex()*DX-1.0f;
      triangle_vertices[i + 1] = cell->yindex()*DY-1.0f;
      //Triangle 1, vertex 2
      triangle_vertices[i + 2] = cell->xindex()*DX-1.0f;
      triangle_vertices[i + 3] = (cell->yindex()+1)*DY-1.0f;
      //Triangle 1, vertex 3
      triangle_vertices[i + 4] = (cell->xindex()+1)*DX-1.0f;
      triangle_vertices[i + 5] = cell->yindex()*DY-1.0f;

      //Triangle 2, vertex 1
      triangle_vertices[i + 6] = (cell->xindex()+1)*DX-1.0f;
      triangle_vertices[i + 7] = cell->yindex()*DY-1.0f;
      //Triangle 2, vertex 2
      triangle_vertices[i + 8] = cell->xindex()*DX-1.0f;
      triangle_vertices[i + 9] = (cell->yindex()+1)*DY-1.0f;
      //Triangle 2, vertex 3
      triangle_vertices[i + 10] = (cell->xindex()+1)*DX-1.0f;
      triangle_vertices[i + 11] = (cell->yindex()+1)*DY-1.0f;
   
      i+=triangles_per_quad * vertices_per_triangle * coordinates_per_vertex;
    }

    glGenBuffers(1, &vbo_triangle);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_triangle);
    glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*n_triangles*vertices_per_triangle*coordinates_per_vertex, triangle_vertices, GL_STATIC_DRAW);

    glGenBuffers(1, &vbo_triangle_colors);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_triangle_colors);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertex_colors), vertex_colors, GL_STATIC_DRAW);
  }

  GLint compile_ok = GL_FALSE, link_ok = GL_FALSE;

  GLuint vs = glCreateShader(GL_VERTEX_SHADER);
  const char *vs_source =
#ifdef GL_ES_VERSION_2_0
    "#version 100\n"  // OpenGL ES 2.0
#else
    "#version 120\n"  // OpenGL 2.1
#endif
    "attribute vec2 coord2d;"
    "attribute vec3 v_color;"
    "varying vec3 f_color;"
    "void main(void) {"
    "  f_color = v_color;"
    "  gl_Position = vec4(coord2d, 0.0, 1.0);"
    "}";

  glShaderSource(vs, 1, &vs_source, NULL);
  glCompileShader(vs);
  glGetShaderiv(vs, GL_COMPILE_STATUS, &compile_ok);
  if (!compile_ok) {
    fprintf(stderr, "Error in vertex shader\n");
    return 0;
  }

  GLuint fs = glCreateShader(GL_FRAGMENT_SHADER);
  const char *fs_source =
#ifdef GL_ES_VERSION_2_0
    "#version 100\n"  // OpenGL ES 2.0
#else
    "#version 120\n"  // OpenGL 2.1
#endif
    "varying vec3 f_color;"
    "void main(void) {"
    "  gl_FragColor = vec4(f_color, 1.0);"
    "}";
  glShaderSource(fs, 1, &fs_source, NULL);
  glCompileShader(fs);
  glGetShaderiv(fs, GL_COMPILE_STATUS, &compile_ok);
  if (!compile_ok) {
    fprintf(stderr, "Error in fragment shader\n");
    return 0;
  }

  program = glCreateProgram();
  glAttachShader(program, vs);
  glAttachShader(program, fs);
  glLinkProgram(program);
  glGetProgramiv(program, GL_LINK_STATUS, &link_ok);
  if (!link_ok) {
    fprintf(stderr, "glLinkProgram:");
    print_log(program);
    return 0;
  }

  const char* attribute_name = "coord2d";
  attribute_coord2d = glGetAttribLocation(program, attribute_name);
  if (attribute_coord2d == -1) {
    fprintf(stderr, "Could not bind attribute %s\n", attribute_name);
    return 0;
  }
  attribute_name = "v_color";
  attribute_v_color = glGetAttribLocation(program, attribute_name);
  if (attribute_v_color == -1) {
    fprintf(stderr, "Could not bind attribute %s\n", attribute_name);
    return 0;
  }

  return 1;
}

void StokesSolver::draw()
{
  const short triangles_per_quad = 2;
  const short vertices_per_triangle = 3;
  const short coordinates_per_vertex = 2;
  const short colors_per_vertex = 3;
  const unsigned long n_triangles = grid.nx * (grid.ny-1) * triangles_per_quad;

  unsigned long i=0;
  for( StaggeredGrid::iterator cell = grid.begin(); !cell->at_top_boundary(); ++cell)
  {
    color c_s = hot(T[cell->self()]);
    color c_u = hot(T[cell->up()]);
    color c_r = hot(T[cell->right()]);
    color c_ur = hot(T[cell->upright()]);

    //Triangle 1, vertex 1
    vertex_colors[i + 0] = c_s.R;
    vertex_colors[i + 1] = c_s.G;
    vertex_colors[i + 2] = c_s.B;
    //Triangle 1, vertex 2
    vertex_colors[i + 3] = c_u.R;
    vertex_colors[i + 4] = c_u.G;
    vertex_colors[i + 5] = c_u.B;
    //Triangle 1, vertex 3
    vertex_colors[i + 6] = c_r.R;
    vertex_colors[i + 7] = c_r.G;
    vertex_colors[i + 8] = c_r.B;
    //Triangle 2, vertex 1
    vertex_colors[i + 9] = c_r.R;
    vertex_colors[i + 10] = c_r.G;
    vertex_colors[i + 11] = c_r.B;
    //Triangle 2, vertex 2
    vertex_colors[i + 12] = c_u.R;
    vertex_colors[i + 13] = c_u.G;
    vertex_colors[i + 14] = c_u.B;
    //Triangle 2, vertex 3
    vertex_colors[i + 15] = c_ur.R;
    vertex_colors[i + 16] = c_ur.G;
    vertex_colors[i + 17] = c_ur.B;

    i += triangles_per_quad*vertices_per_triangle*colors_per_vertex;
  }
    
  glClearColor(1.0, 1.0, 1.0, 1.0);
  glClear(GL_COLOR_BUFFER_BIT);

  glUseProgram(program);
  glEnableVertexAttribArray(attribute_coord2d);
  /* Describe our vertices array to OpenGL (it can't guess its format automatically) */
  glBindBuffer(GL_ARRAY_BUFFER, vbo_triangle);
  glVertexAttribPointer(
    attribute_coord2d, // attribute
    2,                 // number of elements per vertex, here (x,y)
    GL_FLOAT,          // the type of each element
    GL_FALSE,          // take our values as-is
    0,                 // no extra data between each position
    0                  // offset of first element
  );

  glEnableVertexAttribArray(attribute_v_color);
  /* Describe our vertices array to OpenGL (it can't guess its format automatically) */
  glBindBuffer(GL_ARRAY_BUFFER, vbo_triangle_colors);
  glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*n_triangles*vertices_per_triangle*colors_per_vertex, vertex_colors, GL_STATIC_DRAW);
  glVertexAttribPointer(
    attribute_v_color, // attribute
    3,                 // number of elements per vertex, here (x,y)
    GL_FLOAT,          // the type of each element
    GL_FALSE,          // take our values as-is
    0,                 // no extra data between each position
    0);

  /* Push each element in buffer_vertices to the vertex shader */
  glDrawArrays(GL_TRIANGLES, 0, n_triangles*vertices_per_triangle);

  glDisableVertexAttribArray(attribute_coord2d);
  glDisableVertexAttribArray(attribute_v_color);

}
  
