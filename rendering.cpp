#include "convection.h"
#include "color.h"

#ifndef LEGACY_OPENGL
GLuint program;
GLuint vbo_vertices;
GLuint vbo_colors;
GLuint ibo_triangle_vertex_indices;
GLint attribute_coord2d;
GLint attribute_v_color;
#endif //LEGACY_OPENGL

void ConvectionSimulator::setup_opengl()
{
#ifndef LEGACY_OPENGL
  //Setup the vertices, indices, and colors
  {
    GLfloat DX = 2.0/grid.nx;
    GLfloat DY = 2.0/grid.ny;
    const short triangles_per_quad = 2;
    const short vertices_per_triangle = 3;
    const short coordinates_per_vertex = 2;
    const short colors_per_vertex = 3;
    const unsigned long n_triangles = (grid.nx-1) * (grid.ny-1) * triangles_per_quad;
    const unsigned long n_vertices = grid.nx * grid.ny;
    

    vertices = new GLfloat[ n_vertices * coordinates_per_vertex ];
    vertex_colors = new GLfloat[ n_vertices * colors_per_vertex ];
    triangle_vertex_indices = new GLuint[ n_triangles * vertices_per_triangle ];
    
    unsigned long v=0, i=0;
    for( StaggeredGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
    {
      //Vertex for this cell
      vertices[v + 0] = cell->xindex()*DX-1.0f;
      vertices[v + 1] = cell->yindex()*DY-1.0f;

      if ( !cell->at_top_boundary() && !cell->at_right_boundary() )
      {
        //First triangle
        triangle_vertex_indices[i + 0] = cell->self();
        triangle_vertex_indices[i + 1] = cell->self()+grid.nx+1;
        triangle_vertex_indices[i + 2] = cell->self()+grid.nx;

        //Second triangle
        triangle_vertex_indices[i + 3] = cell->self();
        triangle_vertex_indices[i + 4] = cell->self() + 1;
        triangle_vertex_indices[i + 5] = cell->self() + grid.nx + 1;
   
        i += triangles_per_quad * vertices_per_triangle;
      }

      v += coordinates_per_vertex;
    }

    glGenBuffers(1, &ibo_triangle_vertex_indices);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo_triangle_vertex_indices);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint)*n_triangles*vertices_per_triangle, triangle_vertex_indices, GL_STATIC_DRAW);

    glGenBuffers(1, &vbo_vertices);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_vertices);
    glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*n_vertices*coordinates_per_vertex, vertices, GL_STATIC_DRAW);

    glGenBuffers(1, &vbo_colors);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_colors);
    glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*n_vertices*colors_per_vertex, vertex_colors, GL_DYNAMIC_DRAW);
  }

  GLint compile_ok = GL_FALSE, link_ok = GL_FALSE;

  GLuint vs = glCreateShader(GL_VERTEX_SHADER);
  const char *vs_source =
#ifdef GL_ES_VERSION_2_0
    "#version 100\n"  // OpenGL ES 2.0
    "precision mediump float;"
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
    return;
  }

  GLuint fs = glCreateShader(GL_FRAGMENT_SHADER);
  const char *fs_source =
#ifdef GL_ES_VERSION_2_0
    "#version 100\n"  // OpenGL ES 2.0
    "precision mediump float;"
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
    return;
  }

  program = glCreateProgram();
  glAttachShader(program, vs);
  glAttachShader(program, fs);
  glLinkProgram(program);
  glGetProgramiv(program, GL_LINK_STATUS, &link_ok);
  if (!link_ok) {
    fprintf(stderr, "glLinkProgram:");
    return;
  }

  const char* attribute_name = "coord2d";
  attribute_coord2d = glGetAttribLocation(program, attribute_name);
  if (attribute_coord2d == -1) {
    fprintf(stderr, "Could not bind attribute %s\n", attribute_name);
    return;
  }
  attribute_name = "v_color";
  attribute_v_color = glGetAttribLocation(program, attribute_name);
  if (attribute_v_color == -1) {
    fprintf(stderr, "Could not bind attribute %s\n", attribute_name);
    return;
  }

#endif //LEGACY_OPENGL

  return;
}

void ConvectionSimulator::cleanup_opengl()
{
#ifndef LEGACY_OPENGL
  delete[] vertices;
  delete[] vertex_colors;
  delete[] triangle_vertex_indices;

  glDeleteProgram(program);
  glDeleteBuffers(1, &vbo_vertices);
  glDeleteBuffers(1, &vbo_colors);
  glDeleteBuffers(1, &ibo_triangle_vertex_indices);
#endif //LEGACY_OPENGL
}

void ConvectionSimulator::draw( bool draw_composition)
{

#ifndef LEGACY_OPENGL
  const short triangles_per_quad = 2;
  const short vertices_per_triangle = 3;
  const short colors_per_vertex = 3;
  const unsigned long n_triangles = (grid.nx-1) * (grid.ny-1) * triangles_per_quad;
  const unsigned long n_vertices = grid.nx * grid.ny;

  unsigned long v=0;
  for( StaggeredGrid::iterator cell = grid.begin(); !cell->at_top_boundary(); ++cell)
  {
    color c;
    if(draw_composition)
      c = hot(C[cell->self()]);
    else
      c = hot(T[cell->self()]);

    vertex_colors[v + 0] = c.R;
    vertex_colors[v + 1] = c.G;
    vertex_colors[v + 2] = c.B;

    v += colors_per_vertex;
  }
    
  glClearColor(1.0, 1.0, 1.0, 1.0);
  glClear(GL_COLOR_BUFFER_BIT);

  glUseProgram(program);
  glEnableVertexAttribArray(attribute_coord2d);
  /* Describe our vertices array*/
  glBindBuffer(GL_ARRAY_BUFFER, vbo_vertices);
  glVertexAttribPointer(
    attribute_coord2d, // attribute
    2,                 // number of elements per vertex (x,y)
    GL_FLOAT,          // the type of each element
    GL_FALSE,          // take our values as-is
    0,                 // no extra data between each position
    0                  // offset of first element
  );

  glEnableVertexAttribArray(attribute_v_color);
  /* Describe our color array*/
  glBindBuffer(GL_ARRAY_BUFFER, vbo_colors);
  glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*n_vertices*colors_per_vertex, vertex_colors, GL_DYNAMIC_DRAW);
  glVertexAttribPointer(
    attribute_v_color, // attribute
    3,                 // number of elements per verte (r,g,b)
    GL_FLOAT,          // the type of each element
    GL_FALSE,          // take our values as-is
    0,                 // no extra data between each position
    0);

  /* Draw the triangles */
  glDrawElements(GL_TRIANGLES, n_triangles*vertices_per_triangle, GL_UNSIGNED_INT, 0);

  glDisableVertexAttribArray(attribute_coord2d);
  glDisableVertexAttribArray(attribute_v_color);


#else 

  double DX = 2.0/grid.nx;
  double DY = 2.0/grid.ny;

  glClearColor(0.0, 0.0, 0.0, 0.0);
  glClear(GL_COLOR_BUFFER_BIT);
  glOrtho(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0);
  glBegin(GL_TRIANGLE_STRIP);
  for( StaggeredGrid::iterator cell = grid.begin(); !cell->at_top_boundary(); ++cell)
  {
    if (cell->at_left_boundary())
      glBegin(GL_TRIANGLE_STRIP);

    if( !cell->at_right_boundary() )
    {
      color c_s, c_u;
      if (draw_composition)
      {
        c_s = hot(C[cell->self()]);
        c_u = hot(C[cell->up()]);
      }
      else 
      {
        c_s = hot(T[cell->self()]);
        c_u = hot(T[cell->up()]);
      }

      glColor3f(c_s.R, c_s.G, c_s.B);
      glVertex2f(cell->xindex()*DX-1.0, cell->yindex()*DY-1.0);
      glColor3f(c_u.R, c_u.G, c_u.B);
      glVertex2f((cell->xindex())*DX-1.0, (cell->yindex()+1)*DY-1.0);
    }
    else
      glEnd();
  }
  glFlush();

#endif //LEGACY_OPENGL

}
  
