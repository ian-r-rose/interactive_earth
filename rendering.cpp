#include <cmath>
#include "convection.h"
#include "color.h"

#ifndef LEGACY_OPENGL
static GLuint program;
static GLuint vbo_vertices;
static GLuint vbo_colors;
static GLuint ibo_triangle_vertex_indices;
static GLint attribute_coord2d;
static GLint attribute_v_color;
#endif //LEGACY_OPENGL


extern const double flattening;

void ConvectionSimulator::setup_opengl()
{
#ifndef LEGACY_OPENGL
  //Setup the vertices, indices, and colors
  {
    GLfloat DX = 2.0*M_PI/grid.ntheta;
    GLfloat DY = 1.0/(grid.nr-1);
    const short triangles_per_quad = 2;
    const short vertices_per_triangle = 3;
    const short coordinates_per_vertex = 2;
    const short colors_per_vertex = 3;
    const unsigned long n_triangles = grid.ntheta * (grid.nr-1) * triangles_per_quad;
    const unsigned long n_vertices = grid.ntheta * grid.nr;


    vertices = new GLfloat[ n_vertices * coordinates_per_vertex ];
    vertex_colors = new GLfloat[ n_vertices * colors_per_vertex ];
    triangle_vertex_indices = new GLuint[ n_triangles * vertices_per_triangle ];

    unsigned long i=0;
    for( RegularGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
    {
      if ( !cell->at_top_boundary() )
      {
        //First triangle
        triangle_vertex_indices[i + 0] = cell->self();
        triangle_vertex_indices[i + 1] = cell->upright();
        triangle_vertex_indices[i + 2] = cell->up();

        //Second triangle
        triangle_vertex_indices[i + 3] = cell->self();
        triangle_vertex_indices[i + 4] = cell->right();
        triangle_vertex_indices[i + 5] = cell->upright();

        i += triangles_per_quad * vertices_per_triangle;
      }
    }

    glGenBuffers(1, &ibo_triangle_vertex_indices);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo_triangle_vertex_indices);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint)*n_triangles*vertices_per_triangle, triangle_vertex_indices, GL_STATIC_DRAW);

    glGenBuffers(1, &vbo_vertices);
    glGenBuffers(1, &vbo_colors);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
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

void ConvectionSimulator::draw( bool draw_vorticity )
{
  double displacement_factor = 1.0;
#ifndef LEGACY_OPENGL
  const short triangles_per_quad = 2;
  const short vertices_per_triangle = 3;
  const short colors_per_vertex = 3;
  const short coordinates_per_vertex = 2;
  const unsigned long n_triangles = grid.ntheta * (grid.nr-1) * triangles_per_quad;
  const unsigned long n_vertices = grid.ntheta * grid.nr;

  double max_V = 0.;
  double max_S = 0.;
  for( RegularGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
  {
    if(std::abs(V[cell->self()]) > max_V)
      max_V = std::abs(V[cell->self()]);
    if(std::abs(stream[cell->self()]) > max_S)
      max_S = std::abs(stream[cell->self()]);
  }

  unsigned long v=0, i=0;
  for( RegularGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
  {
    //Coordinates for this cell
    const float r = cell->radius();
    const float theta = cell->location().theta;

    //This is kind of cool: add ellipticity to the domain rendering
    //using an epicycle on the normal rendering.
    //TODO: offsetting the epicycle seems to work with 2*angle.
    //Why not 1*angle? Figure this out.
    vertices[v + 0] = r * std::cos(theta);
    vertices[v + 1] = r * std::sin(theta);

    color c;
    if (draw_vorticity)
      c = hot(V[cell->self()]);
    else
      c = hot(T[cell->self()]);

    vertex_colors[i + 0] = c.R;
    vertex_colors[i + 1] = c.G;
    vertex_colors[i + 2] = c.B;

    i += colors_per_vertex;
    v += coordinates_per_vertex;
  }

  glClearColor(0.0, 0.0, 0.0, 1.0);
  glClear(GL_COLOR_BUFFER_BIT);

  glUseProgram(program);

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo_triangle_vertex_indices);

  glEnableVertexAttribArray(attribute_coord2d);
  /* Describe our vertices array*/
  glBindBuffer(GL_ARRAY_BUFFER, vbo_vertices);
  glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*n_vertices*coordinates_per_vertex, vertices, GL_DYNAMIC_DRAW);
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

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
  glBindBuffer(GL_ARRAY_BUFFER, 0);

#else

  GLfloat DX = 2.0*M_PI/(grid.ntheta);
  GLfloat DY = 1.0/(grid.nr-1);

  glClearColor(0.0, 0.0, 0.0, 0.0);
  glClear(GL_COLOR_BUFFER_BIT);
  glOrtho(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0);
  glBegin(GL_TRIANGLE_STRIP);
  for( RegularGrid::iterator cell = grid.begin(); !cell->at_top_boundary(); ++cell)
  {
    if (cell->at_left_boundary() && cell->self() != 0.)
      glBegin(GL_TRIANGLE_STRIP);

    color c_s = hot(T[cell->self()]);
    color c_u = hot(T[cell->up()]);

    const float r = grid.r_inner + (cell->yindex()*DY * (1.0f - grid.r_inner) );
    const float r_up = grid.r_inner + ((cell->yindex()+1)*DY * (1.0f - grid.r_inner) );
    const float theta = cell->xindex()*DX;

    glColor3f(c_s.R, c_s.G, c_s.B);
    glVertex2f(r * std::cos(theta), r * std::sin(theta));
    glColor3f(c_u.R, c_u.G, c_u.B);
    glVertex2f(r_up * std::cos(theta), r_up * std::sin(theta));

    //Add last strip for periodicity
    if (cell->at_right_boundary() )
    {
      color c_s, c_u;
      if (draw_vorticity)
      {
        c_s = hot(V[cell->self()]);
        c_u = hot(V[cell->up()]);
      }
      else
      {
        c_s = hot(T[cell->self()]);
        c_u = hot(T[cell->up()]);
      }

      const float r = grid.r_inner + (cell->yindex()*DY * (1.0f - grid.r_inner) );
      const float r_up = grid.r_inner + ((cell->yindex()+1)*DY * (1.0f - grid.r_inner) );
      const float theta = 0.0;

      glColor3f(c_s.R, c_s.G, c_s.B);
      glVertex2f(r * std::cos(theta), r * std::sin(theta));
      glColor3f(c_u.R, c_u.G, c_u.B);
      glVertex2f(r_up * std::cos(theta), r_up * std::sin(theta));
      glEnd();
    }
  }
  glFlush();

#endif //LEGACY_OPENGL

}

