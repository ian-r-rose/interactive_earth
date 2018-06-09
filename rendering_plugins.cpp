#include <cmath>
#include "rendering_plugins.h"
#include "color.h"

extern color (*colormap)(double);

void Core::setup()
{
  {
    const unsigned long n_triangles = ntheta;
    const unsigned long n_vertices = n_triangles+1;

    vertices = new GLfloat[ n_vertices * coordinates_per_vertex ];
    vertex_colors = new GLfloat[ n_vertices * colors_per_vertex ];
    triangle_vertex_indices = new GLuint[ n_triangles * vertices_per_triangle ];

    const float dtheta = 2.*M_PI/n_triangles;
    const float r = r_inner;
    
    color core_color = colormap(1.0);

    //one vertex at the origin
    vertices[0] = 0.0f;
    vertices[1] = 0.0f - (zoom_factor - 1.0f);
    vertex_colors[0] = core_color.R;
    vertex_colors[1] = core_color.G; 
    vertex_colors[2] = core_color.B; 

    unsigned long v = 2, c=3, i=0; //start at the next vertex index
    for (unsigned long n = 0; n < n_triangles; ++n)
    {
      vertex_colors[c + 0] = core_color.R;
      vertex_colors[c + 1] = core_color.G;
      vertex_colors[c + 2] = core_color.B;

      triangle_vertex_indices[i + 0] = 0;
      triangle_vertex_indices[i + 1] = n+1;
      triangle_vertex_indices[i + 2] = (n == n_triangles-1 ? 1 : n+2);

      v += coordinates_per_vertex;
      c += colors_per_vertex;
      i += vertices_per_triangle;
    }      

    glGenBuffers(1, &plugin_triangle_vertex_indices);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, plugin_triangle_vertex_indices);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint)*n_triangles*vertices_per_triangle, triangle_vertex_indices, GL_STATIC_DRAW);

    glGenBuffers(1, &plugin_vertices);

    glGenBuffers(1, &plugin_vertex_colors);
    glBindBuffer(GL_ARRAY_BUFFER, plugin_vertex_colors);
    glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*n_vertices*colors_per_vertex, vertex_colors, GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
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
    "attribute vec2 plugin_coord2d;"
    "attribute vec3 plugin_v_color;"
    "varying vec3 f_color;"
    "void main(void) {"
    "  f_color = plugin_v_color;"
    "  gl_Position = vec4(plugin_coord2d, 0.0, 1.0);"
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

  plugin_program = glCreateProgram();
  glAttachShader(plugin_program, vs);
  glAttachShader(plugin_program, fs);
  glLinkProgram(plugin_program);
  glGetProgramiv(plugin_program, GL_LINK_STATUS, &link_ok);
  if (!link_ok) {
    fprintf(stderr, "glLinkProgram:");
    return;
  }

  const char* attribute_name = "plugin_coord2d";
  plugin_attribute_coord2d = glGetAttribLocation(plugin_program, attribute_name);
  if (plugin_attribute_coord2d == -1) {
    fprintf(stderr, "Could not bind attribute %s\n", attribute_name);
    return;
  }
  attribute_name = "plugin_v_color";
  plugin_attribute_v_color = glGetAttribLocation(plugin_program, attribute_name);
  if (plugin_attribute_v_color == -1) {
    fprintf(stderr, "Could not bind attribute %s\n", attribute_name);
    return;
  }
  return;
}

void Core::draw()
{
  const unsigned long n_triangles = ntheta;
  const unsigned long n_vertices = n_triangles+1;
  const float dtheta = 2.*M_PI/n_triangles;

  unsigned long v = 2; //start at the next vertex index after the one in the center
  for (unsigned long n = 0; n < n_triangles; ++n)
  {
    const float theta = dtheta*n;
    //Vertices of the inner ellipse
    vertices[v + 0] = zoom_factor*r_inner*std::cos(theta);
    vertices[v + 1] = zoom_factor*r_inner*std::sin(theta) - (zoom_factor-1.0f);

    v += coordinates_per_vertex;
  }

  glUseProgram(plugin_program);

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, plugin_triangle_vertex_indices);

  glEnableVertexAttribArray(plugin_attribute_coord2d);

  glBindBuffer(GL_ARRAY_BUFFER, plugin_vertices);
  glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*n_vertices*coordinates_per_vertex, vertices, GL_DYNAMIC_DRAW);
  glVertexAttribPointer(
    plugin_attribute_coord2d, // attribute
    2,                 // number of elements per vertex (x,y)
    GL_FLOAT,          // the type of each element
    GL_FALSE,          // take our values as-is
    0,                 // no extra data between each position
    0                  // offset of first element
  );

  glEnableVertexAttribArray(plugin_attribute_v_color);
  glBindBuffer(GL_ARRAY_BUFFER, plugin_vertex_colors);
  glVertexAttribPointer(
    plugin_attribute_v_color, // attribute
    3,                 // number of elements per verte (r,g,b)
    GL_FLOAT,          // the type of each element
    GL_FALSE,          // take our values as-is
    0,                 // no extra data between each position
    0);

  glDrawElements(GL_TRIANGLES, n_triangles*vertices_per_triangle, GL_UNSIGNED_INT, 0);

  glDisableVertexAttribArray(plugin_attribute_coord2d);
  glDisableVertexAttribArray(plugin_attribute_v_color);

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
  glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void Core::cleanup()
{
  delete[] vertices;
  delete[] vertex_colors;
  delete[] triangle_vertex_indices;

  glDeleteProgram(plugin_program);
  glDeleteBuffers(1, &plugin_vertices);
  glDeleteBuffers(1, &plugin_vertex_colors);
  glDeleteBuffers(1, &plugin_triangle_vertex_indices);
}
