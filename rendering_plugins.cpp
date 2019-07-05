#include <cmath>
#include <iostream>
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

    color core_color = colormap(1.0);

    //one vertex at the origin
    vertices[0] = 0.0f;
    vertices[1] = 0.0f;
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

  const float angle = sim.spin_angle();
  const float a = 1.0f;
  const float b = 1.0f-flattening;

  unsigned long v = 2; //start at the next vertex index after the one in the center
  for (unsigned long n = 0; n < n_triangles; ++n)
  {
    const float theta = dtheta*n;
    //Vertices of the inner ellipse
    vertices[v + 0] = r_inner*(a+b)/2. * std::cos(theta) +
                      r_inner*(a-b)/2. * std::cos(-theta+ 2.*(angle+M_PI/2.));
    vertices[v + 1] = r_inner*(a+b)/2. * std::sin(theta) +
                      r_inner*(a-b)/2. * std::sin(-theta+ 2.*(angle+M_PI/2.));

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
    3,                 // number of elements per vertex (r,g,b)
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

void Axis::setup()
{
  {
    const unsigned long n_triangles = 2 /*axis*/ + 2 /*arrowheads*/;
    const unsigned long n_vertices = 4 /*axis*/ + 6 /*arrowheads*/;

    vertices = new GLfloat[ n_vertices * coordinates_per_vertex ];
    vertex_colors = new GLfloat[ n_vertices * colors_per_vertex ];
    triangle_vertex_indices = new GLuint[ n_triangles * vertices_per_triangle ];

    color arrow_color;
    arrow_color.R = 0.15; arrow_color.G=0.15; arrow_color.B=0.15;
    unsigned short c = 0;
    for (unsigned short i=0; i<n_vertices; ++i)
    {
      vertex_colors[c + 0] = arrow_color.R;
      vertex_colors[c + 1] = arrow_color.G;
      vertex_colors[c + 2] = arrow_color.B;
      c += colors_per_vertex;
    }

    //Axis
    triangle_vertex_indices[0] = 0;
    triangle_vertex_indices[1] = 1;
    triangle_vertex_indices[2] = 2;
    triangle_vertex_indices[3] = 0;
    triangle_vertex_indices[4] = 2;
    triangle_vertex_indices[5] = 3;
    //Arrowheads
    triangle_vertex_indices[6] = 4;
    triangle_vertex_indices[7] = 5;
    triangle_vertex_indices[8] = 6;
    triangle_vertex_indices[9] = 7;
    triangle_vertex_indices[10] = 8;
    triangle_vertex_indices[11] = 9;

    glGenBuffers(1, &plugin_triangle_vertex_indices);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, plugin_triangle_vertex_indices);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint)*n_triangles*vertices_per_triangle, triangle_vertex_indices, GL_STATIC_DRAW);

    glGenBuffers(1, &plugin_vertices);
    glBindBuffer(GL_ARRAY_BUFFER, plugin_vertices);
    glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*n_vertices*coordinates_per_vertex, vertices, GL_DYNAMIC_DRAW);

    glGenBuffers(1, &plugin_vertex_colors);
    glBindBuffer(GL_ARRAY_BUFFER, plugin_vertex_colors);
    glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*n_vertices*colors_per_vertex, vertex_colors, GL_DYNAMIC_DRAW);

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

void Axis::draw()
{
  const unsigned long n_triangles = 2 /*axis*/ + 2 /*arrowheads*/;
  const unsigned long n_vertices = 4 /*axis*/ + 6 /*arrowheads*/;

  const float d2r = M_PI/180.0;
  const float theta = float( sim.spin_angle() );
  const float r = 0.8 * r_inner * (1.0f-flattening);
  const float arrowhead = 0.1*r_inner;
  const float dtheta = 0.5/r_inner * d2r;

  //Draw the axis
  vertices[0] = r * std::cos(theta - dtheta);
  vertices[1] = r * std::sin(theta - dtheta);
  vertices[2] = -r * std::cos(theta + dtheta);
  vertices[3] = -r * std::sin(theta + dtheta);
  vertices[4] = -r * std::cos(theta - dtheta);
  vertices[5] = -r * std::sin(theta - dtheta);
  vertices[6] = r * std::cos(theta + dtheta);
  vertices[7] = r * std::sin(theta + dtheta);
  //Draw the arrowheads

  //first arrowhead
  vertices[8] = (r + arrowhead) * std::cos(theta);
  vertices[9] = (r + arrowhead) * std::sin(theta);
  vertices[10] = (r - arrowhead) * std::cos(theta + 4.*dtheta);
  vertices[11] = (r - arrowhead) * std::sin(theta + 4.*dtheta);
  vertices[12] = (r - arrowhead) * std::cos(theta - 4.*dtheta);
  vertices[13] = (r - arrowhead) * std::sin(theta - 4.*dtheta);
  //second arrowhead
  vertices[14] = -(r + arrowhead) * std::cos(theta);
  vertices[15] = -(r + arrowhead) * std::sin(theta);
  vertices[16] = -(r - arrowhead) * std::cos(theta + 4.*dtheta);
  vertices[17] = -(r - arrowhead) * std::sin(theta + 4.*dtheta);
  vertices[18] = -(r - arrowhead) * std::cos(theta - 4.*dtheta);
  vertices[19] = -(r - arrowhead) * std::sin(theta - 4.*dtheta);

  glEnable(GL_BLEND);
#ifndef __EMSCRIPTEN__
  glEnable(GL_POLYGON_SMOOTH);
#endif
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

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
    3,                 // number of elements per vertex (r,g,b)
    GL_FLOAT,          // the type of each element
    GL_FALSE,          // take our values as-is
    0,                 // no extra data between each position
    0);

  glDrawElements(GL_TRIANGLES, n_triangles*vertices_per_triangle, GL_UNSIGNED_INT, 0);

  glDisableVertexAttribArray(plugin_attribute_coord2d);
  glDisableVertexAttribArray(plugin_attribute_v_color);

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
  glBindBuffer(GL_ARRAY_BUFFER, 0);

  glDisable(GL_BLEND);
#ifndef __EMSCRIPTEN__
  glDisable(GL_POLYGON_SMOOTH);
#endif
}

void Axis::cleanup()
{
  delete[] vertices;
  delete[] vertex_colors;
  delete[] triangle_vertex_indices;

  glDeleteProgram(plugin_program);
  glDeleteBuffers(1, &plugin_vertices);
  glDeleteBuffers(1, &plugin_vertex_colors);
  glDeleteBuffers(1, &plugin_triangle_vertex_indices);
}

void Seismograph::clear_record()
{
  unsigned int v = 0;
  for (unsigned int i=0; i<n_lines+1; ++i)
  {
    vertices[v + 1] = 0.f;
    v += coordinates_per_vertex;
  }
}

void Seismograph::setup()
{
  {
    vertices = new GLfloat[ n_vertices * coordinates_per_vertex ];
    vertex_colors = new GLfloat[ n_vertices * (colors_per_vertex+1) ];
    line_vertex_indices = new GLuint[ n_lines * vertices_per_line + 1 * vertices_per_triangle /*seismometer*/];

    //Begin with vertices all recording zero
    unsigned int v = 0;
    for (unsigned int i=0; i<n_lines+1; ++i)
    {
      vertices[v + 0] = -0.9f*r_inner*(1.0f-flattening) + float(i)/float(n_vertices) * 1.8f*r_inner*(1.0f-flattening);
      vertices[v + 1] = 0.f;
      v += coordinates_per_vertex;
    }

    color line_color;
    line_color.R = 0.0; line_color.G=0.0; line_color.B=0.0;
    unsigned int c = 0;
    for (unsigned int i=0; i<n_lines+1; ++i)
    {
      vertex_colors[c + 0] = line_color.R;
      vertex_colors[c + 1] = line_color.G;
      vertex_colors[c + 2] = line_color.B;

      //Add transparency to the edges
      if ( i < (n_lines+1)/6)
        vertex_colors[c + 3] = 6*float(i)/float(n_vertices);
      else if ( i > 4 * (n_lines+1)/6 )
        vertex_colors[c + 3] = 6*float(n_lines+1 - i)/float(n_lines+1);
      else
        vertex_colors[c + 3] = 1.0;
      c += colors_per_vertex+1;
    }

    for (unsigned int i =0; i<n_lines; ++i)
    {
      line_vertex_indices[vertices_per_line*i + 0] = i; 
      line_vertex_indices[vertices_per_line*i + 1] = i+1; 
    }

    //Fill the vertex information about the seismometer
    line_vertex_indices[ (n_lines*vertices_per_line) + 0 ] = n_vertices - 3;
    line_vertex_indices[ (n_lines*vertices_per_line) + 1 ] = n_vertices - 2;
    line_vertex_indices[ (n_lines*vertices_per_line) + 2 ] = n_vertices - 1;

    for (int i=vertices_per_triangle; i>0; --i)
    {
      vertex_colors[(n_vertices-i)*(colors_per_vertex+1) + 0] = 0.0;//line_color.R;
      vertex_colors[(n_vertices-i)*(colors_per_vertex+1)+ 1] = 0.3;//line_color.G;
      vertex_colors[(n_vertices-i)*(colors_per_vertex+1) + 2] = 1.0;//line_color.B;
      vertex_colors[(n_vertices-i)*(colors_per_vertex+1) + 3] = 0.8;
    }


    glGenBuffers(1, &plugin_line_vertex_indices);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, plugin_line_vertex_indices);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint)*(n_lines*vertices_per_line + vertices_per_triangle), line_vertex_indices, GL_STATIC_DRAW);

    glGenBuffers(1, &plugin_vertices);
    glBindBuffer(GL_ARRAY_BUFFER, plugin_vertices);
    glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*n_vertices*coordinates_per_vertex, vertices, GL_DYNAMIC_DRAW);

    glGenBuffers(1, &plugin_vertex_colors);
    glBindBuffer(GL_ARRAY_BUFFER, plugin_vertex_colors);
    glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*n_vertices*(colors_per_vertex+1), vertex_colors, GL_STATIC_DRAW);

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
    "attribute vec4 plugin_v_color;"
    "varying vec4 f_color;"
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
    "varying vec4 f_color;"
    "void main(void) {"
    "  gl_FragColor = f_color;"
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

void Seismograph::draw()
{
  //Begin with vertices all recording zero
  unsigned int v = 0;
  for (unsigned int i=0; i<n_lines+1; ++i)
  {
    vertices[v + 1] = vertices[v+3];
    v += coordinates_per_vertex;
  }
  vertices[(n_lines+1)*coordinates_per_vertex - 1] = sim.seismometer_reading() * 0.7*r_inner;

  //Fill the vertex information about the seismometer
  double theta, r;
  sim.seismometer_position(theta, r);
  r += r_inner;

  //Information about ellipticity
  const float angle = sim.spin_angle();
  const float a = 1.0f;
  const float b = 1.0f-flattening;

  //Position of the seismometer in (possibly) elliptical space 
  const float seis_x = r*(a+b)/2. * std::cos(theta) +
                       r*(a-b)/2. * std::cos(-theta+ 2.*(angle+M_PI/2.));
  const float seis_y = r*(a+b)/2. * std::sin(theta) +
                       r*(a-b)/2. * std::sin(-theta+ 2.*(angle+M_PI/2.));

  //Seismometer vertices
  vertices[ (n_lines+1)*coordinates_per_vertex + 0] = seis_x-0.04;
  vertices[ (n_lines+1)*coordinates_per_vertex + 1] = seis_y+0.04;
  vertices[ (n_lines+2)*coordinates_per_vertex + 0] = seis_x+0.04;
  vertices[ (n_lines+2)*coordinates_per_vertex + 1] = seis_y+0.04;
  vertices[ (n_lines+3)*coordinates_per_vertex + 0] = seis_x;
  vertices[ (n_lines+3)*coordinates_per_vertex + 1] = seis_y-0.04;

  glEnable(GL_BLEND);
#ifndef __EMSCRIPTEN__
  glEnable(GL_LINE_SMOOTH);
#endif
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  glLineWidth(2.0);
  glUseProgram(plugin_program);

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, plugin_line_vertex_indices);

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
    4,                 // number of elements per vertex (r,g,b,a)
    GL_FLOAT,          // the type of each element
    GL_FALSE,          // take our values as-is
    0,                 // no extra data between each position
    0);

  glDrawElements(GL_LINES, n_lines*vertices_per_line, GL_UNSIGNED_INT, 0);

  //Draw seismometer
  glDrawElements(GL_TRIANGLES, vertices_per_triangle, GL_UNSIGNED_INT, (void*)(((n_lines)*vertices_per_line)*sizeof(GLuint)));

  glDisableVertexAttribArray(plugin_attribute_coord2d);
  glDisableVertexAttribArray(plugin_attribute_v_color);

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
  glBindBuffer(GL_ARRAY_BUFFER, 0);

  glLineWidth(1.0);
  glDisable(GL_BLEND);
#ifndef __EMSCRIPTEN__
  glDisable(GL_LINE_SMOOTH);
#endif
}

void Seismograph::cleanup()
{
  delete[] vertices;
  delete[] vertex_colors;
  delete[] line_vertex_indices;

  glDeleteProgram(plugin_program);
  glDeleteBuffers(1, &plugin_vertices);
  glDeleteBuffers(1, &plugin_vertex_colors);
  glDeleteBuffers(1, &plugin_line_vertex_indices);
}

void ModeButton::setup()
{
  {
    const unsigned long n_triangles = 100;
    const unsigned long n_lines = 100;
    const unsigned long n_vertices = n_triangles + 1 + n_lines + 1;

    vertices = new GLfloat[ n_vertices * coordinates_per_vertex ];
    vertex_colors = new GLfloat[ n_vertices * colors_per_vertex ];
    vertex_indices = new GLuint[ n_triangles * vertices_per_triangle + n_lines * vertices_per_line ];

    //one vertex at the origin
    vertices[0] = mode_button_x;
    vertices[1] = mode_button_y;

    unsigned long v = 2, c=3, i=0; //start at the next vertex index
    for (unsigned long n = 0; n < n_triangles; ++n)
    {
      vertices[v + 0] = mode_button_x + button_radius * std::cos(float(n)/n_triangles * 2.0f*M_PI);
      vertices[v + 1] = mode_button_y + button_radius * std::sin(float(n)/n_triangles * 2.0f*M_PI);

      vertex_indices[i + 0] = 0;
      vertex_indices[i + 1] = n+1;
      vertex_indices[i + 2] = (n == n_triangles-1 ? 1 : n+2);

      v += coordinates_per_vertex;
      c += colors_per_vertex;
      i += vertices_per_triangle;
    }

    // Make the button icon vertices
    const float padding = 0.01f;
    for (unsigned long n = 0; n <= n_lines; n++)
    {
      const float x = 2.*M_PI*(float(n)/n_lines - 0.5);
      vertices[v + 0] = mode_button_x - button_radius + padding + n * 2.*(button_radius - padding)/n_lines;
      vertices[v + 1] =
        mode_button_y +
        (button_radius - padding) * std::cos(2.0*x/2.0*M_PI) * std::exp(-x*x/2.);

      v += coordinates_per_vertex;
    }

    unsigned long idx = n_triangles * vertices_per_triangle;
    for(unsigned long n = 0; n < n_lines; ++n)
    {
      vertex_indices[idx + 0] = n + n_triangles + 1;
      vertex_indices[idx + 1] = n + n_triangles + 2;

      idx += vertices_per_line;
    }

    for (unsigned long n = 0; n <= n_lines; ++n)
    {
      vertex_colors[c + 0 ] = 0.0;
      vertex_colors[c + 1 ] = 0.0;
      vertex_colors[c + 2 ] = 0.0;

      c += colors_per_vertex;
    }

    glGenBuffers(1, &plugin_vertex_indices);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, plugin_vertex_indices);
    glBufferData(
        GL_ELEMENT_ARRAY_BUFFER,
        sizeof(GLuint)*(n_triangles*vertices_per_triangle + n_lines*vertices_per_line),
        vertex_indices,
        GL_STATIC_DRAW
    );

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

void ModeButton::draw()
{
  const unsigned long n_triangles = 100;
  const unsigned long n_lines = 100;
  const unsigned long n_vertices = n_triangles + 1 + n_lines + 1;

  // Update the color if we are in seismic mode to indicate it is enabled.
  color *button_color;
  color enabled_color = { 1.0, 1.0, 1.0 };
  color disabled_color = { 0.5, 0.5, 0.5 };
  button_color = seismic_mode ? &enabled_color : &disabled_color;

  //one vertex at the origin
  vertex_colors[0] = button_color->R;
  vertex_colors[1] = button_color->G; 
  vertex_colors[2] = button_color->B; 

  unsigned long c=3; //start at the next vertex index
  for (unsigned long n = 0; n < n_triangles; ++n)
  {
    vertex_colors[c + 0] = button_color->R;
    vertex_colors[c + 1] = button_color->G;
    vertex_colors[c + 2] = button_color->B;

    c += colors_per_vertex;
  }

  glEnable(GL_BLEND);
#ifndef __EMSCRIPTEN__
  glEnable(GL_LINE_SMOOTH);
  glEnable(GL_POLYGON_SMOOTH);
#endif
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  glLineWidth(2.0);

  glUseProgram(plugin_program);

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, plugin_vertex_indices);

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
  glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*n_vertices*colors_per_vertex, vertex_colors, GL_STATIC_DRAW);
  glVertexAttribPointer(
    plugin_attribute_v_color, // attribute
    3,                 // number of elements per vertex (r,g,b)
    GL_FLOAT,          // the type of each element
    GL_FALSE,          // take our values as-is
    0,                 // no extra data between each position
    0);

  glDrawElements(GL_TRIANGLES, n_triangles*vertices_per_triangle, GL_UNSIGNED_INT, 0);
  glDrawElements(
      GL_LINES,
      n_lines*vertices_per_line,
      GL_UNSIGNED_INT,
      (void*)(((n_triangles)*vertices_per_triangle)*sizeof(GLuint))
  );

  glDisableVertexAttribArray(plugin_attribute_coord2d);
  glDisableVertexAttribArray(plugin_attribute_v_color);

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
  glBindBuffer(GL_ARRAY_BUFFER, 0);

  glDisable(GL_BLEND);
#ifndef __EMSCRIPTEN__
  glDisable(GL_LINE_SMOOTH);
  glDisable(GL_POLYGON_SMOOTH);
#endif
  glLineWidth(1.0);
}

void ModeButton::cleanup()
{
  delete[] vertices;
  delete[] vertex_colors;
  delete[] vertex_indices;

  glDeleteProgram(plugin_program);
  glDeleteBuffers(1, &plugin_vertices);
  glDeleteBuffers(1, &plugin_vertex_colors);
  glDeleteBuffers(1, &plugin_vertex_indices);
}

void HeatButton::setup()
{
  {
    const unsigned long n_theta = 64;
    const unsigned long n_r = 16;
    const unsigned long n_triangles = n_theta + (n_r - 1)*n_theta*2;
    const unsigned long n_vertices = 1 + n_r * n_theta; // One for the central vertex

    vertices = new GLfloat[ n_vertices * coordinates_per_vertex ];
    vertex_colors = new GLfloat[ n_vertices * colors_per_vertex ];
    vertex_indices = new GLuint[ n_triangles * vertices_per_triangle ];

    //one vertex at the origin
    vertices[0] = heat_button_x;
    vertices[1] = heat_button_y;

    unsigned long n = 1, v = 2, i=0; //start at the next vertex index
    // Set up triangles for the inner fan
    for (unsigned long nt = 0; nt < n_theta; ++nt, ++n)
    {
      const float r = 1.0f/n_r * button_radius;
      const float theta = float(nt) / n_theta * M_PI*2.0f;
      vertices[v + 0] = heat_button_x + r * std::cos(theta);
      vertices[v + 1] = heat_button_y + r * std::sin(theta);

      vertex_indices[i + 0] = 0;
      vertex_indices[i + 1] = n;
      vertex_indices[i + 2] = (nt == n_theta - 1 ? 1 : n + 1);

      v += coordinates_per_vertex;
      i += vertices_per_triangle;
    }
    // Set up the outer triangles
    for (unsigned long nr = 0; nr < n_r-1; ++nr)
      for (unsigned long nt = 0; nt < n_theta; ++nt, ++n)
      {
        const float r = float(nr+1)/n_r * button_radius;
        const float theta = float(nt) / n_theta * M_PI*2.0f;
        vertices[v + 0] = heat_button_x + r * std::cos(theta);
        vertices[v + 1] = heat_button_y + r * std::sin(theta);

        vertex_indices[i + 0] = n - n_theta;
        vertex_indices[i + 1] = n;
        vertex_indices[i + 2] = (nt == n_theta - 1 ? (n - n_theta + 1) : n + 1);
        vertex_indices[i + 3] = n - n_theta;
        vertex_indices[i + 4] = (nt == n_theta - 1 ? (n - n_theta + 1) : n + 1);
        vertex_indices[i + 5] = (nt == n_theta - 1 ? (n - 2*n_theta+1) : n - n_theta + 1);
        
        v += coordinates_per_vertex;
        i += 2*vertices_per_triangle;
      }

    glGenBuffers(1, &plugin_vertex_indices);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, plugin_vertex_indices);
    glBufferData(
        GL_ELEMENT_ARRAY_BUFFER,
        sizeof(GLuint)*(n_triangles*vertices_per_triangle),
        vertex_indices,
        GL_STATIC_DRAW
    );

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

void HeatButton::draw()
{
  const unsigned long n_theta = 64;
  const unsigned long n_r = 16;
  const unsigned long n_triangles = n_theta + (n_r - 1)*n_theta*2;
  const unsigned long n_vertices = 1 + n_r * n_theta; // One for the central vertex

  // Update the button color depending upon whether we are in hot or cold mode
  color hot_color = hot(0.75);
  color cold_color = hot(0.25);
  color disabled_color = { 0.5, 0.5, 0.5 };

  //one vertex at the origin
  const float t_center = alt_press ? 0.0 : 1.0;
  const color center_color = hot(t_center);
  vertex_colors[0] = center_color.R;
  vertex_colors[1] = center_color.G;
  vertex_colors[2] = center_color.B;

  unsigned long c=3; //start at the next vertex index
  for (unsigned long n = 1; n < n_vertices; ++n)
  {
    const unsigned int nr = (n-1)/n_r;
    const float r = float(nr)/n_r;
    const float temp = alt_press ? 0.5 - 0.5*std::exp(-r*r/2.0) : 0.5 + 0.5*std::exp(-r*r/2.0);
    const color vertex_color = hot(temp);
    vertex_colors[c + 0] = vertex_color.R;
    vertex_colors[c + 1] = vertex_color.G;
    vertex_colors[c + 2] = vertex_color.B;

    c += colors_per_vertex;
  }

  glEnable(GL_BLEND);
#ifndef __EMSCRIPTEN__
  glEnable(GL_POLYGON_SMOOTH);
#endif
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  glUseProgram(plugin_program);

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, plugin_vertex_indices);

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
  glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*n_vertices*colors_per_vertex, vertex_colors, GL_STATIC_DRAW);
  glVertexAttribPointer(
    plugin_attribute_v_color, // attribute
    3,                 // number of elements per vertex (r,g,b)
    GL_FLOAT,          // the type of each element
    GL_FALSE,          // take our values as-is
    0,                 // no extra data between each position
    0);

  glDrawElements(GL_TRIANGLES, n_triangles*vertices_per_triangle, GL_UNSIGNED_INT, 0);

  glDisableVertexAttribArray(plugin_attribute_coord2d);
  glDisableVertexAttribArray(plugin_attribute_v_color);

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
  glBindBuffer(GL_ARRAY_BUFFER, 0);

  glDisable(GL_BLEND);
#ifndef __EMSCRIPTEN__
  glDisable(GL_POLYGON_SMOOTH);
#endif
}

void HeatButton::cleanup()
{
  delete[] vertices;
  delete[] vertex_colors;
  delete[] vertex_indices;

  glDeleteProgram(plugin_program);
  glDeleteBuffers(1, &plugin_vertices);
  glDeleteBuffers(1, &plugin_vertex_colors);
  glDeleteBuffers(1, &plugin_vertex_indices);
}
