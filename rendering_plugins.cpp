#include <cmath>
#include "rendering_plugins.h"
#include "color.h"

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
    
    color core_color = hot(1.0);

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

  const float semimajor_axis_angle = (use_geographic_frame ? sim.spin_angle() + M_PI/2. : 0.0f);
  const float a = 1.0f;
  const float b = 1.0f-flattening;

  unsigned long v = 2; //start at the next vertex index after the one in the center
  for (unsigned long n = 0; n < n_triangles; ++n)
  {
    const float theta = dtheta*n;
    //Vertices of the inner ellipse
    vertices[v + 0] = r_inner*(a+b)/2. * std::cos(theta) +
                      r_inner*(a-b)/2. * std::cos(-theta+ 2.*(semimajor_axis_angle));
    vertices[v + 1] = r_inner*(a+b)/2. * std::sin(theta) +
                      r_inner*(a-b)/2. * std::sin(-theta+ 2.*(semimajor_axis_angle));

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
  const float theta = (use_geographic_frame ? sim.spin_angle() : M_PI/2.);
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
  glEnable(GL_POLYGON_SMOOTH);
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

  glDisable(GL_BLEND);
  glDisable(GL_POLYGON_SMOOTH);
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
      vertex_colors[(n_vertices-i)*(colors_per_vertex+1) + 0] = 0.2;//line_color.R;
      vertex_colors[(n_vertices-i)*(colors_per_vertex+1)+ 1] = 0.2;//line_color.G;
      vertex_colors[(n_vertices-i)*(colors_per_vertex+1) + 2] = 0.2;//line_color.B;
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
  theta -= (use_geographic_frame ? 0.0 : sim.spin_angle()+M_PI/2.);
  r += r_inner;

  //Information about ellipticity
  const float semimajor_axis_angle = (use_geographic_frame ? sim.spin_angle() + M_PI/2. : 0.);
  const float a = 1.0f;
  const float b = 1.0f-flattening;

  //Position of the seismometer in (possibly) elliptical space 
  const float seis_x = r*(a+b)/2. * std::cos(theta) +
                       r*(a-b)/2. * std::cos(-theta+ 2.*(semimajor_axis_angle));
  const float seis_y = r*(a+b)/2. * std::sin(theta) +
                       r*(a-b)/2. * std::sin(-theta+ 2.*(semimajor_axis_angle));

  //Seismometer vertices
  vertices[ (n_lines+1)*coordinates_per_vertex + 0] = seis_x-0.02;
  vertices[ (n_lines+1)*coordinates_per_vertex + 1] = seis_y+0.02;
  vertices[ (n_lines+2)*coordinates_per_vertex + 0] = seis_x+0.02;
  vertices[ (n_lines+2)*coordinates_per_vertex + 1] = seis_y+0.02;
  vertices[ (n_lines+3)*coordinates_per_vertex + 0] = seis_x;
  vertices[ (n_lines+3)*coordinates_per_vertex + 1] = seis_y-0.02;

  glEnable(GL_BLEND);
  glEnable(GL_LINE_SMOOTH);
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
    4,                 // number of elements per verte (r,g,b,a)
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
  glDisable(GL_LINE_SMOOTH);
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

