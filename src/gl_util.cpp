#include "gl_util.h"
#include <cmath>
#include <cstdio>
#include <vector>

static GLuint compile_shader(GLenum type, const char *prelude, const char *body, const char *name)
{
  const char *sources[] = {prelude, body};
  GLuint shader = glCreateShader(type);
  glShaderSource(shader, 2, sources, nullptr);
  glCompileShader(shader);

  GLint ok = GL_FALSE;
  glGetShaderiv(shader, GL_COMPILE_STATUS, &ok);
  if (!ok)
  {
    GLint length = 0;
    glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &length);
    std::vector<char> log(length > 1 ? length : 1, '\0');
    glGetShaderInfoLog(shader, log.size(), nullptr, log.data());
    std::fprintf(stderr, "%s: %s shader failed to compile:\n%s\n", name,
                 type == GL_VERTEX_SHADER ? "vertex" : "fragment", log.data());
    glDeleteShader(shader);
    return 0;
  }
  return shader;
}

GLuint create_program(const char *vertex_source, const char *fragment_source, const char *name)
{
  GLuint vs = compile_shader(GL_VERTEX_SHADER, IE_GLSL_VERTEX_PRELUDE, vertex_source, name);
  GLuint fs = compile_shader(GL_FRAGMENT_SHADER, IE_GLSL_FRAGMENT_PRELUDE, fragment_source, name);
  if (!vs || !fs)
  {
    glDeleteShader(vs);
    glDeleteShader(fs);
    return 0;
  }

  GLuint program = glCreateProgram();
  glAttachShader(program, vs);
  glAttachShader(program, fs);
  glLinkProgram(program);
  //The program keeps what it needs; the shader objects can go
  glDeleteShader(vs);
  glDeleteShader(fs);

  GLint ok = GL_FALSE;
  glGetProgramiv(program, GL_LINK_STATUS, &ok);
  if (!ok)
  {
    GLint length = 0;
    glGetProgramiv(program, GL_INFO_LOG_LENGTH, &length);
    std::vector<char> log(length > 1 ? length : 1, '\0');
    glGetProgramInfoLog(program, log.size(), nullptr, log.data());
    std::fprintf(stderr, "%s: shader program failed to link:\n%s\n", name, log.data());
    glDeleteProgram(program);
    return 0;
  }
  return program;
}

static const char *mesh_vertex_source = R"(
layout(location = 0) in vec2 position;
layout(location = 1) in vec4 color;
uniform mat2 transform;
out vec4 vertex_color;
void main()
{
  vertex_color = color;
  gl_Position = vec4(transform * position, 0.0, 1.0);
}
)";

static const char *mesh_fragment_source = R"(
in vec4 vertex_color;
void main()
{
  frag_color = vertex_color;
}
)";

void ColoredMesh::setup(int n_vertices, const GLfloat *positions, const GLfloat *colors,
                        int n_indices, const GLuint *indices, const char *name)
{
  this->n_vertices = n_vertices;
  program = create_program(mesh_vertex_source, mesh_fragment_source, name);
  uniform_transform = glGetUniformLocation(program, "transform");

  glGenVertexArrays(1, &vertex_array);
  glBindVertexArray(vertex_array);

  glGenBuffers(1, &position_buffer);
  glBindBuffer(GL_ARRAY_BUFFER, position_buffer);
  glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat) * 2 * n_vertices, positions, GL_DYNAMIC_DRAW);
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, nullptr);

  glGenBuffers(1, &color_buffer);
  glBindBuffer(GL_ARRAY_BUFFER, color_buffer);
  glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat) * 4 * n_vertices, colors, GL_STATIC_DRAW);
  glEnableVertexAttribArray(1);
  glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, 0, nullptr);

  glGenBuffers(1, &index_buffer);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_buffer);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint) * n_indices, indices, GL_STATIC_DRAW);

  //The vertex array remembers the element buffer, so unbind it first
  glBindVertexArray(0);
  glBindBuffer(GL_ARRAY_BUFFER, 0);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}

void ColoredMesh::update_positions(const GLfloat *positions)
{
  glBindBuffer(GL_ARRAY_BUFFER, position_buffer);
  glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(GLfloat) * 2 * n_vertices, positions);
  glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void ColoredMesh::draw(GLenum mode, int first, int count, const GLfloat transform[4])
{
  glUseProgram(program);
  glUniformMatrix2fv(uniform_transform, 1, GL_FALSE, transform);
  glBindVertexArray(vertex_array);
  //GL takes the byte offset into the index buffer as a pointer
  const size_t offset = sizeof(GLuint) * first;
  glDrawElements(mode, count, GL_UNSIGNED_INT, (void *)offset); // NOLINT(performance-no-int-to-ptr)
  glBindVertexArray(0);
}

void ColoredMesh::cleanup()
{
  glDeleteProgram(program);
  glDeleteBuffers(1, &position_buffer);
  glDeleteBuffers(1, &color_buffer);
  glDeleteBuffers(1, &index_buffer);
  glDeleteVertexArrays(1, &vertex_array);
}

void identity_transform(GLfloat transform[4])
{
  transform[0] = 1.f;
  transform[1] = 0.f;
  transform[2] = 0.f;
  transform[3] = 1.f;
}

void rotation_transform(double angle, GLfloat transform[4])
{
  transform[0] = float(std::cos(angle));
  transform[1] = float(std::sin(angle));
  transform[2] = -float(std::sin(angle));
  transform[3] = float(std::cos(angle));
}
