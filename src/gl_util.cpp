#include "gl_util.h"
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
