/************************************
  Small helpers for working with
  OpenGL shaders.
************************************/
#ifndef GL_UTIL_H
#define GL_UTIL_H

#include "gl_includes.h"

//Compile and link a shader program from vertex and fragment shader bodies,
//prepending the GLSL version preludes from gl_includes.h. Compile and link
//errors are printed with the driver's log. Returns 0 on failure.
GLuint create_program(const char *vertex_source, const char *fragment_source, const char *name);

#endif
