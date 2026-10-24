/************************************
  OpenGL headers for each platform,
  and the GLSL version preludes that
  go with them. We target GL 3.3 core
  on the desktop and OpenGL ES 3.0
  (WebGL2) in the browser.
************************************/
#ifndef GL_INCLUDES_H
#define GL_INCLUDES_H

#if defined(__EMSCRIPTEN__)
#include <GLES3/gl3.h>
#define IE_GLSL_VERSION "#version 300 es\nprecision mediump float;\n"
#elif defined(__APPLE__)
#include <OpenGL/gl3.h>
#define IE_GLSL_VERSION "#version 330 core\n"
#else
#define GL_GLEXT_PROTOTYPES
#include <GL/gl.h>
#include <GL/glext.h>
#define IE_GLSL_VERSION "#version 330 core\n"
#endif

//Prepend to vertex shaders. Inputs are declared with `in`, outputs with `out`.
#define IE_GLSL_VERTEX_PRELUDE IE_GLSL_VERSION
//Prepend to fragment shaders, which write their color to `frag_color`.
#define IE_GLSL_FRAGMENT_PRELUDE IE_GLSL_VERSION "out vec4 frag_color;\n"

#endif
