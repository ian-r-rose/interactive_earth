/************************************
  Small helpers for working with
  OpenGL shaders and meshes.
************************************/
#ifndef GL_UTIL_H
#define GL_UTIL_H

#include "gl_includes.h"

//Compile and link a shader program from vertex and fragment shader bodies,
//prepending the GLSL version preludes from gl_includes.h. Compile and link
//errors are printed with the driver's log. Returns 0 on failure.
GLuint create_program(const char *vertex_source, const char *fragment_source, const char *name);

//Indexed 2D geometry with a color per vertex, drawn with a 2x2 transform
//(e.g. a rotation) applied to the positions. Used for the overlays drawn on
//top of the mantle.
class ColoredMesh
{
public:
  //positions: (x, y) per vertex. colors: (r, g, b, a) per vertex.
  void setup(int n_vertices, const GLfloat *positions, const GLfloat *colors, int n_indices,
             const GLuint *indices, const char *name);
  //Replace all of the positions, e.g. for geometry that changes every frame
  void update_positions(const GLfloat *positions);
  //Draw `count` indices starting at `first`. transform is a column-major 2x2 matrix.
  void draw(GLenum mode, int first, int count, const GLfloat transform[4]);
  void cleanup();

private:
  int n_vertices = 0;
  GLuint program = 0;
  GLuint vertex_array = 0, position_buffer = 0, color_buffer = 0, index_buffer = 0;
  GLint uniform_transform = -1;
};

//Column-major 2x2 matrices for ColoredMesh::draw
void identity_transform(GLfloat transform[4]);
void rotation_transform(double angle, GLfloat transform[4]);

#endif
