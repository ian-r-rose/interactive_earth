#ifndef RENDERING_H
#define RENDERING_H

#include <vector>
#include "gl_includes.h"
#include "color.h"
#include "convection.h"

//The domain is rendered flattened along the spin axis, as if under
//centrifugal forces: positions are R(angle) diag(1, 1 - flattening) R(-angle)
//applied to the circular ones, with angle = spin_angle + pi/2. Returns that
//matrix (or its inverse) in column-major order.
void flattening_shape(double spin_angle, bool inverse, double shape[4]);

//Renders the temperature (or composition) field of the mantle, perturbed by
//any seismic displacement, and the core inside it.
//
//The field is uploaded each frame as a texture indexed by (theta, r) cell,
//and a fragment shader maps every pixel back to polar coordinates, samples
//it, and applies the colormap. No geometry is needed beyond one triangle
//covering the viewport.
class Mantle
{
public:
  Mantle(const ConvectionSimulator &sim) : sim(sim) {}
  void setup();                     //Setup OpenGL data structures
  void draw(bool draw_composition); //Render using OpenGL
  void cleanup();                   //Cleanup OpenGL data structures

private:
  void upload_colormap();

  const ConvectionSimulator &sim;

  GLuint program = 0;
  GLuint vertex_array = 0;     //Empty; the geometry comes from the vertex index
  GLuint field_texture = 0;    //(field, displacement) per grid cell
  GLuint colormap_texture = 0; //256 entry lookup table
  GLint uniform_field = -1, uniform_colormap = -1;
  GLint uniform_inverse_shape = -1, uniform_r_inner = -1, uniform_grid_size = -1;

  std::vector<float> staging; //CPU copy of the field texture
  color (*uploaded_colormap)(double) = nullptr;
};

#endif
