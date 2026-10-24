#ifndef RENDERING_H
#define RENDERING_H

#include <vector>
#include "gl_includes.h"
#include "color.h"
#include "convection.h"

//Renders the temperature (or composition) field of the
//mantle, perturbed by any seismic displacement.
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
  GLuint field_texture = 0;    //(field, displacement) per grid cell
  GLuint colormap_texture = 0; //256 entry lookup table
  GLint uniform_field = -1, uniform_colormap = -1;
  GLint uniform_inverse_shape = -1, uniform_r_inner = -1, uniform_grid_size = -1;

  std::vector<float> staging; //CPU copy of the field texture
  color (*uploaded_colormap)(double) = nullptr;
};

#endif
