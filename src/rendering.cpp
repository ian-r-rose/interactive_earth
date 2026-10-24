#include <cmath>
#include "rendering.h"
#include "gl_util.h"
#include "profiler.h"

//get access to some of the globals
extern const double flattening;
extern color (*colormap)(double);

//Number of entries in the colormap lookup table, matching color.cpp
static const int colormap_size = 256;

//A single triangle that covers the whole viewport, generated from the
//vertex index, so no vertex buffers are needed.
static const char *vertex_source = R"(
out vec2 position;
void main()
{
  vec2 corner = vec2(float((gl_VertexID << 1) & 2), float(gl_VertexID & 2));
  position = 2.0 * corner - 1.0;  // (-1,-1), (3,-1), (-1,3)
  gl_Position = vec4(position, 0.0, 1.0);
}
)";

//Map each pixel back to polar coordinates, sample the field there, and
//color it. The outer radius is 1 in normalized device coordinates.
static const char *fragment_source = R"(
in vec2 position;
uniform sampler2D field;          // (T or C, displacement), one texel per cell
uniform sampler2D colormap;       // lookup table, one row
uniform mat2 inverse_shape;       // undoes the rendered flattening
uniform float r_inner;
uniform vec2 grid_size;           // (ntheta, nr)

const float two_pi = 6.28318530718;

void main()
{
  vec2 p = inverse_shape * position;
  float r = length(p);
  if (r > 1.0 || r < r_inner)
    discard;
  float theta = atan(p.y, p.x);
  if (theta < 0.0)
    theta += two_pi;

  // Cell (i, j) sits at theta = i dtheta and r = r_inner + j dr; put those
  // at texel centers. The texture repeats in theta, which handles the seam.
  float u = theta / two_pi + 0.5 / grid_size.x;
  float v = ((r - r_inner) / (1.0 - r_inner) * (grid_size.y - 1.0) + 0.5) / grid_size.y;
  vec2 value = texture(field, vec2(u, v)).rg;

  // Same lookup as color.cpp: nearest entry to value * size, clamped
  int size = textureSize(colormap, 0).x;
  int index = clamp(int(floor(value.r * float(size) + 0.5)), 0, size - 1);
  vec3 rgb = texelFetch(colormap, ivec2(index, 0), 0).rgb;

  // Seismic displacement brightens or darkens the field
  frag_color = vec4(rgb + value.g, 1.0);
}
)";

void Mantle::setup()
{
  program = create_program(vertex_source, fragment_source, "Mantle");
  uniform_field = glGetUniformLocation(program, "field");
  uniform_colormap = glGetUniformLocation(program, "colormap");
  uniform_inverse_shape = glGetUniformLocation(program, "inverse_shape");
  uniform_r_inner = glGetUniformLocation(program, "r_inner");
  uniform_grid_size = glGetUniformLocation(program, "grid_size");

  const RegularGrid &grid = sim.get_grid();
  staging.resize(2 * grid.ncells);

  //Field texture: half floats are plenty for colors, and are filterable
  //everywhere (including WebGL2), unlike 32 bit floats.
  glGenTextures(1, &field_texture);
  glBindTexture(GL_TEXTURE_2D, field_texture);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RG16F, grid.ntheta, grid.nr, 0, GL_RG, GL_FLOAT, nullptr);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

  //Colormap lookup table, read with texelFetch
  glGenTextures(1, &colormap_texture);
  glBindTexture(GL_TEXTURE_2D, colormap_texture);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, colormap_size, 1, 0, GL_RGBA, GL_FLOAT, nullptr);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glBindTexture(GL_TEXTURE_2D, 0);
  upload_colormap();
}

//Sample the current colormap at the points color.cpp uses for each entry
void Mantle::upload_colormap()
{
  std::vector<float> table(4 * colormap_size);
  for (int i = 0; i < colormap_size; ++i)
  {
    const color c = colormap(double(i) / colormap_size);
    table[4 * i + 0] = c.R;
    table[4 * i + 1] = c.G;
    table[4 * i + 2] = c.B;
    table[4 * i + 3] = 1.0f;
  }
  glBindTexture(GL_TEXTURE_2D, colormap_texture);
  glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, colormap_size, 1, GL_RGBA, GL_FLOAT, table.data());
  glBindTexture(GL_TEXTURE_2D, 0);
  uploaded_colormap = colormap;
}

void Mantle::cleanup()
{
  glDeleteProgram(program);
  glDeleteTextures(1, &field_texture);
  glDeleteTextures(1, &colormap_texture);
}

void Mantle::draw(bool draw_composition)
{
  const RegularGrid &grid = sim.get_grid();

  //Interleave the field and the displacement for upload
  {
    IE_PROFILE(render_prep);
    const double *field = draw_composition ? sim.composition() : sim.temperature();
    const double *D = sim.displacement();
    for (int k = 0; k < grid.ncells; ++k)
    {
      staging[2 * k + 0] = float(field[k]);
      staging[2 * k + 1] = float(D[k]);
    }
  }

  IE_PROFILE(render_gl);
  if (colormap != uploaded_colormap)
    upload_colormap();

  glBindTexture(GL_TEXTURE_2D, field_texture);
  glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, grid.ntheta, grid.nr, GL_RG, GL_FLOAT, staging.data());

  //The domain is rendered flattened along the spin axis, as if under
  //centrifugal forces: positions are R(angle) diag(a, b) R(-angle) applied
  //to the circular ones. The shader needs the inverse of that.
  const double angle = sim.spin_angle() + M_PI / 2.;
  const double a = 1.0, b = 1.0 - flattening;
  const double c = std::cos(angle), s = std::sin(angle);
  const double ia = 1. / a, ib = 1. / b;
  //Symmetric, so column-major order is the same as row-major
  const GLfloat inverse_shape[4] = {float(c * c * ia + s * s * ib), float(c * s * (ia - ib)),
                                    float(c * s * (ia - ib)), float(s * s * ia + c * c * ib)};

  glClearColor(0.0, 0.0, 0.0, 1.0);
  glClear(GL_COLOR_BUFFER_BIT);

  glUseProgram(program);
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, field_texture);
  glActiveTexture(GL_TEXTURE1);
  glBindTexture(GL_TEXTURE_2D, colormap_texture);
  glUniform1i(uniform_field, 0);
  glUniform1i(uniform_colormap, 1);
  glUniformMatrix2fv(uniform_inverse_shape, 1, GL_FALSE, inverse_shape);
  glUniform1f(uniform_r_inner, float(grid.r_inner));
  glUniform2f(uniform_grid_size, float(grid.ntheta), float(grid.nr));

  glDrawArrays(GL_TRIANGLES, 0, 3);

  glActiveTexture(GL_TEXTURE1);
  glBindTexture(GL_TEXTURE_2D, 0);
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, 0);
}
