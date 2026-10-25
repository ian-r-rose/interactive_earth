#include <algorithm>
#include <cmath>
#include "rendering_plugins.h"

/* Axis */

void Axis::setup()
{
  //Geometry for a spin angle of zero; draw() rotates it into place.
  const double r = 0.8 * r_inner * (1.0 - flattening); //Half length of the axis
  const double arrowhead = 0.1 * r_inner;
  const double dtheta = 0.5 / r_inner * M_PI / 180.0; //Angular half width of the axis

  const GLfloat positions[] = {//Axis, a thin quad
                               float(r * std::cos(-dtheta)), float(r * std::sin(-dtheta)),   //
                               float(-r * std::cos(dtheta)), float(-r * std::sin(dtheta)),   //
                               float(-r * std::cos(-dtheta)), float(-r * std::sin(-dtheta)), //
                               float(r * std::cos(dtheta)), float(r * std::sin(dtheta)),     //
                               //Arrowheads at either end
                               float(r + arrowhead), 0.f, //
                               float((r - arrowhead) * std::cos(4. * dtheta)),
                               float((r - arrowhead) * std::sin(4. * dtheta)),
                               float((r - arrowhead) * std::cos(-4. * dtheta)),
                               float((r - arrowhead) * std::sin(-4. * dtheta)),
                               float(-(r + arrowhead)), 0.f, //
                               float(-(r - arrowhead) * std::cos(4. * dtheta)),
                               float(-(r - arrowhead) * std::sin(4. * dtheta)),
                               float(-(r - arrowhead) * std::cos(-4. * dtheta)),
                               float(-(r - arrowhead) * std::sin(-4. * dtheta))};
  const int n_vertices = sizeof(positions) / sizeof(positions[0]) / 2;

  std::vector<GLfloat> colors;
  for (int i = 0; i < n_vertices; ++i)
    colors.insert(colors.end(), {0.15f, 0.15f, 0.15f, 1.0f});

  const GLuint indices[] = {0, 1, 2, 0, 2, 3, 4, 5, 6, 7, 8, 9};
  mesh.setup(n_vertices, positions, colors.data(), 12, indices, "Axis");
}

void Axis::draw()
{
  GLfloat transform[4];
  rotation_transform(sim.spin_angle(), transform);
  mesh.draw(GL_TRIANGLES, 0, 12, transform);
}

void Axis::cleanup()
{
  mesh.cleanup();
}

/* Seismograph */

void Seismograph::clear_record()
{
  std::fill(record.begin(), record.end(), 0.f);
}

void Seismograph::setup()
{
  record.assign(n_points, 0.f);
  vertices.assign(2 * n_vertices, 0.f);

  //The trace is black, fading out at either end. The seismometer is blue.
  std::vector<GLfloat> colors;
  for (int i = 0; i < n_points; ++i)
  {
    float alpha = 1.0f;
    if (i < n_points / 6)
      alpha = 6.f * float(i) / float(n_points);
    else if (i > 5 * n_points / 6)
      alpha = 6.f * float(n_points - i) / float(n_points);
    for (int side = 0; side < 2; ++side)
      colors.insert(colors.end(), {0.f, 0.f, 0.f, alpha});
  }
  for (int i = 0; i < 3; ++i)
    colors.insert(colors.end(), {0.f, 0.3f, 1.0f, 0.8f});

  //Two triangles per segment of the trace, then the seismometer
  std::vector<GLuint> indices;
  for (GLuint i = 0; i + 1 < GLuint(n_points); ++i)
    indices.insert(indices.end(), {2 * i, 2 * i + 1, 2 * i + 2, 2 * i + 1, 2 * i + 3, 2 * i + 2});
  for (GLuint i = 0; i < 3; ++i)
    indices.push_back(2 * n_points + i);

  mesh.setup(n_vertices, vertices.data(), colors.data(), int(indices.size()), indices.data(),
             "Seismograph");
}

void Seismograph::draw(float point_size)
{
  //Scroll the record and add the newest reading
  std::rotate(record.begin(), record.begin() + 1, record.end());
  record.back() = float(sim.seismometer_reading() * 0.7 * r_inner);

  //Lines wider than one pixel are not available in core GL or WebGL,
  //so draw the trace as a strip of triangles, offset to either side of
  //the line, two points wide.
  const double half_width = 1.0 * point_size;
  const double width = 1.8 * r_inner * (1.0 - flattening);
  auto x = [&](int i) { return -0.5 * width + width * i / (n_points - 1); };
  for (int i = 0; i < n_points; ++i)
  {
    const int prev = (i > 0 ? i - 1 : i), next = (i + 1 < n_points ? i + 1 : i);
    const double tx = x(next) - x(prev), ty = record[next] - record[prev];
    const double length = std::sqrt(tx * tx + ty * ty);
    const double nx = -ty / length * half_width, ny = tx / length * half_width;
    vertices[4 * i + 0] = float(x(i) + nx);
    vertices[4 * i + 1] = float(record[i] + ny);
    vertices[4 * i + 2] = float(x(i) - nx);
    vertices[4 * i + 3] = float(record[i] - ny);
  }

  //The seismometer, at its position in the (possibly flattened) domain
  double theta, r;
  sim.seismometer_position(theta, r);
  r += r_inner;
  double shape[4];
  flattening_shape(sim.spin_angle(), false, shape);
  const double px = r * std::cos(theta), py = r * std::sin(theta);
  const double sx = shape[0] * px + shape[2] * py, sy = shape[1] * px + shape[3] * py;
  GLfloat *marker = &vertices[4 * n_points];
  marker[0] = float(sx - 0.04);
  marker[1] = float(sy + 0.04);
  marker[2] = float(sx + 0.04);
  marker[3] = float(sy + 0.04);
  marker[4] = float(sx);
  marker[5] = float(sy - 0.04);

  mesh.update_positions(vertices.data());

  GLfloat transform[4];
  identity_transform(transform);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  mesh.draw(GL_TRIANGLES, 0, 6 * (n_points - 1) + 3, transform);
  glDisable(GL_BLEND);
}

void Seismograph::cleanup()
{
  mesh.cleanup();
}
