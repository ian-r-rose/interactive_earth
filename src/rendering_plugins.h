#ifndef RENDERING_PLUGINS_H
#define RENDERING_PLUGINS_H

#include <vector>
#include "gl_util.h"
#include "rendering.h"

extern const double r_inner;
extern const double flattening;

/******************************************
   Rendering plugins draw on top of the
   mantle, mostly in the space taken by the
   core in the center of the domain, to
   show things of interest or supporting
   data.
******************************************/

class RenderingPlugin
{
public:
  RenderingPlugin(const ConvectionSimulator &sim) : sim(sim) {}

protected:
  const ConvectionSimulator &sim;
};

//The spin axis, for true polar wander
class Axis : RenderingPlugin
{
public:
  Axis(const ConvectionSimulator &sim) : RenderingPlugin(sim) {}
  void setup();
  void draw();
  void cleanup();

private:
  ColoredMesh mesh;
};

//A scrolling trace of the seismometer reading, and a
//marker for where the seismometer is
class Seismograph : RenderingPlugin
{
public:
  Seismograph(const ConvectionSimulator &sim) : RenderingPlugin(sim) {}
  void setup();
  //point_size is the size of one (logical) screen pixel in clip coordinates
  void draw(float point_size);
  void cleanup();

  void clear_record();

private:
  static const int n_points = 1001; //Length of the record
  //The trace is a strip of triangles, two vertices per point,
  //followed by the seismometer marker
  static const int n_vertices = 2 * n_points + 3;

  ColoredMesh mesh;
  std::vector<float> record;     //Recent readings, oldest first
  std::vector<GLfloat> vertices; //(x, y) per vertex
};

#endif
