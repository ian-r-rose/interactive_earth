#ifndef RENDERING_PLUGINS_H
#define RENDERING_PLUGINS_H

#include "convection.h"

extern const double r_inner;
extern const unsigned int nr, ntheta;

//Useful constants
static const short vertices_per_triangle = 3;
static const short vertices_per_line = 2;
static const short coordinates_per_vertex = 2;
static const short colors_per_vertex = 3;
static const short triangles_per_quad = 2;

/******************************************
   Interface for rendering plugins.
   There is space in the center of the
   domain which could be used to display
   various things of interest/supporting
   data. Define a basic plugin architecture
   for accomplishing this.
******************************************/

class RenderingPlugin
{
  public:
    RenderingPlugin( const ConvectionSimulator &sim ) : sim(sim) {}
    void setup();
    void draw();
    void cleanup();

  protected:
    const ConvectionSimulator &sim;
};


class Core : RenderingPlugin
{
  public:
    Core( const ConvectionSimulator &sim ) : RenderingPlugin(sim) {}
    void setup();
    void draw();
    void cleanup();
  private:
    //Data for rendering with OpenGL
    GLfloat* vertices;
    GLfloat* vertex_colors;
    GLuint* triangle_vertex_indices;

    GLuint plugin_program;
    GLuint plugin_vertices;
    GLuint plugin_vertex_colors;
    GLuint plugin_triangle_vertex_indices;
    GLint plugin_attribute_coord2d;
    GLint plugin_attribute_v_color;
};

class Axis : RenderingPlugin
{
  public:
    Axis( const ConvectionSimulator &sim ) : RenderingPlugin(sim) {}
    void setup();
    void draw();
    void cleanup();
  private:
    //Data for rendering with OpenGL
    GLfloat* vertices;
    GLfloat* vertex_colors;
    GLuint* triangle_vertex_indices;

    GLuint plugin_program;
    GLuint plugin_vertices;
    GLuint plugin_vertex_colors;
    GLuint plugin_triangle_vertex_indices;
    GLint plugin_attribute_coord2d;
    GLint plugin_attribute_v_color;
};

class Seismograph : RenderingPlugin
{
  public:
    Seismograph( const ConvectionSimulator &sim ) : RenderingPlugin(sim),
                                                    n_lines(500),
                                                    n_vertices(n_lines+1 /*endpoint*/ + vertices_per_triangle /*seismometer*/)
                                                    {}
    void setup();
    void draw();
    void cleanup();

    void clear_record();

  private:
    const unsigned int n_lines;
    const unsigned int n_vertices;

    //Data for rendering with OpenGL
    GLfloat* vertices;
    GLfloat* vertex_colors;
    GLuint* line_vertex_indices;

    GLuint plugin_program;
    GLuint plugin_vertices;
    GLuint plugin_vertex_colors;
    GLuint plugin_line_vertex_indices;
    GLint plugin_attribute_coord2d;
    GLint plugin_attribute_v_color;
};

#endif
