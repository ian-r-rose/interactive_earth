#ifndef RENDERING_H
#define RENDERING_H

#include "GL/glew.h"
#include "SDL2/SDL.h"
#include "SDL2/SDL_opengl.h"
#include "convection.h"

//Renders the temperature (or composition) field of the
//mantle, perturbed by any seismic displacement.
class Mantle
{
  public:
    Mantle( const ConvectionSimulator &sim ) : sim(sim) {}
    void setup();  //Setup OpenGL data structures
    void draw(bool draw_composition);  //Render using OpenGL
    void cleanup();  //Cleanup OpenGL data structures

  private:
    const ConvectionSimulator &sim;

    //Data for rendering with OpenGL
    GLfloat* vertices;
    GLfloat* vertex_colors;
    GLuint* triangle_vertex_indices;

    GLuint program;
    GLuint vbo_vertices;
    GLuint vbo_colors;
    GLuint ibo_triangle_vertex_indices;
    GLint attribute_coord2d;
    GLint attribute_v_color;
};

#endif
