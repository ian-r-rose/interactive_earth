#include <fftw3.h>
#include "staggered_grid.h"
#include "GL/glew.h"
#include "SDL2/SDL.h"
#include "SDL2/SDL_opengl.h"

#ifndef STOKES_H
#define STOKES_H

class StokesSolver
{
  private:
    
    int nx,ny;
    int ncells;
    double Ra;
    double dt;
    double theta;
    double gamma;

    double heat_source_radius;
    double heat_source;

    StaggeredGrid grid;
    

    double *T;
    double *scratch1, *scratch2;
    double *g, *lux, *luy;
    double *freqs;
    double *vorticity;
    double *stream;
    double *curl_T;
    double *u;
    double *v;

    fftw_plan dst, idst, dft, idft;
    fftw_complex* curl_T_spectral;

    GLfloat* triangle_vertices;  
    GLfloat* vertex_colors;
    GLint* vertex_indices;

    //workhorse functions
    void initialize_temperature();
    double heat(const Point&, const Point&);
    void setup_stokes_problem();
    void setup_diffusion_problem();
    void assemble_curl_T_vector();
   
    //functions for evaluating field at points
    double initial_temperature(const Point&);
    double temperature(const Point&);
    Point velocity(const Point&);

  public:
    StokesSolver( double lx, double ly, int nx, int ny, double Rayleigh);
    ~StokesSolver();

    double rayleigh_number() const;
    double timescale() const;

    void add_heat(double x, double y, bool hot);
    void semi_lagrangian_advect();
    void diffuse_temperature();
    void solve_stokes();
    void draw();
    double nusselt_number();
    void update_state(double rayleigh, double gravity_angle = 0);
    int setup_opengl();
    int cleanup_opengl();
};

#endif
