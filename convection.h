#include <complex>
#include "fftw3.h"
#include "staggered_grid.h"
#include "GL/glew.h"
#include "SDL2/SDL.h"
#include "SDL2/SDL_opengl.h"

#ifndef STOKES_H
#define STOKES_H

class ConvectionSimulator
{
  private:
   
    //Information about the grid
    int nx,ny;
    int ncells;
   
    //Current model parameters
    double Ra;  //Rayleigh number
    double dt;  //timestep
    double theta;  //Direction of gravity

    //Info about how to add heat on clicking
    double heat_source_radius;
    double heat_source;

    //The full grid, which knows how to iterate,
    //get geometric information, etc.
    StaggeredGrid grid;
    
    //Vectors which are used in the solve
    double *T;  //Temperature
    double *D;  //Displacement
    double *Dp;  //Previous displacement
    double *scratch1, *scratch2;  //Scratch vectors for storing temporary information
    double *freqs;  //Frequencies in spectral space for the stokes solve
    double *stream;  //Stream function solution
    double *curl_T; //Curl of temperature for the stream function solution
    double *u; //velocity in the x direction
    double *v; //Velocity in the y direction

    //Auxiliary vectors used in solving the diffusion equation
    double *g, *lux, *luy;
    double gamma; //Auxiliary constant for solving the diffusion equation. 

    //FFTW stuff
    fftw_plan dst, idst, dft, idft; //FFTW plans for doing the forward and inverse transforms
    std::complex<double>* curl_T_spectral;  //Curl of temperature in spectral space

    //Data for rendering with OpenGL
    GLfloat* vertices;  
    GLfloat* vertex_colors;
    GLuint* triangle_vertex_indices;

    //workhorse functions
    void initialize_temperature();  //just like it says
    double heat(const Point&, const Point&);  //heating term at a point, given where the click has happened
    void setup_stokes_problem();  //Setup for spectral solve
    void setup_diffusion_problem(); //Setup for diffusion solve (needs to be called every time the Ra is changed)
    void assemble_curl_T_vector(); //Assembling RHS for spectral solve.  Called every timestep
   
    //functions for evaluating field at points
    double initial_temperature(const Point&);
    double temperature(const Point&);
    Point velocity(const Point&);

  public:
    //Constructor and destructor
    ConvectionSimulator( double lx, double ly, int nx, int ny, double Rayleigh);
    ~ConvectionSimulator();

    //Querying physics information about the solver
    double rayleigh_number() const;  //return Ra
    double timescale() const;  //Characteristic timescale, which is scaled to plume ascent time
    double nusselt_number(); //Calculate nusselt number at a timestep

    void earthquake(double x, double y);  //Add source term for wave equation
    void propagate_seismic_waves(); //evolve the wave equation
    void clear_seismic_waves(); //zero out the displacement vectors

    void add_heat(double x, double y, bool hot); //Add heat at a point, with the bool indicating whether it is hot or cold
    void semi_lagrangian_advect();  //Advect temperature through the velocity field
    void diffuse_temperature(); //Diffuse temperature
    void solve_stokes(); //Solve for velocity field
    void draw();  //Render using OpenGL
    void update_state(double rayleigh, double gravity_angle = 0);  //Update the state of the solver

    void setup_opengl();  //Setup OpenGL data structures
    void cleanup_opengl(); //Cleanup OpenGL data structures
};

#endif
