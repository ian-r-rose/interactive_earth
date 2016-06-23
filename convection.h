#ifndef STOKES_H
#define STOKES_H

#include <complex>
#include "fftw3.h"
#include "regular_grid.h"
#include "GL/glew.h"
#include "SDL2/SDL.h"
#include "SDL2/SDL_opengl.h"
#include "tridiagonal_matrix_solver.h"


class ConvectionSimulator
{
  private:

    enum advection_field { temperature, vorticity };

    //Current model parameters
    double Ra;  //Rayleigh number
    double Ek; //Ekman number
    double Pr; //Prandtl number 
    
    double dt;  //timestep

    //Info about how to add heat/vorticity on clicking
    double heat_source_radius;
    double heat_source;
    double vorticity_source_radius;
    double vorticity_strength;

    //The full grid, which knows how to iterate,
    //get geometric information, etc.
    RegularGrid grid;

    //Vectors which are used in the solve
    double *T;  //Temperature
    double *V;  //Vorticity
    double *scratch;  //Scratch vectors for storing temporary information
    double *stream;  //Stream function solution
    double *vorticity_source; //Curl of temperature for the stream function solution
    double *u; //velocity in the x direction
    double *v; //Velocity in the y direction

    TridiagonalMatrixSolver<std::complex<double> > **poisson_matrices;
    TridiagonalMatrixSolver<std::complex<double> > **temperature_diffusion_matrices;
    TridiagonalMatrixSolver<std::complex<double> > **vorticity_diffusion_matrices;

    //FFTW stuff
    fftw_plan dft_poisson, idft_poisson; //FFTW plans for doing the forward and inverse transforms
    fftw_plan dft_temperature_diffusion, idft_temperature_diffusion; //FFTW plans for doing the forward and inverse transforms
    fftw_plan dft_vorticity_diffusion, idft_vorticity_diffusion; //FFTW plans for doing the forward and inverse transforms
    std::complex<double>* vorticity_source_spectral;  //vorticity source in spectral space
    std::complex<double>* T_spectral;  //Temperature in spectral space
    std::complex<double>* V_spectral;  //Vorticity in spectral space
    std::complex<double> *scratch1_spectral, *scratch2_spectral;  //Scratch complex vectors for temporary stuff


    //Data for rendering with OpenGL
    GLfloat* vertices;
    GLfloat* vertex_colors;
    GLuint* triangle_vertex_indices;

    //workhorse functions
    void initialize_temperature();  //just like it says
    void initialize_vorticity();  //just like it says
    void setup_poisson_problem();
    void setup_vorticity_diffusion_problem();  //Setup for vorticity diffusion solve
    void setup_temperature_diffusion_problem(); //Setup for diffusion solve
    void assemble_vorticity_source_vector(); //Assembling RHS for spectral solve.  Called every timestep
    void semi_lagrangian_advect( const advection_field field );  //Advect vorticity or temperature through the velocity field
    void clip_field( advection_field field, double min, double max); //Clamp the field to be between min and max

    //functions for evaluating field at points
    double initial_temperature(const Point&);
    double initial_vorticity(const Point&);
    double evaluate_temperature(const Point&);
    double evaluate_vorticity(const Point&);
    Point evaluate_velocity(const Point&);

    double heat(const Point &p1, const Point &p2 );
  public:
    //Constructor and destructor
    ConvectionSimulator( double inner_radius, int nx, int ny);
    ~ConvectionSimulator();

    //Querying physics information about the solver
    double rayleigh_number() const;  //return Ra
    double ekman_number() const;  //return Ra
    double prandtl_number() const;  //return Ra
    double nusselt_number(); //Calculate nusselt number at a timestep

    void add_heat(double x, double y, bool hot); //Add heat at a point, with the bool indicating whether it is hot or cold
    void add_vorticity(double x, double y, bool ccw); //Add vorticity at a point, with a bool indicating counterclockwise
    void semi_lagrangian_advect_temperature();  //Advect temperature through the velocity field
    void semi_lagrangian_advect_vorticity();  //Advect vorticity through the velocity field
    void diffuse_temperature(); //Diffuse temperature
    void diffuse_vorticity(); //Diffuse vorticity

    void solve_poisson(); //Solve for velocity field

    void draw(bool);  //Render using OpenGL

    void setup_opengl();  //Setup OpenGL data structures
    void cleanup_opengl(); //Cleanup OpenGL data structures

    void update_state(double rayleigh, double ekman, double prandtl);  //Update the state of the solver
};

#endif
