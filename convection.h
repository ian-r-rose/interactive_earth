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

    enum advection_field { temperature, composition };

    //Current model parameters
    double Ra;  //Rayleigh number
    double buoyancy_number; // Density of chemical field relative to temperature variations
    double dt;  //timestep
    bool include_composition; //whether to include a compositional field

    //Info about how to add heati/composition on clicking
    double heat_source_radius;
    double heat_source;
    double composition_source_radius;
    double composition_source;

    //The full grid, which knows how to iterate,
    //get geometric information, etc.
    RegularGrid grid;

    //Vectors which are used in the solve
    double *T;  //Temperature
    double *C;  //Composition
    double *D;  //Displacement
    double *Dp;  //Previous displacement
    double *scratch;  //Scratch vectors for storing temporary information
    double *stream;  //Stream function solution
    double *curl_density; //Curl of temperature for the stream function solution
    double *u; //velocity in the x direction
    double *v; //Velocity in the y direction

    double spin;
    Point seismometer_location;

    TridiagonalMatrixSolver<std::complex<double> > **stokes_matrices;
    TridiagonalMatrixSolver<std::complex<double> > **diffusion_matrices;

    //FFTW stuff
    fftw_plan dft_stokes, idft_stokes; //FFTW plans for doing the forward and inverse transforms
    fftw_plan dft_diffusion, idft_diffusion; //FFTW plans for doing the forward and inverse transforms
    std::complex<double>* curl_density_spectral;  //Curl of density in spectral space
    std::complex<double>* T_spectral;  //Temperature in spectral space
    std::complex<double> *scratch1_spectral, *scratch2_spectral;  //Scratch complex vectors for temporary stuff


    //Data for rendering with OpenGL
    GLfloat* vertices;
    GLfloat* vertex_colors;
    GLuint* triangle_vertex_indices;

    //workhorse functions
    void initialize_temperature();  //just like it says
    void initialize_composition();  //just like it says
    double heat(const Point&, const Point&);  //heating term at a point, given where the click has happened
    double react(const Point&, const Point&);  //heating term at a point, given where the click has happened
    void setup_stokes_problem();  //Setup for spectral solve
    void setup_diffusion_problem(); //Setup for diffusion solve (needs to be called every time the Ra is changed)
    void assemble_curl_density_vector(); //Assembling RHS for spectral solve.  Called every timestep
    void semi_lagrangian_advect( const advection_field field );  //Advect composition or temperature through the velocity field
    void clip_field( advection_field field, double min, double max); //Clamp the field to be between min and max

    //functions for evaluating field at points
    double initial_temperature(const Point&);
    double initial_composition(const Point&);
    double evaluate_temperature(const Point&);
    double evaluate_composition(const Point&);
    Point evaluate_velocity(const Point&);

  public:
    //Constructor and destructor
    ConvectionSimulator( double inner_radius, int nx, int ny, bool include_composition);
    ~ConvectionSimulator();

    //Querying physics information about the solver
    double rayleigh_number() const;  //return Ra
    double timescale() const;  //Characteristic timescale, which is scaled to plume ascent time
    double nusselt_number(); //Calculate nusselt number at a timestep
    double spin_angle() const; //return angle of spin axis (with respect to x axis) in radians
    double seismometer_reading() const; //Return the displacement at the seismometer
    double timestep() const; //Return the timestep for the simulation.
    void seismometer_position( double &theta, double &r) const; //Query the current location of the seismometer

    void earthquake(double x, double y);  //Add source term for wave equation
    void propagate_seismic_waves(); //evolve the wave equation
    void clear_seismic_waves(); //zero out the displacement vectors
    void place_seismometer( double x, double y); //Add seismometer for monitoring waves at a point

    void add_heat(double x, double y, bool hot); //Add heat at a point, with the bool indicating whether it is hot or cold
    void add_composition(double x, double y); //Add composition at a point
    void semi_lagrangian_advect_temperature();  //Advect temperature through the velocity field
    void semi_lagrangian_advect_composition();  //Advect composition through the velocity field
    void diffuse_temperature(); //Diffuse temperature
    void true_polar_wander(); //Perform TPW

    void solve_stokes(); //Solve for velocity field

    void draw(bool);  //Render using OpenGL

    void setup_opengl();  //Setup OpenGL data structures
    void cleanup_opengl(); //Cleanup OpenGL data structures

    void update_state(double rayleigh);  //Update the state of the solver
};

#endif
