#ifndef STOKES_H
#define STOKES_H

#include <complex>
#include <vector>
#include "fftw3.h"
#include "regular_grid.h"
#include "tridiagonal_matrix_solver.h"

class ConvectionSimulator
{
private:
  enum class advection_field
  {
    temperature,
    composition
  };

  //Current model parameters
  double Ra;                //Rayleigh number
  double buoyancy_number;   // Density of chemical field relative to temperature variations
  double dt;                //timestep
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
  double *T;       //Temperature
  double *C;       //Composition
  double *D;       //Displacement
  double *Dp;      //Previous displacement
  double *scratch; //Scratch vectors for storing temporary information
  std::vector<double> grid_cos_theta, grid_sin_theta; //Of the angle of each grid column
  double *stream;                                     //Stream function solution
  double *curl_density; //Curl of temperature for the stream function solution
  double *u;            //velocity in the x direction
  double *v;            //Velocity in the y direction

  double spin;
  Point seismometer_location;

  //One system per Fourier mode, for the radial direction
  int n_modes() const { return grid.ntheta / 2 + 1; }
  BatchedTridiagonalSolver stokes_matrices;
  BatchedTridiagonalSolver diffusion_matrices;
  //The solves split the modes into blocks of this many, between threads
  static const int mode_block = 32;
  int n_mode_blocks() const { return (n_modes() + mode_block - 1) / mode_block; }

  //FFTW stuff
  fftw_plan row_dft, row_idft; //FFTW plans for the forward and inverse transforms of one row
  std::complex<double> *curl_density_spectral; //Curl of density in spectral space
  std::complex<double> *T_spectral;            //Temperature in spectral space

  //workhorse functions
  void initialize_temperature(); //just like it says
  void initialize_composition(); //just like it says
  double heat(const Point &,
              const Point &); //heating term at a point, given where the click has happened
  double react(const Point &,
               const Point &); //heating term at a point, given where the click has happened
  void setup_fft_plans();      //Plan the transforms, once
  void fft_rows(const double *in, std::complex<double> *out); //Transform each row
  void ifft_rows(std::complex<double> *in, double *out);      //Inverse transform each row
  void setup_stokes_problem();                                //Setup for spectral solve
  void
  setup_diffusion_problem(); //Setup for diffusion solve (needs to be called every time the Ra is changed)
  void assemble_curl_density_vector(); //Assembling RHS for spectral solve.  Called every timestep
  void semi_lagrangian_advect(
      const advection_field field); //Advect composition or temperature through the velocity field
  void clip_field(advection_field field, double min,
                  double max); //Clamp the field to be between min and max

  //functions for evaluating field at points
  double initial_temperature(const Point &);
  double initial_composition(const Point &);

public:
  //Constructor and destructor
  ConvectionSimulator(double inner_radius, int nx, int ny, bool include_composition);
  ~ConvectionSimulator();

  //Querying physics information about the solver
  double rayleigh_number() const; //return Ra
  double timescale() const;       //Characteristic timescale, which is scaled to plume ascent time
  double nusselt_number();        //Calculate nusselt number at a timestep
  double spin_angle() const;      //return angle of spin axis (with respect to x axis) in radians
  double seismometer_reading() const; //Return the displacement at the seismometer
  double timestep() const;            //Return the timestep for the simulation.
  void seismometer_position(double &theta,
                            double &r) const; //Query the current location of the seismometer

  //Read-only access to the solver state, for rendering and diagnostics
  const RegularGrid &get_grid() const { return grid; }
  bool has_composition() const { return include_composition; }
  const double *temperature() const { return T; }
  const double *composition() const { return C; }
  const double *displacement() const { return D; }
  const double *velocity_theta() const { return u; }
  const double *velocity_r() const { return v; }
  const double *stream_function() const { return stream; }

  void set_temperature(const double *values); //Overwrite the temperature field (ncells values)

  void earthquake(double x, double y);        //Add source term for wave equation
  void propagate_seismic_waves();             //evolve the wave equation
  void clear_seismic_waves();                 //zero out the displacement vectors
  void place_seismometer(double x, double y); //Add seismometer for monitoring waves at a point

  void add_heat(double x, double y,
                bool hot); //Add heat at a point, with the bool indicating whether it is hot or cold
  void add_composition(double x, double y);  //Add composition at a point
  void semi_lagrangian_advect_temperature(); //Advect temperature through the velocity field
  void semi_lagrangian_advect_composition(); //Advect composition through the velocity field
  void diffuse_temperature();                //Diffuse temperature
  void true_polar_wander();                  //Perform TPW

  void solve_stokes(); //Solve for velocity field

  void update_state(double rayleigh); //Update the state of the solver
};

#endif
