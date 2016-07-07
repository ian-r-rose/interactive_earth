#include <cmath>
#include "convection.h"
#include <iostream>

//standard math library functions can be unpredictibly slow in
//some implementations, it seems, and not even available in others.
//Here I implement some super basic, unsafe versions of some
//for inlining.
inline double fast_fmod(double x,double y) { return x-((int)(x/y))*y; }
inline double dmin (double x, double y) { return x < y ? x : y; }
inline double dmax (double x, double y) { return x > y ? x : y; }


/*Constructor for the solver. The parameters are in order:
  inner_radius : inner radius, compared to an outer radius of one.
  ntheta : number of cells in theta direction
  nr : number of cells in r direction
  Rayleigh : initial Rayleigh number
*/
ConvectionSimulator::ConvectionSimulator( double inner_radius, int ntheta, int nr):
                          grid(inner_radius, ntheta, nr)
{
  //Allocate memory for data vectors
  T = new double[grid.ncells];
  V = new double[grid.ncells];
  stream = new double[grid.ncells];
  vorticity_source = new double[grid.ncells];
  u = new double[grid.ncells];
  v = new double[grid.ncells];
  scratch = new double[grid.ncells];

  //I actually allocate more memory than necessary for the
  //spectral vector, but this way it makes the indexing simpler
  vorticity_source_spectral = new std::complex<double>[grid.ncells];
  T_spectral = new std::complex<double>[grid.ncells];
  V_spectral = new std::complex<double>[grid.ncells];
  scratch1_spectral = new std::complex<double>[grid.ncells];

  poisson_matrices = new TridiagonalMatrixSolver<std::complex<double> >*[grid.ntheta/2+1];
  temperature_diffusion_matrices = new TridiagonalMatrixSolver<std::complex<double> >*[grid.ntheta/2+1];
  vorticity_diffusion_matrices = new TridiagonalMatrixSolver<std::complex<double> >*[grid.ntheta/2+1];
  for (unsigned int i=0; i<=grid.ntheta/2; ++i)
  {
    poisson_matrices[i] = new TridiagonalMatrixSolver<std::complex<double> >(grid.nr);
    temperature_diffusion_matrices[i] = new TridiagonalMatrixSolver<std::complex<double> >(grid.nr);
    vorticity_diffusion_matrices[i] = new TridiagonalMatrixSolver<std::complex<double> >(grid.nr);
  }

  //Initialize the state
  update_state(1.e6, 1.0e-8, 100.1);
  initialize_temperature();
  initialize_vorticity();

  //Do some setup work for solving stokes and
  //diffustion problems.
  setup_poisson_problem();
  setup_vorticity_diffusion_problem();
  setup_temperature_diffusion_problem();
}

/* Destructor for the solver.*/
ConvectionSimulator::~ConvectionSimulator()
{
  delete[] T;
  delete[] V;
  delete[] stream;
  delete[] vorticity_source;
  delete[] u;
  delete[] v;
  delete[] scratch;
  delete[] vorticity_source_spectral;
  delete[] T_spectral;
  delete[] V_spectral;
  delete[] scratch1_spectral;

  for (unsigned int i=0; i<=grid.ntheta/2; ++i)
  {
    delete poisson_matrices[i];
    delete temperature_diffusion_matrices[i];
    delete vorticity_diffusion_matrices[i];
  }
  delete[] poisson_matrices;
  delete[] vorticity_diffusion_matrices;
  delete[] temperature_diffusion_matrices;

  fftw_destroy_plan(dft_temperature_diffusion);
  fftw_destroy_plan(idft_temperature_diffusion);
  fftw_destroy_plan(dft_vorticity_diffusion);
  fftw_destroy_plan(idft_vorticity_diffusion);
  fftw_destroy_plan(dft_poisson);
  fftw_destroy_plan(idft_poisson);
  fftw_cleanup();
}

/* Functional form of initial temperature field.
   Just start with a constant value based on a guess of the steady state average.*/
double ConvectionSimulator::initial_temperature(const Point &p)
{
  //Guess of internal temperature based on isoviscous boundary layer theory
  //in a cylinder.
  return grid.r_inner/(1.+grid.r_inner);
}

/* Functional form of initial vorticity.
   Just start with zero.*/
double ConvectionSimulator::initial_vorticity(const Point &p)
{
  return 0.0;
}

/* Loop over all the cells and set the initial temperature */
void ConvectionSimulator::initialize_temperature()
{
  for( RegularGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
    T[cell->self()] = initial_temperature(cell->location());
}

/* Loop over all the cells and set the initial vorticity */
void ConvectionSimulator::initialize_vorticity()
{
  for( RegularGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
    V[cell->self()] = initial_vorticity(cell->location());
}

/* Given a heat source centered on p1, calculate the heating at
   point p2, using a simple Gaussian source term */
inline double ConvectionSimulator::heat(const Point &p1, const Point &p2 )
{
  const double x1 = (p1.r+grid.r_inner)*std::cos(p1.theta), y1 = (p1.r+grid.r_inner)*std::sin(p1.theta);
  const double x2 = (p2.r+grid.r_inner)*std::cos(p2.theta), y2 = (p2.r+grid.r_inner)*std::sin(p2.theta);
  const double rsq = (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);
  return heat_source*std::exp( -rsq/2.0/heat_source_radius/heat_source_radius );
}

/* Loop over all the cells and add heat according to where the current heat
   source is. */
void ConvectionSimulator::add_heat(double theta, double r, bool hot)
{
  Point p; p.theta = theta; p.r=r;
  #pragma omp parallel for
  for( unsigned int cell_index=0; cell_index<grid.ncells; ++cell_index)
  {
    RegularGrid::iterator cell(cell_index, grid);
    double temperature = T[cell->self()] + (hot ? 1.0 : -1.0)*heat(p, cell->location())*dt;
    temperature = dmax(dmin(temperature,1.0), 0.0);
    T[cell->self()] = temperature;
  }
}

/* Loop over all the cells and add vorticity, similar to add_heat*/
void ConvectionSimulator::add_vorticity(double theta, double r, bool ccw)
{
  Point p; p.theta = theta; p.r=r;
  #pragma omp parallel for
  for( unsigned int cell_index=0; cell_index<grid.ncells; ++cell_index)
  {
    RegularGrid::iterator cell(cell_index, grid);
    double vorticity = V[cell->self()] + heat(p, cell->location())*dt;
    V[cell->self()] = vorticity;
  }
}

void ConvectionSimulator::generate_vorticity()
{
    #pragma omp parallel for
    for( unsigned int cell_index=0; cell_index<grid.ncells; ++cell_index)
    {
      RegularGrid::iterator cell(cell_index, grid);
      const double r = cell->radius();
      const double curl_T = -(T[cell->right()]-T[cell->left()])/2./grid.dtheta/r;
      V[cell->self()] += Pr*Ra*curl_T*dt;
    }
}

/* Interpolate the velocity onto an arbitrary point
   A bit complicated*/
Point ConvectionSimulator::evaluate_velocity(const Point &p)
{
  Point vel;

  //Get the relevant lower-left cells for the velocities
  RegularGrid::iterator cell = grid.lower_left_cell(p);

  //Determine the local x and y coordinates for the velocity
  double local_theta = fast_fmod(p.theta, grid.dtheta)/grid.dtheta;
  double local_r = (p.r - cell->location().r)/grid.dr;

  unsigned int ul, ur, dl, dr;
  if (cell->at_top_boundary() )
  {
    ul = cell->self();
    ur = cell->right();
  }
  else
  {
    ul = cell->up();
    ur = cell->right();
  }
  dl = cell->self();
  dr = cell->right();

  //get interpolated vx
  vel.theta = linear_interp_2d (local_theta, local_r, u[ul], u[ur],
                            u[dl], u[dr]);
  vel.r = linear_interp_2d (local_theta, local_r, v[ul], v[ur],
                            v[dl], v[dr]);

  return vel;
}

/* Interpolate the temperature onto an arbitrary point. */
inline double ConvectionSimulator::evaluate_temperature(const Point &p)
{
  double temp;
  RegularGrid::iterator cell = grid.lower_left_cell(p);
  double local_theta = fast_fmod(p.theta, grid.dtheta)/grid.dtheta;
  double local_r = ( p.r - cell->location().r )/grid.dr;

  if (cell->at_top_boundary() )
    temp = linear_interp_2d( local_theta, local_r, -T[cell->self()], -T[cell->right()],
                             T[cell->self()], T[cell->right()]);
  else
    temp = linear_interp_2d( local_theta, local_r, T[cell->up()], T[cell->upright()],
                             T[cell->self()], T[cell->right()]);

  return temp;
}

/* Interpolate the vorticity onto an arbitrary point. */
inline double ConvectionSimulator::evaluate_vorticity(const Point &p)
{
  if (p.r < 0.0 || p.r > grid.lr ) return 0.0;

  double vort;
  RegularGrid::iterator cell = grid.lower_left_cell(p);
  double local_theta = fast_fmod(p.theta, grid.dtheta)/grid.dtheta;
  double local_r = ( p.r - cell->location().r )/grid.dr;

  if (cell->at_top_boundary() )
    vort = linear_interp_2d( local_theta, local_r, V[cell->self()], V[cell->right()],
                             V[cell->self()], V[cell->right()]);
  else
    vort = linear_interp_2d( local_theta, local_r, V[cell->up()], V[cell->upright()],
                             V[cell->self()], V[cell->right()]);

  return vort;
}

void ConvectionSimulator::semi_lagrangian_advect_temperature()
{
  advection_field field = temperature;
  semi_lagrangian_advect( field );
}

void ConvectionSimulator::semi_lagrangian_advect_vorticity()
{
  advection_field field = vorticity;
  semi_lagrangian_advect( field );
}

/* Advect the temperature field through the velocity field using
   Semi-lagrangian advection.  This scheme is quite stable, which
   allows me to take VERY large time steps.  The drawback is that it
   is kind of slow.  Here I do it in a very coarse way to make up for
   that.  A more accurate implementation would use better-than-linear
   interpolation and use more iterations. */
void ConvectionSimulator::semi_lagrangian_advect( const advection_field field )
{
  //The goal is to find the temperature at the Lagrangian point which will
  //be advected to the current grid point in one time step.  In general,
  //this point wil not be on the grid.  I find this point using a coarse
  //iterated predictor corrector.

  #pragma omp parallel for
  for( unsigned int cell_index=0; cell_index<grid.ncells; ++cell_index)
  {
    RegularGrid::iterator cell(cell_index, grid);

    Point vel_final, vel_final_cart;  //Velocity at the grid point
    Point vel_takeoff, vel_takeoff_cart; //Velocity at the candidate takeoff point
    Point takeoff_point, takeoff_point_cart; //Candidate takeoff point
    Point final_point, final_point_cart; //grid point

    //These points are known, as they are the grid points in question.  They will
    //not change for this cell.
    vel_final.theta = u[cell_index];
    vel_final.r = v[cell_index];
    final_point = cell->location();
    const double final_r = final_point.r + grid.r_inner;
    const double cos_f = std::cos(final_point.theta);
    const double sin_f = std::sin(final_point.theta);
    final_point_cart.x = final_r * cos_f;
    final_point_cart.y = final_r * sin_f;
    vel_final_cart.x = -vel_final.theta * sin_f + vel_final.r * cos_f;
    vel_final_cart.y =  vel_final.theta * cos_f + vel_final.r * sin_f;


    //Calculate the initial guess for the takeoff point using a basic predictor
    //with a forward euler step
    takeoff_point_cart.x = final_point_cart.x - vel_final_cart.x*dt;
    takeoff_point_cart.y = final_point_cart.y - vel_final_cart.y*dt;
    takeoff_point.theta = std::atan2(takeoff_point_cart.y, takeoff_point_cart.x);
    takeoff_point.r = std::sqrt( takeoff_point_cart.x*takeoff_point_cart.x + takeoff_point_cart.y*takeoff_point_cart.y) - grid.r_inner;
    //Keep it in the domain
    takeoff_point.theta = fast_fmod(takeoff_point.theta + grid.lx, grid.lx);
    takeoff_point.r = dmin( grid.lr, dmax( takeoff_point.r, 0.0) );

    //Iterate on the corrector.  Here I only do one iteration for
    //performance reasons, but in principle we could do more to get
    //a better estimate.
    for(unsigned int i=0; i<1; ++i)
    {
      //Evaluate the velocity at the predictor
      vel_takeoff = evaluate_velocity(takeoff_point);
      const double takeoff_r = takeoff_point.r + grid.r_inner;
      const double cos_t = std::cos(takeoff_point.theta);
      const double sin_t = std::sin(takeoff_point.theta);
      takeoff_point_cart.x = takeoff_r * cos_t;
      takeoff_point_cart.y = takeoff_r * sin_t;
      vel_takeoff_cart.x = -vel_takeoff.theta * sin_t + vel_takeoff.r * cos_t;
      vel_takeoff_cart.y =  vel_takeoff.theta * cos_t + vel_takeoff.r * sin_t;

      //Come up with the corrector using a midpoint rule
      takeoff_point_cart.x = final_point_cart.x - (vel_final_cart.x + vel_takeoff_cart.x)*dt/2.0;
      takeoff_point_cart.y = final_point_cart.y - (vel_final_cart.y + vel_takeoff_cart.y)*dt/2.0;
      takeoff_point.theta = std::atan2(takeoff_point_cart.y, takeoff_point_cart.x);
      takeoff_point.r = std::sqrt( takeoff_point_cart.x*takeoff_point_cart.x + takeoff_point_cart.y*takeoff_point_cart.y) - grid.r_inner;
      //Keep in domain
      takeoff_point.theta = fast_fmod(takeoff_point.theta + grid.lx, grid.lx);
      takeoff_point.r = dmin( grid.lr, dmax( takeoff_point.r, 0.0) );
    }
    if (field == temperature)
      scratch[cell->self()] = evaluate_temperature(takeoff_point);  //Store the temperature we found
    else if (field == vorticity)
      scratch[cell->self()] = evaluate_vorticity(takeoff_point);  //Store the vorticity we found
  }

  //Copy the scratch vector into the temperature or vorticity vector.
  if (field == temperature)
    for( RegularGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
      T[cell->self()] = scratch[cell->self()];
  else if (field == vorticity)
    for( RegularGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
      V[cell->self()] = scratch[cell->self()];
}

/* Solve the temperature diffusion equation implicitly with reverse Euler.
   I do the x and y direction separately, with the x direction
   done spectrally using FFTs, and the y direction with finite
   differences, solving the resulting tridiagonal system using
   the O(N) Thomas algorithm. */
void ConvectionSimulator::diffuse_temperature()
{
  //Apply boundary conditions
  for( RegularGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
  {
    if (cell->at_top_boundary())
      T[cell->self()] = 0.0;
    else if (cell->at_bottom_boundary())
      T[cell->self()] = 1.0;
  }

  //Execute the forward fourier transform
  fftw_execute(dft_temperature_diffusion);  //X direction

  #pragma omp parallel for
  for (unsigned int l = 0; l <= grid.ntheta/2; ++l)
    temperature_diffusion_matrices[l]->solve( T_spectral+l, grid.ntheta, scratch1_spectral+l);

  //Execute the inverse Fourier transform
  fftw_execute(idft_temperature_diffusion);  //X direction

  //Renormalize
  for( RegularGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
    T[cell->self()] = scratch[cell->self()]/grid.ntheta;
}

/* Solve the vorticity diffusion equation implicitly with reverse Euler.
   I do the x and y direction separately, with the x direction
   done spectrally using FFTs, and the y direction with finite
   differences, solving the resulting tridiagonal system using
   the O(N) Thomas algorithm. */
void ConvectionSimulator::diffuse_vorticity()
{
  //Apply boundary conditions
  for( RegularGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
  {
    if (cell->at_top_boundary())
      V[cell->self()] = 0.0;
    else if (cell->at_bottom_boundary())
      V[cell->self()] = 0.0;
  }

  //Execute the forward fourier transform
  fftw_execute(dft_vorticity_diffusion);  //X direction

  #pragma omp parallel for
  for (unsigned int l = 0; l <= grid.ntheta/2; ++l)
    vorticity_diffusion_matrices[l]->solve( V_spectral+l, grid.ntheta, scratch1_spectral+l);

  //Execute the inverse Fourier transform
  fftw_execute(idft_vorticity_diffusion);  //X direction

  //Renormalize
  for( RegularGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
    V[cell->self()] = scratch[cell->self()]/grid.ntheta;
}

void ConvectionSimulator::setup_temperature_diffusion_problem()
{
  //Setup the tridiagonal matrices in the y direction
  double *upper_diag = new double[grid.nr];
  double *lower_diag = new double[grid.nr];

  const double dr = grid.dr;
  //Upper and lower diagonals are the same for each frequency
  for( unsigned int i=0; i < grid.nr; ++i )
  {
    const double r = grid.r_inner + i*dr;
    upper_diag[i] = -dt/dr/dr * (1. + 0.5*dr/r);
    lower_diag[i] = -dt/dr/dr * (1. - 0.5*dr/r);
  }
  upper_diag[0] = 0.0;  //fix for lower B.C.
  lower_diag[grid.nr-1] = 0.0; //fix for upper B.C.

  //Diagonals include a term for the frequency in the x-direction
  #pragma omp parallel for
  for (unsigned int l = 0; l <= grid.ntheta/2; ++l)
  {
     double factor = (l*l);
     double *diag = new double[grid.nr];
     diag[0] = 1.0; diag[grid.nr-1] = 1.0;
     for (unsigned int i=1; i<grid.nr-1; ++i)
     {
       const double r = grid.r_inner + i*dr;
       diag[i] = 1. + 2.*dt/dr/dr + factor*dt/r/r;
     }

     temperature_diffusion_matrices[l]->initialize(lower_diag, diag, upper_diag);
     delete[] diag;
  }

  delete[] upper_diag;
  delete[] lower_diag;

  int n[1];
  int stride, dist, howmany;

  //transform in the x direction with a full dft
  n[0] = grid.ntheta; stride = 1; dist = grid.ntheta; howmany = grid.nr;
  dft_temperature_diffusion = fftw_plan_many_dft_r2c(1, n, howmany, T, NULL, stride, dist,
                                         reinterpret_cast<fftw_complex*>(T_spectral),
                                         NULL, stride, dist, FFTW_ESTIMATE);
  idft_temperature_diffusion = fftw_plan_many_dft_c2r(1, n, howmany,
                                          reinterpret_cast<fftw_complex*>(scratch1_spectral),
                                          NULL, stride, dist, scratch, NULL,
                                          stride, dist, FFTW_ESTIMATE);

}

void ConvectionSimulator::setup_vorticity_diffusion_problem()
{
  //Setup the tridiagonal matrices in the y direction
  double *upper_diag = new double[grid.nr];
  double *lower_diag = new double[grid.nr];

  const double dr = grid.dr;
  //Upper and lower diagonals are the same for each frequency
  for( unsigned int i=0; i < grid.nr; ++i )
  {
    const double r = grid.r_inner + i*dr;
    upper_diag[i] = -dt/dr/dr * (1. + 0.5*dr/r) * Pr;
    lower_diag[i] = -dt/dr/dr * (1. - 0.5*dr/r) * Pr;
  }
  upper_diag[0] = 0.0;  //fix for lower B.C.
  lower_diag[grid.nr-1] = 0.0; //fix for upper B.C.

  //Diagonals include a term for the frequency in the x-direction
  #pragma omp parallel for
  for (unsigned int l = 0; l <= grid.ntheta/2; ++l)
  {
     double factor = (l*l);
     double *diag = new double[grid.nr];
     diag[0] = 1.0; diag[grid.nr-1] = 1.0;
     for (unsigned int i=1; i<grid.nr-1; ++i)
     {
       const double r = grid.r_inner + i*dr;
       diag[i] = (1. + 2.*dt/dr/dr*Pr + factor*dt/r/r*Pr);
     }

     vorticity_diffusion_matrices[l]->initialize(lower_diag, diag, upper_diag);
     delete[] diag;
  }

  delete[] upper_diag;
  delete[] lower_diag;

  int n[1];
  int stride, dist, howmany;

  //transform in the x direction with a full dft
  n[0] = grid.ntheta; stride = 1; dist = grid.ntheta; howmany = grid.nr;
  dft_vorticity_diffusion = fftw_plan_many_dft_r2c(1, n, howmany, V, NULL, stride, dist,
                                         reinterpret_cast<fftw_complex*>(V_spectral),
                                         NULL, stride, dist, FFTW_ESTIMATE);
  idft_vorticity_diffusion = fftw_plan_many_dft_c2r(1, n, howmany,
                                          reinterpret_cast<fftw_complex*>(scratch1_spectral),
                                          NULL, stride, dist, scratch, NULL,
                                          stride, dist, FFTW_ESTIMATE);

}

/* Setup the stokes solve with the stream function formulation.
   We assemble the vector that stores the frequencies of the
   eigenmodes, as well as tell FFTW how to do the transforms */
void ConvectionSimulator::setup_poisson_problem()
{
  double *upper_diag = new double[grid.nr];
  double *diag = new double[grid.nr];
  double *lower_diag = new double[grid.nr];

  const double dr = grid.dr;
  for( unsigned int i=0; i < grid.nr; ++i )
  {
    const double r = grid.r_inner + i*dr;
    upper_diag[i] = 1./dr/dr + 0.5/dr/r;
    lower_diag[i] = 1./dr/dr - 0.5/dr/r;
  }
  upper_diag[0] = 0.0;  //fix for lower B.C.
  lower_diag[grid.nr-1] = 0.0; //fix for upper B.C.

  for (unsigned int l = 0; l <= grid.ntheta/2; ++l)
  {
     double factor = (l*l);
     diag[0] = 1.0; diag[grid.nr-1] = 1.0;
     for (unsigned int i=1; i<grid.nr-1; ++i)
     {
       const double r = grid.r_inner + i*dr;
       diag[i] = -2./dr/dr - factor/r/r;
     }

     poisson_matrices[l]->initialize(lower_diag, diag, upper_diag);
  }

  delete[] upper_diag;
  delete[] diag;
  delete[] lower_diag;

  int n[1];
  int stride, dist, howmany;

  //transform in the x direction with a full dft
  n[0] = grid.ntheta; stride = 1; dist = grid.ntheta; howmany = grid.ny;
  dft_poisson = fftw_plan_many_dft_r2c(1, n, howmany, V, NULL, stride, dist,
                                      reinterpret_cast<fftw_complex*>(V_spectral),
                                      NULL, stride, dist, FFTW_ESTIMATE);
  idft_poisson = fftw_plan_many_dft_c2r(1, n, howmany,
                                       reinterpret_cast<fftw_complex*>(scratch1_spectral),
                                       NULL, stride, dist, scratch, NULL,
                                       stride, dist, FFTW_ESTIMATE);

}


/*Calculate the Nusselt number for the given temperature field.
  In principle I could calculate the heat flux through the bottom
  boundary as well, and take the average, but as it is, I am only
  doing the heat flux through the top boundary*/
double ConvectionSimulator::nusselt_number()
{
  double heat_flux = 0;

  for ( RegularGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
  {
    if( cell->at_top_boundary() )
    {
      double grad_T = (T[cell->self()] - T[cell->down()])/grid.dr;
      heat_flux += grad_T;
    }
  }
  return -heat_flux * grid.lr/grid.ntheta;
}


void ConvectionSimulator::assemble_vorticity_source_vector()
{
  #pragma omp parallel for
  for( unsigned int cell_index=0; cell_index<grid.ncells; ++cell_index)
  {
    RegularGrid::iterator cell(cell_index, grid);
    const double r = cell->radius();
  }
}


/* Solve the Poisson equation for the stream function.*/
void ConvectionSimulator::solve_poisson()
{
  //Execute the forward fourier transform
  fftw_execute(dft_poisson);  //X direction

  #pragma omp parallel for
  for (unsigned int l = 0; l <= grid.ntheta/2; ++l)
     poisson_matrices[l]->solve( V_spectral+l, grid.ntheta, scratch1_spectral+l);

  //Execute the inverse Fourier transform
  fftw_execute(idft_poisson);  //X direction

  //Renormalize
  for( RegularGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
    stream[cell->self()] = scratch[cell->self()]/grid.ntheta;

  //Come up with the velocities by taking finite differences of the stream function.
  //I could also take these derivatives in spectral space, but that would mean
  //more fourier transforms, so this should be considerably cheaper.
  #pragma omp parallel for
  for( unsigned int cell_index=0; cell_index<grid.ncells; ++cell_index)
  {
    RegularGrid::iterator cell(cell_index, grid);
    const double r = cell->radius();
    if( cell->at_top_boundary() )
      u[cell->self()] = (stream[cell->self()] - stream[cell->down()])/grid.dr;
    else if ( cell->at_bottom_boundary() )
      u[cell->self()] = (stream[cell->up()] - stream[cell->self()])/grid.dr;
    else
      u[cell->self()] = (stream[cell->up()] - stream[cell->down()])/2.0/grid.dr;

    v[cell->self()] = -(stream[cell->right()] - stream[cell->left()])/2.0/grid.dtheta/r;
  }
}

void ConvectionSimulator::update_state(double rayleigh, double ekman, double prandtl)
{
  const double Ra_c = 657.;

  Ra = rayleigh;
  Ek = ekman;
  Pr = prandtl;

  const double length_scale = std::pow(Ra/2./Ra_c, -1./3.)*grid.lr;
  const double Nu = std::pow(Ra/Ra_c/2., 1./3.) / 2.;  //Nusselt
  const double velocity_scale = std::sqrt( Ra * Nu );
  const double cfl = dmin(grid.dr, grid.r_inner*grid.dtheta)/velocity_scale;

  //Estimate other state properties based on simple isoviscous scalings
  dt = cfl * 20.0; //Roughly 10x CFL, thanks to semi-lagrangian
  heat_source_radius = length_scale*0.5;  //Radius of order the boundary layer thickness
  heat_source = 1.*velocity_scale/grid.lr*2.; //Heat a blob of order the ascent time for thta blob
  vorticity_source_radius = grid.lr/10.; //Size of vorticity blob
  vorticity_strength=heat_source; //Strength of vorticity addition

  setup_temperature_diffusion_problem();  //need to recompute the auxiliary vectors for the diffusion problem
  setup_vorticity_diffusion_problem();  //need to recompute the auxiliary vectors for the diffusion problem

}

double ConvectionSimulator::rayleigh_number() const
{
  return Ra;
}

double ConvectionSimulator::ekman_number() const
{
  return Ek;
}

double ConvectionSimulator::prandtl_number() const
{
  return Pr;
}
