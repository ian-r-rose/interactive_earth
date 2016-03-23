#include <cmath>
#include "convection.h"

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
ConvectionSimulator::ConvectionSimulator( double inner_radius, int ntheta, int nr, double Rayleigh, bool include_composition):
                          Ra(Rayleigh),
                          grid(inner_radius, ntheta, nr),
                          include_composition (include_composition)
{
  //Allocate memory for data vectors
  T = new double[grid.ncells];
  D = new double[grid.ncells];
  Dp = new double[grid.ncells];
  stream = new double[grid.ncells];
  curl_density = new double[grid.ncells];
  u = new double[grid.ncells];
  v = new double[grid.ncells];
  scratch = new double[grid.ncells];

  if(include_composition)
    C = new double[grid.ncells];

  //I actually allocate more memory than necessary for the
  //spectral vector, but this way it makes the indexing simpler
  curl_density_spectral = new std::complex<double>[grid.ncells];
  T_spectral = new std::complex<double>[grid.ncells];
  scratch1_spectral = new std::complex<double>[grid.ncells];
  scratch2_spectral = new std::complex<double>[grid.ncells];

  stokes_matrices = new TridiagonalMatrixSolver<std::complex<double> >*[grid.ntheta/2+1];
  diffusion_matrices = new TridiagonalMatrixSolver<std::complex<double> >*[grid.ntheta/2+1];
  for (unsigned int i=0; i<=grid.ntheta/2; ++i)
  {
    stokes_matrices[i] = new TridiagonalMatrixSolver<std::complex<double> >(grid.nr);
    diffusion_matrices[i] = new TridiagonalMatrixSolver<std::complex<double> >(grid.nr);
  }

  //Initialize the state
  buoyancy_number = 1.0;
  update_state(Rayleigh);
  initialize_temperature();
  if(include_composition)
    initialize_composition();

  //Initialize displacement to zero
  clear_seismic_waves();

  //Do some setup work for solving stokes and
  //diffustion problems.
  setup_stokes_problem();
  setup_diffusion_problem();
}

/* Destructor for the solver.*/
ConvectionSimulator::~ConvectionSimulator()
{
  delete[] T;
  delete[] D;
  delete[] Dp;
  delete[] stream;
  delete[] curl_density;
  delete[] u;
  delete[] v;
  delete[] scratch;
  delete[] curl_density_spectral;
  delete[] T_spectral;
  delete[] scratch1_spectral;
  delete[] scratch2_spectral;

  if(include_composition)
    delete[] C;

  for (unsigned int i=0; i<=grid.ntheta/2; ++i)
  {
    delete stokes_matrices[i];
    delete diffusion_matrices[i];
  }
  delete[] stokes_matrices;
  delete[] diffusion_matrices;

  fftw_destroy_plan(dft_diffusion);
  fftw_destroy_plan(idft_diffusion);
  fftw_destroy_plan(dft_stokes);
  fftw_destroy_plan(idft_stokes);
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

/* Functional form of initial composition.
   Just start with zero.*/
double ConvectionSimulator::initial_composition(const Point &p)
{
  return 0.0;
}

/* Loop over all the cells and set the initial temperature */
void ConvectionSimulator::initialize_temperature()
{
  for( RegularGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
    T[cell->self()] = initial_temperature(cell->location());
}

/* Loop over all the cells and set the initial composition */
void ConvectionSimulator::initialize_composition()
{
  for( RegularGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
    C[cell->self()] = initial_composition(cell->location());
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

/* Given a heat source centered on p1, calculate the composition addition at point p2*/ 
inline double ConvectionSimulator::react(const Point &p1, const Point &p2 )
{
  const double x1 = (p1.r+grid.r_inner)*std::cos(p1.theta), y1 = (p1.r+grid.r_inner)*std::sin(p1.theta);
  const double x2 = (p2.r+grid.r_inner)*std::cos(p2.theta), y2 = (p2.r+grid.r_inner)*std::sin(p2.theta);
  const double rsq = (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);
  return rsq < composition_source_radius*composition_source_radius ? composition_source : 0.0;
}

/* Loop over all the cells and add heat according to where the current heat
   source is. */
void ConvectionSimulator::add_heat(double theta, double r, bool hot)
{
  Point p; p.theta = theta; p.r=r;
  for( RegularGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
    T[cell->self()] = T[cell->self()] + (hot ? 1.0 : -1.0)*heat(p, cell->location())*dt;
}

/* Loop over all the cells and add composition, similar to add_heat*/
void ConvectionSimulator::add_composition(double theta, double r)
{
  Point p; p.theta = theta; p.r=r;
  for( RegularGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
    C[cell->self()] = C[cell->self()] + react(p, cell->location())*dt;
}

/* Loop over the cells and zero out the displacement vectors */
void ConvectionSimulator::clear_seismic_waves()
{
  for( RegularGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
  {
    D[cell->self()] = 0.;
    Dp[cell->self()] = 0.;
  }
}


/* Add an earthquake at the point x,y.  I use a functional form of a
  Mexican hat wavelet.  At some point it may be good to revisit what
  functional form is best.  */
void ConvectionSimulator::earthquake(double theta, double r)
{
  const double x1 = (r+grid.r_inner)*std::cos(theta), y1 = (r+grid.r_inner)*std::sin(theta);

  const double earthquake_radius = grid.dr*4.;  //Somewhat arbitrary radius
  const double prefactor = 2. / std::sqrt( 3. * earthquake_radius * M_PI * M_PI);

  for( RegularGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
  {
    Point p2 = cell->location();
    const double x2 = (p2.r+grid.r_inner)*std::cos(p2.theta), y2 = (p2.r+grid.r_inner)*std::sin(p2.theta);
    const double rsq = (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);
    const double dist = rsq/earthquake_radius/earthquake_radius;
    D[cell->self()] += prefactor * (1.0 - dist) * std::exp( -dist/2. );
  }
}

/* Explicitly propagate the wave equation. */
void ConvectionSimulator::propagate_seismic_waves()
{
  const double reference_speed = 1.0;  //dummy wavespeed
  const double tstep = 0.6*dmin(grid.r_inner*grid.dtheta,grid.dr)/reference_speed; // Timestep to satisfy cfl
  double dissipation = 0.5;  //Empirically chosen dissipation

  //We can get away with an explicit timestepping scheme for the wave equation,
  //so no complicated tridiagonal matrix inversion or any such nonsense here
  for( RegularGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
  {
    double laplacian;
    const double r = cell->location().r + grid.r_inner;
    const double dr = grid.dr;
    const double dtheta = grid.dtheta;

    //Make the wavespeed temperature dependent.  This is a HUGE
    //temperature dependence so it is pretty obvious.
    double speed = reference_speed*(1.0 - 0.7*T[cell->self()]);
    if(speed < 0.0) speed = 0.0;

    //Five point laplacian stencil
    if( cell->at_top_boundary() )
      laplacian = (1./dr/dr + 0.5/r/dr)*D[cell->down()] + (1./dr/dr - 0.5/r/dr)*D[cell->down()]
                  + 1./r/r/dtheta/dtheta * ( D[cell->right()] + D[cell->left()] ) -
                  (2./r/r/dtheta/dtheta + 2./dr/dr)*D[cell->self()];
    else if ( cell->at_bottom_boundary() )
      laplacian = (1./dr/dr + 0.5/r/dr)*D[cell->up()] + (1./dr/dr - 0.5/r/dr)*D[cell->up()]
                  + 1./r/r/dtheta/dtheta * ( D[cell->right()] + D[cell->left()] ) -
                  (2./r/r/dtheta/dtheta + 2./dr/dr)*D[cell->self()];
    else
      laplacian = (1./dr/dr + 0.5/r/dr)*D[cell->up()] + (1./dr/dr - 0.5/r/dr)*D[cell->down()]
                  + 1./r/r/dtheta/dtheta * ( D[cell->right()] + D[cell->left()] ) -
                  (2./r/r/dtheta/dtheta + 2./dr/dr)*D[cell->self()];

    //Explicitly propagate wave
    scratch[cell->self()] = (2.0 * D[cell->self()] - Dp[cell->self()]
                             + tstep*tstep*speed*speed*laplacian + tstep*dissipation*D[cell->self()])
                             /(1.0 + dissipation*tstep);
  }

  for( RegularGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
  {
    //Current displacement becomes previous displacement
    Dp[cell->self()] = D[cell->self()];
    //Copy over new displacement
    D[cell->self()] = scratch[cell->self()];

    //clip field to tamp down on instabilities and keep it in reasonable values
    if(D[cell->self()] > 1.0) D[cell->self()] = 1.0;
    if(D[cell->self()] < -1.0) D[cell->self()] = -1.0;
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

/* Interpolate the composition onto an arbitrary point. */
inline double ConvectionSimulator::evaluate_composition(const Point &p)
{
  if (p.r < 0.0 || p.r > grid.lr ) return 0.0;

  double comp;
  RegularGrid::iterator cell = grid.lower_left_cell(p);
  double local_theta = fast_fmod(p.theta, grid.dtheta)/grid.dtheta;
  double local_r = ( p.r - cell->location().r )/grid.dr;

  if (cell->at_top_boundary() )
    comp = linear_interp_2d( local_theta, local_r, C[cell->self()], C[cell->right()],
                             C[cell->self()], C[cell->right()]);
  else
    comp = linear_interp_2d( local_theta, local_r, C[cell->up()], C[cell->upright()],
                             C[cell->self()], C[cell->right()]);

  return comp;
}

void ConvectionSimulator::semi_lagrangian_advect_temperature()
{
  advection_field field = temperature;
  semi_lagrangian_advect( field );
}

void ConvectionSimulator::semi_lagrangian_advect_composition()
{
  advection_field field = composition;
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

  Point vel_final, vel_final_cart;  //Velocity at the grid point
  Point vel_takeoff, vel_takeoff_cart; //Velocity at the candidate takeoff point
  Point takeoff_point, takeoff_point_cart; //Candidate takeoff point
  Point final_point, final_point_cart; //grid point


  for( RegularGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
  {
    //These points are known, as they are the grid points in question.  They will
    //not change for this cell.
    unsigned int cell_index = cell->self();
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
    else if (field == composition)
      scratch[cell->self()] = evaluate_composition(takeoff_point);  //Store the composition we found
  }

  //Copy the scratch vector into the temperature or composition vector.
  if (field == temperature)
    for( RegularGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
      T[cell->self()] = scratch[cell->self()];
  else if (field == composition)
    for( RegularGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
      C[cell->self()] = scratch[cell->self()];
}

/* Solve the diffusion equation implicitly with reverse Euler.
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
  fftw_execute(dft_diffusion);  //X direction

  for (unsigned int l = 0; l <= grid.ntheta/2.; ++l)
    diffusion_matrices[l]->solve( T_spectral+l, grid.ntheta, scratch1_spectral+l);

  //Execute the inverse Fourier transform
  fftw_execute(idft_diffusion);  //X direction

  //Renormalize
  for( RegularGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
    T[cell->self()] = scratch[cell->self()]/grid.ntheta;
}


void ConvectionSimulator::setup_diffusion_problem()
{
  //Setup the tridiagonal matrices in the y direction
  double *upper_diag = new double[grid.nr];
  double *diag = new double[grid.nr];
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
  for (unsigned int l = 0; l <= grid.ntheta/2.; ++l)
  {
     double factor = (l*l);
     diag[0] = 1.0; diag[grid.nr-1] = 1.0;
     for (unsigned int i=1; i<grid.nr-1; ++i)
     {
       const double r = grid.r_inner + i*dr;
       diag[i] = 1. + 2.*dt/dr/dr + factor*dt/r/r;
     }

     diffusion_matrices[l]->initialize(lower_diag, diag, upper_diag);
  }

  delete[] upper_diag;
  delete[] diag;
  delete[] lower_diag;

  int n[1];
  int stride, dist, howmany;

  //transform in the x direction with a full dft
  n[0] = grid.ntheta; stride = 1; dist = grid.ntheta; howmany = grid.nr;
  dft_diffusion = fftw_plan_many_dft_r2c(1, n, howmany, T, NULL, stride, dist,
                                         reinterpret_cast<fftw_complex*>(T_spectral),
                                         NULL, stride, dist, FFTW_ESTIMATE);
  idft_diffusion = fftw_plan_many_dft_c2r(1, n, howmany,
                                          reinterpret_cast<fftw_complex*>(scratch1_spectral),
                                          NULL, stride, dist, scratch, NULL,
                                          stride, dist, FFTW_ESTIMATE);

}

/* Setup the stokes solve with the stream function formulation.
   We assemble the vector that stores the frequencies of the
   eigenmodes, as well as tell FFTW how to do the transforms */
void ConvectionSimulator::setup_stokes_problem()
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

  for (unsigned int l = 0; l <= grid.ntheta/2.; ++l)
  {
     double factor = (l*l);
     diag[0] = 1.0; diag[grid.nr-1] = 1.0;
     for (unsigned int i=1; i<grid.nr-1; ++i)
     {
       const double r = grid.r_inner + i*dr;
       diag[i] = -2./dr/dr - factor/r/r;
     }

     stokes_matrices[l]->initialize(lower_diag, diag, upper_diag);
  }

  delete[] upper_diag;
  delete[] diag;
  delete[] lower_diag;

  int n[1];
  int stride, dist, howmany;

  //transform in the x direction with a full dft
  n[0] = grid.ntheta; stride = 1; dist = grid.ntheta; howmany = grid.ny;
  dft_stokes = fftw_plan_many_dft_r2c(1, n, howmany, curl_density, NULL, stride, dist,
                                      reinterpret_cast<fftw_complex*>(curl_density_spectral),
                                      NULL, stride, dist, FFTW_ESTIMATE);
  idft_stokes = fftw_plan_many_dft_c2r(1, n, howmany,
                                       reinterpret_cast<fftw_complex*>(scratch2_spectral),
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


//The curl of the temperature is what is relevant for the stream
//function calculation.  This calculates that curl.
void ConvectionSimulator::assemble_curl_density_vector()
{
  //Assemble curl_density vector
  if (include_composition)
  {
    for( RegularGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
    {
      const double r = cell->location().r + grid.r_inner;
      double dTdtheta, dCdtheta;

      if(cell->at_bottom_boundary() || cell->at_top_boundary())
      {
        dTdtheta = 0.0;
        dCdtheta = 0.0;
      }
      else
      {
        dTdtheta = (T[cell->right()]-T[cell->left()])/2./grid.dtheta;
        dCdtheta = (C[cell->right()]-C[cell->left()])/2./grid.dtheta;
      }

      curl_density[cell->self()] = Ra*( dTdtheta - buoyancy_number * dCdtheta )/r;
    }
  }
  else
  {
    for( RegularGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
    {
      const double r = cell->location().r + grid.r_inner;
      if(cell->at_bottom_boundary() || cell-> at_top_boundary())
        curl_density[cell->self()] = 0.0;
      else
        curl_density[cell->self()] = Ra*(T[cell->right()] - T[cell->left()])
                                   /r/2.0/grid.dtheta;
    }
  }
}


/* Actually solve the biharmonic equation for the Stokes system.*/
void ConvectionSimulator::solve_stokes()
{
  //Come up with the RHS of the spectral solve
  assemble_curl_density_vector();

  //Execute the forward fourier transform
  fftw_execute(dft_stokes);  //X direction

  for (unsigned int l = 0; l <= grid.ntheta/2.; ++l)
  {
     stokes_matrices[l]->solve( curl_density_spectral+l, grid.ntheta, scratch1_spectral+l);
     stokes_matrices[l]->solve( scratch1_spectral+l, grid.ntheta, scratch2_spectral+l);
  }

  //Execute the inverse Fourier transform
  fftw_execute(idft_stokes);  //X direction

  //Renormalize
  for( RegularGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
    stream[cell->self()] = scratch[cell->self()]/grid.ntheta;

  //Come up with the velocities by taking finite differences of the stream function.
  //I could also take these derivatives in spectral space, but that would mean
  //more fourier transforms, so this should be considerably cheaper.
  for( RegularGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
  {
    const double r = cell->location().r + grid.r_inner;
    if( cell->at_top_boundary() )
      u[cell->self()] = (stream[cell->self()] - stream[cell->down()])/grid.dr;
    else if ( cell->at_bottom_boundary() )
      u[cell->self()] = (stream[cell->up()] - stream[cell->self()])/grid.dr;
    else
      u[cell->self()] = (stream[cell->up()] - stream[cell->down()])/2.0/grid.dr;

    v[cell->self()] = -(stream[cell->right()] - stream[cell->left()])/2.0/grid.dtheta/r;
  }


}


void ConvectionSimulator::update_state(double rayleigh)
{
  double Ra_c = 657.;  //critical rayleigh number
  double length_scale = std::pow(rayleigh/2./Ra_c, -1./3.)*grid.lr;  //calculate a provisional length scale

  //The resolution basically sets the maximum Ra we can use.  Estimate the minimum length scale,
  //and if that is smaller than the resolution, cap the Rayleigh number.
  const double boundary_layer_cells = 4.;  //grid cells per boundary layer
  if (length_scale < boundary_layer_cells *grid.dr)
    Ra = 2.*Ra_c*std::pow( grid.lr/boundary_layer_cells/grid.dr, 3.0);
  else Ra = rayleigh;

  length_scale = std::pow(Ra/2./Ra_c, -1./3.)*grid.lr;
  const double Nu = std::pow(Ra/Ra_c/2., 1./3.) / 2.;  //Nusselt
  const double velocity_scale = std::sqrt( Ra * Nu );
  const double cfl = dmin(grid.dr, grid.r_inner*grid.dtheta)/velocity_scale;

  //Estimate other state properties based on simple isoviscous scalings
  dt = cfl * 20.0; //Roughly 10x CFL, thanks to semi-lagrangian
  heat_source_radius = length_scale*0.5;  //Radius of order the boundary layer thickness
  heat_source = velocity_scale/grid.lr*2.; //Heat a blob of order the ascent time for thta blob
  composition_source_radius = grid.lr/10.; //Size of compositional blob
  composition_source=heat_source; //Strength of compositional reaction

  setup_diffusion_problem();  //need to recompute the auxiliary vectors for the diffusion problem

}


double ConvectionSimulator::rayleigh_number() const
{
  return Ra;
}

double ConvectionSimulator::timescale() const
{
  double Ra_c = 657.;  //critical rayleigh number
  const double Nu = std::pow(Ra/Ra_c/2., 1./3.) / 2.;  //Nusselt
  const double velocity_scale = std::sqrt( Ra * Nu );
  return grid.lr/velocity_scale; //Approximately the ascent time for a plume
}
