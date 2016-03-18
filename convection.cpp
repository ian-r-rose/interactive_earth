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
  lx : domain size in x direction
  ly : domain size in y direction
  nx : number of cells  in x direction
  ny : number of cells  in y direction
  Rayleigh : initial Rayleigh number
*/
ConvectionSimulator::ConvectionSimulator( double inner_radius, int nx, int ny, double Rayleigh):
                          nx(nx), ny(ny), ncells(nx*ny), Ra(Rayleigh),
                          grid(inner_radius, nx, ny)
{
  //Allocate memory for data vectors
  T = new double[ncells];
  D = new double[ncells];
  Dp = new double[ncells];
  stream = new double[ncells];
  curl_T = new double[ncells];
  u = new double[ncells];
  v = new double[ncells];
  scratch = new double[ncells];

  //I actually allocate more memory than necessary for the 
  //spectral vector, but this way it makes the indexing simpler
  curl_T_spectral = new std::complex<double>[ncells];
  T_spectral = new std::complex<double>[ncells];
  scratch1_spectral = new std::complex<double>[ncells];
  scratch2_spectral = new std::complex<double>[ncells];

  stokes_matrices = new TridiagonalMatrixSolver<std::complex<double> >*[grid.nx/2];
  diffusion_matrices = new TridiagonalMatrixSolver<std::complex<double> >*[grid.nx/2];
  for (unsigned int i=0; i<=grid.nx/2; ++i)
  {
    stokes_matrices[i] = new TridiagonalMatrixSolver<std::complex<double> >(grid.ny);
    diffusion_matrices[i] = new TridiagonalMatrixSolver<std::complex<double> >(grid.ny);
  }

  //Initialize the state
  update_state(Rayleigh);
  initialize_temperature();

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
  delete[] curl_T;
  delete[] u;
  delete[] v;
  delete[] scratch;
  delete[] curl_T_spectral;
  delete[] T_spectral;
  delete[] scratch1_spectral;
  delete[] scratch2_spectral;

  for (unsigned int i=0; i<=grid.nx/2; ++i)
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
}
 
/* Functional form of initial temperature field.  
   Just start with a constant value of one half.*/
double ConvectionSimulator::initial_temperature(const Point &p)
{
//  if (std::sqrt( (0.35-p.x)*(0.35-p.x)+(0.5-p.y)*(0.5-p.y))  < 0.05 ) return 1.0;
//  else if (std::sqrt( (1.65-p.x)*(1.65-p.x)+(0.5-p.y)*(0.5-p.y))  < 0.05 ) return 0.0;
//  else return 0.5;
  return 0.5;
}

/* Loop over all the cells and set the initial temperature */
void ConvectionSimulator::initialize_temperature()
{
  for( RegularGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
    T[cell->self()] = initial_temperature(cell->location());
}
 
/* Given a heat source centered on p1, calculate the heating at
   point p2, using a simple Gaussian source term */ 
inline double ConvectionSimulator::heat(const Point &p1, const Point &p2 )
{
  const double x1 = (p1.y+grid.r_inner)*std::cos(p1.x), y1 = (p1.y+grid.r_inner)*std::sin(p1.x);
  const double x2 = (p2.y+grid.r_inner)*std::cos(p2.x), y2 = (p2.y+grid.r_inner)*std::sin(p2.x);
  const double rsq = (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);
  return heat_source*std::exp( -rsq/2.0/heat_source_radius/heat_source_radius );
}

/* Loop over all the cells and add heat according to where the current heat
   source is. */
void ConvectionSimulator::add_heat(double x, double y, bool hot)
{
  Point p; p.x = x; p.y=y;
  for( RegularGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
    T[cell->self()] = T[cell->self()] + (hot ? 1.0 : -1.0)*heat(p, cell->location())*dt;
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

  const double earthquake_radius = grid.dy*4.;  //Somewhat arbitrary radius
  const double prefactor = 2. / std::sqrt( 3. * earthquake_radius * M_PI * M_PI);

  for( RegularGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
  {
    Point p2 = cell->location();
    const double x2 = (p2.y+grid.r_inner)*std::cos(p2.x), y2 = (p2.y+grid.r_inner)*std::sin(p2.x);
    const double rsq = (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);
    const double dist = rsq/earthquake_radius/earthquake_radius;
    D[cell->self()] += prefactor * (1.0 - dist) * std::exp( -dist/2. );
  }
}

/* Explicitly propagate the wave equation. */
void ConvectionSimulator::propagate_seismic_waves()
{
  const double reference_speed = 1.0;  //dummy wavespeed
  const double dt = 0.6*dmin(grid.dx,grid.dy)/reference_speed; // Timestep to satisfy cfl
  double dissipation = 0.5;  //Empirically chosen dissipation

  //We can get away with an explicit timestepping scheme for the wave equation,
  //so no complicated tridiagonal matrix inversion or any such nonsense here
  for( RegularGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
  {
    double laplacian;
    const double r = cell->location().y + grid.r_inner;
    const double dr = grid.dy;
    const double dtheta = grid.dx;

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
                             + dt*dt*speed*speed*laplacian + dt*dissipation*D[cell->self()])
                             /(1.0 + dissipation*dt); 
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
Point ConvectionSimulator::velocity(const Point &p)
{
  Point vel;

  //Get the relevant lower-left cells for the velocities
  RegularGrid::iterator cell = grid.lower_left_cell(p); 

  //Determine the local x and y coordinates for the velocity
  double local_x = fast_fmod(p.x, grid.dx)/grid.dx;
  double local_y = (p.y - cell->location().y)/grid.dy;

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
  vel.x = linear_interp_2d (local_x, local_y, u[ul], u[ur],
                            u[dl], u[dr]);
  vel.y = linear_interp_2d (local_x, local_y, v[ul], v[ur],
                            v[dl], v[dr]);

  return vel;
} 

/* Interpolate the temperature onto an arbitrary point. */
inline double ConvectionSimulator::temperature(const Point &p)
{
  double temp;
  RegularGrid::iterator cell = grid.lower_left_cell(p); 
  double local_x = fast_fmod(p.x, grid.dx)/grid.dx;
  double local_y = ( p.y - cell->location().y )/grid.dy;

  if (cell->at_top_boundary() )
    temp = linear_interp_2d( local_x, local_y, -T[cell->self()], -T[cell->right()],  
                             T[cell->self()], T[cell->right()]);
  else 
    temp = linear_interp_2d( local_x, local_y, T[cell->up()], T[cell->upright()],
                             T[cell->self()], T[cell->right()]);

  return temp;
}

/* Advect the temperature field through the velocity field using
   Semi-lagrangian advection.  This scheme is quite stable, which
   allows me to take VERY large time steps.  The drawback is that it
   is kind of slow.  Here I do it in a very coarse way to make up for
   that.  A more accurate implementation would use better-than-linear
   interpolation and use more iterations. */
void ConvectionSimulator::semi_lagrangian_advect()
{
  //The goal is to find the temperature at the Lagrangian point which will 
  //be advected to the current grid point in one time step.  In general,
  //this point wil not be on the grid.  I find this point using a coarse
  //iterated predictor corrector.  

  Point vel_final;  //Velocity at the grid point
  Point vel_takeoff; //Velocity at the candidate takeoff point
  Point takeoff_point; //Candidate takeoff point
  Point final_point; //grid point

  for( RegularGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
  {
    //These points are known, as they are the grid points in question.  They will
    //not change for this cell.
    unsigned int cell_index = cell->self();
    vel_final.x = u[cell_index];
    vel_final.y = v[cell_index];
    final_point = cell->location();

  
    //Calculate the initial guess for the takeoff point using a basic predictor
    //with a forward euler step
    takeoff_point.x = final_point.x - vel_final.x*dt;
    takeoff_point.y = final_point.y - vel_final.y*dt;
    //Keep it in the domain
    takeoff_point.x = fast_fmod(takeoff_point.x + grid.lx, grid.lx);
    takeoff_point.y = dmin( grid.ly, dmax( takeoff_point.y, 0.0) );

    //Iterate on the corrector.  Here I only do one iteration for
    //performance reasons, but in principle we could do more to get
    //a better estimate.
    for(unsigned int i=0; i<1; ++i)
    {
      //Evaluate the velocity at the predictor
      vel_takeoff = velocity(takeoff_point);
      //Come up with the corrector using a midpoint rule
      takeoff_point.x = final_point.x - (vel_final.x + vel_takeoff.x)*dt/2.0;
      takeoff_point.y = final_point.y - (vel_final.y + vel_takeoff.y)*dt/2.0;
      //Keep in domain
      takeoff_point.x = fast_fmod(takeoff_point.x + grid.lx, grid.lx);
      takeoff_point.y = dmin( grid.ly, dmax( takeoff_point.y, 0.0) );
    } 
    scratch[cell_index] = temperature(takeoff_point);  //Store the temperature we found
  }
   
  //Copy the scratch vector into the temperature vector. 
  for( RegularGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
    T[cell->self()] = scratch[cell->self()];
}

/* Solve the diffusion equation implicitly with reverse Euler.
   I do the x and y direction separately, so in all cases I am
   solving a series of tridiagonal matrix equations.  This allows
   me to solve it in O(N).  The tridiagonal solves are done using
   the Thomas algorithm.  The x-direction is more complicated, 
   because the periodicity makes it not strictly tridiagonal. To
   get around this, I condense the matrix with one step of Gaussian
   elimination, so that it is separated into the tridiagonal and
   non-tridiagonal parts.  This is kind of complicated and needs 
   to be better documented.  The good news is that it is fast and
   unconditionally stable.  */
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

  for (unsigned int l = 0; l <= grid.nx/2.; ++l)
    diffusion_matrices[l]->solve( T_spectral+l, grid.nx, scratch1_spectral+l);
     
  //Execute the inverse Fourier transform
  fftw_execute(idft_diffusion);  //X direction

  //Renormalize
  for( RegularGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
    T[cell->self()] = scratch[cell->self()]/grid.nx;
}

  
void ConvectionSimulator::setup_diffusion_problem()
{
  //Setup the tridiagonal matrices in the y direction
  double *upper_diag = new double[grid.ny];
  double *diag = new double[grid.ny];
  double *lower_diag = new double[grid.ny];

  const double dr = grid.dy;
  //Upper and lower diagonals are the same for each frequency
  for( unsigned int i=0; i < grid.ny; ++i )
  {
    const double r = grid.r_inner + i*dr;
    upper_diag[i] = -dt/dr/dr * (1. + 0.5*dr/r);
    lower_diag[i] = -dt/dr/dr * (1. - 0.5*dr/r);
  }
  upper_diag[0] = 0.0;  //fix for lower B.C.
  lower_diag[grid.ny-1] = 0.0; //fix for upper B.C.

  //Diagonals include a term for the frequency in the x-direction
  for (unsigned int l = 0; l <= grid.nx/2.; ++l)
  {
     double factor = (l*l);
     diag[0] = 1.0; diag[grid.ny-1] = 1.0;
     for (unsigned int i=1; i<grid.ny-1; ++i)
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
  n[0] = grid.nx; stride = 1; dist = grid.nx; howmany = grid.ny; 
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
  double *upper_diag = new double[grid.ny];
  double *diag = new double[grid.ny];
  double *lower_diag = new double[grid.ny];

  const double dr = grid.dy;
  for( unsigned int i=0; i < grid.ny; ++i )
  {
    const double r = grid.r_inner + i*dr;
    upper_diag[i] = 1./dr/dr + 0.5/dr/r;
    lower_diag[i] = 1./dr/dr - 0.5/dr/r;
  }
  upper_diag[0] = 0.0;  //fix for lower B.C.
  lower_diag[grid.ny-1] = 0.0; //fix for upper B.C.

  for (unsigned int l = 0; l <= grid.nx/2.; ++l)
  {
     double factor = (l*l);
     diag[0] = 1.0; diag[grid.ny-1] = 1.0;
     for (unsigned int i=1; i<grid.ny-1; ++i)
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
  n[0] = grid.nx; stride = 1; dist = grid.nx; howmany = grid.ny; 
  dft_stokes = fftw_plan_many_dft_r2c(1, n, howmany, curl_T, NULL, stride, dist,
                                      reinterpret_cast<fftw_complex*>(curl_T_spectral), 
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
      double grad_T = (T[cell->self()] - T[cell->down()])/grid.dy;
      heat_flux += grad_T;
    }
  }
  return -heat_flux * grid.ly/grid.nx;
}
  
      
//The curl of the temperature is what is relevant for the stream
//function calculation.  This calculates that curl.
void ConvectionSimulator::assemble_curl_T_vector()
{
  //Assemble curl_T vector
  for( RegularGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
  {
    const double r = cell->location().y + grid.r_inner;
    if(cell->at_bottom_boundary())
      curl_T[cell->self()] = 0.0; 
    else if (cell->at_top_boundary())
      curl_T[cell->self()] = 0.0; 
    else
      curl_T[cell->self()] = Ra*(T[cell->right()] - T[cell->left()])
                                 /r/2.0/grid.dx;
  }
}


/* Actually solve the biharmonic equation for the Stokes system.*/
void ConvectionSimulator::solve_stokes()
{ 
  //Come up with the RHS of the spectral solve
  assemble_curl_T_vector();

  //Execute the forward fourier transform
  fftw_execute(dft_stokes);  //X direction

  for (unsigned int l = 0; l <= grid.nx/2.; ++l)
  {
     stokes_matrices[l]->solve( curl_T_spectral+l, grid.nx, scratch1_spectral+l);
     stokes_matrices[l]->solve( scratch1_spectral+l, grid.nx, scratch2_spectral+l);
  }
     
  //Execute the inverse Fourier transform
  fftw_execute(idft_stokes);  //X direction

  //Renormalize
  for( RegularGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
    stream[cell->self()] = scratch[cell->self()]/grid.nx;

  //Come up with the velocities by taking finite differences of the stream function.
  //I could also take these derivatives in spectral space, but that would mean 
  //more fourier transforms, so this should be considerably cheaper.
  for( RegularGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
  {
    const double r = cell->location().y + grid.r_inner;
    if( cell->at_top_boundary() )
      u[cell->self()] = (stream[cell->self()] - stream[cell->down()])/grid.dy;
    else if ( cell->at_bottom_boundary() )
      u[cell->self()] = (stream[cell->up()] - stream[cell->self()])/grid.dy;
    else
      u[cell->self()] = (stream[cell->up()] - stream[cell->down()])/2.0/grid.dy;

    v[cell->self()] = -(stream[cell->right()] - stream[cell->left()])/2.0/grid.dx/r;
  }


}


void ConvectionSimulator::update_state(double rayleigh)
{
  double Ra_c = 657.;  //critical rayleigh number
  double length_scale = std::pow(rayleigh/2./Ra_c, -1./3.)*grid.ly;  //calculate a provisional length scale

  //The resolution basically sets the maximum Ra we can use.  Estimate the minimum length scale,
  //and if that is smaller than the resolution, cap the Rayleigh number.
  const double boundary_layer_cells = 4.;  //grid cells per boundary layer
  if (length_scale < boundary_layer_cells *grid.dy)
    Ra = 2.*Ra_c*std::pow( grid.ly/boundary_layer_cells/grid.dy, 3.0);
  else Ra = rayleigh;

  length_scale = std::pow(Ra/2./Ra_c, -1./3.)*grid.ly;  
  const double Nu = std::pow(Ra/Ra_c/2., 1./3.) / 2.;  //Nusselt
  const double velocity_scale = std::sqrt( Ra * Nu );
  const double cfl = grid.dy/velocity_scale;

  //Estimate other state properties based on simple isoviscous scalings
  dt = cfl * 10.0; //Roughly 10x CFL, thanks to semi-lagrangian
  heat_source_radius = length_scale*0.5;  //Radius of order the boundary layer thickness
  heat_source = velocity_scale/grid.ly*2.; //Heat a blob of order the ascent time for thta blob

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
  return grid.ly/velocity_scale; //Approximately the ascent time for a plume
}
