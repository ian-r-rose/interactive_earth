#include "convection.h"
#include <cmath>
#include <iostream>


inline double fast_fmod(double x,double y) { return x-((int)(x/y))*y; }


/*Constructor for the solver. The parameters are in order:
  lx : domain size in x direction
  ly : domain size in y direction
  nx : number of cells  in x direction
  ny : number of cells  in y direction
  Rayleigh : initial Rayleigh number
*/
ConvectionSimulator::ConvectionSimulator( double lx, double ly, int nx, int ny, double Rayleigh):
                          nx(nx), ny(ny), ncells(nx*ny), Ra(Rayleigh),
                          grid(lx, ly, nx, ny), 
                          theta(0.0)
{
  //Allocate memory for data vectors
  T = new double[ncells];
  C = new double[ncells];
  stream = new double[ncells];
  curl_density = new double[ncells];
  lux = new double[ncells];
  luy = new double[ncells];
  u = new double[ncells];
  v = new double[ncells];
  g = new double[ncells];
  freqs = new double[ncells];
  scratch1 = new double[ncells];
  scratch2 = new double[ncells];

  //I actually allocate more memory than necessary for the 
  //spectral vector, but this way it makes the indexing simpler
  curl_density_spectral = new fftw_complex[ncells];

  //Initialize the state
  buoyancy_number = 1.0;
  update_state(Rayleigh, theta);
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
  delete[] C;
  delete[] stream;
  delete[] curl_density;
  delete[] lux;
  delete[] luy;
  delete[] u;
  delete[] v;
  delete[] g;
  delete[] freqs;
  delete[] scratch1;
  delete[] scratch2;
  delete[] curl_density_spectral;
 
  fftw_destroy_plan(dst);
  fftw_destroy_plan(idst);
  fftw_destroy_plan(dft);
  fftw_destroy_plan(idft);
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
  for( StaggeredGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
    T[cell->self()] = initial_temperature(cell->center());
}
 
/* Given a heat source centered on p1, calculate the heating at
   point p2, using a simple Gaussian source term */ 
inline double ConvectionSimulator::heat(const Point &p1, const Point &p2 )
{
  const double rsq = (p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y);
  return heat_source*std::exp( -rsq/2.0/heat_source_radius/heat_source_radius );
}

/* Given a heat source centered on p1, calculate the heating at
   point p2, using a simple Gaussian source term */ 
inline double ConvectionSimulator::react(const Point &p1, const Point &p2 )
{
  const double rsq = (p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y);
  return rsq < 0.1*0.1 ? heat_source : 0.0;
}

/* Loop over all the cells and add heat according to where the current heat
   source is. */
void ConvectionSimulator::add_heat(double x, double y, bool hot)
{
  Point p; p.x = x; p.y=y;
  for( StaggeredGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
    T[cell->self()] = T[cell->self()] + (hot ? 1.0 : -1.0)*heat(p, cell->center())*dt;
}
  
/* Loop over all the cells and add composition, similar to add_heat*/
void ConvectionSimulator::add_composition(double x, double y)
{
  Point p; p.x = x; p.y=y;
  for( StaggeredGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
    C[cell->self()] = C[cell->self()] + react(p, cell->center())*dt;
}



/* Interpolate the velocity onto an arbitrary point
   A bit complicated due to the staggered nature of the grid */
Point ConvectionSimulator::velocity(const Point &p)
{
  Point vel;

  //Get the relevant lower-left cells for the x and y velocities
  StaggeredGrid::iterator x_cell = grid.lower_left_vface_cell(p); 
  StaggeredGrid::iterator y_cell = grid.lower_left_hface_cell(p); 

  //Determine the local x and y coordinates for the x velocity
  double vx_local_x = fast_fmod(p.x, grid.dx)/grid.dx;
  double vx_local_y = (p.y - x_cell->vface().y)/grid.dy;

  //Determine the local x and y coordinates for the y velocity
  double vy_local_x = fast_fmod(p.x + grid.dx*0.5, grid.dx)/grid.dx;
  double vy_local_y = (p.y - y_cell->hface().y)/grid.dy;

  //get interpolated vx
  if(x_cell->at_top_boundary())
    vel.x = linear_interp_2d (vx_local_x, vx_local_y, u[x_cell->self()], 
                              u[x_cell->right()], u[x_cell->self()], u[x_cell->right()]);
  else if(x_cell->at_bottom_boundary() && vx_local_y < 0.0)
    vel.x = linear_interp_2d (vx_local_x, 1.0+vx_local_y, u[x_cell->self()], 
                              u[x_cell->right()], u[x_cell->self()], u[x_cell->right()]);
  else
    vel.x = linear_interp_2d (vx_local_x, vx_local_y, u[x_cell->up()], u[x_cell->upright()],
                              u[x_cell->self()], u[x_cell->right()]);

  //get interpolated vy
  if(y_cell->at_top_boundary())
    vel.y = linear_interp_2d (vy_local_x, vy_local_y, 0.0, 0.0 ,
                              v[y_cell->self()], v[y_cell->right()]);
  else
    vel.y = linear_interp_2d (vy_local_x, vy_local_y, v[y_cell->up()], v[y_cell->upright()],
                              v[y_cell->self()], v[y_cell->right()]);


  return vel;
} 

/* Interpolate the temperature onto an arbitrary point. */
inline double ConvectionSimulator::evaluate_temperature(const Point &p)
{
  double value;
  StaggeredGrid::iterator cell = grid.lower_left_center_cell(p); 
  double local_x = fast_fmod(p.x + grid.dx*0.5, grid.dx)/grid.dx;
  double local_y = ( p.y - cell->center().y )/grid.dy;

  if (cell->at_top_boundary() )
    value = linear_interp_2d( local_x, local_y, -T[cell->self()], -T[cell->right()],  
                             T[cell->self()], T[cell->right()]);
  else if (cell->at_bottom_boundary() && local_y < 0.0)
    value = linear_interp_2d( local_x, local_y-grid.dy, 
                             T[cell->self()], T[cell->right()], 2.0-T[cell->self()], 2.0-T[cell->right()]);
  else
    value = linear_interp_2d( local_x, local_y, T[cell->up()], T[cell->upright()],
                             T[cell->self()], T[cell->right()]);

  return value;
}

/* Interpolate the temperature onto an arbitrary point. */
inline double ConvectionSimulator::evaluate_composition(const Point &p)
{
  if (p.y < 0.0 || p.y > 1.0 ) return 0.0;

  double value;
  StaggeredGrid::iterator cell = grid.lower_left_center_cell(p); 
  double local_x = (p.x - cell->center().x)/grid.dx;
  double local_y = (p.y - cell->center().y)/grid.dy;

  if (cell->at_top_boundary() )
    value = linear_interp_2d( local_x, local_y, C[cell->self()], C[cell->right()],  
                             C[cell->self()], C[cell->right()]);
  else if (cell->at_bottom_boundary() && local_y < 0.0)
   value = linear_interp_2d( local_x, local_y-grid.dy, 
                             C[cell->self()], C[cell->right()], C[cell->self()], C[cell->right()]);
  else
    value = linear_interp_2d( local_x, local_y, C[cell->up()], C[cell->upright()],
                             C[cell->self()], C[cell->right()]);

  return value;
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
void ConvectionSimulator::semi_lagrangian_advect( advection_field field)
{
  //The goal is to find the temperature at the Lagrangian point which will 
  //be advected to the current grid point in one time step.  In general,
  //this point wil not be on the grid.  I find this point using a coarse
  //iterated predictor corrector.  

  Point vel_final;  //Velocity at the grid point
  Point vel_takeoff; //Velocity at the candidate takeoff point
  Point takeoff_point; //Candidate takeoff point
  Point final_point; //grid point

  double *F;
  if( field == temperature ) F = T;
  else if( field == composition ) F = C;

  for( StaggeredGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
  {
    //These points are known, as they are the grid points in question.  They will
    //not change for this cell.
    vel_final = velocity (cell->center() );
    final_point = cell->center();

  
    //Calculate the initial guess for the takeoff point using a basic predictor
    //with a forward euler step
    takeoff_point.x = final_point.x - vel_final.x*dt;
    takeoff_point.y = final_point.y - vel_final.y*dt;
    //Keep it in the domain
    takeoff_point.x = fast_fmod(fast_fmod(takeoff_point.x, grid.lx) + grid.lx, grid.lx);
    takeoff_point.y = std::min( grid.ly, std::max( takeoff_point.y, 0.0) );

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
      takeoff_point.x = fast_fmod(fast_fmod(takeoff_point.x, grid.lx) + grid.lx, grid.lx);
      takeoff_point.y = std::min( grid.ly, std::max( takeoff_point.y, 0.0) );
    } 
    if (field == temperature ) 
      scratch1[cell->self()] = evaluate_temperature(takeoff_point);  //Store the temperature we found
    else if (field == composition ) 
      scratch1[cell->self()] = evaluate_composition(takeoff_point);  //Store the composition we found
  }
   
  //Copy the scratch vector into the field vector. 
  for( StaggeredGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
    F[cell->self()] = scratch1[cell->self()];
 
  clip_field(field, 0.0, 1.0);
}

void ConvectionSimulator::clip_field( advection_field field, double min, double max)
{
  double *F;
  if( field == temperature ) F = T;
  else if( field == composition ) F = C;

  for( StaggeredGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
  {
    double val = F[cell->self()];
    F[cell->self()] = ( val > max ? max : (val < min ? min: val ) );
  }
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

  //First do diffusion in the y direction
  double eta_y = dt/grid.dy/grid.dy;
  double alpha_y = -eta_y;
  double lambda_y = (1.0 + 2.0*eta_y);
  //Forward sweep
  for( StaggeredGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
  {
    if(cell->at_bottom_boundary())
      scratch1[cell->self()] = 1.0;
    else if (cell->at_top_boundary())
      scratch1[cell->self()] = 0.0;
    else
      scratch1[cell->self()] = (T[cell->self()] - scratch1[cell->down()]*alpha_y)/(lambda_y-luy[cell->down()]*alpha_y);
  }
  //Do reverse sweep for y direction
  for( StaggeredGrid::reverse_iterator cell = grid.rbegin(); cell != grid.rend(); ++cell)
  {
    if(cell->at_top_boundary())
      T[cell->self()] = 0.0;
    if(cell->at_bottom_boundary())
      T[cell->self()] = 1.0;
    else
      T[cell->self()] = scratch1[cell->self()] - luy[cell->self()]*T[cell->up()];
  }


  //Now do the x direction, which is more complicated due to the periodicity.
  double eta_x = dt/grid.dx/grid.dx;
  double alpha_x = -eta_x;
  double lambda_x = (1.0 + 2.0*eta_x);

  //First sweep!
  for( StaggeredGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
  {
    if(cell->at_left_boundary())
      scratch1[cell->self()] = T[cell->self()]/lambda_x;
    else if (cell->at_right_boundary())
      scratch1[cell->self()] = 0.0; 
    else
      scratch1[cell->self()] = (T[cell->self()]-scratch1[cell->left()]*alpha_x)/(lambda_x-lux[cell->left()]*alpha_x);
  }
  //Do reverse sweep for x direction
  StaggeredGrid::reverse_iterator trailer = grid.rbegin();
  for( StaggeredGrid::reverse_iterator cell = grid.rbegin(); cell != grid.rend(); ++cell)
  {
    if(cell->at_right_boundary())
    {
      ++cell;
      scratch2[cell->self()] = scratch1[cell->self()];
    }
    else
      scratch2[cell->self()] = scratch1[cell->self()] - lux[cell->self()]*scratch2[cell->right()];
 
    if(cell->at_left_boundary())
    {
      double xn = (T[trailer->self()] - alpha_x*(scratch2[cell->self()] + scratch2[trailer->self()-1]))/(lambda_x-gamma);
      T[trailer->self()] = xn;
      ++trailer;
      for(; (trailer->at_right_boundary()==false && trailer != grid.rend()); ++trailer)
        T[trailer->self()] = scratch2[trailer->self()] - xn*g[trailer->self()];
    }
  }
}

  
void ConvectionSimulator::setup_diffusion_problem()
{

  double eta_x = dt/grid.dx/grid.dx;
  double eta_y = dt/grid.dy/grid.dy;
  double alpha_x = -eta_x;
  double lambda_x = (1.0 + 2.0*eta_x);
  double alpha_y = -eta_y;
  double lambda_y = (1.0 + 2.0*eta_y);

  for( StaggeredGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
  {
    //assemble condensed lux vector
    if(cell->at_left_boundary())
    {
      lux[cell->self()] = alpha_x/lambda_x;
      scratch2[cell->self()] = alpha_x/lambda_x;
    }
    else if (cell->at_right_boundary())
    {
      lux[cell->self()] = 0.0; lux[cell->left()] = 0.0; //solving condensed system, so no last DOF
      scratch2[cell->self()] = 0.0; 
      scratch2[cell->left()] = (alpha_x-scratch2[cell->left()-1]*alpha_x)/(lambda_x-lux[cell->left()-1]*alpha_x);
    }
    else
    {
      lux[cell->self()] = alpha_x/(lambda_x - lux[cell->left()]*alpha_x);
      scratch2[cell->self()] = (-scratch2[cell->left()]*alpha_x)/(lambda_x-lux[cell->left()]*alpha_x);
    }

    //Assemble luy vector
    if(cell->at_bottom_boundary())
      luy[cell->self()] = 0.0;
    else if (cell->at_top_boundary())
      luy[cell->self()] = 0.0;
    else
      luy[cell->self()] = alpha_y/(lambda_y - luy[cell->down()]*alpha_y);
  }
  //Do reverse sweep for x direction
  for( StaggeredGrid::reverse_iterator cell = grid.rbegin(); cell != grid.rend(); ++cell)
  {
    if(cell->at_right_boundary())
    {
      ++cell;
      g[cell->self()] = scratch2[cell->self()];
    }
    else
      g[cell->self()] = scratch2[cell->self()] - lux[cell->self()]*g[cell->right()];
 
    if(cell->at_left_boundary())
      gamma = alpha_x*(g[cell->self()] + g[cell->left()-1]);
  }
}

/* Setup the stokes solve with the stream function formulation.
   We assemble the vector that stores the frequencies of the
   eigenmodes, as well as tell FFTW how to do the transforms */
void ConvectionSimulator::setup_stokes_problem()
{

  for( StaggeredGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
  {
    int l = cell->xindex();
    int m = cell->yindex()+1;
    double factor = 1.0/(4.0*M_PI*M_PI*l*l/grid.lx/grid.lx + M_PI*M_PI*m*m/grid.ly/grid.ly);
    factor = factor*factor;
    freqs[cell->self()]=factor;
  }

  int n[1];
  int stride, dist, howmany;

  //transform in the y direction with a dst
  n[0] = grid.ny; stride = grid.nx; dist = 1; howmany = grid.nx; 
  fftw_r2r_kind f[1] = {FFTW_RODFT00};
  dst = fftw_plan_many_r2r(1, n, howmany, curl_density, NULL, stride, dist,
                                     scratch1, NULL, stride, dist,
                                     f, FFTW_ESTIMATE);
  idst = fftw_plan_many_r2r(1, n, howmany, scratch1, NULL, stride, dist,
                                     scratch2, NULL, stride, dist,
                                     f, FFTW_ESTIMATE);

  //transform in the x direction with a full dft
  n[0] = grid.nx; stride = 1; dist = grid.nx; howmany = grid.ny; 
  dft = fftw_plan_many_dft_r2c(1, n, howmany, scratch1, NULL, stride, dist,
                     curl_density_spectral, NULL, stride, dist, FFTW_ESTIMATE);
  idft = fftw_plan_many_dft_c2r(1, n, howmany, curl_density_spectral, NULL, stride, dist,
                     scratch1, NULL, stride, dist, FFTW_ESTIMATE);

}

  
/*Calculate the Nusselt number for the given temperature field.
  In principle I could calculate the heat flux through the bottom
  boundary as well, and take the average, but as it is, I am only
  doing the heat flux through the top boundary*/
double ConvectionSimulator::nusselt_number() 
{
  double heat_flux;

  for ( StaggeredGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
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
void ConvectionSimulator::assemble_curl_density_vector()
{
  //Assemble curl_density vector
  for( StaggeredGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
  {
    double dTdx, dTdy, dCdx, dCdy, curl_T, curl_C;
    dCdx = C[cell->self()] - C[cell->left()];

    if(!cell->at_bottom_boundary())
    {
      dTdx = (T[cell->self()] - T[cell->left()] + T[cell->down()] - T[cell->downleft()])/2.0/grid.dx;
      dCdx = (C[cell->self()] - C[cell->left()] + C[cell->down()] - C[cell->downleft()])/2.0/grid.dx;
      dTdy = (T[cell->left()] - T[cell->downleft()] + T[cell->self()] - T[cell->down()])/2.0/grid.dy;
      dCdy = (C[cell->left()] - C[cell->downleft()] + C[cell->self()] - C[cell->down()])/2.0/grid.dy;
    }
    else
    {
      dTdx = 0.0;
      dCdx = (C[cell->self()] - C[cell->left()])/grid.dx;
      dTdy = (T[cell->left()] + T[cell->self()] - 0.0 )/grid.dy/2.0/2.0;
      dCdy = 0.0;
    }

    curl_T = (std::cos(theta*M_PI/180.0) * dTdx - std::sin(theta*M_PI/180.0) * dTdy);
    curl_C = (std::cos(theta*M_PI/180.0) * dCdx - std::sin(theta*M_PI/180.0) * dCdy);

    curl_density[cell->self()] = Ra * (curl_T - buoyancy_number * curl_C);
  }
}


/* Actually solve the biharmonic equation for the Stokes system.*/
void ConvectionSimulator::solve_stokes()
{ 
  //Come up with the RHS of the spectral solve
  assemble_curl_density_vector();

  //Execute the forward fourier transform
  fftw_execute(dst);  //Y direction
  fftw_execute(dft);  //X direction

  //Multiply by the relevant frequencies.
  for( StaggeredGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
  {
    if(cell->xindex() <= grid.nx/2)
    {
      curl_density_spectral[cell->self()][0]*=freqs[cell->self()];
      curl_density_spectral[cell->self()][1]*=freqs[cell->self()];
    }
  }

  //Execute the inverse Fourier transform
  fftw_execute(idft);  //Y direction
  fftw_execute(idst);  //X direction

  //Renormalize
  for( StaggeredGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
    stream[cell->self()] = scratch2[cell->self()]/(grid.ny+1)/2.0/grid.nx;

  //Come up with the velocities by taking finite differences of the stream function.
  //I could also take these derivatives in spectral space, but that would mean 
  //more fourier transforms, so this should be considerably cheaper.
  for( StaggeredGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
  {
    if( cell->at_top_boundary() == false)
      u[cell->self()] = (stream[cell->up()] - stream[cell->self()])/grid.dy;
    else
      u[cell->self()] = (stream[cell->self()] - stream[cell->down()])/grid.dy;

    v[cell->self()] = -(stream[cell->right()] - stream[cell->self()])/grid.dx;
  }


}


void ConvectionSimulator::update_state(double rayleigh, double gravity_angle)
{
  theta = gravity_angle; //update the angle of gravity
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
  dt = cfl * 5.0; //Roughly 10x CFL, thanks to semi-lagrangian
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
