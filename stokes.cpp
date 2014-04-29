#include "stokes.h"
#include <cmath>

StokesSolver::StokesSolver( double lx, double ly, int nx, int ny, double Rayleigh):
                          nx(nx), ny(ny), ncells(nx*ny), Ra(Rayleigh),
                          grid(lx, ly, nx, ny), 
                          theta(0.0)
{
  T = new double[ncells];
  vorticity = new double[ncells];
  stream = new double[ncells];
  curl_T = new double[ncells];
  lux = new double[ncells];
  luy = new double[ncells];
  u = new double[ncells];
  v = new double[ncells];
  g = new double[ncells];
  freqs = new double[ncells];
  scratch1 = new double[ncells];
  scratch2 = new double[ncells];

  curl_T_spectral = new fftw_complex[nx*ny];

  update_state(Rayleigh, theta);
  
  initialize_temperature();
  assemble_curl_T_vector();
  setup_stokes_problem();
  setup_diffusion_problem();

}

StokesSolver::~StokesSolver()
{
  delete[] T;
  delete[] vorticity;
  delete[] stream;
  delete[] curl_T;
  delete[] lux;
  delete[] luy;
  delete[] u;
  delete[] v;
  delete[] g;
  delete[] freqs;
  delete[] scratch1;
  delete[] scratch2;
  delete[] curl_T_spectral;
 
  fftw_destroy_plan(dst);
  fftw_destroy_plan(idst);
  fftw_destroy_plan(dft);
  fftw_destroy_plan(idft);
}
 

double StokesSolver::initial_temperature(const Point &p)
{
//  if (std::sqrt( (0.35-p.x)*(0.35-p.x)+(0.5-p.y)*(0.5-p.y))  < 0.05 ) return 1.0;
//  else if (std::sqrt( (1.65-p.x)*(1.65-p.x)+(0.5-p.y)*(0.5-p.y))  < 0.05 ) return 0.0;
//  else return 0.5;
  return 0.5;
}

void StokesSolver::initialize_temperature()
{
  for( StaggeredGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
    T[cell->self()] = initial_temperature(cell->center());
}
 
inline double StokesSolver::heat(const Point &p1, const Point &p2 )
{
  const double rsq = (p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y);
  return heat_source*std::exp( -rsq/2.0/heat_source_radius/heat_source_radius );
}

void StokesSolver::add_heat(double x, double y, bool hot)
{
  Point p; p.x = x; p.y=y;
  for( StaggeredGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
    T[cell->self()] = T[cell->self()] + (hot ? 1.0 : -1.0)*heat(p, cell->center())*dt;
}
  



//interpolate the velocity onto an arbitrary point
//Kind of a mess...
Point StokesSolver::velocity(const Point &p)
{
  Point vel;
  StaggeredGrid::iterator x_cell = grid.lower_left_vface_cell(p); 
  StaggeredGrid::iterator y_cell = grid.lower_left_hface_cell(p); 

  double vx_local_x = (p.x - x_cell->vface().x)/grid.dx;
  double vx_local_y = (p.y - x_cell->vface().y)/grid.dy;

  double vy_local_x = (p.x - y_cell->hface().x)/grid.dx;
  double vy_local_y = (p.y - y_cell->hface().y)/grid.dy;

  //get interpolated vx
  if(x_cell->at_top_boundary())
    vel.x = linear_interp_2d (vx_local_x, vx_local_y, u[x_cell->self()], 
                              u[x_cell->right()], u[x_cell->self()], u[x_cell->right()]);
  else if(x_cell->at_bottom_boundary() && vx_local_x < 0.0)
    vel.x = linear_interp_2d (vx_local_x, vx_local_y, u[x_cell->self()], 
                              u[x_cell->right()], u[x_cell->self()], u[x_cell->right()]);
  else
    vel.x = linear_interp_2d (vx_local_x, vx_local_y, u[x_cell->up()], u[x_cell->upright()],
                              u[x_cell->self()], u[x_cell->right()]);

  //get interpolated vy
  if(y_cell->at_top_boundary())
    vel.y = linear_interp_2d (vy_local_x, vy_local_y, v[y_cell->self()], 
                              v[y_cell->right()], v[y_cell->self()], v[y_cell->right()]);
  else
    vel.y = linear_interp_2d (vy_local_x, vy_local_y, v[y_cell->up()], v[y_cell->upright()],
                              v[y_cell->self()], v[y_cell->right()]);


  return vel;
} 

inline double StokesSolver::temperature(const Point &p)
{
  double temp;
  StaggeredGrid::iterator cell = grid.lower_left_center_cell(p); 
  double local_x = (p.x - cell->center().x)/grid.dx;
  double local_y = (p.y - cell->center().y)/grid.dy;

  if (cell->at_top_boundary() )
    temp = linear_interp_2d( local_x, local_y, 0.0, 0.0,  
                             T[cell->self()], T[cell->right()]);
  else if (cell->at_bottom_boundary() && local_y < 0.0)
    temp = linear_interp_2d( local_x, local_y-grid.dy, 
                             T[cell->self()], T[cell->right()], 1.0, 1.0);
  else
    temp = linear_interp_2d( local_x, local_y, T[cell->up()], T[cell->upright()],
                             T[cell->self()], T[cell->right()]);

  return temp;
}
  

void StokesSolver::semi_lagrangian_advect()
{

  Point vel_final;
  Point vel_takeoff;
  Point takeoff_point;
  Point final_point;

  for( StaggeredGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
  {
    vel_final.x = (u[cell->self()]+u[cell->right()])/2.0;
    vel_final.y = (v[cell->up()]+v[cell->self()])/2.0;
    final_point = cell->center();
  
    takeoff_point.x = final_point.x - vel_final.x*dt;
    takeoff_point.y = final_point.y - vel_final.y*dt;
    takeoff_point.x = (takeoff_point.x < 0.0 ? takeoff_point.x+grid.lx : 
                      (takeoff_point.x >= grid.lx ? takeoff_point.x-grid.lx : takeoff_point.x));
    takeoff_point.y = (takeoff_point.y < 0.0 ? takeoff_point.y+grid.ly : 
                      (takeoff_point.y >= grid.ly ? takeoff_point.y-grid.ly : takeoff_point.y));
    for(unsigned int i=0; i<1; ++i)
    {
      vel_takeoff = velocity(takeoff_point);
      takeoff_point.x = final_point.x - (vel_final.x + vel_takeoff.x)*dt/2.0;
      takeoff_point.y = final_point.y - (vel_final.y + vel_takeoff.y)*dt/2.0;
      takeoff_point.x = (takeoff_point.x < 0.0 ? takeoff_point.x+grid.lx : 
                        (takeoff_point.x >= grid.lx ? takeoff_point.x-grid.lx : takeoff_point.x));
      takeoff_point.y = (takeoff_point.y < 0.0 ? takeoff_point.y+grid.ly : 
                        (takeoff_point.y >= grid.ly ? takeoff_point.y-grid.ly : takeoff_point.y));
    } 
    scratch1[cell->self()] = temperature(takeoff_point);
  }
  
  for( StaggeredGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
    T[cell->self()] = scratch1[cell->self()];
}

void StokesSolver::diffuse_temperature()
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

  
void StokesSolver::setup_diffusion_problem()
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

void StokesSolver::setup_stokes_problem()
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
  dst = fftw_plan_many_r2r(1, n, howmany, curl_T, NULL, stride, dist,
                                     scratch1, NULL, stride, dist,
                                     f, FFTW_ESTIMATE);
  idst = fftw_plan_many_r2r(1, n, howmany, scratch1, NULL, stride, dist,
                                     scratch2, NULL, stride, dist,
                                     f, FFTW_ESTIMATE);

  //transform in the x direction with a full dft
  n[0] = grid.nx; stride = 1; dist = grid.nx; howmany = grid.ny; 
  dft = fftw_plan_many_dft_r2c(1, n, howmany, scratch1, NULL, stride, dist,
                     curl_T_spectral, NULL, stride, dist, FFTW_ESTIMATE);
  idft = fftw_plan_many_dft_c2r(1, n, howmany, curl_T_spectral, NULL, stride, dist,
                     scratch1, NULL, stride, dist, FFTW_ESTIMATE);

}

  

void StokesSolver::assemble_curl_T_vector()
{
  //Assemble curl_T vector
  for( StaggeredGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
  {
    if(cell->at_bottom_boundary())
      curl_T[cell->self()] = 0.0; 
    else if (cell->at_top_boundary())
      curl_T[cell->self()] = 0.0; 
    else
      curl_T[cell->self()] = Ra*std::cos(theta*M_PI/180.0)*(T[cell->self()] - T[cell->left()]
                            + T[cell->down()] - T[cell->downleft()])/2.0/grid.dx
                            - Ra*std::sin(theta*M_PI/180.0)*(T[cell->left()] - T[cell->downleft()]
                            + T[cell->self()] - T[cell->down()])/2.0/grid.dy;
  }
}

void StokesSolver::solve_stokes()
{ 
  assemble_curl_T_vector();

  fftw_execute(dst);
  fftw_execute(dft);

  for( StaggeredGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
  {
    if(cell->xindex() <= grid.nx/2)
    {
      curl_T_spectral[cell->self()][0]*=freqs[cell->self()];
      curl_T_spectral[cell->self()][1]*=freqs[cell->self()];
    }
  }

  fftw_execute(idft);
  fftw_execute(idst);

  
  for( StaggeredGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
    stream[cell->self()] = scratch2[cell->self()]/(grid.ny+1)/2.0/grid.nx;

  //Come up with the velocities
  for( StaggeredGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
  {
    if( cell->at_top_boundary() == false)
      u[cell->self()] = (stream[cell->up()] - stream[cell->self()])/grid.dy;
    else
      u[cell->self()] = (stream[cell->self()] - stream[cell->down()])/grid.dy;

    v[cell->self()] = -(stream[cell->right()] - stream[cell->self()])/grid.dx;
  }


}


void StokesSolver::update_state(double rayleigh, double gravity_angle)
{
  theta = gravity_angle;
  double length_scale = std::pow(rayleigh, -1./3.)*grid.ly;  //calculate a provisional length scale

  if (length_scale < grid.ly/grid.ny/2.0) Ra = std::pow( grid.ny*2.0, 3.0);
  else Ra = rayleigh;

  dt = grid.ly/grid.ny * std::pow(Ra,-2./3.) * 20.0; //Roughly 20x CFL, thanks to semi-lagrangian
  heat_source_radius = length_scale*5.0;  //Radius of order the boundary layer thickness
  heat_source = std::pow(Ra, 2./3.)/grid.ly*1.e0; //Heat a blob of order the ascent time for thta blob

  setup_diffusion_problem();  //need to recompute the auxiliary vectors for the diffusion problem

}


double StokesSolver::rayleigh_number() const
{
  return Ra;
}
