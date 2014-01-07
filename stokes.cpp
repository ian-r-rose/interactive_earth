#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>
#include <Amesos.h>
#include <AztecOO.h>
#include <Teuchos_TimeMonitor.hpp>

#include "stokes.h"
#include "color.h"
#include <vector>
#include <cmath>
#include <fstream>



Teuchos::RCP<Teuchos::Time> Diffusion = Teuchos::TimeMonitor::getNewCounter("Diffusion time");
Teuchos::RCP<Teuchos::Time> Advection = Teuchos::TimeMonitor::getNewCounter("Advection time");
Teuchos::RCP<Teuchos::Time> Stokes = Teuchos::TimeMonitor::getNewCounter("Stokes time");
Teuchos::RCP<Teuchos::Time> Draw = Teuchos::TimeMonitor::getNewCounter("Draw time");

StokesSolver::StokesSolver( double lx, double ly, int nx, int ny):
#ifdef EPETRA_MPI
                          Comm(MPI_COMM_WORLD),
#endif
                          nx(nx), ny(ny), ncells(nx*ny), Ra(1.0e7),
                          grid(lx, ly, nx, ny), dt(.7e-5),
                          map(ncells, 0, Comm), theta(0.0),
                          T(map), vorticity(map), stream(map), curl_T(map), 
                          u(map), v(map), g(map), lux(map), luy(map), freqs(map),
                          scratch1(map), scratch2(map), scratch3(map), scratch4(map),
                          poisson_matrix(Copy, map, 5),
                          diffusion_updown(Copy, map, 3), diffusion_leftright(Copy,map,3)
{
  initialize_temperature();
  assemble_curl_T_vector();
  assemble_stokes_matrix();
  poisson_problem.SetOperator(&poisson_matrix);

  MLPrec = new ML_Epetra::MultiLevelPreconditioner(poisson_matrix, true);
  ifpack_precon = new Ifpack_ILUT(&poisson_matrix);
  ifpack_precon->Compute();
  preconditioner = MLPrec;

  
  aztec_solver.SetAztecOption(AZ_solver, AZ_gmres);
  aztec_solver.SetAztecOption(AZ_output, AZ_none);
  aztec_solver.SetAztecOption(AZ_keep_info, 1);

  poisson_problem.SetLHS(&vorticity);
  poisson_problem.SetRHS(&curl_T);
//  aztec_solver.SetPrecOperator(preconditioner);
  aztec_solver.SetProblem(poisson_problem);

  aztec_solver.Iterate(100, 1.e-2);

  poisson_problem.SetLHS(&stream);
  poisson_problem.SetRHS(&vorticity);
  aztec_solver.SetProblem(poisson_problem);
//  aztec_solver.SetPrecOperator(preconditioner);
  aztec_solver.Iterate(100, 1.e-2);

  assemble_diffusion_matrix();

}

double StokesSolver::initial_temperature(const Point &p)
{
  if (std::sqrt( (0.35-p.x)*(0.35-p.x)+(0.5-p.y)*(0.5-p.y))  < 0.05 ) return 1.0;
  else if (std::sqrt( (1.65-p.x)*(1.65-p.x)+(0.5-p.y)*(0.5-p.y))  < 0.05 ) return 0.0;
  else return 0.5;
}

void StokesSolver::initialize_temperature()
{
  for( StaggeredGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
    T[cell->self()] = initial_temperature(cell->center());
}

//interpolate the velocity onto an arbitrary point
//Kind of a mess...
Point StokesSolver::velocity(const Point &p)
{
//  Teuchos::TimeMonitor LocalTimer(*VInterp);
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
  Teuchos::TimeMonitor LocalTimer(*Advection);

  Point vel_final;
  Point vel_takeoff;
  Point takeoff_point;
  Point final_point;
  scratch1.PutScalar(0.0);

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
  T = scratch1;
}

void StokesSolver::upwind_advect()
{
  Teuchos::TimeMonitor LocalTimer(*Advection);
 
  double dTdx;
  double dTdy;
  double vx, vy;
  scratch1.PutScalar(0.0);

  for( StaggeredGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
  {
    vx = (u[cell->self()]+u[cell->right()])/2.0;
    if (vx <= 0.0) 
      dTdx = (T[cell->right()]-T[cell->self()])/grid.dx;
    else
      dTdx = (T[cell->self()]-T[cell->left()])/grid.dx;
    if( cell->at_top_boundary())
    {
      vy = 0;
      dTdy = (0.0 - T[cell->self()])/grid.dy/2.0;
    }
    else if (cell->at_bottom_boundary())
    {
      vy = 0;
      dTdy = (T[cell->self()] - 1.0)/grid.dy/2.0;
    }
    else
    {
      vy = (v[cell->up()]+v[cell->self()])/2.0;
      if (vy <= 0.0) 
        dTdy = (T[cell->up()]-T[cell->self()])/grid.dy;
      else
        dTdy = (T[cell->self()]-T[cell->down()])/grid.dy;
    }
    
    scratch1[cell->self()] = T[cell->self()] - dt*(vx*dTdx+vy*dTdy);
  }
  T = scratch1;
}


void StokesSolver::diffuse_temperature()
{
  Teuchos::TimeMonitor LocalTimer(*Diffusion);

  //First do diffusion in the y direction
  double eta_y = dt/grid.dy/grid.dy;
  double alpha_y = -eta_y;
  double lambda_y = (1.0 + 2.0*eta_y);
  scratch1.PutScalar(0.0);
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
  scratch1.PutScalar(0.0); //c prime
  scratch2.PutScalar(0.0); //d prime

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

/*void StokesSolver::diffuse_temperature()
{
  Teuchos::TimeMonitor LocalTimer(*Diffusion);

  diffusion_ud_problem.SetRHS(&T);
  diffusion_ud_problem.SetLHS(&scratch1);
  amesos_ud_solver->Solve();
 
  diffusion_lr_problem.SetRHS(&scratch1);
  diffusion_lr_problem.SetLHS(&T);
  amesos_lr_solver->Solve();
  
  //enforce top and bottom bcs
  for( StaggeredGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
    if(cell->at_top_boundary())
      T[cell->self()] = 0.0;
    else if(cell->at_bottom_boundary())
      T[cell->self()] = 1.0;
}*/
   
  
/*void StokesSolver::assemble_diffusion_matrix()
{
  std::vector<double> values;
  std::vector<int> indices; //5 point stencil

  //Assemble poisson matrix 
  for( StaggeredGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
  {
    // Do the up/down diffusion matrix
    indices.clear();
    values.clear();
    if(cell->at_top_boundary() || cell->at_bottom_boundary())
    {
      indices.push_back(cell->self());
      values.push_back(1.0);
    }
    else
    {
      indices.push_back(cell->self()); values.push_back(1.0+dt*2.0/grid.dx/grid.dx);
      indices.push_back( cell->up() ); values.push_back(-dt*1.0/grid.dx/grid.dx);
      indices.push_back( cell->down() ); values.push_back(-dt*1.0/grid.dx/grid.dx);
    }
    int ierr =  diffusion_updown.InsertGlobalValues(cell->self(), indices.size(), &values[0], &indices[0]);
    // Do the left/right diffusion matrix
    indices.clear();
    values.clear();
    if(cell->at_top_boundary() || cell->at_bottom_boundary())
    {
      indices.push_back(cell->self());
      values.push_back(1.0);
    }
    else
    {
      indices.push_back(cell->self()); values.push_back(1.0+dt*2.0/grid.dx/grid.dx);
      indices.push_back( cell->left() ); values.push_back(-dt*1.0/grid.dx/grid.dx);
      indices.push_back( cell->right() ); values.push_back(-dt*1.0/grid.dx/grid.dx);
    }
    ierr =  diffusion_leftright.InsertGlobalValues(cell->self(), indices.size(), &values[0], &indices[0]);
  }
  diffusion_updown.FillComplete();
  diffusion_leftright.FillComplete();
  diffusion_ud_problem.SetOperator(&diffusion_updown);
  diffusion_lr_problem.SetOperator(&diffusion_leftright);

  Amesos Amesos_Factory;
  amesos_ud_solver = Amesos_Factory.Create("Amesos_Klu", diffusion_ud_problem);
  amesos_lr_solver = Amesos_Factory.Create("Amesos_Klu", diffusion_lr_problem);
  
  amesos_ud_solver->SymbolicFactorization();
  amesos_ud_solver->NumericFactorization();
  amesos_lr_solver->SymbolicFactorization();
  amesos_lr_solver->NumericFactorization();
}*/

void StokesSolver::assemble_diffusion_matrix()
{

  double eta_x = dt/grid.dx/grid.dx;
  double eta_y = dt/grid.dy/grid.dy;
  double alpha_x = -eta_x;
  double lambda_x = (1.0 + 2.0*eta_x);
  double alpha_y = -eta_y;
  double lambda_y = (1.0 + 2.0*eta_y);
  g.PutScalar(0.0);
  lux.PutScalar(0.0);
  luy.PutScalar(0.0);
  scratch1.PutScalar(0.0); //c prime
  scratch2.PutScalar(0.0); //d prime

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

void StokesSolver::assemble_stokes_matrix()
{
  std::vector<double> values;
  std::vector<int> indices; //5 point stencil

  //Assemble poisson matrix 
  for( StaggeredGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
  {
    StaggeredGrid::iterator up_cell(cell->up(), grid);
    StaggeredGrid::iterator down_cell(cell->down(), grid);

    indices.clear();
    values.clear();

    if(cell->at_top_boundary() || cell->at_bottom_boundary())
    {
      indices.push_back(cell->self());
      values.push_back(1.0);
    }
    else
    {
      indices.push_back(cell->self()); values.push_back(4.0/grid.dx/grid.dy);
      indices.push_back( cell->left() ); values.push_back(-1.0/grid.dx/grid.dy);
      indices.push_back( cell->right() ); values.push_back(-1.0/grid.dx/grid.dy);
      if(!up_cell->at_top_boundary())
        {indices.push_back( cell->up() ); values.push_back(-1.0/grid.dx/grid.dy);}
      if(!down_cell->at_bottom_boundary())
        {indices.push_back( cell->down() ); values.push_back(-1.0/grid.dx/grid.dy);}
    }
    int ierr =  poisson_matrix.InsertGlobalValues(cell->self(), indices.size(), &values[0], &indices[0]);
  }
  poisson_matrix.FillComplete();

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

  curl_T_spectral = new fftw_complex[ncells];

  //transform in the y direction with a dst
  n[0] = grid.ny; stride = grid.nx; dist = 1; howmany = grid.nx; 
  fftw_r2r_kind f[1] = {FFTW_RODFT00};
  dst = fftw_plan_many_r2r(1, n, howmany, curl_T.Values(), NULL, stride, dist,
                                     scratch1.Values(), NULL, stride, dist,
                                     f, FFTW_ESTIMATE);
  idst = fftw_plan_many_r2r(1, n, howmany, scratch1.Values(), NULL, stride, dist,
                                     scratch2.Values(), NULL, stride, dist,
                                     f, FFTW_ESTIMATE);

  //transform in the x direction with a full dft
  n[0] = grid.nx; stride = 1; dist = grid.nx; howmany = grid.ny; 
  dft = fftw_plan_many_dft_r2c(1, n, howmany, scratch1.Values(), NULL, stride, dist,
                     curl_T_spectral, NULL, stride, dist, FFTW_ESTIMATE);
  idft = fftw_plan_many_dft_c2r(1, n, howmany, curl_T_spectral, NULL, stride, dist,
                     scratch1.Values(), NULL, stride, dist, FFTW_ESTIMATE);
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
  Teuchos::TimeMonitor LocalTimer(*Stokes);
  
  assemble_curl_T_vector();


  fftw_execute(dst);
  fftw_execute(dft);

  for( StaggeredGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
  {
    curl_T_spectral[cell->self()][0]*=freqs[cell->self()];
    curl_T_spectral[cell->self()][1]*=freqs[cell->self()];
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

/*void StokesSolver::solve_stokes()
{  
  Teuchos::TimeMonitor LocalTimer(*Stokes);
  
  assemble_curl_T_vector();

  poisson_problem.SetLHS(&vorticity);
  poisson_problem.SetRHS(&curl_T);
  aztec_solver.SetProblem(poisson_problem);
//  aztec_solver.SetPrecOperator(preconditioner);
  aztec_solver.SetAztecOption(AZ_pre_calc, AZ_reuse);
  aztec_solver.Iterate(100, 1.e-1);
  poisson_problem.SetLHS(&stream);
  poisson_problem.SetRHS(&vorticity);
  aztec_solver.SetProblem(poisson_problem);
//  aztec_solver.SetPrecOperator(preconditioner);
  aztec_solver.SetAztecOption(AZ_pre_calc, AZ_reuse);
  aztec_solver.Iterate(100, 1.e-1);


  //Come up with the velocities
  for( StaggeredGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
  {
    if( cell->at_top_boundary() == false)
      u[cell->self()] = (stream[cell->up()] - stream[cell->self()])/grid.dy;
    else
      u[cell->self()] = (stream[cell->self()] - stream[cell->down()])/grid.dy;

    v[cell->self()] = -(stream[cell->right()] - stream[cell->self()])/grid.dx;
  }

} */


void StokesSolver::draw()
{
  Teuchos::TimeMonitor LocalTimer(*Draw);
  double DX = 2.0/grid.nx;
  double DY = 2.0/grid.ny;

  glClearColor(0.0, 0.0, 0.0, 0.0);
  glClear(GL_COLOR_BUFFER_BIT);
  glOrtho(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0);
  glBegin(GL_TRIANGLE_STRIP);
  for( StaggeredGrid::iterator cell = grid.begin(); !cell->at_top_boundary(); ++cell)
  {
    if (cell->at_left_boundary())
      glBegin(GL_TRIANGLE_STRIP);

    if( !cell->at_right_boundary() )
    {
      color c_s = hot(T[cell->self()]);
      color c_u = hot(T[cell->up()]);

      glColor3f(c_s.R, c_s.G, c_s.B);
      glVertex2f(cell->xindex()*DX-1.0, cell->yindex()*DY-1.0);
      glColor3f(c_u.R, c_u.G, c_u.B);
      glVertex2f((cell->xindex())*DX-1.0, (cell->yindex()+1)*DY-1.0);
    }
    else
      glEnd();
  }
  glFlush();
}
  
