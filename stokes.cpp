#include <Epetra_SerialComm.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>
#include <Amesos.h>
#include <AztecOO.h>
#include <Teuchos_TimeMonitor.hpp>

#include "stokes.h"
#include <vector>
#include <cmath>
#include <fstream>


Teuchos::RCP<Teuchos::Time> Diffusion = Teuchos::TimeMonitor::getNewCounter("Diffusion time");
Teuchos::RCP<Teuchos::Time> Advection = Teuchos::TimeMonitor::getNewCounter("Advection time");
Teuchos::RCP<Teuchos::Time> Stokes = Teuchos::TimeMonitor::getNewCounter("Stokes time");
Teuchos::RCP<Teuchos::Time> Draw = Teuchos::TimeMonitor::getNewCounter("Draw time");

StokesSolver::StokesSolver( double lx, double ly, int nx, int ny):
                          nx(nx), ny(ny), ncells(nx*ny), Ra(1.0e7),
                          grid(lx, ly, nx, ny), dt(5.e-7),
                          map(ncells, 0, Comm),
                          T(map), vorticity(map), stream(map), dTdx(map), 
                          u(map), v(map), scratch1(map), scratch2(map),
                          poisson_matrix(Copy, map, 5),
                          diffusion_updown(Copy, map, 3), diffusion_leftright(Copy,map,3)
{
  initialize_temperature();
  assemble_dTdx_vector();
  assemble_stokes_matrix();

  poisson_problem.SetOperator(&poisson_matrix);
  assemble_dTdx_vector();
  poisson_problem.SetLHS(&vorticity);
  poisson_problem.SetRHS(&dTdx);
  aztec_solver.SetProblem(poisson_problem);
  aztec_solver.SetAztecOption(AZ_solver, AZ_gmres);
//  aztec_solver.SetAztecOption(AZ_precond, AZ_Jacobi);
  aztec_solver.SetAztecOption(AZ_output, AZ_none);
  aztec_solver.SetAztecOption(AZ_keep_info, 1);
  aztec_solver.Iterate(100, 1.e-2);
  poisson_problem.SetLHS(&stream);
  poisson_problem.SetRHS(&vorticity);
  aztec_solver.SetProblem(poisson_problem);
  aztec_solver.Iterate(100, 1.e-2);

  assemble_diffusion_matrix();
}

double StokesSolver::initial_temperature(const Point &p)
{
//  return -std::exp( (-(0.5-p.x)*(0.5-p.x) - (0.5-p.y)*(0.5-p.y) )/.05);
//  return ( std::sqrt( (0.5-p.x)*(0.5-p.x)+(0.5-p.y)*(0.5-p.y)) < 0.15 ? 1.0 : 0.0);
//  double temp = 1.0-p.y/grid.ly + 0.1*std::sin( 2.0*3.14159*p.x/grid.lx )*(p.y)*(p.y-grid.ly)/grid.ly/grid.ly;
  return 0.5;
//  return (temp > 0.0 ? ( temp < 1.0 ? temp : 1.0 ) : 0.0);
}

void StokesSolver::initialize_temperature()
{
  for( StaggeredGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
    T[cell->self()] = initial_temperature(cell->center());
}

int StokesSolver::cell_id(const Point&p)
{
  int xindex = p.x/grid.dx;
  int yindex = p.y/grid.dy;
  return nx*yindex + xindex;
}

void StokesSolver::upwind_advect()
{
  Teuchos::TimeMonitor LocalTimer(*Advection);
 
  double dTdx;
  double dTdy;
  double vx, vy;
  scratch1.PutScalar(1.0);

  for( StaggeredGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
  {
    if( cell->at_right_boundary())
    {
      vx = 0;
      dTdx = 0;
    }
    else if (cell->at_left_boundary())
    {
      vx = 0;
      dTdx = 0;
    }
    else
    {
      vx = (u[cell->self()]+u[cell->right()])/2.0;
      if (vx <= 0.0) 
        dTdx = (T[cell->right()]-T[cell->self()])/grid.dx;
      else
        dTdx = (T[cell->self()]-T[cell->left()])/grid.dx;
    }
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
  
}
   
  
void StokesSolver::assemble_diffusion_matrix()
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
    else if(cell->at_left_boundary())
    {
      indices.push_back(cell->self()); values.push_back(1.0);
      indices.push_back(cell->right()); values.push_back(-1.0);
    } 
    else if(cell->at_right_boundary())
    {
      indices.push_back(cell->self()); values.push_back(1.0);
      indices.push_back(cell->left()); values.push_back(-1.0);
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
    else if(cell->at_left_boundary())
    {
      indices.push_back(cell->self()); values.push_back(1.0);
      indices.push_back(cell->right()); values.push_back(1.0);
    } 
    else if(cell->at_right_boundary())
    {
      indices.push_back(cell->self()); values.push_back(1.0);
      indices.push_back(cell->left()); values.push_back(1.0);
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

}

void StokesSolver::assemble_stokes_matrix()
{
  std::vector<double> values;
  std::vector<int> indices; //5 point stencil

  //Assemble poisson matrix 
  for( StaggeredGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
  {
    indices.clear();
    values.clear();
    if(cell->at_boundary())
    {
      indices.push_back(cell->self());
      values.push_back(1.0);
    }
    else
    {
      indices.push_back(cell->self()); values.push_back(-4.0/grid.dx/grid.dy);
      indices.push_back( cell->left() ); values.push_back(1.0/grid.dx/grid.dy);
      indices.push_back( cell->right() ); values.push_back(1.0/grid.dx/grid.dy);
      indices.push_back( cell->up() ); values.push_back(1.0/grid.dx/grid.dy);
      indices.push_back( cell->down() ); values.push_back(1.0/grid.dx/grid.dy);
    }

    int ierr =  poisson_matrix.InsertGlobalValues(cell->self(), indices.size(), &values[0], &indices[0]);
  }
  poisson_matrix.FillComplete();
//  new ML_Epetra::MultiLevelPreconditioner(poisson_matrix, true);
//  if(!MLPrec) std::cout<<"ERROR"<<std::endl;
//  aztec_solver.SetPrecOperator(MLPrec);

}

  

void StokesSolver::assemble_dTdx_vector()
{
  //Assemble dTdx vector
  for( StaggeredGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
    if(cell->at_boundary())
      dTdx[cell->self()] = 0.0;
    else
      dTdx[cell->self()] = Ra*(T[cell->self()] - T[cell->left()]
                            + T[cell->down()] - T[cell->downleft()])/2.0/grid.dx;
}

 
void StokesSolver::solve_stokes()
{  
  Teuchos::TimeMonitor LocalTimer(*Stokes);
  
  assemble_dTdx_vector();
  poisson_problem.SetLHS(&vorticity);
  poisson_problem.SetRHS(&dTdx);
  aztec_solver.SetProblem(poisson_problem);
  aztec_solver.SetAztecOption(AZ_pre_calc, AZ_reuse);
  aztec_solver.Iterate(100, 1.e-1);
  poisson_problem.SetLHS(&stream);
  poisson_problem.SetRHS(&vorticity);
  aztec_solver.SetProblem(poisson_problem);
  aztec_solver.SetAztecOption(AZ_pre_calc, AZ_reuse);
  aztec_solver.Iterate(100, 1.e-1);

  //Come up with the velocities
  for( StaggeredGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
  {
    if( cell->at_top_boundary() == false)
      u[cell->self()] = (stream[cell->up()] - stream[cell->self()])/grid.dy;
    else
      u[cell->self()] = (stream[cell->self()] - stream[cell->down()])/grid.dy;
    if ( cell->at_right_boundary() == false)
      v[cell->self()] = -(stream[cell->right()] - stream[cell->self()])/grid.dx;
    else 
      v[cell->self()] = -(stream[cell->self()] - stream[cell->left()])/grid.dx;
  }

} 


void StokesSolver::draw()
{
  Teuchos::TimeMonitor LocalTimer(*Draw);
  double DX = 2.0/grid.nx;
  double DY = 2.0/grid.ny;

  glClearColor(0.0, 0.0, 0.0, 0.0);
  glClear(GL_COLOR_BUFFER_BIT);
  glOrtho(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0);
  glBegin(GL_TRIANGLES);
  for( StaggeredGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
    if( cell->at_right_boundary() == false && cell->at_top_boundary() == false)
    {
      glColor3f(T[cell->self()], .5, 1.0-T[cell->self()]);
      glVertex2f(cell->xindex()*DX-1.0, cell->yindex()*DY-1.0);
      glColor3f(T[cell->upright()], .5, 1.0-T[cell->upright()]);
      glVertex2f((cell->xindex()+1)*DX-1.0, (cell->yindex()+1)*DY-1.0);
      glColor3f(T[cell->up()], .5, 1.0-T[cell->up()]);
      glVertex2f((cell->xindex())*DX-1.0, (cell->yindex()+1)*DY-1.0);

      glColor3f(T[cell->self()], .5, 1.0-T[cell->self()]);
      glVertex2f(cell->xindex()*DX-1.0, cell->yindex()*DY-1.0);
      glColor3f(T[cell->right()], .5, 1.0-T[cell->right()]);
      glVertex2f((cell->xindex()+1)*DX-1.0, cell->yindex()*DY-1.0);
      glColor3f(T[cell->upright()], .5, 1.0-T[cell->upright()]);
      glVertex2f((cell->xindex()+1)*DX-1.0, (cell->yindex()+1)*DY-1.0);

    }
  glEnd();
  glFlush();
}
  
