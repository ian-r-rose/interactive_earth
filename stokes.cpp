#include <Epetra_SerialComm.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>
#include <Amesos.h>
#include <AztecOO.h>

#include "stokes.h"
#include <vector>
#include <cmath>
#include <fstream>


StokesSolver::StokesSolver( double lx, double ly, int nx, int ny):
                          nx(nx), ny(ny), ncells(nx*ny), Ra(10.0),
                          grid(lx, ly, nx, ny),
                          map(ncells, 0, Comm),
                          T(map), vorticity(map), stream(map),
                          dTdx(map), u(map), v(map), Tnew(map),
                          poisson_matrix(Copy, map, 5)
{
  initialize_temperature();
  assemble_stokes_matrix();

  Amesos Amesos_Factory;
  poisson_problem.SetOperator(&poisson_matrix);
  poisson_solver = Amesos_Factory.Create("Amesos_Klu", poisson_problem);
  poisson_solver->SymbolicFactorization(); 
  poisson_solver->NumericFactorization();
}

double StokesSolver::initial_temperature(const Point &p)
{
//  return -std::exp( (-(0.5-p.x)*(0.5-p.x) - (0.5-p.y)*(0.5-p.y) )/.05);
  return ( std::sqrt( (0.5-p.x)*(0.5-p.x)+(0.5-p.y)*(0.5-p.y)) < 0.15 ? 0.0 : 1.0);
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
  double dTdx;
  double dTdy;
  double vx, vy;

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
    
    Tnew[cell->self()] = T[cell->self()] - .1*(vx*dTdx+vy*dTdy);
  }
  T = Tnew;
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

}

void StokesSolver::assemble_diffusion_rhs()
{
  for( StaggeredGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
    if(cell->at_boundary())
      dTdx[cell->self()] = 0.0;
    else
      dTdx[cell->self()] = Ra*(T[cell->self()] - T[cell->left()]
                            + T[cell->down()] - T[cell->downleft()])/2.0/grid.dx;
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
  
  assemble_dTdx_vector();
  poisson_problem.SetLHS(&vorticity);
  poisson_problem.SetRHS(&dTdx);
  poisson_solver->Solve();

  poisson_problem.SetLHS(&stream);
  poisson_problem.SetRHS(&vorticity);
  poisson_solver->Solve();

  //Come up with the velocities
  for( StaggeredGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
  {
    if( cell->at_top_boundary() == false)
      u[cell->self()] = -(stream[cell->up()] - stream[cell->self()])/grid.dy;
    else
      u[cell->self()] = -(stream[cell->self()] - stream[cell->down()])/grid.dy;
    if ( cell->at_right_boundary() == false)
      v[cell->self()] = (stream[cell->right()] - stream[cell->self()])/grid.dx;
    else 
      v[cell->self()] = (stream[cell->self()] - stream[cell->left()])/grid.dx;
  }

  std::ofstream temperature("temperature.txt");
  std::ofstream vort("vorticity.txt");
  std::ofstream psi("stream.txt");
  std::ofstream vel("velocity.txt");
  for( StaggeredGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
  {
    vort<<cell->corner().x<<"\t"<<cell->corner().y<<"\t"<<vorticity[cell->self()]<<std::endl;
    psi<<cell->corner().x<<"\t"<<cell->corner().y<<"\t"<<stream[cell->self()]<<std::endl;
    temperature<<cell->center().x<<"\t"<<cell->center().y<<"\t"<<T[cell->self()]<<std::endl;
    vel<<cell->corner().x<<"\t"<<cell->corner().y<<"\t"<<u[cell->self()]<<"\t"<<v[cell->self()]<<std::endl;
  }
} 


void StokesSolver::draw()
{
  double DX = 2.0/grid.nx;
  double DY = 2.0/grid.ny;

  glClearColor(0.0, 0.0, 0.0, 0.0);
  glClear(GL_COLOR_BUFFER_BIT);
  glOrtho(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0);
  glBegin(GL_TRIANGLES);
  for( StaggeredGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
    if( cell->at_right_boundary() == false && cell->at_top_boundary() == false)
    {
      glVertex2f(cell->xindex()*DX-1.0, cell->yindex()*DY-1.0);
      glColor3f(T[cell->self()], 0.0, 1.0-T[cell->self()]);
      glVertex2f((cell->xindex()+1)*DX-1.0, (cell->yindex()+1)*DY-1.0);
      glColor3f(T[cell->upright()], 0.0, 1.0-T[cell->upright()]);
      glVertex2f((cell->xindex())*DX-1.0, (cell->yindex()+1)*DY-1.0);
      glColor3f(T[cell->up()], 0.0, 1.0-T[cell->up()]);

      glVertex2f(cell->xindex()*DX-1.0, cell->yindex()*DY-1.0);
      glColor3f(T[cell->self()], 0.0, 1.0-T[cell->self()]);
      glVertex2f((cell->xindex()+1)*DX-1.0, cell->yindex()*DY-1.0);
      glColor3f(T[cell->right()], 0.0, 1.0-T[cell->right()]);
      glVertex2f((cell->xindex()+1)*DX-1.0, (cell->yindex()+1)*DY-1.0);
      glColor3f(T[cell->upright()], 0.0, 1.0-T[cell->upright()]);

    }
  glEnd();
  glFlush();
}
  
