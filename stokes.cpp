#include <Epetra_SerialComm.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>
#include <Amesos.h>

#include "stokes.h"
#include <vector>
#include <cmath>
#include <fstream>


StokesSolver::StokesSolver( double lx, double ly, int nx, int ny):
                          nx(nx), ny(ny), ncells(nx*ny), Ra(10.0),
                          grid(lx, ly, nx, ny),
                          map(ncells, 0, Comm),
                          T(map), vorticity(map), stream(map),
                          dTdx(map), u(map), v(map),
                          poisson_matrix(Copy, map, 5)
{
  initialize_temperature();
  assemble_stokes_matrix();
}

double StokesSolver::initial_temperature(const Point &p)
{
//  return -std::exp( (-(0.5-p.x)*(0.5-p.x) - (0.5-p.y)*(0.5-p.y) )/.05);
  return ( std::sqrt( (0.5-p.x)*(0.5-p.x)+(0.5-p.y)*(0.5-p.y)) < 0.15 ? 1.0 : 2.0);
}

void StokesSolver::initialize_temperature()
{
  for( StaggeredGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
    T[cell->self()] = initial_temperature(cell->center());
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
  Epetra_LinearProblem poisson_problem;
  poisson_problem.SetOperator(&poisson_matrix);

  Amesos Amesos_Factory;
  Amesos_BaseSolver *Solver;
  Solver = Amesos_Factory.Create("Amesos_Klu", poisson_problem);
  Solver->SymbolicFactorization(); 
  Solver->NumericFactorization();

  
  assemble_dTdx_vector();
  poisson_problem.SetLHS(&vorticity);
  poisson_problem.SetRHS(&dTdx);
  Solver->Solve();

  poisson_problem.SetLHS(&stream);
  poisson_problem.SetRHS(&vorticity);
  Solver->Solve();

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


