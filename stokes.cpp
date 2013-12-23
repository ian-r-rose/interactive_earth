#include <Epetra_ConfigDefs.h>
#include <Epetra_SerialComm.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_Time.h>
#include <AztecOO.h>
#include <Amesos.h>
#include <ml_include.h>
#include <Epetra_LinearProblem.h>
#include <ml_MultiLevelOperator.h>
#include <ml_epetra_utils.h>

#include "staggered_grid.h"
#include <vector>
#include <cmath>
#include <fstream>

double temperature(const Point &p)
{
//  return -std::exp( (-(0.5-p.x)*(0.5-p.x) - (0.5-p.y)*(0.5-p.y) )/.05);
  return ( std::sqrt( (0.5-p.x)*(0.5-p.x)+(0.5-p.y)*(0.5-p.y)) < 0.15 ? 1.0 : 2.0);
}

int main(int argc, char **argv)
{
  Epetra_SerialComm Comm;

  int nx = 100, ny = 100;
  int ncells = nx*ny;
  double Ra = 10.e5; 

  Epetra_Map map(ncells, 0, Comm);

  Epetra_Vector T(map);
  Epetra_Vector vorticity(map);
  Epetra_Vector stream(map);
  Epetra_Vector dTdx(map);
  Epetra_Vector u(map);
  Epetra_Vector v(map);

  StaggeredGrid grid(1.0, 1.0, nx, ny);
  Epetra_CrsMatrix mat(Copy, map, 5);
  std::vector<int> indices; //5 point stencil
  std::vector<double> values;

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
      indices.push_back(cell->self()); values.push_back(-4.0);
      indices.push_back( cell->left() ); values.push_back(1.0);
      indices.push_back( cell->right() ); values.push_back(1.0);
      indices.push_back( cell->up() ); values.push_back(1.0);
      indices.push_back( cell->down() ); values.push_back(1.0);
    }

    int ierr =  mat.InsertGlobalValues(cell->self(), indices.size(), &values[0], &indices[0]);
  }
  mat.FillComplete();

  //Assemble T vector
  for( StaggeredGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
    T[cell->self()] = temperature(cell->center());

  //Assemble dTdx vector
  for( StaggeredGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
    if(cell->at_boundary())
      dTdx[cell->self()] = 0.0;
    else
      dTdx[cell->self()] = Ra*(T[cell->self()] - T[cell->left()]
                            + T[cell->down()] - T[cell->downleft()])/2.0/grid.dx;
 
  
  Epetra_LinearProblem poisson_problem;
  poisson_problem.SetOperator(&mat);

  Amesos Amesos_Factory;
  Amesos_BaseSolver *Solver;
  Solver = Amesos_Factory.Create("Amesos_Klu", poisson_problem);
  Solver->SymbolicFactorization(); 
  Solver->NumericFactorization();

  
  poisson_problem.SetLHS(&vorticity);
  poisson_problem.SetRHS(&dTdx);
  Solver->Solve();

  poisson_problem.SetLHS(&stream);
  poisson_problem.SetRHS(&vorticity);
  Solver->Solve();

  //Come up with the velocities
  for( StaggeredGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
  {
    u[cell->self()] = -(stream[cell->up()] - stream[cell->self()])/grid.dy;
    v[cell->self()] = (stream[cell->right()] - stream[cell->self()])/grid.dx;
  }

  std::ofstream temperature("temperature.txt");
  std::ofstream vort("vorticity.txt");
  std::ofstream psi("stream.txt");
  std::ofstream vel("velocity.txt");
  for( StaggeredGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
  {
    vort<<cell->corner().x<<"\t"<<cell->corner().y<<"\t"<<vorticity[cell->self()]<<std::endl;
    psi<<cell->corner().x<<"\t"<<cell->corner().y<<"\t"<<stream[cell->self()]<<std::endl;
    temperature<<cell->corner().x<<"\t"<<cell->corner().y<<"\t"<<T[cell->self()]<<std::endl;
    vel<<cell->corner().x<<"\t"<<cell->corner().y<<"\t"<<u[cell->self()]<<"\t"<<v[cell->self()]<<std::endl;
  }
} 
