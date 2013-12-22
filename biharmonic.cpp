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

double heat(const Point &p)
{
  return -std::exp( (-(0.5-p.x)*(0.5-p.x) - (0.5-p.y)*(0.5-p.y) )/.05);
}

int main(int argc, char **argv)
{
  Epetra_SerialComm Comm;

  int nx = 100, ny = 100;
  int ncells = nx*ny;

  Epetra_Map map(ncells, 0, Comm);

  Epetra_Vector rhs(map);
  Epetra_Vector T(map);

  StaggeredGrid grid(1.0, 1.0, nx, ny);
  Epetra_CrsMatrix mat(Copy, map, 5);
  std::vector<int> indices; //5 point stencil
  std::vector<double> values;
  for( StaggeredGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
  {
    indices.clear();
    values.clear();
    if(cell->at_boundary())
    {
      indices.push_back(cell->self());
      values.push_back(1.0);
      rhs[cell->self()] = 0.0;     
    }
    else
    {
      indices.push_back(cell->self()); values.push_back(-4.0);
      indices.push_back( cell->left() ); values.push_back(1.0);
      indices.push_back( cell->right() ); values.push_back(1.0);
      indices.push_back( cell->up() ); values.push_back(1.0);
      indices.push_back( cell->down() ); values.push_back(1.0);
      rhs[cell->self()] = heat(cell->corner());
    }

    int ierr =  mat.InsertGlobalValues(cell->self(), indices.size(), &values[0], &indices[0]);
  }
  mat.FillComplete();
 
  
  Epetra_LinearProblem laplace_problem;
  laplace_problem.SetOperator(&mat);
  laplace_problem.SetLHS(&T);
  laplace_problem.SetRHS(&rhs);

  Amesos Amesos_Factory;
  Amesos_BaseSolver *Solver;
  Solver = Amesos_Factory.Create("Amesos_Klu", laplace_problem);
  Solver->SymbolicFactorization(); 
  Solver->NumericFactorization();
  Solver->Solve();

  std::ofstream file("out.txt");
  for( StaggeredGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
    file<<cell->corner().x<<"\t"<<cell->corner().y<<"\t"<<T[cell->self()]<<std::endl;
  file.close(); 
} 
