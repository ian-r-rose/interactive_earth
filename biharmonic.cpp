#include <Epetra_ConfigDefs.h>
#include <Epetra_SerialComm.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_Time.h>
#include <AztecOO.h>
#include <ml_include.h>
#include <Epetra_LinearProblem.h>
#include <ml_MultiLevelOperator.h>
#include <ml_epetra_utils.h>

#include "staggered_grid.h"
#include <vector>


int main(int argc, char **argv)
{
  Epetra_SerialComm Comm;

  int nx = 10, ny = 5;
  int ncells = nx*ny;

  Epetra_Map map(ncells, 0, Comm);

  Epetra_Vector stream(map);
  Epetra_Vector vorticity(map);
  Epetra_Vector u(map);
  Epetra_Vector v(map);
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
    std::cout<<ierr<<std::endl;
  }
  mat.FillComplete();

  
 
} 
