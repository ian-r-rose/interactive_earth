#include <Epetra_SerialComm.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>
#include "staggered_grid.h"

#ifndef STOKES_H
#define STOKES_H

class StokesSolver
{
  private:
    
    int nx,ny;
    int ncells;
    double Ra;

    StaggeredGrid grid;
    
    Epetra_SerialComm Comm;
    Epetra_Map map;

    Epetra_Vector T;
    Epetra_Vector vorticity;
    Epetra_Vector stream;
    Epetra_Vector dTdx;
    Epetra_Vector u;
    Epetra_Vector v;

    Epetra_CrsMatrix poisson_matrix;

 
    void initialize_temperature();
    double initial_temperature(const Point&);
    void assemble_stokes_matrix();
    void assemble_dTdx_vector();

  public:
    StokesSolver( double lx, double ly, int nx, int ny);
    void solve_stokes();
};

#endif
