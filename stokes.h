#include <Epetra_SerialComm.h>
#include <Amesos.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>
#include "staggered_grid.h"
#include <GL/gl.h>

#ifndef STOKES_H
#define STOKES_H

class StokesSolver
{
  private:
    
    int nx,ny;
    int ncells;
    double Ra;

    StaggeredGrid grid;
    Amesos_BaseSolver *poisson_solver;
    Epetra_LinearProblem poisson_problem;
    
    Epetra_SerialComm Comm;
    Epetra_Map map;

    Epetra_Vector T;
    Epetra_Vector Tnew;
    Epetra_Vector vorticity;
    Epetra_Vector stream;
    Epetra_Vector dTdx;
    Epetra_Vector u;
    Epetra_Vector v;

    Epetra_CrsMatrix poisson_matrix;

 
    //workhorse functions
    void initialize_temperature();
    void assemble_stokes_matrix();
    void assemble_diffusion_rhs();
    void assemble_dTdx_vector();
   
    //functions for evaluating field at points
    double initial_temperature(const Point&);
    double temperature(const Point&);
    Point velocity(const Point&);
    int cell_id(const Point&);

  public:
    StokesSolver( double lx, double ly, int nx, int ny);
    void upwind_advect();
    void diffuse_temperature();
    void solve_stokes();
    void draw();
};

#endif
