#include <AztecOO.h>
#include <Amesos.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>
#include <ml_MultiLevelPreconditioner.h>

#ifdef EPETRA_MPI
#  include <mpi.h>
#  include <Epetra_MpiComm.h>
#else
#  include <Epetra_SerialComm.h>
#endif

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
    double dt;

    StaggeredGrid grid;
    AztecOO aztec_solver;
    Amesos_BaseSolver* amesos_ud_solver;
    Amesos_BaseSolver* amesos_lr_solver;
    Epetra_LinearProblem diffusion_ud_problem, diffusion_lr_problem;
    Epetra_LinearProblem poisson_problem;
    
#ifdef EPETRA_MPI
    Epetra_MpiComm Comm;
#else
    Epetra_SerialComm Comm;
#endif

    Epetra_Map map;

    Epetra_Vector T;
    Epetra_Vector scratch1, scratch2;
    Epetra_Vector vorticity;
    Epetra_Vector stream;
    Epetra_Vector dTdx;
    Epetra_Vector u;
    Epetra_Vector v;

    Epetra_CrsMatrix poisson_matrix;
    Epetra_CrsMatrix diffusion_updown;
    Epetra_CrsMatrix diffusion_leftright;
    
    ML_Epetra::MultiLevelPreconditioner * MLPrec;

 
    //workhorse functions
    void initialize_temperature();
    void assemble_stokes_matrix();
    void assemble_diffusion_rhs();
    void assemble_diffusion_matrix();
    void assemble_dTdx_vector();
   
    //functions for evaluating field at points
    double initial_temperature(const Point&);
    double temperature(const Point&);
    Point velocity(const Point&);

  public:
    StokesSolver( double lx, double ly, int nx, int ny);
    void upwind_advect();
    void semi_lagrangian_advect();
    void diffuse_temperature();
    void solve_stokes();
    void draw();
};

#endif
