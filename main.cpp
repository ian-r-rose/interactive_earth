#include <GL/freeglut.h>
#include <GL/gl.h>
#include "stokes.h"
#include <Teuchos_TimeMonitor.hpp>

const unsigned int nx = 150;
const unsigned int ny = 75;
StokesSolver* handle = NULL;

void renderFunction()
{
  static int i=0;
  if(i%1==0)
    handle->draw();
  if(i%2== 0)
    handle->solve_stokes();
//  handle->upwind_advect();
  handle->semi_lagrangian_advect();
  handle->diffuse_temperature();
  ++i;
  std::cout<<"Step "<<i<<std::endl;
  glutPostRedisplay();
  if (i%100 == 0)
    Teuchos::TimeMonitor::summarize();
}

/* Main method - main entry point of application
the freeglut library does the window creation work for us, 
regardless of the platform. */
int main(int argc, char** argv)
{
#ifdef EPETRA_MPI
    MPI_Init(&argc,&argv);
#endif

    StokesSolver stokes(2.0, 1.0, nx,ny);
    handle = &stokes;

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_SINGLE);
    glutInitWindowSize(5*nx,5*ny);
    glutInitWindowPosition(10,10);
    glutCreateWindow("Convection");
    glutDisplayFunc(renderFunction);
    glutMainLoop();    
    return 0;
}
