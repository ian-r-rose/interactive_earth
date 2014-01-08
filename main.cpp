#include <GL/freeglut.h>
#include <GL/gl.h>
#include "stokes.h"

const unsigned int nx = 200;
const unsigned int ny = 100;
StokesSolver* handle = NULL;

void renderFunction( int ms)
{
  static int i=0;
  if(i%1==0)
    handle->draw();
  if(i%2== 0)
    handle->solve_stokes();
  handle->semi_lagrangian_advect();
  handle->diffuse_temperature();
  ++i;
  std::cout<<"Step "<<i<<std::endl;

  glutTimerFunc(ms, renderFunction, 0);
  glutPostRedisplay();

}

/* Main method - main entry point of application
the freeglut library does the window creation work for us, 
regardless of the platform. */
int main(int argc, char** argv)
{
    StokesSolver stokes(2.0, 1.0, nx,ny);
    handle = &stokes;

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_SINGLE);
    glutInitWindowSize(5*nx,5*ny);
    glutInitWindowPosition(10,10);
    glutCreateWindow("Convection");
//    glutDisplayFunc(renderFunction);
    glutTimerFunc(50, renderFunction, 0);
    glutMainLoop();    
    return 0;
}
