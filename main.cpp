#include <GL/freeglut.h>
#include <GL/gl.h>
#include "stokes.h"
#include <Teuchos_TimeMonitor.hpp>

const unsigned int nx = 200;
const unsigned int ny = 100;
StokesSolver stokes(2.0, 1.0, nx,ny);

/* display function - code from:
     http://fly.cc.fer.hr/~unreal/theredbook/chapter01.html
This is the actual usage of the OpenGL library. 
The following code is the same for any platform */
void renderFunction()
{
  static int i=0;
  if(i%10== 0)
    stokes.solve_stokes();
  stokes.upwind_advect();
  stokes.diffuse_temperature();
  if(i%10==0)
    stokes.draw();
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
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_SINGLE);
    glutInitWindowSize(5*nx,5*ny);
    glutInitWindowPosition(50,50);
    glutCreateWindow("OpenGL - First window demo");
    glutDisplayFunc(renderFunction);
    glutMainLoop();    
    return 0;
}
