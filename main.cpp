#include "GL/freeglut.h"
#include "GL/gl.h"
#include "stokes.h"
StokesSolver stokes(1.0, 1.0, 100,100);

/* display function - code from:
     http://fly.cc.fer.hr/~unreal/theredbook/chapter01.html
This is the actual usage of the OpenGL library. 
The following code is the same for any platform */
void renderFunction()
{
  static int i=0;
  if(i%10 == 0)
    stokes.solve_stokes();
  stokes.upwind_advect();
  stokes.draw();
  ++i;
  glutPostRedisplay();
}

/* Main method - main entry point of application
the freeglut library does the window creation work for us, 
regardless of the platform. */
int main(int argc, char** argv)
{
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_SINGLE);
    glutInitWindowSize(500,500);
    glutInitWindowPosition(100,100);
    glutCreateWindow("OpenGL - First window demo");
    glutDisplayFunc(renderFunction);
    glutMainLoop();    
    return 0;
}
