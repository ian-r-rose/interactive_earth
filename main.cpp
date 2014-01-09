#include <GL/freeglut.h>
#include <GL/gl.h>
#include "stokes.h"

const unsigned int nx = 350;
const unsigned int ny = 175;
const double lx = 2.0;
const double ly = 1.0;
const unsigned int scale = 5;

const unsigned int xpix = nx*scale;
const unsigned int ypix = ny*scale;

int heat_state = 0;
double hx, hy;

StokesSolver* handle = NULL;

void keyboardFunction(unsigned char key, int x, int y)
{
  if(key == 27)
    exit(0);
}

void motionFunction( int x, int y)
{
  hx = lx*(double(x)/double(xpix));
  hy = ly*(1.0-double(y)/double(ypix));
}

void mouseFunction(int button, int state, int x, int y)
{
  if(state==GLUT_DOWN)
  {
     if(button == GLUT_LEFT_BUTTON)
       heat_state = 1;
     if(button == GLUT_RIGHT_BUTTON)
       heat_state = -1;
     hx = lx*(double(x)/double(xpix));
     hy = ly*(1.0-double(y)/double(ypix));
  }
  else heat_state=0;
}

void renderFunction(/* int ms*/)
{
  static int i=0;
  if(i%1==0)
    handle->draw();
  if(i%2== 0)
    handle->solve_stokes();

  if(heat_state != 0) handle->add_heat(hx, hy, (heat_state==1 ? true : false));
  handle->semi_lagrangian_advect();
  handle->diffuse_temperature();

  ++i;
  std::cout<<"Step "<<i<<std::endl;

//  glutTimerFunc(ms, renderFunction, 0);
  glutPostRedisplay();

}

/* Main method - main entry point of application
the freeglut library does the window creation work for us, 
regardless of the platform. */
int main(int argc, char** argv)
{
    StokesSolver stokes(lx, ly, nx,ny, 1.e7);
    handle = &stokes;

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_SINGLE);
    glutInitWindowSize(scale*nx,scale*ny);
    glutInitWindowPosition(10,10);
    glutCreateWindow("Convection");

    glutDisplayFunc(renderFunction);
    glutMotionFunc(motionFunction);
    glutMouseFunc(mouseFunction);
    glutKeyboardFunc(keyboardFunction); 

//    glutTimerFunc(50, renderFunction, 0);
    glutMainLoop();    
    return 0;
}
