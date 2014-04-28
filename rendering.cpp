#include "stokes.h"
#include "color.h"


void StokesSolver::draw()
{
  double DX = 2.0/grid.nx;
  double DY = 2.0/grid.ny;

  glClearColor(0.0, 0.0, 0.0, 0.0);
  glClear(GL_COLOR_BUFFER_BIT);
  glOrtho(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0);
  glBegin(GL_TRIANGLE_STRIP);
  for( StaggeredGrid::iterator cell = grid.begin(); !cell->at_top_boundary(); ++cell)
  {
    if (cell->at_left_boundary())
      glBegin(GL_TRIANGLE_STRIP);

    if( !cell->at_right_boundary() )
    {
      color c_s = hot(T[cell->self()]);
      color c_u = hot(T[cell->up()]);

      glColor3f(c_s.R, c_s.G, c_s.B);
      glVertex2f(cell->xindex()*DX-1.0, cell->yindex()*DY-1.0);
      glColor3f(c_u.R, c_u.G, c_u.B);
      glVertex2f((cell->xindex())*DX-1.0, (cell->yindex()+1)*DY-1.0);
    }
    else
      glEnd();
  }
  glFlush();
}
  
