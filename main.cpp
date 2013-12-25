#include "stokes.h"


int main()
{

  StokesSolver stokes(1.0, 1.0, 100,100);
  for(int i = 0; i<20; ++i)
  {
    stokes.solve_stokes();
    for(int j=0; j<10; ++j)
      stokes.upwind_advect();
  }
  return 0;
}

