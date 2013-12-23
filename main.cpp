#include "stokes.h"


int main()
{

  StokesSolver stokes(1.0, 1.0, 100,100);
  stokes.solve_stokes();
  return 0;
}

