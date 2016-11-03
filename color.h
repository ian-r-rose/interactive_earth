#ifndef COLOR_H
#define COLOR_H

#include <cmath>

//Basic utilities for trying out colormaps

struct color
{
  double R, G, B;
};



/*******************************
  Matplotlib colormaps thanks
  to Stefan van der Walt,
  Nathaniel J. Smith and Eric
  Firing.
*******************************/

color hot(double x);
color magma(double x);
color viridis(double x);
color plasma(double x);
color inferno(double x);

#endif
