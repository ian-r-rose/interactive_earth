#ifndef COLOR_H
#define COLOR_H

//Basic utilities for trying out colormaps

struct color
{
  double R, G, B;
};

inline color hot(double x)
{
  color c;
  double red_cutoff = 0.3;
  double green_cutoff = 0.7;
  c.R = (x < red_cutoff ?  (x/red_cutoff) : 1.0);
  c.G = (x < red_cutoff ? 0.0 : (x < green_cutoff ? (x-red_cutoff)/(green_cutoff-red_cutoff) : 1.0));
  c.B = (x < green_cutoff ? 0.0 : (x-green_cutoff)/(1.0-green_cutoff));
  return c;
}

#endif
