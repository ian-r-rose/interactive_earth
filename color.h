//Basic utilities for tryin out colormaps

struct color
{
  double R, G, B;
};

color hot(double x)
{
  color c;
  double red_cutoff = 0.4;
  double green_cutoff = 0.8;
  c.R = (x < red_cutoff ?  (x/red_cutoff) : 1.0);
  c.G = (x < red_cutoff ? 0.0 : (x < green_cutoff ? (x-red_cutoff)/(green_cutoff-red_cutoff) : 1.0));
  c.B = (x < green_cutoff ? 0.0 : (x-green_cutoff)/(1.0-green_cutoff));
  return c;
}

color hot_with_composition( double x, double y)
{
  color c;
  double red_cutoff = 0.4;
  double green_cutoff = 0.8;
  c.R = (x < red_cutoff ?  (x/red_cutoff) : 1.0);
  c.G = (x < red_cutoff ? 0.0 : (x < green_cutoff ? (x-red_cutoff)/(green_cutoff-red_cutoff) : 1.0));
  c.B = (x < green_cutoff ? 0.0 : (x-green_cutoff)/(1.0-green_cutoff));

  if(y > 0.8)
  {
    c.R = 0.2*c.R;
    c.G = 0.2*c.G;
    c.B = 0.8;
  }

  return c;
}
