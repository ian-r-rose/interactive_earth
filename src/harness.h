/************************************
  Shared pieces for driving the solver
  without user input: a deterministic
  set of heat sources, the timestep
  sequence used by the interactive
  app, and summary diagnostics. Used
  by the benchmark, the tests, and the
  app's screenshot mode.
************************************/
#ifndef HARNESS_H
#define HARNESS_H

#include <algorithm>
#include <cmath>
#include "convection.h"
#include "profiler.h"

//Number of initial steps for which the fixed heat sources are applied
const int seed_steps = 30;

//Four hot plumes near the base and two cold blobs near the top,
//standing in for user clicks.
inline void add_seed_heat_sources(ConvectionSimulator &sim)
{
  for (int k = 0; k < 4; ++k)
    sim.add_heat(0.3 + k * M_PI / 2., 0.05, true);
  for (int k = 0; k < 2; ++k)
    sim.add_heat(1.0 + k * M_PI, 0.45, false);
}

//One convection timestep, in the same order as the interactive app.
inline void advance(ConvectionSimulator &sim, int step_number)
{
  const double start = now_ms();
  sim.solve_stokes();
  if (step_number < seed_steps)
    add_seed_heat_sources(sim);
  sim.semi_lagrangian_advect_temperature();
  if (sim.has_composition())
    sim.semi_lagrangian_advect_composition();
  sim.diffuse_temperature();
  profiler().accumulate(Phase::step, now_ms() - start);
}

struct Diagnostics
{
  double nusselt;
  double mean_T, rms_T, min_T, max_T;
  double v_rms, v_max;
  bool finite;
};

inline Diagnostics compute_diagnostics(ConvectionSimulator &sim)
{
  const int n = sim.get_grid().ncells;
  const double *T = sim.temperature();
  const double *u = sim.velocity_theta();
  const double *v = sim.velocity_r();
  Diagnostics d = {sim.nusselt_number(), 0., 0., T[0], T[0], 0., 0., true};
  for (int i = 0; i < n; ++i)
  {
    const double vv = u[i] * u[i] + v[i] * v[i];
    d.mean_T += T[i];
    d.rms_T += T[i] * T[i];
    d.min_T = std::min(d.min_T, T[i]);
    d.max_T = std::max(d.max_T, T[i]);
    d.v_rms += vv;
    d.v_max = std::max(d.v_max, std::sqrt(vv));
    d.finite = d.finite && std::isfinite(T[i]) && std::isfinite(vv);
  }
  d.mean_T /= n;
  d.rms_T = std::sqrt(d.rms_T / n);
  d.v_rms = std::sqrt(d.v_rms / n);
  d.finite = d.finite && std::isfinite(d.nusselt);
  return d;
}

#endif
