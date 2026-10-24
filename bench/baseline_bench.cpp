// Temporary baseline harness for the *unmodified* solver.
// Superseded by ie_bench once the physics/render split lands.
// Reads private solver state for diagnostics; do not imitate.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <iterator>
#include <string>
#include <vector>
#include "GL/glew.h"
#include "SDL2/SDL.h"
#include "SDL2/SDL_opengl.h"
#include "fftw3.h"

#define private public
#include "../convection.h"
#undef private

#ifdef _OPENMP
#include <omp.h>
#endif

// Globals referenced by rendering code are not linked here.
const double flattening = 0.0;

using clk = std::chrono::steady_clock;

struct Series
{
  std::vector<double> ms;
  void add(clk::time_point a, clk::time_point b)
  {
    ms.push_back(std::chrono::duration<double, std::milli>(b - a).count());
  }
  double pct(double p) const
  {
    std::vector<double> s(ms);
    std::sort(s.begin(), s.end());
    return s[std::min(s.size() - 1, size_t(p * (s.size() - 1) + 0.5))];
  }
};

int main(int argc, char **argv)
{
  int ntheta = 1024, nr = 128, steps = 1000, warmup = 50, seed_steps = 30;
  double ra = 1.e6;
  for (int i = 1; i + 1 < argc; i += 2)
  {
    std::string a = argv[i];
    if (a == "--ntheta") ntheta = atoi(argv[i + 1]);
    else if (a == "--nr") nr = atoi(argv[i + 1]);
    else if (a == "--steps") steps = atoi(argv[i + 1]);
    else if (a == "--warmup") warmup = atoi(argv[i + 1]);
    else if (a == "--ra") ra = atof(argv[i + 1]);
  }
  int threads = 1;
#ifdef _OPENMP
  threads = omp_get_max_threads();
#endif

  ConvectionSimulator sim(0.5, ntheta, nr, false);
  sim.update_state(ra);

  Series stokes, advect, diffuse, step, waves;
  std::vector<double> nu_series;
  for (int s = 0; s < warmup + steps; ++s)
  {
    auto t0 = clk::now();
    sim.solve_stokes();
    auto t1 = clk::now();
    // Deterministic "clicks": four hot plumes near the base, two cold
    // blobs near the top, applied for the first few steps.
    if (s < seed_steps)
    {
      for (int k = 0; k < 4; ++k)
        sim.add_heat(0.3 + k * M_PI / 2., 0.05, true);
      for (int k = 0; k < 2; ++k)
        sim.add_heat(1.0 + k * M_PI, 0.45, false);
    }
    auto t2 = clk::now();
    sim.semi_lagrangian_advect_temperature();
    auto t3 = clk::now();
    sim.diffuse_temperature();
    auto t4 = clk::now();
    if (s >= warmup)
    {
      stokes.add(t0, t1);
      advect.add(t2, t3);
      diffuse.add(t3, t4);
      step.add(t0, t4);
    }
    if (s % 50 == 0) nu_series.push_back(sim.nusselt_number());
  }

  // Seismic phase timing, from a fixed source.
  sim.clear_seismic_waves();
  sim.earthquake(1.0, 0.25);
  for (int s = 0; s < 200; ++s)
  {
    auto a = clk::now();
    sim.propagate_seismic_waves();
    waves.add(a, clk::now());
  }
  const double seis = sim.seismometer_reading();

  // Diagnostics.
  const int n = sim.grid.ncells;
  double sumT = 0, sumT2 = 0, v2 = 0, vmax = 0;
  bool finite = true;
  for (int i = 0; i < n; ++i)
  {
    sumT += sim.T[i];
    sumT2 += sim.T[i] * sim.T[i];
    const double vv = sim.u[i] * sim.u[i] + sim.v[i] * sim.v[i];
    v2 += vv;
    vmax = std::max(vmax, std::sqrt(vv));
    finite = finite && std::isfinite(sim.T[i]) && std::isfinite(vv);
  }

  printf("{\n  \"ntheta\": %d, \"nr\": %d, \"ra\": %g, \"ra_effective\": %g,\n", ntheta, nr, ra,
         sim.rayleigh_number());
  printf("  \"steps\": %d, \"warmup\": %d, \"threads\": %d,\n", steps, warmup, threads);
  auto ser = [](const char *name, const Series &x, bool comma) {
    printf("  \"%s_ms\": {\"median\": %.4f, \"p95\": %.4f}%s\n", name, x.pct(0.5), x.pct(0.95),
           comma ? "," : "");
  };
  ser("stokes", stokes, true);
  ser("advect", advect, true);
  ser("diffuse", diffuse, true);
  ser("waves", waves, true);
  ser("step", step, true);
  printf("  \"steps_per_s\": %.1f,\n", 1000. / step.pct(0.5));
  printf("  \"diag\": {\"nu\": %.10g, \"mean_T\": %.10g, \"rms_T\": %.10g, \"v_rms\": %.10g, "
         "\"v_max\": %.10g, \"seismometer\": %.10g, \"finite\": %s},\n",
         sim.nusselt_number(), sumT / n, std::sqrt(sumT2 / n), std::sqrt(v2 / n), vmax, seis,
         finite ? "true" : "false");
  printf("  \"nu_every_50\": [");
  for (size_t i = 0; i < nu_series.size(); ++i)
    printf("%s%.6g", i ? ", " : "", nu_series[i]);
  printf("]\n}\n");
  return finite ? 0 : 1;
}
