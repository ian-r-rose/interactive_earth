/************************************
  Headless benchmark for the solver.
  Runs the same timestep sequence as
  the interactive app from a fixed
  set of heat sources, and reports
  per-phase timings and diagnostics.

  ie_bench [--ntheta N] [--nr N] [--ra RA] [--steps N] [--warmup N]
           [--wave-steps N] [--threads N] [--composition] [--json]
           [--dump out.ppm]
************************************/
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "color.h"
#include "convection.h"
#include "harness.h"
#include "profiler.h"

static void usage()
{
  std::fprintf(stderr, "usage: ie_bench [--ntheta N] [--nr N] [--ra RA] [--steps N] [--warmup N]\n"
                       "                [--wave-steps N] [--threads N] [--composition] [--json]\n"
                       "                [--dump out.ppm]\n");
  std::exit(2);
}

//Render the temperature field as an annulus with the "hot" colormap,
//nearest-neighbor sampled, and write it as a binary PPM.
static bool dump_temperature(const ConvectionSimulator &sim, const char *path, int size = 768)
{
  const RegularGrid &grid = sim.get_grid();
  const double *T = sim.temperature();
  std::vector<unsigned char> img(3 * size * size, 0);
  for (int py = 0; py < size; ++py)
    for (int px = 0; px < size; ++px)
    {
      //Pixel center in [-1,1], y up
      const double x = 2. * (px + 0.5) / size - 1.;
      const double y = 1. - 2. * (py + 0.5) / size;
      const double r = std::sqrt(x * x + y * y);
      color c;
      if (r > grid.r_outer)
        continue;
      else if (r < grid.r_inner)
        c = hot(1.0);
      else
      {
        double theta = std::atan2(y, x);
        if (theta < 0.)
          theta += 2. * M_PI;
        const int i = int(theta / grid.dtheta + 0.5) % grid.ntheta;
        const int j = std::min(grid.nr - 1, int((r - grid.r_inner) / grid.dr + 0.5));
        c = hot(T[j * grid.ntheta + i]);
      }
      unsigned char *p = &img[3 * (py * size + px)];
      p[0] = (unsigned char)(255. * std::min(1., std::max(0., c.R)) + 0.5);
      p[1] = (unsigned char)(255. * std::min(1., std::max(0., c.G)) + 0.5);
      p[2] = (unsigned char)(255. * std::min(1., std::max(0., c.B)) + 0.5);
    }
  FILE *f = std::fopen(path, "wb");
  if (!f)
    return false;
  std::fprintf(f, "P6\n%d %d\n255\n", size, size);
  std::fwrite(img.data(), 1, img.size(), f);
  std::fclose(f);
  return true;
}

int main(int argc, char **argv)
{
  int ntheta = 1024, nr = 128, steps = 1000, warmup = 50, wave_steps = 200;
  double ra = 1.e6;
  bool composition = false, json = false;
  const char *dump_path = nullptr;

  for (int i = 1; i < argc; ++i)
  {
    const std::string a = argv[i];
    const bool has_value = i + 1 < argc;
    if (a == "--composition")
      composition = true;
    else if (a == "--json")
      json = true;
    else if (!has_value)
      usage();
    else if (a == "--ntheta")
      ntheta = std::atoi(argv[++i]);
    else if (a == "--nr")
      nr = std::atoi(argv[++i]);
    else if (a == "--ra")
      ra = std::atof(argv[++i]);
    else if (a == "--steps")
      steps = std::atoi(argv[++i]);
    else if (a == "--warmup")
      warmup = std::atoi(argv[++i]);
    else if (a == "--wave-steps")
      wave_steps = std::atoi(argv[++i]);
    else if (a == "--dump")
      dump_path = argv[++i];
    else if (a == "--threads")
    {
      const int t = std::atoi(argv[++i]);
#ifdef _OPENMP
      omp_set_num_threads(t);
#else
      if (t != 1)
        std::fprintf(stderr, "Built without OpenMP; ignoring --threads\n");
#endif
    }
    else
      usage();
  }
  if (steps < 1)
    usage();

  int threads = 1;
#ifdef _OPENMP
  threads = omp_get_max_threads();
#endif

  //Keep every sample, so the statistics cover the whole run
  profiler().set_window(0);

  ConvectionSimulator sim(0.5, ntheta, nr, composition);
  sim.update_state(ra);

  for (int s = 0; s < warmup + steps; ++s)
  {
    advance(sim, s);
    if (s < warmup)
      profiler().reset();
    else
      profiler().commit();
  }
  const Diagnostics d = compute_diagnostics(sim);

  //Seismic waves from a fixed source, timed separately
  if (wave_steps > 0)
  {
    sim.clear_seismic_waves();
    sim.earthquake(1.0, 0.25);
    for (int s = 0; s < wave_steps; ++s)
    {
      sim.propagate_seismic_waves();
      profiler().commit();
    }
  }

  if (dump_path && !dump_temperature(sim, dump_path))
    std::fprintf(stderr, "Could not write %s\n", dump_path);

  const Profiler &prof = profiler();
  const double step_ms = prof.percentile(Phase::step, 0.5);
  if (json)
  {
    std::printf("{\"ntheta\": %d, \"nr\": %d, \"ra\": %g, \"ra_effective\": %.10g, "
                "\"composition\": %s, \"steps\": %d, \"warmup\": %d, \"threads\": %d,\n",
                ntheta, nr, ra, sim.rayleigh_number(), composition ? "true" : "false", steps,
                warmup, threads);
    std::printf(" \"steps_per_s\": %.2f,\n \"phases_ms\": %s,\n", 1000. / step_ms,
                prof.json().c_str());
    std::printf(" \"diag\": {\"nu\": %.10g, \"mean_T\": %.10g, \"rms_T\": %.10g, \"min_T\": %.10g, "
                "\"max_T\": %.10g, \"v_rms\": %.10g, \"v_max\": %.10g, \"finite\": %s}}\n",
                d.nusselt, d.mean_T, d.rms_T, d.min_T, d.max_T, d.v_rms, d.v_max,
                d.finite ? "true" : "false");
  }
  else
  {
    std::printf("%dx%d  Ra=%.3g  threads=%d  steps=%d (+%d warmup)\n", ntheta, nr,
                sim.rayleigh_number(), threads, steps, warmup);
    std::printf("%-16s %9s %9s %9s %7s\n", "phase", "median", "p95", "mean", "%step");
    for (int i = 0; i < static_cast<int>(Phase::count); ++i)
    {
      const Phase p = static_cast<Phase>(i);
      if (!prof.has(p))
        continue;
      std::printf("%-16s %9.3f %9.3f %9.3f %6.1f%%\n", phase_name(p), prof.percentile(p, 0.5),
                  prof.percentile(p, 0.95), prof.mean(p),
                  p == Phase::waves ? 0. : 100. * prof.mean(p) / prof.mean(Phase::step));
    }
    std::printf("%.1f steps/s\n", 1000. / step_ms);
    std::printf("Nu=%.6g  T mean=%.4g rms=%.4g [%.3g, %.3g]  v rms=%.4g max=%.4g  %s\n", d.nusselt,
                d.mean_T, d.rms_T, d.min_T, d.max_T, d.v_rms, d.v_max,
                d.finite ? "finite" : "NON-FINITE");
  }
  return d.finite ? 0 : 1;
}
