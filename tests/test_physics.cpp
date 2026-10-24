/************************************
  Physics regression tests for the
  solver. Each test checks a property
  that optimizations must preserve.
  Run all tests, or name some:

    test_physics [conduction] [stokes] ...
************************************/
#include <cmath>
#include <cstdio>
#include <cstring>
#include <functional>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "../convection.h"
#include "../harness.h"

static int failures = 0;

#define CHECK(cond, ...)                                       \
  do                                                           \
  {                                                            \
    if (!(cond))                                               \
    {                                                          \
      std::printf("    FAIL %s:%d: %s\n      ", __FILE__, __LINE__, #cond); \
      std::printf(__VA_ARGS__);                                \
      std::printf("\n");                                       \
      ++failures;                                              \
    }                                                          \
  } while (0)

static bool close(double value, double expected, double rel_tol)
{
  return std::fabs(value - expected) <= rel_tol*std::fabs(expected);
}

const double r_inner = 0.5;
const int ntheta = 1024, nr = 128;

/* Below the critical Rayleigh number there is no convection, and the
   temperature relaxes to the conductive profile in an annulus,
   T(r) = ln(r/r_o)/ln(r_i/r_o), with a Nusselt number of one. */
static void test_conduction()
{
  ConvectionSimulator sim(r_inner, ntheta, nr, false);
  sim.update_state(100.);
  for (int s = 0; s < 100; ++s)
  {
    sim.solve_stokes();
    sim.semi_lagrangian_advect_temperature();
    sim.diffuse_temperature();
  }

  const RegularGrid &grid = sim.get_grid();
  double max_error = 0.;
  for (int j = 0; j < nr; ++j)
    for (int i = 0; i < ntheta; ++i)
    {
      const double r = grid.r_inner + j*grid.dr;
      const double exact = std::log(r/grid.r_outer)/std::log(grid.r_inner/grid.r_outer);
      max_error = std::max(max_error, std::fabs(sim.temperature()[j*ntheta + i] - exact));
    }
  const Diagnostics d = compute_diagnostics(sim);
  CHECK(max_error < 1.e-5, "max |T - T_conductive| = %g", max_error);
  CHECK(close(d.nusselt, 1.0, 0.01), "Nu = %.6f, expected 1", d.nusselt);
  CHECK(d.v_max < 1.e-12, "v_max = %g, expected no flow", d.v_max);
}

/* Solve a dense linear system with partial pivoting. Deliberately
   independent of the solver's Thomas algorithm. */
static std::vector<double> dense_solve(std::vector<std::vector<double>> A, std::vector<double> b)
{
  const int n = b.size();
  for (int k = 0; k < n; ++k)
  {
    int p = k;
    for (int i = k+1; i < n; ++i)
      if (std::fabs(A[i][k]) > std::fabs(A[p][k])) p = i;
    std::swap(A[k], A[p]);
    std::swap(b[k], b[p]);
    for (int i = k+1; i < n; ++i)
    {
      const double f = A[i][k]/A[k][k];
      for (int j = k; j < n; ++j) A[i][j] -= f*A[k][j];
      b[i] -= f*b[k];
    }
  }
  std::vector<double> x(n);
  for (int i = n-1; i >= 0; --i)
  {
    double sum = b[i];
    for (int j = i+1; j < n; ++j) sum -= A[i][j]*x[j];
    x[i] = sum/A[i][i];
  }
  return x;
}

/* For a temperature field that is a single Fourier mode in theta,
   T = g(r) sin(l theta), the discrete Stokes problem decouples: the
   stream function is p(r) cos(l theta), where p solves the radial
   biharmonic problem for wavenumber l. Compare the solver's stream
   function and velocities with an independent dense solve of the
   same discretization. This exercises the FFTs, mode indexing,
   normalization, tridiagonal solves, and velocity differencing. */
static void test_stokes_single_mode()
{
  ConvectionSimulator sim(r_inner, ntheta, nr, false);
  sim.update_state(1.e4);
  const RegularGrid &grid = sim.get_grid();
  const double Ra = sim.rayleigh_number();
  const int l = 3;
  const double dr = grid.dr, dtheta = grid.dtheta;

  std::vector<double> T(grid.ncells);
  std::vector<double> g(nr), r(nr);
  for (int j = 0; j < nr; ++j)
  {
    r[j] = grid.r_inner + j*dr;
    g[j] = std::sin(M_PI*j/(nr-1.));
    for (int i = 0; i < ntheta; ++i)
      T[j*ntheta + i] = g[j]*std::sin(l*i*dtheta);
  }
  sim.set_temperature(T.data());
  sim.solve_stokes();

  //Centered differences in theta turn d/dtheta of sin(l theta) into
  //s_l cos(l theta), with s_l = sin(l dtheta)/dtheta.
  const double s_l = std::sin(l*dtheta)/dtheta;
  std::vector<double> rhs(nr, 0.);
  for (int j = 1; j < nr-1; ++j)
    rhs[j] = Ra*g[j]*s_l/r[j];

  //Radial operator d^2/dr^2 + 1/r d/dr - l^2/r^2 with zero Dirichlet rows
  std::vector<std::vector<double>> A(nr, std::vector<double>(nr, 0.));
  A[0][0] = A[nr-1][nr-1] = 1.;
  for (int j = 1; j < nr-1; ++j)
  {
    A[j][j-1] = 1./dr/dr - 0.5/dr/r[j];
    A[j][j] = -2./dr/dr - double(l*l)/r[j]/r[j];
    A[j][j+1] = 1./dr/dr + 0.5/dr/r[j];
  }
  const std::vector<double> p = dense_solve(A, dense_solve(A, rhs));

  double p_max = 0.;
  for (double x : p) p_max = std::max(p_max, std::fabs(x));

  double psi_err = 0., u_err = 0., v_err = 0., u_max = 0., v_max = 0.;
  for (int j = 0; j < nr; ++j)
    for (int i = 0; i < ntheta; ++i)
    {
      const int k = j*ntheta + i;
      const double c = std::cos(l*i*dtheta), s = std::sin(l*i*dtheta);
      psi_err = std::max(psi_err, std::fabs(sim.stream_function()[k] - p[j]*c));

      double dpdr;
      if (j == 0) dpdr = (p[1] - p[0])/dr;
      else if (j == nr-1) dpdr = (p[nr-1] - p[nr-2])/dr;
      else dpdr = (p[j+1] - p[j-1])/2./dr;
      const double u = dpdr*c;
      const double v = p[j]*s_l*s/r[j];
      u_err = std::max(u_err, std::fabs(sim.velocity_theta()[k] - u));
      v_err = std::max(v_err, std::fabs(sim.velocity_r()[k] - v));
      u_max = std::max(u_max, std::fabs(u));
      v_max = std::max(v_max, std::fabs(v));
    }
  CHECK(p_max > 0., "trivial stream function");
  CHECK(psi_err < 1.e-9*p_max, "stream function error %g (scale %g)", psi_err, p_max);
  CHECK(u_err < 1.e-9*u_max, "u_theta error %g (scale %g)", u_err, u_max);
  CHECK(v_err < 1.e-9*v_max, "u_r error %g (scale %g)", v_err, v_max);
}

/* A source in a uniform-temperature mantle should reach a seismometer
   at the expected travel time, and the wave energy should decay. */
static void test_wave_arrival()
{
  ConvectionSimulator sim(r_inner, ntheta, nr, false);
  const RegularGrid &grid = sim.get_grid();
  //Must match propagate_seismic_waves()
  const double tstep = 0.6*std::min(grid.r_inner*grid.dtheta, grid.dr);
  const double speed = 1.0 - 0.7*sim.temperature()[0];
  const double r = 0.25;  //depth coordinate, mid-mantle

  for (double angle : {0.3, 0.6})
  {
    sim.clear_seismic_waves();
    sim.earthquake(0.0, r);
    sim.place_seismometer(angle, r);
    const double distance = 2.*(r + grid.r_inner)*std::sin(angle/2.);
    const double expected_step = distance/speed/tstep;

    double peak = 0., early_energy = 0., late_energy = 0.;
    int peak_step = -1;
    bool finite = true;
    for (int s = 0; s < 1200; ++s)
    {
      sim.propagate_seismic_waves();
      const double x = std::fabs(sim.seismometer_reading());
      if (x > peak) { peak = x; peak_step = s; }
      if (s == 300 || s == 1199)
      {
        double e = 0.;
        for (int k = 0; k < grid.ncells; ++k)
        {
          e += sim.displacement()[k]*sim.displacement()[k];
          finite = finite && std::isfinite(sim.displacement()[k]);
        }
        (s == 300 ? early_energy : late_energy) = e;
      }
    }
    CHECK(finite, "non-finite displacement");
    CHECK(close(peak_step, expected_step, 0.05),
          "angle %.1f: peak at step %d, expected %.1f", angle, peak_step, expected_step);
    CHECK(late_energy < 0.8*early_energy, "energy %g -> %g did not decay", early_energy, late_energy);
  }
}

/* Time-averaged Nusselt number and RMS velocity from the fixed heat
   sources, compared with values recorded from the original code
   (247ae78). At these Rayleigh numbers 1000 steps is not yet chaotic,
   so the tolerances allow roundoff-level changes (FFT library, loop
   order) but catch real changes to the physics. */
static void check_convection(double ra, double nu_ref, double vrms_ref, double rel_tol)
{
  ConvectionSimulator sim(r_inner, ntheta, nr, false);
  sim.update_state(ra);
  double nu = 0., vrms = 0., min_T = 0., max_T = 1.;
  bool finite = true;
  int n = 0;
  for (int s = 0; s < 1000; ++s)
  {
    advance(sim, s);
    if (s >= 500)
    {
      const Diagnostics d = compute_diagnostics(sim);
      nu += d.nusselt;
      vrms += d.v_rms;
      min_T = std::min(min_T, d.min_T);
      max_T = std::max(max_T, d.max_T);
      finite = finite && d.finite;
      ++n;
    }
  }
  nu /= n;
  vrms /= n;
  CHECK(finite, "Ra=%g: non-finite state", ra);
  CHECK(close(nu, nu_ref, rel_tol), "Ra=%g: mean Nu %.6f, reference %.6f", ra, nu, nu_ref);
  CHECK(close(vrms, vrms_ref, rel_tol), "Ra=%g: mean v_rms %.4f, reference %.4f", ra, vrms, vrms_ref);
  CHECK(min_T >= 0. && max_T <= 1., "Ra=%g: T outside [0,1]: [%g, %g]", ra, min_T, max_T);
}

static void test_convection_regression()
{
  check_convection(1.e5, 5.130927, 103.9836, 0.01);
  check_convection(1.e6, 10.494390, 453.0153, 0.03);
}

/* Results should not depend on the number of threads. This holds
   bitwise today, which makes refactors easy to verify exactly. If a
   future change (e.g. a parallel reduction) breaks it on purpose,
   relax this to a tolerance. */
static void test_thread_independence()
{
#ifdef _OPENMP
  const int max_threads = omp_get_max_threads();
  if (max_threads < 2)
  {
    std::printf("    (skipped: only one thread available)\n");
    return;
  }
  std::vector<std::vector<double>> results;
  for (int threads : {1, max_threads})
  {
    omp_set_num_threads(threads);
    ConvectionSimulator sim(r_inner, ntheta, nr, false);
    sim.update_state(1.e6);
    for (int s = 0; s < 100; ++s)
      advance(sim, s);
    results.emplace_back(sim.temperature(), sim.temperature() + sim.get_grid().ncells);
  }
  omp_set_num_threads(max_threads);
  CHECK(results[0] == results[1], "temperature differs between 1 and %d threads", max_threads);
#else
  std::printf("    (skipped: built without OpenMP)\n");
#endif
}

int main(int argc, char **argv)
{
  const struct { const char *name; std::function<void()> run; } tests[] = {
    {"conduction", test_conduction},
    {"stokes", test_stokes_single_mode},
    {"waves", test_wave_arrival},
    {"convection", test_convection_regression},
    {"threads", test_thread_independence},
  };

  int run = 0;
  for (const auto &t : tests)
  {
    bool selected = argc == 1;
    for (int i = 1; i < argc; ++i)
      selected = selected || std::strcmp(argv[i], t.name) == 0;
    if (!selected) continue;

    const int before = failures;
    const double start = now_ms();
    std::printf("%-12s ...\n", t.name);
    t.run();
    std::printf("%-12s %s (%.1f s)\n", t.name, failures == before ? "ok" : "FAILED",
                (now_ms() - start)/1000.);
    ++run;
  }
  if (run == 0)
  {
    std::printf("No tests matched\n");
    return 2;
  }
  std::printf("%d failure(s)\n", failures);
  return failures == 0 ? 0 : 1;
}
