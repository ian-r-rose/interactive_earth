#include <algorithm>
#include <cmath>
#include <vector>
#include "convection.h"
#include "profiler.h"

//standard math library functions can be unpredictibly slow in
//some implementations, it seems, and not even available in others.
//Here I implement some super basic, unsafe versions of some
//for inlining.
inline double fast_fmod(double x, double y)
{
  return x - ((int)(x / y)) * y;
}
inline double dmin(double x, double y)
{
  return x < y ? x : y;
}
inline double dmax(double x, double y)
{
  return x > y ? x : y;
}

/*Constructor for the solver. The parameters are in order:
  inner_radius : inner radius, compared to an outer radius of one.
  ntheta : number of cells in theta direction
  nr : number of cells in r direction
  Rayleigh : initial Rayleigh number
*/
ConvectionSimulator::ConvectionSimulator(double inner_radius, int ntheta, int nr,
                                         bool include_composition)
    : include_composition(include_composition), grid(inner_radius, ntheta, nr),
      stokes_matrices(grid.nr, n_modes()), diffusion_matrices(grid.nr, n_modes())
{
  //Sines and cosines of the grid angles, for advection
  grid_cos_theta.resize(grid.ntheta);
  grid_sin_theta.resize(grid.ntheta);
  for (int i = 0; i < grid.ntheta; ++i)
  {
    grid_cos_theta[i] = std::cos(i * grid.dtheta);
    grid_sin_theta[i] = std::sin(i * grid.dtheta);
  }

  //Allocate memory for data vectors
  T = fftw_alloc_real(grid.ncells);
  D = new double[grid.ncells];
  Dp = new double[grid.ncells];
  stream = fftw_alloc_real(grid.ncells);
  curl_density = fftw_alloc_real(grid.ncells);
  u = new double[grid.ncells];
  v = new double[grid.ncells];
  scratch = fftw_alloc_real(grid.ncells);

  if (include_composition)
    C = new double[grid.ncells];

  //I actually allocate more memory than necessary for the
  //spectral vector, but this way it makes the indexing simpler
  curl_density_spectral = reinterpret_cast<std::complex<double> *>(fftw_alloc_complex(grid.ncells));
  T_spectral = reinterpret_cast<std::complex<double> *>(fftw_alloc_complex(grid.ncells));

  //Plan the transforms before initializing anything, since planning may
  //overwrite the arrays it is given.
  setup_fft_plans();

  //Initialize the state
  buoyancy_number = 1.0;
  update_state(1.e99);
  initialize_temperature();
  if (include_composition)
    initialize_composition();
  spin = M_PI / 2.;
  seismometer_location.theta = M_PI / 2.;
  seismometer_location.r = 0.5;

  //Initialize displacement to zero
  clear_seismic_waves();

  //The Stokes matrices do not depend on the state; the diffusion
  //matrices are set up by update_state.
  setup_stokes_problem();
}

/* Destructor for the solver.*/
ConvectionSimulator::~ConvectionSimulator()
{
  fftw_free(T);
  delete[] D;
  delete[] Dp;
  fftw_free(stream);
  fftw_free(curl_density);
  delete[] u;
  delete[] v;
  fftw_free(scratch);
  fftw_free(curl_density_spectral);
  fftw_free(T_spectral);

  if (include_composition)
    delete[] C;

  fftw_destroy_plan(row_dft);
  fftw_destroy_plan(row_idft);
}

/* Functional form of initial temperature field.
   Just start with a constant value based on a guess of the steady state average.*/
double ConvectionSimulator::initial_temperature(const Point & /*p*/)
{
  //Guess of internal temperature based on isoviscous boundary layer theory
  //in a cylinder.
  return grid.r_inner / (1. + grid.r_inner);
}

/* Functional form of initial composition.
   Just start with zero.*/
double ConvectionSimulator::initial_composition(const Point & /*p*/)
{
  return 0.0;
}

/* Loop over all the cells and set the initial temperature */
void ConvectionSimulator::initialize_temperature()
{
  for (RegularGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
    T[cell->self()] = initial_temperature(cell->location());
}

/* Loop over all the cells and set the initial composition */
void ConvectionSimulator::initialize_composition()
{
  for (RegularGrid::iterator cell = grid.begin(); cell != grid.end(); ++cell)
    C[cell->self()] = initial_composition(cell->location());
}

/* Call f(k, rsq) for every cell k, with rsq its squared distance from
   the point (theta, r), split between threads. */
template <typename F>
void ConvectionSimulator::for_each_cell_distance(double theta, double r, F f) const
{
  const double x1 = (r + grid.r_inner) * std::cos(theta), y1 = (r + grid.r_inner) * std::sin(theta);
#pragma omp parallel for
  for (int j = 0; j < grid.nr; ++j)
  {
    const double radius = grid.r_inner + j * grid.dr;
    for (int i = 0; i < grid.ntheta; ++i)
    {
      const double dx = x1 - radius * grid_cos_theta[i], dy = y1 - radius * grid_sin_theta[i];
      f(j * grid.ntheta + i, dx * dx + dy * dy);
    }
  }
}

//Gaussians are not evaluated beyond exp(-gaussian_cutoff), about 4e-18,
//which is too small to change the fields they are added to.
static const double gaussian_cutoff = 40.;

/* Loop over all the cells and add heat according to where the current heat
   source is. */
void ConvectionSimulator::add_heat(double theta, double r, bool hot)
{
  IE_PROFILE(heat);
  //A Gaussian source term
  const double amplitude = (hot ? 1.0 : -1.0) * heat_source * dt;
  const double inv_two_variance = 1. / (2. * heat_source_radius * heat_source_radius);
  for_each_cell_distance(theta, r,
                         [&](int k, double rsq)
                         {
                           const double arg = rsq * inv_two_variance;
                           const double heat = (arg < gaussian_cutoff ? std::exp(-arg) : 0.);
                           T[k] = dmax(dmin(T[k] + amplitude * heat, 1.0), 0.0);
                         });
}

/* Loop over all the cells and add composition, similar to add_heat*/
void ConvectionSimulator::add_composition(double theta, double r)
{
  IE_PROFILE(heat);
  //A uniform source within composition_source_radius
  const double radius_sq = composition_source_radius * composition_source_radius;
  const double amount = composition_source * dt;
  for_each_cell_distance(theta, r,
                         [&](int k, double rsq)
                         {
                           const double composition = C[k] + (rsq < radius_sq ? amount : 0.0);
                           C[k] = dmax(dmin(composition, 1.0), 0.0);
                         });
}

/* Loop over the cells and zero out the displacement vectors */
void ConvectionSimulator::clear_seismic_waves()
{
  std::fill(D, D + grid.ncells, 0.);
  std::fill(Dp, Dp + grid.ncells, 0.);
}

/* Place the seismometer at a point (theta,r)*/
void ConvectionSimulator::place_seismometer(double theta, double r)
{
  seismometer_location.theta = theta;
  seismometer_location.r = r;
}

/* Add an earthquake at the point x,y.  I use a functional form of a
  Mexican hat wavelet.  At some point it may be good to revisit what
  functional form is best.  */
void ConvectionSimulator::earthquake(double theta, double r)
{
  const double earthquake_radius = grid.dr * 4.; //Somewhat arbitrary radius
  const double prefactor = 2.e-1 / M_E;

  for_each_cell_distance(theta, r,
                         [&](int k, double rsq)
                         {
                           const double dist = rsq / 2. / earthquake_radius / earthquake_radius;
                           D[k] += prefactor * (1. - dist) * std::exp(-dist);
                         });
}

/* Explicitly propagate the wave equation. */
void ConvectionSimulator::propagate_seismic_waves()
{
  IE_PROFILE(waves);
  const double reference_speed = 1.0; //dummy wavespeed
  const double tstep =
      0.6 * dmin(grid.r_inner * grid.dtheta, grid.dr) / reference_speed; // Timestep to satisfy cfl
  double dissipation = 0.5; //Empirically chosen dissipation
  double dispersion =
      1.0e-1; //Used to damp out any total displacment that is caused by lots of earthquake sources
  const int ntheta = grid.ntheta, nr = grid.nr;

//We can get away with an explicit timestepping scheme for the wave equation,
//so no complicated tridiagonal matrix inversion or any such nonsense here
#pragma omp parallel for
  for (int j = 0; j < nr; ++j)
    for (int i = 0; i < ntheta; ++i)
    {
      const int self = j * ntheta + i;
      const int right = j * ntheta + (i == ntheta - 1 ? 0 : i + 1);
      const int left = j * ntheta + (i == 0 ? ntheta - 1 : i - 1);
      const int up = self + ntheta, down = self - ntheta;

      double laplacian;
      const double r = grid.r_inner + j * grid.dr;
      const double dr = grid.dr;
      const double dtheta = grid.dtheta;

      //Make the wavespeed temperature dependent.  This is a HUGE
      //temperature dependence so it is pretty obvious.
      double speed = reference_speed * (1.0 - 0.7 * T[self]);
      speed = (speed < 0.1 ? 0.1 : (speed > 1.0 ? 1.0 : speed));

      //Five point laplacian stencil
      if (j == nr - 1) //Top boundary
        laplacian = (1. / dr / dr + 0.5 / r / dr) * D[down] +
                    (1. / dr / dr - 0.5 / r / dr) * D[down] +
                    1. / r / r / dtheta / dtheta * (D[right] + D[left]) -
                    (2. / r / r / dtheta / dtheta + 2. / dr / dr) * D[self];
      else if (j == 0) //Bottom boundary
        laplacian = (1. / dr / dr + 0.5 / r / dr) * D[up] + (1. / dr / dr - 0.5 / r / dr) * D[up] +
                    1. / r / r / dtheta / dtheta * (D[right] + D[left]) -
                    (2. / r / r / dtheta / dtheta + 2. / dr / dr) * D[self];
      else
        laplacian = (1. / dr / dr + 0.5 / r / dr) * D[up] +
                    (1. / dr / dr - 0.5 / r / dr) * D[down] +
                    1. / r / r / dtheta / dtheta * (D[right] + D[left]) -
                    (2. / r / r / dtheta / dtheta + 2. / dr / dr) * D[self];

      //Explicitly propagate wave
      scratch[self] = (2.0 * D[self] - Dp[self] + tstep * tstep * speed * speed * laplacian +
                       tstep * dissipation * D[self] - tstep * tstep * dispersion * D[self]) /
                      (1.0 + dissipation * tstep);
    }

#pragma omp parallel for
  for (int k = 0; k < grid.ncells; ++k)
  {
    //Current displacement becomes previous displacement
    Dp[k] = D[k];
    //Copy over new displacement
    D[k] = scratch[k];

    //clip field to tamp down on instabilities and keep it in reasonable values
    if (D[k] > 1.0)
      D[k] = 1.0;
    if (D[k] < -1.0)
      D[k] = -1.0;
  }
}

void ConvectionSimulator::semi_lagrangian_advect_temperature()
{
  semi_lagrangian_advect(advection_field::temperature);
}

void ConvectionSimulator::semi_lagrangian_advect_composition()
{
  semi_lagrangian_advect(advection_field::composition);
}

/* Advect the temperature field through the velocity field using
   Semi-lagrangian advection.  This scheme is quite stable, which
   allows me to take VERY large time steps.  The drawback is that it
   is kind of slow.  Here I do it in a very coarse way to make up for
   that.  A more accurate implementation would use better-than-linear
   interpolation and use more iterations. */
void ConvectionSimulator::semi_lagrangian_advect(const advection_field field)
{
  IE_PROFILE(advect);
  //The goal is to find the temperature at the Lagrangian point which will
  //be advected to the current grid point in one time step.  In general,
  //this point wil not be on the grid.  I find this point using a coarse
  //iterated predictor corrector.
  //
  //This is the innermost loop of the simulation, so it avoids the cell
  //iterators and trigonometry: cell indices are computed directly, the
  //grid angles' sines and cosines are precomputed, and those of the
  //takeoff points come from their Cartesian coordinates.

  const int ntheta = grid.ntheta, nr = grid.nr;
  const double dtheta = grid.dtheta, dr = grid.dr;
  const double r_inner = grid.r_inner, lr = grid.lr;
  const double half_dt = dt / 2.0;
  const double inv_dtheta = 1.0 / dtheta, inv_dr = 1.0 / dr;
  const double *cos_theta = grid_cos_theta.data(), *sin_theta = grid_sin_theta.data();
  const double *field_values = (field == advection_field::temperature ? T : C);

  //Locate a point (r measured from the inner radius, already clamped to
  //[0, lr]) on the grid: the cell at its lower left, the cells to the right
  //and above (wrapping in theta), and the local coordinates in that cell.
  struct Location
  {
    int self, right, up, upright;
    bool at_top;
    double local_theta, local_r;
  };
  //
  //The takeoff point's angle is found relative to the grid point's, which
  //is almost always a small angle, where a short series is far cheaper
  //than std::atan2 and agrees with it to about 1e-11 (|z| <= tan(pi/12),
  //where the first omitted term is z^17/17). Larger angles fall back to
  //std::atan2. The branch is almost never taken, so it predicts well.
  const double tan_pi_12 = 0.2679491924311227;
  auto locate = [&](double x, double y, double r, double theta_f, double cos_f, double sin_f)
  {
    //Coordinates in a frame rotated to the grid point's angle
    const double along = x * cos_f + y * sin_f, across = y * cos_f - x * sin_f;
    double offset;
    if (along > 0. && std::fabs(across) <= tan_pi_12 * along)
    {
      const double z = across / along, z2 = z * z;
      offset =
          z *
          (1. + z2 * (-1. / 3. +
                      z2 * (1. / 5. +
                            z2 * (-1. / 7. +
                                  z2 * (1. / 9. +
                                        z2 * (-1. / 11. + z2 * (1. / 13. + z2 * (-1. / 15.))))))));
    }
    else
      offset = std::atan2(across, along);

    //Angle in [0, 2 pi)
    double theta = theta_f + offset;
    if (theta < 0.)
      theta += grid.lx;
    else if (theta >= grid.lx)
      theta -= grid.lx;
    const double theta_cells = theta * inv_dtheta;
    int i = fast_floor(theta_cells);
    const double local_theta = theta_cells - i;
    if (i >= ntheta)
      i -= ntheta;
    else if (i < 0)
      i += ntheta;
    const double r_cells = r * inv_dr;
    int j = fast_floor(r_cells);
    j = (j < 0 ? 0 : (j >= nr ? nr - 1 : j));

    Location loc;
    loc.self = j * ntheta + i;
    loc.right = j * ntheta + (i + 1 == ntheta ? 0 : i + 1);
    loc.at_top = (j == nr - 1);
    loc.up = loc.self + ntheta;
    loc.upright = loc.right + ntheta;
    loc.local_theta = local_theta;
    loc.local_r = r_cells - j;
    return loc;
  };

#pragma omp parallel for collapse(2)
  for (int j = 0; j < nr; ++j)
    for (int i = 0; i < ntheta; ++i)
    {
      const int k = j * ntheta + i;

      //These points are known, as they are the grid points in question.  They will
      //not change for this cell.
      const double cos_f = cos_theta[i], sin_f = sin_theta[i];
      const double final_r = j * dr + r_inner;
      const double final_x = final_r * cos_f, final_y = final_r * sin_f;
      const double vel_final_x = -u[k] * sin_f + v[k] * cos_f;
      const double vel_final_y = u[k] * cos_f + v[k] * sin_f;

      //Calculate the initial guess for the takeoff point using a basic predictor
      //with a forward euler step
      double x = final_x - vel_final_x * dt;
      double y = final_y - vel_final_y * dt;
      double rho = std::sqrt(x * x + y * y);
      //Keep it in the domain
      double r = dmin(lr, dmax(rho - r_inner, 0.0));
      Location loc = locate(x, y, r, i * dtheta, cos_f, sin_f);

      //Iterate on the corrector.  Here I only do one iteration for
      //performance reasons, but in principle we could do more to get
      //a better estimate.
      {
        //Evaluate the velocity at the predictor, interpolating linearly.
        //At the top row there is nothing above, so use the top values.
        const int ul = (loc.at_top ? loc.self : loc.up);
        const int ur = (loc.at_top ? loc.right : loc.upright);
        const int dl = loc.self, dr_ = loc.right;
        const double vel_theta =
            linear_interp_2d(loc.local_theta, loc.local_r, u[ul], u[ur], u[dl], u[dr_]);
        const double vel_r =
            linear_interp_2d(loc.local_theta, loc.local_r, v[ul], v[ur], v[dl], v[dr_]);

        //The takeoff angle's cosine and sine, from its Cartesian position
        const double inv_rho = 1.0 / rho;
        const double cos_t = x * inv_rho, sin_t = y * inv_rho;
        const double vel_takeoff_x = -vel_theta * sin_t + vel_r * cos_t;
        const double vel_takeoff_y = vel_theta * cos_t + vel_r * sin_t;

        //Come up with the corrector using a midpoint rule
        x = final_x - (vel_final_x + vel_takeoff_x) * half_dt;
        y = final_y - (vel_final_y + vel_takeoff_y) * half_dt;
        rho = std::sqrt(x * x + y * y);
        //Keep in domain
        r = dmin(lr, dmax(rho - r_inner, 0.0));
        loc = locate(x, y, r, i * dtheta, cos_f, sin_f);
      }

      //Interpolate the field at the takeoff point
      const double *f = field_values;
      double value;
      if (field == advection_field::temperature)
      {
        //Above the top boundary the temperature is reflected about zero
        if (loc.at_top)
          value = linear_interp_2d(loc.local_theta, loc.local_r, -f[loc.self], -f[loc.right],
                                   f[loc.self], f[loc.right]);
        else
          value = linear_interp_2d(loc.local_theta, loc.local_r, f[loc.up], f[loc.upright],
                                   f[loc.self], f[loc.right]);
      }
      else
      {
        if (loc.at_top)
          value = linear_interp_2d(loc.local_theta, loc.local_r, f[loc.self], f[loc.right],
                                   f[loc.self], f[loc.right]);
        else
          value = linear_interp_2d(loc.local_theta, loc.local_r, f[loc.up], f[loc.upright],
                                   f[loc.self], f[loc.right]);
      }
      scratch[k] = value;
    }

  //Copy the scratch vector into the temperature or composition vector,
  //each thread the part it just computed
  double *destination = (field == advection_field::temperature ? T : C);
#pragma omp parallel for
  for (int k = 0; k < grid.ncells; ++k)
    destination[k] = scratch[k];
}

/* Solve the diffusion equation implicitly with reverse Euler.
   I do the x and y direction separately, with the x direction
   done spectrally using FFTs, and the y direction with finite
   differences, solving the resulting tridiagonal system using
   the O(N) Thomas algorithm. */
void ConvectionSimulator::diffuse_temperature()
{
  //Apply boundary conditions
  {
    IE_PROFILE(diffuse_misc);
    std::fill(T, T + grid.ntheta, 1.0);                               //Bottom row
    std::fill(T + (grid.nr - 1) * grid.ntheta, T + grid.ncells, 0.0); //Top row
  }

  //Execute the forward fourier transform
  {
    IE_PROFILE(diffuse_fft);
    fft_rows(T, T_spectral); //X direction
  }

  //Solve in place, normalizing for the unnormalized inverse transform
  {
    IE_PROFILE(diffuse_tridiag);
    const double normalization = 1.0 / grid.ntheta;
#pragma omp parallel for
    for (int b = 0; b < n_mode_blocks(); ++b)
      diffusion_matrices.solve(T_spectral, T_spectral, grid.ntheta, b * mode_block,
                               std::min((b + 1) * mode_block, n_modes()), normalization);
  }

  //Execute the inverse Fourier transform
  {
    IE_PROFILE(diffuse_fft);
    ifft_rows(T_spectral, T); //X direction
  }
}

void ConvectionSimulator::setup_diffusion_problem()
{
  //Setup the tridiagonal matrices in the y direction
  std::vector<double> upper_diag(grid.nr);
  std::vector<double> lower_diag(grid.nr);

  const double dr = grid.dr;
  //Upper and lower diagonals are the same for each frequency
  for (int i = 0; i < grid.nr; ++i)
  {
    const double r = grid.r_inner + i * dr;
    upper_diag[i] = -dt / dr / dr * (1. + 0.5 * dr / r);
    lower_diag[i] = -dt / dr / dr * (1. - 0.5 * dr / r);
  }
  upper_diag[0] = 0.0;           //fix for lower B.C.
  lower_diag[grid.nr - 1] = 0.0; //fix for upper B.C.

//Diagonals include a term for the frequency in the x-direction
#pragma omp parallel for
  for (int l = 0; l <= grid.ntheta / 2; ++l)
  {
    double factor = (l * l);
    std::vector<double> diag(grid.nr);
    diag[0] = 1.0;
    diag[grid.nr - 1] = 1.0;
    for (int i = 1; i < grid.nr - 1; ++i)
    {
      const double r = grid.r_inner + i * dr;
      diag[i] = 1. + 2. * dt / dr / dr + factor * dt / r / r;
    }

    diffusion_matrices.initialize(l, lower_diag.data(), diag.data(), upper_diag.data());
  }
}

/* Setup the stokes solve with the stream function formulation.
   We assemble the vector that stores the frequencies of the
   eigenmodes, as well as tell FFTW how to do the transforms */
void ConvectionSimulator::setup_stokes_problem()
{
  std::vector<double> upper_diag(grid.nr);
  std::vector<double> diag(grid.nr);
  std::vector<double> lower_diag(grid.nr);

  const double dr = grid.dr;
  for (int i = 0; i < grid.nr; ++i)
  {
    const double r = grid.r_inner + i * dr;
    upper_diag[i] = 1. / dr / dr + 0.5 / dr / r;
    lower_diag[i] = 1. / dr / dr - 0.5 / dr / r;
  }
  upper_diag[0] = 0.0;           //fix for lower B.C.
  lower_diag[grid.nr - 1] = 0.0; //fix for upper B.C.

  for (int l = 0; l <= grid.ntheta / 2; ++l)
  {
    double factor = (l * l);
    diag[0] = 1.0;
    diag[grid.nr - 1] = 1.0;
    for (int i = 1; i < grid.nr - 1; ++i)
    {
      const double r = grid.r_inner + i * dr;
      diag[i] = -2. / dr / dr - factor / r / r;
    }

    stokes_matrices.initialize(l, lower_diag.data(), diag.data(), upper_diag.data());
  }
}

/* Tell FFTW how to do the transforms in the x direction: a real DFT of
   one row of the grid. The same plans serve every row of every array,
   through FFTW's new-array execute functions (see fft_rows), which are
   safe to call from several threads at once. They are made once, since
   planning is expensive (and with anything other than FFTW_ESTIMATE,
   overwrites the arrays). */
void ConvectionSimulator::setup_fft_plans()
{
  //New-array execution needs every row to have the alignment the plans
  //were made with. That holds for the usual power of two ntheta, since the
  //arrays come from fftw_alloc; otherwise, give up on aligned SIMD.
  //FFTW_MEASURE finds plans about 10% faster, but that is under 1% of a
  //step, costs about 0.3 s of planning, and makes results depend on timing.
  const bool rows_aligned = (grid.ntheta * sizeof(double)) % 64 == 0;
  const unsigned flags = FFTW_ESTIMATE | (rows_aligned ? 0u : unsigned(FFTW_UNALIGNED));

  row_dft =
      fftw_plan_dft_r2c_1d(grid.ntheta, T, reinterpret_cast<fftw_complex *>(T_spectral), flags);
  row_idft = fftw_plan_dft_c2r_1d(grid.ntheta, reinterpret_cast<fftw_complex *>(T_spectral),
                                  scratch, flags);
}

/* Forward and inverse transforms of every row. The rows are independent,
   so they are split between threads; each one is computed the same way
   regardless of the number of threads. Like all of FFTW's real inverse
   transforms, the inverse overwrites its input. */
void ConvectionSimulator::fft_rows(const double *in, std::complex<double> *out)
{
#pragma omp parallel for
  for (int j = 0; j < grid.nr; ++j)
    fftw_execute_dft_r2c(row_dft, const_cast<double *>(in) + j * grid.ntheta,
                         reinterpret_cast<fftw_complex *>(out + j * grid.ntheta));
}

void ConvectionSimulator::ifft_rows(std::complex<double> *in, double *out)
{
#pragma omp parallel for
  for (int j = 0; j < grid.nr; ++j)
    fftw_execute_dft_c2r(row_idft, reinterpret_cast<fftw_complex *>(in + j * grid.ntheta),
                         out + j * grid.ntheta);
}

/*Calculate the Nusselt number for the given temperature field.
  In principle I could calculate the heat flux through the bottom
  boundary as well, and take the average, but as it is, I am only
  doing the heat flux through the top boundary*/
double ConvectionSimulator::nusselt_number()
{
  double heat_flux = 0;
  const double *top = T + (grid.nr - 1) * grid.ntheta, *below = top - grid.ntheta;
  for (int i = 0; i < grid.ntheta; ++i)
    heat_flux += (top[i] - below[i]) / grid.dr;

  return -heat_flux *                           //sum of fluxes
         grid.r_outer / grid.ntheta *           //metric terms
         std::log(grid.r_outer / grid.r_inner); //geometry factor.
}

//The curl of the temperature is what is relevant for the stream
//function calculation.  This calculates that curl.
void ConvectionSimulator::assemble_curl_density_vector()
{
  IE_PROFILE(stokes_rhs);
  const int ntheta = grid.ntheta;
  //Centered differences in theta of each row of the field, scaled by factor
  //and added to out
  auto add_dtheta = [ntheta](const double *field, double factor, double *out)
  {
    out[0] += factor * (field[1] - field[ntheta - 1]);
    for (int i = 1; i < ntheta - 1; ++i)
      out[i] += factor * (field[i + 1] - field[i - 1]);
    out[ntheta - 1] += factor * (field[0] - field[ntheta - 2]);
  };

#pragma omp parallel for
  for (int j = 0; j < grid.nr; ++j)
  {
    double *out = curl_density + j * ntheta;
    std::fill(out, out + ntheta, 0.0);
    if (j == 0 || j == grid.nr - 1)
      continue; //Zero at the boundaries

    //Ra / r times the theta derivative of the density
    const double r = grid.r_inner + j * grid.dr;
    const double factor = Ra / r / 2. / grid.dtheta;
    add_dtheta(T + j * ntheta, factor, out);
    if (include_composition)
      add_dtheta(C + j * ntheta, -buoyancy_number * factor, out);
  }
}

/* Actually solve the biharmonic equation for the Stokes system.*/
void ConvectionSimulator::solve_stokes()
{
  //Come up with the RHS of the spectral solve
  assemble_curl_density_vector();

  //Execute the forward fourier transform
  {
    IE_PROFILE(stokes_fft);
    fft_rows(curl_density, curl_density_spectral); //X direction
  }

  //The biharmonic operator is the square of the Laplacian, so solve with
  //the Laplacian twice, each block of modes while it is in cache. Solve in
  //place, normalizing for the unnormalized inverse transform.
  {
    IE_PROFILE(stokes_tridiag);
    const double normalization = 1.0 / grid.ntheta;
#pragma omp parallel for
    for (int b = 0; b < n_mode_blocks(); ++b)
    {
      const int first = b * mode_block, last = std::min((b + 1) * mode_block, n_modes());
      stokes_matrices.solve(curl_density_spectral, curl_density_spectral, grid.ntheta, first, last);
      stokes_matrices.solve(curl_density_spectral, curl_density_spectral, grid.ntheta, first, last,
                            normalization);
    }
  }

  //Execute the inverse Fourier transform
  {
    IE_PROFILE(stokes_fft);
    ifft_rows(curl_density_spectral, stream); //X direction
  }

  //Come up with the velocities by taking finite differences of the stream function.
  //I could also take these derivatives in spectral space, but that would mean
  //more fourier transforms, so this should be considerably cheaper.
  //At the boundaries, the radial derivative is one sided.
  IE_PROFILE(stokes_velocity);
  const int ntheta = grid.ntheta;
#pragma omp parallel for
  for (int j = 0; j < grid.nr; ++j)
  {
    const int above = (j == grid.nr - 1 ? j : j + 1), below = (j == 0 ? j : j - 1);
    const double u_factor = 1.0 / ((above - below) * grid.dr);
    const double r = grid.r_inner + j * grid.dr;
    const double v_factor = -1.0 / (2.0 * grid.dtheta * r);
    const double *s_row = stream + j * ntheta;
    const double *s_above = stream + above * ntheta, *s_below = stream + below * ntheta;
    double *u_row = u + j * ntheta, *v_row = v + j * ntheta;
    for (int i = 0; i < ntheta; ++i)
    {
      const int left = (i == 0 ? ntheta - 1 : i - 1), right = (i == ntheta - 1 ? 0 : i + 1);
      u_row[i] = (s_above[i] - s_below[i]) * u_factor;
      v_row[i] = (s_row[right] - s_row[left]) * v_factor;
    }
  }
}

void ConvectionSimulator::true_polar_wander()
{
  IE_PROFILE(tpw);

  double moment_of_inertia[3] = {0., 0., 0.};
  double omega[2] = {std::cos(spin), std::sin(spin)};

  //First, calculate the moment of inertia tensor. Each row is summed by
  //one thread, then the rows in order, so the result does not depend on
  //the number of threads.
  std::vector<double> row_moments(3 * grid.nr);
#pragma omp parallel for
  for (int j = 0; j < grid.nr; ++j)
  {
    const double r = grid.r_inner + j * grid.dr;
    const double volume = r * grid.dr * grid.dtheta;
    double yy = 0., xx = 0., xy = 0.;
    for (int i = 0; i < grid.ntheta; ++i)
    {
      const int k = j * grid.ntheta + i;
      double density = -T[k];
      if (include_composition)
        density = -(T[k] - buoyancy_number * C[k]);

      //The TPW calculation can go wonky when the temperature gets
      //cold enough that the sign changes. Not really sure what is
      //going on at the moment, but this fixes it
      density = dmin(density, 0.0);

      const double x = r * grid_cos_theta[i];
      const double y = r * grid_sin_theta[i];
      yy += density * y * y;
      xx += density * x * x;
      xy += density * x * y;
    }
    row_moments[3 * j + 0] = yy * volume;
    row_moments[3 * j + 1] = xx * volume;
    row_moments[3 * j + 2] = xy * volume;
  }
  for (int j = 0; j < grid.nr; ++j)
  {
    moment_of_inertia[0] += row_moments[3 * j + 0]; // E_11
    moment_of_inertia[1] += row_moments[3 * j + 1]; // E_22
    moment_of_inertia[2] -= row_moments[3 * j + 2]; // E_21 E12
  }

  //Not strictly necessary, but doesn't hurt to subtract the trace
  double reference_moment = (moment_of_inertia[0] + moment_of_inertia[1]) / 2.;
  moment_of_inertia[0] -= reference_moment;
  moment_of_inertia[1] -= reference_moment;

  double omega_dot[2];
  double E_dot_omega[2];
  double omega_dot_E_dot_omega_omega[2];

  //Break up the advectin into a bunch of small forward Euler steps.
  //This is way inefficient, but overall not the most expensive thing
  //in the model, so probably not a big deal.
  //TODO: figure out a better way to timestep the nonlinear ODE
  unsigned int n_euler_steps = 100;
  double small_dt = dt / n_euler_steps;
  for (unsigned int i = 0; i < n_euler_steps; ++i)
  {
    //Compute rotational stress
    E_dot_omega[0] = moment_of_inertia[0] * omega[0] + moment_of_inertia[2] * omega[1];
    E_dot_omega[1] = moment_of_inertia[2] * omega[0] + moment_of_inertia[1] * omega[1];
    omega_dot_E_dot_omega_omega[0] =
        (omega[0] * E_dot_omega[0] + omega[1] * E_dot_omega[1]) * omega[0];
    omega_dot_E_dot_omega_omega[1] =
        (omega[0] * E_dot_omega[0] + omega[1] * E_dot_omega[1]) * omega[1];

    double Fr = 0.003;
    double factor = 2. / 19.; //Tuning factor to make the rate look good.

    omega_dot[0] = factor * Ra / Fr * (E_dot_omega[0] - omega_dot_E_dot_omega_omega[0]);
    omega_dot[1] = factor * Ra / Fr * (E_dot_omega[1] - omega_dot_E_dot_omega_omega[1]);

    omega[0] += omega_dot[0] * small_dt;
    omega[1] += omega_dot[1] * small_dt;
    double omega_norm = std::sqrt(omega[0] * omega[0] + omega[1] * omega[1]);
    omega[0] /= omega_norm;
    omega[1] /= omega_norm;
  }

  spin = std::atan2(omega[1], omega[0]);
}

void ConvectionSimulator::update_state(double rayleigh)
{
  double Ra_c = 657.; //critical rayleigh number
  double length_scale =
      std::pow(rayleigh / 2. / Ra_c, -1. / 3.) * grid.lr; //calculate a provisional length scale

  //The resolution basically sets the maximum Ra we can use.  Estimate the minimum length scale,
  //and if that is smaller than the resolution, cap the Rayleigh number.
  const double boundary_layer_cells = 4.; //grid cells per boundary layer
  if (length_scale < boundary_layer_cells * grid.dr)
    Ra = 2. * Ra_c * std::pow(grid.lr / boundary_layer_cells / grid.dr, 3.0);
  else
    Ra = rayleigh;

  length_scale = std::pow(Ra / 2. / Ra_c, -1. / 3.) * grid.lr;
  const double Nu = std::pow(Ra / Ra_c / 2., 1. / 3.) / 2.; //Nusselt
  const double velocity_scale = std::sqrt(Ra * Nu);
  const double cfl = dmin(grid.dr, grid.r_inner * grid.dtheta) / velocity_scale;

  //Estimate other state properties based on simple isoviscous scalings
  dt = cfl * 20.0;                             //Roughly 10x CFL, thanks to semi-lagrangian
  heat_source_radius = length_scale * 0.5;     //Radius of order the boundary layer thickness
  heat_source = velocity_scale / grid.lr * 2.; //Heat a blob of order the ascent time for thta blob
  composition_source_radius = grid.lr / 10.;   //Size of compositional blob
  composition_source = heat_source;            //Strength of compositional reaction

  setup_diffusion_problem(); //The diffusion matrices depend on dt
}

void ConvectionSimulator::set_temperature(const double *values)
{
  for (int i = 0; i < grid.ncells; ++i)
    T[i] = values[i];
}

double ConvectionSimulator::spin_angle() const
{
  return spin;
}

double ConvectionSimulator::rayleigh_number() const
{
  return Ra;
}

double ConvectionSimulator::seismometer_reading() const
{
  return D[grid.cell_id(seismometer_location)];
}

void ConvectionSimulator::seismometer_position(double &theta, double &r) const
{
  theta = seismometer_location.theta;
  r = seismometer_location.r;
}

double ConvectionSimulator::timestep() const
{
  return dt;
}

double ConvectionSimulator::timescale() const
{
  double Ra_c = 657.;                                       //critical rayleigh number
  const double Nu = std::pow(Ra / Ra_c / 2., 1. / 3.) / 2.; //Nusselt
  const double velocity_scale = std::sqrt(Ra * Nu);
  return grid.lr / velocity_scale; //Approximately the ascent time for a plume
}
