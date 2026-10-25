#include "rendering.h"
#include "rendering_plugins.h"
#include "color.h"
#include "profiler.h"
#include "harness.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <vector>

#include "sokol_app.h"

#ifdef __EMSCRIPTEN__
#include <emscripten/bind.h>
#endif

/*********************************************
    MODIFY THSE PARAMETERS TO CHANGE THE
    BEHAVIOR OF THE SIMULATION.
*********************************************/

//Whether to include a chemical field.
//The simulation will be faster without an additional advected field.
bool include_composition = false;

//Whether to do TPW calculation
bool include_tpw = false;

//Number of cells in the theta and r directions.
//This is the primary control on resolution,
//as well as performance, and sets the maximum
//Rayleigh number (which goes as nr^3).
//Override with IE_GRID=NTHETAxNR, e.g. IE_GRID=1536x192
//(on the web, the grid parameter of the page URL).
struct GridSize
{
  int ntheta, nr;
};
GridSize grid_size()
{
  GridSize size = {1024, 128};
  const char *setting = std::getenv("IE_GRID");
  if (setting && *setting)
  {
    int ntheta = 0, nr = 0;
    if (std::sscanf(setting, "%dx%d", &ntheta, &nr) == 2 && ntheta >= 8 && nr >= 4)
      size = {ntheta, nr};
    else
      std::fprintf(stderr, "Ignoring IE_GRID=%s; expected NTHETAxNR, such as 1536x192\n", setting);
  }
  return size;
}
const GridSize grid = grid_size();

//Aspect ratio of the computational domain
//is set by the inner radius, where the outer
//radius is assumed to be 1.0
const double r_inner = 0.5;

//Render the simulation flattened, as if
//acting under centrifugal forces.
const double flattening = 0.0;

/*********************************************
    PROBABLY DON'T MODIFY THE REST
    UNLESS YOU KNOW WHAT YOU ARE DOING.
*********************************************/

//Size of computational domain.  The ratio of
//ntheta to nr should be roughly the same
//as ltheta to lr to avoid a funny shaped grid,
//which will result in a strange looking simulation.
const double ltheta = 2. * M_PI;
const double lr = 1.0 - r_inner;

//The simulation is drawn in a square viewport centered in the window.
//Its lower left corner (in GL coordinates) and size, in pixels.
int view_x = 0, view_y = 0, view_size = 1;
//Size of one logical (unscaled) screen pixel in clip coordinates
float point_size = 2.f;

//Whether to add heat to the simulation on mouse click.
int click_state = 0;
//Location of heat-adding, chemistry adding, or earthquake generation
double click_theta, click_r;

//Whether we are in earthquake mode
bool seismic_mode = false;

//Whether to solve the advection-diffusion equation
bool advection_diffusion = true;

//Whether to draw composition or temperature fields
bool draw_composition = false;

//Function pointer to which color scale we are using
color (*colormap)(double) = &hot;

//The simulation time.
double simulation_time = 0.0;

//Global solver
ConvectionSimulator simulator(r_inner, grid.ntheta, grid.nr, include_composition);
Mantle mantle(simulator);
Axis axis(simulator);
Seismograph seismograph(simulator);

//Deterministic screenshot mode for regression checks. If IE_SCREENSHOT
//is set, seed convection with fixed heat sources, run IE_SCREENSHOT_STEPS
//timesteps, write the framebuffer to that path as a PPM, and quit.
const char *screenshot_path = std::getenv("IE_SCREENSHOT");
const int screenshot_steps =
    std::getenv("IE_SCREENSHOT_STEPS") ? std::atoi(std::getenv("IE_SCREENSHOT_STEPS")) : 200;

void save_screenshot(const char *path)
{
  const int n = view_size;
  std::vector<unsigned char> pixels(3 * n * n);
  glPixelStorei(GL_PACK_ALIGNMENT, 1);
  glReadPixels(view_x, view_y, n, n, GL_RGB, GL_UNSIGNED_BYTE, pixels.data());
  FILE *f = std::fopen(path, "wb");
  if (!f)
  {
    std::cerr << "Could not open " << path << '\n';
    return;
  }
  std::fprintf(f, "P6\n%d %d\n255\n", n, n);
  //OpenGL rows start at the bottom
  for (int row = n - 1; row >= 0; --row)
    std::fwrite(&pixels[3 * n * row], 1, 3 * n, f);
  std::fclose(f);
}

//Convert window coordinates (in pixels, from the top left) into
//coordinates in [-0.5, 0.5] relative to the square viewport, y up.
inline void window_to_view(const float wx, const float wy, float *x, float *y)
{
  const float top = sapp_heightf() - float(view_y + view_size);
  *x = (wx - float(view_x)) / float(view_size) - 0.5f;
  *y = 1.0f - (wy - top) / float(view_size) - 0.5f;
}

//Given an x and y location, compute the r, theta location of the
//simulator domain. Can handle the case of ellipticity.
inline void compute_simulator_location(const float x, const float y, float *theta, float *r)
{
  //Scale the click coordinates so that they are
  //in the correct place if the domain is elliptical
  float xpp = x, ypp = y;
  if (flattening > 1.e-2)
  {
    //Rotate the click position so that the ellipse
    //corresponds to the Cartesian axes.
    float rot_angle = simulator.spin_angle() - M_PI / 2.;
    float xp = x * std::cos(-rot_angle) - y * std::sin(-rot_angle);
    float yp = x * std::sin(-rot_angle) + y * std::cos(-rot_angle);

    //Do the inverse flattening (only along y)
    yp = yp / (1. - flattening);

    //Rotate back to the initial angle
    xpp = xp * std::cos(rot_angle) - yp * std::sin(rot_angle);
    ypp = xp * std::sin(rot_angle) + yp * std::cos(rot_angle);
  }
  *theta = std::atan2(ypp, xpp);
  *theta = (*theta < 0. ? *theta + 2. * M_PI : *theta);
  *r = 2. * std::sqrt(xpp * xpp + ypp * ypp);
}

//Update where to add heat
inline void handle_mouse_motion(const sapp_event *event)
{
  float x, y;
  window_to_view(event->mouse_x, event->mouse_y, &x, &y);

  float theta, r;
  compute_simulator_location(x, y, &theta, &r);

  click_theta = ltheta * theta / 2. / M_PI;
  click_r = lr * (r - r_inner) / (1.0f - r_inner);
}

//Change the Rayleigh number on scrolling
inline void handle_mouse_wheel(const sapp_event *event)
{
  if (event->scroll_y == 0.f)
    return;
  double rayleigh = simulator.rayleigh_number();
  double factor = std::pow(10.0, 1. / 100. * (event->scroll_y < 0 ? -1. : 1.0));
  simulator.update_state(rayleigh * factor);
}

//Toggle whether to add heat, and whether it should
//be positive or negative
inline void handle_mouse_button(const sapp_event *event)
{
  if (event->type == SAPP_EVENTTYPE_MOUSE_DOWN)
  {
    //Ctrl-click acts as a right click, for one-button mice and trackpads
    const bool ctrl = event->modifiers & SAPP_MODIFIER_CTRL;
    if (event->mouse_button == SAPP_MOUSEBUTTON_LEFT)
      click_state = ctrl ? -1 : 1;
    if (event->mouse_button == SAPP_MOUSEBUTTON_RIGHT)
      click_state = -1;

    float x, y;
    window_to_view(event->mouse_x, event->mouse_y, &x, &y);

    float theta, r;
    compute_simulator_location(x, y, &theta, &r);

    click_theta = theta;
    click_r = r - r_inner;
  }
  else
  {
    click_state = 0;
  }
}

inline bool in_domain(const float /*theta*/, const float r)
{
  return (r + r_inner < 1.) && (r > 0.);
}

void cycle_colorscale()
{
  const unsigned int n_colormaps = 5;
  static color (*maps[n_colormaps])(double) = {&hot, &viridis, &inferno, &seismic, &gist_earth};
  static unsigned short i = 0;
  i = (i + 1) % n_colormaps;
  colormap = maps[i];
}

void text_output()
{
  //rough timescale conversion for Earth's mantle:
  //L^2/kappa / seconds/Gyr
  const double timescale = 2800.e3 * 2800.e3 / 1.e-6 / M_PI / 1.e7 / 1.e9;
  //Unicode characters for exponentiation
  const char *exponents[10] = {"\u2070", "\u00B9", "\u00B2", "\u00b3", "\u2074",
                               "\u2075", "\u2076", "\u2077", "\u2078", "\u2079"};
  //Output Rayleigh information
  const double Ra = simulator.rayleigh_number();
  if (Ra < 1000)
    std::cout << "Ra: " << std::setprecision(3) << Ra << "\t";
  else
  {
    const int upper = int(std::floor(std::log10(Ra)));
    const int lower = int(rint(Ra / std::pow(10, upper)));
    if (lower == 1)
      std::cout << "Ra: " << "  10" << exponents[upper] << "\t";
    else
      std::cout << "Ra: " << lower << "x10" << exponents[upper] << "\t";
  }

  //Output Nusselt information
  std::cout << "Nu: " << std::setprecision(2) << simulator.nusselt_number() << "\t";

  //Output time information
  if (simulation_time * timescale < 1.)
    std::cout << "\tTime: " << int(rint(simulation_time * timescale * 1000.)) << " million years\t";
  else
    std::cout << "\tTime: " << std::setprecision(2) << simulation_time * timescale
              << " billion years\t";

  //Output performance information, averaged over recent frames
  const Profiler &prof = profiler();
  const double diffuse = prof.mean(Phase::diffuse_fft) + prof.mean(Phase::diffuse_tridiag) +
                         prof.mean(Phase::diffuse_misc);
  const double stokes = prof.mean(Phase::stokes_rhs) + prof.mean(Phase::stokes_fft) +
                        prof.mean(Phase::stokes_tridiag) + prof.mean(Phase::stokes_velocity);
  std::cout << std::fixed << std::setprecision(1) << "| step " << prof.mean(Phase::step)
            << " ms (stokes " << stokes << ", advect " << prof.mean(Phase::advect) << ", diffuse "
            << diffuse << ")  x " << prof.mean(Phase::simulate) / prof.mean(Phase::step)
            << " per frame  | render " << prof.mean(Phase::render_prep) << " + "
            << prof.mean(Phase::render_gl) << " ms  | " << 1000. / prof.mean(Phase::frame)
            << " fps";
  std::cout.unsetf(std::ios_base::floatfield);

  //End the line, and flush so the output appears even when piped
  std::cout << '\n' << std::flush;
}

//Number of convection timesteps taken
int step_count = 0;

//One convection timestep
void convection_step()
{
  const double step_start = now_ms();

  //At the moment, the stokes solve is not the limiting factor,
  //so it does not hurt to do it every timestep
  simulator.solve_stokes();

  //Add heat if the user is clicking
  if (click_state != 0 && ((include_composition && !draw_composition) || (!include_composition)) &&
      in_domain(click_theta, click_r))
    simulator.add_heat(click_theta, click_r, (click_state == 1 ? true : false));
  if (click_state != 0 && include_composition && draw_composition &&
      in_domain(click_theta, click_r))
    simulator.add_composition(click_theta, click_r);

  //Same fixed heat sources as the headless benchmark
  if (screenshot_path && step_count < seed_steps)
    add_seed_heat_sources(simulator);

  //The user can do some neat painting by turning off advection and diffusion
  if (advection_diffusion)
  {
    //Advect temperature and composition fields
    simulator.semi_lagrangian_advect_temperature();
    if (include_composition)
      simulator.semi_lagrangian_advect_composition();

    //Diffuse temperature
    simulator.diffuse_temperature();
  }

  //Do TPW
  if (include_tpw)
    simulator.true_polar_wander();

  //One sample per step for the solver phases
  profiler().accumulate(Phase::step, now_ms() - step_start);
  profiler().commit();

  ++step_count;
  simulation_time += simulator.timestep();
}

//One step of the seismic waves
void seismic_step()
{
  //Make earthquakes
  if (click_state == 1 && in_domain(click_theta, click_r))
    simulator.earthquake(click_theta, click_r);
  else if (click_state == -1 && in_domain(click_theta, click_r))
    simulator.place_seismometer(click_theta, click_r);
  //Propagate waves
  simulator.propagate_seismic_waves();
}

/* Frame pacing. Drawing waits for the display's refresh (vsync), so with
   one step per frame a fast machine sits idle for most of each frame.
   Instead, take as many convection steps as fit in a budget, leaving the
   rest of the frame for rendering and the system. When a step takes
   longer than that (slow machines, large grids), this is one step per
   frame, as before.

   The budget is for 60 frames per second, whatever the display: on faster
   displays, simulation gets the time that more frames would have taken.
   (Estimating the refresh rate from frame intervals is unreliable, since
   compositors make the intervals jitter by several milliseconds.) The
   budget is measured after drawing, which is where the wait for vsync
   happens, so most of the frame is still ahead. */
const double step_budget_ms = 0.8 * 1000. / 60.;
const int max_steps_per_frame = 8;

void simulate()
{
  if (seismic_mode)
  {
    //The waves and the seismograph are paced by the frame
    seismic_step();
    return;
  }

  //Screenshots must be deterministic, so take one step per frame
  const int max_steps = (screenshot_path ? 1 : max_steps_per_frame);
  const double budget = step_budget_ms;
  const double start = now_ms();
  int steps = 0;
  double last_step = 0.;
  do
  {
    const double step_start = now_ms();
    convection_step();
    last_step = now_ms() - step_start;
    ++steps;
  } while (steps < max_steps && (now_ms() - start) + last_step < budget);
  profiler().accumulate(Phase::simulate, now_ms() - start);
}

//Set up all the OpenGL stuff, once the window and context exist
void init()
{
  mantle.setup();
  if (include_tpw)
    axis.setup();
  seismograph.setup();
}

//Cleanup
void cleanup()
{
  mantle.cleanup();
  if (include_tpw)
    axis.cleanup();
  seismograph.cleanup();
}

void frame()
{
  //Fit the largest square that we can in the window. Sizes are in
  //framebuffer pixels, which may be smaller than logical pixels on
  //high density displays.
  const int width = sapp_width(), height = sapp_height();
  view_size = (width < height ? width : height);
  view_x = (width - view_size) / 2;
  view_y = (height - view_size) / 2;
  glViewport(view_x, view_y, view_size, view_size);
  point_size = 2.f * sapp_dpi_scale() / float(view_size);

  //Draw the current state
  mantle.draw(include_composition && draw_composition);
  {
    IE_PROFILE(render_gl);
    if (include_tpw && !seismic_mode)
      axis.draw();
    if (seismic_mode)
      seismograph.draw(point_size);
  }

  if (screenshot_path && step_count == screenshot_steps)
  {
    save_screenshot(screenshot_path);
    sapp_quit();
    return;
  }

  simulate();

  //Time between frames, including any wait for vsync
  static double last_frame = now_ms();
  const double frame_start = now_ms();
  profiler().accumulate(Phase::frame, frame_start - last_frame);
  last_frame = frame_start;
  profiler().commit();

#ifndef __EMSCRIPTEN__
  //Output scaling and performance information about once a second
  static double last_output = frame_start;
  if (!seismic_mode && frame_start - last_output >= 1000.)
  {
    text_output();
    last_output = frame_start;
  }
#endif
}

void event(const sapp_event *event)
{
  //Screenshot mode must be deterministic, so ignore user input
  if (screenshot_path)
    return;

  switch (event->type)
  {
  case SAPP_EVENTTYPE_KEY_DOWN:
    //Some systems seem to have a hair trigger on the
    //keypress repeats. Filter those out.
    if (event->key_repeat)
      break;
    switch (event->key_code)
    {
    case SAPP_KEYCODE_SPACE:
      seismic_mode = !seismic_mode;
      simulator.clear_seismic_waves();
      seismograph.clear_record();
      break;

#ifndef __EMSCRIPTEN__
    //Quitting should be handled by navigating to another
    //webpage or closing the browser when using, emscripten,
    //so just disable quitting on escape
    case SAPP_KEYCODE_ESCAPE:
      sapp_request_quit();
      break;
#endif
    case SAPP_KEYCODE_TAB:
      draw_composition = !draw_composition;
      sapp_consume_event(); //Don't move the browser's focus
      break;

    case SAPP_KEYCODE_BACKSPACE:
      advection_diffusion = !advection_diffusion;
      sapp_consume_event(); //Don't navigate back in the browser
      break;

    case SAPP_KEYCODE_C:
      cycle_colorscale();
      break;

    default:
      break;
    }
    break;
  case SAPP_EVENTTYPE_MOUSE_DOWN:
  case SAPP_EVENTTYPE_MOUSE_UP:
    handle_mouse_button(event);
    break;
  case SAPP_EVENTTYPE_MOUSE_MOVE:
    handle_mouse_motion(event);
    break;
  case SAPP_EVENTTYPE_MOUSE_SCROLL:
    handle_mouse_wheel(event);
    sapp_consume_event(); //Don't scroll the page
    break;
  default:
    break;
  }
}

sapp_desc sokol_main(int /*argc*/, char ** /*argv*/)
{
  std::printf("Grid %dx%d, maximum Ra %.2g\n", grid.ntheta, grid.nr, simulator.rayleigh_number());

  sapp_desc desc = {};
  desc.init_cb = init;
  desc.frame_cb = frame;
  desc.cleanup_cb = cleanup;
  desc.event_cb = event;
  desc.window_title = "Convection";
  desc.gl.major_version = 3;
  desc.gl.minor_version = 3;
  desc.high_dpi = true; //Render at full resolution on high density displays
  return desc;
}

#ifdef __EMSCRIPTEN__
//Emscripten bindings so that the Javascript on the page can
//query the simulator for its state.

//Embind requires a raw function pointer, so wrap the simulator
//member function call in a function
double emscripten_rayleigh()
{
  return simulator.rayleigh_number();
}
std::string emscripten_profile()
{
  return profiler().json();
}
double emscripten_nusselt()
{
  return simulator.nusselt_number();
}
EMSCRIPTEN_BINDINGS(my_module)
{
  emscripten::function("rayleigh", &emscripten_rayleigh);
  emscripten::function("nusselt", &emscripten_nusselt);
  emscripten::function("profile", &emscripten_profile);
}
#endif
