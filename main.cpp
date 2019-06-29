#include "GL/glew.h"
#include "SDL2/SDL_opengl.h"
#include "SDL2/SDL.h"
#include "convection.h"
#include "rendering_plugins.h"
#include "color.h"
#include <cmath>
#include <iostream>
#include <iomanip>

#ifdef __EMSCRIPTEN__
#include <emscripten/emscripten.h>
#include <emscripten/bind.h>
#include <emscripten/html5.h>
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
//as well as performance.
const unsigned int ntheta = 512;
const unsigned int nr = 64;

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
const double ltheta = 2.*M_PI;
const double lr = 1.0-r_inner;

//Total number of pixels in x and y directions
int xpix;
int ypix;

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

// Button geometries
const float mode_button_x = -0.9;
const float mode_button_y = -0.9;
const float mode_button_radius = 0.1;



//Global solver
ConvectionSimulator simulator(r_inner, ntheta,nr, include_composition);
Axis axis(simulator);
Core core(simulator);
ModeButton modebutton(simulator);
Seismograph seismograph(simulator);

//Structures for initializing a window and OpenGL conext
SDL_GLContext context;
SDL_Window *window=NULL;

//Given an x and y location, compute the r, theta location of the 
//simulator domain. Can handle the case of ellipticity.
//x and y are assumed to range from -0.5 to 0.5, not -1 to -1.
//Why? I'm not sure and it is lost to the sands of time.
inline void compute_simulator_location( const float x, const float y, float *theta, float *r )
{
  //Scale the click coordinates so that they are
  //in the correct place if the domain is elliptical
  float xpp = x, ypp = y;
  if (flattening > 1.e-2)
  {
    //Rotate the click position so that the ellipse
    //corresponds to the Cartesian axes.
    float rot_angle = simulator.spin_angle() - M_PI/2.;
    float xp = x * std::cos(-rot_angle) - y * std::sin(-rot_angle);
    float yp = x * std::sin(-rot_angle) + y * std::cos(-rot_angle);

    //Do the inverse flattening
    yp = yp / (1.-flattening);
    xp = xp;

    //Rotate back to the initial angle
    xpp = xp * std::cos(rot_angle) - yp * std::sin(rot_angle);
    ypp = xp * std::sin(rot_angle) + yp * std::cos(rot_angle);
  }
  *theta = std::atan2( ypp, xpp );
  *theta = (*theta < 0. ? *theta + 2.*M_PI : *theta );
  *r = 2.*std::sqrt( xpp*xpp + ypp*ypp );
}

void toggle_seismic_mode()
{
  seismic_mode = !seismic_mode;
  simulator.clear_seismic_waves();
  seismograph.clear_record();
  if (seismic_mode)
  {
    colormap = &seismic;
  }
  else
  {
    colormap = &hot;
  }
}

inline bool check_buttons(float x, float y)
{
  // Handle the weird choice of -0.5 to 0.5
  float xp = x*2.0f;
  float yp = y*2.0f;
  if (std::sqrt((xp-mode_button_x)*(xp-mode_button_x) + (yp-mode_button_y)*(yp-mode_button_y)) < mode_button_radius)
  {
    toggle_seismic_mode();
    return true;
  }
  return false;
}

//Given x,y from -0.5 to 0.5, from a mouse button
//or finger event, handle it in the simulator.
inline void handle_mouse_or_finger_motion(float x, float y)
{
  float theta, r;

  compute_simulator_location( x, y, &theta, &r);

  click_theta = ltheta * theta / 2. / M_PI;
  click_r = lr*(r-r_inner)/(1.0f-r_inner);
}

//Update where to add heat
inline void handle_mouse_motion(SDL_MouseMotionEvent *event)
{
  // We handle toucn events separately.
  if (event->which == SDL_TOUCH_MOUSEID) {
    return;
  }
  const float x = float(event->x)/float(xpix)-0.5f;
  const float y = 1.0f-float(event->y)/float(ypix)-0.5f;
  handle_mouse_or_finger_motion(x, y);
}

//Update where to add heat
inline void handle_finger_motion(SDL_TouchFingerEvent *event)
{
  float x = event->x - 0.5f;
  float y = 0.5f - event->y;
  handle_mouse_or_finger_motion(x, y);
}

//Change the Rayleigh number on scrolling
inline void handle_mouse_wheel(SDL_MouseWheelEvent *event)
{
  double rayleigh = simulator.rayleigh_number();
  double factor = std::pow(10.0, 1./100.* (event->y < 0 ? -1. : 1.0) );
  simulator.update_state( rayleigh * factor );
}

//Toggle whether to add heat, and whether it should
//be positive or negative
inline void handle_mouse_button(SDL_MouseButtonEvent *event)
{
  // We handle toucn events separately.
  if (event->which == SDL_TOUCH_MOUSEID) {
    return;
  }
  if(event->state==SDL_PRESSED)
  {
    if(event->button == SDL_BUTTON_LEFT)
      click_state = 1;
    if(event->button == SDL_BUTTON_RIGHT)
      click_state = -1;

    float x = float(event->x)/float(xpix)-0.5f;
    float y = 1.0f-float(event->y)/float(ypix)-0.5f;

    check_buttons(x, y);

    float theta, r;
    compute_simulator_location(x, y, &theta, &r);

    click_theta = theta;
    click_r = r-r_inner;
  }
  else
  {
    click_state=0;
  }
}

//Toggle whether to add heat, and whether it should
//be positive or negative
inline void handle_finger_down(SDL_TouchFingerEvent *event)
{
  if(event->type==SDL_FINGERDOWN)
  {
    click_state = 1;
    float x = event->x - 0.5f;
    float y = 0.5f - event->y;

    check_buttons(x, y);
    float theta, r;
    compute_simulator_location(x, y, &theta, &r);

    click_theta = theta;
    click_r = r-r_inner;
  }
  else
  {
    click_state=0;
  }
}

inline bool in_domain( const float theta, const float r )
{
  return (r+r_inner < 1.) && ( r > 0.);
}

void cycle_colorscale()
{
  const unsigned int n_colormaps = 5;
  static color (*maps[n_colormaps])(double) =
    { &hot, &viridis, &inferno, &seismic, &gist_earth };
  static unsigned short i = 0;
  i = (i+1)%n_colormaps;
  colormap = maps[i];
}

void text_output()
{
  //rough timescale conversion for Earth's mantle:
  //L^2/kappa / seconds/Gyr
  const double timescale = 2800.e3*2800.e3/1.e-6 /M_PI/1.e7/1.e9;
  //Unicode characters for exponentiation
  const char* exponents[10] = { "\u2070", "\u00B9", "\u00B2",
                                "\u00b3", "\u2074", "\u2075",
                                "\u2076", "\u2077", "\u2078",
                                "\u2079"};
  //Output Rayleigh information
  const double Ra = simulator.rayleigh_number();
  if( Ra < 1000 )
    std::cout<<"Ra: "<<std::setprecision(3)<<Ra<<"\t";
  else
  {
    const int upper = int(std::floor(std::log10(Ra)));
    const int lower = int(rint(Ra/std::pow(10, upper)));
    if (lower == 1)
      std::cout<<"Ra: "<<"  10"<<exponents[upper]<<"\t";
    else
      std::cout<<"Ra: "<<lower<<"x10"<<exponents[upper]<<"\t";
  }

  //Output Nusselt information
  std::cout<<"Nu: "<<std::setprecision(2)<<simulator.nusselt_number()<<"\t";

  //Output time information
  if (simulation_time*timescale < 1.)
    std::cout<<"\tTime: "<<int(rint(simulation_time*timescale*1000.))<<" million years\t";
  else
    std::cout<<"\tTime: "<<std::setprecision(2)<<simulation_time*timescale<<" billion years\t";

  //clear the line
  std::cout<<std::endl;
}

//Actually perform the timestep
void timestep()
{
  static int i=0;  //Keep track of timestep number
  simulator.draw( include_composition && draw_composition );  //Draw to screen
  core.draw();
  if (include_tpw && !seismic_mode)
    axis.draw();
  if (seismic_mode)
    seismograph.draw();
  modebutton.draw();

  //Do the convection problem if not in seismic mode
  if( !seismic_mode )
  {
    //At the moment, the stokes solve is not the limiting factor,
    //so it does not hurt to do it every timestep
    simulator.solve_stokes();

    //Add heat if the user is clicking
    if(click_state != 0 && ( (include_composition && !draw_composition)||(!include_composition)) && in_domain(click_theta, click_r) )
      simulator.add_heat(click_theta, click_r, (click_state==1 ? true : false));
    if(click_state != 0 && include_composition && draw_composition && in_domain(click_theta, click_r) )
      simulator.add_composition(click_theta, click_r);

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

#ifndef __EMSCRIPTEN__
    //Output scaling information
    text_output();
#endif

    //increment timestep
    ++i;
    simulation_time += simulator.timestep();
  }
  //If we are in seismic mode
  else
  {
    //Make earthquakes
    if(click_state == 1 && in_domain(click_theta,click_r) )
      simulator.earthquake(click_theta, click_r);
    else if(click_state == -1 && in_domain(click_theta,click_r) )
      simulator.place_seismometer(click_theta, click_r);
    //Propagate waves
    simulator.propagate_seismic_waves();
  }
}


//Initialize the windows, set up all the OpenGL stuff
void init()
{
    SDL_Init(SDL_INIT_VIDEO);
    SDL_SetHint(SDL_HINT_MAC_CTRL_CLICK_EMULATE_RIGHT_CLICK, "1");

#ifndef __EMSCRIPTEN__
    //Get the number of displays
    unsigned int n_displays = SDL_GetNumVideoDisplays();
    if (n_displays < 1)
       std::cerr<<"SDL : Could not find displays"<<std::endl;

    //Identify the largest display
    int max_display_dimension = 0;
    SDL_Rect draw_rect;
    float screen_fraction = 0.8;
    for (unsigned int i=0; i < n_displays; ++i)
    {
      SDL_Rect r;
      SDL_GetDisplayBounds( i, &r );
      if ( r.w > max_display_dimension )
      {
        max_display_dimension = r.w;
        draw_rect = r;
      }
      if ( r.h > max_display_dimension )
      {
        max_display_dimension = r.h;
        draw_rect = r;
      }
    }
    //set the size of the window
    xpix = int( screen_fraction * float( (draw_rect.w < draw_rect.h ? draw_rect.w : draw_rect.h )) );
    ypix = xpix;

    window = SDL_CreateWindow(
       "Convection",
        draw_rect.x + draw_rect.w/2 - xpix/2 , draw_rect.y + draw_rect.h/2 - ypix/2,
        xpix, ypix,
        SDL_WINDOW_OPENGL);
#else
    int width, height;
    emscripten_get_canvas_element_size("#canvas", &width, &height);
    xpix = (width < height ? width : height);
    ypix = xpix;

    window = SDL_CreateWindow(
       "Convection",
        SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED,
        xpix, ypix,
        SDL_WINDOW_OPENGL);
#endif

    context = SDL_GL_CreateContext(window);
    if (!context)
        SDL_LogError(SDL_LOG_CATEGORY_APPLICATION, "SDL_GL_CreateContext(): %s\n", SDL_GetError());

    glewExperimental=GL_TRUE;
    GLenum glew_status = glewInit();
    if (glew_status != GLEW_OK) {
      fprintf(stderr, "Error: %s\n", glewGetErrorString(glew_status));
    }


    simulator.setup_opengl();
    core.setup();
    if (include_tpw)
      axis.setup();
    seismograph.setup();
    modebutton.setup();
}

//Cleanup
void quit()
{
    simulator.cleanup_opengl();
    core.cleanup();
    if(include_tpw)
      axis.cleanup();
    seismograph.cleanup();
    modebutton.cleanup();

    SDL_GL_DeleteContext(context);
    SDL_Quit();
    exit(0);
}

void loop()
{
  SDL_Event event;
  while(SDL_PollEvent(&event))
  {
    switch(event.type)
    {
      case SDL_QUIT:
        quit();
      case SDL_KEYDOWN:
        //Some systems seem to have a hair trigger on the
        //keypress repeats. Filter those out.
        if(!event.key.repeat)
        {
          switch(event.key.keysym.sym)
          {
            case SDLK_SPACE:
              toggle_seismic_mode();
              break;

#ifndef __EMSCRIPTEN__
            //Quitting should be handled by navigating to another
            //webpage or closing the browser when using, emscripten,
            //so just disable quitting on escape
            case SDLK_ESCAPE:
              quit();
              break;
#endif
            case SDLK_TAB:
              draw_composition = !draw_composition;
              break;

            case SDLK_BACKSPACE:
              advection_diffusion = !advection_diffusion;
              break;

            case SDLK_c:
              cycle_colorscale();
              break;

            default:
              break;
          }
        }
        break;
      case SDL_MOUSEBUTTONDOWN:
      case SDL_MOUSEBUTTONUP:
        handle_mouse_button(&event.button);
        break;
      case SDL_MOUSEMOTION:
        handle_mouse_motion(&event.motion);
        break;
      case SDL_MOUSEWHEEL:
        handle_mouse_wheel(&event.wheel);
        break;
      case SDL_FINGERDOWN:
      case SDL_FINGERUP:
        handle_finger_down(&event.tfinger);
        break;
      case SDL_FINGERMOTION:
        handle_finger_motion(&event.tfinger);
        break;
      default:
        break;
    }
  }
  timestep();
  SDL_GL_SwapWindow(window);
}

int main(int argc, char** argv)
{
    init();

    //The browser needs to handle the event loop if we are using
    //Emscripten, so in that case we need to give the event loop
    //away.  Otherwise we do it.
#ifdef __EMSCRIPTEN__
    emscripten_set_main_loop(loop, 0, 1);
#else
    while (true) {
        loop();
    }
#endif

    quit();
    return 0;
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
double emscripten_nusselt()
{
  return simulator.nusselt_number();
}
EMSCRIPTEN_BINDINGS(my_module) {
  emscripten::function("rayleigh", &emscripten_rayleigh );
  emscripten::function("nusselt", &emscripten_nusselt );
}
#endif
