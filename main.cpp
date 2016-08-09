#include "GL/glew.h"
#include "SDL2/SDL_opengl.h"
#include "SDL2/SDL.h"
#include "convection.h"
#include "rendering_plugins.h"
#include <cmath>
#include <iostream>
#include <iomanip>

#ifdef __EMSCRIPTEN__
#include <emscripten/emscripten.h>
#endif

/*********************************************
    MODIFY THSE PARAMETERS TO CHANGE THE
    BEHAVIOR OF THE SIMULATION.
*********************************************/

//Number of cells in the theta and r directions.
//This is the primary control on resolution,
//as well as performance.
const unsigned int ntheta = 1024;
const unsigned int nr = 128;

//Aspect ratio of the computational domain
//is set by the inner radius, where the outer
//radius is assumed to be 1.0
const double r_inner = 0.35;

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

//Whether to draw composition or temperature fields
bool draw_vorticity = false;

//Global solver
ConvectionSimulator simulator(r_inner, ntheta,nr);
Core core(simulator);

//Structures for initializing a window and OpenGL conext
SDL_GLContext context;
SDL_Window *window=NULL;

//Given an x and y location, compute the r, theta location of the 
//simulator domain. Can handle the case of ellipticity.
inline void compute_simulator_location( const float x, const float y, float *theta, float *r )
{
  //Scale the click coordinates so that they are
  //in the correct place if the domain is elliptical
  float xpp = x, ypp = y;
  *theta = std::atan2( ypp, xpp );
  *theta = (*theta < 0. ? *theta + 2.*M_PI : *theta );
  *r = 2.*std::sqrt( xpp*xpp + ypp*ypp );
}

//Update where to add heat
inline void handle_mouse_motion(SDL_MouseMotionEvent *event)
{
  float x = float(event->x)/float(xpix)-0.5f;
  float y = 1.0f-float(event->y)/float(ypix)-0.5f;

  float theta, r;
  compute_simulator_location( x, y, &theta, &r);

  click_theta = ltheta * theta / 2. / M_PI;
  click_r = lr*(r-r_inner)/(1.0f-r_inner);
}

//Change the Rayleigh number on scrolling
inline void handle_mouse_wheel(SDL_MouseWheelEvent *event)
{
  double rayleigh = simulator.rayleigh_number();
  double factor = std::pow(10.0, 1./100.* (event->y < 0 ? -1. : 1.0) );
  //simulator.update_state( rayleigh * factor, 1000., 1.0 );
}

//Toggle whether to add heat, and whether it should
//be positive or negative
inline void handle_mouse_button(SDL_MouseButtonEvent *event)
{
  if(event->state==SDL_PRESSED)
  {
    if(event->button == SDL_BUTTON_LEFT)
      click_state = 1;
    if(event->button == SDL_BUTTON_RIGHT)
      click_state = -1;

    float x = float(event->x)/float(xpix)-0.5f;
    float y = 1.0f-float(event->y)/float(ypix)-0.5f;

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

//Actually perform the timestep
void timestep()
{
  static int i=0;  //Keep track of timestep number
  simulator.draw( draw_vorticity );  //Draw to screen
  core.draw();

  simulator.solve_poisson();

  //Add heat if the user is clicking
  if(click_state != 0 && !draw_vorticity && in_domain(click_theta, click_r) )
    simulator.add_heat(click_theta, click_r, (click_state==1 ? true : false));
  if(click_state != 0 && draw_vorticity && in_domain(click_theta, click_r) )
    simulator.add_vorticity(click_theta, click_r, (click_state==1 ? true : false));
  simulator.generate_vorticity();
  //Advect temperature and composition fields
  simulator.semi_lagrangian_advect_temperature();
  simulator.semi_lagrangian_advect_vorticity();

  //Diffuse temperature and vorticity
  simulator.diffuse_temperature();
  simulator.diffuse_vorticity();

  //increment timestep
  ++i;
  std::cout<<"Ra: "<<simulator.rayleigh_number()<<std::endl;
}


//Initialize the windows, set up all the OpenGL stuff
void init()
{
    SDL_Init(SDL_INIT_VIDEO);
    SDL_SetHint(SDL_HINT_MAC_CTRL_CLICK_EMULATE_RIGHT_CLICK, "1");

#ifndef __EMSCRIPTEN__
    //Get the number of displays
    int n_displays = SDL_GetNumVideoDisplays();
    if (n_displays < 1)
       std::cerr<<"SDL : Could not find displays"<<std::endl;

    //Identify the largest display
    unsigned int max_display_dimension = 0;
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
    int width, height, isFullscreen;
    emscripten_get_canvas_size( &width, &height, &isFullscreen );
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
}

//Cleanup
void quit()
{
    simulator.cleanup_opengl();
    core.cleanup();

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
        if(event.key.keysym.sym == SDLK_TAB)
          draw_vorticity = ! draw_vorticity;
        //Quitting should be handled by navigating to another
        //webpage or closing the browser when using, emscripten,
        //so just disable quitting on escape
#ifndef __EMSCRIPTEN__
        else if(event.key.keysym.sym == SDLK_ESCAPE)
          quit();
#endif
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
