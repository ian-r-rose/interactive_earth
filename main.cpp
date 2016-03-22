#include "GL/glew.h"
#include "SDL2/SDL_opengl.h"
#include "SDL2/SDL.h"
#include "convection.h"
#include <cmath>
#include <iostream>
#include <iomanip>

#ifdef __EMSCRIPTEN__
#include <emscripten/emscripten.h>
#endif

//Number of cells in the x and y directions.
//This is the primary control on resolution,
//as well as performance/
const unsigned int nx = 1024;
const unsigned int ny = 128;

//Size of computational domain.  The ratio of
//lx to ly should be the same of nx to ny,
//otherwise the convective features will
//look kind of squashed and funny.
const double r_inner = 0.5;
const double lx = 2.*M_PI;
const double ly = 1.0-r_inner;

//Initial Rayleigh number of simulation
const double Ra = 1.e7;

//Factor for how much to blow up the rendered
//triangles so that they are bigger on screen
const unsigned int scale = 1;

//Total number of pixels in x and y directions
const unsigned int xpix = nx*scale;
const unsigned int ypix = nx*scale;

//Whether to add heat to the simulation on mouse click.
int click_state = 0;
//Location of heat-adding
double hx, hy;

//Whether to solve the stokes equation for a given timestep.
bool solve_stokes = true;

//Whether we are in earthquake mode
bool seismic_mode = false;

//Global solver
ConvectionSimulator simulator(r_inner, nx,ny, Ra);

//Structures for initializing a window and OpenGL conext
SDL_GLContext context;
SDL_Window *window=NULL;

//Update where to add heat
inline void handle_mouse_motion(SDL_MouseMotionEvent *event)
{
  const float xx = float(event->x)/float(xpix);
  const float yy = 1.0f-float(event->y)/float(ypix);
  float theta = std::atan2( yy-0.5f, xx-0.5f );
  theta = (theta < 0. ? theta + 2.*M_PI : theta );
  const float r = 2.*std::sqrt( (xx-0.5f)*(xx-0.5f) + (yy-0.5f)*(yy-0.5f) );

  hx = lx * theta / 2. / M_PI;
  hy = ly*(r-r_inner)/(1.0f-r_inner);
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
  if(event->state==SDL_PRESSED)
  {
    if(event->button == SDL_BUTTON_LEFT)
      click_state = 1;
    if(event->button == SDL_BUTTON_RIGHT)
      click_state = -1;

    const float xx = float(event->x)/float(xpix);
    const float yy = 1.0f-float(event->y)/float(ypix);
    float theta = std::atan2( yy-0.5f, xx-0.5f );
    theta = (theta < 0. ? theta + 2.*M_PI : theta );
    const float r = 2.*std::sqrt( (xx-0.5f)*(xx-0.5f) + (yy-0.5f)*(yy-0.5f) );

    hx = lx * theta / 2. / M_PI;
    hy = ly*(r-r_inner)/(1.0f-r_inner);

  }
  else
  {
    click_state=0;
  }
}

inline bool in_domain( const float x, const float y )
{
  const float xx = float(x)/float(xpix);
  const float yy = 1.0f-float(y)/float(ypix);
  const float r = 2.*std::sqrt( (xx-0.5f)*(xx-0.5f) + (yy-0.5f)*(yy-0.5f) );
  //return (r < 1.) && ( r > r_inner);
  return true;
}

//Actually perform the timestep
void timestep()
{
  static int i=0;  //Keep track of timestep number
  simulator.draw();  //Draw to screen

  //Do the convection problem if not in seismic mode
  if( !seismic_mode )
  {
    //I have found that it usually looks okay if we only update the velocity solution
    //every other timestep.  Any more delay and it starts to look funny.
    if(solve_stokes && i%2==0)
      simulator.solve_stokes();

    //Add heat if the user is clicking
    if(click_state != 0 && in_domain(hx, hy) ) simulator.add_heat(hx, hy, (click_state==1 ? true : false));

    //Advect temperature field
    simulator.semi_lagrangian_advect();

    //Diffuse temperature
    simulator.diffuse_temperature();

    //Output scaling information
    std::cout<<"Ra: "<<std::setprecision(3)<<simulator.rayleigh_number()
             <<"\tNu: "<<simulator.nusselt_number()<<std::endl;
    //increment timestep
    ++i;
  }
  //If we are in seismic mode
  else
  {
    //Make earthquakes
    if(click_state != 0 && in_domain(hx,hy) ) simulator.earthquake(hx, hy);
    //Propagate waves
    simulator.propagate_seismic_waves();
  }
}


//Initialize the windows, set up all the OpenGL stuff
void init()
{
    SDL_Init(SDL_INIT_VIDEO);


    window = SDL_CreateWindow(
       "Convection",
        SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED,
        xpix, ypix,
        SDL_WINDOW_OPENGL);
    context = SDL_GL_CreateContext(window);
    if (!context)
        SDL_LogError(SDL_LOG_CATEGORY_APPLICATION, "SDL_GL_CreateContext(): %s\n", SDL_GetError());

    glewExperimental=GL_TRUE;
    GLenum glew_status = glewInit();
    if (glew_status != GLEW_OK) {
      fprintf(stderr, "Error: %s\n", glewGetErrorString(glew_status));
    }


    simulator.setup_opengl();
}

//Cleanup
void quit()
{
    simulator.cleanup_opengl();
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
        if(event.key.keysym.sym == SDLK_SPACE)
        {
          seismic_mode = !seismic_mode;
          simulator.clear_seismic_waves();
        }
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
