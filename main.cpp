#include "GL/glew.h"
#include "SDL2/SDL_opengl.h"
#include "SDL2/SDL.h"
#include "stokes.h"
#include <cmath>
#include <iostream>
#include <iomanip>

#ifdef __EMSCRIPTEN__
#include <emscripten/emscripten.h>
#endif

//Number of cells in the x and y directions.
//This is the primary control on resolution,
//as well as performance/
const unsigned int nx = 400;
const unsigned int ny = 100;

//Size of computational domain.  The ratio of
//lx to ly should be the same of nx to ny,
//otherwise the convective features will
//look kind of squashed and funny.
const double lx = 4.0;
const double ly = 1.0;

//Initial Rayleigh number of simulation
const double Ra = 1.e7;

//Factor for how much to blow up the rendered 
//triangles so that they are bigger on screen
const unsigned int scale = 3;

//Total number of pixels in x and y directions
const unsigned int xpix = nx*scale;
const unsigned int ypix = ny*scale;

//Whether to add heat to the simulation on mouse click.
int heat_state = 0;
//Location of heat-adding
double hx, hy;

//Whether to solve the stokes equation for a given timestep.
bool solve_stokes = true;

//Pointer for the solver so that the various event handlers
//can access it.  Did it this way due to the way GLUT is 
//organized, not really necessary now that I am using SDL.
StokesSolver* handle = NULL;

//Structures for initializing a window and OpenGL conext
SDL_GLContext context;
SDL_Window *window=NULL;

//Update where to add heat
inline void handle_mouse_motion(SDL_MouseMotionEvent *event)
{
  hx = lx*(double(event->x)/double(xpix));
  hy = ly*(1.0-double(event->y)/double(ypix));
}

//Change the Rayleigh number on scrolling
inline void handle_mouse_wheel(SDL_MouseWheelEvent *event)
{
  double rayleigh = handle->rayleigh_number();
  double factor = std::pow(10.0, 1./20.* (event->y < 0 ? -1. : 1.0) );
  handle->update_state( rayleigh * factor );
}
  
//Toggle whether to add heat, and whether it should
//be positive or negative
inline void handle_mouse_button(SDL_MouseButtonEvent *event)
{
  if(event->state==SDL_PRESSED)
  {
     if(event->button == SDL_BUTTON_LEFT)
       heat_state = 1;
     if(event->button == SDL_BUTTON_RIGHT)
       heat_state = -1;
     hx = lx*(double(event->x)/double(xpix));
     hy = ly*(1.0-double(event->y)/double(ypix));
  }
  else heat_state=0;
}

//Actually perform the timestep
void timestep()
{
  static int i=0;  //Keep track of timestep number
  handle->draw();  //Draw to screen

  //The stokes solve is typically the most expensive part.  I have found
  //That it usually looks okay if we only update the velocity solution 
  //every other timestep.  Any more delay and it starts to look funny.
  if(solve_stokes && i%2==0)
    handle->solve_stokes();

  //Add heat if the user is clicking
  if(heat_state != 0) handle->add_heat(hx, hy, (heat_state==1 ? true : false));

  //Advect temperature field
  handle->semi_lagrangian_advect();

  //Diffuse temperature
  handle->diffuse_temperature();

  //Output scaling information
  std::cout<<"Ra: "<<std::setprecision(3)<<handle->rayleigh_number()
           <<"\tNu: "<<handle->nusselt_number()<<std::endl;
  //increment timestep
  ++i;
}


//Initialize the windows, set up all the OpenGL stuff
void init()
{
    SDL_Init(SDL_INIT_VIDEO);


    window = SDL_CreateWindow(
       "Convection", 10, 10, scale*nx, scale*ny, 
        SDL_WINDOW_OPENGL|SDL_WINDOW_RESIZABLE);
    context = SDL_GL_CreateContext(window);
    if (!context)
        SDL_LogError(SDL_LOG_CATEGORY_APPLICATION, "SDL_GL_CreateContext(): %s\n", SDL_GetError());

    glewExperimental=GL_TRUE;
    GLenum glew_status = glewInit();
    if (glew_status != GLEW_OK) {
      fprintf(stderr, "Error: %s\n", glewGetErrorString(glew_status));
    }


    handle->setup_opengl();
}

//Cleanup
void quit()
{
    handle->cleanup_opengl();
    SDL_GL_DeleteContext(context);
    SDL_Quit();
    exit(0);
}

void loop()
{
  SDL_Event event;
  SDL_Delay(1);
  while(SDL_PollEvent(&event))
  {
    switch(event.type)
    {
      case SDL_QUIT:
        quit();
      case SDL_KEYDOWN:
        //Neat feature of only advecting and diffusing, but not updating the
        //Stokes solution.  Not realistic, just fun to play with.
        if(event.key.keysym.sym == SDLK_SPACE)
          solve_stokes = !solve_stokes;
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
    StokesSolver stokes(lx, ly, nx,ny, Ra);
    handle = &stokes;
 
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
