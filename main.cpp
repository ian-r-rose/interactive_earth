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

const unsigned int nx = 200;
const unsigned int ny = 50;
const double lx = 4.0;
const double ly = 1.0;
const double Ra = 1.e7;
const unsigned int scale = 8;

const unsigned int xpix = nx*scale;
const unsigned int ypix = ny*scale;

int heat_state = 0;
double hx, hy;
bool solve_stokes = true;

StokesSolver* handle = NULL;

SDL_GLContext context;
SDL_Window *window=NULL;

inline void handle_mouse_motion(SDL_MouseMotionEvent *event)
{
  hx = lx*(double(event->x)/double(xpix));
  hy = ly*(1.0-double(event->y)/double(ypix));
}

inline void handle_mouse_wheel(SDL_MouseWheelEvent *event)
{
  double rayleigh = handle->rayleigh_number();
  double factor = std::pow(10.0, 1./20.* (event->y < 0 ? -1. : 1.0) );
  handle->update_state( rayleigh * factor );
}
  

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

void timestep()
{
  static int i=0;
  if(i%1==0)
    handle->draw();
  if(solve_stokes && i%2==0)
    handle->solve_stokes();

  if(heat_state != 0) handle->add_heat(hx, hy, (heat_state==1 ? true : false));
  handle->semi_lagrangian_advect();
  handle->diffuse_temperature();
  std::cout<<"Ra: "<<std::setprecision(3)<<handle->rayleigh_number()
           <<"\tNu: "<<handle->nusselt_number()<<std::endl;
  ++i;
}


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
        if(event.key.keysym.sym == SDLK_ESCAPE)
          quit();
        else if(event.key.keysym.sym == SDLK_SPACE)
          solve_stokes = !solve_stokes;
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
  SDL_RenderPresent(renderer);
}

int main(int argc, char** argv)
{
    StokesSolver stokes(lx, ly, nx,ny, Ra);
    handle = &stokes;
 
    init();

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
