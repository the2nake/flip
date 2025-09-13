#include <SDL3/SDL_events.h>
#include <SDL3/SDL_hints.h>
#include <SDL3/SDL_oldnames.h>
#include <SDL3/SDL_render.h>
#include <SDL3/SDL_video.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <time.h>

#include <SDL3/SDL.h>
#include <SDL3/SDL_main.h>

#define WINDOW_W 400
#define WINDOW_H 400

int main() {
  printf("\n -- running flip fluid simulator -- \n\n");

  SDL_Init(SDL_INIT_VIDEO);

  SDL_Window *window = SDL_CreateWindow("ffs", WINDOW_W, WINDOW_H,
                                        SDL_WINDOW_HIGH_PIXEL_DENSITY);
  SDL_Renderer *renderer = SDL_CreateRenderer(window, NULL);

  printf("on wayland, try SDL_VIDEODRIVER=wayland\n\n");
  
  printf("content dpi scaling: %.2f\n",
         SDL_GetDisplayContentScale(SDL_GetDisplayForWindow(window)));
  printf("window dpi scaling: %.2f\n", SDL_GetWindowDisplayScale(window));
  printf("window pixel density: %.2f\n", SDL_GetWindowPixelDensity(window));

  bool running = true;

  double update_time;
  long long cycles = 0;

  while (running) {
    double update0 = clock() / (double)CLOCKS_PER_SEC;

    SDL_Event event;
    while (SDL_PollEvent(&event)) {
      switch (event.type) {
      case SDL_EVENT_QUIT:
        running = false;
        break;
      default:
        break;
      }
    }

    double update1 = clock() / (double)CLOCKS_PER_SEC;
    update_time += update1 - update0;

    SDL_RenderPresent(renderer);
    SDL_Delay(20);
    ++cycles;
  }

  printf("\n -- finished -- \n\n");

  printf("time per cycle: %f ms\n\n", 1000.0 * update_time / cycles);

  return 0;
}