#include <SDL3/SDL.h>
#include <SDL3/SDL_main.h>

#include <array>
#include <chrono>
#include <cstdio>

const int grid_x = 150;
const int grid_y = 150;
const int pix_w = 10;
const int pix_h = 10;

const int frame_ms = 1000 / 60;

int main() {
  printf("\n -- running flip fluid simulator -- \n\n");

  std::array<std::array<bool, 10>, 10> grid{{
      {0, 1, 0, 1, 0, 1, 0, 1, 0, 1},
      {1, 0, 1, 0, 1, 0, 1, 0, 1, 0},
      {0, 1, 0, 1, 0, 1, 0, 1, 0, 1},
      {1, 0, 1, 0, 1, 0, 1, 0, 1, 0},
      {0, 1, 0, 1, 0, 1, 0, 1, 0, 1},
      {1, 0, 1, 0, 1, 0, 1, 0, 1, 0},
      {0, 1, 0, 1, 0, 1, 0, 1, 0, 1},
      {1, 0, 1, 0, 1, 0, 1, 0, 1, 0},
      {0, 1, 0, 1, 0, 1, 0, 1, 0, 1},
      {1, 0, 1, 0, 1, 0, 1, 0, 1, 0},
  }};

  SDL_Init(SDL_INIT_VIDEO);

  SDL_Window *window = SDL_CreateWindow("ffs", 400, 400, 0);
  SDL_Renderer *renderer = SDL_CreateRenderer(window, nullptr);

  bool running = true;
  long long cycle = 0;
  long long render_duration = 0;

  while (running) {
    auto start = std::chrono::high_resolution_clock::now();

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

    if (cycle % 15 == 1) {
      for (int i = 0; i < grid.size(); ++i) {
        for (int j = 0; j < grid[i].size(); ++j) {
          grid[i][j] = !grid[i][j];
        }
      }
    }

    auto render_start = std::chrono::high_resolution_clock::now();

    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
    SDL_RenderFillRect(renderer, nullptr);

    SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
    for (int i = 0; i < grid.size(); ++i) {
      for (int j = 0; j < grid[i].size(); ++j) {
        if (grid[i][j]) {
          SDL_FRect rect{(float)grid_x + pix_w * j, (float)grid_y + pix_h * i,
                         pix_w, pix_h};
          SDL_RenderFillRect(renderer, &rect);
        }
      }
    }

    SDL_RenderPresent(renderer);

    long render_time =
        std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::high_resolution_clock::now() - render_start)
            .count();
    render_duration += render_time;

    long time = std::chrono::duration_cast<std::chrono::milliseconds>(
                    std::chrono::high_resolution_clock::now() - start)
                    .count();
    SDL_Delay(std::max(frame_ms - time, 0L));
    ++cycle;
  }
  printf("avg render time: %lld μs / frame\n", render_duration / cycle);

  SDL_Quit();
  return 0;
}