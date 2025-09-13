#include <SDL3/SDL.h>
#include <SDL3/SDL_main.h>

#include <array>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdio>

const int grid_x = 100;
const int grid_y = 100;
const int pix_w = 5;
const int pix_h = 5;

const int frame_ms = 1000 / 60;

const int n = 20;                          // rows
const int m = 40;                           // columns
const float spacing = 10.0f;               // grid spacing
const float timestep = frame_ms / 1000.0f; // seconds

const float relaxation =
    1.8f; // from 1 to 2, scales incompressibility calculation

std::array<std::array<int, m>, n> states{};
std::array<std::array<float, m + 1>, 2 * n + 1> vels{};
std::array<std::array<int, m>, n> avels{};

int cell_state(int i, int j) {
  if (0 <= i && i < n && 0 <= j && j < m) {
    return states[i][j];
  } else {
    return 0;
  };
};

void clamp(int &value, int low, int hi) {
  if (value < low) {
    value = low;
  } else if (value > hi) {
    value = hi;
  }
}

float vel(int i, int j) {
  clamp(i, 0, vels.size() - 1);
  clamp(j, 0, m + i % 2 - 1);
  return vels[i][j];

  // skip the lower
  if (0 <= i && i < vels.size() && 0 <= j && j < m + i % 2) {
    return vels[i][j];
  } else {
    // printf("%d %d\n", i, j);
    return 0.0f;
  };
};

void print_vels() {
  printf("---\n");
  for (int i = 0; i < vels.size(); i += 1) {
    printf("|");
    if (i % 2 == 0) {
      printf("     ");
    }
    for (int j = 0; j < m + i % 2; ++j) {
      // vels[i][j] += 9.0f; // gravity
      if (i % 2 && j + 1 < m + i % 2) {
        printf("%+3.1f  .  ", vels[i][j]);
      } else {
        printf("%+3.1f     ", vels[i][j]);
      }
    }
    printf("\n");
  }
};

int main() {
  printf("\n -- running flip fluid simulator -- \n\n");

  // initial state
  for (auto &row : states) {
    for (auto &entry : row) {
      entry = 1;
    }
  }
  for (auto &row : vels) {
    for (auto &entry : row) {
      entry = 0.0f;
    }
  }

  // // left wall fan
  // for (int i = 0; i < n; ++i) {
  //   states[i][0] = 0;
  // }

  for (int i = 8; i < 16; ++i) {
    vels[2 * i + 1][0] = 4;
  }

  // rules
  // vertical (y, even indices) elements are m + 1 long
  // horizontal (x, odd indices) elements are m long
  // extra element should be 0

  SDL_Init(SDL_INIT_VIDEO);

  SDL_Window *window = SDL_CreateWindow("ffs", 400, 400, 0);
  SDL_Renderer *renderer = SDL_CreateRenderer(window, nullptr);

  bool running = true;
  long long cycle = 0;
  long long render_duration = 0;

  float max_v = 0;

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

    // simulation

    // 1. external forces (every x velocity)
    for (int i = 2; i < vels.size() - 2; i += 2) {
      for (int j = 0; j < m; ++j) {
        vels[i][j] += 9.0f; // gravity
      }
    }

    // 2. incompressibility (every cell)
    for (int i = 0; i < 100; ++i) {
      for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
          if (states[i][j] == 0) {
            continue;
          }

          float &c = vels[2 * i + 1][j + 1]; // right flow (out)
          float &a = vels[2 * i + 1][j];     // left flow (in)
          float &d = vels[2 * i + 2][j];     // bottom flow (out)
          float &b = vels[2 * i][j];         // top flow (in)
          const float net_outflow = c - a + d - b;

          const int sc = cell_state(i, j + 1); // right state
          const int sa = cell_state(i, j - 1); // left state
          const int sd = cell_state(i + 1, j); // bottom state
          const int sb = cell_state(i - 1, j); // top state
          const int s = sc + sa + sd + sb;
          if (s == 0) {
            continue;
          }

          // printf("%f %.0f\n", net_outflow, s);

          c -= net_outflow * relaxation * (sc / (float)s);
          a += net_outflow * relaxation * (sa / (float)s);
          d -= net_outflow * relaxation * (sd / (float)s);
          b += net_outflow * relaxation * (sb / (float)s);
        }
      }
    }

    // print_vels();

    auto new_vels = vels;
    // 3. semi-lagrangian advection UNCHECKED
    for (int i = 0; i < vels.size(); ++i) {
      for (int j = 0; j < m + i % 2; ++j) {
        float vf_x = 0;
        float vf_y = 0;

        bool horizontal = i % 2;

        if (horizontal) {
          if (!cell_state(i / 2, j - 1) || !cell_state(i / 2, j)) {
            continue;
          }
          vf_x = vels[i][j];
          const float v1 = vel(i - 1, j);     // top right
          const float v2 = vel(i + 1, j);     // bottom right
          const float v3 = vel(i + 1, j - 1); // bottom left
          const float v4 = vel(i - 1, j - 1); // top left
          vf_y = (v1 + v2 + v3 + v4) / 4.0f;
        } else {
          if (!cell_state(i / 2 - 1, j) || !cell_state(i / 2, j)) {
            continue;
          }
          const float v1 = vel(i - 1, j + 1); // top right
          const float v2 = vel(i + 1, j + 1); // bottom right
          const float v3 = vel(i + 1, j);     // bottom left
          const float v4 = vel(i - 1, j);     // top left
          vf_x = (v1 + v2 + v3 + v4) / 4.0f;
          vf_y = vels[i][j];
        }

        const float xf_x = spacing * (j + (horizontal ? 0.0f : 0.5f));
        const float xf_y = spacing * i * 0.5f;

        // work backwards to check what (approximate) point would've been there
        const float xi_x = xf_x - timestep * vf_x;
        const float xi_y = xf_y - timestep * vf_y;

        // printf("(%d %d) vf: %.3f %.3f\n", i, j, vf_x, vf_y);
        // printf("(%d %d) xf: %.0f %.0f\n", i, j, xf_x, xf_y);
        // printf("(%d %d) xi: %.3f %.3f\n", i, j, xi_x, xi_y);

        // weighted average based on area opposite the known velocities
        auto interpolate_vel = [](int cell_i, int cell_j, float xi_x,
                                  float xi_y) -> float {
          const float cell_x = spacing * (cell_j + (cell_i % 2 ? 0.0f : 0.5f));
          const float cell_y = spacing * cell_i * 0.5f;
          // printf("(%d %d) cell: %.0f %.0f\n", cell_i, cell_j, cell_x,
          // cell_y);

          const float rel_x = (xi_x - cell_x) / spacing;
          const float rel_y = (xi_y - cell_y) / spacing;
          // printf("---> rel: %.4f %.4f\n", rel_x, rel_y);

          assert(rel_x >= -0.0001);
          assert(rel_y >= -0.0001);
          assert(rel_x < 1.0001);
          assert(rel_y < 1.0001);

          const float w1 = rel_x * (1.0f - rel_y);          // top right
          const float w2 = rel_x * rel_y;                   // bottom right
          const float w3 = (1.0f - rel_x) * rel_y;          // bottom left
          const float w4 = (1.0f - rel_x) * (1.0f - rel_y); // top left

          // ignore wall velocities? or set to zero? i choose to set to zero
          return w1 * vel(cell_i, cell_j + 1) +
                 w2 * vel(cell_i + 2, cell_j + 1) +
                 w3 * vel(cell_i + 2, cell_j) + w4 * vel(cell_i, cell_j);
        };

        if (horizontal) {
          const int cell_i = 2 * std::floor(xi_y / spacing - 0.5) + 1;
          const int cell_j = std::floor(xi_x / spacing);
          new_vels[i][j] = interpolate_vel(cell_i, cell_j, xi_x, xi_y);
        } else {
          const int cell_i = 2 * std::floor(xi_y / spacing);
          const int cell_j = std::floor(xi_x / spacing - 0.5);
          new_vels[i][j] = interpolate_vel(cell_i, cell_j, xi_x, xi_y);
        }
      }
    }

    vels = new_vels;

    // drawing
    auto render_start = std::chrono::high_resolution_clock::now();

    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
    SDL_RenderFillRect(renderer, nullptr);

    for (int i = 0; i < states.size(); ++i) {
      for (int j = 0; j < states[i].size(); ++j) {
        if (states[i][j]) {
          // vel_at, scale hue under max_vel,
          uint8_t col =
              std::min(255, (int)((10 / 255.0) * std::abs(vels[2 * i + 1][j])));
          SDL_SetRenderDrawColor(renderer, 100, col, 155, 255);
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