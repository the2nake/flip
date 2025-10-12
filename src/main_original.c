#include <SDL3/SDL_events.h>
#include <SDL3/SDL_hints.h>
#include <SDL3/SDL_pixels.h>
#include <SDL3/SDL_render.h>
#include <SDL3/SDL_surface.h>
#include <SDL3/SDL_video.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <time.h>

#include <SDL3/SDL.h>
#include <SDL3/SDL_main.h>

#define SIM_W 288
#define SIM_H 102
#define CELL_W 4
#define CELL_H CELL_W
#define PARTICLE_COUNT 10000
#define BACKTRACK_PRECISION 4

#define WINDOW_W SIM_W *CELL_W
#define WINDOW_H SIM_H *CELL_H

const float frametime = 1.0 / 60.0;
const float k_timestep = 1.0 / 60.0;
const float k_gravity = 9.81;
const float k_relax = 1.9;
const int k_iters = 40;

float now() { return (float)clock() / CLOCKS_PER_SEC; }

typedef enum { WALL = 0, WATER = 1, AIR = 2 } state_e;
typedef struct {
  float x1;
  float x2;
  float v1;
  float v2;
} particle_s;

typedef struct {
  state_e s[SIM_H][SIM_W];       // state
  float v1[(SIM_W + 1) * SIM_H]; // horizontal velocity
  float v2[SIM_W * (SIM_H + 1)]; // vertical velocity
  int x1n;
  int x2n;

  particle_s particles[PARTICLE_COUNT];
} simulator_s;

float *x_vel(simulator_s *sim, int i) {
  if (0 <= i && i < sim->x1n) {
    return &sim->v1[i];
  } else {
    return NULL;
  }
}

float *y_vel(simulator_s *sim, int i) {
  if (0 <= i && i < sim->x2n) {
    return &sim->v2[i];
  } else {
    return NULL;
  }
}

bool in_bounds(int i, int j) {
  return 0 <= i && i < SIM_H && 0 <= j && j < SIM_W;
}

bool on_edge(int i, int j) {
  if (i == 0 || i == SIM_H - 1) {
    return true;
  }
  if (j == 0 || j == SIM_W - 1) {
    return true;
  }
  return false;
}

// use in_wall for particles
bool wall_at(simulator_s *sim, int i, int j) {
  if (in_bounds(i, j)) {
    return sim->s[i][j] == WALL;
  }
  return false;
}

// use wall_at for cells
bool in_wall(simulator_s *sim, particle_s particle) {
  int i = particle.x1 / CELL_W;
  int j = particle.x2 / CELL_H;
  return in_bounds(i, j) &&
         wall_at(sim, i, j); // check bounds to not lose particles
}

simulator_s initialise() {
  simulator_s sim;

  for (int i = 0; i < SIM_H; ++i) {
    for (int j = 0; j < SIM_W; ++j) {
      // if (i > SIM_H / 2) {
      sim.s[i][j] = WATER;
      // } else {
      //   sim.s[i][j] = AIR;
      // }
    }
  }

  for (int i = 0; i < SIM_H; ++i) {
    for (int j = 0; j < SIM_W; ++j) {
      if (on_edge(i, j)) {
        sim.s[i][j] = WALL;
      }
    }
  }

  sim.x1n = (SIM_W + 1) * SIM_H;
  sim.x2n = (SIM_H + 1) * SIM_W;

  for (int i = 0; i < sim.x1n; ++i) {
    sim.v1[i] = 0.f;
  }

  for (int i = 0; i < sim.x2n; ++i) {
    sim.v2[i] = 0.f;
  }

  // TODO: distribute particles inside water cells
  particle_s default_particle = {.x1 = 0, .x2 = 0, .v1 = 0, .v2 = 0};
  for (int i = 0; i < PARTICLE_COUNT; ++i) {
    sim.particles[i] = default_particle;
  }

  return sim;
}

void move_particles(simulator_s *sim) {
  // TODO: separate particles using LUT method, radix sorting?

  for (int i = 0; i < PARTICLE_COUNT; ++i) {
    sim->particles[i].v2 += k_gravity * k_timestep;
    float dx1 = sim->particles[i].v1 * k_timestep;
    float dx2 = sim->particles[i].v2 * k_timestep;

    sim->particles[i].x1 += dx1;
    sim->particles[i].x2 += dx2;

    // TODO: particles bounce off of walls with raycasting
    for (int tries = 0; tries < 10 && in_wall(sim, sim->particles[i]);
         ++tries) {
      sim->particles[i].x1 -= dx1 / BACKTRACK_PRECISION;
      sim->particles[i].x2 -= dx1 / BACKTRACK_PRECISION;
    }
  }
}

void velocity_to_grid(simulator_s *sim) {
  float weight_sum;
  for (int i = 0; i < sim->x1n; ++i) {
  }
  for (int i = 0; i < PARTICLE_COUNT; ++i) {
    particle_s *p = sim->particles;
    int c_x1 = p[i].x1 / CELL_W;
    int c_x2 = p[i].x2 / CELL_H;

    float l_x1 = p[i].x1 - c_x1;
    float l_x2 = p[i].x2 - c_x2;
  }
}

void incompress(simulator_s *sim, int iters) {
  // TODO: compensate for high density areas
  for (int n = 0; n < iters; ++n) {
    for (int i = 1; i < SIM_H - 1; ++i) {
      for (int j = 1; j < SIM_W - 1; ++j) {
        if (wall_at(sim, i, j)) {
          continue;
        }

        bool sl = sim->s[i][j - 1] != WALL;
        bool sr = sim->s[i][j + 1] != WALL;
        bool su = sim->s[i - 1][j] != WALL;
        bool sd = sim->s[i + 1][j] != WALL;
        float s = sl + sr + su + sd;
        if (sl & sr & su & sd == 0) {
          continue;
        }

        const int v1_i = i * (SIM_W + 1) + j; // left
        const int v2_i = i * SIM_W + j;       // top

        float *vl = &sim->v1[v1_i];
        float *vr = &sim->v1[v1_i + 1];
        float *vu = &sim->v2[v2_i];
        float *vd = &sim->v2[v2_i + (SIM_W + 1)];

        float flow = k_relax * (-*vl + *vr + -*vu + *vd);

        *vl += flow * (sl / s);
        *vr -= flow * (sr / s);
        *vu += flow * (su / s);
        *vd -= flow * (sd / s);
      }
    }
  }
}

void distribute_to_particles(simulator_s *sim) {}

void render_simulation(SDL_Renderer *renderer, simulator_s *sim) {

  const SDL_Color c_wall = {155, 155, 155};
  const SDL_Color c_fluid = {200, 220, 255};
  const SDL_Color c_air = {255, 255, 255};

  for (int i = 0; i < SIM_H; ++i) {
    for (int j = 0; j < SIM_W; ++j) {
      SDL_FRect rect = {
          .x = CELL_W * j, .y = CELL_H * i, .w = CELL_W, .h = CELL_H};
      switch (sim->s[i][j]) {
      case WATER:
        SDL_SetRenderDrawColor(renderer, c_fluid.r, c_fluid.g, c_fluid.b,
                               SDL_ALPHA_OPAQUE);
        break;
      case WALL:
        SDL_SetRenderDrawColor(renderer, c_wall.r, c_wall.g, c_wall.b,
                               SDL_ALPHA_OPAQUE);
        break;
      case AIR:
        SDL_SetRenderDrawColor(renderer, c_air.r, c_air.g, c_air.b,
                               SDL_ALPHA_OPAQUE);
        break;
      }
      SDL_RenderFillRect(renderer, &rect);
    }
  }
}

int main() {
  printf("\n");

  SDL_Init(SDL_INIT_VIDEO);

  SDL_Window *window = SDL_CreateWindow("ffs", WINDOW_W, WINDOW_H,
                                        SDL_WINDOW_HIGH_PIXEL_DENSITY);
  SDL_Renderer *renderer = SDL_CreateRenderer(window, NULL);

  printf("on wayland, try SDL_VIDEODRIVER=wayland\n\n");

  printf("content dpi scaling: %.2f\n",
         SDL_GetDisplayContentScale(SDL_GetDisplayForWindow(window)));
  printf("window dpi scaling: %.2f\n", SDL_GetWindowDisplayScale(window));
  printf("window pixel density: %.2f\n\n", SDL_GetWindowPixelDensity(window));

  float dpi = SDL_GetWindowPixelDensity(window);
  SDL_SetRenderScale(renderer, dpi, dpi);

  printf(" -- running flip fluid simulator -- \n\n");

  bool running = true;

  float update_time = 0.0;
  long long cycles = 0;

  simulator_s sim = initialise();

  while (running) {
    float t0 = now();

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

    // simulation code

    move_particles(&sim);
    velocity_to_grid(&sim);
    incompress(&sim, k_iters);
    distribute_to_particles(&sim);

    float update1 = now();
    update_time += update1 - t0;

    // rendering code

    render_simulation(renderer, &sim);

    SDL_RenderPresent(renderer);

    float t1 = now();
    SDL_Delay(1000 * fmax(frametime - (now() - t0), 0.f));
    ++cycles;
  }

  printf("\n -- finished -- \n\n");

  printf("time per cycle: %f ms\n\n", 1000.0 * update_time / cycles);

  return 0;
}