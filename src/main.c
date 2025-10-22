#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <SDL3/SDL.h>
#include <SDL3/SDL_main.h>

/*
TODO: future considerations
cache locality
more efficient data type for copy away from grid velocities
using restrict keyword on pointers?
*/

//==================
//    PARAMETERS
//==================

constexpr int SIM_W = 20;
constexpr int SIM_H = 20;

constexpr int CELL_W = 10;
constexpr int CELL_H = CELL_W;

constexpr int DENSITY = 2;
constexpr int PARTICLES_PER_CELL = DENSITY * DENSITY;
constexpr int MAX_PARTICLES = PARTICLES_PER_CELL * SIM_W * SIM_H;

constexpr int BACKTRACK_PRECISION = 4;

// TODO! convert window_scale to internal difference between graphical and
// internal scale

constexpr int WINDOW_W = SIM_W * CELL_W;
constexpr int WINDOW_H = SIM_H * CELL_H;
constexpr float WINDOW_SCALE = 2.f;

const float k_frametime = 1.0 / 60.0;
const float k_timestep = 1.0 / 60.0;
const float k_gravity = 9.81;
const float k_relax = 1.9;
const int k_iters = 40;

float now() { return (float)clock() / CLOCKS_PER_SEC; }

typedef enum { solid_e = 0, water_e = 1, air_e = 2 } state_e_t;

typedef struct {
  float x1; // x position
  float x2; // y position
  float v1; // x velocity
  float v2; // y velocity
} particle_t;

constexpr int V1N = (SIM_W + 1) * SIM_H;
constexpr int V2N = SIM_W * (SIM_H + 1);

state_e_t s[SIM_H][SIM_W]; // states
float v1[V1N];             // horizontal velocity
float v2[V2N];             // vertical velocity
float w1[V1N];             // fixed point weights
float w2[V2N];             // fixed point weights

particle_t *particles;
int n_particles = MAX_PARTICLES;

//=============
//    DEBUG
//=============

void print_v2s() {
  for (int i = 0; i < SIM_H + 1; ++i) {
    for (int j = 0; j < SIM_H; ++j) {
      printf("%4.1f ", v2[SIM_W * i + j]);
    }
    printf("\n");
  }
  printf("========================================\n");
}

void print_w2s() {
  for (int i = 0; i < SIM_H + 1; ++i) {
    for (int j = 0; j < SIM_H; ++j) {
      printf("%.1f ", w2[SIM_W * i + j]);
    }
    printf("\n");
  }
  printf("========================================\n");
}

//=============
//    UTILS
//=============

float *x_vel(int i) {
  if (0 <= i && i < V1N) {
    return &v1[i];
  } else {
    return NULL;
  }
}

float *y_vel(int i) {
  if (0 <= i && i < V2N) {
    return &v2[i];
  } else {
    return NULL;
  }
}

int in_rangei(int value, int lo, int hi) { return lo <= value && value <= hi; }

float in_rangef(float value, float lo, float hi) {
  return lo <= value && value <= hi;
}

float clamp(float value, const float lo, const float hi) {
  if (value > hi) {
    return hi;
  }
  if (value < lo) {
    return lo;
  }
  return value;
}

bool cell_in_bounds(int i, int j) {
  return 0 <= i && i < SIM_H && 0 <= j && j < SIM_W;
}

bool particle_in_bounds(particle_t *p) {
  return 0.f <= p->x1 && p->x1 <= SIM_W * CELL_W && 0.f <= p->x2 &&
         p->x2 <= SIM_H * CELL_H;
}

void particle_enforce_bounds(particle_t *p) {
  if (!particle_in_bounds(p)) {
    p->x1 = clamp(p->x1, 0.f, SIM_W * CELL_W);
    p->x2 = clamp(p->x2, 0.f, SIM_H * CELL_H);
  }
}

bool cell_on_edge(int i, int j) {
  if (i == 0 || i == SIM_H - 1) {
    return true;
  }
  if (j == 0 || j == SIM_W - 1) {
    return true;
  }
  return false;
}

bool cell_is(int i, int j, state_e_t state) {
  if (cell_in_bounds(i, j)) {
    return s[i][j] == state;
  }
  return false;
}

void set_cell_at(particle_t particle, state_e_t state) {
  int j = particle.x1 / CELL_W;
  int i = particle.x2 / CELL_H;
  if (cell_in_bounds(i, j)) {
    s[i][j] = state;
  }
}

bool particle_in(particle_t particle, state_e_t state) {
  int j = particle.x1 / CELL_W;
  int i = particle.x2 / CELL_H;
  return cell_is(i, j, state); // check bounds to not lose particles
}

float particle_velocity(particle_t *p) {
  return sqrtf(p->v1 * p->v1 + p->v2 * p->v2);
}

//=============
//    MAIN
//=============

void setup_scaling(SDL_Window *window, SDL_Renderer *renderer);

void initialise();
void set_state_half_water_box();
void reset_velocity_field();
void distribute_particles();
int count_water_cells();

void advection();

void v_particles_to_grid();
void add_v1_weight(int idx, float v, float w);
void add_v2_weight(int idx, float v, float w);

void projection(int iters);

void v_grid_to_particles();

void render_simulation(SDL_Renderer *renderer);

int main() {
  printf("\n");

  SDL_Init(SDL_INIT_VIDEO);

  SDL_Window *window =
      SDL_CreateWindow("ffs", WINDOW_W * WINDOW_SCALE, WINDOW_H * WINDOW_SCALE,
                       SDL_WINDOW_HIGH_PIXEL_DENSITY);
  SDL_Renderer *renderer = SDL_CreateRenderer(window, NULL);

  setup_scaling(window, renderer);

  printf(" -- running flip fluid simulator -- \n\n");

  bool running = true;

  float update_time = 0.0;
  long long cycles = 0;

  initialise();

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

    advection();
    v_particles_to_grid();
    // projection(k_iters);
    v_grid_to_particles();

    update_time += now() - t0;

    // rendering code

    render_simulation(renderer);

    SDL_RenderPresent(renderer);

    SDL_Delay(1000 * fmax(k_frametime - (now() - t0), 0.f));
    ++cycles;
  }

  free(particles);

  printf("\n -- finished -- \n\n");

  printf("time per cycle: %f ms\n\n", 1000.0 * update_time / cycles);

  return 0;
}

void setup_scaling(SDL_Window *w, SDL_Renderer *r) {
  printf("on wayland, try SDL_VIDEODRIVER=wayland\n\n");
  float content_dpi = SDL_GetDisplayContentScale(SDL_GetDisplayForWindow(w));
  float window_dpi = SDL_GetWindowDisplayScale(w);
  float pixel_density = SDL_GetWindowPixelDensity(w);

  printf("content dpi scaling: %.2f\n", content_dpi);
  printf("window dpi scaling: %.2f\n", window_dpi);
  printf("window pixel density: %.2f\n\n", pixel_density);

  float scale = pixel_density * WINDOW_SCALE;
  SDL_SetRenderScale(r, scale, scale);
}

void initialise() {
  set_state_half_water_box();
  reset_velocity_field();
  distribute_particles();
}

void reset_velocity_field() {
  for (int i = 0; i < V1N; ++i) {
    v1[i] = 0.f;
    w1[i] = 1.f;
  }

  for (int i = 0; i < V2N; ++i) {
    v2[i] = 0.f;
    w1[i] = 1.f;
  }
}

void set_state_half_water_box() {
  for (int i = 0; i < SIM_H; ++i) {
    for (int j = 0; j < SIM_W; ++j) {
      if (cell_on_edge(i, j)) {
        s[i][j] = solid_e;
      } else if (j < SIM_W / 2) {
        s[i][j] = water_e;
      } else {
        s[i][j] = air_e;
      }
    }
  }
}

void distribute_particles() {
  particles = malloc(n_particles * sizeof(particle_t));
  n_particles = PARTICLES_PER_CELL * count_water_cells();

  constexpr float x_gap = (float)CELL_W / DENSITY;
  constexpr float y_gap = (float)CELL_H / DENSITY;
  constexpr float left_padding = x_gap / 2.f;
  constexpr float top_padding = y_gap / 2.f;

  int idx = 0;
  for (int i = 0; i < SIM_H; ++i) {
    for (int j = 0; j < SIM_W; ++j) {
      if (!cell_is(i, j, water_e))
        continue;

      for (int k = 0; k < PARTICLES_PER_CELL; ++k) {
        particles[idx] = (particle_t){
            .x1 = j * CELL_W + left_padding + x_gap * (k % DENSITY),
            .x2 = i * CELL_H + top_padding + y_gap * (int)(k / DENSITY),
            .v1 = 0.f,
            .v2 = 0.f};
        ++idx;
      }
    }
  }
}

int count_water_cells() {
  int res = 0;
  for (int i = 0; i < SIM_H; ++i) {
    for (int j = 0; j < SIM_W; ++j) {
      if (cell_is(i, j, water_e))
        ++res;
    }
  }
  return res;
}

void advection() {
  // TODO? separate particles using LUT method, radix sorting
  // set fluid cells to air
  for (int i = 0; i < SIM_H; ++i) {
    for (int j = 0; j < SIM_W; ++j) {
      if (s[i][j] != solid_e) {
        s[i][j] = air_e;
      }
    }
  }

  // update particle locations and cell states
  for (int i = 0; i < n_particles; ++i) {
    particles[i].v2 += k_gravity * k_timestep;
    float dx1 = particles[i].v1 * k_timestep;
    float dx2 = particles[i].v2 * k_timestep;

    particles[i].x1 += dx1;
    particles[i].x2 += dx2;

    particle_enforce_bounds(&particles[i]);
    // TODO? particles bounce off of walls with raycasting
    for (int tries = 0; tries < 10 && particle_in(particles[i], solid_e);
         ++tries) {
      particles[i].x1 -= dx1 / BACKTRACK_PRECISION;
      particles[i].x2 -= dx2 / BACKTRACK_PRECISION;
    }
    set_cell_at(particles[i], water_e);
  }
}

typedef struct {
  int idx;
  float weight;
} cell_weight_t;

void v_particles_to_grid() {
  particle_t *p = particles;
  reset_velocity_field();
  for (int i = 0; i < n_particles; ++i) {
    // x velocities: no change on x, staggered upwards by CELL_H / 2
    int v1_col = p[i].x1 / CELL_W;
    int v1_row = p[i].x2 / CELL_H - 0.5;
    int v1_i = v1_row * (SIM_W + 1) + v1_col; // index of top-left x vel

    float v1_dx = p[i].x1 - v1_col * CELL_W;
    float v1_dy = p[i].x2 - (v1_row + 0.5) * CELL_H;

    // y velocities: no change on y, staggered left by CELL_W / 2
    int v2_col = p[i].x1 / CELL_W - 0.5;
    int v2_row = p[i].x2 / CELL_H;
    int v2_i = v2_row * SIM_W + v2_col; // index of top-left y vel

    float v2_dx = p[i].x1 - (v2_col + 0.5) * CELL_W;
    float v2_dy = p[i].x2 - v2_row * CELL_H;

    // clang-format off
    cell_weight_t v1_weights[4] = {
        {v1_i                  , v1_dx            * v1_dy           },
        {v1_i + (SIM_W + 1)    , v1_dx            * (CELL_H - v1_dy)},
        {v1_i               + 1, (CELL_W - v1_dx) * v1_dy           },
        {v1_i + (SIM_W + 1) + 1, (CELL_W - v1_dx) * (CELL_H - v1_dy)}};

    cell_weight_t v2_weights[4] = {
        {v2_i            ,           v2_dx  *           v2_dy },
        {v2_i + SIM_W    ,           v2_dx  * (CELL_H - v2_dy)},
        {v2_i         + 1, (CELL_W - v2_dx) *           v2_dy },
        {v2_i + SIM_W + 1, (CELL_W - v2_dx) * (CELL_H - v2_dy)}};
    // clang-format on

    // add weights
    for (int i = 0; i < 4; ++i) {
      add_v1_weight(v1_weights[i].idx, p[i].v1, v1_weights[i].weight);
      add_v2_weight(v2_weights[i].idx, p[i].v2, v2_weights[i].weight);
    }

    // somehow this looped list pairing thing is faster by 0.1 ms and easier to
    // read. vectorization? i don't even know
  }

  for (int i = 0; i < V1N; ++i) {
    if (isfinite(v1[i]) && w1[i])
      v1[i] /= w1[i];
  }
  for (int i = 0; i < V2N; ++i) {
    if (isfinite(v2[i]) && w2[i])
      v2[i] /= w2[i];
  }
}

void add_v1_weight(int idx, float v, float w) {
  if (!in_rangei(idx, 0, V1N))
    return;

  int i = idx / (SIM_W + 1);
  int j = idx % (SIM_W + 1);

  bool v1_is_water = cell_is(i, j - 1, water_e) || cell_is(i, j, water_e);
  if (v1_is_water) {
    v1[idx] += w * v;
    w1[idx] += w;
  } else {
    v1[idx] = NAN;
    w1[idx] = 0;
  }
}

void add_v2_weight(int idx, float v, float w) {
  if (!in_rangei(idx, 0, V2N))
    return;

  int i = idx / SIM_W;
  int j = idx % SIM_W;

  bool v2_is_water = cell_is(i - 1, j, water_e) || cell_is(i, j, water_e);
  if (v2_is_water) {
    v2[idx] += w * v;
    w2[idx] += w;
  } else {
    v2[idx] = NAN;
    w2[idx] = 0;
  }
}

void projection(int iters) {
  // TODO! fix for handling nan cells
  // TODO: compensate for high density areas
  // TODO? do you handle incompressibility including air/water boundaries?
  for (int n = 0; n < iters; ++n) {
    for (int i = 1; i < SIM_H - 1; ++i) {
      for (int j = 1; j < SIM_W - 1; ++j) {
        if (!cell_is(i, j, water_e))
          continue;

        bool sl = s[i][j - 1] == water_e;
        bool sr = s[i][j + 1] == water_e;
        bool su = s[i - 1][j] == water_e;
        bool sd = s[i + 1][j] == water_e;
        float s = sl + sr + su + sd;
        if (sl & sr & su & sd == 0) {
          continue;
        }

        const int v1_i = i * (SIM_W + 1) + j; // left
        const int v2_i = i * SIM_W + j;       // top

        float *vl = &v1[v1_i];
        float *vr = &v1[v1_i + 1];
        float *vu = &v2[v2_i];
        float *vd = &v2[v2_i + (SIM_W + 1)];

        float flow = k_relax * (-*vl + *vr + -*vu + *vd);

        *vl += flow * (sl / s);
        *vr -= flow * (sr / s);
        *vu += flow * (su / s);
        *vd -= flow * (sd / s);
      }
    }
  }
}

void v_grid_to_particles() {}

void render_simulation(SDL_Renderer *renderer) {
  // draw cells
  const SDL_Color c_wall = {155, 155, 155};
  const SDL_Color c_fluid = {200, 220, 255};
  const SDL_Color c_air = {255, 255, 255};

  for (int i = 0; i < SIM_H; ++i) {
    for (int j = 0; j < SIM_W; ++j) {
      SDL_FRect rect = {
          .x = CELL_W * j, .y = CELL_H * i, .w = CELL_W, .h = CELL_H};
      switch (s[i][j]) {
      case water_e:
        SDL_SetRenderDrawColor(renderer, c_fluid.r, c_fluid.g, c_fluid.b,
                               SDL_ALPHA_OPAQUE);
        break;
      case solid_e:
        SDL_SetRenderDrawColor(renderer, c_wall.r, c_wall.g, c_wall.b,
                               SDL_ALPHA_OPAQUE);
        break;
      case air_e:
        SDL_SetRenderDrawColor(renderer, c_air.r, c_air.g, c_air.b,
                               SDL_ALPHA_OPAQUE);
        break;
      }
      SDL_RenderFillRect(renderer, &rect);
    }
  }
  // draw particles
  SDL_SetRenderDrawColor(renderer, 0, 0, 0, SDL_ALPHA_OPAQUE);
  for (int i = 0; i < n_particles; ++i) {
    particle_t p = particles[i];
    SDL_RenderPoint(renderer, p.x1, p.x2);
  }
}
