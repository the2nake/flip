#include <SDL3/SDL_pixels.h>
#include <SDL3/SDL_render.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <SDL3/SDL.h>
#include <SDL3/SDL_main.h>

#include "flip.h"

/*
TODO: future considerations
cache locality
more efficient data type for copy away from grid velocities
using restrict keyword on pointers?
*/

//=============
//    MAIN
//=============

void setup_scaling(SDL_Window *window, SDL_Renderer *renderer);

void initialise();

void advection();
void v_particles_to_grid();
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

  float sum_t1 = 0.f, sum_t2 = 0.f, sum_t3 = 0.f, sum_t4 = 0.f;

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

    sum_t1 += time_of(advection());
    sum_t2 += time_of(v_particles_to_grid());
    sum_t3 += time_of(projection(k_iters));
    sum_t4 += time_of(v_grid_to_particles());

    update_time += now() - t0;

    // rendering code

    render_simulation(renderer);

    SDL_RenderPresent(renderer);

    SDL_Delay(1000 * fmax(k_frametime - (now() - t0), 0.f));
    ++cycles;
  }

  free(particles);

  printf("\n -- finished -- \n\n");

  printf("time per cycle: %f ms\n", 1000.0 * update_time / cycles);
  printf("   advect: %f ms\n", 1000.0 * sum_t1 / cycles);
  printf("  to_grid: %f ms\n", 1000.0 * sum_t2 / cycles);
  printf("  project: %f ms\n", 1000.0 * sum_t3 / cycles);
  printf("  to_part: %f ms\n", 1000.0 * sum_t4 / cycles);

  printf("\n");

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

void set_state_half_water_box();
void reset_velocity_field();
void distribute_particles();
int count_water_cells();
void update_prior_velocities();

void initialise() {
  set_state_half_water_box();
  distribute_particles();
  reset_velocity_field();
  update_prior_velocities();
}

void set_state_half_water_box() {
  for (int i = 0; i < SIM_H; ++i) {
    for (int j = 0; j < SIM_W; ++j) {
      if (cell_on_edge(i, j)) {
        states[i][j] = solid_e;
      } else if (j < SIM_W / 2) {
        states[i][j] = water_e;
      } else {
        states[i][j] = air_e;
      }
    }
  }
}

void distribute_particles() {
  n_particles = PARTICLES_PER_CELL * count_water_cells();
  particles = malloc(n_particles * sizeof(particle_t));
  particles_w = malloc(n_particles * 8 * sizeof(cell_weight_t));

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

void reset_velocity_field() {
  for (int i = 0; i < V1N; ++i) {
    v1[i] = 0.f;
    w1[i] = 1.f;
  }

  for (int i = 0; i < V2N; ++i) {
    v2[i] = 0.f;
    w2[i] = 1.f;
  }
}

void update_prior_velocities() {
  memcpy(v1_prior, v1, V1N * sizeof(float));
  memcpy(v2_prior, v2, V2N * sizeof(float));
}

void advection() {
  // TODO? separate particles using LUT method, radix sorting
  // set fluid cells to air
  for (int i = 0; i < SIM_H; ++i) {
    for (int j = 0; j < SIM_W; ++j) {
      if (states[i][j] != solid_e) {
        states[i][j] = air_e;
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
    for (int tries = 0; tries < 10 && particle_in(&particles[i], solid_e);
         ++tries) {
      particles[i].x1 -= dx1 / BACKTRACK_PRECISION;
      particles[i].x2 -= dx2 / BACKTRACK_PRECISION;
    }
    set_cell_at(&particles[i], water_e);
  }
}

void compute_weights(particle_t *p, cell_weight_t *w);
void add_v1_weight(cell_weight_t *cell, float v); // returns true if successful
void add_v2_weight(cell_weight_t *cell, float v); // returns true if successful
void enforce_solid_velocity_field() {}            // TODO!

void v_particles_to_grid() {
  particle_t *p = particles;
  reset_velocity_field();
  for (int i = 0; i < n_particles; ++i) {
    // somehow this looped list pairing thing is faster by 0.1 ms and easier to
    // read. vectorization? i don't even know

    compute_weights(&p[i], &particles_w[8 * i]);

    // store weights
    for (int j = 0; j < 4; ++j) {
      cell_weight_t *v1_w = &particles_w[8 * i + j];
      cell_weight_t *v2_w = &particles_w[8 * i + 4 + j];
      add_v1_weight(v1_w, p[i].v1);
      add_v2_weight(v2_w, p[i].v2);
    }
  }

  for (int i = 0; i < V1N; ++i) {
    if (isfinite(v1[i]) && w1[i])
      v1[i] /= w1[i];
  }
  for (int i = 0; i < V2N; ++i) {
    if (isfinite(v2[i]) && w2[i])
      v2[i] /= w2[i];
  }

  enforce_solid_velocity_field();
  update_prior_velocities();
}

void compute_weights(particle_t *p, cell_weight_t *w) {
  int v1_i = particle_v1_index(p); // index of top-left x vel
  int v2_i = particle_v2_index(p); // index of top-left y vel

  int v1_row, v1_col, v2_row, v2_col;
  coordinates_from_index(v1_i, SIM_W + 1, &v1_row, &v1_col);
  coordinates_from_index(v2_i, SIM_W, &v2_row, &v2_col);

  // x velocities: no change on x, staggered upwards by CELL_H / 2
  // y velocities: no change on y, staggered left by CELL_W / 2
  float v1_dx, v1_dy, v2_dx, v2_dy;
  position_in_v1_grid(p, v1_row, v1_col, &v1_dx, &v1_dy);
  position_in_v2_grid(p, v2_row, v2_col, &v2_dx, &v2_dy);

  // clang-format off
  w[0] = (cell_weight_t){v1_i                  , v1_dx            * v1_dy           };
  w[1] = (cell_weight_t){v1_i + (SIM_W + 1)    , v1_dx            * (CELL_H - v1_dy)};
  w[2] = (cell_weight_t){v1_i               + 1, (CELL_W - v1_dx) * v1_dy           };
  w[3] = (cell_weight_t){v1_i + (SIM_W + 1) + 1, (CELL_W - v1_dx) * (CELL_H - v1_dy)};

  w[4] = (cell_weight_t){v2_i            ,           v2_dx  *           v2_dy };
  w[5] = (cell_weight_t){v2_i + SIM_W    ,           v2_dx  * (CELL_H - v2_dy)};
  w[6] = (cell_weight_t){v2_i         + 1, (CELL_W - v2_dx) *           v2_dy };
  w[7] = (cell_weight_t){v2_i + SIM_W + 1, (CELL_W - v2_dx) * (CELL_H - v2_dy)};
  // clang-format on
}

void add_v1_weight(cell_weight_t *c, float v) {
  if (!in_rangei(c->i, 0, V1N)) {
    c->w = 0.f;
    return;
  }

  int i = c->i / (SIM_W + 1);
  int j = c->i % (SIM_W + 1);

  bool v1_is_water = cell_is(i, j - 1, water_e) || cell_is(i, j, water_e);
  bool v1_touches_solid = cell_is(i, j - 1, solid_e) || cell_is(i, j, solid_e);
  if (v1_is_water && !v1_touches_solid) {
    v1[c->i] += c->w * v;
    w1[c->i] += c->w;
  } else {
    v1[c->i] = v1_touches_solid ? 0.f : NAN;
    w1[c->i] = 0.f;
    c->w = 0.f;
  }
}

void add_v2_weight(cell_weight_t *c, float v) {
  if (!in_rangei(c->i, 0, V2N)) {
    c->w = 0.f;
    return;
  }

  int i = c->i / SIM_W;
  int j = c->i % SIM_W;

  bool v2_is_water = cell_is(i - 1, j, water_e) || cell_is(i, j, water_e);
  bool v2_touches_solid = cell_is(i - 1, j, solid_e) || cell_is(i, j, solid_e);
  if (v2_is_water && !v2_touches_solid) {
    v2[c->i] += c->w * v;
    w2[c->i] += c->w;
  } else {
    v2[c->i] = v2_touches_solid ? 0.f : NAN;
    w2[c->i] = 0.f;
    c->w = 0.f;
  }
}

// enforce inflow = outflow iteratively
void projection(int iters) {
  // TODO! fix particle lagging in the air, check its velocity field
  // TODO: compensate for high density areas
  for (int n = 0; n < iters; ++n) {
    for (int i = 1; i < SIM_H - 1; ++i) {
      for (int j = 1; j < SIM_W - 1; ++j) {
        if (states[i][j] != water_e)
          continue;

        bool sl = states[i][j - 1] != solid_e;
        bool sr = states[i][j + 1] != solid_e;
        bool su = states[i - 1][j] != solid_e;
        bool sd = states[i + 1][j] != solid_e;
        int s = sl + sr + su + sd;
        if (!s) {
          continue;
        }

        const int v1_i = i * (SIM_W + 1) + j; // left
        const int v2_i = i * SIM_W + j;       // top

        float *vl = &v1[v1_i];
        float *vr = &v1[v1_i + 1];
        float *vu = &v2[v2_i];
        float *vd = &v2[v2_i + SIM_W];

        float net = (*vr + *vd - *vl - *vu);
        float flow = k_relax * net / s;

        *vl += sl * flow;
        *vr -= sr * flow;
        *vu += su * flow;
        *vd -= sd * flow;
      }
    }
  }
  // assign velocities between non-water cells to nan
  for (int i = 0; i < SIM_H; ++i) {
    for (int j = 0; j < SIM_W + 1; ++j) {
      bool borders_water = cell_is(i, j - 1, water_e) || cell_is(i, j, water_e);
      if (!borders_water) {
        v1[i * (SIM_W + 1) + j] = NAN;
      }
    }
  }

  for (int i = 0; i < SIM_H + 1; ++i) {
    for (int j = 0; j < SIM_W; ++j) {
      bool borders_water = cell_is(i - 1, j, water_e) || cell_is(i, j, water_e);
      if (!borders_water) {
        v2[i * SIM_W + j] = NAN;
      }
    }
  }
}

typedef enum { v1_e, v2_e } field_e_t;

void update_particle(int i, field_e_t field, bool pic);

void v_grid_to_particles() {
  bool pic = false;

  for (int i = 0; i < n_particles; ++i) {
    update_particle(i, v1_e, pic);
    update_particle(i, v2_e, pic);
  }
}

void update_particle(int i, field_e_t field, bool pic) {
  assert(0 <= i && i < n_particles);

  // clang-format off
  cell_weight_t *c = &particles_w[8 * i + (field == v2_e) * 4];
  float *v         = field == v1_e ? v1               : v2;
  float *v_prior   = field == v1_e ? v1_prior         : v2_prior;
  float *v_out     = field == v1_e ? &particles[i].v1 : &particles[i].v2;
  // clang-format on

  float dv = 0.f;
  float w = 0.f;

  for (int j = 0; j < 4; ++j, ++c) {
    if (c->i < 0 || !c->w || !isfinite(v[c->i]))
      continue;

    float change = v[c->i] - (pic ? 0.f : v_prior[c->i]);
    dv += change * c->w;
    w += c->w;
  }

  *v_out = dv / w + (pic ? 0.f : *v_out);
}

void render_particles(SDL_Renderer *renderer);
void render_velocities(SDL_Renderer *renderer);

void render_simulation(SDL_Renderer *renderer) {
  // draw cells
  const SDL_Color c_wall = {155, 155, 155};
  const SDL_Color c_fluid = {200, 220, 255};
  const SDL_Color c_air = {255, 255, 255};

  for (int i = 0; i < SIM_H; ++i) {
    for (int j = 0; j < SIM_W; ++j) {
      SDL_FRect rect = {
          .x = CELL_W * j, .y = CELL_H * i, .w = CELL_W, .h = CELL_H};
      switch (states[i][j]) {
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

  render_particles(renderer);
  // render_velocity_field(renderer);
}

void render_particles(SDL_Renderer *renderer) {
  SDL_SetRenderDrawColor(renderer, 0, 0, 0, SDL_ALPHA_OPAQUE);
  for (int i = 0; i < n_particles; ++i) {
    particle_t p = particles[i];
    SDL_RenderPoint(renderer, p.x1, p.x2);
  }
}

void render_velocities(SDL_Renderer *renderer) {
  SDL_SetRenderDrawColor(renderer, 180, 255, 0, SDL_ALPHA_OPAQUE);
  // render x velocities
  for (int i = 0; i < V1N; ++i) {
    int row = 0;
    int col = 0;
    coordinates_from_index(i, SIM_W + 1, &row, &col);
    int origin_x = col * CELL_W;
    int origin_y = (row + 0.5) * CELL_H;
    SDL_RenderLine(renderer, origin_x, origin_y, origin_x + v1[i], origin_y);
  }
  // render y velocities
  for (int i = 0; i < V2N; ++i) {
    int row = 0, col = 0;
    coordinates_from_index(i, SIM_W, &row, &col);
    int origin_x = (col + 0.5) * CELL_W;
    int origin_y = row * CELL_H;
    SDL_RenderLine(renderer, origin_x, origin_y, origin_x, origin_y + v2[i]);
  }
}