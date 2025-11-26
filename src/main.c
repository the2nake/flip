#include <SDL3/SDL.h>
#include <SDL3/SDL_events.h>
#include <SDL3/SDL_keycode.h>
#include <SDL3/SDL_log.h>
#include <SDL3/SDL_main.h>
#include <SDL3/SDL_oldnames.h>
#include <SDL3/SDL_pixels.h>
#include <SDL3/SDL_render.h>
#include <SDL3/SDL_stdinc.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "flip.h"

//=============
//    MAIN
//=============

void setup_scaling(SDL_Window *window, SDL_Renderer *renderer);

void initialise();

void compute_density();
void advect();
void v_to_grid();
void project(int iters);
void v_to_particles();

void render_simulation(SDL_Renderer *renderer);

bool paused = true;
bool show_particles = false;
bool show_velocities = false;

void set_colora(SDL_Renderer *renderer, const SDL_Color *color,
                const uint8_t a) {
  SDL_SetRenderDrawColor(renderer, color->r, color->g, color->b, a);
}

void set_color(SDL_Renderer *renderer, const SDL_Color *color) {
  SDL_SetRenderDrawColor(renderer, color->r, color->g, color->b, color->a);
}

int main() {
  printf("\n");

  SDL_Init(SDL_INIT_VIDEO);

  const int window_pixel_w = WINDOW_W * WINDOW_SCALE;
  const int window_pixel_h = WINDOW_H * WINDOW_SCALE;

  SDL_Window *window = SDL_CreateWindow("ffs", window_pixel_w, window_pixel_h,
                                        SDL_WINDOW_HIGH_PIXEL_DENSITY);
  SDL_Renderer *renderer = SDL_CreateRenderer(window, NULL);
  SDL_SetRenderDrawBlendMode(renderer, SDL_BLENDMODE_BLEND);

  setup_scaling(window, renderer);

  printf(" -- running flip fluid simulator -- \n\n");

  bool running = true;

  float update_time = 0.0;
  long long cycles = 0;

  float sum_t1 = 0.f, sum_t2 = 0.f, sum_t3 = 0.f, sum_t4 = 0.f, sum_t5 = 0.f;

  initialise();

  while (running) {
    SDL_Event event;

    while (SDL_PollEvent(&event)) {
      switch (event.type) {
        case SDL_EVENT_QUIT:
          running = false;
          break;
        case (SDL_EVENT_KEY_UP):
          switch (event.key.key) {
            case SDLK_P:
              show_particles = !show_particles;
              break;
            case SDLK_V:
              show_velocities = !show_velocities;
              break;
            case SDLK_SPACE:
              paused = !paused;
              break;
            default:
              break;
          }
        default:
          break;
      }
    }

    float t0 = now();

    if (!paused) {
      // simulation code
      sum_t1 += time_of(compute_density());
      sum_t2 += time_of(advect());
      sum_t3 += time_of(v_to_grid());
      sum_t4 += time_of(project(k_iters));
      sum_t5 += time_of(v_to_particles());

      update_time += now() - t0;
      ++cycles;
    }

    render_simulation(renderer);

    if (paused) {
      set_color(renderer, &(SDL_Color){0, 0, 0, SDL_ALPHA_OPAQUE});
      SDL_RenderDebugText(renderer, window_pixel_w * 0.1, window_pixel_h * 0.1,
                          "PAUSED: [space] to resume");
    }

    SDL_RenderPresent(renderer);

    SDL_Delay(1000 * fmax(k_frametime - (now() - t0), 0.f));
  }

  free(particles);

  printf("\n -- finished -- \n\n");

  printf("time per cycle: %f ms\n", 1000.f * update_time / cycles);
  printf("  density: %f ms\n", 1000.f * sum_t1 / cycles);
  printf("   advect: %f ms\n", 1000.f * sum_t2 / cycles);
  printf("  to_grid: %f ms\n", 1000.f * sum_t3 / cycles);
  printf("  project: %f ms\n", 1000.f * sum_t4 / cycles);
  printf("  to_part: %f ms\n", 1000.f * sum_t5 / cycles);

  printf("\n");

  return 0;
}

void setup_scaling(SDL_Window *w, SDL_Renderer *r) {
  printf("on wayland, try SDL_VIDEODRIVER=wayland\n\n");
  float content_dpi = SDL_GetDisplayContentScale(SDL_GetDisplayForWindow(w));
  float window_dpi = SDL_GetWindowDisplayScale(w);
  float pixel_density = SDL_GetWindowPixelDensity(w);

  printf(" content dpi scaling: %.2f\n", content_dpi);
  printf("  window dpi scaling: %.2f\n", window_dpi);
  printf("window pixel density: %.2f\n\n", pixel_density);

  SDL_SetRenderScale(r, pixel_density, pixel_density);
}

void set_state_half_water_box();
void distribute_particles();
void reset_velocity_field();
void update_prior_velocities();

void initialise() {
  set_state_half_water_box();
  distribute_particles();
  reset_velocity_field();

  update_prior_velocities();
  compute_density();
}

void set_state_half_water_box() {
  for (int i = 0; i < SIM_H; ++i) {
    for (int j = 0; j < SIM_W; ++j) {
      if (cell_on_edge(i, j)) {
        states[i][j] = solid_e;
      } else if (j * 2 < SIM_W && i > SIM_H / 5) {
        states[i][j] = water_e;
      } else {
        states[i][j] = air_e;
      }
    }
  }
}

int count_water_cells();

void distribute_particles() {
  n_particles = PARTICLES_PER_CELL * count_water_cells();
  particles = malloc(n_particles * sizeof(particle_t));
  vel_ws = malloc(n_particles * 8 * sizeof(vel_weight_t));

  constexpr float x_gap = PARTICLE_SIZE;
  constexpr float y_gap = PARTICLE_SIZE;
  constexpr float left_pad = x_gap / 2.f;
  constexpr float top_pad = y_gap / 2.f;

  int idx = 0;

  for (int i = 0; i < SIM_H; ++i) {
    for (int j = 0; j < SIM_W; ++j) {
      if (!cell_is(i, j, water_e)) continue;

      for (int k = 0; k < PARTICLES_PER_CELL; ++k) {
        particles[idx] = (particle_t){
            .x1 = j * CELL_W + left_pad + x_gap * (k % PARTICLE_PACKING),
            .x2 = i * CELL_H + top_pad + y_gap * (int)(k / PARTICLE_PACKING),
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
      if (cell_is(i, j, water_e)) ++res;
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

void compute_density() {
  constexpr float cell_area = CELL_W * CELL_H;

  for (int i = 0; i < SIM_H; ++i) {
    for (int j = 0; j < SIM_W; ++j) {
      densities[i][j] = states[i][j] == solid_e ? NAN : 0.f;
    }
  }

  for (int i = 0; i < n_particles; ++i) {
    int c_i = particles[i].x2 / CELL_H - 0.5;
    int c_j = particles[i].x1 / CELL_W - 0.5;

    if (c_i < 0 || c_j < 0 || c_i >= SIM_H - 1 || c_j >= SIM_W - 1) continue;

    float c_x1 = (c_j + 0.5) * CELL_W;
    float c_x2 = (c_i + 0.5) * CELL_H;
    float dx1 = particles[i].x1 - c_x1;
    float dx2 = particles[i].x2 - c_x2;

    // clang-format off
    densities[c_i][c_j]         += (CELL_W - dx1) * (CELL_H - dx2) / cell_area;
    densities[c_i][c_j + 1]     +=           dx1  * (CELL_H - dx2) / cell_area;
    densities[c_i + 1][c_j]     += (CELL_W - dx1) *           dx2  / cell_area;
    densities[c_i + 1][c_j + 1] +=           dx1  *           dx2  / cell_area;
    // clang-format on
  }
}

// TODO! fix viscosity
void advect() {
  // TODO? separate particles using LUT method, radix sorting
  // set all non-solid cells to air
  for (int i = 0; i < SIM_H; ++i) {
    for (int j = 0; j < SIM_W; ++j) {
      if (states[i][j] != solid_e) { states[i][j] = air_e; }
    }
  }

  // move particles and update cell states
  for (int i = 0; i < n_particles; ++i) {
    particle_t *p = &particles[i];
    p->v2 += k_gravity * k_timestep;
    float dx1 = p->v1 * k_timestep;
    float dx2 = p->v2 * k_timestep;

    p->x1 += dx1;
    p->x2 += dx2;

    particle_enforce_bounds(p);
    // TODO? particles bounce off of walls with raycasting

    bool in_cell = particle_in(p, solid_e);
    if (in_cell) {
      int c_i = 0, c_j = 0;
      get_particle_cell(p, &c_i, &c_j);
      get_cell_normal(c_i, c_j, &dx1, &dx2);

      int tries = 0;
      for (; tries < BACKTRACK_ATTEMPTS && in_cell; ++tries) {
        p->x1 += BACKTRACK_RANGE * dx1 * CELL_W / BACKTRACK_ATTEMPTS;
        p->x2 += BACKTRACK_RANGE * dx2 * CELL_H / BACKTRACK_ATTEMPTS;
        // particle_enforce_bounds(p);
        in_cell = particle_in(p, solid_e);
      }

      if (in_cell) {
        SDL_LogWarn(SDL_LOG_CATEGORY_APPLICATION, "particle inside a solid");
      } else if (!particle_in_bounds(p)) {
        SDL_LogWarn(SDL_LOG_CATEGORY_APPLICATION, "particle out of bounds");
      }
    }
    set_cell_at(&particles[i], water_e);
  }
}

void compute_weights(particle_t *p, vel_weight_t *vel_w);
void add_weight(field_e_t field, vel_weight_t *vel_w, float vel);

void v_to_grid() {
  reset_velocity_field();

  for (int i = 0; i < n_particles; ++i) {
    // somehow this looped list pairing thing is faster by 0.1 ms and
    // easier to read. vectorization? i don't even know
    compute_weights(&particles[i], &vel_ws[8 * i]);

    for (int j = 0; j < 4; ++j) {
      add_weight(v1_e, &vel_ws[8 * i + j], particles[i].v1);
      add_weight(v2_e, &vel_ws[8 * i + j + 4], particles[i].v2);
    }
  }

  for (int i = 0; i < V1N; ++i) {
    if (isfinite(v1[i]) && w1[i]) v1[i] /= w1[i];
  }
  for (int i = 0; i < V2N; ++i) {
    if (isfinite(v2[i]) && w2[i]) v2[i] /= w2[i];
  }
}

void compute_weights(particle_t *p, vel_weight_t *vel_w) {
  int v1_i = particle_v1_index(p);  // index of top-left x vel
  int v2_i = particle_v2_index(p);  // index of top-left y vel

  int v1_row, v1_col, v2_row, v2_col;
  coordinates_from_index(v1_i, SIM_W + 1, &v1_row, &v1_col);
  coordinates_from_index(v2_i, SIM_W, &v2_row, &v2_col);

  // x velocities: no change on x, staggered upwards by CELL_H / 2
  // y velocities: no change on y, staggered left by CELL_W / 2
  float v1_dx, v1_dy, v2_dx, v2_dy;
  position_in_v1_grid(p, v1_row, v1_col, &v1_dx, &v1_dy);
  position_in_v2_grid(p, v2_row, v2_col, &v2_dx, &v2_dy);

  // clang-format off
  vel_w[0] = (vel_weight_t){v1_i                  , v1_dx            * v1_dy           };
  vel_w[1] = (vel_weight_t){v1_i + (SIM_W + 1)    , v1_dx            * (CELL_H - v1_dy)};
  vel_w[2] = (vel_weight_t){v1_i               + 1, (CELL_W - v1_dx) * v1_dy           };
  vel_w[3] = (vel_weight_t){v1_i + (SIM_W + 1) + 1, (CELL_W - v1_dx) * (CELL_H - v1_dy)};

  vel_w[4] = (vel_weight_t){v2_i            ,           v2_dx  *           v2_dy };
  vel_w[5] = (vel_weight_t){v2_i + SIM_W    ,           v2_dx  * (CELL_H - v2_dy)};
  vel_w[6] = (vel_weight_t){v2_i         + 1, (CELL_W - v2_dx) *           v2_dy };
  vel_w[7] = (vel_weight_t){v2_i + SIM_W + 1, (CELL_W - v2_dx) * (CELL_H - v2_dy)};
  // clang-format on
}

void add_weight(field_e_t field, vel_weight_t *vel_w, float vel) {
  if (!in_rangei(vel_w->i, 0, field == v1_e ? V1N : V2N)) {
    vel_w->w = 0.f;
    return;
  }

  float *vf = field == v1_e ? v1 : v2;
  float *wf = field == v1_e ? w1 : w2;
  int i = vel_w->i / (SIM_W + (field == v1_e));
  int j = vel_w->i % (SIM_W + (field == v1_e));

  bool is_water, touches_solid;

  if (field == v1_e) {
    is_water = cell_is(i, j - 1, water_e) || cell_is(i, j, water_e);
    touches_solid = cell_is(i, j - 1, solid_e) || cell_is(i, j, solid_e);
  } else {
    is_water = cell_is(i - 1, j, water_e) || cell_is(i, j, water_e);
    touches_solid = cell_is(i - 1, j, solid_e) || cell_is(i, j, solid_e);
  }

  if (is_water && !touches_solid) {
    vf[vel_w->i] += vel_w->w * vel;
    wf[vel_w->i] += vel_w->w;
  } else {
    vf[vel_w->i] = touches_solid ? 0.f : NAN;
    wf[vel_w->i] = 0.f;
    vel_w->w = 0.f;
  }
}

// enforce inflow = outflow iteratively
void project(int iters) {
  update_prior_velocities();

  // TODO! fix particle lagging in the air, check its velocity field
  // TODO: compensate for high density areas
  for (int n = 0; n < iters; ++n) {
    for (int i = 1; i < SIM_H - 1; ++i) {
      for (int j = 1; j < SIM_W - 1; ++j) {
        if (states[i][j] != water_e) continue;

        bool sl = states[i][j - 1] != solid_e;
        bool sr = states[i][j + 1] != solid_e;
        bool su = states[i - 1][j] != solid_e;
        bool sd = states[i + 1][j] != solid_e;
        int s = sl + sr + su + sd;
        if (!s) continue;

        const int v1_i = i * (SIM_W + 1) + j;  // left
        const int v2_i = i * SIM_W + j;        // top

        float *vl = &v1[v1_i];
        float *vr = &v1[v1_i + 1];
        float *vu = &v2[v2_i];
        float *vd = &v2[v2_i + SIM_W];

        float net = (*vr + *vd - *vl - *vu) -
                    k_stiffness * (densities[i][j] - PARTICLES_PER_CELL);
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
      if (!borders_water) v1[i * (SIM_W + 1) + j] = NAN;
    }
  }
  for (int i = 0; i < SIM_H + 1; ++i) {
    for (int j = 0; j < SIM_W; ++j) {
      bool borders_water = cell_is(i - 1, j, water_e) || cell_is(i, j, water_e);
      if (!borders_water) v2[i * SIM_W + j] = NAN;
    }
  }
}

void update_particle(int i, field_e_t field, float pic);

void v_to_particles() {
  for (int i = 0; i < n_particles; ++i) {
    update_particle(i, v1_e, k_flip);
    update_particle(i, v2_e, k_flip);
  }
}

void update_particle(int i, field_e_t field, float flip) {
  assert(0 <= i && i < n_particles);
  assert(-0.01 <= flip && flip <= 1.01);

  // clang-format off
  vel_weight_t *vel_w = &vel_ws[8 * i + (field == v2_e) * 4];
  float *vf           = field == v1_e ? v1               : v2;
  float *v_prior      = field == v1_e ? v1_prior         : v2_prior;
  float *v_out        = field == v1_e ? &particles[i].v1 : &particles[i].v2;
  // clang-format on

  float v_pic = 0.f;
  float v_flip = 0.f;
  float w = 0.f;

  for (int j = 0; j < 4; ++j, ++vel_w) {
    if (vel_w->i < 0 || !vel_w->w || !isfinite(vf[vel_w->i])) continue;

    v_pic += vf[vel_w->i] * vel_w->w;
    v_flip += (vf[vel_w->i] - v_prior[vel_w->i]) * vel_w->w;
    w += vel_w->w;
  }

  v_flip = v_flip / w + *v_out;
  v_pic = v_pic / w;

  *v_out = lerp(v_pic, v_flip, flip);
}

// clang-format off
const SDL_Color c_wall  = {155, 155, 155, SDL_ALPHA_OPAQUE};
const SDL_Color c_fluid = {200, 220, 255, SDL_ALPHA_OPAQUE};
const SDL_Color c_air   = {255, 254, 255, SDL_ALPHA_OPAQUE};
// clang-format on

void render_cells(SDL_Renderer *renderer);
void render_particles(SDL_Renderer *renderer);
void render_velocities(SDL_Renderer *renderer);

void render_simulation(SDL_Renderer *renderer) {
  set_color(renderer, &c_air);
  SDL_RenderClear(renderer);

  render_cells(renderer);
  if (show_particles) { render_particles(renderer); }
  if (show_velocities) { render_velocities(renderer); }
}

int opacity(float density) {
  return iclamp(2 * 255 * density / PARTICLES_PER_CELL, 0, 255);
}

void render_cells(SDL_Renderer *renderer) {
  for (int i = 0; i < SIM_H; ++i) {
    for (int j = 0; j < SIM_W; ++j) {
      SDL_FRect rect = {WINDOW_SCALE * CELL_W * j, WINDOW_SCALE * CELL_H * i,
                        WINDOW_SCALE * CELL_W, WINDOW_SCALE * CELL_H};
      switch (states[i][j]) {
        case water_e:
          set_colora(renderer, &c_fluid, opacity(densities[i][j]));
          break;
        case solid_e:
          set_color(renderer, &c_wall);
          break;
        case air_e:
          set_color(renderer, &c_air);
          break;
      }
      SDL_RenderFillRect(renderer, &rect);
    }
  }
}

void render_particles(SDL_Renderer *renderer) {
  set_color(renderer, &(SDL_Color){4, 2, 0, SDL_ALPHA_OPAQUE});
  for (int i = 0; i < n_particles; ++i) {
    particle_t p = particles[i];
    SDL_RenderPoint(renderer, WINDOW_SCALE * p.x1, WINDOW_SCALE * p.x2);
  }
}

void render_velocities(SDL_Renderer *renderer) {
  const float scale = 0.1 * WINDOW_SCALE;
  set_color(renderer, &(SDL_Color){160, 255, 0, SDL_ALPHA_OPAQUE});

  // x velocities
  for (int i = 0; i < V1N; ++i) {
    int row = 0, col = 0;
    coordinates_from_index(i, SIM_W + 1, &row, &col);
    float x = WINDOW_SCALE * col * CELL_W;
    float y = WINDOW_SCALE * (row + 0.5) * CELL_H;
    SDL_RenderLine(renderer, x, y, x + v1[i] * scale, y);
  }

  // y velocities
  for (int i = 0; i < V2N; ++i) {
    int row = 0, col = 0;
    coordinates_from_index(i, SIM_W, &row, &col);
    float x = WINDOW_SCALE * (col + 0.5) * CELL_W;
    float y = WINDOW_SCALE * row * CELL_H;
    SDL_RenderLine(renderer, x, y, x, y + v2[i] * scale);
  }
}