#include <SDL3/SDL.h>
#include <SDL3/SDL_main.h>
#include <SDL3/SDL_pixels.h>
#include <SDL3/SDL_render.h>
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

void advect();
void v_to_grid();
void project(int iters);
void v_to_particles();

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

    sum_t1 += time_of(advect());
    sum_t2 += time_of(v_to_grid());
    sum_t3 += time_of(project(k_iters));
    sum_t4 += time_of(v_to_particles());

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
      } else if (j * 2 < SIM_W && i > SIM_H / 5) {
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
  constexpr float left_pad = x_gap / 2.f;
  constexpr float top_pad = y_gap / 2.f;

  int idx = 0;

  for (int i = 0; i < SIM_H; ++i) {
    for (int j = 0; j < SIM_W; ++j) {
      if (!cell_is(i, j, water_e)) continue;

      for (int k = 0; k < PARTICLES_PER_CELL; ++k) {
        particles[idx] = (particle_t){
            .x1 = j * CELL_W + left_pad + x_gap * (k % DENSITY),
            .x2 = i * CELL_H + top_pad + y_gap * (int)(k / DENSITY),
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
      int j = p->x1 / CELL_W, i = p->x2 / CELL_H;
      float cell_x = (j + 0.5) * CELL_W;
      float cell_y = (i + 0.5) * CELL_H;

      int tries = 0;
      for (; tries < BACKTRACK && in_cell; ++tries) {
        p->x1 -= 2 * (dx1 > 0 ? 1 : -1) * fmax(fabs(dx1) / BACKTRACK, 1.5f);
        p->x2 -= 2 * (dx2 > 0 ? 1 : -1) * fmax(fabs(dx2) / BACKTRACK, 1.5f);
        in_cell = particle_in(p, solid_e);
      }

      // hacky way to avoid surface normals
      // if (fabs(i - SIM_W * 0.5) > fabs(j - SIM_H * 0.5)) {
      //   if ((p->v1 > 0) == (cell_x - p->x1 > 0)) { p->v1 *= -.95; }
      // } else {
      //   if ((p->v2 > 0) == (cell_y - p->x2 > 0)) { p->v2 *= -.95; }
      // }
    }
    set_cell_at(&particles[i], water_e);
  }
}

void compute_weights(particle_t *p, cell_weight_t *w);
void add_weight(field_e_t field, cell_weight_t *cell, float v);
void enforce_solid_velocity_field() {
}  // TODO! make sure solid block-caused velocities are not overridden

void v_to_grid() {
  particle_t *p = particles;
  reset_velocity_field();
  for (int i = 0; i < n_particles; ++i) {
    // somehow this looped list pairing thing is faster by 0.1 ms and
    // easier to read. vectorization? i don't even know

    compute_weights(&p[i], &particles_w[8 * i]);

    for (int j = 0; j < 4; ++j) {
      add_weight(v1_e, &particles_w[8 * i + j], p[i].v1);
      add_weight(v2_e, &particles_w[8 * i + j + 4], p[i].v2);
    }
  }

  for (int i = 0; i < V1N; ++i) {
    if (isfinite(v1[i]) && w1[i]) v1[i] /= w1[i];
  }
  for (int i = 0; i < V2N; ++i) {
    if (isfinite(v2[i]) && w2[i]) v2[i] /= w2[i];
  }

  enforce_solid_velocity_field();
  update_prior_velocities();
}

void compute_weights(particle_t *p, cell_weight_t *w) {
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

void add_weight(field_e_t field, cell_weight_t *c, float v) {
  if (!in_rangei(c->i, 0, field == v1_e ? V1N : V2N)) {
    c->w = 0.f;
    return;
  }

  float *vf = field == v1_e ? v1 : v2;
  float *wf = field == v1_e ? w1 : w2;
  int i = c->i / (SIM_W + (field == v1_e));
  int j = c->i % (SIM_W + (field == v1_e));

  bool is_water, touches_solid;

  if (field == v1_e) {
    is_water = cell_is(i, j - 1, water_e) || cell_is(i, j, water_e);
    touches_solid = cell_is(i, j - 1, solid_e) || cell_is(i, j, solid_e);
  } else {
    is_water = cell_is(i - 1, j, water_e) || cell_is(i, j, water_e);
    touches_solid = cell_is(i - 1, j, solid_e) || cell_is(i, j, solid_e);
  }

  if (is_water && !touches_solid) {
    vf[c->i] += c->w * v;
    wf[c->i] += c->w;
  } else {
    vf[c->i] = touches_solid ? 0.f : NAN;
    wf[c->i] = 0.f;
    c->w = 0.f;
  }
}

// enforce inflow = outflow iteratively
void project(int iters) {
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
  float pic = 0.02;

  for (int i = 0; i < n_particles; ++i) {
    update_particle(i, v1_e, k_flip);
    update_particle(i, v2_e, k_flip);
  }
}

void update_particle(int i, field_e_t field, float flip) {
  assert(0 <= i && i < n_particles);
  assert(-0.01 <= flip && flip <= 1.01);

  // clang-format off
  cell_weight_t *c = &particles_w[8 * i + (field == v2_e) * 4];
  float *vf         = field == v1_e ? v1               : v2;
  float *v_prior   = field == v1_e ? v1_prior         : v2_prior;
  float *v_out     = field == v1_e ? &particles[i].v1 : &particles[i].v2;
  // clang-format on

  float v_pic = 0.f;
  float v_flip = 0.f;
  float w = 0.f;

  for (int j = 0; j < 4; ++j, ++c) {
    if (c->i < 0 || !c->w || !isfinite(vf[c->i])) continue;

    v_pic += vf[c->i] * c->w;
    v_flip += (vf[c->i] - v_prior[c->i]) * c->w;
    w += c->w;
  }

  v_flip = v_flip / w + *v_out;
  v_pic = v_pic / w;

  *v_out = lerp(v_pic, v_flip, flip);
}

void render_cells(SDL_Renderer *renderer);
void render_particles(SDL_Renderer *renderer);
void render_velocities(SDL_Renderer *renderer);

void render_simulation(SDL_Renderer *renderer) {
  render_cells(renderer);
  render_particles(renderer);
  // render_velocities(renderer);
}

void set_color(SDL_Renderer *renderer, const SDL_Color *color) {
  SDL_SetRenderDrawColor(renderer, color->r, color->g, color->b, color->a);
}

void render_cells(SDL_Renderer *renderer) {
  // clang-format off
  const SDL_Color c_wall  = {155, 155, 155, SDL_ALPHA_OPAQUE};
  const SDL_Color c_fluid = {200, 220, 255, SDL_ALPHA_OPAQUE};
  const SDL_Color c_air   = {255, 255, 255, SDL_ALPHA_OPAQUE};
  // clang-format on

  for (int i = 0; i < SIM_H; ++i) {
    for (int j = 0; j < SIM_W; ++j) {
      SDL_FRect rect = {CELL_W * j, CELL_H * i, CELL_W, CELL_H};
      switch (states[i][j]) {
        case water_e:
          set_color(renderer, &c_fluid);
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
    SDL_RenderPoint(renderer, p.x1, p.x2);
  }
}

void render_velocities(SDL_Renderer *renderer) {
  const float scale = 0.25;
  set_color(renderer, &(SDL_Color){180, 255, 0, SDL_ALPHA_OPAQUE});

  // x velocities
  for (int i = 0; i < V1N; ++i) {
    int row = 0, col = 0;
    coordinates_from_index(i, SIM_W + 1, &row, &col);
    int x = col * CELL_W;
    int y = (row + 0.5) * CELL_H;
    SDL_RenderLine(renderer, x, y, x + v1[i] * scale, y);
  }

  // y velocities
  for (int i = 0; i < V2N; ++i) {
    int row = 0, col = 0;
    coordinates_from_index(i, SIM_W, &row, &col);
    int mid_x = (col + 0.5) * CELL_W;
    int mid_y = row * CELL_H;
    SDL_RenderLine(renderer, mid_x, mid_y, mid_x, mid_y + v2[i] * scale);
  }
}