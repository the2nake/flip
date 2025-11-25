#include "flip.h"

#include <math.h>
#include <stdio.h>
#include <time.h>

const float k_frametime = 1.0 / 60.0;
const float k_timestep = 1.0 / 60.0;
const float k_gravity = 9.81 * SIM_H;
const float k_relax = 1.9;
const float k_flip = 0.93;
const int k_iters = 100;

state_e_t states[SIM_H][SIM_W];  // states
float densities[SIM_H][SIM_W];   // densities

float v1[V1N];  // horizontal velocity
float v2[V2N];  // vertical velocity
float w1[V1N];  // velocity field weights
float w2[V2N];  // velocity field weights
float v1_prior[V1N];
float v2_prior[V1N];

particle_t *particles;
vel_weight_t *vel_ws;
int n_particles = MAX_PARTICLES;

//=============
//    DEBUG
//=============

void print_field(float *arr, int w, int h, const char *name) {
  printf("-- %s --\n", name);
  for (int i = 0; i < h; ++i) {
    for (int j = 0; j < w; ++j) { printf("%4.1f ", arr[w * i + j]); }
    printf("\n");
  }
  printf("---\n");
}

void print_v1() { print_field(v1, SIM_W + 1, SIM_H, "x vel"); }
void print_v2() { print_field(v2, SIM_W, SIM_H + 1, "y vel"); }

void print_w1() { print_field(w1, SIM_W + 1, SIM_H, "x weights"); }
void print_w2() { print_field(w2, SIM_W, SIM_H + 1, "y weights"); }

//=============
//    UTILS
//=============

float now() { return (float)clock() / CLOCKS_PER_SEC; }

bool in_rangei(int val, int lo, int hi) { return lo <= val && val <= hi; }
bool in_rangef(float val, float lo, float hi) { return lo <= val && val <= hi; }

int iclamp(int val, const int lo, const int hi) {
  if (val > hi) { return hi; }
  if (val < lo) { return lo; }
  return val;
}

float fclamp(float val, const float lo, const float hi) {
  if (val > hi) { return hi; }
  if (val < lo) { return lo; }
  return val;
}

float lerp(float a, float b, float t) { return a + t * (b - a); }

bool cell_in_bounds(int i, int j) {
  return 0 <= i && i < SIM_H && 0 <= j && j < SIM_W;
}

bool particle_in_bounds(particle_t *p) {
  return 0.f <= p->x1 && p->x1 <= SIM_W * CELL_W && 0.f <= p->x2 &&
         p->x2 <= SIM_H * CELL_H;
}

void particle_enforce_bounds(particle_t *p) {
  if (!particle_in_bounds(p)) {
    p->x1 = fclamp(p->x1, 0.f, SIM_W * CELL_W);
    p->x2 = fclamp(p->x2, 0.f, SIM_H * CELL_H);
  }
}

// index of top-left x vel
int particle_v1_index(particle_t *p) {
  int v1_col = p->x1 / CELL_W;
  int v1_row = p->x2 / CELL_H - 0.5;
  return v1_row * (SIM_W + 1) + v1_col;
}

// index of top-left y vel
int particle_v2_index(particle_t *p) {
  int v2_col = p->x1 / CELL_W - 0.5;
  int v2_row = p->x2 / CELL_H;
  return v2_row * SIM_W + v2_col;
}

void coordinates_from_index(int index, int row_width, int *row, int *col) {
  *row = index / row_width;
  *col = index % row_width;
}

void position_in_v1_grid(particle_t *p, int row, int col, float *dx,
                         float *dy) {
  *dx = p->x1 - col * CELL_W;
  *dy = p->x2 - (row + 0.5) * CELL_H;
}

void position_in_v2_grid(particle_t *p, int row, int col, float *dx,
                         float *dy) {
  *dx = p->x1 - (col + 0.5) * CELL_W;
  *dy = p->x2 - row * CELL_H;
}

bool cell_on_edge(int i, int j) {
  if (i == 0 || i == SIM_H - 1) { return true; }
  if (j == 0 || j == SIM_W - 1) { return true; }
  return false;
}

bool cell_is(int i, int j, state_e_t state) {
  if (cell_in_bounds(i, j)) { return states[i][j] == state; }
  return false;
}

void set_cell_at(particle_t *particle, state_e_t state) {
  int j = particle->x1 / CELL_W;
  int i = particle->x2 / CELL_H;
  if (cell_in_bounds(i, j)) { states[i][j] = state; }
}

bool particle_in(particle_t *particle, state_e_t state) {
  int j = particle->x1 / CELL_W;
  int i = particle->x2 / CELL_H;
  return cell_is(i, j, state);  // check bounds to not lose particles
}

float particle_velocity(particle_t *p) { return hypot(p->v1, p->v2); }