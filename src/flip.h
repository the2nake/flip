#pragma once

#include "hash_grid.h"

//==================
//   PARAMETERS
//==================

constexpr int SIM_W = 60;
constexpr int SIM_H = 60;

constexpr int CELL_W = 2;
constexpr int CELL_H = CELL_W;

constexpr int PARTICLE_PACKING = 3;
constexpr int PARTICLES_PER_CELL = 9;

constexpr float PARTICLE_SIZE = (float)CELL_W / PARTICLE_PACKING;
constexpr float DENSITY_0 = PARTICLES_PER_CELL;

constexpr float BACKTRACK_RANGE = 2.f;
constexpr int BACKTRACK_ATTEMPTS = 10;

constexpr int WINDOW_W = SIM_W * CELL_W;
constexpr int WINDOW_H = SIM_H * CELL_H;
constexpr float WINDOW_SCALE = 6.0f;

//==================
//   MACROS
//==================

#define time_of(x)           \
  ({                         \
    float time_of_s = now(); \
    x;                       \
    float time_of_e = now(); \
    time_of_e - time_of_s;   \
  })

extern const float k_frametime;
extern const float k_timestep;
extern const float k_gravity;

extern const float k_relax;
extern const float k_flip;
extern const float k_stiffness;

extern const int k_iters;

typedef enum { solid_e = 0, water_e = 1, air_e = 2 } state_e_t;
typedef enum { v1_e, v2_e } field_e_t;

typedef struct {
  float x1;  // x position
  float x2;  // y position
  float v1;  // x velocity
  float v2;  // y velocity
} particle_t;

typedef struct {
  int i;
  int j;
} cell_t;

constexpr int V1N = (SIM_W + 1) * SIM_H;
constexpr int V2N = SIM_W * (SIM_H + 1);

extern state_e_t states[SIM_H][SIM_W];  // states
extern float densities[SIM_H][SIM_W];   // densities

extern float v1[V1N];  // horizontal velocity
extern float v2[V2N];  // vertical velocity
extern float w1[V1N];  // velocity field weights
extern float w2[V2N];  // velocity field weights
extern float v1_prior[V1N];
extern float v2_prior[V1N];

typedef struct {
  int i;
  float w;
} vel_weight_t;

extern particle_t *particles;
extern vel_weight_t *vel_ws;
extern int n_particles;
extern hash_grid_t particle_grid;

float now();
void print_field(float *arr, int w, int h, const char *name);
void print_v1();
void print_v2();
void print_w1();
void print_w2();
void inspect(particle_t* p);

bool in_rangei(int val, int lo, int hi);
bool in_rangef(float val, float lo, float hi);

int iclamp(int val, const int lo, const int hi);
float fclamp(float val, const float lo, const float hi);
float lerp(float a, float b, float t);
void normalise(float *v1, float *v2);

bool cell_in_bounds(int i, int j);
bool particle_in_bounds(particle_t *p);

void particle_enforce_bounds(particle_t *p);

int particle_v1_index(particle_t *p);  // index of top-left x vel
int particle_v2_index(particle_t *p);  // index of top-left y vel

void coordinates_from_index(int index, int row_width, int *row, int *col);

void position_in_v1_grid(particle_t *p, int row, int col, float *dx, float *dy);
void position_in_v2_grid(particle_t *p, int row, int col, float *dx, float *dy);

bool cell_on_edge(int i, int j);
bool cell_is(int i, int j, state_e_t state);
void get_cell_normal(int i, int j, float *v1, float *v2);
void set_cell_at(particle_t *particle, state_e_t state);

void get_particle_cell(particle_t *particle, int *i, int *j);
bool particle_in(particle_t *particle, state_e_t state);
float particle_velocity(particle_t *p);
