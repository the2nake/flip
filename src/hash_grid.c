#include "hash_grid.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "flip.h"

hash_grid_t hg_init(int lookup_size, int item_count) {
  return (hash_grid_t){.lookup = malloc(lookup_size * sizeof(int)),
                       .indices = malloc(item_count * sizeof(int)),
                       .n_lookup = lookup_size,
                       .n_indices = item_count};
}

void hg_compute(hash_grid_t* hg) {
  // clear hash grid
  for (int i = 0; i < hg->n_lookup; ++i) { hg->lookup[i] = 0; }
  for (int i = 0; i < hg->n_indices; ++i) { hg->indices[i] = -1; }

  // count particles in each cell
  for (int i = 0; i < hg->n_indices; ++i) {
    int c_i = 0, c_j = 0;
    get_particle_cell(&particles[i], &c_i, &c_j);
    if (isnan(particles[i].v1) || isnan(particles[i].v2)) {
      inspect(&particles[i]);
    }
    hg->lookup[c_i * SIM_W + c_j] += 1;
  }

  // compute prefix sum
  int sum = 0;
  for (int i = 0; i < hg->n_lookup; ++i) {
    int temp = sum;
    sum += hg->lookup[i];
    hg->lookup[i] = temp;
  }

  for (int i = 0; i < hg->n_indices; ++i) {
    int c_i = 0, c_j = 0;
    get_particle_cell(&particles[i], &c_i, &c_j);

    int pos = hg->lookup[c_i * SIM_W + c_j];
    const int fpos = pos;
    assert(pos >= 0);  // otherwise the lookup was constructed incorrectly

    // look for an empty spot in `indices`
    while (pos < hg->n_indices && hg->indices[pos] != -1) { ++pos; }
    if (pos >= hg->n_indices) {
      printf("%d, %d, (%d, %d)\n", i, fpos, c_i, c_j);
      printf("%f %f\n", particles[i].x1, particles[i].x2);
    }
    assert(pos < hg->n_indices);
    hg->indices[pos] = i;
  }
}

void hg_free(hash_grid_t* hg) {
  if (hg->lookup != nullptr) { free(hg->lookup); }
  if (hg->indices != nullptr) { free(hg->indices); }
}