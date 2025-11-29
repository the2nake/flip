#include "hash_grid.h"

#include <assert.h>
#include <stdlib.h>

#include "flip.h"

hash_grid_t hg_init(int lookup_size, int item_count) {
  hash_grid_t hg = {.lookup = nullptr,
                    .indices = nullptr,
                    .n_lookup = lookup_size + 1,  //  `+1` for guard cell at end
                    .n_indices = item_count};

  hg.lookup = malloc(hg.n_lookup * sizeof(int));
  hg.indices = malloc(hg.n_indices * sizeof(int));
  return hg;
}

void hg_compute(hash_grid_t* hg) {
  // clear hash grid
  for (int i = 0; i < hg->n_lookup; ++i) { hg->lookup[i] = 0; }
  for (int i = 0; i < hg->n_indices; ++i) { hg->indices[i] = -1; }

  const int cols = PARTICLE_PACKING * SIM_W;

  // count particles in each cell
  for (int i = 0; i < hg->n_indices; ++i) {
    int c_j = PARTICLE_PACKING * particles[i].x1 / CELL_W;
    int c_i = PARTICLE_PACKING * particles[i].x2 / CELL_H;
    hg->lookup[c_i * cols + c_j] += 1;
  }

  // compute prefix sum
  int sum = 0;
  for (int i = 0; i < hg->n_lookup; ++i) {
    int temp = sum;
    sum += hg->lookup[i];
    hg->lookup[i] = temp;
  }

  for (int i = 0; i < hg->n_indices; ++i) {
    int c_j = PARTICLE_PACKING * particles[i].x1 / CELL_W;
    int c_i = PARTICLE_PACKING * particles[i].x2 / CELL_H;

    int pos = hg->lookup[c_i * cols + c_j];
    assert(pos >= 0);  // otherwise the lookup was constructed incorrectly

    // look for an empty spot in `indices`
    while (pos < hg->n_indices && hg->indices[pos] != -1) { ++pos; }
    assert(pos < hg->n_indices);
    hg->indices[pos] = i;
  }
}

void hg_free(hash_grid_t* hg) {
  if (hg->lookup != nullptr) { free(hg->lookup); }
  if (hg->indices != nullptr) { free(hg->indices); }
}
