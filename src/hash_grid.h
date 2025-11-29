#pragma once

typedef struct {
  int* lookup;
  int* indices;
  int n_lookup;
  int n_indices;
} hash_grid_t;

hash_grid_t hg_init(int lookup_size, int item_count);
void hg_compute(hash_grid_t* hg);
void hg_free(hash_grid_t* hg);
