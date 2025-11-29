#include "hash_grid.h"

#include <stdlib.h>

hash_grid_t hg_init(int lookup_size, int item_count) {
  return (hash_grid_t){.lookup = malloc(lookup_size * sizeof(int)),
                       .indices = malloc(item_count * sizeof(int)),
                       .n_lookup = lookup_size,
                       .n_indices = item_count};
}

void hg_free(hash_grid_t hg) {
  if (hg.lookup != nullptr) { free(hg.lookup); }
  if (hg.indices != nullptr) { free(hg.indices); }
}