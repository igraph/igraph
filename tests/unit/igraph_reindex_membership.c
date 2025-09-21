/*
   igraph library.
   Copyright (C) 2025  The igraph development team <igraph@igraph.org>

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include <igraph.h>

#include "test_utilities.h"

void check(igraph_int_t n, igraph_int_t type_index_count, igraph_int_t type_count) {
    igraph_vector_int_t membership;
    igraph_vector_int_t mapping;
    igraph_vector_int_t new_to_old;
    igraph_int_t nb_clusters;

    igraph_vector_int_init(&new_to_old, 0);

    igraph_vector_int_init(&mapping, type_count);
    for (igraph_int_t i=0; i < type_count; i++) {
        VECTOR(mapping)[i] = RNG_INTEGER(0, type_index_count-1);
    }

    igraph_vector_int_init(&membership, n);
    for (igraph_int_t i=0; i < n; i++) {
        VECTOR(membership)[i] = VECTOR(mapping)[ RNG_INTEGER(0, type_count-1) ];
    }

    igraph_int_t min_original_type = n > 0 ? igraph_vector_int_min(&membership) : 0;
    igraph_int_t max_original_type = n > 0 ? igraph_vector_int_max(&membership) : 0;

    printf("Old: "); print_vector_int(&membership);
    igraph_reindex_membership(&membership, &new_to_old, &nb_clusters);
    printf("New: "); print_vector_int(&membership);
    printf("New to old:  "); print_vector_int(&new_to_old);
    printf("nb_clusters: %" IGRAPH_PRId "\n\n", nb_clusters);

    IGRAPH_ASSERT(igraph_vector_int_size(&membership) == n);
    if (n == 0) {
        IGRAPH_ASSERT(nb_clusters == 0);
        IGRAPH_ASSERT(igraph_vector_int_size(&new_to_old) == 0);
    } else {
        IGRAPH_ASSERT(nb_clusters == igraph_vector_int_max(&membership) + 1);
        IGRAPH_ASSERT(nb_clusters <= type_count);
        IGRAPH_ASSERT(igraph_vector_int_min(&membership) == 0);

        IGRAPH_ASSERT(igraph_vector_int_size(&new_to_old) == nb_clusters);
        IGRAPH_ASSERT(igraph_vector_int_min(&new_to_old) == min_original_type);
        IGRAPH_ASSERT(igraph_vector_int_max(&new_to_old) == max_original_type);
    }

    igraph_vector_int_destroy(&membership);
    igraph_vector_int_destroy(&mapping);
    igraph_vector_int_destroy(&new_to_old);
}

int main(void) {
    igraph_rng_seed(igraph_rng_default(), 42);

    check(0, 10, 7);

    check(1, 1, 1);
    check(1, 10, 2);

    check(10, 10, 7);
    check(10, 100, 7);

    check(7, 7, 7);
    check(7, 100, 7);

    VERIFY_FINALLY_STACK();

    return 0;
}
