/*
   IGraph library.
   Copyright (C) 2021  The igraph development team <igraph@igraph.org>

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
#include "test_utilities.inc"

void init_vm(igraph_vector_t *type_dist,
             int v0, int v1,
             igraph_matrix_t *pref_matrix,
             int m00, int m10, int m01, int m11) {
    igraph_vector_init_int_end(type_dist, -1, v0, v1, -1);
    igraph_matrix_init(pref_matrix, 2, 2);
    MATRIX(*pref_matrix, 0, 0) = m00;
    MATRIX(*pref_matrix, 1, 0) = m10;
    MATRIX(*pref_matrix, 0, 1) = m01;
    MATRIX(*pref_matrix, 1, 1) = m11;
}

#define DESTROY_GVM() do {               \
    igraph_destroy(&g);                  \
    igraph_vector_destroy(&type_dist);   \
    igraph_matrix_destroy(&pref_matrix); \
    } while(0)

int main() {
    igraph_t g;
    igraph_vector_t type_dist;
    igraph_vector_t node_type_vec;
    igraph_matrix_t pref_matrix;
    igraph_bool_t bipartite;

    igraph_rng_seed(igraph_rng_default(), 42);

    igraph_vector_init(&node_type_vec, 0);
    /*No vertices*/
    init_vm(&type_dist, 1, 0, &pref_matrix, 0, 0, 0, 1);
    IGRAPH_ASSERT(igraph_establishment_game(&g, /*nodes*/ 0, /*types*/ 2, /*edges_per_step*/ 0, &type_dist, &pref_matrix, /*directed*/ 0, /*node_type_vec*/ NULL) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_vcount(&g) == 0);
    IGRAPH_ASSERT(!igraph_is_directed(&g));
    DESTROY_GVM();

    /*Zero matrix elements for only possible vertex type means no edges*/
    init_vm(&type_dist, 1, 0, &pref_matrix, 0, 0, 0, 1);
    IGRAPH_ASSERT(igraph_establishment_game(&g, /*nodes*/ 20, /*types*/ 2, /*edges_per_step*/ 5, &type_dist, &pref_matrix, /*directed*/ 0, /*node_type_vec*/ NULL) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_ecount(&g) == 0);
    IGRAPH_ASSERT(igraph_vcount(&g) == 20);
    DESTROY_GVM();

    /*Two types with only cross terms makes a bipartite graph*/
    init_vm(&type_dist, 1, 1, &pref_matrix, 0, 1, 1, 0);
    IGRAPH_ASSERT(igraph_establishment_game(&g, /*nodes*/ 20, /*types*/ 2, /*edges_per_step*/ 5, &type_dist, &pref_matrix, /*directed*/ 1, &node_type_vec) == IGRAPH_SUCCESS);
    igraph_is_bipartite(&g, &bipartite, NULL);
    IGRAPH_ASSERT(bipartite);
    IGRAPH_ASSERT(igraph_is_directed(&g));
    IGRAPH_ASSERT(igraph_vector_min(&node_type_vec) == 0);
    IGRAPH_ASSERT(igraph_vector_max(&node_type_vec) == 1);
    IGRAPH_ASSERT(igraph_vector_size(&node_type_vec) == 20);
    DESTROY_GVM();

    /*Distribution of types should have at least one possible value*/
    igraph_vector_init(&type_dist, 0);
    igraph_matrix_init(&pref_matrix, 0, 0);
    igraph_set_error_handler(igraph_error_handler_ignore);
    IGRAPH_ASSERT(igraph_establishment_game(&g, /*nodes*/ 20, /* types*/ 2, /*edges_per_step*/ 5, &type_dist, &pref_matrix, /*directed*/ 0, /*node_type_vec*/ NULL) == IGRAPH_EINVAL);
    DESTROY_GVM();

    igraph_vector_destroy(&node_type_vec);

    VERIFY_FINALLY_STACK();
    return 0;
}
