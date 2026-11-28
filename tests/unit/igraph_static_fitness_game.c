/*
   igraph library.
   Copyright (C) 2026  The igraph development team <igraph@igraph.org>

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

int main(void) {
    igraph_t graph;
    igraph_vector_t good, bad;
    const igraph_real_t bad_values[] = { -1.0, IGRAPH_NAN, IGRAPH_INFINITY };

    igraph_rng_seed(igraph_rng_default(), 42);

    igraph_vector_init_real(&good, 3, 1.0, 2.0, 1.0);
    igraph_vector_init_real(&bad, 3, 1.0, 2.0, 1.0);

    /* Valid input */
    igraph_static_fitness_game(&graph, 10, &good, &good, IGRAPH_LOOPS_SW | IGRAPH_MULTI_SW);
    IGRAPH_ASSERT(igraph_vcount(&graph) == 3);
    IGRAPH_ASSERT(igraph_ecount(&graph) == 10);
    igraph_destroy(&graph);

    /* Invalid fitness values are rejected in each fitness vector. */
    for (igraph_int_t i = 0; i < (igraph_int_t) (sizeof(bad_values) / sizeof(bad_values[0])); i++) {
        VECTOR(bad)[1] = bad_values[i];

        /* Request a single edge so that the "too many edges" check cannot fail first. */
        CHECK_ERROR(igraph_static_fitness_game(&graph, 1, &bad, NULL, IGRAPH_SIMPLE_SW), IGRAPH_EINVAL);
        CHECK_ERROR(igraph_static_fitness_game(&graph, 10, &bad, NULL, IGRAPH_LOOPS_SW | IGRAPH_MULTI_SW), IGRAPH_EINVAL);
        CHECK_ERROR(igraph_static_fitness_game(&graph, 10, &bad, &good, IGRAPH_LOOPS_SW | IGRAPH_MULTI_SW), IGRAPH_EINVAL);
        CHECK_ERROR(igraph_static_fitness_game(&graph, 10, &good, &bad, IGRAPH_LOOPS_SW | IGRAPH_MULTI_SW), IGRAPH_EINVAL);
    }

    igraph_vector_destroy(&good);
    igraph_vector_destroy(&bad);

    VERIFY_FINALLY_STACK();
    return 0;
}
