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

void check_and_destroy(igraph_matrix_t *result, igraph_real_t half_size) {
    igraph_real_t min, max;
    igraph_matrix_minmax(result, &min, &max);
    IGRAPH_ASSERT(min >= -half_size);
    IGRAPH_ASSERT(max <= half_size);
    igraph_matrix_destroy(result);
}

int main() {
    igraph_t g;
    igraph_matrix_t result;

    igraph_rng_seed(igraph_rng_default(), 42);

    printf("No vertices:\n");
    igraph_small(&g, 0, 0, -1);
    igraph_matrix_init(&result, 0, 0);
    IGRAPH_ASSERT(igraph_layout_graphopt(&g, &result, /*niter*/ 500, /*node_charge*/ 0.001, /*node_mass*/ 30, /*spring_length*/ 0, /*spring constant*/ 1, /*max_sa_movement*/ 5, /*use seed*/ 0) == IGRAPH_SUCCESS);
    print_matrix(&result);
    igraph_matrix_destroy(&result);
    igraph_destroy(&g);

    printf("One vertex.\n");
    igraph_small(&g, 1, 0, -1);
    igraph_matrix_init(&result, 0, 0);
    IGRAPH_ASSERT(igraph_layout_graphopt(&g, &result, /*niter*/ 500, /*node_charge*/ 0.001, /*node_mass*/ 30, /*spring_length*/ 0, /*spring constant*/ 1, /*max_sa_movement*/ 5, /*use seed*/ 0) == IGRAPH_SUCCESS);
    check_and_destroy(&result, 1.0);
    igraph_destroy(&g);

    printf("Full graph of 4 vertices, no loops.\n");
    igraph_full(&g, 4, 0, 0);
    igraph_matrix_init(&result, 0, 0);
    IGRAPH_ASSERT(igraph_layout_graphopt(&g, &result, /*niter*/ 500, /*node_charge*/ 0.001, /*node_mass*/ 30, /*spring_length*/ 0, /*spring constant*/ 1, /*max_sa_movement*/ 5, /*use seed*/ 0) == IGRAPH_SUCCESS);
    check_and_destroy(&result, 20.0);
    igraph_destroy(&g);

    printf("4 vertices, disconnected, with loops and multi-edges.\n");
    igraph_small(&g, 4, 0, 0,0, 0,0, 0,0, 1,2, 1,2, 1,2, -1);
    igraph_matrix_init(&result, 0, 0);
    IGRAPH_ASSERT(igraph_layout_graphopt(&g, &result, /*niter*/ 500, /*node_charge*/ 0.001, /*node_mass*/ 30, /*spring_length*/ 0, /*spring constant*/ 1, /*max_sa_movement*/ 5, /*use seed*/ 0) == IGRAPH_SUCCESS);
    check_and_destroy(&result, 100.0);
    igraph_destroy(&g);

    printf("Full graph of 4 vertices, no loops with no repulsion.\n");
    igraph_full(&g, 4, 0, 0);
    igraph_matrix_init(&result, 0, 0);
    IGRAPH_ASSERT(igraph_layout_graphopt(&g, &result, /*niter*/ 500, /*node_charge*/ 0.0, /*node_mass*/ 30, /*spring_length*/ 0, /*spring constant*/ 1, /*max_sa_movement*/ 5, /*use seed*/ 0) == IGRAPH_SUCCESS);
    check_and_destroy(&result, 1.0);
    igraph_destroy(&g);

    printf("4 vertices in a line, with no repulsion, spring length 1 and a seed:\n");
    igraph_small(&g, 4, 0, 0,1, 1,2, 2,3, -1);
    igraph_real_t seed[] = {0.15, -0.15, 0.05, -0.05, -0.05, 0.05, -0.15, 0.15};
    matrix_init_real_row_major(&result, 4, 2, seed);
    IGRAPH_ASSERT(igraph_layout_graphopt(&g, &result, /*niter*/ 500, /*node_charge*/ 0.0, /*node_mass*/ 30, /*spring_length*/ 1, /*spring constant*/ 10, /*max_sa_movement*/ 5, /*use seed*/ 1) == IGRAPH_SUCCESS);
    print_matrix(&result);
    igraph_matrix_destroy(&result);
    igraph_destroy(&g);

    printf("4 vertices in a line, with no repulsion, spring length 1 and a seed, no sa_movement:\n");
    igraph_small(&g, 4, 0, 0,1, 1,2, 2,3, -1);
    matrix_init_real_row_major(&result, 4, 2, seed);
    IGRAPH_ASSERT(igraph_layout_graphopt(&g, &result, /*niter*/ 500, /*node_charge*/ 0.0, /*node_mass*/ 30, /*spring_length*/ 1, /*spring constant*/ 10, /*max_sa_movement*/ 0, /*use seed*/ 1) == IGRAPH_SUCCESS);
    print_matrix(&result);
    igraph_matrix_destroy(&result);
    igraph_destroy(&g);

    printf("Wrong size seed.\n");
    igraph_small(&g, 4, 0, 0,1, 1,2, 2,3, -1);
    matrix_init_real_row_major(&result, 3, 2, seed);
    IGRAPH_ASSERT(igraph_layout_graphopt(&g, &result, /*niter*/ 500, /*node_charge*/ 0.0, /*node_mass*/ 30, /*spring_length*/ 1, /*spring constant*/ 10, /*max_sa_movement*/ 0, /*use seed*/ 1) == IGRAPH_SUCCESS);
    check_and_destroy(&result, 1.0);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();
    return 0;
}
