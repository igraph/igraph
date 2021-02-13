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

void make_box(int vertices, float half_size, igraph_vector_t *bounds) {
    for (int i = 0; i < 4; i++) {
        igraph_vector_init(&bounds[i], vertices);
    }
    igraph_vector_fill(&bounds[0], -half_size);
    igraph_vector_fill(&bounds[1], half_size);
    igraph_vector_fill(&bounds[2], -half_size);
    igraph_vector_fill(&bounds[3], half_size);
}

void destroy_bounds(igraph_vector_t *bounds) {
    for (int i = 0; i < 4; i++) {
        igraph_vector_destroy(&bounds[i]);
    }
}

int main() {
    igraph_t g;
    igraph_matrix_t result;
    igraph_vector_t bounds[4];
    igraph_vector_t weights;
    igraph_rng_seed(igraph_rng_default(), 42);

    printf("Empty graph:\n");
    igraph_small(&g, 0, 0, -1);
    igraph_matrix_init(&result, 0, 0);
    IGRAPH_ASSERT(igraph_layout_fruchterman_reingold(&g, &result, /*use_seed*/ 0,
                  /*niter*/ 0, /*start_temp*/ 0.0, IGRAPH_LAYOUT_NOGRID,
                  /*weight*/ NULL, /*minx*/ NULL, /*maxx*/ NULL, /*miny*/ NULL,
                  /*maxy*/ NULL) == IGRAPH_SUCCESS);
    print_matrix(&result);
    igraph_matrix_destroy(&result);
    igraph_destroy(&g);

    printf("Full graph of 10 vertices with loops and bounds:\n");
    igraph_full(&g, 10, 0, 1);
    igraph_matrix_init(&result, 0, 0);
    make_box(10, 1.0, bounds);
    IGRAPH_ASSERT(igraph_layout_fruchterman_reingold(&g, &result, /*use_seed*/ 0,
                  /*niter*/ 100, /*start_temp*/ 1.0, IGRAPH_LAYOUT_AUTOGRID,
                  /*weight*/ NULL, &bounds[0], &bounds[1], &bounds[2], &bounds[3]) == IGRAPH_SUCCESS);
    print_matrix(&result);
    igraph_matrix_destroy(&result);
    destroy_bounds(bounds);

    printf("Full graph of 10 vertices with loops and weights:\n");
    igraph_matrix_init(&result, 0, 0);
    igraph_vector_init(&weights, 55);
    igraph_vector_fill(&weights, 100);
    IGRAPH_ASSERT(igraph_layout_fruchterman_reingold(&g, &result, /*use_seed*/ 0,
                  /*niter*/ 100, /*start_temp*/ 1.0, IGRAPH_LAYOUT_AUTOGRID,
                  &weights, /*minx*/ NULL, /*maxx*/ NULL, /*miny*/ NULL,
                  /*maxy*/ NULL) == IGRAPH_SUCCESS);
    print_matrix(&result);
    igraph_matrix_destroy(&result);
    igraph_destroy(&g);
    igraph_vector_destroy(&weights);

    printf("Full graph of 5 vertices, seed and no iterations:\n");
    igraph_full(&g, 5, 0, 0);
    igraph_real_t seed[10] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
    matrix_init_real_row_major(&result, 5, 2, seed);
    IGRAPH_ASSERT(igraph_layout_fruchterman_reingold(&g, &result, /*use_seed*/ 1,
                  /*niter*/ 0, /*start_temp*/ 100, IGRAPH_LAYOUT_GRID,
                  /*weight*/ NULL, /*minx*/ NULL, /*maxx*/ NULL, /*miny*/ NULL,
                  /*maxy*/ NULL) == IGRAPH_SUCCESS);
    print_matrix(&result);
    igraph_matrix_destroy(&result);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();
    return 0;
}
