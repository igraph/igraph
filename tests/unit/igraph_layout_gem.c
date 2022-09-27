/* IGraph library.
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
#include "test_utilities.h"

int main(void) {
    igraph_t graph;
    igraph_matrix_t res;
    igraph_bool_t use_seed;
    igraph_bool_t do_not_use_seed;
    igraph_integer_t maxiter;
    igraph_real_t temp_max;
    igraph_real_t temp_min;
    igraph_real_t temp_init;

    igraph_rng_seed(igraph_rng_default(), 42);
    use_seed = 1;
    do_not_use_seed = 0;
    maxiter = 10;
    temp_max = 1;
    temp_min = 0.1;
    temp_init = 1;

    printf("Check if 0 vertex graph crashes.\n");
    igraph_small(&graph, 0, IGRAPH_UNDIRECTED, -1);
    igraph_matrix_init(&res, 0, 2);

    igraph_layout_gem(&graph, &res, use_seed, maxiter, temp_max, temp_min, temp_init);
    igraph_matrix_print(&res);
    igraph_destroy(&graph);
    igraph_matrix_destroy(&res);

    printf("Check seeded 1-vertex graph, no iterations.\n");
    maxiter = 0;
    igraph_small(&graph, 1, IGRAPH_UNDIRECTED, -1);
    {
        int elem[] = {3,4};
        matrix_init_int_row_major(&res, 1, 2, elem);
    }

    igraph_layout_gem(&graph, &res, use_seed, maxiter, temp_max, temp_min, temp_init);
    igraph_matrix_print(&res);
    igraph_destroy(&graph);
    igraph_matrix_destroy(&res);
    maxiter = 40;

    printf("Check if 0 vertex graph crashes without a seed.\n");
    igraph_small(&graph, 0, IGRAPH_UNDIRECTED, -1);
    igraph_matrix_init(&res, 0, 2);

    igraph_layout_gem(&graph, &res, do_not_use_seed, maxiter, temp_max, temp_min, temp_init);
    igraph_matrix_print(&res);
    igraph_destroy(&graph);
    igraph_matrix_destroy(&res);

    printf("Check unseeded 1-vertex graph.\n");
    maxiter = 40;
    igraph_small(&graph, 1, IGRAPH_UNDIRECTED, -1);
    igraph_matrix_init(&res, 0, 2);

    igraph_layout_gem(&graph, &res, do_not_use_seed, maxiter, temp_max, temp_min, temp_init);
    igraph_destroy(&graph);
    igraph_matrix_destroy(&res);

    printf("Check 2-vertex graph.\n");
    maxiter = 100000;
    igraph_small(&graph, 2, IGRAPH_UNDIRECTED, 0,1, -1);
    {
        int elem[] = {3, 4, 5, 6};
        matrix_init_int_row_major(&res, 2, 2, elem);
    }

    igraph_layout_gem(&graph, &res, use_seed, maxiter, temp_max, temp_min, temp_init);
    IGRAPH_ASSERT(igraph_matrix_min(&res) > -1000);
    IGRAPH_ASSERT(igraph_matrix_max(&res) < 1000);
    igraph_destroy(&graph);

    VERIFY_FINALLY_STACK();

    maxiter = 3;
    printf("Check negative maxiter.\n");
    CHECK_ERROR(igraph_layout_gem(&graph, &res, use_seed, -1, temp_max, temp_min, temp_init), IGRAPH_EINVAL);
    printf("Check negative temp_max.\n");
    CHECK_ERROR(igraph_layout_gem(&graph, &res, use_seed, maxiter, -1, temp_min, temp_init), IGRAPH_EINVAL);
    printf("Check negative temp_min.\n");
    CHECK_ERROR(igraph_layout_gem(&graph, &res, use_seed, maxiter, temp_max, -1, temp_init), IGRAPH_EINVAL);
    printf("Check negative temp_init.\n");
    CHECK_ERROR(igraph_layout_gem(&graph, &res, use_seed, maxiter, temp_max, temp_min, -1), IGRAPH_EINVAL);
    printf("Check temp_init not between min and max.\n");
    CHECK_ERROR(igraph_layout_gem(&graph, &res, use_seed, maxiter, 1, 2, 1.5), IGRAPH_EINVAL);
    printf("Check too many rows in seed.\n");
    igraph_matrix_destroy(&res);
    igraph_matrix_init(&res, 10, 2);
    CHECK_ERROR(igraph_layout_gem(&graph, &res, use_seed, maxiter, temp_max, temp_min, temp_init), IGRAPH_EINVAL);
    printf("Check too many columns in seed.\n");
    igraph_matrix_destroy(&res);
    igraph_matrix_init(&res, 2, 5);
    CHECK_ERROR(igraph_layout_gem(&graph, &res, use_seed, maxiter, temp_max, temp_min, temp_init), IGRAPH_EINVAL);

    igraph_matrix_destroy(&res);
    igraph_destroy(&graph);

    VERIFY_FINALLY_STACK();
    return 0;
}
