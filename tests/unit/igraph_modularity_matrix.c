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

void test_print_destroy(igraph_t *g, igraph_vector_t *weights, float resolution, igraph_matrix_t *modmat, igraph_bool_t directed) {
    int i, j;
    IGRAPH_ASSERT(igraph_modularity_matrix(g, weights, resolution, modmat, directed) == IGRAPH_SUCCESS);
    for (i = 0; i < igraph_matrix_nrow(modmat); i++) {
        for (j = 0; j < igraph_matrix_ncol(modmat); j++) {
            if (fabs(MATRIX(*modmat, i, j)) < 1e-10) {
                MATRIX(*modmat, i, j) = 0;
            }
        }
    }
    print_matrix(modmat);
    igraph_destroy(g);
    igraph_matrix_destroy(modmat);
    if (weights) {
        igraph_vector_destroy(weights);
    }
}

int main() {
    igraph_t g;
    igraph_vector_t weights, membership;
    igraph_matrix_t modmat;
    igraph_real_t modularity, test_modularity;
    int i, j;

    printf("No vertices:\n");
    igraph_small(&g, 0, /*directed*/0, -1);
    igraph_matrix_init(&modmat, 0, 0);
    test_print_destroy(&g, NULL, 1.0, &modmat, 0);

    printf("No edges:\n");
    igraph_small(&g, 3, /*directed*/0, -1);
    igraph_matrix_init(&modmat, 0, 0);
    test_print_destroy(&g, NULL, 1.0, &modmat, 0);

    printf("Triangle with no resolution should give the adjacency matrix:\n");
    igraph_small(&g, 3, /*directed*/0, 0,1, 0,2, 1,2, -1);
    igraph_matrix_init(&modmat, 0, 0);
    test_print_destroy(&g, NULL, 0.0, &modmat, 0);

    printf("Triangle and point with self-loop, undirected :\n");
    igraph_small(&g, 4, /*directed*/0, 0,1, 0,2, 1,2, 3,3, -1);
    igraph_matrix_init(&modmat, 0, 0);
    test_print_destroy(&g, NULL, 1.0, &modmat, 0);

    printf("Triangle and point with self-loop, directed, but direction ignored:\n");
    igraph_small(&g, 4, /*directed*/1, 0,1, 0,2, 1,2, 3,3, -1);
    igraph_matrix_init(&modmat, 0, 0);
    test_print_destroy(&g, NULL, 1.0, &modmat, 0);

    printf("Triangle and point with self-loop, directed:\n");
    igraph_small(&g, 4, /*directed*/1, 0,1, 0,2, 1,2, 3,3, -1);
    igraph_matrix_init(&modmat, 0, 0);
    test_print_destroy(&g, NULL, 1.0, &modmat, 1);

    printf("Triangle with weights 0, 1, 2:\n");
    igraph_small(&g, 3, /*directed*/0, 0,1, 0,2, 1,2, -1);
    igraph_vector_init_int(&weights, 3, 0, 1, 2);
    igraph_matrix_init(&modmat, 0, 0);
    test_print_destroy(&g, &weights, 1.0, &modmat, 0);

    printf("Triangle with weights 0, -1, -2:\n");
    igraph_small(&g, 3, /*directed*/0, 0,1, 0,2, 1,2, -1);
    igraph_vector_init_int(&weights, 3, 0, -1, -2);
    igraph_matrix_init(&modmat, 0, 0);
    test_print_destroy(&g, &weights, 1.0, &modmat, 0);

    printf("Directed triangle with weights 0, 1, 2:\n");
    igraph_small(&g, 3, /*directed*/1, 0,1, 0,2, 1,2, -1);
    igraph_vector_init_int(&weights, 3, 0, 1, 2);
    igraph_matrix_init(&modmat, 0, 0);
    test_print_destroy(&g, &weights, 1.0, &modmat, 1);

    printf("Triangle with weights -1, 0, 1 will cause divisions by zero:\n");
    igraph_small(&g, 3, /*directed*/0, 0,1, 0,2, 1,2, -1);
    igraph_vector_init_int(&weights, 3, -1, 0, 1);
    igraph_matrix_init(&modmat, 0, 0);
    test_print_destroy(&g, &weights, 1.0, &modmat, 0);

    printf("Comparison with modularity:\n");
    igraph_small(&g, 5, /*directed*/1, 0,1, 0,2, 1,2, 3,4, 4,0, -1);
    igraph_vector_init_int(&weights, 5, 1, 2, 3, 4, 5);
    igraph_vector_init_int(&membership, 5, 0, 0, 0, 1, 1);
    igraph_matrix_init(&modmat, 0, 0);
    IGRAPH_ASSERT(igraph_modularity_matrix(&g, &weights, 0.7, &modmat, 1) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_modularity(&g, &membership, &weights, 0.7, 1, &modularity) == IGRAPH_SUCCESS);
    print_matrix(&modmat);
    test_modularity = 0;
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            test_modularity += MATRIX(modmat, i, j);
        }
    }
    for (i = 3; i < 5; i++) {
        for (j = 3; j < 5; j++) {
            test_modularity += MATRIX(modmat, i, j);
        }
    }
    printf("Modularity: %g, modularity via matrix: %g\n", modularity, test_modularity / igraph_vector_sum(&weights));
    igraph_destroy(&g);
    igraph_vector_destroy(&membership);
    igraph_vector_destroy(&weights);
    igraph_matrix_destroy(&modmat);

    VERIFY_FINALLY_STACK();
    return 0;
}
