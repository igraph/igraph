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
    IGRAPH_ASSERT(igraph_modularity_matrix(g, weights, resolution, modmat, directed) == IGRAPH_SUCCESS);
    print_matrix(modmat);
    igraph_destroy(g);
    igraph_matrix_destroy(modmat);
    if (weights) {
        igraph_vector_destroy(weights);
    }
}

int main() {
    igraph_t g;
    igraph_vector_t weights;
    igraph_matrix_t modmat;

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

    VERIFY_FINALLY_STACK();
    return 0;
}
