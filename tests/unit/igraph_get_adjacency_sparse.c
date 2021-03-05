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

void print_and_destroy(igraph_t *g, igraph_get_adjacency_t type) {
    igraph_spmatrix_t result;
    igraph_spmatrix_init(&result, 0, 0);
    IGRAPH_ASSERT(igraph_get_adjacency_sparse(g, &result, type) == IGRAPH_SUCCESS);
    print_spmatrix(&result);
    igraph_spmatrix_destroy(&result);
    printf("\n");
}

int main() {
    igraph_t g_null, g_lm, g_lm_undir, g_empty;

    igraph_small(&g_null, 0, 0, -1);
    igraph_small(&g_empty, 3, 0, -1);
    igraph_small(&g_lm, 6, 1, 0,1, 0,2, 1,1, 1,3, 2,0, 2,3, 3,4, 3,4, -1);
    igraph_small(&g_lm_undir, 6, 0, 0,1, 0,2, 1,1, 1,3, 2,0, 2,3, 3,4, 3,4, -1);

    printf("No vertices:\n");
    print_and_destroy(&g_null, IGRAPH_GET_ADJACENCY_BOTH);

    printf("No edges:\n");
    print_and_destroy(&g_empty, IGRAPH_GET_ADJACENCY_BOTH);

    printf("Disconnected graph with loops and multiple edges:\n");
    print_and_destroy(&g_lm, IGRAPH_GET_ADJACENCY_BOTH);

    printf("Same graph, undirected, symmetric matrix:\n");
    print_and_destroy(&g_lm_undir, IGRAPH_GET_ADJACENCY_BOTH);

    printf("Same graph, undirected, upper triangular matrix:\n");
    print_and_destroy(&g_lm_undir, IGRAPH_GET_ADJACENCY_UPPER);

    printf("Same graph, undirected, lower triangular matrix:\n");
    print_and_destroy(&g_lm_undir, IGRAPH_GET_ADJACENCY_LOWER);

    igraph_destroy(&g_null);
    igraph_destroy(&g_empty);
    igraph_destroy(&g_lm);
    igraph_destroy(&g_lm_undir);

    VERIFY_FINALLY_STACK();
    return 0;
}
