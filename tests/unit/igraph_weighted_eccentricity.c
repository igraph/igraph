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

#include "test_utilities.h"

void print_and_destroy(igraph_t *graph, igraph_neimode_t mode) {
    igraph_vector_t ecc;
    igraph_vector_t weights;
    igraph_vector_init(&ecc, 0);
    igraph_vector_init(&weights, igraph_ecount(graph));
    for (int i = 0; i < igraph_ecount(graph); i++) {
        VECTOR(weights)[i] = 1;
    }
    igraph_weighted_eccentricity(graph, &weights, &ecc, igraph_vss_all(), mode);
    print_vector(&ecc);
    igraph_destroy(graph);
    igraph_vector_destroy(&ecc);
    igraph_vector_destroy(&weights);
}

int main() {
    igraph_t g;

    printf("Null graph:\n");
    igraph_empty(&g, 0, IGRAPH_UNDIRECTED);
    print_and_destroy(&g, IGRAPH_OUT);

    printf("\nSingleton graph:\n");
    igraph_empty(&g, 1, IGRAPH_UNDIRECTED);
    print_and_destroy(&g, IGRAPH_OUT);

    printf("\nPath with isolated vertex:\n");
    igraph_small(&g, 3, IGRAPH_UNDIRECTED,
                 0,2,
                 -1);
    print_and_destroy(&g, IGRAPH_OUT);

    printf("\nUndirected path graph:\n");
    igraph_ring(&g, 5, IGRAPH_UNDIRECTED, /* mutual */ 0, /* circular */ 0);
    print_and_destroy(&g, IGRAPH_OUT);

    printf("\nDirected path graph:\n");
    igraph_ring(&g, 5, IGRAPH_DIRECTED, /* mutual */ 0, /* circular */ 0);
    print_and_destroy(&g, IGRAPH_OUT);

    printf("\nUndirected star:\n");
    igraph_star(&g, 10, IGRAPH_STAR_UNDIRECTED, 0);
    print_and_destroy(&g, IGRAPH_OUT);

    printf("\nOut-star:\n");
    igraph_star(&g, 10, IGRAPH_STAR_OUT, 0);
    print_and_destroy(&g, IGRAPH_ALL);

    printf("\nIn-star:\n");
    igraph_star(&g, 10, IGRAPH_STAR_OUT, 0);
    print_and_destroy(&g, IGRAPH_OUT);

    VERIFY_FINALLY_STACK();

    return 0;
}
