/*
   igraph library.
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

void print_and_destroy_weighted(igraph_t *graph, igraph_neimode_t mode, igraph_vector_t *weights) {
    igraph_vector_t ecc;
    igraph_vector_init(&ecc, 0);
    igraph_eccentricity(graph, weights, &ecc, igraph_vss_all(), mode);
    print_vector(&ecc);
    igraph_destroy(graph);
    igraph_vector_destroy(&ecc);
    if (weights) {
        igraph_vector_destroy(weights);
    }
}

void print_and_destroy(igraph_t *graph, igraph_neimode_t mode) {
    igraph_vector_t weights;
    igraph_vector_init(&weights, igraph_ecount(graph));
    for (int i = 0; i < igraph_ecount(graph); i++) {
        VECTOR(weights)[i] = 1;
    }
    print_and_destroy_weighted(graph, mode, &weights);
}

int main(void) {
    igraph_t g;
    igraph_vector_t weights;
    igraph_vector_t ecc;
    igraph_vector_init(&ecc, 0);

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

    printf("\nOut-star, IGRAPH_OUT:\n");
    igraph_star(&g, 10, IGRAPH_STAR_OUT, 0);
    print_and_destroy(&g, IGRAPH_OUT);

    printf("\nOut-star with NULL weights:\n");
    igraph_star(&g, 10, IGRAPH_STAR_OUT, 0);
    print_and_destroy_weighted(&g, IGRAPH_ALL, NULL);

    printf("\nOut-star with weights:\n");
    igraph_vector_init_real(&weights, 9, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9);
    igraph_small(&g, 10, IGRAPH_DIRECTED, 0,1, 0,2, 0,3, 0,4, 0,5, 0,6, 0,7, 0,8, 0,9, -1);
    print_and_destroy_weighted(&g, IGRAPH_ALL, &weights);

    VERIFY_FINALLY_STACK();

    printf("\nCheck wrong number of weights error.\n");
    igraph_vector_init(&weights, 1);
    igraph_empty(&g, 1, IGRAPH_UNDIRECTED);
    CHECK_ERROR(igraph_eccentricity(&g, &weights, &ecc, igraph_vss_all(), IGRAPH_OUT), IGRAPH_EINVAL);

    printf("Check NaN weight error.\n");
    igraph_destroy(&g);
    igraph_small(&g, 2, IGRAPH_UNDIRECTED, 0,1, -1);
    VECTOR(weights)[0] = IGRAPH_NAN;
    CHECK_ERROR(igraph_eccentricity(&g, &weights, &ecc, igraph_vss_all(), IGRAPH_OUT), IGRAPH_EINVAL);

    printf("Check negative weight error.\n");
    VECTOR(weights)[0] = -1;
    CHECK_ERROR(igraph_eccentricity(&g, &weights, &ecc, igraph_vss_all(), IGRAPH_OUT), IGRAPH_EINVAL);

    igraph_vector_destroy(&ecc);
    igraph_vector_destroy(&weights);
    igraph_destroy(&g);
    VERIFY_FINALLY_STACK();

    return 0;
}
