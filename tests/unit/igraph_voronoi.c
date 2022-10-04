/*
   IGraph library.
   Copyright (C) 2022  The igraph development team <igraph@igraph.org>

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
    igraph_t g;
    igraph_vector_int_t generators, membership;
    igraph_vector_t distances, weights;

    igraph_rng_seed(igraph_rng_default(), 42);

    /* Init result vars */

    igraph_vector_int_init(&membership, 0);
    igraph_vector_init(&distances, 0);

    /* Edge cases */

    printf("Singleton graph.\n");

    igraph_empty(&g, 1, IGRAPH_DIRECTED);
    igraph_vector_int_init_int(&generators, 1,
                               0);

    igraph_voronoi(&g, &membership, &distances, &generators, NULL, IGRAPH_ALL, IGRAPH_VORONOI_RANDOM);
    print_vector_int(&membership);
    print_vector(&distances);

    igraph_vector_int_destroy(&generators);
    igraph_destroy(&g);

    /* Disconnected directed multigraph */

    printf("\n\nDisconnected directed multigraph.\n");

    igraph_small(&g, 0, IGRAPH_DIRECTED,
                 0,2, 1,2, 2,3, 3,4, 5,4, 6,4, 2,3, 1,1,
                 -1);
    igraph_vector_int_init_int(&generators, 2,
                               0, 1);

    printf("\nTiebreaking: 'first'\n");
    igraph_voronoi(&g, &membership, &distances, &generators, NULL, IGRAPH_OUT, IGRAPH_VORONOI_FIRST);
    print_vector_int(&membership);
    print_vector(&distances);

    printf("\nTiebreaking: 'last'\n");
    igraph_voronoi(&g, &membership, &distances, &generators, NULL, IGRAPH_OUT, IGRAPH_VORONOI_LAST);
    print_vector_int(&membership);
    print_vector(&distances);

    printf("\nTiebreaking: 'random'\n");
    igraph_voronoi(&g, &membership, &distances, &generators, NULL, IGRAPH_OUT, IGRAPH_VORONOI_RANDOM);
    print_vector_int(&membership);
    print_vector(&distances);

    igraph_vector_int_destroy(&generators);
    igraph_destroy(&g);

    /* Karate club network */

    printf("\n\nKarate club, unweighted.\n");

    igraph_famous(&g, "Zachary");

    igraph_vector_int_init_int(&generators, 3,
                               0, 32, 24);

    printf("\nTiebreaking: 'first'\n");
    igraph_voronoi(&g, &membership, &distances, &generators, NULL, IGRAPH_ALL, IGRAPH_VORONOI_FIRST);
    print_vector_int(&membership);
    print_vector(&distances);

    printf("\nTiebreaking: 'last'\n");
    igraph_voronoi(&g, &membership, &distances, &generators, NULL, IGRAPH_ALL, IGRAPH_VORONOI_LAST);
    print_vector_int(&membership);
    print_vector(&distances);

    printf("\nTiebreaking: 'random'\n");
    igraph_voronoi(&g, &membership, &distances, &generators, NULL, IGRAPH_ALL, IGRAPH_VORONOI_RANDOM);
    print_vector_int(&membership);
    print_vector(&distances);

    printf("\n\nKarate club, betweenness weighted.\n");

    igraph_vector_init(&weights, 0);
    igraph_edge_betweenness(&g, &weights, IGRAPH_UNDIRECTED, NULL);

    printf("\nTiebreaking: 'first'\n");
    igraph_voronoi(&g, &membership, &distances, &generators, &weights, IGRAPH_ALL, IGRAPH_VORONOI_FIRST);
    print_vector_int(&membership);
    print_vector(&distances);

    printf("\nTiebreaking: 'last'\n");
    igraph_voronoi(&g, &membership, &distances, &generators, &weights, IGRAPH_ALL, IGRAPH_VORONOI_LAST);
    print_vector_int(&membership);
    print_vector(&distances);

    printf("\nTiebreaking: 'random'\n");
    igraph_voronoi(&g, &membership, &distances, &generators, &weights, IGRAPH_ALL, IGRAPH_VORONOI_RANDOM);
    print_vector_int(&membership);
    print_vector(&distances);

    igraph_vector_destroy(&weights);

    igraph_vector_int_destroy(&generators);
    igraph_destroy(&g);

    /* Destroy result vars */

    igraph_vector_destroy(&distances);
    igraph_vector_int_destroy(&membership);

    VERIFY_FINALLY_STACK();

    return 0;
}
