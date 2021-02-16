/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2009-2021  The igraph development team

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include <igraph.h>

#include "test_utilities.inc"

igraph_bool_t bfs_callback(const igraph_t *graph,
                           igraph_integer_t vid,
                           igraph_integer_t pred,
                           igraph_integer_t succ,
                           igraph_integer_t rank,
                           igraph_integer_t dist,
                           void *extra) {
    printf(" %li", (long int) vid);
    return 0;
}

int main() {

    igraph_t graph, ring;
    igraph_vector_t order, rank, father, pred, succ, dist;
    igraph_vector_t restricted;
    igraph_vector_t roots;
    long int i;

    igraph_ring(&ring, 10, /*directed=*/ 0, /*mutual=*/ 0, /*circular=*/ 1);
    igraph_disjoint_union(&graph, &ring, &ring);
    igraph_destroy(&ring);

    igraph_vector_init(&order, 0);
    igraph_vector_init(&rank, 0);
    igraph_vector_init(&father, 0);
    igraph_vector_init(&pred, 0);
    igraph_vector_init(&succ, 0);
    igraph_vector_init(&dist, 0);

    igraph_bfs(&graph, /*root=*/0, /*roots=*/ 0, /*neimode=*/ IGRAPH_OUT,
               /*unreachable=*/ 1, /*restricted=*/ 0,
               &order, &rank, &father, &pred, &succ, &dist,
               /*callback=*/ 0, /*extra=*/ 0);

    print_vector_round(&order);
    print_vector_round(&rank);
    print_vector_round(&father);
    print_vector_round(&pred);
    print_vector_round(&succ);
    print_vector_round(&dist);

    igraph_vector_destroy(&order);
    igraph_vector_destroy(&rank);
    igraph_vector_destroy(&father);
    igraph_vector_destroy(&pred);
    igraph_vector_destroy(&succ);
    igraph_vector_destroy(&dist);

    /* Test the callback */

    printf("(");
    igraph_bfs(&graph, /*root=*/ 0, /*roots=*/ 0, /*neimode=*/ IGRAPH_OUT,
               /*unreachable=*/ 1, /*restricted=*/ 0,
               0, 0, 0, 0, 0, 0, &bfs_callback, 0);
    printf(" )\n");

    /* Test different roots */

    printf("(");
    igraph_bfs(&graph, /*root=*/ 2, /*roots=*/ 0, /*neimode=*/ IGRAPH_OUT,
               /*unreachable=*/ 1, /*restricted=*/ 0,
               0, 0, 0, 0, 0, 0, &bfs_callback, 0);
    printf(" )\n");

    /* Test restricted */

    igraph_vector_init(&restricted, 0);
    for (i = 5; i < igraph_vcount(&graph); i++) {
        igraph_vector_push_back(&restricted, i);
    }
    printf("(");
    igraph_bfs(&graph, /*root=*/ 5, /*roots=*/ 0, /*neimode=*/ IGRAPH_OUT,
               /*unreachable=*/ 1, &restricted,
               0, 0, 0, 0, 0, 0, &bfs_callback, 0);
    printf(" )\n");

    /* Root not in restricted set */

    printf("(");
    igraph_bfs(&graph, /*root=*/ 4, /*roots=*/ 0, /*neimode=*/ IGRAPH_OUT,
               /*unreachable=*/ 1, &restricted,
               0, 0, 0, 0, 0, 0, &bfs_callback, 0);
    printf(" )\n");

    printf("(");
    igraph_bfs(&graph, /*root=*/ 3, /*roots=*/ 0, /*neimode=*/ IGRAPH_OUT,
               /*unreachable=*/ 0, &restricted,
               0, 0, 0, 0, 0, 0, &bfs_callback, 0);
    printf(" )\n");

    /* Multiple root vertices */

    igraph_vector_init(&roots, 3);
    VECTOR(roots)[0] = 3;
    VECTOR(roots)[1] = 4;
    VECTOR(roots)[2] = 6;
    printf("(");
    igraph_bfs(&graph, /*root=*/ -1, &roots, /*neimode=*/ IGRAPH_OUT,
               /*unreachable=*/ 0, &restricted,
               0, 0, 0, 0, 0, 0, &bfs_callback, 0);
    printf(" )\n");

    igraph_vector_destroy(&roots);
    igraph_vector_destroy(&restricted);
    igraph_destroy(&graph);

    VERIFY_FINALLY_STACK();

    return 0;
}
