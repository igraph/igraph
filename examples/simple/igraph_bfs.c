/*
   igraph library.
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

int main(void) {

    igraph_t graph, ring;
    igraph_vector_int_t order, rank, father, pred, succ, dist;

    /* Initialize the library. */
    igraph_setup();

    /* Create a disjoint union of two rings */
    igraph_ring(&ring, 10, /*directed=*/ 0, /*mutual=*/ 0, /*circular=*/ 1);
    igraph_disjoint_union(&graph, &ring, &ring);
    igraph_destroy(&ring);

    /* Initialize the vectors where the result will be stored. Any of these
     * can be omitted and replaced with a null pointer when calling
     * igraph_bfs() */
    igraph_vector_int_init(&order, 0);
    igraph_vector_int_init(&rank, 0);
    igraph_vector_int_init(&father, 0);
    igraph_vector_int_init(&pred, 0);
    igraph_vector_int_init(&succ, 0);
    igraph_vector_int_init(&dist, 0);

    /* Now call the BFS function */
    igraph_bfs(&graph, /*root=*/0, /*roots=*/ NULL, /*neimode=*/ IGRAPH_OUT,
               /*unreachable=*/ 1, /*restricted=*/ NULL,
               &order, &rank, &father, &pred, &succ, &dist,
               /*callback=*/ NULL, /*extra=*/ NULL);

    /* Print the results */
    igraph_vector_int_print(&order);
    igraph_vector_int_print(&rank);
    igraph_vector_int_print(&father);
    igraph_vector_int_print(&pred);
    igraph_vector_int_print(&succ);
    igraph_vector_int_print(&dist);

    /* Cleam up after ourselves */
    igraph_vector_int_destroy(&order);
    igraph_vector_int_destroy(&rank);
    igraph_vector_int_destroy(&father);
    igraph_vector_int_destroy(&pred);
    igraph_vector_int_destroy(&succ);
    igraph_vector_int_destroy(&dist);

    igraph_destroy(&graph);

    return 0;
}
