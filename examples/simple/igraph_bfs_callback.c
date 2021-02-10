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

    /* Create a disjoint union of two rings */
    igraph_ring(&ring, 10, /*directed=*/ 0, /*mutual=*/ 0, /*circular=*/ 1);
    igraph_disjoint_union(&graph, &ring, &ring);
    igraph_destroy(&ring);

    /* Now call the BFS function */
    printf("(");
    igraph_bfs(&graph, /*root=*/0, /*roots=*/ 0, /*neimode=*/ IGRAPH_OUT,
               /*unreachable=*/ 1, /*restricted=*/ 0,
               /*order=*/ 0, /*rank=*/ 0, /*father=*/ 0, /*pred=*/ 0,
               /*succ=*/ 0, /*dist=*/ 0,
               /*callback=*/ bfs_callback, /*extra=*/ 0);
    printf(" )\n");

    /* Cleam up after ourselves */
    igraph_destroy(&graph);

    return 0;
}
