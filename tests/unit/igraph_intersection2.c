/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2013  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard st, Cambridge MA, 02139 USA

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
#include <stdio.h>

#include "test_utilities.inc"

int main() {

    igraph_t star, ring, uni, result;
    igraph_vector_ptr_t glist;

    igraph_star(&star, 11, IGRAPH_STAR_UNDIRECTED, /*center=*/ 10);
    igraph_ring(&ring, 10, IGRAPH_UNDIRECTED, /*mutual=*/ 0, /*circular=*/ 1);

    igraph_union(&uni, &star, &ring, /*edge_map1=*/ 0, /*edge_map2=*/ 0);

    igraph_intersection(&result, &uni, &star, /*edge_map1*/ 0,
                        /*edge_map2=*/ 0);
    igraph_write_graph_edgelist(&result, stdout);

    igraph_destroy(&result);

    /* ---------------------------- */

    igraph_vector_ptr_init(&glist, 2);
    VECTOR(glist)[0] = &uni;
    VECTOR(glist)[1] = &star;

    igraph_intersection_many(&result, &glist, /*edgemaps=*/ 0);
    printf("--\n");
    igraph_write_graph_edgelist(&result, stdout);

    igraph_vector_ptr_destroy(&glist);
    igraph_destroy(&result);
    igraph_destroy(&uni);
    igraph_destroy(&ring);
    igraph_destroy(&star);

    VERIFY_FINALLY_STACK();

    return 0;
}
