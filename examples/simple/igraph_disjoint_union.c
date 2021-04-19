/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2006-2020  The igraph development team <igraph@igraph.org>

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
#include <stdio.h>

int main() {
    igraph_t left, right, uni;
    igraph_vector_ptr_t glist;
    long int i, n;

    igraph_small(&left, 4, IGRAPH_UNDIRECTED, 0,1, 1,2, 2,2, 2,3, -1);
    igraph_small(&right, 5, IGRAPH_UNDIRECTED, 0,1, 1,2, 2,2, 2,4, -1);

    igraph_disjoint_union(&uni, &left, &right);
    igraph_write_graph_edgelist(&uni, stdout);
    printf("\n");

    igraph_destroy(&left);
    igraph_destroy(&right);
    igraph_destroy(&uni);

    /* Empty graph list; the result is the directed null graph. */
    igraph_vector_ptr_init(&glist, 0);
    igraph_disjoint_union_many(&uni, &glist);
    if (!igraph_is_directed(&uni) || igraph_vcount(&uni) != 0) {
        return 1;
    }
    igraph_vector_ptr_destroy(&glist);
    igraph_destroy(&uni);

    /* Non-empty graph list. */
    igraph_vector_ptr_init(&glist, 10);
    n = igraph_vector_ptr_size(&glist);
    for (i = 0; i < n; i++) {
        VECTOR(glist)[i] = calloc(1, sizeof(igraph_t));
        igraph_small(VECTOR(glist)[i], 2, IGRAPH_DIRECTED, 0,1, 1,0, -1);
    }
    if (!igraph_is_directed(&uni)) {
        return 2;
    }

    igraph_disjoint_union_many(&uni, &glist);
    igraph_write_graph_edgelist(&uni, stdout);
    printf("\n");

    /* Destroy and free the graph list. */
    n = igraph_vector_ptr_size(&glist);
    for (i = 0; i < n; i++) {
        igraph_destroy(VECTOR(glist)[i]);
        free(VECTOR(glist)[i]);
    }
    igraph_vector_ptr_destroy(&glist);
    igraph_destroy(&uni);

    return 0;
}
