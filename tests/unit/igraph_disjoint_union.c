/*
   igraph library.
   Copyright (C) 2006-2024  The igraph development team <igraph@igraph.org>

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

#include "test_utilities.h"

int main(void) {
    igraph_t left, right, uni;
    igraph_vector_ptr_t glist;
    igraph_int_t i, n;

    igraph_small(&left, 4, IGRAPH_UNDIRECTED, 0,1, 1,2, 2,2, 2,3, -1);
    igraph_small(&right, 5, IGRAPH_UNDIRECTED, 0,1, 1,2, 2,2, 2,4, -1);

    igraph_disjoint_union(&uni, &left, &right);
    print_graph(&uni);
    printf("\n");

    igraph_destroy(&left);
    igraph_destroy(&right);
    igraph_destroy(&uni);

    /* Empty graph list; the result is the directed null graph. */
    igraph_vector_ptr_init(&glist, 0);
    igraph_disjoint_union_many(&uni, &glist);
    print_graph(&uni);
    igraph_vector_ptr_destroy(&glist);
    igraph_destroy(&uni);

    /* Non-empty graph list. */
    igraph_vector_ptr_init(&glist, 10);
    n = igraph_vector_ptr_size(&glist);
    for (i = 0; i < n; i++) {
        VECTOR(glist)[i] = IGRAPH_CALLOC(1, igraph_t);
        igraph_small(VECTOR(glist)[i], 2, IGRAPH_DIRECTED, 0,1, 1,0, -1);
    }
    igraph_disjoint_union_many(&uni, &glist);
    print_graph(&uni);
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
