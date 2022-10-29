/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#include "test_utilities.h"

int main(void) {
    const igraph_integer_t edges[] = { 0, 1, 0, 2, 1, 6, 2, 6, 1, 3, 1, 4, 1, 5,
                                    3, 2, 4, 2, 5, 2
                                  };
    igraph_t g;
    igraph_vector_int_t edgev;
    igraph_vector_int_list_t resvertices, resedges;
    igraph_vector_int_t parents, inbound_edges;
    igraph_integer_t vcount, i;

    igraph_vector_int_view(&edgev, edges, sizeof(edges) / sizeof(edges[0]));
    vcount = igraph_vector_int_max(&edgev) + 1;
    igraph_create(&g, &edgev, vcount, IGRAPH_DIRECTED);

    igraph_vector_int_list_init(&resvertices, 0);
    igraph_vector_int_list_init(&resedges, 0);
    igraph_vector_int_init(&parents, 0);
    igraph_vector_int_init(&inbound_edges, 0);

    igraph_get_shortest_paths(&g, &resvertices, &resedges, /*from=*/ 0,
                              /*to=*/ igraph_vss_all(), /*mode=*/ IGRAPH_OUT,
                              &parents, &inbound_edges);

    for (i = 0; i < vcount; i++) {
        igraph_vector_int_t *v1 = igraph_vector_int_list_get_ptr(&resvertices, i);
        igraph_vector_int_t *v2 = igraph_vector_int_list_get_ptr(&resedges, i);
        printf("%" IGRAPH_PRId " V: ", i);
        igraph_vector_int_print(v1);
        printf("%" IGRAPH_PRId " E: ", i);
        igraph_vector_int_print(v2);
    }
    printf("pred: ");
    igraph_vector_int_print(&parents);
    printf("inbe: ");
    igraph_vector_int_print(&inbound_edges);

    igraph_vector_int_destroy(&inbound_edges);
    igraph_vector_int_destroy(&parents);
    igraph_vector_int_list_destroy(&resedges);
    igraph_vector_int_list_destroy(&resvertices);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();

    return 0;
}
