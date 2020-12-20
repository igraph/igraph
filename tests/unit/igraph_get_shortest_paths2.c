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

#include "test_utilities.inc"

int main() {
    const igraph_real_t edges[] = { 0, 1, 0, 2, 1, 6, 2, 6, 1, 3, 1, 4, 1, 5,
                                    3, 2, 4, 2, 5, 2
                                  };
    igraph_t g;
    igraph_vector_t edgev;
    igraph_vector_ptr_t resvertices, resedges;
    igraph_vector_long_t predecessors, inbound_edges;
    int vcount, i;

    igraph_vector_view(&edgev, edges, sizeof(edges) / sizeof(igraph_real_t));
    vcount = igraph_vector_max(&edgev) + 1;
    igraph_create(&g, &edgev, vcount, IGRAPH_DIRECTED);

    igraph_vector_ptr_init(&resvertices, vcount);
    igraph_vector_ptr_init(&resedges, vcount);
    igraph_vector_long_init(&predecessors, 0);
    igraph_vector_long_init(&inbound_edges, 0);

    for (i = 0; i < vcount; i++) {
        igraph_vector_t *v1 = malloc(sizeof(igraph_vector_t));
        igraph_vector_t *v2 = malloc(sizeof(igraph_vector_t));
        if (!v1 || !v2) {
            exit(2);
        }
        igraph_vector_init(v1, 0);
        igraph_vector_init(v2, 0);
        VECTOR(resvertices)[i] = v1;
        VECTOR(resedges)[i] = v2;
    }

    igraph_get_shortest_paths(&g, &resvertices, &resedges, /*from=*/ 0,
                              /*to=*/ igraph_vss_all(), /*mode=*/ IGRAPH_OUT,
                              &predecessors, &inbound_edges);

    for (i = 0; i < vcount; i++) {
        igraph_vector_t *v1 = VECTOR(resvertices)[i];
        igraph_vector_t *v2 = VECTOR(resedges)[i];
        printf("%i V: ", i);
        igraph_vector_print(v1);
        printf("%i E: ", i);
        igraph_vector_print(v2);
    }
    printf("pred: ");
    igraph_vector_long_print(&predecessors);
    printf("inbe: ");
    igraph_vector_long_print(&inbound_edges);

    igraph_vector_long_destroy(&inbound_edges);
    igraph_vector_long_destroy(&predecessors);
    for (i = 0; i < vcount; i++) {
        igraph_vector_t *v1 = VECTOR(resvertices)[i];
        igraph_vector_t *v2 = VECTOR(resedges)[i];
        igraph_vector_destroy(v1);
        igraph_vector_destroy(v2);
        igraph_free(v1);
        igraph_free(v2);
    }
    igraph_vector_ptr_destroy(&resedges);
    igraph_vector_ptr_destroy(&resvertices);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();

    return 0;
}
