/*
   igraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
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
#include <stdlib.h>

int main(void) {

    igraph_t g;
    igraph_vector_int_list_t vecs, evecs;
    igraph_vector_int_t parents, inbound;
    igraph_int_t i;
    igraph_real_t weights[] = { 0, 2, 1, 0, 5, 2, 1, 1, 0, 2, 2, 8, 1, 1, 3, 1, 1, 4, 2, 1 };
    igraph_vector_t weights_vec = igraph_vector_view(weights, sizeof(weights) / sizeof(weights[0]));
    igraph_vs_t vs;

    /* Initialize the library. */
    igraph_setup();

    igraph_vector_int_list_init(&vecs, 0);
    igraph_vector_int_list_init(&evecs, 0);
    igraph_vector_int_init(&parents, 0);
    igraph_vector_int_init(&inbound, 0);

    igraph_vs_vector_small(&vs, 0, 1, 3, 5, 2, 1, -1);
    igraph_small(&g, 10, IGRAPH_DIRECTED,
                 0, 1, 0, 2, 0, 3,   1, 2, 1, 4, 1, 5,
                 2, 3, 2, 6,         3, 2, 3, 6,
                 4, 5, 4, 7,         5, 6, 5, 8, 5, 9,
                 7, 5, 7, 8,         8, 9,
                 5, 2,
                 2, 1,
                 -1);

    igraph_get_shortest_paths_dijkstra(&g, /*vertices=*/ &vecs,
                                       /*edges=*/ &evecs, /*from=*/ 0, /*to=*/ vs,
                                       &weights_vec, IGRAPH_OUT,
                                       &parents,
                                       /*inbound_edges=*/ &inbound);
    printf("Vertices:\n");
    for (i = 0; i < igraph_vector_int_list_size(&vecs); i++) {
        igraph_vector_int_print(igraph_vector_int_list_get_ptr(&vecs, i));
    }

    printf("\nEdges:\n");
    for (i = 0; i < igraph_vector_int_list_size(&evecs); i++) {
        igraph_vector_int_print(igraph_vector_int_list_get_ptr(&evecs, i));
    }

    printf("\nParents:\n");
    igraph_vector_int_print(&parents);

    printf("\nInbound:\n");
    igraph_vector_int_print(&inbound);

    igraph_vector_int_list_destroy(&vecs);
    igraph_vector_int_list_destroy(&evecs);
    igraph_vector_int_destroy(&parents);
    igraph_vector_int_destroy(&inbound);

    igraph_vs_destroy(&vs);
    igraph_destroy(&g);

    return 0;
}
