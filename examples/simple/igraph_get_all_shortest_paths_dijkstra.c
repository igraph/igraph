/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
/*
   IGraph library.
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

void print_and_destroy_items(igraph_vector_int_list_t* vec) {
    igraph_integer_t i;

    /* Sort the paths in a deterministic manner to avoid problems with
     * different qsort() implementations on different platforms */
    igraph_vector_int_list_sort(vec, igraph_vector_int_colex_cmp);
    for (i = 0; i < igraph_vector_int_list_size(vec); i++) {
        igraph_vector_int_print(igraph_vector_int_list_get_ptr(vec, i));
    }

    igraph_vector_int_list_destroy(vec);
}

int main(void) {

    igraph_t g;
    igraph_vector_int_list_t vertices, edges;

    igraph_real_t weights[] = { 0, 2, 1, 0, 5, 2, 1, 1, 0, 2, 2, 8, 1, 1, 3, 1, 1, 4, 2, 1 };

    igraph_vector_t weights_vec;
    igraph_vector_int_t nrgeo;
    igraph_vs_t vs;

    igraph_vector_int_list_init(&vertices, 0);
    igraph_vector_int_list_init(&edges, 0);
    igraph_vs_vector_small(&vs, 1, 3, 4, 5, 2, 1, -1);
    igraph_vector_int_init(&nrgeo, 0);
    igraph_small(&g, 10, IGRAPH_DIRECTED,
                 0, 1, 0, 2, 0, 3,   1, 2, 1, 4, 1, 5,
                 2, 3, 2, 6,         3, 2, 3, 6,
                 4, 5, 4, 7,         5, 6, 5, 8, 5, 9,
                 7, 5, 7, 8,         8, 9,
                 5, 2,
                 2, 1,
                 -1);

    igraph_vector_view(&weights_vec, weights, sizeof(weights) / sizeof(weights[0]));
    igraph_get_all_shortest_paths_dijkstra(
                &g,
                /*vertices=*/ &vertices, /*edges=*/ &edges, /*nrgeo=*/ &nrgeo,
                /*from=*/ 0, /*to=*/ vs,
                /*weights=*/ &weights_vec, /*mode=*/ IGRAPH_OUT);

    printf("Vertices:\n");
    print_and_destroy_items(&vertices);
    printf("\nEdges:\n");
    print_and_destroy_items(&edges);
    printf("\nNumber of geodesics:\n");
    igraph_vector_int_print(&nrgeo);

    igraph_vector_int_destroy(&nrgeo);
    igraph_vs_destroy(&vs);
    igraph_destroy(&g);

    return 0;
}
