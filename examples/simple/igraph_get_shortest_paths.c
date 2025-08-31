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

int check_evecs(const igraph_t *graph, const igraph_vector_int_list_t *vecs,
                const igraph_vector_int_list_t *evecs, int error_code) {

    igraph_bool_t directed = igraph_is_directed(graph);
    igraph_int_t i, n = igraph_vector_int_list_size(vecs);
    if (igraph_vector_int_list_size(evecs) != n) {
        exit(error_code + 1);
    }

    for (i = 0; i < n; i++) {
        igraph_vector_int_t *vvec = igraph_vector_int_list_get_ptr(vecs, i);
        igraph_vector_int_t *evec = igraph_vector_int_list_get_ptr(evecs, i);
        igraph_int_t j, n2 = igraph_vector_int_size(evec);
        if (igraph_vector_int_size(vvec) == 0 && n2 == 0) {
            continue;
        }
        if (igraph_vector_int_size(vvec) != n2 + 1) {
            exit(error_code + 2);
        }
        for (j = 0; j < n2; j++) {
            igraph_int_t edge = VECTOR(*evec)[j];
            igraph_int_t from = VECTOR(*vvec)[j];
            igraph_int_t to = VECTOR(*vvec)[j + 1];
            if (directed) {
                if (from != IGRAPH_FROM(graph, edge) ||
                    to   != IGRAPH_TO  (graph, edge)) {
                    exit(error_code);
                }
            } else {
                igraph_int_t from2 = IGRAPH_FROM(graph, edge);
                igraph_int_t to2 = IGRAPH_TO(graph, edge);
                igraph_int_t min1 = from < to ? from : to;
                igraph_int_t max1 = from < to ? to : from;
                igraph_int_t min2 = from2 < to2 ? from2 : to2;
                igraph_int_t max2 = from2 < to2 ? to2 : from2;
                if (min1 != min2 || max1 != max2) {
                    exit(error_code + 3);
                }
            }
        }
    }

    return 0;
}

int main(void) {

    igraph_t g;
    igraph_vector_int_list_t vecs, evecs;
    igraph_vector_int_t parents, inbound;
    igraph_int_t i;
    igraph_vs_t vs;

    /* Initialize the library. */
    igraph_setup();

    igraph_ring(&g, 10, IGRAPH_DIRECTED, 0, 1);

    igraph_vector_int_list_init(&vecs, 0);
    igraph_vector_int_list_init(&evecs, 0);
    igraph_vector_int_init(&parents, 0);
    igraph_vector_int_init(&inbound, 0);

    igraph_vs_vector_small(&vs, 1, 3, 5, 2, 1,  -1);

    igraph_get_shortest_paths(&g, NULL, &vecs, &evecs, 0, vs, IGRAPH_OUT, &parents, &inbound);

    check_evecs(&g, &vecs, &evecs, 10);

    for (i = 0; i < igraph_vector_int_list_size(&vecs); i++) {
        igraph_vector_int_print(igraph_vector_int_list_get_ptr(&vecs, i));
    }

    igraph_vector_int_print(&parents);
    igraph_vector_int_print(&inbound);

    igraph_vector_int_list_destroy(&vecs);
    igraph_vector_int_list_destroy(&evecs);
    igraph_vector_int_destroy(&parents);
    igraph_vector_int_destroy(&inbound);

    igraph_vs_destroy(&vs);
    igraph_destroy(&g);

    if (!IGRAPH_FINALLY_STACK_EMPTY) {
        return 1;
    }

    return 0;
}
