/* -*- mode: C -*-  */
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

void print_vector(igraph_vector_t *v) {
    long int i, l = igraph_vector_size(v);
    for (i = 0; i < l; i++) {
        printf(" %li", (long int) VECTOR(*v)[i]);
    }
    printf("\n");
}

int check_evecs(const igraph_t *graph, const igraph_vector_ptr_t *vecs,
                const igraph_vector_ptr_t *evecs, int error_code) {

    igraph_bool_t directed = igraph_is_directed(graph);
    long int i, n = igraph_vector_ptr_size(vecs);
    if (igraph_vector_ptr_size(evecs) != n) {
        exit(error_code + 1);
    }

    for (i = 0; i < n; i++) {
        igraph_vector_t *vvec = VECTOR(*vecs)[i];
        igraph_vector_t *evec = VECTOR(*evecs)[i];
        long int j, n2 = igraph_vector_size(evec);
        if (igraph_vector_size(vvec) == 0 && n2 == 0) {
            continue;
        }
        if (igraph_vector_size(vvec) != n2 + 1) {
            exit(error_code + 2);
        }
        for (j = 0; j < n2; j++) {
            long int edge = VECTOR(*evec)[j];
            long int from = VECTOR(*vvec)[j];
            long int to = VECTOR(*vvec)[j + 1];
            if (directed) {
                if (from != IGRAPH_FROM(graph, edge) ||
                    to   != IGRAPH_TO  (graph, edge)) {
                    exit(error_code);
                }
            } else {
                long int from2 = IGRAPH_FROM(graph, edge);
                long int to2 = IGRAPH_TO(graph, edge);
                long int min1 = from < to ? from : to;
                long int max1 = from < to ? to : from;
                long int min2 = from2 < to2 ? from2 : to2;
                long int max2 = from2 < to2 ? to2 : from2;
                if (min1 != min2 || max1 != max2) {
                    exit(error_code + 3);
                }
            }
        }
    }

    return 0;
}

int main() {

    igraph_t g;
    igraph_vector_ptr_t vecs, evecs;
    igraph_vector_long_t pred, inbound;
    long int i;
    igraph_vs_t vs;

    igraph_ring(&g, 10, IGRAPH_DIRECTED, 0, 1);

    igraph_vector_ptr_init(&vecs, 5);
    igraph_vector_ptr_init(&evecs, 5);
    igraph_vector_long_init(&pred, 0);
    igraph_vector_long_init(&inbound, 0);

    for (i = 0; i < igraph_vector_ptr_size(&vecs); i++) {
        VECTOR(vecs)[i] = calloc(1, sizeof(igraph_vector_t));
        igraph_vector_init(VECTOR(vecs)[i], 0);
        VECTOR(evecs)[i] = calloc(1, sizeof(igraph_vector_t));
        igraph_vector_init(VECTOR(evecs)[i], 0);
    }
    igraph_vs_vector_small(&vs, 1, 3, 5, 2, 1,  -1);

    igraph_get_shortest_paths(&g, &vecs, &evecs, 0, vs, IGRAPH_OUT, &pred, &inbound);

    check_evecs(&g, &vecs, &evecs, 10);

    for (i = 0; i < igraph_vector_ptr_size(&vecs); i++) {
        print_vector(VECTOR(vecs)[i]);
        igraph_vector_destroy(VECTOR(vecs)[i]);
        free(VECTOR(vecs)[i]);
        igraph_vector_destroy(VECTOR(evecs)[i]);
        free(VECTOR(evecs)[i]);
    }

    igraph_vector_long_print(&pred);
    igraph_vector_long_print(&inbound);

    igraph_vector_ptr_destroy(&vecs);
    igraph_vector_ptr_destroy(&evecs);
    igraph_vector_long_destroy(&pred);
    igraph_vector_long_destroy(&inbound);
    igraph_vs_destroy(&vs);
    igraph_destroy(&g);

    if (!IGRAPH_FINALLY_STACK_EMPTY) {
        return 1;
    }

    return 0;
}
