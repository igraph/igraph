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

int check_evecs(const igraph_t *graph, const igraph_vector_int_list_t *vecs,
                const igraph_vector_int_list_t *evecs, int error_code) {

    igraph_bool_t directed = igraph_is_directed(graph);
    igraph_integer_t i, n = igraph_vector_int_list_size(vecs);
    if (igraph_vector_int_list_size(evecs) != n) {
        exit(error_code + 1);
    }

    for (i = 0; i < n; i++) {
        igraph_vector_int_t *vvec = igraph_vector_int_list_get_ptr(vecs, i);
        igraph_vector_int_t *evec = igraph_vector_int_list_get_ptr(evecs, i);
        igraph_integer_t j, n2 = igraph_vector_int_size(evec);
        if (igraph_vector_int_size(vvec) == 0 && n2 == 0) {
            continue;
        }
        if (igraph_vector_int_size(vvec) != n2 + 1) {
            exit(error_code + 2);
        }
        for (j = 0; j < n2; j++) {
            igraph_integer_t edge = VECTOR(*evec)[j];
            igraph_integer_t from = VECTOR(*vvec)[j];
            igraph_integer_t to = VECTOR(*vvec)[j + 1];
            if (directed) {
                if (from != IGRAPH_FROM(graph, edge) ||
                    to   != IGRAPH_TO  (graph, edge)) {
                    exit(error_code);
                }
            } else {
                igraph_integer_t from2 = IGRAPH_FROM(graph, edge);
                igraph_integer_t to2 = IGRAPH_TO(graph, edge);
                igraph_integer_t min1 = from < to ? from : to;
                igraph_integer_t max1 = from < to ? to : from;
                igraph_integer_t min2 = from2 < to2 ? from2 : to2;
                igraph_integer_t max2 = from2 < to2 ? to2 : from2;
                if (min1 != min2 || max1 != max2) {
                    exit(error_code + 3);
                }
            }
        }
    }

    return 0;
}

int check_parents_inbound(const igraph_t* graph, const igraph_vector_int_t* parents,
                       const igraph_vector_int_t* inbound, int start, int error_code) {
    igraph_integer_t i, n = igraph_vcount(graph);

    if (igraph_vector_int_size(parents) != n ||
        igraph_vector_int_size(inbound) != n) {
        exit(error_code);
    }

    if (VECTOR(*parents)[start] != -1 || VECTOR(*inbound)[start] != -1) {
        printf("%" IGRAPH_PRId "\n", VECTOR(*parents)[start]);
        printf("%" IGRAPH_PRId "\n", VECTOR(*inbound)[start]);
        exit(error_code + 1);
    }

    for (i = 0; i < n; i++) {
        if (VECTOR(*parents)[i] == -2) {
            if (VECTOR(*inbound)[i] != -1) {
                exit(error_code + 2);
            }
        } else if (VECTOR(*parents)[i] == -1) {
            if (i != start) {
                exit(error_code + 3);
            }
            if (VECTOR(*inbound)[i] != -1) {
                exit(error_code + 4);
            }
        } else {
            igraph_integer_t eid = VECTOR(*inbound)[i];
            igraph_integer_t u = IGRAPH_FROM(graph, eid), v = IGRAPH_TO(graph, eid);
            if (v != i && !igraph_is_directed(graph)) {
                igraph_integer_t dummy = u;
                u = v;
                v = dummy;
            }
            if (v != i) {
                exit(error_code + 5);
            } else if (u != VECTOR(*parents)[i]) {
                exit(error_code + 6);
            }
        }
    }

    return 0;
}

int main() {

    igraph_t g;
    igraph_vector_int_list_t vecs, evecs;
    igraph_vector_int_t parents, inbound;
    igraph_integer_t i;
    igraph_real_t weights[] = { 1, 2, 3, 4, 5, 1, 1, 1, 1, 1 };
    igraph_real_t weights2[] = { 0, 2, 1, 0, 5, 2, 1, 1, 0, 2, 2, 8, 1, 1, 3, 1, 1, 4, 2, 1 };
    igraph_vector_t weights_vec;
    igraph_vs_t vs;

    /* Simple ring graph without weights */

    igraph_ring(&g, 10, IGRAPH_UNDIRECTED, 0, 1);

    igraph_vector_int_list_init(&vecs, 0);
    igraph_vector_int_list_init(&evecs, 0);
    igraph_vector_int_init(&parents, 0);
    igraph_vector_int_init(&inbound, 0);

    igraph_vs_vector_small(&vs, 0, 1, 3, 5, 2, 1,  -1);

    igraph_get_shortest_paths_dijkstra(&g, /*vertices=*/ &vecs,
                                       /*edges=*/ &evecs, /*from=*/ 0, /*to=*/ vs,
                                       /*weights=*/ 0, /*mode=*/ IGRAPH_OUT,
                                       &parents,
                                       /*inbound_edges=*/ &inbound);

    check_evecs(&g, &vecs, &evecs, 10);
    check_parents_inbound(&g, &parents, &inbound, /* from= */ 0, 40);

    for (i = 0; i < igraph_vector_int_list_size(&vecs); i++) {
        igraph_vector_int_print(igraph_vector_int_list_get_ptr(&vecs, i));
    }

    /* Same ring, but with weights */

    igraph_vector_view(&weights_vec, weights, sizeof(weights) / sizeof(weights[0]));
    igraph_get_shortest_paths_dijkstra(&g, /*vertices=*/ &vecs,
                                       /*edges=*/ &evecs, /*from=*/ 0, /*to=*/ vs,
                                       &weights_vec, IGRAPH_OUT,
                                       &parents,
                                       /*inbound_edges=*/ &inbound);

    check_evecs(&g, &vecs, &evecs, 20);
    check_parents_inbound(&g, &parents, &inbound, /* from= */ 0, 50);

    for (i = 0; i < igraph_vector_int_list_size(&vecs); i++) {
        igraph_vector_int_print(igraph_vector_int_list_get_ptr(&vecs, i));
    }

    igraph_destroy(&g);

    /* More complicated example */

    igraph_small(&g, 10, IGRAPH_DIRECTED,
                 0, 1, 0, 2, 0, 3,    1, 2, 1, 4, 1, 5,
                 2, 3, 2, 6,         3, 2, 3, 6,
                 4, 5, 4, 7,         5, 6, 5, 8, 5, 9,
                 7, 5, 7, 8,         8, 9,
                 5, 2,
                 2, 1,
                 -1);

    igraph_vector_view(&weights_vec, weights2, sizeof(weights2) / sizeof(weights2[0]));
    igraph_get_shortest_paths_dijkstra(&g, /*vertices=*/ &vecs,
                                       /*edges=*/ &evecs, /*from=*/ 0, /*to=*/ vs,
                                       &weights_vec, IGRAPH_OUT,
                                       &parents,
                                       /*inbound_edges=*/ &inbound);

    check_evecs(&g, &vecs, &evecs, 30);
    check_parents_inbound(&g, &parents, &inbound, /* from= */ 0, 60);

    for (i = 0; i < igraph_vector_int_list_size(&vecs); i++) {
        igraph_vector_int_print(igraph_vector_int_list_get_ptr(&vecs, i));
    }

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
