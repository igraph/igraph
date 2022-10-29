/* -*- mode: C -*-  */
/*
   IGraph library.
   (C) 2006-2021 The igraph development team  Gabor Csardi <csardi.gabor@gmail.com>

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
#include <stdlib.h>

void check_evecs(const igraph_t *graph, const igraph_vector_int_list_t *vecs,
                 const igraph_vector_int_list_t *evecs) {

    igraph_bool_t directed = igraph_is_directed(graph);
    igraph_integer_t i, n = igraph_vector_int_list_size(vecs);

    IGRAPH_ASSERT(igraph_vector_int_list_size(evecs) == n);

    for (i = 0; i < n; i++) {
        igraph_vector_int_t *vvec = igraph_vector_int_list_get_ptr(vecs, i);
        igraph_vector_int_t *evec = igraph_vector_int_list_get_ptr(evecs, i);
        igraph_integer_t j, n2 = igraph_vector_int_size(evec);
        if (igraph_vector_int_size(vvec) == 0 && n2 == 0) {
            continue;
        }
        IGRAPH_ASSERT(igraph_vector_int_size(vvec) == n2 + 1);

        for (j = 0; j < n2; j++) {
            igraph_integer_t edge = VECTOR(*evec)[j];
            igraph_integer_t from = VECTOR(*vvec)[j];
            igraph_integer_t to = VECTOR(*vvec)[j + 1];
            if (directed) {
                IGRAPH_ASSERT(from == IGRAPH_FROM(graph, edge) &&
                              to == IGRAPH_TO(graph, edge));
            } else {
                igraph_integer_t from2 = IGRAPH_FROM(graph, edge);
                igraph_integer_t to2 = IGRAPH_TO(graph, edge);
                igraph_integer_t min1 = from < to ? from : to;
                igraph_integer_t max1 = from < to ? to : from;
                igraph_integer_t min2 = from2 < to2 ? from2 : to2;
                igraph_integer_t max2 = from2 < to2 ? to2 : from2;
                IGRAPH_ASSERT(min1 == min2 && max1 == max2);
            }
        }
    }
}

void check_parents_inbound(const igraph_t* graph, const igraph_vector_int_t* parents,
                        const igraph_vector_int_t* inbound, int start) {

    igraph_integer_t i, n = igraph_vcount(graph);

    IGRAPH_ASSERT(igraph_vector_int_size(parents) == n);
    IGRAPH_ASSERT(igraph_vector_int_size(inbound) == n);

    IGRAPH_ASSERT(VECTOR(*parents)[start] == -1 && VECTOR(*inbound)[start] == -1);

    for (i = 0; i < n; i++) {
        if (VECTOR(*parents)[i] == -2) {
            IGRAPH_ASSERT(VECTOR(*inbound)[i] == -1);

        } else if (VECTOR(*parents)[i] == -1) {
            IGRAPH_ASSERT(i == start);

            IGRAPH_ASSERT(VECTOR(*inbound)[i] == -1);

        } else {
            igraph_integer_t eid = VECTOR(*inbound)[i];
            igraph_integer_t u = IGRAPH_FROM(graph, eid), v = IGRAPH_TO(graph, eid);
            if (v != i && !igraph_is_directed(graph)) {
                igraph_integer_t dummy = u;
                u = v;
                v = dummy;
            }
            IGRAPH_ASSERT(v == i);
            IGRAPH_ASSERT(u == VECTOR(*parents)[i]);
        }
    }
}

int main(void) {
    igraph_t g;
    igraph_vector_int_list_t vecs, evecs;
    igraph_vector_int_t parents, inbound;
    igraph_real_t weights_data_0[] = { 0, 2, 1, 0, 5, 2, 1, 1, 0, 2, 2, 8, 1, 1, 3, 1, 1, 4, 2, 1 };
    igraph_real_t weights_data_1[] = { 6, 7, 8, -4, -2, -3, 9, 2, 7 };
    igraph_real_t weights_data_2[] = { 6, 7, 2, -4, -2, -3, 9, 2, 7 };
    igraph_vector_t weights_vec;
    igraph_vs_t vs;
    igraph_integer_t vs_size;

    igraph_small(&g, 10, IGRAPH_DIRECTED,
                0, 1, 0, 2, 0, 3,    1, 2, 1, 4, 1, 5,
                2, 3, 2, 6,         3, 2, 3, 6,
                4, 5, 4, 7,         5, 6, 5, 8, 5, 9,
                7, 5, 7, 8,         8, 9,
                5, 2,
                2, 1,
                -1);

    igraph_vector_int_init(&parents, 0);
    igraph_vector_int_init(&inbound, 0);

    printf("Paths to only some vertices\n");

    igraph_vs_vector_small(&vs, 0, 1, 3, 5, 2, 1,  -1);
    igraph_vs_size(&g, &vs, &vs_size);

    igraph_vector_int_list_init(&vecs, 0);
    igraph_vector_int_list_init(&evecs, 0);

    igraph_vector_view(&weights_vec, weights_data_0, sizeof(weights_data_0) / sizeof(weights_data_0[0]));
    igraph_get_shortest_paths_bellman_ford(&g, /*vertices=*/ &vecs, /*edges=*/ &evecs,
                                           /*from=*/ 0, /*to=*/ vs,
                                           &weights_vec, IGRAPH_OUT,
                                           &parents,
                                           /*inbound_edges=*/ &inbound);

    check_evecs(&g, &vecs, &evecs);
    check_parents_inbound(&g, &parents, &inbound, /* from= */ 0);

    print_vector_int_list(&vecs);

    printf("\nPaths to all vertices\n");

    vs_size = igraph_vcount(&g);

    igraph_get_shortest_paths_bellman_ford(&g, /*vertices=*/ &vecs, /*edges=*/ &evecs,
                                           /*from=*/ 0, /*to=*/ igraph_vss_all(),
                                           &weights_vec, IGRAPH_OUT,
                                           &parents,
                                           /*inbound_edges=*/ &inbound);

    check_evecs(&g, &vecs, &evecs);
    check_parents_inbound(&g, &parents, &inbound, /* from= */ 0);

    print_vector_int_list(&vecs);

    igraph_vector_int_list_destroy(&vecs);
    igraph_vector_int_list_destroy(&evecs);

    igraph_vector_int_destroy(&parents);
    igraph_vector_int_destroy(&inbound);

    igraph_vs_destroy(&vs);
    igraph_destroy(&g);


    printf("\nGraph with negative weights\n");

    /***************************************/

    /* Graph with negative weights */

    igraph_vector_int_list_init(&vecs, 0);
    igraph_vector_int_list_init(&evecs, 0);
    igraph_vector_int_init(&parents, 0);
    igraph_vector_int_init(&inbound, 0);

    igraph_vs_vector_small(&vs, 0, 1, 3, 2, 1,  -1);
    igraph_small(&g, 5, IGRAPH_DIRECTED,
                 0, 1, 0, 3, 1, 3, 1, 4, 2, 1, 3, 2, 3, 4, 4, 0, 4, 2,
                 -1);

    igraph_vector_view(&weights_vec, weights_data_1, sizeof(weights_data_1) / sizeof(weights_data_1[0]));
    igraph_get_shortest_paths_bellman_ford(&g, /*vertices=*/ &vecs, /*edges=*/ &evecs,
                                           /*from=*/ 0, /*to=*/ vs,
                                           &weights_vec, IGRAPH_OUT,
                                           &parents,
                                           /*inbound_edges=*/ &inbound);

    check_evecs(&g, &vecs, &evecs);
    check_parents_inbound(&g, &parents, &inbound, /* from= */ 0);

    print_vector_int_list(&vecs);

    /***************************************/

    /* Same graph with negative loop */
    igraph_set_error_handler(igraph_error_handler_ignore);
    igraph_vector_view(&weights_vec, weights_data_2,
                       sizeof(weights_data_2) / sizeof(weights_data_2[0]));
    IGRAPH_ASSERT(igraph_get_shortest_paths_bellman_ford(&g, /*vertices=*/ &vecs, /*edges=*/ &evecs,
                                                         /*from=*/ 0, /*to=*/ vs,
                                                         &weights_vec, IGRAPH_OUT,
                                                         &parents,
                                                         /*inbound_edges=*/ &inbound) == IGRAPH_ENEGLOOP);

    igraph_vector_int_list_destroy(&vecs);
    igraph_vector_int_list_destroy(&evecs);
    igraph_vector_int_destroy(&parents);
    igraph_vector_int_destroy(&inbound);

    igraph_vs_destroy(&vs);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();
    return 0;
}
