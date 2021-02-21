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
#include "test_utilities.inc"
#include <stdlib.h>

void check_evecs(const igraph_t *graph, const igraph_vector_ptr_t *vecs,
                 const igraph_vector_ptr_t *evecs) {

    igraph_bool_t directed = igraph_is_directed(graph);
    long int i, n = igraph_vector_ptr_size(vecs);

    IGRAPH_ASSERT(igraph_vector_ptr_size(evecs) == n);

    for (i = 0; i < n; i++) {
        igraph_vector_t *vvec = VECTOR(*vecs)[i];
        igraph_vector_t *evec = VECTOR(*evecs)[i];
        long int j, n2 = igraph_vector_size(evec);
        if (igraph_vector_size(vvec) == 0 && n2 == 0) {
            continue;
        }
        IGRAPH_ASSERT(igraph_vector_size(vvec) == n2 + 1);

        for (j = 0; j < n2; j++) {
            long int edge = VECTOR(*evec)[j];
            long int from = VECTOR(*vvec)[j];
            long int to = VECTOR(*vvec)[j + 1];
            if (directed) {
                IGRAPH_ASSERT(from == IGRAPH_FROM(graph, edge) &&
                              to == IGRAPH_TO(graph, edge));
            } else {
                long int from2 = IGRAPH_FROM(graph, edge);
                long int to2 = IGRAPH_TO(graph, edge);
                long int min1 = from < to ? from : to;
                long int max1 = from < to ? to : from;
                long int min2 = from2 < to2 ? from2 : to2;
                long int max2 = from2 < to2 ? to2 : from2;
                IGRAPH_ASSERT(min1 == min2 && max1 == max2);
            }
        }
    }
}

void check_pred_inbound(const igraph_t* graph, const igraph_vector_long_t* pred,
                        const igraph_vector_long_t* inbound, int start) {

    long int i, n = igraph_vcount(graph);

    IGRAPH_ASSERT(igraph_vector_long_size(pred) == n);
    IGRAPH_ASSERT(igraph_vector_long_size(inbound) == n);

    IGRAPH_ASSERT(VECTOR(*pred)[start] == start && VECTOR(*inbound)[start] == -1);

    for (i = 0; i < n; i++) {
        if (VECTOR(*pred)[i] == -1) {
            IGRAPH_ASSERT(VECTOR(*inbound)[i] == -1);

        } else if (VECTOR(*pred)[i] == i) {
            IGRAPH_ASSERT(i == start);

            IGRAPH_ASSERT(VECTOR(*inbound)[i] == -1);

        } else {
            long int eid = VECTOR(*inbound)[i];
            long int u = IGRAPH_FROM(graph, eid), v = IGRAPH_TO(graph, eid);
            if (v != i && !igraph_is_directed(graph)) {
                long int dummy = u;
                u = v;
                v = dummy;
            }
            IGRAPH_ASSERT(v == i);
            IGRAPH_ASSERT(u == VECTOR(*pred)[i]);
        }
    }
}

int main() {
    igraph_t g;
    igraph_vector_ptr_t vecs, evecs;
    igraph_vector_long_t pred, inbound;
    long int i;
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

    igraph_vector_long_init(&pred, 0);
    igraph_vector_long_init(&inbound, 0);

    printf("Paths to only some vertices\n");

    igraph_vs_vector_small(&vs, 0, 1, 3, 5, 2, 1,  -1);
    igraph_vs_size(&g, &vs, &vs_size);

    igraph_vector_ptr_init(&vecs, vs_size);
    igraph_vector_ptr_init(&evecs, vs_size);

    for (i = 0; i < igraph_vector_ptr_size(&vecs); i++) {
        VECTOR(vecs)[i] = calloc(1, sizeof(igraph_vector_t));
        igraph_vector_init(VECTOR(vecs)[i], 0);
        VECTOR(evecs)[i] = calloc(1, sizeof(igraph_vector_t));
        igraph_vector_init(VECTOR(evecs)[i], 0);
    }

    igraph_vector_view(&weights_vec, weights_data_0, sizeof(weights_data_0) / sizeof(igraph_real_t));
    igraph_get_shortest_paths_bellman_ford(&g, /*vertices=*/ &vecs, /*edges=*/ &evecs,
                                           /*from=*/ 0, /*to=*/ vs,
                                           &weights_vec, IGRAPH_OUT,
                                           /*predecessors=*/ &pred,
                                           /*inbound_edges=*/ &inbound);

    check_evecs(&g, &vecs, &evecs);
    check_pred_inbound(&g, &pred, &inbound, /* from= */ 0);

    for (i = 0; i < igraph_vector_ptr_size(&vecs); i++) {
        print_vector_round(VECTOR(vecs)[i]);
        igraph_vector_destroy(VECTOR(vecs)[i]);
        free(VECTOR(vecs)[i]);
        igraph_vector_destroy(VECTOR(evecs)[i]);
        free(VECTOR(evecs)[i]);
    }

    printf("\nPaths to all vertices\n");

    vs_size = igraph_vcount(&g);

    igraph_vector_ptr_resize(&vecs, vs_size);
    igraph_vector_ptr_resize(&evecs, vs_size);

    for (i = 0; i < igraph_vector_ptr_size(&vecs); i++) {
        VECTOR(vecs)[i] = calloc(1, sizeof(igraph_vector_t));
        igraph_vector_init(VECTOR(vecs)[i], 0);
        VECTOR(evecs)[i] = calloc(1, sizeof(igraph_vector_t));
        igraph_vector_init(VECTOR(evecs)[i], 0);
    }

    igraph_get_shortest_paths_bellman_ford(&g, /*vertices=*/ &vecs, /*edges=*/ &evecs,
                                           /*from=*/ 0, /*to=*/ igraph_vss_all(),
                                           &weights_vec, IGRAPH_OUT,
                                           /*predecessors=*/ &pred,
                                           /*inbound_edges=*/ &inbound);

    check_evecs(&g, &vecs, &evecs);
    check_pred_inbound(&g, &pred, &inbound, /* from= */ 0);

    for (i = 0; i < igraph_vector_ptr_size(&vecs); i++) {
        print_vector_round(VECTOR(vecs)[i]);
        igraph_vector_destroy(VECTOR(vecs)[i]);
        free(VECTOR(vecs)[i]);
        igraph_vector_destroy(VECTOR(evecs)[i]);
        free(VECTOR(evecs)[i]);
    }

    igraph_vector_ptr_destroy(&vecs);
    igraph_vector_ptr_destroy(&evecs);

    igraph_vector_long_destroy(&pred);
    igraph_vector_long_destroy(&inbound);

    igraph_vs_destroy(&vs);
    igraph_destroy(&g);


    printf("\nGraph with negative weights\n");

    /***************************************/

    /* Graph with negative weights */

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
    igraph_vs_vector_small(&vs, 0, 1, 3, 2, 1,  -1);
    igraph_small(&g, 5, IGRAPH_DIRECTED,
                 0, 1, 0, 3, 1, 3, 1, 4, 2, 1, 3, 2, 3, 4, 4, 0, 4, 2,
                 -1);

    igraph_vector_view(&weights_vec, weights_data_1, sizeof(weights_data_1) / sizeof(igraph_real_t));
    igraph_get_shortest_paths_bellman_ford(&g, /*vertices=*/ &vecs, /*edges=*/ &evecs,
                                           /*from=*/ 0, /*to=*/ vs,
                                           &weights_vec, IGRAPH_OUT,
                                           /*predecessors=*/ &pred,
                                           /*inbound_edges=*/ &inbound);

    check_evecs(&g, &vecs, &evecs);
    check_pred_inbound(&g, &pred, &inbound, /* from= */ 0);

    for (i = 0; i < igraph_vector_ptr_size(&vecs); i++) {
        print_vector_round(VECTOR(vecs)[i]);
        igraph_vector_destroy(VECTOR(vecs)[i]);
        free(VECTOR(vecs)[i]);
        igraph_vector_destroy(VECTOR(evecs)[i]);
        free(VECTOR(evecs)[i]);
    }

    /***************************************/

    /* Same graph with negative loop */
    igraph_set_error_handler(igraph_error_handler_ignore);
    igraph_vector_view(&weights_vec, weights_data_2,
                       sizeof(weights_data_2) / sizeof(igraph_real_t));
    IGRAPH_ASSERT(igraph_get_shortest_paths_bellman_ford(&g, /*vertices=*/ &vecs, /*edges=*/ &evecs,
                                                         /*from=*/ 0, /*to=*/ vs,
                                                         &weights_vec, IGRAPH_OUT,
                                                         /*predecessors=*/ &pred,
                                                         /*inbound_edges=*/ &inbound) == IGRAPH_ENEGLOOP);

    igraph_vector_ptr_destroy(&vecs);
    igraph_vector_ptr_destroy(&evecs);
    igraph_vector_long_destroy(&pred);
    igraph_vector_long_destroy(&inbound);

    igraph_vs_destroy(&vs);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();
    return 0;
}
