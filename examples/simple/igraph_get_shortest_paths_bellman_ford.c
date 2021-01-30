#include <igraph.h>

#include <stdlib.h>


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

int check_pred_inbound(const igraph_t* graph, const igraph_vector_long_t* pred,
                       const igraph_vector_long_t* inbound, int start, int error_code) {
    long int i, n = igraph_vcount(graph);

    if (igraph_vector_long_size(pred) != n ||
        igraph_vector_long_size(inbound) != n) {
        exit(error_code);
    }

    if (VECTOR(*pred)[start] != start || VECTOR(*inbound)[start] != -1) {
        exit(error_code + 1);
    }

    for (i = 0; i < n; i++) {
        if (VECTOR(*pred)[i] == -1) {
            if (VECTOR(*inbound)[i] != -1) {
                exit(error_code + 2);
            }
        } else if (VECTOR(*pred)[i] == i) {
            if (i != start) {
                exit(error_code + 3);
            }
            if (VECTOR(*inbound)[i] != -1) {
                exit(error_code + 4);
            }
        } else {
            long int eid = VECTOR(*inbound)[i];
            long int u = IGRAPH_FROM(graph, eid), v = IGRAPH_TO(graph, eid);
            if (v != i && !igraph_is_directed(graph)) {
                long int dummy = u;
                u = v;
                v = dummy;
            }
            if (v != i) {
                exit(error_code + 5);
            } else if (u != VECTOR(*pred)[i]) {
                exit(error_code + 6);
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
    igraph_real_t weights_data_0[] = { 0, 2, 1, 0, 5, 2, 1, 1, 0, 2, 2, 8, 1, 1, 3, 1, 1, 4, 2, 1 };
    igraph_real_t weights_data_1[] = { 6, 7, 8, -4, -2, -3, 9, 2, 7 };
    igraph_real_t weights_data_2[] = { 6, 7, 2, -4, -2, -3, 9, 2, 7 };
    igraph_vector_t weights_vec;
    igraph_vs_t vs;
    igraph_vector_ptr_init(&vecs, 6);
    igraph_vector_ptr_init(&evecs, 6);
    igraph_vector_long_init(&pred, 0);
    igraph_vector_long_init(&inbound, 0);

    for (i = 0; i < igraph_vector_ptr_size(&vecs); i++) {
        VECTOR(vecs)[i] = calloc(1, sizeof(igraph_vector_t));
        igraph_vector_init(VECTOR(vecs)[i], 0);VECTOR(evecs)[i] = calloc(1, sizeof(igraph_vector_t));
        igraph_vector_init(VECTOR(evecs)[i], 0);
    }
    igraph_vs_vector_small(&vs, 0, 1, 3, 5, 2, 1,  -1);
    igraph_small(&g, 10, IGRAPH_DIRECTED,
                0, 1, 0, 2, 0, 3,    1, 2, 1, 4, 1, 5,
                2, 3, 2, 6,         3, 2, 3, 6,
                4, 5, 4, 7,         5, 6, 5, 8, 5, 9,
                7, 5, 7, 8,         8, 9,
                5, 2,
                2, 1,
                -1);
    igraph_vector_view(&weights_vec, weights_data_0, sizeof(weights_data_0) / sizeof(igraph_real_t));
    igraph_get_shortest_paths_bellman_ford(&g, /*vertices=*/ &vecs,
                                    /*edges=*/ &evecs, /*from=*/ 0, /*to=*/ vs,
                                    &weights_vec, IGRAPH_OUT,
                                    /*predecessors=*/ &pred,
                                    /*inbound_edges=*/ &inbound);

    check_evecs(&g, &vecs, &evecs, 30);
    check_pred_inbound(&g, &pred, &inbound, /* from= */ 0, 60);

    for (i = 0; i < igraph_vector_ptr_size(&vecs); i++) {
            print_vector(VECTOR(vecs)[i]);
            igraph_vector_destroy(VECTOR(vecs)[i]);
            free(VECTOR(vecs)[i]);
            igraph_vector_destroy(VECTOR(evecs)[i]);
            free(VECTOR(evecs)[i]);
    }

    printf("\n");

    /***************************************/

    /* Graph with negative weights */

    igraph_vector_ptr_init(&vecs, 5);
    igraph_vector_ptr_init(&evecs, 5);
    igraph_vector_long_init(&pred, 0);
    igraph_vector_long_init(&inbound, 0);

    for (i = 0; i < igraph_vector_ptr_size(&vecs); i++) {
        VECTOR(vecs)[i] = calloc(1, sizeof(igraph_vector_t));
        igraph_vector_init(VECTOR(vecs)[i], 0);VECTOR(evecs)[i] = calloc(1, sizeof(igraph_vector_t));
        igraph_vector_init(VECTOR(evecs)[i], 0);
    }
    igraph_vs_vector_small(&vs, 0, 1, 3, 2, 1,  -1);
    igraph_small(&g, 5, IGRAPH_DIRECTED,
                 0, 1, 0, 3, 1, 3, 1, 4, 2, 1, 3, 2, 3, 4, 4, 0, 4, 2, -1);

    igraph_vector_view(&weights_vec, weights_data_1, sizeof(weights_data_1) / sizeof(igraph_real_t));
    igraph_get_shortest_paths_bellman_ford(&g, /*vertices=*/ &vecs,
                                    /*edges=*/ &evecs, /*from=*/ 0, /*to=*/ vs,
                                    &weights_vec, IGRAPH_OUT,
                                    /*predecessors=*/ &pred,
                                    /*inbound_edges=*/ &inbound);

    check_evecs(&g, &vecs, &evecs, 30);
    check_pred_inbound(&g, &pred, &inbound, /* from= */ 0, 60);

    for (i = 0; i < igraph_vector_ptr_size(&vecs); i++) {
            print_vector(VECTOR(vecs)[i]);
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
    if (igraph_get_shortest_paths_bellman_ford(&g, /*vertices=*/ &vecs,
                                    /*edges=*/ &evecs, /*from=*/ 0, /*to=*/ vs,
                                    &weights_vec, IGRAPH_OUT,
                                    /*predecessors=*/ &pred,
                                    /*inbound_edges=*/ &inbound) != IGRAPH_ENEGLOOP) {
        return 1;
    }

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