/* IGraph library.
   Copyright (C) 2022  The igraph development team <igraph@igraph.org>

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include <igraph.h>
#include "test_utilities.h"

#define LENGTH 10

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


igraph_error_t no_heuristic(igraph_real_t *result, igraph_integer_t vertex_id, void *extra) {
    (void) vertex_id;
    (void) extra;
    *result = 0;
    return IGRAPH_SUCCESS;
}

igraph_error_t lattice_heuristic(igraph_real_t *result, igraph_integer_t vertex_id, void *extra) {
    int x[4];
    for (int i = 0; i < 4; i++) {
        x[i] = vertex_id % LENGTH;
        vertex_id /= LENGTH;
    }
    *result = abs(LENGTH - 1 - x[0]) + x[1] + x[2] + x[3];
    return IGRAPH_SUCCESS;
}

struct xyt {
    igraph_vector_t x;
    igraph_vector_t y;
    igraph_integer_t target;
};

igraph_error_t euclidean_heuristic(igraph_real_t *result, igraph_integer_t vertex_id, void *extra) {
    struct xyt *xyp = extra;
    igraph_real_t xt, xf, yt, yf;
    xt = VECTOR(xyp->x)[xyp->target];
    yt = VECTOR(xyp->y)[xyp->target];
    xf = VECTOR(xyp->x)[vertex_id];
    yf = VECTOR(xyp->y)[vertex_id];
    *result = sqrt((xt-xf)*(xt-xf) + (yt-yf)*(yt-yf));
    return IGRAPH_SUCCESS;
}
int main(void) {

    igraph_t g;
    igraph_vector_int_list_t vecs, evecs;
    igraph_vector_int_t parents, inbound;
    igraph_integer_t i;
    igraph_real_t weights[] = { 1, 2, 3, 4, 5, 1, 1, 1, 1, 1 };
    igraph_vector_t weights_vec;
    igraph_vs_t vs;
    struct xyt xy;
    igraph_vector_int_t dimvector;
    igraph_integer_t dims[] = {LENGTH, LENGTH, LENGTH, LENGTH};

    //set seed for grg random graph generation
    igraph_rng_seed(igraph_rng_default(), 42);

    /* Simple ring graph without weights */

    igraph_ring(&g, 10, IGRAPH_UNDIRECTED, 0, 1);

    igraph_vector_int_list_init(&vecs, 0);
    igraph_vector_int_list_init(&evecs, 0);
    igraph_vector_int_init(&parents, 0);
    igraph_vector_int_init(&inbound, 0);

    igraph_vs_vector_small(&vs, 0, 1, 3, 5, 2, 1,  -1);

    printf("Astar, unweighted, no heuristic:\n");
    igraph_get_shortest_paths_astar(&g, /*vertices=*/ &vecs,
                                       /*edges=*/ &evecs, /*from=*/ 0, /*to=*/ vs,
                                       /*weights=*/ 0, /*mode=*/ IGRAPH_OUT,
                                       &parents,
                                       /*inbound_edges=*/ &inbound, no_heuristic, NULL);

    check_evecs(&g, &vecs, &evecs, /*error code*/10);
    check_parents_inbound(&g, &parents, &inbound, /* from= */ 0, /*error code*/40);

    for (i = 0; i < igraph_vector_int_list_size(&vecs); i++) {
        igraph_vector_int_print(igraph_vector_int_list_get_ptr(&vecs, i));
    }

    /* Same ring, but with weights */

    printf("Astar, weighted, no heuristic:\n");
    igraph_vector_view(&weights_vec, weights, sizeof(weights) / sizeof(weights[0]));
    igraph_get_shortest_paths_astar(&g, /*vertices=*/ &vecs,
                                       /*edges=*/ &evecs, /*from=*/ 0, /*to=*/ vs,
                                       &weights_vec, IGRAPH_OUT,
                                       &parents,
                                       /*inbound_edges=*/ &inbound, no_heuristic, NULL);

    check_evecs(&g, &vecs, &evecs, 20);
    check_parents_inbound(&g, &parents, &inbound, /* from= */ 0, /*error code*/50);

    for (i = 0; i < igraph_vector_int_list_size(&vecs); i++) {
        igraph_vector_int_print(igraph_vector_int_list_get_ptr(&vecs, i));
    }

    igraph_destroy(&g);
    igraph_vs_destroy(&vs);

    printf("Astar, unweighted, lattice with manhattan distance heuristic:\n");
    igraph_vector_int_view(&dimvector, dims, sizeof(dims)/sizeof(dims[0]));
    igraph_vs_vector_small(&vs, LENGTH-1, -1);

    igraph_square_lattice(&g, &dimvector, /*nei*/ 1, IGRAPH_UNDIRECTED, /*mutual*/ false, /*periodic*/NULL);

    //print_graph_canon(&g);
    igraph_get_shortest_paths_astar(&g, /*vertices=*/ &vecs,
                                       /*edges=*/ &evecs, /*from=*/ 0, /*to=*/ vs,
                                       /*weights*/NULL, IGRAPH_OUT,
                                       &parents,
                                       /*inbound_edges=*/ &inbound, lattice_heuristic, NULL);
    igraph_vector_int_print(&(VECTOR(vecs)[0]));
    check_evecs(&g, &vecs, &evecs, /*error code*/60);
    check_parents_inbound(&g, &parents, &inbound, /* from= */ 0, /*error code*/70);

    igraph_destroy(&g);
 
    printf("Astar, unweighted, grg with euclidean distance heuristic:\n");
    xy.target = LENGTH-1; //just because that was the end vertex last test
    igraph_vector_init(&xy.x, 0);
    igraph_vector_init(&xy.y, 0);

    igraph_grg_game(&g, /*nodes*/100, /*radius*/0.2, /*torus*/ false, &xy.x, &xy.y); 
    igraph_vector_init(&weights_vec, igraph_ecount(&g));

    for (int i = 0; i < igraph_ecount(&g); i++) {
        igraph_real_t xt, xf, yt, yf;
        xt = VECTOR(xy.x)[IGRAPH_TO(&g, i)];
        xf = VECTOR(xy.x)[IGRAPH_FROM(&g, i)];
        yt = VECTOR(xy.y)[IGRAPH_TO(&g, i)];
        yf = VECTOR(xy.y)[IGRAPH_FROM(&g, i)];
        VECTOR(weights_vec)[i] = sqrt((xt-xf)*(xt-xf) + (yt-yf)*(yt-yf));
    }

    igraph_get_shortest_paths_astar(&g, /*vertices=*/ &vecs,
                                       /*edges=*/ &evecs, /*from=*/ 0, /*to=*/ vs,
                                       /*weights*/&weights_vec, IGRAPH_OUT,
                                       &parents,
                                       /*inbound_edges=*/ &inbound, euclidean_heuristic, &xy);

    igraph_vector_int_print(&(VECTOR(vecs)[0]));
    check_evecs(&g, &vecs, &evecs, /*error code*/80);
    check_parents_inbound(&g, &parents, &inbound, /* from= */ 0, /*error code*/90);

    printf("Check with dijkstra:\n");
    igraph_get_shortest_paths_dijkstra(&g, /*vertices=*/ &vecs,
                                       /*edges=*/ &evecs, /*from=*/ 0, /*to=*/ vs,
                                       /*weights*/&weights_vec, IGRAPH_OUT,
                                       &parents,
                                       /*inbound_edges=*/ &inbound);

    igraph_vector_int_print(&(VECTOR(vecs)[0]));

    igraph_vector_int_list_destroy(&vecs);
    igraph_vector_int_list_destroy(&evecs);
    igraph_vector_int_destroy(&parents);
    igraph_vector_int_destroy(&inbound);
    igraph_vector_destroy(&weights_vec);
    igraph_vector_destroy(&xy.x);
    igraph_vector_destroy(&xy.y);

    igraph_vs_destroy(&vs);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();
    return 0;
}
