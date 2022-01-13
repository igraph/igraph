/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
/*
   IGraph library.
   Copyright (C) 2005-2021 The igraph development team

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
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/


#include "igraph_paths.h"

#include "igraph_adjlist.h"
#include "igraph_interface.h"
#include "igraph_memory.h"
#include "igraph_qsort.h"

#include "core/indheap.h"
#include "core/interruption.h"

igraph_error_t igraph_get_widest_paths(const igraph_t *graph,
                                       igraph_vector_ptr_t *vertices,
                                       igraph_vector_ptr_t *edges,
                                       igraph_integer_t from,
                                       igraph_vs_t to,
                                       const igraph_vector_t *weights,
                                       igraph_neimode_t mode,
                                       igraph_vector_int_t *predecessors,
                                       igraph_vector_int_t *inbound_edges) {

    /* Implementation details: This is a Dijkstra algorithm with a
    binary heap, modified to support widest paths. The heap is indexed,
    so it stores both the widest path to a node, as well as it's index. We
    use a 2 way heap so that we can query indexes directly in the heap.

    To adapt a Dijkstra to handle widest path, instead of prioritising candidate
    nodes with the minimum distance, we prioritise those with the maximum
    width instead. When adding a node into our set of 'completed' nodes, we
    update all neighbouring nodes with a width that is equal to the min of the
    width to the current node and the width of the edge.

    We denote the widest path from a node to itself as infinity, and the widest
    path from a node to a node it cannot reach as negative infinity.
    */

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_real_t my_posinfinity = IGRAPH_POSINFINITY;
    igraph_real_t my_neginfinity = IGRAPH_NEGINFINITY;
    igraph_vit_t vit;
    igraph_2wheap_t Q;
    igraph_lazy_inclist_t inclist;
    igraph_vector_t widths;
    igraph_integer_t *parents;
    igraph_bool_t *is_target;
    igraph_integer_t i, to_reach;

    if (!weights) {
        IGRAPH_ERROR("Weight vector is required.", IGRAPH_EINVAL);
    }

    if (igraph_vector_size(weights) != no_of_edges) {
        IGRAPH_ERROR("Weight vector length does not match", IGRAPH_EINVAL);
    }

    if (igraph_vector_is_any_nan(weights)) {
        IGRAPH_ERROR("Weight vector must not contain NaN values.", IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_vit_create(graph, to, &vit));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit);

    if (vertices && IGRAPH_VIT_SIZE(vit) != igraph_vector_ptr_size(vertices)) {
        IGRAPH_ERROR("Size of `vertices' and `to' should match.", IGRAPH_EINVAL);
    }
    if (edges && IGRAPH_VIT_SIZE(vit) != igraph_vector_ptr_size(edges)) {
        IGRAPH_ERROR("Size of `edges' and `to' should match.", IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_2wheap_init(&Q, no_of_nodes));
    IGRAPH_FINALLY(igraph_2wheap_destroy, &Q);
    IGRAPH_CHECK(igraph_lazy_inclist_init(graph, &inclist, mode, IGRAPH_LOOPS));
    IGRAPH_FINALLY(igraph_lazy_inclist_destroy, &inclist);

    IGRAPH_VECTOR_INIT_FINALLY(&widths, no_of_nodes);
    igraph_vector_fill(&widths, my_neginfinity);

    parents = IGRAPH_CALLOC(no_of_nodes, igraph_integer_t);
    if (parents == 0) {
        IGRAPH_ERROR("Can't calculate widest paths.", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, parents);
    is_target = IGRAPH_CALLOC(no_of_nodes, igraph_bool_t);
    if (is_target == 0) {
        IGRAPH_ERROR("Can't calculate widest paths.", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, is_target);

    /* Mark the vertices we need to reach */
    to_reach = IGRAPH_VIT_SIZE(vit);
    for (IGRAPH_VIT_RESET(vit); !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit)) {
        if (!is_target[ IGRAPH_VIT_GET(vit) ]) {
            is_target[ IGRAPH_VIT_GET(vit) ] = 1;
        } else {
            to_reach--;       /* this node was given multiple times */
        }
    }

    VECTOR(widths)[from] = my_posinfinity;
    parents[from] = 0;
    igraph_2wheap_push_with_index(&Q, from, my_posinfinity);

    while (!igraph_2wheap_empty(&Q) && to_reach > 0) {
        igraph_integer_t nlen, maxnei = igraph_2wheap_max_index(&Q);
        igraph_real_t maxwidth = igraph_2wheap_delete_max(&Q);
        igraph_vector_int_t *neis;

        IGRAPH_ALLOW_INTERRUPTION();

        if (is_target[maxnei]) {
            is_target[maxnei] = 0;
            to_reach--;
        }

        /* Now check all neighbors of 'maxnei' for a wider path */
        neis = igraph_lazy_inclist_get(&inclist, maxnei);
        nlen = igraph_vector_int_size(neis);
        for (i = 0; i < nlen; i++) {
            igraph_integer_t edge = VECTOR(*neis)[i];
            igraph_integer_t tto = IGRAPH_OTHER(graph, edge, maxnei);
            igraph_real_t edgewidth = VECTOR(*weights)[edge];
            igraph_real_t altwidth = maxwidth < edgewidth ? maxwidth : edgewidth;
            igraph_real_t curdist = VECTOR(widths)[tto];
            if (curdist < 0) {
                /* This is the first assigning a width to this vertex */
                VECTOR(widths)[tto] = altwidth;
                parents[tto] = edge + 1;
                IGRAPH_CHECK(igraph_2wheap_push_with_index(&Q, tto, altwidth));
            } else if (altwidth > curdist) {
                /* This is a wider path */
                VECTOR(widths)[tto] = altwidth;
                parents[tto] = edge + 1;
                IGRAPH_CHECK(igraph_2wheap_modify(&Q, tto, altwidth));
            }
        }
    } /* !igraph_2wheap_empty(&Q) */


    if (to_reach > 0) {
        IGRAPH_WARNING("Couldn't reach some vertices.");
    }

    /* Create `predecessors' if needed */
    if (predecessors) {
        IGRAPH_CHECK(igraph_vector_int_resize(predecessors, no_of_nodes));

        for (i = 0; i < no_of_nodes; i++) {
            if (i == from) {
                /* i is the start vertex */
                VECTOR(*predecessors)[i] = i;
            } else if (parents[i] <= 0) {
                /* i was not reached */
                VECTOR(*predecessors)[i] = -1;
            } else {
                /* i was reached via the edge with ID = parents[i] - 1 */
                VECTOR(*predecessors)[i] = IGRAPH_OTHER(graph, parents[i] - 1, i);
            }
        }
    }

    /* Create `inbound_edges' if needed */
    if (inbound_edges) {
        IGRAPH_CHECK(igraph_vector_int_resize(inbound_edges, no_of_nodes));

        for (i = 0; i < no_of_nodes; i++) {
            if (parents[i] <= 0) {
                /* i was not reached */
                VECTOR(*inbound_edges)[i] = -1;
            } else {
                /* i was reached via the edge with ID = parents[i] - 1 */
                VECTOR(*inbound_edges)[i] = parents[i] - 1;
            }
        }
    }
    /* Reconstruct the widest paths based on vertex and/or edge IDs */
    if (vertices || edges) {
        for (IGRAPH_VIT_RESET(vit), i = 0; !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit), i++) {
            igraph_integer_t node = IGRAPH_VIT_GET(vit);
            igraph_integer_t size, act, edge;
            igraph_vector_int_t *vvec = 0, *evec = 0;
            if (vertices) {
                vvec = VECTOR(*vertices)[i];
                igraph_vector_int_clear(vvec);
            }
            if (edges) {
                evec = VECTOR(*edges)[i];
                igraph_vector_int_clear(evec);
            };
            IGRAPH_ALLOW_INTERRUPTION();

            size = 0;
            act = node;
            while (parents[act]) {
                size++;
                edge = parents[act] - 1;
                act = IGRAPH_OTHER(graph, edge, act);
            }
            if (vvec && (size > 0 || node == from)) {
                IGRAPH_CHECK(igraph_vector_int_resize(vvec, size + 1));
                VECTOR(*vvec)[size] = node;
            }
            if (evec) {
                IGRAPH_CHECK(igraph_vector_int_resize(evec, size));
            }
            act = node;
            while (parents[act]) {
                edge = parents[act] - 1;
                act = IGRAPH_OTHER(graph, edge, act);
                size--;
                if (vvec) {
                    VECTOR(*vvec)[size] = act;
                }
                if (evec) {
                    VECTOR(*evec)[size] = edge;
                }
            }
        }
    }

    igraph_lazy_inclist_destroy(&inclist);
    igraph_2wheap_destroy(&Q);
    igraph_vector_destroy(&widths);
    IGRAPH_FREE(is_target);
    IGRAPH_FREE(parents);
    igraph_vit_destroy(&vit);
    IGRAPH_FINALLY_CLEAN(6);

    return IGRAPH_SUCCESS;
}

igraph_error_t igraph_get_widest_path(const igraph_t *graph,
                                      igraph_vector_int_t *vertices,
                                      igraph_vector_int_t *edges,
                                      igraph_integer_t from,
                                      igraph_integer_t to,
                                      const igraph_vector_t *weights,
                                      igraph_neimode_t mode) {

    igraph_vector_ptr_t vertices2, *vp = &vertices2;
    igraph_vector_ptr_t edges2, *ep = &edges2;

    if (vertices) {
        IGRAPH_CHECK(igraph_vector_ptr_init(&vertices2, 1));
        IGRAPH_FINALLY(igraph_vector_ptr_destroy, &vertices2);
        VECTOR(vertices2)[0] = vertices;
    } else {
        vp = 0;
    }
    if (edges) {
        IGRAPH_CHECK(igraph_vector_ptr_init(&edges2, 1));
        IGRAPH_FINALLY(igraph_vector_ptr_destroy, &edges2);
        VECTOR(edges2)[0] = edges;
    } else {
        ep = 0;
    }

    IGRAPH_CHECK(igraph_get_widest_paths(graph, vp, ep,
                 from, igraph_vss_1(to),
                 weights, mode, 0, 0));

    if (edges) {
        igraph_vector_ptr_destroy(&edges2);
        IGRAPH_FINALLY_CLEAN(1);
    }
    if (vertices) {
        igraph_vector_ptr_destroy(&vertices2);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return IGRAPH_SUCCESS;
}

igraph_error_t igraph_widest_paths_floyd_warshall(const igraph_t *graph,
                                   igraph_matrix_t *res,
                                   const igraph_vs_t from,
                                   const igraph_vs_t to,
                                   const igraph_vector_t *weights,
                                   igraph_neimode_t mode) {

    /* Implementation Details: This is a modified Floyd Warshall algorithm
    which computes the widest path between every pair of nodes. The key
    difference between this and the regular Floyd Warshall is that instead
    of updating the distance between two nodes to be the minimum of itself
    and the distance through an intermediate node, we instead set the width
    to be the maximum of itself and the width through the intermediate node.

    We denote the widest path from a node to itself as infinity, and the widest
    path from a node to a node it cannot reach as negative infinity.
    */

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_lazy_inclist_t inclist;
    igraph_matrix_t adj;
    igraph_integer_t i, j, k;
    igraph_real_t my_posinfinity = IGRAPH_POSINFINITY;
    igraph_real_t my_neginfinity = IGRAPH_NEGINFINITY;
    igraph_vit_t fromvit, tovit;
    igraph_integer_t no_of_from, no_of_to;

    if (!weights) {
        IGRAPH_ERROR("Weight vector is required.", IGRAPH_EINVAL);
    }

    if (igraph_vector_size(weights) != no_of_edges) {
        IGRAPH_ERRORF("Weight vector length (%ld) does not match number "
                      " of edges (%ld).", IGRAPH_EINVAL,
                      igraph_vector_size(weights), no_of_edges);
    }

    if (igraph_vector_is_any_nan(weights)) {
        IGRAPH_ERROR("Weight vector must not contain NaN values.", IGRAPH_EINVAL);
    }

    /* Construct adjacency matrix */
    IGRAPH_CHECK(igraph_matrix_init(&adj, no_of_nodes, no_of_nodes));
    IGRAPH_FINALLY(igraph_matrix_destroy, &adj);

    IGRAPH_CHECK(igraph_lazy_inclist_init(graph, &inclist, mode, IGRAPH_LOOPS));
    IGRAPH_FINALLY(igraph_lazy_inclist_destroy, &inclist);

    for (i = 0; i < no_of_nodes; i++) {
        for (j = 0; j < no_of_nodes; j++) {
            if (i == j) {
                MATRIX(adj, i, j) = my_posinfinity;
            } else {
                MATRIX(adj, i, j) = my_neginfinity;
            }
        }
    }

    for (i = 0; i < no_of_nodes; i++) {
        igraph_vector_int_t *neis = igraph_lazy_inclist_get(&inclist, i);
        igraph_integer_t nlen = igraph_vector_int_size(neis);
        for (j = 0; j < nlen; j++) {
            igraph_integer_t edge = VECTOR(*neis)[j];
            igraph_integer_t tto = IGRAPH_OTHER(graph, edge, i);
            if (VECTOR(*weights)[edge] >  MATRIX(adj, i, tto)) {
                MATRIX(adj, i, tto) = VECTOR(*weights)[edge];
            }
        }
    }

    /* Run modified Floyd Warshall */
    for (k = 0; k < no_of_nodes; k++) {
        for (i = 0; i < no_of_nodes; i++) {
            if (i == k) continue;
            for (j = 0; j < no_of_nodes; j++) {
                if (i == j || j == k) continue;

                IGRAPH_ALLOW_INTERRUPTION();

                igraph_real_t altwidth = MATRIX(adj, i, k);
                if (MATRIX(adj, k, j) < altwidth) {
                    altwidth = MATRIX(adj, k, j);
                }
                if (altwidth > MATRIX(adj, i, j)) {
                    MATRIX(adj, i, j) = altwidth;
                }
            }
        }
    }

    /* Write into results matrix */
    IGRAPH_CHECK(igraph_vit_create(graph, from, &fromvit));
    IGRAPH_FINALLY(igraph_vit_destroy, &fromvit);
    no_of_from = IGRAPH_VIT_SIZE(fromvit);

    IGRAPH_CHECK(igraph_vit_create(graph, to, &tovit));
    IGRAPH_FINALLY(igraph_vit_destroy, &tovit);
    no_of_to = IGRAPH_VIT_SIZE(tovit);

    IGRAPH_CHECK(igraph_matrix_resize(res, no_of_from, no_of_to));

    for (IGRAPH_VIT_RESET(fromvit), i = 0; !IGRAPH_VIT_END(fromvit); IGRAPH_VIT_NEXT(fromvit), i++) {
        igraph_integer_t u = IGRAPH_VIT_GET(fromvit);
        for (IGRAPH_VIT_RESET(tovit),j = 0; !IGRAPH_VIT_END(tovit); IGRAPH_VIT_NEXT(tovit), j++) {
            igraph_integer_t v = IGRAPH_VIT_GET(tovit);
            MATRIX(*res, i, j) = MATRIX(adj, u, v);
        }
    }

    /* Final deletion stuff */
    igraph_vit_destroy(&tovit);
    igraph_vit_destroy(&fromvit);
    igraph_lazy_inclist_destroy(&inclist);
    igraph_matrix_destroy(&adj);
    IGRAPH_FINALLY_CLEAN(4);

    return IGRAPH_SUCCESS;
}

igraph_error_t igraph_widest_paths_dijkstra(const igraph_t *graph,
                                   igraph_matrix_t *res,
                                   const igraph_vs_t from,
                                   const igraph_vs_t to,
                                   const igraph_vector_t *weights,
                                   igraph_neimode_t mode) {

    /* Implementation details: This is a Dijkstra algorithm with a
    binary heap, modified to support widest paths. The heap is indexed,
    so it stores both the widest path to a node, as well as it's index. We
    use a 2 way heap so that we can query indexes directly in the heap.

    To adapt a Dijkstra to handle widest path, instead of prioritising candidate
    nodes with the minimum distance, we prioritise those with the maximum
    width instead. When adding a node into our set of 'completed' nodes, we
    update all neighbouring nodes with a width that is equal to the min of the
    width to the current node and the width of the edge.

    We denote the widest path from a node to itself as infinity, and the widest
    path from a node to a node it cannot reach as negative infinity.
    */

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_2wheap_t Q;
    igraph_vit_t fromvit, tovit;
    igraph_integer_t no_of_from, no_of_to;
    igraph_lazy_inclist_t inclist;
    igraph_integer_t i, j;
    igraph_real_t my_posinfinity = IGRAPH_POSINFINITY;
    igraph_real_t my_neginfinity = IGRAPH_NEGINFINITY;
    igraph_bool_t all_to;
    igraph_vector_int_t indexv;

    if (!weights) {
        IGRAPH_ERROR("Weight vector is required.", IGRAPH_EINVAL);
    }

    if (igraph_vector_size(weights) != no_of_edges) {
        IGRAPH_ERRORF("Weight vector length (%ld) does not match number "
                      " of edges (%ld).", IGRAPH_EINVAL,
                      igraph_vector_size(weights), no_of_edges);
    }

    if (igraph_vector_is_any_nan(weights)) {
        IGRAPH_ERROR("Weight vector must not contain NaN values.", IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_vit_create(graph, from, &fromvit));
    IGRAPH_FINALLY(igraph_vit_destroy, &fromvit);
    no_of_from = IGRAPH_VIT_SIZE(fromvit);

    IGRAPH_CHECK(igraph_2wheap_init(&Q, no_of_nodes));
    IGRAPH_FINALLY(igraph_2wheap_destroy, &Q);
    IGRAPH_CHECK(igraph_lazy_inclist_init(graph, &inclist, mode, IGRAPH_LOOPS));
    IGRAPH_FINALLY(igraph_lazy_inclist_destroy, &inclist);

    all_to = igraph_vs_is_all(&to);
    if (all_to) {
        no_of_to = no_of_nodes;
    } else {
        IGRAPH_VECTOR_INT_INIT_FINALLY(&indexv, no_of_nodes);
        IGRAPH_CHECK(igraph_vit_create(graph, to, &tovit));
        IGRAPH_FINALLY(igraph_vit_destroy, &tovit);
        no_of_to = IGRAPH_VIT_SIZE(tovit);
        for (i = 0; !IGRAPH_VIT_END(tovit); IGRAPH_VIT_NEXT(tovit)) {
            igraph_integer_t v = IGRAPH_VIT_GET(tovit);
            if (VECTOR(indexv)[v]) {
                IGRAPH_ERROR("Duplicate vertices in `to', this is not allowed.",
                             IGRAPH_EINVAL);
            }
            VECTOR(indexv)[v] = ++i;
        }
    }

    IGRAPH_CHECK(igraph_matrix_resize(res, no_of_from, no_of_to));
    igraph_matrix_fill(res, my_neginfinity);

    for (IGRAPH_VIT_RESET(fromvit), i = 0;
         !IGRAPH_VIT_END(fromvit);
         IGRAPH_VIT_NEXT(fromvit), i++) {

        igraph_integer_t reached = 0;
        igraph_integer_t source = IGRAPH_VIT_GET(fromvit);
        igraph_2wheap_clear(&Q);
        igraph_2wheap_push_with_index(&Q, source, my_posinfinity);

        while (!igraph_2wheap_empty(&Q)) {
            igraph_integer_t maxnei = igraph_2wheap_max_index(&Q);
            igraph_real_t maxwidth = igraph_2wheap_deactivate_max(&Q);
            igraph_vector_int_t *neis;
            igraph_integer_t nlen;

            IGRAPH_ALLOW_INTERRUPTION();

            if (all_to) {
                MATRIX(*res, i, maxnei) = maxwidth;
            } else {
                if (VECTOR(indexv)[maxnei]) {
                    MATRIX(*res, i, VECTOR(indexv)[maxnei] - 1) = maxwidth;
                    reached++;
                    if (reached == no_of_to) {
                        igraph_2wheap_clear(&Q);
                        break;
                    }
                }
            }

            /* Now check all neighbors of 'maxnei' for a wider path*/
            neis = igraph_lazy_inclist_get(&inclist, maxnei);
            nlen = igraph_vector_int_size(neis);
            for (j = 0; j < nlen; j++) {
                igraph_integer_t edge = VECTOR(*neis)[j];
                igraph_integer_t tto = IGRAPH_OTHER(graph, edge, maxnei);
                igraph_real_t edgewidth = VECTOR(*weights)[edge];
                igraph_real_t altwidth = maxwidth < edgewidth ? maxwidth : edgewidth;
                igraph_bool_t active = igraph_2wheap_has_active(&Q, tto);
                igraph_bool_t has = igraph_2wheap_has_elem(&Q, tto);
                igraph_real_t curwidth = active ? igraph_2wheap_get(&Q, tto) : my_posinfinity;
                if (!has) {
                    /* This is the first time assigning a width to this vertex */
                    IGRAPH_CHECK(igraph_2wheap_push_with_index(&Q, tto, altwidth));
                } else if (altwidth > curwidth) {
                    /* This is a wider path */
                    IGRAPH_CHECK(igraph_2wheap_modify(&Q, tto, altwidth));
                }
            }

        } /* !igraph_2wheap_empty(&Q) */

    } /* !IGRAPH_VIT_END(fromvit) */

    if (!all_to) {
        igraph_vit_destroy(&tovit);
        igraph_vector_int_destroy(&indexv);
        IGRAPH_FINALLY_CLEAN(2);
    }

    igraph_lazy_inclist_destroy(&inclist);
    igraph_2wheap_destroy(&Q);
    igraph_vit_destroy(&fromvit);
    IGRAPH_FINALLY_CLEAN(3);

    return IGRAPH_SUCCESS;
}
