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

    // TODO: implement a modified dijkstra to calculate widest paths from a source to every other node

    return IGRAPH_UNIMPLEMENTED;
}

igraph_error_t igraph_get_widest_path(const igraph_t *graph,
                                      igraph_vector_int_t *vertices,
                                      igraph_vector_int_t *edges,
                                      igraph_integer_t from,
                                      igraph_integer_t to,
                                      const igraph_vector_t *weights,
                                      igraph_neimode_t mode) {

    // TODO: calls igraph_get_widest_paths and returns the one specific path to "to"

    return IGRAPH_UNIMPLEMENTED;
}

igraph_error_t igraph_get_all_widest_paths(const igraph_t *graph,
                                            igraph_vector_ptr_t *vertices,
                                            igraph_vector_ptr_t *edges,
                                            igraph_vector_int_t *nrgeo,
                                            igraph_integer_t from, igraph_vs_t to,
                                            const igraph_vector_t *weights,
                                            igraph_neimode_t mode) {

    // TODO: implement a modified dijkstra that stores all paths

    return IGRAPH_UNIMPLEMENTED;
}



igraph_error_t igraph_widest_paths_floyd_warshall(const igraph_t *graph,
                                   igraph_matrix_t *res,
                                   const igraph_vs_t from,
                                   const igraph_vs_t to,
                                   const igraph_vector_t *weights,
                                   igraph_neimode_t mode) {

    // TODO: implement floyd warshalls to calculate widest paths between all nodes

    return IGRAPH_UNIMPLEMENTED;
}

igraph_error_t igraph_widest_paths_dijkstra(const igraph_t *graph,
                                   igraph_matrix_t *res,
                                   const igraph_vs_t from,
                                   const igraph_vs_t to,
                                   const igraph_vector_t *weights,
                                   igraph_neimode_t mode) {

    /* Implementation details: This is a copy of the igraph_shortest_paths_dijkstra() algorithm,
    with some small modifications. Namely:

        - We prioritise nodes with the widest path so far in the heap.
        - When adding nodes into the heap, we set the value to be equal to the min of the current
          widest path and the weight of the edge.
        - We allow negative weights.
        - The widest path from a node to itself has width infinity.
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
        return igraph_shortest_paths(graph, res, from, to, mode);
    }

    if (igraph_vector_size(weights) != no_of_edges) {
        IGRAPH_ERRORF("Weight vector length (%ld) does not match number "
                      " of edges (%ld).", IGRAPH_EINVAL,
                      igraph_vector_size(weights), no_of_edges);
    }

    if (no_of_edges > 0) {
        igraph_real_t min = igraph_vector_min(weights);
        if (igraph_is_nan(min)) {
            IGRAPH_ERROR("Weight vector must not contain NaN values.", IGRAPH_EINVAL);
        }
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
                IGRAPH_ERROR("Duplicate vertices in `to', this is not allowed",
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
                igraph_real_t edgeweight = VECTOR(*weights)[edge];
                igraph_real_t altwidth = maxwidth < edgeweight ? maxwidth : edgeweight;
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
