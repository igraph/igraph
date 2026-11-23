/*
   igraph library.
   Copyright (C) 2026  The igraph development team <igraph@igraph.org>

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

#include "igraph_path_cover.h"

#include "igraph_adjlist.h"
#include "igraph_dqueue.h"
#include "igraph_error.h"
#include "igraph_interface.h"

#include "flow/path_cover_internal.h"

/**
 * \function igraph_maximum_antichain
 * \brief Maximum antichain of a DAG, extracted from a solved minimum path cover flow.
 *
 * By Dilworth's theorem, the width of a minimum path cover of a DAG equals
 * the size of a maximum antichain of its reachability order (a set of
 * vertices no two of which are connected by a directed path). This
 * function extracts such a maximum antichain from the auxiliary flow
 * network produced by \ref igraph_minimum_path_cover(), without solving a
 * separate flow problem.
 *
 * </para><para>
 * The \p flow and \p demand vectors must come from a \em solved (i.e.
 * minimum, not merely feasible) flow, such as the ones returned by
 * \ref igraph_minimum_path_cover() through its \c flow_network_flow and
 * \c flow_network_demand parameters. Passing a flow that has not been
 * fully minimized results in undefined behavior.
 *
 * \param flow_network The auxiliary flow network, as returned by
 *        \ref igraph_minimum_path_cover() through its \c flow_network
 *        parameter.
 * \param flow The flow on each edge of \p flow_network, as returned by
 *        \ref igraph_minimum_path_cover() through its
 *        \c flow_network_flow parameter.
 * \param demand The demand (lower bound) on each edge of \p flow_network,
 *        as returned by \ref igraph_minimum_path_cover() through its
 *        \c flow_network_demand parameter.
 * \param antichain Pointer to an initialized vector. The vertex IDs of a
 *        maximum antichain of the original graph are stored here.
 * \return Error code.
 *
 * Time complexity: O(|V|+|E|), where |V| and |E| are the number of
 * vertices and edges in \p flow_network.
 *
 * </para><para>
 * The two-phase residual-reachability extraction below is ported from
 * maxantichain_from_minflow() in
 * https://github.com/algbio/PerformanceMPC/blob/main/src/mpc/antichain.cpp.
 *
 * </para><para>
 * Reference:
 *
 * </para><para>
 * Dilworth RP: A Decomposition Theorem for Partially Ordered Sets.
 * Annals of Mathematics, Second Series, 51(1):161-166, 1950.
 *
 * \sa \ref igraph_minimum_path_cover().
 */
igraph_error_t igraph_maximum_antichain(
        const igraph_t *flow_network,
        const igraph_vector_int_t *flow,
        const igraph_vector_int_t *demand,
        igraph_vector_int_t *antichain) {

    igraph_int_t network_nodes = igraph_vcount(flow_network);
    igraph_int_t no_of_nodes;
    igraph_int_t source, sink;
    igraph_inclist_t out_inc, in_inc;
    igraph_vector_bool_t reachable, claimed;
    igraph_dqueue_int_t queue;
    igraph_int_t i;

    if (network_nodes < 2 || (network_nodes - 2) % 2 != 0) {
        IGRAPH_ERROR("Invalid flow network; it must have been produced by "
                     "igraph_minimum_path_cover().", IGRAPH_EINVAL);
    }
    if (igraph_vector_int_size(flow) != igraph_ecount(flow_network) ||
        igraph_vector_int_size(demand) != igraph_ecount(flow_network)) {
        IGRAPH_ERROR("The flow and demand vectors must have one entry per edge "
                     "of the flow network.", IGRAPH_EINVAL);
    }

    no_of_nodes = (network_nodes - 2) / 2;
    source = 2 * no_of_nodes;
    sink = 2 * no_of_nodes + 1;

    IGRAPH_CHECK(igraph_inclist_init(flow_network, &out_inc, IGRAPH_OUT, IGRAPH_LOOPS_TWICE));
    IGRAPH_FINALLY(igraph_inclist_destroy, &out_inc);
    IGRAPH_CHECK(igraph_inclist_init(flow_network, &in_inc, IGRAPH_IN, IGRAPH_LOOPS_TWICE));
    IGRAPH_FINALLY(igraph_inclist_destroy, &in_inc);

    IGRAPH_VECTOR_BOOL_INIT_FINALLY(&reachable, network_nodes);
    IGRAPH_VECTOR_BOOL_INIT_FINALLY(&claimed, network_nodes);

    IGRAPH_DQUEUE_INT_INIT_FINALLY(&queue, 0);
    IGRAPH_CHECK(igraph_dqueue_int_push(&queue, source));
    VECTOR(reachable)[source] = true;

    /* Phase 1: mark everything reachable from the source in the residual
     * sense: forward along edges with slack, backward along any edge. */
    while (!igraph_dqueue_int_empty(&queue)) {
        igraph_int_t s = igraph_dqueue_int_pop(&queue);
        igraph_vector_int_t *neis;
        igraph_int_t nlen, j;

        neis = igraph_inclist_get(&out_inc, s);
        nlen = igraph_vector_int_size(neis);
        for (j = 0; j < nlen; j++) {
            igraph_int_t eid = VECTOR(*neis)[j];
            igraph_int_t u = IGRAPH_TO(flow_network, eid);
            if (!VECTOR(reachable)[u] && VECTOR(*flow)[eid] > VECTOR(*demand)[eid]) {
                VECTOR(reachable)[u] = true;
                IGRAPH_CHECK(igraph_dqueue_int_push(&queue, u));
            }
        }

        neis = igraph_inclist_get(&in_inc, s);
        nlen = igraph_vector_int_size(neis);
        for (j = 0; j < nlen; j++) {
            igraph_int_t eid = VECTOR(*neis)[j];
            igraph_int_t u = IGRAPH_FROM(flow_network, eid);
            if (!VECTOR(reachable)[u]) {
                VECTOR(reachable)[u] = true;
                IGRAPH_CHECK(igraph_dqueue_int_push(&queue, u));
            }
        }
    }

    /* This must hold for a genuinely minimum flow: the sink can never be
     * reachable from the source in the residual graph, by the usual
     * max-flow/min-cut termination argument. */
    IGRAPH_ASSERT(!VECTOR(reachable)[sink]);

    /* Phase 2: the antichain consists of the source-side endpoints of
     * tight, positive-demand edges (i.e. saturated vertex-demand edges)
     * crossing from the reachable side to the unreachable side. */
    igraph_vector_int_clear(antichain);
    for (i = 0; i < network_nodes; i++) {
        igraph_vector_int_t *neis;
        igraph_int_t nlen, j;

        if (!VECTOR(reachable)[i]) {
            continue;
        }

        neis = igraph_inclist_get(&out_inc, i);
        nlen = igraph_vector_int_size(neis);
        for (j = 0; j < nlen; j++) {
            igraph_int_t eid = VECTOR(*neis)[j];
            igraph_int_t u = IGRAPH_TO(flow_network, eid);
            if (VECTOR(*flow)[eid] == VECTOR(*demand)[eid] && VECTOR(*demand)[eid] >= 1 &&
                !VECTOR(reachable)[u] && !VECTOR(claimed)[u]) {
                IGRAPH_CHECK(igraph_vector_int_push_back(antichain, igraph_i_mpc_v_r(i, no_of_nodes)));
                VECTOR(claimed)[u] = true;
            }
        }
    }

    igraph_dqueue_int_destroy(&queue);
    igraph_vector_bool_destroy(&claimed);
    igraph_vector_bool_destroy(&reachable);
    igraph_inclist_destroy(&in_inc);
    igraph_inclist_destroy(&out_inc);
    IGRAPH_FINALLY_CLEAN(5);

    return IGRAPH_SUCCESS;
}
