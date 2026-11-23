/*
   igraph library.
   Copyright (C) 2025  The igraph development team <igraph@igraph.org>

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

#include "igraph_structural.h"

#include "igraph_interface.h"
#include "igraph_isomorphism.h"

/* Returns the total number of possible edges given the number of vertices in a graph,
 * whether it is directed, and whether loops should be assumed possible.
 *
 * No loops, undirected: n * (n-1) / 2
 * No loops, directed:   n * (n-1)
 * Loops, undirected:    n * (n+1) / 2
 * Loops, directed:      n^2
 */
static igraph_real_t total_possible_edges(igraph_int_t vcount,
                                          igraph_bool_t directed,
                                          igraph_bool_t loops) {
    igraph_real_t nv = (igraph_real_t) vcount;
    if (!loops) {
        return (directed ? nv * (nv - 1) : nv * (nv - 1) / 2);
    } else {
        return (directed ? nv * nv : nv * (nv + 1) / 2);
    }
}

/**
 * \function igraph_rich_club_sequence
 * \brief Density sequence of subgraphs formed by sequential vertex removal.
 *
 * \experimental
 *
 * This function takes a graph and a vertex ordering as input, sequentially
 * removes the vertices in the given order, and calculates the density of the
 * remaining subgraph after each removal.
 *
 * </para><para>
 * Density is calculated as the ratio of the number of edges (or total edge
 * weight, if weighted) to the number of total possible edges in the graph.
 * The latter is dependent on whether the graph is directed and whether
 * self-loops are assumed to be possible: for undirected graphs without
 * self-loops, this total is given by <code>n(n-1)/2</code>,
 * and for directed graphs by <code>n(n-1)</code>.
 * When self-loops are allowed, these are adjusted to <code>n(n+1)/2</code>
 * for undirected and <code>n^2</code> for directed graphs.
 *
 * </para><para>
 * Vertex order can be sorted by degree so that the resulting density sequence
 * helps reveal how interconnected a graph is across different degree levels,
 * or the presence of a "rich-club" effect.
 *
 * \param graph The graph object to analyze.
 * \param weights Vector of edge weights. If \c NULL all weights are
 *    assumed to be 1.
 * \param res Initialized vector, the result will be written here. <code>res[i]</code>
 *    contains the density of the remaining graph after \c i vertices have been
 *    removed. If \p normalized is set to \c false, it contains the remaining
 *    edge count (or remaining total edge weights if weights were given).
 * \param vertex_order Vector giving the order in which vertices are removed.
 * \param normalized If \c false, return edge counts (or total edge weights).
 *    If \c true, divide by the largest possible edge count to obtain densities.
 * \param loops Whether self-loops are assumed to be possible. Ignored when
 *    normalized is not requested.
 * \param directed If false, directed graphs will be treated as undirected.
 *    Ignored with undirected graphs.
 *
 * \return Error code: \c IGRAPH_EINVAL: invalid vertex_order vector and/or weight vector
 * lengths
 *
 * Time complexity: O(|V| + |E|)
 * where |V| is the number of vertices and |E| the number of edges in the graph given.
 *
 * \sa \ref igraph_density(), which uses the same calculation of total possible edges.
 */
igraph_error_t igraph_rich_club_sequence(
        const igraph_t *graph, const igraph_vector_t *weights,
        igraph_vector_t *res,
        const igraph_vector_int_t *vertex_order,
        igraph_bool_t normalized,
        igraph_bool_t loops, igraph_bool_t directed) {

    const igraph_int_t vcount = igraph_vcount(graph);
    const igraph_int_t ecount = igraph_ecount(graph);
    igraph_vector_int_t order_of;

    // Error handling: invalid vertex_order and weights sizes
    if (igraph_vector_int_size(vertex_order) != vcount) {
        IGRAPH_ERRORF("Vertex order vector length (%" IGRAPH_PRId ") does not match "
                      "number of vertices (%" IGRAPH_PRId ").",
                      IGRAPH_EINVAL,
                      igraph_vector_int_size(vertex_order), vcount);
    }
    if (weights && igraph_vector_size(weights) != ecount) {
        IGRAPH_ERRORF("Weight vector length (%" IGRAPH_PRId ") does not match "
                      "number of edges (%" IGRAPH_PRId ").",
                      IGRAPH_EINVAL,
                      igraph_vector_size(weights), ecount);
    }

    if (! igraph_is_directed(graph)) {
        directed = false;
    }

    IGRAPH_CHECK(igraph_vector_resize(res, vcount));
    igraph_vector_null(res);

    IGRAPH_VECTOR_INT_INIT_FINALLY(&order_of, vcount);
    IGRAPH_CHECK(igraph_invert_permutation(vertex_order, &order_of));

    igraph_bool_t warning_issued = false;

    // remaining_total_weight vector: number of edges (or total edge weight) removed by index
    for (igraph_int_t eid = 0; eid < ecount; eid++) {
        igraph_int_t v1 = IGRAPH_FROM(graph, eid);
        igraph_int_t v2 = IGRAPH_TO(graph, eid);
        if (!loops && normalized && !warning_issued && v1 == v2) {
            IGRAPH_WARNING("Self-loops were requested to be assumed absent, "
                           "but a self-loop was encountered. Density calculations will proceed "
                           "with the assumption of no loops.");
            warning_issued = true;
        }
        igraph_int_t order_v1 = VECTOR(order_of)[v1]; // order of endpoints
        igraph_int_t order_v2 = VECTOR(order_of)[v2];

        igraph_int_t edge_removal_index = (order_v1 < order_v2 ? order_v1 : order_v2);
        VECTOR(*res)[edge_removal_index] += (weights ? VECTOR(*weights)[eid] : 1);
    }

    // remaining_total_weight vector: edges (or total edge weight) remaining after i removals
    igraph_real_t total = 0;
    for (igraph_int_t i = vcount - 1; i >= 0; i--) {
        total += VECTOR(*res)[i];
        VECTOR(*res)[i] = total;
    }

    // Normalize edge counts to densities
    if (normalized) {
        for (igraph_int_t i = 0; i < vcount; i++) {
            // (vcount - i) = the number of vertices left in this loop
            VECTOR(*res)[i] =
                    VECTOR(*res)[i] /
                    total_possible_edges(vcount - i, directed, loops);
        }
    }

    igraph_vector_int_destroy(&order_of);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}


/**
 * \function igraph_rich_club_coefficient
 * \brief Degree-based rich-club coefficient of a graph.
 *
 * \experimental
 *
 * The rich-club coefficient is a metric used to evaluate whether well-connected
 * vertices (hubs) tend to connect to other well-connected vertices.
 *
 * </para><para>
 * For a given degree threshold \p k, the rich-club coefficient \phi(k) is defined
 * as the density of the subgraph induced by vertices with degree strictly greater
 * than \p k.
 *
 * </para><para>
 * The output vector \p res is resized to length <code>max(degree) + 1</code>,
 * where <code>res[k]</code> contains \phi(k). When no vertices have degree
 * greater than \p k, or fewer than two vertices remain (when \p loops is false),
 * <code>res[k]</code> evaluates to \c NaN.
 *
 * \param graph The graph object to analyze.
 * \param weights Optional edge weights. If \c NULL, all edge weights are
 *    assumed to be 1.
 * \param res Initialized vector, the result will be written here. <code>res[k]</code>
 *    contains the rich-club coefficient for degree threshold \p k.
 * \param mode For directed graphs, selects whether out-degree (\c IGRAPH_OUT),
 *    in-degree (\c IGRAPH_IN), or total degree (\c IGRAPH_ALL) determines vertex
 *    degree for ranking. Ignored for undirected graphs.
 * \param loops Whether self-loops are considered when calculating degree and whether
 *    they are assumed to be possible in the theoretical denominator.
 * \param directed Whether directed density is calculated. Ignored for undirected graphs.
 *
 * \return Error code:
 *    \c IGRAPH_EINVAL: invalid weight vector length.
 *    \c IGRAPH_EINVMODE: invalid mode for directed graph.
 *
 * Time complexity: O(|V| log |V| + |E|), where |V| is the number of vertices and |E|
 * is the number of edges.
 *
 * \sa \ref igraph_rich_club_sequence(), \ref igraph_density().
 */
igraph_error_t igraph_rich_club_coefficient(
        const igraph_t *graph,
        const igraph_vector_t *weights,
        igraph_vector_t *res,
        igraph_neimode_t mode,
        igraph_bool_t loops,
        igraph_bool_t directed) {

    const igraph_int_t vcount = igraph_vcount(graph);
    const igraph_int_t ecount = igraph_ecount(graph);

    if (weights && igraph_vector_size(weights) != ecount) {
        IGRAPH_ERRORF("Weight vector length (%" IGRAPH_PRId ") does not match "
                      "number of edges (%" IGRAPH_PRId ").",
                      IGRAPH_EINVAL,
                      igraph_vector_size(weights), ecount);
    }

    if (igraph_is_directed(graph)) {
        if (mode != IGRAPH_OUT && mode != IGRAPH_IN && mode != IGRAPH_ALL) {
            IGRAPH_ERROR("Invalid mode for directed graph.", IGRAPH_EINVMODE);
        }
    } else {
        directed = false;
        mode = IGRAPH_ALL;
    }

    /* Empty graph case */
    if (vcount == 0) {
        IGRAPH_CHECK(igraph_vector_resize(res, 0));
        return IGRAPH_SUCCESS;
    }

    igraph_vector_int_t degrees;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&degrees, vcount);
    igraph_loops_t deg_loops = loops ? (igraph_is_directed(graph) ? IGRAPH_LOOPS_ONCE : IGRAPH_LOOPS_TWICE)
                                     : IGRAPH_NO_LOOPS;
    IGRAPH_CHECK(igraph_degree(graph, &degrees, igraph_vss_all(), mode, deg_loops));

    igraph_int_t maxdeg = igraph_vector_int_max(&degrees);
    IGRAPH_CHECK(igraph_vector_resize(res, maxdeg + 1));
    igraph_vector_fill(res, IGRAPH_NAN);

    igraph_vector_int_t order;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&order, vcount);
    IGRAPH_CHECK(igraph_vector_int_sort_ind(&degrees, &order, IGRAPH_ASCENDING));

    igraph_vector_t seq_res;
    IGRAPH_VECTOR_INIT_FINALLY(&seq_res, vcount);
    IGRAPH_CHECK(igraph_rich_club_sequence(graph, weights, &seq_res, &order, true, loops, directed));

    igraph_int_t m_k = 0;
    for (igraph_int_t k = 0; k <= maxdeg; k++) {
        while (m_k < vcount && VECTOR(degrees)[VECTOR(order)[m_k]] <= k) {
            m_k++;
        }
        igraph_int_t remaining = vcount - m_k;
        if (remaining >= 2 || (loops && remaining >= 1)) {
            VECTOR(*res)[k] = VECTOR(seq_res)[m_k];
        } else {
            VECTOR(*res)[k] = IGRAPH_NAN;
        }
    }

    igraph_vector_destroy(&seq_res);
    igraph_vector_int_destroy(&order);
    igraph_vector_int_destroy(&degrees);
    IGRAPH_FINALLY_CLEAN(3);

    return IGRAPH_SUCCESS;
}
