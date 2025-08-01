/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
/*
   IGraph library.
   Copyright (C) 2007-2020 The igraph development team

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

#include "igraph_community.h"

#include "igraph_interface.h"
#include "igraph_structural.h"

/**
 * \function igraph_modularity
 * \brief Calculates the modularity of a graph with respect to some clusters or vertex types.
 *
 * The modularity of a graph with respect to some clustering of the vertices
 * (or assignment of vertex types)
 * measures how strongly separated the different clusters are from each
 * other compared to a random null model. It is defined as
 *
 * </para><para>
 * <code>Q = 1/(2m) sum_ij (A_ij - γ k_i k_j / (2m)) δ(c_i,c_j)</code>,
 *
 * </para><para>
 * where \c m is the number of edges, <code>A_ij</code> is the adjacency matrix,
 * \c k_i is the degree of vertex \c i, \c c_i is the cluster that vertex \c i belongs to
 * (or its vertex type), <code>δ(i,j)=1</code> if <code>i=j</code> and 0 otherwise,
 * and the sum goes over all \c i, \c j pairs of vertices. Note that in this formula,
 * the diagonal of the adjacency matrix contains twice the number of self-loops.
 *
 * </para><para>
 * The resolution parameter \c γ allows weighting the random null model, which
 * might be useful when finding partitions with a high modularity. Maximizing modularity
 * with higher values of the resolution parameter typically results in more, smaller clusters
 * when finding partitions with a high modularity. Lower values typically results in
 * fewer, larger clusters. The original definition of modularity is retrieved
 * when setting <code>γ = 1</code>.
 *
 * </para><para>
 * Modularity can also be calculated on directed graphs. This only requires a relatively
 * modest change,
 *
 * </para><para>
 * <code>Q = 1/m sum_ij (A_ij - γ k^out_i k^in_j / m) δ(c_i,c_j)</code>,
 *
 * </para><para>
 * where \c k^out_i is the out-degree of node \c i and \c k^in_j is the in-degree of node \c j.
 *
 * </para><para>
 * Modularity on weighted graphs is also meaningful. When taking
 * edge weights into account, \c A_ij equals the weight of the corresponding edge
 * (or 0 if there is no edge), \c k_i is the strength (i.e. the weighted degree) of
 * vertex \c i, with similar counterparts for a directed graph, and \c m is the total
 * weight of all edges.
 *
 * </para><para>
 * Note that the modularity is not well-defined for graphs with no edges.
 * igraph returns \c NaN for graphs with no edges; see
 * https://github.com/igraph/igraph/issues/1539 for
 * a detailed discussion.
 *
 * </para><para>
 * For the original definition of modularity, see Newman, M. E. J., and Girvan, M.
 * (2004). Finding and evaluating community structure in networks.
 * Physical Review E 69, 026113. https://doi.org/10.1103/PhysRevE.69.026113
 *
 * </para><para>
 * For the directed definition of modularity, see Leicht, E. A., and Newman, M. E.
 * J. (2008). Community Structure in Directed Networks. Physical Review Letters 100,
 * 118703. https://doi.org/10.1103/PhysRevLett.100.118703
 *
 * </para><para>
 * For the introduction of the resolution parameter \c γ, see Reichardt, J., and
 * Bornholdt, S. (2006). Statistical mechanics of community detection. Physical
 * Review E 74, 016110. https://doi.org/10.1103/PhysRevE.74.016110
 *
 * \param graph      The input graph.
 * \param membership Numeric vector of integer values which gives the type of each
 *                   vertex, i.e. the cluster to which it belongs.
 *                   It does not have to be consecutive, i.e. empty communities
 *                   are allowed.
 * \param weights    Weight vector or \c NULL if no weights are specified.
 * \param resolution The resolution parameter \c γ. Must not be negative.
 *                   Set it to 1 to use the classical definition of modularity.
 * \param directed   Whether to use the directed or undirected version of modularity.
 *                   Ignored for undirected graphs.
 * \param modularity Pointer to a real number, the result will be
 *                   stored here.
 * \return Error code.
 *
 * \sa \ref igraph_modularity_matrix()
 *
 * Time complexity: O(|V|+|E|), the number of vertices plus the number
 * of edges.
 */
igraph_error_t igraph_modularity(const igraph_t *graph,
                      const igraph_vector_int_t *membership,
                      const igraph_vector_t *weights,
                      const igraph_real_t resolution,
                      const igraph_bool_t directed,
                      igraph_real_t *modularity) {

    igraph_vector_t k_out, k_in;
    igraph_integer_t no_of_partitions;
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_real_t e; /* count/fraction of edges/weights within partitions */
    igraph_real_t m; /* edge count / weight sum */
    igraph_integer_t c1, c2;
    /* Only consider the graph as directed if it actually is directed */
    igraph_bool_t use_directed = directed && igraph_is_directed(graph);
    igraph_real_t directed_multiplier = (use_directed ? 1 : 2);

    if (igraph_vector_int_size(membership) != igraph_vcount(graph)) {
        IGRAPH_ERROR("Membership vector size differs from number of vertices.",
                     IGRAPH_EINVAL);
    }
    if (resolution < 0.0) {
      IGRAPH_ERROR("The resolution parameter must not be negative.", IGRAPH_EINVAL);
    }

    if (no_of_edges == 0) {
        /* Special case: the modularity of graphs with no edges is not
         * well-defined */
        *modularity = IGRAPH_NAN;
        return IGRAPH_SUCCESS;
    }

    /* At this point, the 'membership' vector does not have length zero,
       thus it is safe to call igraph_vector_max() and min(). */

    no_of_partitions = igraph_vector_int_max(membership) + 1;

    if (igraph_vector_int_min(membership) < 0) {
        IGRAPH_ERROR("Invalid membership vector: negative entry.", IGRAPH_EINVAL);
    }

    IGRAPH_VECTOR_INIT_FINALLY(&k_out, no_of_partitions);
    IGRAPH_VECTOR_INIT_FINALLY(&k_in, no_of_partitions);

    e = 0.0;
    if (weights) {
        if (igraph_vector_size(weights) != no_of_edges)
            IGRAPH_ERROR("Weight vector size differs from number of edges.",
                         IGRAPH_EINVAL);
        m = 0.0;
        for (igraph_integer_t i = 0; i < no_of_edges; i++) {
            igraph_real_t w = VECTOR(*weights)[i];
            if (w < 0) {
                IGRAPH_ERROR("Negative weight in weight vector.", IGRAPH_EINVAL);
            }
            c1 = VECTOR(*membership)[ IGRAPH_FROM(graph, i) ];
            c2 = VECTOR(*membership)[ IGRAPH_TO(graph, i) ];
            if (c1 == c2) {
                e += directed_multiplier * w;
            }
            VECTOR(k_out)[c1] += w;
            VECTOR(k_in)[c2]  += w;
            m += w;
        }
    } else {
        m = no_of_edges;
        for (igraph_integer_t i = 0; i < no_of_edges; i++) {
            c1 = VECTOR(*membership)[ IGRAPH_FROM(graph, i) ];
            c2 = VECTOR(*membership)[ IGRAPH_TO(graph, i) ];
            if (c1 == c2) {
                e += directed_multiplier;
            }
            VECTOR(k_out)[c1] += 1;
            VECTOR(k_in)[c2]  += 1;
        }
    }

    if (!use_directed) {
        /* Graph is undirected, simply add vectors */
        igraph_vector_add(&k_out, &k_in);
        igraph_vector_update(&k_in, &k_out);
    }

    /* Divide all vectors by total weight. */
    igraph_vector_scale(&k_out, 1.0/( directed_multiplier * m ) );
    igraph_vector_scale(&k_in, 1.0/( directed_multiplier * m ) );
    e /= directed_multiplier * m;

    if (m > 0) {
        *modularity = e;
        for (igraph_integer_t i = 0; i < no_of_partitions; i++) {
            *modularity -= resolution * VECTOR(k_out)[i] * VECTOR(k_in)[i];
        }
    } else {
        *modularity = IGRAPH_NAN;
    }

    igraph_vector_destroy(&k_out);
    igraph_vector_destroy(&k_in);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_modularity_matrix_get_adjacency(
           const igraph_t *graph, igraph_matrix_t *res,
           const igraph_vector_t *weights, igraph_bool_t directed) {

    /* Specifically used to handle weights and/or ignore direction */
    igraph_eit_t edgeit;
    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t from, to;

    IGRAPH_CHECK(igraph_matrix_resize(res, no_of_nodes, no_of_nodes));
    igraph_matrix_null(res);
    IGRAPH_CHECK(igraph_eit_create(graph, igraph_ess_all(IGRAPH_EDGEORDER_ID), &edgeit));
    IGRAPH_FINALLY(igraph_eit_destroy, &edgeit);

    if (weights) {
        for (; !IGRAPH_EIT_END(edgeit); IGRAPH_EIT_NEXT(edgeit)) {
            igraph_integer_t edge = IGRAPH_EIT_GET(edgeit);
            from = IGRAPH_FROM(graph, edge);
            to = IGRAPH_TO(graph, edge);
            MATRIX(*res, from, to) += VECTOR(*weights)[edge];
            if (!directed) {
                MATRIX(*res, to, from) += VECTOR(*weights)[edge];
            }
        }
    } else {
        for (; !IGRAPH_EIT_END(edgeit); IGRAPH_EIT_NEXT(edgeit)) {
            igraph_integer_t edge = IGRAPH_EIT_GET(edgeit);
            from = IGRAPH_FROM(graph, edge);
            to = IGRAPH_TO(graph, edge);
            MATRIX(*res, from, to) += 1;
            if (!directed) {
                MATRIX(*res, to, from) += 1;
            }
        }
    }

    igraph_eit_destroy(&edgeit);
    IGRAPH_FINALLY_CLEAN(1);
    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_modularity_matrix
 * \brief Calculates the modularity matrix.
 *
 * This function returns the modularity matrix, which is defined as
 *
 * </para><para>
 * <code>B_ij = A_ij - γ k_i k_j / (2m)</code>
 *
 * </para><para>
 * for undirected graphs, where \c A_ij is the adjacency matrix, \c γ is the
 * resolution parameter, \c k_i is the degree of vertex \c i, and \c m is the
 * number of edges in the graph. When there are no edges, or the weights add up
 * to zero, the result is undefined.
 *
 * </para><para>
 * For directed graphs the modularity matrix is changed to
 *
 * </para><para>
 * <code>B_ij = A_ij - γ k^out_i k^in_j / m</code>
 *
 * </para><para>
 * where <code>k^out_i</code> is the out-degree of node \c i and <code>k^in_j</code> is the
 * in-degree of node \c j.
 *
 * </para><para>
 * Note that self-loops in undirected graphs are multiplied by 2 in this
 * implementation. If weights are specified, the weighted counterparts of the adjacency
 * matrix and degrees are used.
 *
 * \param graph      The input graph.
 * \param weights    Edge weights, pointer to a vector. If this is a null pointer
 *                   then every edge is assumed to have a weight of 1.
 * \param resolution The resolution parameter \c γ. Must not be negative.
 *                   Default is 1. Lower values favor fewer, larger communities;
 *                   higher values favor more, smaller communities.
 * \param modmat     Pointer to an initialized matrix in which the modularity
 *                   matrix is stored.
 * \param directed   For directed graphs: if the edges should be treated as
 *                   undirected. For undirected graphs this is ignored.
 * \return Error code.
 *
 * \sa \ref igraph_modularity()
 */
igraph_error_t igraph_modularity_matrix(const igraph_t *graph,
                             const igraph_vector_t *weights,
                             const igraph_real_t resolution,
                             igraph_matrix_t *modmat,
                             igraph_bool_t directed) {

    const igraph_integer_t no_of_nodes = igraph_vcount(graph);
    const igraph_integer_t no_of_edges = igraph_ecount(graph);
    const igraph_real_t sw = weights ? igraph_vector_sum(weights) : no_of_edges;
    igraph_vector_t deg, in_deg, out_deg;
    igraph_real_t scaling_factor;

    if (weights && igraph_vector_size(weights) != no_of_edges) {
        IGRAPH_ERROR("Invalid weight vector length.", IGRAPH_EINVAL);
    }

    if (resolution < 0.0) {
        IGRAPH_ERROR("The resolution parameter must not be negative.", IGRAPH_EINVAL);
    }

    if (!igraph_is_directed(graph)) {
        directed = false;
    }
    IGRAPH_CHECK(igraph_i_modularity_matrix_get_adjacency(graph, modmat, weights, directed));

    /* Performance notes:
     *  - Iterating in column-major order makes a large difference.
     *  - Applying the scaling_factor to in_deg (or out_deg) first to reduce the
     *    number of multiplications does not make an appreciable performance
     *    difference. However, doing this in the undirected case causes the result
     *    matrix to sometimes not be strictly symmetric due to the non-associativity
     *    of floating point multiplication.
     */

    if (directed) {
        IGRAPH_VECTOR_INIT_FINALLY(&in_deg, no_of_nodes);
        IGRAPH_VECTOR_INIT_FINALLY(&out_deg, no_of_nodes);

        IGRAPH_CHECK(igraph_strength(graph, &in_deg, igraph_vss_all(), IGRAPH_IN,
                                     IGRAPH_LOOPS, weights));
        IGRAPH_CHECK(igraph_strength(graph, &out_deg, igraph_vss_all(), IGRAPH_OUT,
                                     IGRAPH_LOOPS, weights));

        scaling_factor = resolution / sw;

        for (igraph_integer_t j = 0; j < no_of_nodes; j++) {
            for (igraph_integer_t i = 0; i < no_of_nodes; i++) {
                MATRIX(*modmat, i, j) -= VECTOR(out_deg)[i] * VECTOR(in_deg)[j] * scaling_factor;
            }
        }
        igraph_vector_destroy(&in_deg);
        igraph_vector_destroy(&out_deg);
        IGRAPH_FINALLY_CLEAN(2);
    } else {
        IGRAPH_VECTOR_INIT_FINALLY(&deg, no_of_nodes);

        IGRAPH_CHECK(igraph_strength(graph, &deg, igraph_vss_all(), IGRAPH_ALL,
                                     IGRAPH_LOOPS, weights));

        scaling_factor = resolution / 2.0 / sw;

        for (igraph_integer_t j = 0; j < no_of_nodes; j++) {
            for (igraph_integer_t i = 0; i < no_of_nodes; i++) {
                MATRIX(*modmat, i, j) -= VECTOR(deg)[i] * VECTOR(deg)[j] * scaling_factor;
            }
        }
        igraph_vector_destroy(&deg);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return IGRAPH_SUCCESS;
}
