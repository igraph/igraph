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

#include "igraph_conversion.h"
#include "igraph_interface.h"
#include "igraph_structural.h"

/**
 * \function igraph_modularity
 * \brief Calculate the modularity of a graph with respect to some clusters or vertex types
 *
 * The modularity of a graph with respect to some clustering of the vertices
 * (or assignment of vertex types)
 * measures how strongly separated the different clusters are from each
 * other compared to a random null model. It is defined as
 *
 * </para><para>
 * <code>Q = 1/(2m) sum_ij (A_ij - gamma * k_i * k_j / (2m)) * d(c_i,c_j)</code>,
 *
 * </para><para>
 * where \c m is the number of edges, <code>A_ij</code> is the adjacency matrix,
 * \c k_i is the degree of vertex \c i, \c c_i is the cluster that vertex \c i belongs to
 * (or its vertex type), <code>d(i,j)=1</code> if <code>i=j</code> and 0 otherwise,
 * and the sum goes over all <code>i, j</code> pairs of vertices.
 *
 * </para><para>
 * The resolution parameter \c gamma allows weighting the random null model, which
 * might be useful when finding partitions with a high modularity. Maximizing modularity
 * with higher values of the resolution parameter typically result in more, smaller clusters
 * when finding partitions with a high modularity. Lower values typically results in
 * fewer, larger clusters. The original definition of modularity is retrieved
 * when setting <code>gamma=1</code>.
 *
 * </para><para>
 * Modularity can also be calculated on directed graphs. This only requires a relatively
 * modest change
 *
 * </para><para>
 * <code>Q = 1/(m) sum_ij (A_ij - gamma * k^out_i * k^in_j / m) * d(c_i,c_j)</code>,
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
 * igraph assumes that the modularity of graphs with no edges is \c NaN; see
 * the discussion in https://github.com/igraph/igraph/issues/1539 for
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
 * For the introduction of the resolution parameter, see Reichardt, J., and
 * Bornholdt, S. (2006). Statistical mechanics of community detection. Physical
 * Review E 74, 016110. https://doi.org/10.1103/PhysRevE.74.016110
 *
 * \param graph      The input graph.
 * \param membership Numeric vector of integer values which gives the type of each
 *                   vertex, i.e. the cluster to which it belongs.
 *                   It does not have to be consecutive, i.e. empty communities
 *                   are allowed.
 * \param weights    Weight vector or \c NULL if no weights are specified.
 * \param resolution Resolution parameter. Must be greater than or equal to 0.
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
                      const igraph_vector_t *membership,
                      const igraph_vector_t *weights,
                      const igraph_real_t resolution,
                      const igraph_bool_t directed,
                      igraph_real_t *modularity) {

    igraph_vector_t e, k_out, k_in;
    igraph_long_t types = (igraph_long_t) igraph_vector_max(membership) + 1;
    igraph_long_t no_of_edges = igraph_ecount(graph);
    igraph_long_t i;
    igraph_long_t from, to;
    igraph_real_t m;
    igraph_long_t c1, c2;
    /* Only consider the graph as directed if it actually is directed */
    igraph_bool_t use_directed = directed && igraph_is_directed(graph);
    igraph_real_t directed_multiplier = (use_directed ? 1 : 2);

    if (igraph_vector_size(membership) != igraph_vcount(graph)) {
        IGRAPH_ERROR("Membership vector size differs from number of vertices.",
                     IGRAPH_EINVAL);
    }
    if (resolution < 0.0) {
      IGRAPH_ERROR("The resolution parameter must be non-negative.", IGRAPH_EINVAL);
    }

    if (no_of_edges == 0) {
        /* Special case: the modularity of graphs with no edges is not
         * well-defined */
        if (modularity) {
            *modularity = IGRAPH_NAN;
        }
        return IGRAPH_SUCCESS;
    }

    /* At this point, the 'membership' vector does not have length zero,
       thus it is safe to call igraph_vector_min(). */
    if (igraph_vector_min(membership) < 0) {
        IGRAPH_ERROR("Invalid membership vector: negative entry.", IGRAPH_EINVAL);
    }

    IGRAPH_VECTOR_INIT_FINALLY(&e, types);
    IGRAPH_VECTOR_INIT_FINALLY(&k_out, types);
    IGRAPH_VECTOR_INIT_FINALLY(&k_in, types);

    if (weights) {
        if (igraph_vector_size(weights) != no_of_edges)
            IGRAPH_ERROR("Vector size differs from number of edges.",
                         IGRAPH_EINVAL);
        m = 0.0;
        for (i = 0; i < no_of_edges; i++) {
            igraph_real_t w = VECTOR(*weights)[i];
            if (w < 0) {
                IGRAPH_ERROR("Negative weight in weight vector.", IGRAPH_EINVAL);
            }
            igraph_edge(graph, (igraph_long_t) i, &from, &to);
            c1 = (igraph_long_t) VECTOR(*membership)[from];
            c2 = (igraph_long_t) VECTOR(*membership)[to];
            if (c1 == c2) {
                VECTOR(e)[c1] += directed_multiplier * w;
            }
            VECTOR(k_out)[c1] += w;
            VECTOR(k_in)[c2]  += w;
            m += w;
        }
    } else {
        m = no_of_edges;
        for (i = 0; i < no_of_edges; i++) {
            igraph_edge(graph, (igraph_long_t) i, &from, &to);
            c1 = (igraph_long_t) VECTOR(*membership)[from];
            c2 = (igraph_long_t) VECTOR(*membership)[to];
            if (c1 == c2) {
                VECTOR(e)[c1] += directed_multiplier;
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
    igraph_vector_scale(&e, 1.0/( directed_multiplier * m ) );

    *modularity = 0.0;
    if (m > 0) {
        for (i = 0; i < types; i++) {
            *modularity += VECTOR(e)[i];
            *modularity -= resolution * VECTOR(k_out)[i] * VECTOR(k_in)[i];
        }
    }

    igraph_vector_destroy(&e);
    igraph_vector_destroy(&k_out);
    igraph_vector_destroy(&k_in);
    IGRAPH_FINALLY_CLEAN(3);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_modularity_matrix
 * \brief Calculate the modularity matrix
 *
 * This function returns the modularity matrix defined as
 *
 * </para><para>
 * <code>B_ij = A_ij - gamma * k_i * k_j / (2m)</code>
 *
 * </para><para>
 * where \c A_ij is the adjacency matrix, \c gamma is the resolution parameter,
 * \c k_i is the degree of vertex \c i, and \c m is the number of edges in the graph.
 *
 * Note that self-loops are multiplied by 2 in this
 * implementation. If weights are specified, the weighted counterparts are used.
 *
 * \param graph      The input graph.
 * \param weights    Edge weights, pointer to a vector. If this is a null pointer
 *                   then every edge is assumed to have a weight of 1.
 * \param resolution Resolution parameter. Must be greater than or equal to 0.
 *                   Default is 1. Lower values favor fewer, larger communities;
 *                   higher values favor more, smaller communities.
 * \param modmat     Pointer to an initialized matrix in which the modularity
 *                   matrix is stored.
 *
 * \sa \ref igraph_modularity()
 */
igraph_error_t igraph_modularity_matrix(const igraph_t *graph,
                             const igraph_vector_t *weights,
                             const igraph_real_t resolution,
                             igraph_matrix_t *modmat) {

    igraph_long_t no_of_nodes = igraph_vcount(graph);
    igraph_long_t no_of_edges = igraph_ecount(graph);
    igraph_real_t sw = weights ? igraph_vector_sum(weights) : no_of_edges;
    igraph_vector_t deg;
    igraph_long_t i, j;

    if (weights && igraph_vector_size(weights) != no_of_edges) {
        IGRAPH_ERROR("Invalid weight vector length", IGRAPH_EINVAL);
    }

    if (resolution < 0.0) {
        IGRAPH_ERROR("The resolution parameter must be non-negative", IGRAPH_EINVAL);
    }

    IGRAPH_VECTOR_INIT_FINALLY(&deg, no_of_nodes);
    if (!weights) {
        IGRAPH_CHECK(igraph_degree(graph, &deg, igraph_vss_all(), IGRAPH_ALL,
                                   IGRAPH_LOOPS));
    } else {
        IGRAPH_CHECK(igraph_strength(graph, &deg, igraph_vss_all(), IGRAPH_ALL,
                                     IGRAPH_LOOPS, weights));
    }
    IGRAPH_CHECK(igraph_get_adjacency(graph, modmat, IGRAPH_GET_ADJACENCY_BOTH,
                                      /*eids=*/ 0));

    for (i = 0; i < no_of_nodes; i++) {
        MATRIX(*modmat, i, i) *= 2;
    }
    for (i = 0; i < no_of_nodes; i++) {
        for (j = 0; j < no_of_nodes; j++) {
            MATRIX(*modmat, i, j) -= resolution * VECTOR(deg)[i] * VECTOR(deg)[j] / 2.0 / sw;
        }
    }

    igraph_vector_destroy(&deg);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}
