/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
/*
   IGraph library.
   Copyright (C) 2005-2023  The igraph development team <igraph@igraph.org>

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

/**
 * \function igraph_maxdegree
 * \brief The maximum degree in a graph (or set of vertices).
 *
 * The largest in-, out- or total degree of the specified vertices is
 * calculated. If the graph has no vertices, or \p vids is empty,
 * 0 is returned, as this is the smallest possible value for degrees.
 *
 * \param graph The input graph.
 * \param res Pointer to an integer (\c igraph_integer_t), the result
 *        will be stored here.
 * \param vids Vector giving the vertex IDs for which the maximum degree will
 *        be calculated.
 * \param mode Defines the type of the degree.
 *        \c IGRAPH_OUT, out-degree,
 *        \c IGRAPH_IN, in-degree,
 *        \c IGRAPH_ALL, total degree (sum of the
 *        in- and out-degree).
 *        This parameter is ignored for undirected graphs.
 * \param loops Boolean, gives whether the self-loops should be
 *        counted.
 * \return Error code:
 *         \c IGRAPH_EINVVID: invalid vertex ID.
 *         \c IGRAPH_EINVMODE: invalid mode argument.
 *
 * Time complexity: O(v) if \p loops is \c true, and O(v*d) otherwise. v is the number
 * of vertices for which the degree will be calculated, and d is their
 * (average) degree.
 *
 * \sa \ref igraph_degree() to retrieve the degrees for several vertices.
 */
igraph_error_t igraph_maxdegree(const igraph_t *graph, igraph_integer_t *res,
                     igraph_vs_t vids, igraph_neimode_t mode,
                     igraph_bool_t loops) {

    igraph_vector_int_t tmp;

    IGRAPH_VECTOR_INT_INIT_FINALLY(&tmp, 0);

    IGRAPH_CHECK(igraph_degree(graph, &tmp, vids, mode, loops));
    if (igraph_vector_int_size(&tmp) == 0) {
        *res = 0;
    } else {
        *res = igraph_vector_int_max(&tmp);
    }

    igraph_vector_int_destroy(&tmp);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_avg_nearest_neighbor_degree_weighted(const igraph_t *graph,
        igraph_vs_t vids,
        igraph_neimode_t mode,
        igraph_neimode_t neighbor_degree_mode,
        igraph_vector_t *knn,
        igraph_vector_t *knnk,
        const igraph_vector_t *weights) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_vector_int_t neis, edge_neis;
    igraph_integer_t no_vids;
    igraph_vit_t vit;
    igraph_vector_t my_knn_v, *my_knn = knn;
    igraph_vector_t strength;
    igraph_vector_int_t deg;
    igraph_integer_t maxdeg;
    igraph_vector_t deghist;

    if (igraph_vector_size(weights) != igraph_ecount(graph)) {
        IGRAPH_ERROR("Invalid weight vector size.", IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit);
    no_vids = IGRAPH_VIT_SIZE(vit);

    if (!knn) {
        IGRAPH_VECTOR_INIT_FINALLY(&my_knn_v, no_vids);
        my_knn = &my_knn_v;
    } else {
        IGRAPH_CHECK(igraph_vector_resize(knn, no_vids));
    }

    /* Get degree of neighbours */
    IGRAPH_VECTOR_INT_INIT_FINALLY(&deg, no_of_nodes);
    IGRAPH_CHECK(igraph_degree(graph, &deg, igraph_vss_all(),
                               neighbor_degree_mode, IGRAPH_LOOPS));
    IGRAPH_VECTOR_INIT_FINALLY(&strength, no_of_nodes);

    /* Get strength of all nodes */
    IGRAPH_CHECK(igraph_strength(graph, &strength, igraph_vss_all(),
                                 mode, IGRAPH_LOOPS, weights));

    /* Get maximum degree for initialization */
    IGRAPH_CHECK(igraph_maxdegree(graph, &maxdeg, igraph_vss_all(),
                                  mode, IGRAPH_LOOPS));
    IGRAPH_VECTOR_INT_INIT_FINALLY(&neis, maxdeg);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&edge_neis, maxdeg);
    igraph_vector_int_clear(&neis);
    igraph_vector_int_clear(&edge_neis);

    if (knnk) {
        IGRAPH_CHECK(igraph_vector_resize(knnk, maxdeg));
        igraph_vector_null(knnk);
        IGRAPH_VECTOR_INIT_FINALLY(&deghist, maxdeg);
    }

    for (igraph_integer_t i = 0; !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit), i++) {
        igraph_real_t sum = 0.0;
        igraph_integer_t v = IGRAPH_VIT_GET(vit);
        igraph_integer_t nv;
        igraph_real_t str = VECTOR(strength)[v];
        /* Get neighbours and incident edges */
        IGRAPH_CHECK(igraph_neighbors(graph, &neis, v, mode));
        IGRAPH_CHECK(igraph_incident(graph, &edge_neis, v, mode));
        nv = igraph_vector_int_size(&neis);
        for (igraph_integer_t j = 0; j < nv; j++) {
            igraph_integer_t nei = VECTOR(neis)[j];
            igraph_integer_t e = VECTOR(edge_neis)[j];
            igraph_real_t w = VECTOR(*weights)[e];
            sum += w * VECTOR(deg)[nei];
        }
        if (str != 0.0) {
            VECTOR(*my_knn)[i] = sum / str;
        } else {
            VECTOR(*my_knn)[i] = IGRAPH_NAN;
        }
        if (knnk && nv > 0) {
            VECTOR(*knnk)[nv - 1] += sum;
            VECTOR(deghist)[nv - 1] += str;
        }
    }

    igraph_vector_int_destroy(&edge_neis);
    igraph_vector_int_destroy(&neis);
    IGRAPH_FINALLY_CLEAN(2);

    if (knnk) {
        for (igraph_integer_t i = 0; i < maxdeg; i++) {
            igraph_real_t dh = VECTOR(deghist)[i];
            if (dh != 0) {
                VECTOR(*knnk)[i] /= dh;
            } else {
                VECTOR(*knnk)[i] = IGRAPH_NAN;
            }
        }

        igraph_vector_destroy(&deghist);
        IGRAPH_FINALLY_CLEAN(1);
    }

    igraph_vector_destroy(&strength);
    igraph_vector_int_destroy(&deg);
    IGRAPH_FINALLY_CLEAN(2);

    if (!knn) {
        igraph_vector_destroy(&my_knn_v);
        IGRAPH_FINALLY_CLEAN(1);
    }

    igraph_vit_destroy(&vit);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_avg_nearest_neighbor_degree
 * \brief Average neighbor degree.
 *
 * Calculates the average degree of the neighbors for each vertex (\p knn), and
 * optionally, the same quantity as a function of the vertex degree (\p knnk).
 *
 * </para><para>
 * For isolated vertices \p knn is set to NaN. The same is done in \p knnk for
 * vertex degrees that don't appear in the graph.
 *
 * </para><para>
 * The weighted version computes a weighted average of the neighbor degrees as
 *
 * </para><para>
 * <code>k_nn_u = 1/s_u sum_v w_uv k_v</code>,
 *
 * </para><para>
 * where <code>s_u = sum_v w_uv</code> is the sum of the incident edge weights
 * of vertex \c u, i.e. its strength.
 * The sum runs over the neighbors \c v of vertex \c u
 * as indicated by \p mode. <code>w_uv</code> denotes the weighted adjacency matrix
 * and <code>k_v</code> is the neighbors' degree, specified by \p neighbor_degree_mode.
 * This is equation (6) in the reference below.
 *
 * </para><para>
 * When only the <code>k_nn(k)</code> degree correlation function is needed,
 * \ref igraph_degree_correlation_vector() can be used as well. This function provides
 * more flexible control over how degree at each end of directed edges are computed.
 *
 * </para><para>
 * Reference:
 *
 * </para><para>
 * A. Barrat, M. Barthélemy, R. Pastor-Satorras, and A. Vespignani,
 * The architecture of complex weighted networks,
 * Proc. Natl. Acad. Sci. USA 101, 3747 (2004).
 * https://dx.doi.org/10.1073/pnas.0400087101
 *
 * \param graph The input graph. It may be directed.
 * \param vids The vertices for which the calculation is performed.
 * \param mode The type of neighbors to consider in directed graphs.
 *   \c IGRAPH_OUT considers out-neighbors, \c IGRAPH_IN in-neighbors
 *   and \c IGRAPH_ALL ignores edge directions.
 * \param neighbor_degree_mode The type of degree to average in directed graphs.
 *   \c IGRAPH_OUT averages out-degrees, \c IGRAPH_IN averages in-degrees
 *   and \c IGRAPH_ALL ignores edge directions for the degree calculation.
 * \param vids The vertices for which the calculation is performed.
 * \param knn Pointer to an initialized vector, the result will be
 *   stored here. It will be resized as needed. Supply a \c NULL pointer
 *   here if you only want to calculate \c knnk.
 * \param knnk Pointer to an initialized vector, the average
 *   neighbor degree as a function of the vertex degree is stored
 *   here. This is sometimes referred to as the <code>k_nn(k)</code>
 *   degree correlation function. The first (zeroth) element is for degree
 *   one vertices, etc. The calculation is done based only on the vertices
 *   \p vids. Supply a \c NULL pointer here if you don't want to calculate this.
 * \param weights Optional edge weights. Supply a null pointer here
 *   for the non-weighted version.
 *
 * \return Error code.
 *
 * \sa \ref igraph_degree_correlation_vector() for computing only the degree correlation function,
 * with more flexible control over degree computations.
 *
 * Time complexity: O(|V|+|E|), linear in the number of vertices and
 * edges.
 *
 * \example examples/simple/igraph_avg_nearest_neighbor_degree.c
 */
igraph_error_t igraph_avg_nearest_neighbor_degree(const igraph_t *graph,
                                       igraph_vs_t vids,
                                       igraph_neimode_t mode,
                                       igraph_neimode_t neighbor_degree_mode,
                                       igraph_vector_t *knn,
                                       igraph_vector_t *knnk,
                                       const igraph_vector_t *weights) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_vector_int_t neis;
    igraph_integer_t no_vids;
    igraph_vit_t vit;
    igraph_vector_t my_knn_v, *my_knn = knn;
    igraph_vector_int_t deg;
    igraph_integer_t maxdeg;
    igraph_vector_int_t deghist;

    if (weights) {
        return igraph_i_avg_nearest_neighbor_degree_weighted(graph, vids,
                mode, neighbor_degree_mode, knn, knnk, weights);
    }

    IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit);
    no_vids = IGRAPH_VIT_SIZE(vit);

    if (!knn) {
        IGRAPH_VECTOR_INIT_FINALLY(&my_knn_v, no_vids);
        my_knn = &my_knn_v;
    } else {
        IGRAPH_CHECK(igraph_vector_resize(knn, no_vids));
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&deg, no_of_nodes);
    IGRAPH_CHECK(igraph_degree(graph, &deg, igraph_vss_all(),
                               neighbor_degree_mode, IGRAPH_LOOPS));
    IGRAPH_CHECK(igraph_maxdegree(graph, &maxdeg, igraph_vss_all(), mode, IGRAPH_LOOPS));
    IGRAPH_VECTOR_INT_INIT_FINALLY(&neis, maxdeg);
    igraph_vector_int_clear(&neis);

    if (knnk) {
        IGRAPH_CHECK(igraph_vector_resize(knnk, maxdeg));
        igraph_vector_null(knnk);
        IGRAPH_VECTOR_INT_INIT_FINALLY(&deghist, maxdeg);
    }

    for (igraph_integer_t i = 0; !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit), i++) {
        igraph_real_t sum = 0.0;
        igraph_integer_t v = IGRAPH_VIT_GET(vit);
        igraph_integer_t nv;
        IGRAPH_CHECK(igraph_neighbors(graph, &neis, v, mode));
        nv = igraph_vector_int_size(&neis);
        for (igraph_integer_t j = 0; j < nv; j++) {
            igraph_integer_t nei = VECTOR(neis)[j];
            sum += VECTOR(deg)[nei];
        }
        if (nv != 0) {
            VECTOR(*my_knn)[i] = sum / nv;
        } else {
            VECTOR(*my_knn)[i] = IGRAPH_NAN;
        }
        if (knnk && nv > 0) {
            VECTOR(*knnk)[nv - 1] += VECTOR(*my_knn)[i];
            VECTOR(deghist)[nv - 1] += 1;
        }
    }

    if (knnk) {
        for (igraph_integer_t i = 0; i < maxdeg; i++) {
            igraph_integer_t dh = VECTOR(deghist)[i];
            if (dh != 0) {
                VECTOR(*knnk)[i] /= dh;
            } else {
                VECTOR(*knnk)[i] = IGRAPH_NAN;
            }
        }
        igraph_vector_int_destroy(&deghist);
        IGRAPH_FINALLY_CLEAN(1);
    }

    igraph_vector_int_destroy(&neis);
    igraph_vector_int_destroy(&deg);
    igraph_vit_destroy(&vit);
    IGRAPH_FINALLY_CLEAN(3);

    if (!knn) {
        igraph_vector_destroy(&my_knn_v);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_degree_correlation_vector
 * \brief Degree correlation function.
 *
 * \experimental
 *
 * Computes the degree correlation function <code>k_nn(k)</code>, defined as the
 * mean degree of the targets of directed edges whose source has degree \c k.
 * The averaging is done over all directed edges. The \p from_mode and \p to_mode
 * parameters control how the source and target vertex degrees are computed.
 * This way the out-in, out-out, in-in and in-out degree correlation functions
 * can all be computed.
 *
 * </para><para>
 * In undirected graphs, edges are treated as if they were a pair of reciprocal directed
 * ones.
 *
 * </para><para>
 * If P_ij is the joint degree distribution of the graph, computable with
 * \ref igraph_joint_degree_distribution(), then
 * <code>k_nn(k) = (sum_j j P_kj) / (sum_j P_kj)</code>.
 *
 * </para><para>
 * The function \ref igraph_avg_nearest_neighbor_degree(), whose main purpose is to
 * calculate the average neighbor degree for each vertex separately, can also compute
 * <code>k_nn(k)</code>. It differs from this function in that it can take a subset
 * of vertices to base the calculation on, but it does not allow the same fine-grained
 * control over how degrees are computed.
 *
 * </para><para>
 * References:
 *
 * </para><para>
 * R. Pastor-Satorras, A. Vazquez, A. Vespignani:
 * Dynamical and Correlation Properties of the Internet,
 * Phys. Rev. Lett., vol. 87, pp. 258701 (2001).
 * https://doi.org/10.1103/PhysRevLett.87.258701
 *
 * </para><para>
 * A. Vazquez, R. Pastor-Satorras, A. Vespignani:
 * Large-scale topological and dynamical properties of the Internet,
 * Phys. Rev. E, vol. 65, pp. 066130 (2002).
 * https://doi.org/10.1103/PhysRevE.65.066130
 *
 * </para><para>
 * A. Barrat, M. Barthélemy, R. Pastor-Satorras, and A. Vespignani,
 * The architecture of complex weighted networks,
 * Proc. Natl. Acad. Sci. USA 101, 3747 (2004).
 * https://dx.doi.org/10.1073/pnas.0400087101
 *
 * \param graph The input graph.
 * \param weights An optional weight vector. If not \c NULL, weighted averages will be computed.
 * \param knnk An initialized vector, the result will be written here.
 *    <code>knnk[d]</code> will contain the mean degree of vertices connected to
 *    by vertices of degree \c d. Note that in contrast to
 *    \ref igraph_avg_nearest_neighbor_degree(), <code>d=0</code> is also
 *    included.
 * \param from_mode How to compute the degree of sources? Can be \c IGRAPH_OUT
 *    for out-degree, \c IGRAPH_IN for in-degree, or \c IGRAPH_ALL for total degree.
 *    Ignored in undirected graphs.
 * \param to_mode How to compute the degree of sources? Can be \c IGRAPH_OUT
 *    for out-degree, \c IGRAPH_IN for in-degree, or \c IGRAPH_ALL for total degree.
 *    Ignored in undirected graphs.
 * \param directed_neighbors Whether to consider <code>u -> v</code> connections
 *    to be directed. Undirected connections are treated as reciprocal directed ones,
 *    i.e. both <code>u -> v</code> and <code>v -> u</code> will be considered.
 *    Ignored in undirected graphs.
 * \return Error code.
 *
 * \sa \ref igraph_avg_nearest_neighbor_degree() for computing the average neighbour
 * degree of a set of vertices, \ref igraph_joint_degree_distribution() to get the
 * complete joint degree distribution, and \ref igraph_assortativity_degree()
 * to compute the degree assortativity.
 *
 * Time complexity: O(|E| + |V|)
 */
igraph_error_t igraph_degree_correlation_vector(
        const igraph_t *graph, const igraph_vector_t *weights,
        igraph_vector_t *knnk,
        igraph_neimode_t from_mode, igraph_neimode_t to_mode,
        igraph_bool_t directed_neighbors) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_integer_t maxdeg;
    igraph_vector_t weight_sums;
    igraph_vector_int_t *deg_from, *deg_to, deg_out, deg_in, deg_all;

    if (weights && igraph_vector_size(weights) != no_of_edges) {
        IGRAPH_ERRORF("Weight vector length (%" IGRAPH_PRId ") does not match number of edges (%" IGRAPH_PRId ").",
                      IGRAPH_EINVAL,
                      igraph_vector_size(weights), no_of_edges);
    }

    if (! igraph_is_directed(graph)) {
        from_mode = to_mode = IGRAPH_ALL;
        directed_neighbors = false;
    }

    igraph_bool_t have_out = from_mode == IGRAPH_OUT || to_mode == IGRAPH_OUT;
    igraph_bool_t have_in  = from_mode == IGRAPH_IN  || to_mode == IGRAPH_IN;
    igraph_bool_t have_all = from_mode == IGRAPH_ALL || to_mode == IGRAPH_ALL;

    if (have_out) {
        IGRAPH_VECTOR_INT_INIT_FINALLY(&deg_out, no_of_nodes);
        IGRAPH_CHECK(igraph_degree(graph, &deg_out, igraph_vss_all(), IGRAPH_OUT, /* loops */ true));
    }

    if (have_in) {
        IGRAPH_VECTOR_INT_INIT_FINALLY(&deg_in, no_of_nodes);
        IGRAPH_CHECK(igraph_degree(graph, &deg_in, igraph_vss_all(), IGRAPH_IN, /* loops */ true));
    }

    if (have_all) {
        IGRAPH_VECTOR_INT_INIT_FINALLY(&deg_all, no_of_nodes);
        IGRAPH_CHECK(igraph_degree(graph, &deg_all, igraph_vss_all(), IGRAPH_ALL, /* loops */ true));
    }

    switch (from_mode) {
    case IGRAPH_OUT: deg_from = &deg_out; break;
    case IGRAPH_IN:  deg_from = &deg_in;  break;
    case IGRAPH_ALL: deg_from = &deg_all; break;
    default:
        IGRAPH_ERROR("Invalid 'from' mode.", IGRAPH_EINVMODE);
    }

    switch (to_mode) {
    case IGRAPH_OUT: deg_to = &deg_out; break;
    case IGRAPH_IN:  deg_to = &deg_in;  break;
    case IGRAPH_ALL: deg_to = &deg_all; break;
    default:
        IGRAPH_ERROR("Invalid 'to' mode.", IGRAPH_EINVMODE);
    }

    maxdeg = no_of_edges > 0 ? igraph_vector_int_max(deg_from) : 0;

    IGRAPH_VECTOR_INIT_FINALLY(&weight_sums, maxdeg+1);

    IGRAPH_CHECK(igraph_vector_resize(knnk, maxdeg+1));
    igraph_vector_null(knnk);

    for (igraph_integer_t eid=0; eid < no_of_edges; eid++) {
        igraph_integer_t from = IGRAPH_FROM(graph, eid);
        igraph_integer_t to   = IGRAPH_TO(graph, eid);
        igraph_integer_t fromdeg = VECTOR(*deg_from)[from];
        igraph_integer_t todeg   = VECTOR(*deg_to)[to];
        igraph_real_t w = weights ? VECTOR(*weights)[eid] : 1;

        VECTOR(weight_sums)[fromdeg] += w;
        VECTOR(*knnk)[fromdeg] += w * todeg;

        /* Treat undirected edges as reciprocal directed ones */
        if (! directed_neighbors) {
            VECTOR(weight_sums)[todeg] += w;
            VECTOR(*knnk)[todeg] += w * fromdeg;
        }
    }

    IGRAPH_CHECK(igraph_vector_div(knnk, &weight_sums));

    igraph_vector_destroy(&weight_sums);
    IGRAPH_FINALLY_CLEAN(1);

    /* In reverse order of initialization: */

    if (have_all) {
        igraph_vector_int_destroy(&deg_all);
        IGRAPH_FINALLY_CLEAN(1);
    }

    if (have_in) {
        igraph_vector_int_destroy(&deg_in);
        IGRAPH_FINALLY_CLEAN(1);
    }

    if (have_out) {
        igraph_vector_int_destroy(&deg_out);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return IGRAPH_SUCCESS;
}

igraph_error_t igraph_i_strength_all(
        const igraph_t *graph, igraph_vector_t *res,
        igraph_neimode_t mode, igraph_bool_t loops,
        const igraph_vector_t *weights) {

    // When calculating strength for all vertices, iterating over edges is faster
    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);

    IGRAPH_CHECK(igraph_vector_resize(res, no_of_nodes));
    igraph_vector_null(res);

    if (!igraph_is_directed(graph)) {
        mode = IGRAPH_ALL;
    }

    if (loops) {
        if (mode & IGRAPH_OUT) {
            for (igraph_integer_t edge = 0; edge < no_of_edges; ++edge) {
                VECTOR(*res)[IGRAPH_FROM(graph, edge)] += VECTOR(*weights)[edge];
            }
        }
        if (mode & IGRAPH_IN) {
            for (igraph_integer_t edge = 0; edge < no_of_edges; ++edge) {
                VECTOR(*res)[IGRAPH_TO(graph, edge)] += VECTOR(*weights)[edge];
            }
        }
    } else {
        if (mode & IGRAPH_OUT) {
            for (igraph_integer_t edge = 0; edge < no_of_edges; ++edge) {
                igraph_integer_t from = IGRAPH_FROM(graph, edge);
                if (from != IGRAPH_TO(graph, edge)) {
                   VECTOR(*res)[from] += VECTOR(*weights)[edge];
                }
            }
        }
        if (mode & IGRAPH_IN) {
            for (igraph_integer_t edge = 0; edge < no_of_edges; ++edge) {
                igraph_integer_t to = IGRAPH_TO(graph, edge);
                if (IGRAPH_FROM(graph, edge) != to) {
                    VECTOR(*res)[to] += VECTOR(*weights)[edge];
                }
            }
        }
    }

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_strength
 * \brief Strength of the vertices, also called weighted vertex degree.
 *
 * In a weighted network the strength of a vertex is the sum of the
 * weights of all incident edges. In a non-weighted network this is
 * exactly the vertex degree.
 *
 * \param graph The input graph.
 * \param res Pointer to an initialized vector, the result is stored
 *   here. It will be resized as needed.
 * \param vids The vertices for which the calculation is performed.
 * \param mode Gives whether to count only outgoing (\c IGRAPH_OUT),
 *   incoming (\c IGRAPH_IN) edges or both (\c IGRAPH_ALL).
 *   This parameter is ignored for undirected graphs.
 * \param loops A logical scalar, whether to count loop edges as well.
 * \param weights A vector giving the edge weights. If this is a \c NULL
 *   pointer, then \ref igraph_degree() is called to perform the
 *   calculation.
 * \return Error code.
 *
 * Time complexity: O(|V|+|E|), linear in the number vertices and
 * edges.
 *
 * \sa \ref igraph_degree() for the traditional, non-weighted version.
 */
igraph_error_t igraph_strength(const igraph_t *graph, igraph_vector_t *res,
                    const igraph_vs_t vids, igraph_neimode_t mode,
                    igraph_bool_t loops, const igraph_vector_t *weights) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_vit_t vit;
    igraph_integer_t no_vids;
    igraph_vector_int_t degrees;
    igraph_vector_int_t neis;

    if (! weights) {
        IGRAPH_VECTOR_INT_INIT_FINALLY(&degrees, no_of_nodes);
        IGRAPH_CHECK(igraph_vector_resize(res, no_of_nodes));
        IGRAPH_CHECK(igraph_degree(graph, &degrees, vids, mode, loops));
        for (igraph_integer_t i = 0; i < no_of_nodes; i++) {
            VECTOR(*res)[i] = VECTOR(degrees)[i];
        }
        igraph_vector_int_destroy(&degrees);
        IGRAPH_FINALLY_CLEAN(1);
        return IGRAPH_SUCCESS;
    }

    if (igraph_vector_size(weights) != igraph_ecount(graph)) {
        IGRAPH_ERROR("Invalid weight vector length.", IGRAPH_EINVAL);
    }

    if (mode != IGRAPH_OUT && mode != IGRAPH_IN && mode != IGRAPH_ALL) {
        IGRAPH_ERROR("Invalid mode for vertex strength calculation.", IGRAPH_EINVMODE);
    }

    if (igraph_vs_is_all(&vids)) {
        return igraph_i_strength_all(graph, res, mode, loops, weights);
    }

    IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit);
    no_vids = IGRAPH_VIT_SIZE(vit);

    IGRAPH_VECTOR_INT_INIT_FINALLY(&neis, 0);
    IGRAPH_CHECK(igraph_vector_int_reserve(&neis, no_of_nodes));
    IGRAPH_CHECK(igraph_vector_resize(res, no_vids));
    igraph_vector_null(res);

    if (loops) {
        for (igraph_integer_t i = 0; !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit), i++) {
            IGRAPH_CHECK(igraph_incident(graph, &neis, IGRAPH_VIT_GET(vit), mode));
            const igraph_integer_t n = igraph_vector_int_size(&neis);
            for (igraph_integer_t j = 0; j < n; j++) {
                igraph_integer_t edge = VECTOR(neis)[j];
                VECTOR(*res)[i] += VECTOR(*weights)[edge];
            }
        }
    } else {
        for (igraph_integer_t i = 0; !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit), i++) {
            IGRAPH_CHECK(igraph_incident(graph, &neis, IGRAPH_VIT_GET(vit), mode));
            const igraph_integer_t n = igraph_vector_int_size(&neis);
            for (igraph_integer_t j = 0; j < n; j++) {
                igraph_integer_t edge = VECTOR(neis)[j];
                igraph_integer_t from = IGRAPH_FROM(graph, edge);
                igraph_integer_t to = IGRAPH_TO(graph, edge);
                if (from != to) {
                    VECTOR(*res)[i] += VECTOR(*weights)[edge];
                }
            }
        }
    }

    igraph_vit_destroy(&vit);
    igraph_vector_int_destroy(&neis);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}


/**
 * \function igraph_sort_vertex_ids_by_degree
 * \brief Calculate a list of vertex IDs sorted by degree of the corresponding vertex.
 *
 * The list of vertex IDs is returned in a vector that is sorted
 * in ascending or descending order of vertex degree.
 *
 * \param graph The input graph.
 * \param outvids Pointer to an initialized vector that will be
 *        resized and will contain the ordered vertex IDs.
 * \param vids Input vertex selector of vertex IDs to include in
 *        calculation.
 * \param mode Defines the type of the degree.
 *        \c IGRAPH_OUT, out-degree,
 *        \c IGRAPH_IN, in-degree,
 *        \c IGRAPH_ALL, total degree (sum of the
 *        in- and out-degree).
 *        This parameter is ignored for undirected graphs.
 * \param loops Boolean, gives whether the self-loops should be
 *        counted.
 * \param order Specifies whether the ordering should be ascending
 *        (\c IGRAPH_ASCENDING) or descending (\c IGRAPH_DESCENDING).
 * \param only_indices If true, then return a sorted list of indices
 *        into a vector corresponding to \c vids, rather than a list
 *        of vertex IDs. This parameter is ignored if \c vids is set
 *        to all vertices via \ref igraph_vs_all() or \ref igraph_vss_all(),
 *        because in this case the indices and vertex IDs are the
 *        same.
 * \return Error code:
 *         \c IGRAPH_EINVVID: invalid vertex ID.
 *         \c IGRAPH_EINVMODE: invalid mode argument.
 *
 */
igraph_error_t igraph_sort_vertex_ids_by_degree(const igraph_t *graph,
                                     igraph_vector_int_t *outvids,
                                     igraph_vs_t vids,
                                     igraph_neimode_t mode,
                                     igraph_bool_t loops,
                                     igraph_order_t order,
                                     igraph_bool_t only_indices) {
    igraph_integer_t i, n;
    igraph_vector_int_t degrees;
    igraph_vector_int_t vs_vec;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&degrees, 0);
    IGRAPH_CHECK(igraph_degree(graph, &degrees, vids, mode, loops));
    IGRAPH_CHECK(igraph_vector_int_qsort_ind(&degrees, outvids, order));
    if (only_indices || igraph_vs_is_all(&vids) ) {
        igraph_vector_int_destroy(&degrees);
        IGRAPH_FINALLY_CLEAN(1);
    } else {
        IGRAPH_VECTOR_INT_INIT_FINALLY(&vs_vec, 0);
        IGRAPH_CHECK(igraph_vs_as_vector(graph, vids, &vs_vec));
        n = igraph_vector_int_size(outvids);
        for (i = 0; i < n; i++) {
            VECTOR(*outvids)[i] = VECTOR(vs_vec)[VECTOR(*outvids)[i]];
        }
        igraph_vector_int_destroy(&vs_vec);
        igraph_vector_int_destroy(&degrees);
        IGRAPH_FINALLY_CLEAN(2);
    }
    return IGRAPH_SUCCESS;
}
