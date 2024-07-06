/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2009-2023  The igraph development team <igraph@igraph.org>

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

#include "igraph_mixing.h"

#include "igraph_interface.h"
#include "igraph_structural.h"

/**
 * \function igraph_assortativity_nominal
 * \brief Assortativity of a graph based on vertex categories.
 *
 * Assuming the vertices of the input graph belong to different
 * categories, this function calculates the assortativity coefficient of
 * the graph. The assortativity coefficient is between minus one and one
 * and it is one if all connections stay within categories, it is
 * minus one, if the network is perfectly disassortative. For a
 * randomly connected network it is (asymptotically) zero.
 *
 * </para><para>
 * The unnormalized version, computed when \p normalized is set to false,
 * is identical to the modularity, and is defined as follows for
 * directed networks:
 * </para><para>
 * <code>1/m sum_ij (A_ij - k^out_i k^in_j / m) d(i,j)</code>,
 * </para><para>
 * where \c m denotes the number of edges, \c A_ij is the adjacency matrix,
 * <code>k^out</code> and <code>k^in</code> are the out- and in-degrees,
 * and <code>d(i,j)</code> is one if vertices \c i and \c j are in the same
 * category and zero otherwise.
 *
 * </para><para>
 * The normalized assortativity coefficient is obtained by dividing the
 * previous expression by
 * </para><para>
 * <code>1/m sum_ij (m - k^out_i k^in_j d(i,j) / m)</code>.
 * </para><para>
 * It can take any value within the interval [-1, 1].
 *
 * </para><para>
 * Undirected graphs are effectively treated as directed ones with all-reciprocal
 * edges. Thus, self-loops are taken into account twice in undirected graphs.
 *
 * </para><para>
 * References:
 *
 * </para><para>
 * M. E. J. Newman: Mixing patterns in networks,
 * Phys. Rev. E 67, 026126 (2003)
 * https://doi.org/10.1103/PhysRevE.67.026126.
 * See section II and equation (2) for the definition of the concept.
 *
 * </para><para>
 * For an educational overview of assortativity, see
 * M. E. J. Newman,
 * Networks: An Introduction, Oxford University Press (2010).
 * https://doi.org/10.1093/acprof%3Aoso/9780199206650.001.0001.
 *
 * \param graph The input graph, it can be directed or undirected.
 * \param types Integer vector giving the vertex categories. The types
 *    are represented by integers starting at zero.
 * \param res Pointer to a real variable, the result is stored here.
 * \param directed Boolean, it gives whether to consider edge
 *    directions in a directed graph. It is ignored for undirected
 *    graphs.
 * \param normalized Boolean, whether to compute the usual normalized
 *    assortativity. The unnormalized version is identical to
 *    modularity. Supply true here to compute the standard assortativity.
 * \return Error code.
 *
 * Time complexity: O(|E|+t), |E| is the number of edges, t is the
 * number of vertex types.
 *
 * \sa \ref igraph_assortativity() for computing the assortativity
 * based on continuous vertex values instead of discrete categories.
 * \ref igraph_modularity() to compute generalized modularity.
 * \ref igraph_joint_type_distribution() to obtain the mixing matrix.
 *
 * \example examples/simple/igraph_assortativity_nominal.c
 */

igraph_error_t igraph_assortativity_nominal(const igraph_t *graph,
                                            const igraph_vector_int_t *types,
                                            igraph_real_t *res,
                                            igraph_bool_t directed,
                                            igraph_bool_t normalized) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_real_t no_of_edges_real = no_of_edges;   /* for divisions */
    igraph_integer_t no_of_types;
    igraph_vector_int_t ai, bi, eii;
    igraph_real_t sumaibi = 0.0, sumeii = 0.0;

    if (igraph_vector_int_size(types) != no_of_nodes) {
        IGRAPH_ERROR("Invalid types vector length.", IGRAPH_EINVAL);
    }

    if (no_of_nodes == 0) {
        *res = IGRAPH_NAN;
        return IGRAPH_SUCCESS;
    }

    /* 'types' length > 0 here, safe to call vector_min() */
    if (igraph_vector_int_min(types) < 0) {
        IGRAPH_ERROR("Vertex types must not be negative.", IGRAPH_EINVAL);
    }

    directed = directed && igraph_is_directed(graph);

    no_of_types = igraph_vector_int_max(types) + 1;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&ai, no_of_types);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&bi, no_of_types);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&eii, no_of_types);

    for (igraph_integer_t e = 0; e < no_of_edges; e++) {
        igraph_integer_t from = IGRAPH_FROM(graph, e);
        igraph_integer_t to = IGRAPH_TO(graph, e);
        igraph_integer_t from_type = VECTOR(*types)[from];
        igraph_integer_t to_type = VECTOR(*types)[to];

        VECTOR(ai)[from_type] += 1;
        VECTOR(bi)[to_type] += 1;
        if (from_type == to_type) {
            VECTOR(eii)[from_type] += 1;
        }
        if (!directed) {
            if (from_type == to_type) {
                VECTOR(eii)[from_type] += 1;
            }
            VECTOR(ai)[to_type] += 1;
            VECTOR(bi)[from_type] += 1;
        }
    }

    for (igraph_integer_t i = 0; i < no_of_types; i++) {
        sumaibi += (VECTOR(ai)[i] / no_of_edges_real) * (VECTOR(bi)[i] / no_of_edges_real);
        sumeii  += (VECTOR(eii)[i] / no_of_edges_real);
    }

    if (!directed) {
        sumaibi /= 4.0;
        sumeii  /= 2.0;
    }

    if (normalized) {
        *res = (sumeii - sumaibi) / (1.0 - sumaibi);
    } else {
        *res = (sumeii - sumaibi);
    }

    igraph_vector_int_destroy(&eii);
    igraph_vector_int_destroy(&bi);
    igraph_vector_int_destroy(&ai);
    IGRAPH_FINALLY_CLEAN(3);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_assortativity
 * \brief Assortativity based on numeric properties of vertices.
 *
 * This function calculates the assortativity coefficient of a
 * graph based on given values \c x_i for each vertex \c i. This type of
 * assortativity coefficient equals the Pearson correlation of the values
 * at the two ends of the edges.
 *
 * </para><para>
 * The unnormalized covariance of values, computed when \p normalized is
 * set to false, is defined as follows in a directed graph:
 * </para><para>
 * <code>cov(x_out, x_in) = 1/m sum_ij (A_ij - k^out_i k^in_j / m) x_i x_j</code>,
 * </para><para>
 * where \c m denotes the number of edges, \c A_ij is the adjacency matrix, and
 * <code>k^out</code> and <code>k^in</code> are the out- and in-degrees.
 * \c x_out and \c x_in refer to the sets of vertex values at the start and end of
 * the directed edges.
 *
 * </para><para>
 * The normalized covariance, i.e. Pearson correlation, is obtained by dividing
 * the previous expression by
 * <code>sqrt(var(x_out)) sqrt(var(x_in))</code>, where
 * </para><para>
 * <code>var(x_out) = 1/m sum_i k^out_i x_i^2 - (1/m sum_i k^out_i x_i^2)^2</code>
 * </para><para>
 * <code>var(x_in)  = 1/m sum_j k^in_j x_j^2 - (1/m sum_j k^in_j x_j^2)^2</code>
 *
 * </para><para>
 * Undirected graphs are effectively treated as directed graphs where all edges
 * are reciprocal. Therefore, self-loops are effectively considered twice in
 * undirected graphs.
 *
 * </para><para>
 * References:
 *
 * </para><para>
 * M. E. J. Newman: Mixing patterns
 * in networks, Phys. Rev. E 67, 026126 (2003)
 * https://doi.org/10.1103/PhysRevE.67.026126.
 * See section III and equation (21) for the definition, and equation (26) for
 * performing the calculation in directed graphs with the degrees as values.
 *
 * </para><para>
 * M. E. J. Newman: Assortative mixing in networks,
 * Phys. Rev. Lett. 89, 208701 (2002)
 * https://doi.org/10.1103/PhysRevLett.89.208701.
 * See equation (4) for performing the calculation in undirected
 * graphs with the degrees as values.
 *
 * </para><para>
 * For an educational overview of the concept of assortativity, see
 * M. E. J. Newman,
 * Networks: An Introduction, Oxford University Press (2010).
 * https://doi.org/10.1093/acprof%3Aoso/9780199206650.001.0001.
 *
 * \param graph The input graph, it can be directed or undirected.
 * \param values The vertex values, these can be arbitrary numeric
 *     values.
 * \param values_in A second value vector to be used for the incoming
 *     edges when calculating assortativity for a directed graph.
 *     Supply \c NULL here if you want to use the same values
 *     for outgoing and incoming edges. This argument is ignored
 *     (with a warning) if it is not a null pointer and the undirected
 *     assortativity coefficient is being calculated.
 * \param res Pointer to a real variable, the result is stored here.
 * \param directed Boolean, whether to consider edge directions for
 *     directed graphs. It is ignored for undirected graphs.
 * \param normalized Boolean, whether to compute the normalized
 *     covariance, i.e. Pearson correlation. Supply true here to
 *     compute the standard assortativity.
 * \return Error code.
 *
 * Time complexity: O(|E|), linear in the number of edges of the
 * graph.
 *
 * \sa \ref igraph_assortativity_nominal() if you have discrete vertex
 * categories instead of numeric labels, and \ref
 * igraph_assortativity_degree() for the special case of assortativity
 * based on vertex degrees.
 */

igraph_error_t igraph_assortativity(const igraph_t *graph,
                         const igraph_vector_t *values,
                         const igraph_vector_t *values_in,
                         igraph_real_t *res,
                         igraph_bool_t directed,
                         igraph_bool_t normalized) {

    const igraph_integer_t no_of_nodes = igraph_vcount(graph);
    const igraph_integer_t no_of_edges = igraph_ecount(graph);

    directed = directed && igraph_is_directed(graph);

    if (!directed && values_in) {
        IGRAPH_WARNING("Incoming vertex values ignored when calculating undirected assortativity.");
    }

    if (igraph_vector_size(values) != no_of_nodes) {
        IGRAPH_ERROR("Invalid vertex values vector length.", IGRAPH_EINVAL);
    }

    if (values_in && igraph_vector_size(values_in) != no_of_nodes) {
        IGRAPH_ERROR("Invalid incoming vertex values vector length.", IGRAPH_EINVAL);
    }

    if (!directed) {
        igraph_real_t num1 = 0.0, num2 = 0.0, den1 = 0.0;

        for (igraph_integer_t e = 0; e < no_of_edges; e++) {
            igraph_integer_t from = IGRAPH_FROM(graph, e);
            igraph_integer_t to = IGRAPH_TO(graph, e);
            igraph_real_t from_value = VECTOR(*values)[from];
            igraph_real_t to_value = VECTOR(*values)[to];

            num1 += from_value * to_value;
            num2 += from_value + to_value;
            if (normalized) {
                den1 += from_value * from_value + to_value * to_value;
            }
        }

        num1 /= no_of_edges;
        if (normalized) {
            den1 /= no_of_edges * 2.0;
        }
        num2 /= no_of_edges * 2.0;
        num2 = num2 * num2;

        if (normalized) {
            *res = (num1 - num2) / (den1 - num2);
        } else {
            *res = (num1 - num2);
        }

    } else {
        igraph_real_t num1 = 0.0, num2 = 0.0, num3 = 0.0,
                      den1 = 0.0, den2 = 0.0;
        igraph_real_t num, den;

        if (!values_in) {
            values_in = values;
        }

        for (igraph_integer_t e = 0; e < no_of_edges; e++) {
            igraph_integer_t from = IGRAPH_FROM(graph, e);
            igraph_integer_t to = IGRAPH_TO(graph, e);
            igraph_real_t from_value = VECTOR(*values)[from];
            igraph_real_t to_value = VECTOR(*values_in)[to];

            num1 += from_value * to_value;
            num2 += from_value;
            num3 += to_value;
            if (normalized) {
                den1 += from_value * from_value;
                den2 += to_value * to_value;
            }
        }

        num = num1 - num2 * num3 / no_of_edges;
        if (normalized) {
            den = sqrt(den1 - num2 * num2 / no_of_edges) *
                    sqrt(den2 - num3 * num3 / no_of_edges);

            *res = num / den;
        } else {
            *res = num / no_of_edges;
        }
    }

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_assortativity_degree
 * \brief Assortativity of a graph based on vertex degree.
 *
 * Assortativity based on vertex degree, please see the discussion at
 * the documentation of \ref igraph_assortativity() for details.
 * This function simply calls \ref igraph_assortativity() with
 * the degrees as the vertex values and normalization enabled.
 * In the directed case, it uses out-degrees as out-values and
 * in-degrees as in-values.
 *
 * </para><para>
 * For regular graphs, i.e. graphs in which all vertices have the
 * same degree, computing degree correlations is not meaningful,
 * and this function returns NaN.
 *
 * \param graph The input graph, it can be directed or undirected.
 * \param res Pointer to a real variable, the result is stored here.
 * \param directed Boolean, whether to consider edge directions for
 *     directed graphs. This argument is ignored for undirected
 *     graphs. Supply true here to do the natural thing, i.e. use
 *     directed version of the measure for directed graphs and the
 *     undirected version for undirected graphs.
 * \return Error code.
 *
 * Time complexity: O(|E|+|V|), |E| is the number of edges, |V| is
 * the number of vertices.
 *
 * \sa \ref igraph_assortativity() for the general function
 * calculating assortativity for any kind of numeric vertex values,
 * and \ref igraph_joint_degree_distribution() to get the complete
 * joint degree distribution.
 *
 * \example examples/simple/igraph_assortativity_degree.c
 */

igraph_error_t igraph_assortativity_degree(const igraph_t *graph,
                                igraph_real_t *res,
                                igraph_bool_t directed) {

    directed = directed && igraph_is_directed(graph);
    igraph_integer_t no_of_nodes = igraph_vcount(graph);

    /* This function uses igraph_strength() instead of igraph_degree() in order to obtain
     * a vector of reals instead of a vector of integers. */
    if (directed) {
        igraph_vector_t indegree, outdegree;
        IGRAPH_VECTOR_INIT_FINALLY(&indegree, no_of_nodes);
        IGRAPH_VECTOR_INIT_FINALLY(&outdegree, no_of_nodes);
        IGRAPH_CHECK(igraph_strength(graph, &indegree, igraph_vss_all(), IGRAPH_IN, IGRAPH_LOOPS, NULL));
        IGRAPH_CHECK(igraph_strength(graph, &outdegree, igraph_vss_all(), IGRAPH_OUT, IGRAPH_LOOPS, NULL));
        IGRAPH_CHECK(igraph_assortativity(graph, &outdegree, &indegree, res, /* directed */ true, /* normalized */ true));
        igraph_vector_destroy(&indegree);
        igraph_vector_destroy(&outdegree);
        IGRAPH_FINALLY_CLEAN(2);
    } else {
        igraph_vector_t degree;
        IGRAPH_VECTOR_INIT_FINALLY(&degree, no_of_nodes);
        IGRAPH_CHECK(igraph_strength(graph, &degree, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS, NULL));
        IGRAPH_CHECK(igraph_assortativity(graph, &degree, 0, res, /* directed */ false, /* normalized */ true));
        igraph_vector_destroy(&degree);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_joint_degree_matrix
 * \brief The joint degree matrix of a graph.
 *
 * \experimental
 *
 * In graph theory, the joint degree matrix \c J_ij of a graph gives the number
 * of edges, or sum of edge weights, between vertices of degree \c i and degree
 * \c j. This function stores \c J_ij into <code>jdm[i-1, j-1]</code>.
 * Each edge, including self-loops, is counted precisely once, both in undirected
 * and directed graphs.
 *
 * </para><para>
 * <code>sum_(i,j) J_ij</code> is the total number of edges (or total edge weight)
 * \c m in the graph, where <code>(i,j)</code> refers to ordered or unordered
 * pairs in directed and undirected graphs, respectively. Thus <code>J_ij / m</code>
 * is the probability that an edge chosen at random (with probability proportional
 * to its weight) connects vertices with degrees \c i and \c j.
 *
 * </para><para>
 * Note that \c J_ij is similar, but not identical to the joint degree
 * \em distribution, computed by \ref igraph_joint_degree_distribution(),
 * which is defined for \em ordered <code>(i, j)</code> degree
 * pairs even in the undirected case. When considering undirected graphs, the
 * diagonal of the joint degree distribution is twice that of the joint
 * degree matrix.
 *
 * </para><para>
 * References:
 *
 * </para><para>
 * Isabelle Stanton and Ali Pinar:
 * Constructing and sampling graphs with a prescribed joint degree distribution.
 * ACM J. Exp. Algorithmics 17, Article 3.5 (2012).
 * https://doi.org/10.1145/2133803.2330086
 *
 * \param graph A pointer to an initialized graph object.
 * \param weights A vector containing the weights of the edges. If passing a
 *        \c NULL pointer, edges will be assumed to have unit weights, i.e.
 *        the matrix entries will be connection counts.
 * \param jdm A pointer to an initialized matrix that will be resized. The values
 *        will be written here.
 * \param max_out_degree Number of rows in the result, i.e. the largest (out-)degree
 *        to consider. If negative, the largest (out-)degree of the graph will
 *        be used.
 * \param max_in_degree Number of columns in the result, i.e. the largest (in-)degree
 *        to consider. If negative, the largest (in-)degree of the graph will
 *        be used.
 * \return Error code.
 *
 * \sa \ref igraph_joint_degree_distribution() to count ordered vertex pairs instead of
 * edges, or to obtain a normalized matrix.
 *
 * Time complexity: O(E), where E is the number of edges in input graph.
 */

igraph_error_t igraph_joint_degree_matrix(
        const igraph_t *graph, const igraph_vector_t *weights,
        igraph_matrix_t *jdm,
        igraph_integer_t max_out_degree, igraph_integer_t max_in_degree) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_eit_t eit;
    igraph_integer_t eid;
    igraph_integer_t v1id;
    igraph_integer_t v2id;
    igraph_integer_t v1deg;
    igraph_integer_t v2deg;

    if (weights && igraph_vector_size(weights) != no_of_edges) {
        IGRAPH_ERRORF("Weight vector length (%" IGRAPH_PRId ") does not match number of edges (%" IGRAPH_PRId ").",
                      IGRAPH_EINVAL,
                      igraph_vector_size(weights), no_of_edges);
    }

    if (igraph_is_directed(graph)) {
        igraph_vector_int_t out_degrees;
        igraph_vector_int_t in_degrees;

        // Compute max degrees
        IGRAPH_VECTOR_INT_INIT_FINALLY(&out_degrees, no_of_nodes);
        IGRAPH_VECTOR_INT_INIT_FINALLY(&in_degrees, no_of_nodes);
        IGRAPH_CHECK(igraph_degree(graph, &out_degrees, igraph_vss_all(), IGRAPH_OUT, true));
        IGRAPH_CHECK(igraph_degree(graph, &in_degrees, igraph_vss_all(), IGRAPH_IN, true));

        if (max_out_degree < 0) {
            max_out_degree = no_of_nodes > 0 ? igraph_vector_int_max(&out_degrees) : 0;
        }

        if (max_in_degree < 0) {
            max_in_degree = no_of_nodes > 0 ? igraph_vector_int_max(&in_degrees) : 0;
        }

        IGRAPH_CHECK(igraph_matrix_resize(jdm, max_out_degree, max_in_degree));
        igraph_matrix_null(jdm);

        IGRAPH_CHECK(igraph_eit_create(graph, igraph_ess_all(IGRAPH_EDGEORDER_ID), &eit));
        IGRAPH_FINALLY(igraph_eit_destroy, &eit);
        for (; !IGRAPH_EIT_END(eit); IGRAPH_EIT_NEXT(eit)) {
            eid = IGRAPH_EIT_GET(eit);
            v1id = IGRAPH_FROM(graph, eid);
            v2id = IGRAPH_TO(graph, eid);
            v1deg = VECTOR(out_degrees)[v1id];
            v2deg = VECTOR(in_degrees)[v2id];
            if (v1deg <= max_out_degree && v2deg <= max_in_degree) {
                MATRIX(*jdm, v1deg-1, v2deg-1) += weights ? VECTOR(*weights)[eid] : 1;
            }
        }

        igraph_eit_destroy(&eit);
        igraph_vector_int_destroy(&in_degrees);
        igraph_vector_int_destroy(&out_degrees);
        IGRAPH_FINALLY_CLEAN(3);

    } else {
        igraph_vector_int_t degrees;
        igraph_integer_t maxdeg;

        IGRAPH_VECTOR_INT_INIT_FINALLY(&degrees, no_of_nodes);
        IGRAPH_CHECK(igraph_degree(graph, &degrees, igraph_vss_all(), IGRAPH_ALL, true));

        // Compute max degree of the graph only if needed
        if (max_out_degree < 0 || max_in_degree < 0) {
            maxdeg = no_of_nodes > 0 ? igraph_vector_int_max(&degrees) : 0;
        }

        if (max_out_degree < 0) {
            max_out_degree = maxdeg;
        }
        if (max_in_degree < 0) {
            max_in_degree = maxdeg;
        }

        IGRAPH_CHECK(igraph_matrix_resize(jdm, max_out_degree, max_in_degree));
        igraph_matrix_null(jdm);

        IGRAPH_CHECK(igraph_eit_create(graph, igraph_ess_all(IGRAPH_EDGEORDER_ID), &eit));
        IGRAPH_FINALLY(igraph_eit_destroy, &eit);
        while (!IGRAPH_EIT_END(eit)) {
            eid = IGRAPH_EIT_GET(eit);
            v1id = IGRAPH_FROM(graph, eid);
            v2id = IGRAPH_TO(graph, eid);
            v1deg = VECTOR(degrees)[v1id];
            v2deg = VECTOR(degrees)[v2id];

            // Undirected JDMs are symmetrical, needs to be accounted for this when indexing.
            if (v1deg <= max_out_degree && v2deg <= max_in_degree) {
                MATRIX(*jdm, v1deg-1, v2deg-1) += weights ? VECTOR(*weights)[eid] : 1;
            }
            // Do not double-count connections between same-degree vertices.
            if (v1deg != v2deg && v2deg <= max_out_degree && v1deg <= max_in_degree) {
                MATRIX(*jdm, v2deg-1, v1deg-1) += weights ? VECTOR(*weights)[eid] : 1;
            }

            IGRAPH_EIT_NEXT(eit);
        }

        igraph_eit_destroy(&eit);
        igraph_vector_int_destroy(&degrees);
        IGRAPH_FINALLY_CLEAN(2);
    }

    return IGRAPH_SUCCESS;
}

/**
 * Common implementation for igraph_joint_type_distribution() and igraph_joint_degree_distribution()
 *
 * For the max_from/to_type parameters, negative values mean "automatic". These are used
 * only with igraph_joint_degree_distribution().
 *
 * check_types controls whether types should be validated to be non-negative. Validation
 * is only necessary with igraph_joint_type_distribution() but not with igraph_joint_degree_distribution().
 *
 * directed_neighbors must NOT be true when the graph is undirected.
 */
static igraph_error_t mixing_matrix(
        const igraph_t *graph, const igraph_vector_t *weights,
        igraph_matrix_t *p,
        const igraph_vector_int_t *from_types, const igraph_vector_int_t *to_types,
        igraph_bool_t directed_neighbors, igraph_bool_t normalized,
        igraph_integer_t max_from_type, igraph_integer_t max_to_type,
        igraph_bool_t check_types) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_integer_t nrow, ncol;
    igraph_real_t sum;
    igraph_bool_t negative_weight;

    if (igraph_vector_int_size(from_types) != no_of_nodes) {
        IGRAPH_ERROR("Length of 'from' type vector must agree with vertex count.", IGRAPH_EINVAL);
    }

    if (igraph_vector_int_size(to_types) != no_of_nodes) {
        IGRAPH_ERROR("Length of 'to' type vector must agree with vertex count.", IGRAPH_EINVAL);
    }

    if (weights && igraph_vector_size(weights) != no_of_edges) {
        IGRAPH_ERRORF("Weight vector length (%" IGRAPH_PRId ") does not match number of edges (%" IGRAPH_PRId ").",
                      IGRAPH_EINVAL,
                      igraph_vector_size(weights), no_of_edges);
    }

    if (max_from_type < 0) {
        if (no_of_nodes == 0) {
            nrow = 0;
        } else {
            nrow = igraph_vector_int_max(from_types) + 1;
        }
    } else {
        nrow = max_from_type + 1;
    }

    if (max_to_type < 0) {
        if (no_of_nodes == 0) {
            ncol = 0;
        } else if (to_types == from_types) {
            /* Avoid computing the maximum again if target vertex types
             * are the same as source vertex types. */
            ncol = nrow;
        } else {
            ncol = igraph_vector_int_max(to_types) + 1;
        }
    } else {
        ncol = max_to_type + 1;
    }

    if (check_types && no_of_nodes > 0) {
        igraph_integer_t min;

        min = igraph_vector_int_min(from_types);
        if (min < 0) {
            IGRAPH_ERROR("Invalid source vertex type.", IGRAPH_EINVAL);
        }

        if (to_types != from_types) {
            min = igraph_vector_int_min(from_types);
            if (min < 0) {
                IGRAPH_ERROR("Invalid target vertex type.", IGRAPH_EINVAL);
            }
        }
    }

    IGRAPH_CHECK(igraph_matrix_resize(p, nrow, ncol));
    igraph_matrix_null(p);

    sum = 0;
    negative_weight = false;
    for (igraph_integer_t eid=0; eid < no_of_edges; eid++) {
        igraph_integer_t from = IGRAPH_FROM(graph, eid);
        igraph_integer_t to   = IGRAPH_TO(graph, eid);
        igraph_integer_t from_type = VECTOR(*from_types)[from];
        igraph_integer_t to_type   = VECTOR(*to_types)[to];
        igraph_real_t w = weights ? VECTOR(*weights)[eid] : 1;

        if (from_type >= nrow || to_type >= ncol) {
            continue;
        }

        MATRIX(*p, from_type, to_type) += w;
        sum += w;

        if (! directed_neighbors) {
            MATRIX(*p, to_type, from_type) += w;
            sum += w;
        }

        if (w < 0) {
            negative_weight = true;
        }
    }

    if (normalized) {
        if (negative_weight) {
            /* When some edge weights are negative, they cannot be interpreted as sampling weights,
             * and the sum of weights may be zero, potentially leading to Inf/NaN results. */
            IGRAPH_WARNING("Negative edge weights are present. Normalization may not be meaningful.");
        }
        if (no_of_edges > 0) {
            /* Scale only when there are some edges, thus 'sum' can be non-zero. */
            igraph_matrix_scale(p, 1.0 / sum);
        }
    }

    return IGRAPH_SUCCESS;
}


/**
 * \function igraph_joint_degree_distribution
 * \brief The joint degree distribution of a graph.
 *
 * \experimental
 *
 * Computes the joint degree distribution \c P_ij of a graph, used in the
 * study of degree correlations. \c P_ij is the probability that a randomly
 * chosen ordered pair of \em connected vertices have degrees \c i and \c j.
 *
 * </para><para>
 * In directed graphs, directionally connected <code>u -> v</code> pairs
 * are considered. The joint degree distribution of an undirected graph is the
 * same as that of the corresponding directed graph in which all connection are
 * bidirectional, assuming that \p from_mode is \c IGRAPH_OUT, \p to_mode is
 * \c IGRAPH_IN and \p directed_neighbors is true.
 *
 * </para><para>
 * When \p normalized is false, <code>sum_ij P_ij</code> gives the total
 * number of connections in a directed graph, or twice that value in an
 * undirected graph. The sum is taken over ordered <code>(i,j)</code> degree
 * pairs.
 *
 * </para><para>
 * The joint degree distribution relates to other concepts used in the study of
 * degree correlations. If \c P_ij is normalized then the degree correlation
 * function <code>k_nn(k)</code> is obtained as
 *
 * </para><para>
 * <code>k_nn(k) = (sum_j j P_kj) / (sum_j P_kj)</code>.
 *
 * </para><para>
 * The non-normalized degree assortativity is obtained as
 *
 * </para><para>
 * <code>a = sum_ij i j (P_ij - q_i r_j)</code>,
 *
 * </para><para>
 * where <code>q_i = sum_k P_ik</code> and <code>r_j = sum_k P_kj</code>.
 *
 * </para><para>
 * Note that the joint degree distribution \c P_ij is similar, but not identical
 * to the joint degree matrix \c J_ij computed by \ref igraph_joint_degree_matrix().
 * If the graph is undirected, then the diagonal entries of an unnormalized \c P_ij
 * are double that of \c J_ij, as any undirected connection between same-degree vertices
 * is counted in both directions. In contrast to \ref igraph_joint_degree_matrix(),
 * this function returns matrices which include the row and column corresponding
 * to zero degrees. In directed graphs, this row and column is not necessarily
 * zero when \p from_mode is different from \c IGRAPH_OUT or \p to_mode is different
 * from \c IGRAPH_IN.
 *
 * </para><para>
 * References:
 *
 * </para><para>
 * M. E. J. Newman: Mixing patterns in networks,
 * Phys. Rev. E 67, 026126 (2003)
 * https://doi.org/10.1103/PhysRevE.67.026126.
 *
 * \param graph A pointer to an initialized graph object.
 * \param weights A vector containing the weights of the edges. If passing a
 *    \c NULL pointer, edges will be assumed to have unit weights.
 * \param p A pointer to an initialized matrix that will be resized. The \c P_ij
 *    value will be written into <code>p[i,j]</code>.
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
 * \param normalized Whether to normalize the matrix so that entries sum to 1.0.
 *    If false, matrix entries will be connection counts. Normalization is not
 *    meaningful if some edge weights are negative.
 * \param max_from_degree The largest source vertex degree to consider. If negative,
 *    the largest source degree will be used. The row count of the result matrix
 *    is one larger than this value.
 * \param max_to_degree The largest target vertex degree to consider. If negative,
 *    the largest target degree will be used. The column count of the result matrix
 *    is one larger than this value.
 * \return Error code.
 *
 * \sa \ref igraph_joint_degree_matrix() for computing the joint degree matrix;
 * \ref igraph_assortativity_degree() and \ref igraph_assortativity() for
 * degree correlations coefficients, and \ref igraph_degree_correlation_vector()
 * for the degree correlation function.
 *
 * Time complexity: O(E), where E is the number of edges in the input graph.
 */
igraph_error_t igraph_joint_degree_distribution(
        const igraph_t *graph, const igraph_vector_t *weights,
        igraph_matrix_t *p,
        igraph_neimode_t from_mode, igraph_neimode_t to_mode,
        igraph_bool_t directed_neighbors,
        igraph_bool_t normalized,
        igraph_integer_t max_from_degree, igraph_integer_t max_to_degree) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_vector_int_t *deg_from, *deg_to, deg_out, deg_in, deg_all;

    /* Make sure directionality parameters are consistent for undirected graphs. */
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
        IGRAPH_ERROR("Invalid 'from' degree mode.", IGRAPH_EINVMODE);
    }

    switch (to_mode) {
    case IGRAPH_OUT: deg_to = &deg_out; break;
    case IGRAPH_IN:  deg_to = &deg_in;  break;
    case IGRAPH_ALL: deg_to = &deg_all; break;
    default:
        IGRAPH_ERROR("Invalid 'to' degree mode.", IGRAPH_EINVMODE);
    }

    IGRAPH_CHECK(mixing_matrix(graph,
                               weights, p,
                               deg_from, deg_to,
                               directed_neighbors, normalized,
                               max_from_degree, max_to_degree,
                               /*check_types=*/ false));

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


/**
 * \function igraph_joint_type_distribution
 * \brief Mixing matrix for vertex categories.
 *
 * \experimental
 *
 * Computes the mixing matrix M_ij, i.e. the joint distribution of vertex types
 * at the endpoints directed of edges. Categories are represented by non-negative integer
 * indices, passed in \p from_types and \p to_types. The row and column counts of \p m
 * will be one larger than the largest source and target type, respectively. Re-index type
 * vectors using \ref igraph_reindex_membership() if they are not contiguous integers,
 * to avoid producing a very large matrix.
 *
 * </para><para>
 * M_ij is proportional to the probability that a randomly chosen ordered pair of vertices
 * have types \c i and \c j.
 *
 * </para><para>
 * When there is a single categorization of vertices, i.e. \p from_types and \p to_types
 * are the same, M_ij is related to the modularity (\ref igraph_modularity()) and nominal
 * assortativity (\ref igraph_assortativity_nominal()). Let <code>a_i = sum_j M_ij</code> and
 * <code>b_j = sum_i M_ij</code>. If M_ij is normalized, i.e. <code>sum_ij M_ij = 1</code>,
 * and the types represent membership in vertex partitions, then the modularity of the
 * partitioning can be computed as
 *
 * </para><para>
 * <code>Q = sum_ii M_ii - sum_i a_i b_i</code>
 *
 * </para><para>
 * The normalized nominal assortativity is
 *
 * </para><para>
 * <code>Q / (1 - sum_i a_i b_i)</code>
 *
 * </para><para>
 * \ref igraph_joint_degree_distribution() is a special case of this function, with
 * categories consisting vertices of the same degree.
 *
 * </para><para>
 * References:
 *
 * </para><para>
 * M. E. J. Newman: Mixing patterns in networks,
 * Phys. Rev. E 67, 026126 (2003)
 * https://doi.org/10.1103/PhysRevE.67.026126.
 *
 * \param graph The input graph.
 * \param p The mixing matrix M_ij will be stored here.
 * \param weights A vector containing the weights of the edges. If passing a
 *    \c NULL pointer, edges will be assumed to have unit weights.
 * \param from_types Vertex types for source vertices. These must be non-negative integers.
 * \param to_types Vertex types for target vertices. These must be non-negative integers.
 *    If \c NULL, it is assumed to be the same as \p from_types.
 * \param directed Whether to treat edges are directed. Ignored for undirected graphs.
 * \param normalized Whether to normalize the matrix so that entries sum to 1.0.
 *    If false, matrix entries will be connection counts. Normalization is not
 *    meaningful if some edge weights are negative.
 * \return Error code.
 *
 * \sa \ref igraph_joint_degree_distribution() to compute the joint distribution
 * of vertex degrees; \ref igraph_modularity() to compute the modularity of
 * a vertex partitioning; \ref igraph_assortativity_nominal() to compute
 * assortativity based on vertex categories.
 *
 * Time complexity: O(E), where E is the number of edges in the input graph.
 */
igraph_error_t igraph_joint_type_distribution(
        const igraph_t *graph, const igraph_vector_t *weights,
        igraph_matrix_t *p,
        const igraph_vector_int_t *from_types, const igraph_vector_int_t *to_types,
        igraph_bool_t directed, igraph_bool_t normalized) {

    IGRAPH_ASSERT(from_types != NULL);
    if (to_types == NULL) {
        to_types = from_types;
    }
    if (! igraph_is_directed(graph)) {
        directed = false;
    }
    return mixing_matrix(graph, weights, p, from_types, to_types, directed, normalized, -1, -1, true);
}
