/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2009-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge, MA 02139 USA

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
 * http://doi.org/10.1103/PhysRevLett.89.208701.
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
 * calculating assortativity for any kind of numeric vertex values.
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
