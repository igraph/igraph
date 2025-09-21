/*
   igraph library.
   Copyright (C) 2007-2020  The igraph development team <igraph@igraph.org>

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

#include "igraph_centrality.h"

#include "igraph_interface.h"
#include "igraph_structural.h"
#include "igraph_vector.h"

#include "core/math.h"

/**
 * \function igraph_centralization
 * \brief Calculate the centralization score from the node level scores.
 *
 * For a centrality score defined on the vertices of a graph, it is
 * possible to define a graph level centralization index, by
 * calculating the sum of the deviations from the maximum centrality
 * score. Consequently, the higher the centralization index of the
 * graph, the more centralized the structure is.
 *
 * </para><para>
 * In order to make graphs of different sizes comparable,
 * the centralization index is usually normalized to a number between
 * zero and one, by dividing the (unnormalized) centralization score
 * of the most centralized structure with the same number of vertices.
 *
 * </para><para>
 * For most centrality indices, the most centralized structure is the
 * star graph, a single center connected to all other nodes in the network.
 * There is some variation depending on whether the graph is directed or not,
 * whether loop edges are allowed, etc.
 *
 * </para><para>
 * This function simply calculates the graph level index, if the node
 * level scores and the theoretical maximum are given. It is called by
 * all the measure-specific centralization functions. It uses the calculation
 *
 * </para><para>
 * <code>C = sum_v ((max_u c_u) - c_v)</code>
 *
 * </para><para>
 * where \c c are the centrality scores passed in \p scores. If \p normalized
 * is \c true, then <code>C/theoretical_max</code> is returned.
 *
 * \param scores A vector containing the node-level centrality scores.
 * \param theoretical_max The graph level centrality score of the most
 *     centralized graph with the same number of vertices. Only used
 *     if \c normalized set to true.
 * \param normalized Boolean, whether to normalize the centralization
 *     by dividing the supplied theoretical maximum.
 * \return The graph level index.
 *
 * \sa \ref igraph_centralization_degree(), \ref
 * igraph_centralization_betweenness(), \ref
 * igraph_centralization_closeness(), and \ref
 * igraph_centralization_eigenvector_centrality() for specific
 * centralization functions.
 *
 * Time complexity: O(n), the length of the score vector.
 *
 * \example examples/simple/centralization.c
 */

igraph_real_t igraph_centralization(const igraph_vector_t *scores,
                                    igraph_real_t theoretical_max,
                                    igraph_bool_t normalized) {

    igraph_int_t no_of_nodes = igraph_vector_size(scores);
    igraph_real_t cent;

    if (no_of_nodes != 0) {
        igraph_real_t maxscore = igraph_vector_max(scores);
        cent = no_of_nodes * maxscore - igraph_vector_sum(scores);
        if (normalized) {
            cent = cent / theoretical_max;
        }
    } else {
        cent = IGRAPH_NAN;
    }

    return cent;
}

/**
 * \function igraph_centralization_degree
 * \brief Calculate vertex degree and graph centralization.
 *
 * This function calculates the degree of the vertices by passing its
 * arguments to \ref igraph_degree(); and it calculates the graph
 * level centralization index based on the results by calling \ref
 * igraph_centralization().
 *
 * \param graph The input graph.
 * \param res A vector if you need the node-level degree scores, or a
 *     null pointer otherwise.
 * \param mode Constant the specifies the type of degree for directed
 *     graphs. Possible values: \c IGRAPH_IN, \c IGRAPH_OUT and \c
 *     IGRAPH_ALL. This argument is ignored for undirected graphs.
 * \param loops Specifies how to treat loop edges when calculating the
 *     degree (and the centralization). \c IGRAPH_NO_LOOPS ignores loop
 *     edges; \c IGRAPH_LOOPS_ONCE counts each loop edge only once;
 *     \c IGRAPH_LOOPS_TWICE counts each loop edge twice in undirected
 *     graphs and once in directed graphs.
 * \param centralization Pointer to a real number, the centralization
 *     score is placed here.
 * \param theoretical_max Pointer to real number or a null pointer. If
 *     not a null pointer, then the theoretical maximum graph
 *     centrality score for a graph with the same number vertices is
 *     stored here.
 * \param normalized Boolean, whether to calculate a normalized
 *     centralization score. See \ref igraph_centralization() for how
 *     the normalization is done.
 * \return Error code.
 *
 * \sa \ref igraph_centralization(), \ref igraph_degree().
 *
 * Time complexity: the complexity of \ref igraph_degree() plus O(n),
 * the number of vertices queried, for calculating the centralization
 * score.
 */

igraph_error_t igraph_centralization_degree(
    const igraph_t *graph, igraph_vector_t *res, igraph_neimode_t mode,
    igraph_loops_t loops, igraph_real_t *centralization,
    igraph_real_t *theoretical_max, igraph_bool_t normalized
) {

    igraph_vector_t myscores;
    igraph_vector_t *scores = res;
    igraph_real_t *tmax = theoretical_max, mytmax;

    if (!tmax) {
        tmax = &mytmax;
    }

    if (!res) {
        scores = &myscores;
        IGRAPH_VECTOR_INIT_FINALLY(scores, 0);
    }

    IGRAPH_CHECK(igraph_strength(graph, scores, igraph_vss_all(), mode, loops, 0));
    IGRAPH_CHECK(igraph_centralization_degree_tmax(graph, 0, mode, loops, tmax));

    *centralization = igraph_centralization(scores, *tmax, normalized);

    if (!res) {
        igraph_vector_destroy(scores);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_centralization_degree_tmax
 * \brief Theoretical maximum for graph centralization based on degree.
 *
 * This function returns the theoretical maximum graph centrality
 * based on vertex degree.
 *
 * </para><para>
 * There are two ways to call this function, the first is to supply a
 * graph as the \p graph argument, and then the number of
 * vertices is taken from this object, and its directedness is
 * considered as well. The \p nodes argument is ignored in
 * this case. The \p mode argument is also ignored if the
 * supplied graph is undirected.
 *
 * </para><para>
 * The other way is to supply a null pointer as the \p graph
 * argument. In this case the \p nodes and \p mode
 * arguments are considered.
 *
 * </para><para>
 * The most centralized structure is the star. More specifically, for
 * undirected graphs it is the star, for directed graphs it is the
 * in-star or the out-star.
 *
 * \param graph A graph object or a null pointer, see the description
 *     above.
 * \param nodes The number of nodes. This is ignored if the
 *     \p graph argument is not a null pointer.
 * \param mode Constant, whether the calculation is based on in-degree
 *     (\c IGRAPH_IN), out-degree (\c IGRAPH_OUT)
 *     or total degree (\c IGRAPH_ALL). This is ignored if
 *     the \p graph argument is not a null pointer and the
 *     given graph is undirected.
 * \param loops Specifies how to treat loop edges when calculating the
 *     degree (and the centralization). \c IGRAPH_NO_LOOPS ignores loop
 *     edges; \c IGRAPH_LOOPS_ONCE counts each loop edge only once;
 *     \c IGRAPH_LOOPS_TWICE counts each loop edge twice in undirected
 *     graphs and once in directed graphs.
 * \param res Pointer to a real variable, the result is stored here.
 * \return Error code.
 *
 * Time complexity: O(1).
 *
 * \sa \ref igraph_centralization_degree() and \ref
 * igraph_centralization().
 */

igraph_error_t igraph_centralization_degree_tmax(
    const igraph_t *graph, igraph_int_t nodes, igraph_neimode_t mode,
    igraph_loops_t loops, igraph_real_t *res
) {

    igraph_bool_t directed = (mode != IGRAPH_ALL);
    igraph_real_t real_nodes;

    if (graph) {
        directed = igraph_is_directed(graph);
        nodes = igraph_vcount(graph);
    } else {
        if (nodes < 0) {
            IGRAPH_ERROR("Number of vertices must not be negative.", IGRAPH_EINVAL);
        }
    }

    if (nodes == 0) {
        *res = IGRAPH_NAN;
        return IGRAPH_SUCCESS;
    }

    real_nodes = nodes;    /* implicit cast to igraph_real_t */

    if (directed) {
        switch (mode) {
        case IGRAPH_IN:
        case IGRAPH_OUT:
            if (loops == IGRAPH_NO_LOOPS) {
                *res = (real_nodes - 1) * (real_nodes - 1);
            } else {
                *res = (real_nodes - 1) * real_nodes;
            }
            break;
        case IGRAPH_ALL:
            if (loops == IGRAPH_NO_LOOPS) {
                *res = 2 * (real_nodes - 1) * (real_nodes - 2);
            } else {
                *res = 2 * (real_nodes - 1) * (real_nodes - 1);
            }
            break;
        }
    } else {
        if (loops == IGRAPH_NO_LOOPS) {
            *res = (real_nodes - 1) * (real_nodes - 2);
        } else if (loops == IGRAPH_LOOPS_ONCE) {
            *res = (real_nodes - 1) * (real_nodes - 1);
        } else {
            *res = (real_nodes - 1) * real_nodes;
        }
    }

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_centralization_betweenness
 * \brief Calculate vertex betweenness and graph centralization.
 *
 * This function calculates the betweenness centrality of the vertices
 * by passing its arguments to \ref igraph_betweenness(); and it
 * calculates the graph level centralization index based on the
 * results by calling \ref igraph_centralization().
 *
 * \param graph The input graph.
 * \param res A vector if you need the node-level betweenness scores, or a
 *     null pointer otherwise.
 * \param directed Boolean, whether to consider directed paths when
 *     calculating betweenness.
 * \param centralization Pointer to a real number, the centralization
 *     score is placed here.
 * \param theoretical_max Pointer to real number or a null pointer. If
 *     not a null pointer, then the theoretical maximum graph
 *     centrality score for a graph with the same number vertices is
 *     stored here.
 * \param normalized Boolean, whether to calculate a normalized
 *     centralization score. See \ref igraph_centralization() for how
 *     the normalization is done.
 * \return Error code.
 *
 * \sa \ref igraph_centralization(), \ref igraph_betweenness().
 *
 * Time complexity: the complexity of \ref igraph_betweenness() plus
 * O(n), the number of vertices queried, for calculating the
 * centralization score.
 */

igraph_error_t igraph_centralization_betweenness(const igraph_t *graph,
                                      igraph_vector_t *res,
                                      igraph_bool_t directed,
                                      igraph_real_t *centralization,
                                      igraph_real_t *theoretical_max,
                                      igraph_bool_t normalized) {

    igraph_vector_t myscores;
    igraph_vector_t *scores = res;
    igraph_real_t *tmax = theoretical_max, mytmax;

    if (!tmax) {
        tmax = &mytmax;
    }

    if (!res) {
        scores = &myscores;
        IGRAPH_VECTOR_INIT_FINALLY(scores, 0);
    }

    IGRAPH_CHECK(igraph_betweenness(graph, /*weights=*/ 0, scores, igraph_vss_all(), directed, false));

    IGRAPH_CHECK(igraph_centralization_betweenness_tmax(graph, 0, directed, tmax));

    *centralization = igraph_centralization(scores, *tmax, normalized);

    if (!res) {
        igraph_vector_destroy(scores);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_centralization_betweenness_tmax
 * \brief Theoretical maximum for graph centralization based on betweenness.
 *
 * This function returns the theoretical maximum graph centrality
 * based on vertex betweenness.
 *
 * </para><para>
 * There are two ways to call this function, the first is to supply a
 * graph as the \p graph argument, and then the number of
 * vertices is taken from this object, and its directedness is
 * considered as well. The \p nodes argument is ignored in
 * this case. The \p directed argument is also ignored if the
 * supplied graph is undirected.
 *
 * </para><para>
 * The other way is to supply a null pointer as the \p graph
 * argument. In this case the \p nodes and \p directed
 * arguments are considered.
 *
 * </para><para>
 * The most centralized structure is the star.
 *
 * \param graph A graph object or a null pointer, see the description
 *     above.
 * \param nodes The number of nodes. This is ignored if the
 *     \p graph argument is not a null pointer.
 * \param directed Boolean, whether to use directed paths in
 *     the betweenness calculation. This argument is ignored if
 *     \p graph is not a null pointer and it is undirected.
 * \param res Pointer to a real variable, the result is stored here.
 * \return Error code.
 *
 * Time complexity: O(1).
 *
 * \sa \ref igraph_centralization_betweenness() and \ref
 * igraph_centralization().
 */

igraph_error_t igraph_centralization_betweenness_tmax(const igraph_t *graph,
        igraph_int_t nodes,
        igraph_bool_t directed,
        igraph_real_t *res) {
    igraph_real_t real_nodes;

    if (graph) {
        directed = directed && igraph_is_directed(graph);
        nodes = igraph_vcount(graph);
    } else {
        if (nodes < 0) {
            IGRAPH_ERROR("Number of vertices must not be negative.", IGRAPH_EINVAL);
        }
    }

    if (nodes == 0) {
        *res = IGRAPH_NAN;
        return IGRAPH_SUCCESS;
    }

    real_nodes = nodes;    /* implicit cast to igraph_real_t */

    if (directed) {
        *res = (real_nodes - 1) * (real_nodes - 1) * (real_nodes - 2);
    } else {
        *res = (real_nodes - 1) * (real_nodes - 1) * (real_nodes - 2) / 2.0;
    }

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_centralization_closeness
 * \brief Calculate vertex closeness and graph centralization.
 *
 * This function calculates the closeness centrality of the vertices
 * by passing its arguments to \ref igraph_closeness(); and it
 * calculates the graph level centralization index based on the
 * results by calling \ref igraph_centralization().
 *
 * \param graph The input graph.
 * \param res A vector if you need the node-level closeness scores, or a
 *     null pointer otherwise.
 * \param mode Constant the specifies the type of closeness for directed
 *     graphs. Possible values: \c IGRAPH_IN, \c IGRAPH_OUT and \c
 *     IGRAPH_ALL. This argument is ignored for undirected graphs. See
 *     \ref igraph_closeness() argument with the same name for more.
 * \param centralization Pointer to a real number, the centralization
 *     score is placed here.
 * \param theoretical_max Pointer to real number or a null pointer. If
 *     not a null pointer, then the theoretical maximum graph
 *     centrality score for a graph with the same number vertices is
 *     stored here.
 * \param normalized Boolean, whether to calculate a normalized
 *     centralization score. See \ref igraph_centralization() for how
 *     the normalization is done.
 * \return Error code.
 *
 * \sa \ref igraph_centralization(), \ref igraph_closeness().
 *
 * Time complexity: the complexity of \ref igraph_closeness() plus
 * O(n), the number of vertices queried, for calculating the
 * centralization score.
 */

igraph_error_t igraph_centralization_closeness(const igraph_t *graph,
                                    igraph_vector_t *res,
                                    igraph_neimode_t mode,
                                    igraph_real_t *centralization,
                                    igraph_real_t *theoretical_max,
                                    igraph_bool_t normalized) {

    igraph_vector_t myscores;
    igraph_vector_t *scores = res;
    igraph_real_t *tmax = theoretical_max, mytmax;

    if (!tmax) {
        tmax = &mytmax;
    }

    if (!res) {
        scores = &myscores;
        IGRAPH_VECTOR_INIT_FINALLY(scores, 0);
    }

    IGRAPH_CHECK(igraph_closeness(graph, scores, NULL, NULL, igraph_vss_all(), mode,
                                  /*weights=*/ 0, /*normalized=*/ 1));

    IGRAPH_CHECK(igraph_centralization_closeness_tmax(graph, 0, mode,
                 tmax));

    *centralization = igraph_centralization(scores, *tmax, normalized);

    if (!res) {
        igraph_vector_destroy(scores);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_centralization_closeness_tmax
 * \brief Theoretical maximum for graph centralization based on closeness.
 *
 * This function returns the theoretical maximum graph centrality
 * based on vertex closeness.
 *
 * </para><para>
 * There are two ways to call this function, the first is to supply a
 * graph as the \p graph argument, and then the number of
 * vertices is taken from this object, and its directedness is
 * considered as well. The \p nodes argument is ignored in
 * this case. The \p mode argument is also ignored if the
 * supplied graph is undirected.
 *
 * </para><para>
 * The other way is to supply a null pointer as the \p graph
 * argument. In this case the \p nodes and \p mode
 * arguments are considered.
 *
 * </para><para>
 * The most centralized structure is the star.
 *
 * \param graph A graph object or a null pointer, see the description
 *     above.
 * \param nodes The number of nodes. This is ignored if the
 *     \p graph argument is not a null pointer.
 * \param mode Constant, specifies what kind of distances to consider
 *     to calculate closeness. See the \p mode argument of
 *     \ref igraph_closeness() for details. This argument is ignored
 *     if \p graph is not a null pointer and it is undirected.
 * \param res Pointer to a real variable, the result is stored here.
 * \return Error code.
 *
 * Time complexity: O(1).
 *
 * \sa \ref igraph_centralization_closeness() and \ref
 * igraph_centralization().
 */

igraph_error_t igraph_centralization_closeness_tmax(const igraph_t *graph,
        igraph_int_t nodes,
        igraph_neimode_t mode,
        igraph_real_t *res) {
    igraph_real_t real_nodes;

    if (graph) {
        nodes = igraph_vcount(graph);
        if (!igraph_is_directed(graph)) {
            mode = IGRAPH_ALL;
        }
    } else {
        if (nodes < 0) {
            IGRAPH_ERROR("Number of vertices must not be negative.", IGRAPH_EINVAL);
        }
    }

    if (nodes == 0) {
        *res = IGRAPH_NAN;
        return IGRAPH_SUCCESS;
    }

    real_nodes = nodes;    /* implicit cast to igraph_real_t */

    if (mode != IGRAPH_ALL) {
        *res = (real_nodes - 1) * (1.0 - 1.0 / real_nodes);
    } else {
        *res = (real_nodes - 1) * (real_nodes - 2) / (2.0 * real_nodes - 3);
    }

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_centralization_eigenvector_centrality
 * \brief Calculate eigenvector centrality scores and graph centralization.
 *
 * This function calculates the eigenvector centrality of the vertices
 * by passing its arguments to \ref igraph_eigenvector_centrality();
 * and it calculates the graph level centralization index based on the
 * results by calling \ref igraph_centralization().
 *
 * </para><para>
 * Note that vertex-level eigenvector centrality scores do not have
 * a natural scale. As with any eigenvector, their interpretation is
 * invariant to scaling by a constant factor. However, due to how
 * graph-level \em centralization is defined, its value depends on the
 * specific scale/normalization used for vertex-level scores. Which of
 * two graphs will have a higher eigenvector \em centralization depends
 * on the choice of normalization for centralities. This function makes
 * the specific choice of scaling vertex-level centrality scores by their
 * maximum (i.e. it uses the ∞-norm). Other normalization choices, such
 * as the 1-norm or 2-norm are not currently implemented.
 *
 * \param graph The input graph.
 * \param vector A vector if you need the node-level eigenvector
 *      centrality scores, or a null pointer otherwise.
 * \param value If not a null pointer, then the leading eigenvalue is
 *      stored here.
 * \param mode How to consider edge directions in directed graphs.
 *     See \ref igraph_eigenvector_centrality() for details. Ignored
 *     for directed graphs.
 * \param options Options to ARPACK. See \ref igraph_arpack_options_t
 *    for details. Note that the function overwrites the
 *    <code>n</code> (number of vertices) parameter and
 *    it always starts the calculation from a non-random vector
 *    calculated based on the degree of the vertices.
 * \param centralization Pointer to a real number, the centralization
 *     score is placed here.
 * \param theoretical_max Pointer to real number or a null pointer. If
 *     not a null pointer, then the theoretical maximum graph
 *     centrality score for a graph with the same number vertices is
 *     stored here.
 * \param normalized Boolean, whether to calculate a normalized
 *     centralization score. See \ref igraph_centralization() for how
 *     the normalization is done.
 * \return Error code.
 *
 * \sa \ref igraph_centralization(), \ref igraph_eigenvector_centrality().
 *
 * Time complexity: the complexity of \ref
 * igraph_eigenvector_centrality() plus O(|V|), the number of vertices
 * for the calculating the centralization.
 */

igraph_error_t igraph_centralization_eigenvector_centrality(
    const igraph_t *graph,
    igraph_vector_t *vector,
    igraph_real_t *value,
    igraph_neimode_t mode,
    igraph_arpack_options_t *options,
    igraph_real_t *centralization,
    igraph_real_t *theoretical_max,
    igraph_bool_t normalized) {

    igraph_vector_t myscores;
    igraph_vector_t *scores = vector;
    igraph_real_t realvalue, *myvalue = value;
    igraph_real_t *tmax = theoretical_max, mytmax;

    if (!tmax) {
        tmax = &mytmax;
    }

    if (!vector) {
        scores = &myscores;
        IGRAPH_VECTOR_INIT_FINALLY(scores, 0);
    }
    if (!value) {
        myvalue = &realvalue;
    }

    IGRAPH_CHECK(igraph_eigenvector_centrality(graph, scores, myvalue, mode,
                 /*weights=*/ NULL,
                 options));

    IGRAPH_CHECK(igraph_centralization_eigenvector_centrality_tmax(
                     graph, 0, mode, tmax));

    *centralization = igraph_centralization(scores, *tmax, normalized);

    if (!vector) {
        igraph_vector_destroy(scores);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_centralization_eigenvector_centrality_tmax
 * \brief Theoretical maximum centralization for eigenvector centrality.
 *
 * This function returns the theoretical maximum graph centrality
 * based on vertex eigenvector centrality.
 *
 * </para><para>
 * There are two ways to call this function, the first is to supply a
 * graph as the \p graph argument, and then the number of
 * vertices is taken from this object, and its directedness is
 * considered as well. The \p nodes argument is ignored in
 * this case. The \p mode argument is also ignored if the
 * supplied graph is undirected.
 *
 * </para><para>
 * The other way is to supply a null pointer as the \p graph. argument.
 * In this case the \p nodes and \p mode arguments are considered.
 *
 * </para><para>
 * The most centralized directed structure is the in-star with \p mode
 * set to \c IGRAPH_OUT, and the out-star with \p mode set to \c IGRAPH_IN.
 * The most centralized undirected structure is the graph with a single edge.
 *
 * \param graph A graph object or a null pointer, see the description
 *     above.
 * \param nodes The number of nodes. This is ignored if the
 *     \p graph argument is not a null pointer.
 * \param mode How to consider edge directions in directed graphs.
 *     See \ref igraph_eigenvector_centrality() for details. This argument
 *     is ignored if \p graph is not a null pointer and it is
 *     undirected.
 * \param res Pointer to a real variable, the result is stored here.
 * \return Error code.
 *
 * Time complexity: O(1).
 *
 * \sa \ref igraph_centralization_closeness() and \ref
 * igraph_centralization().
 */

igraph_error_t igraph_centralization_eigenvector_centrality_tmax(
    const igraph_t *graph,
    igraph_int_t nodes,
    igraph_neimode_t mode,
    igraph_real_t *res) {

    if (graph) {
        nodes = igraph_vcount(graph);
        if (! igraph_is_directed(graph)) {
            mode = IGRAPH_ALL;
        }
    } else {
        if (nodes < 0) {
            IGRAPH_ERROR("Number of vertices must not be negative.", IGRAPH_EINVAL);
        }
    }

    if (nodes == 0) {
        *res = IGRAPH_NAN;
        return IGRAPH_SUCCESS;
    }

    if (nodes == 1) {
        *res = 0;
        return IGRAPH_SUCCESS;
    }

    if (mode != IGRAPH_ALL) {
        *res = nodes - 1;
    } else {
        *res = nodes - 2;
    }

    return IGRAPH_SUCCESS;
}
