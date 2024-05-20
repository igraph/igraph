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

#include "igraph_cocitation.h"

#include "igraph_adjlist.h"
#include "igraph_interface.h"

#include "core/interruption.h"

#include <math.h>

static igraph_error_t igraph_i_cocitation_real(const igraph_t *graph, igraph_matrix_t *res,
                           igraph_vs_t vids, igraph_neimode_t mode,
                           igraph_vector_t *weights);

/**
 * \ingroup structural
 * \function igraph_cocitation
 * \brief Cocitation coupling.
 *
 * Two vertices are cocited if there is another vertex citing both of
 * them. \ref igraph_cocitation() simply counts how many times two vertices are
 * cocited.
 * The cocitation score for each given vertex and all other vertices
 * in the graph will be calculated.
 *
 * \param graph The graph object to analyze.
 * \param res Pointer to a matrix, the result of the calculation will
 *        be stored here. The number of its rows is the same as the
 *        number of vertex IDs in \p vids, the number of
 *        columns is the number of vertices in the graph.
 * \param vids The vertex IDs of the vertices for which the
 *        calculation will be done.
 * \return Error code:
 *         \c IGRAPH_EINVVID: invalid vertex ID.
 *
 * Time complexity: O(|V|d^2), |V| is
 * the number of vertices in the graph,
 * d is the (maximum) degree of
 * the vertices in the graph.
 *
 * \sa \ref igraph_bibcoupling()
 *
 * \example examples/simple/igraph_cocitation.c
 */

igraph_error_t igraph_cocitation(const igraph_t *graph, igraph_matrix_t *res,
                      const igraph_vs_t vids) {
    return igraph_i_cocitation_real(graph, res, vids, IGRAPH_OUT, NULL);
}

/**
 * \ingroup structural
 * \function igraph_bibcoupling
 * \brief Bibliographic coupling.
 *
 * The bibliographic coupling of two vertices is the number
 * of other vertices they both cite, \ref igraph_bibcoupling() calculates
 * this.
 * The bibliographic coupling  score for each given vertex and all
 * other vertices in the graph will be calculated.
 *
 * \param graph The graph object to analyze.
 * \param res Pointer to a matrix, the result of the calculation will
 *        be stored here. The number of its rows is the same as the
 *        number of vertex IDs in \p vids, the number of
 *        columns is the number of vertices in the graph.
 * \param vids The vertex IDs of the vertices for which the
 *        calculation will be done.
 * \return Error code:
 *         \c IGRAPH_EINVVID: invalid vertex ID.
 *
 * Time complexity: O(|V|d^2),
 * |V| is the number of vertices in
 * the graph, d is the (maximum)
 * degree of the vertices in the graph.
 *
 * \sa \ref igraph_cocitation()
 *
 * \example examples/simple/igraph_cocitation.c
 */

igraph_error_t igraph_bibcoupling(const igraph_t *graph, igraph_matrix_t *res,
                       const igraph_vs_t vids) {
    return igraph_i_cocitation_real(graph, res, vids, IGRAPH_IN, NULL);
}

/**
 * \ingroup structural
 * \function igraph_similarity_inverse_log_weighted
 * \brief Vertex similarity based on the inverse logarithm of vertex degrees.
 *
 * The inverse log-weighted similarity of two vertices is the number of
 * their common neighbors, weighted by the inverse logarithm of their degrees.
 * It is based on the assumption that two vertices should be considered
 * more similar if they share a low-degree common neighbor, since high-degree
 * common neighbors are more likely to appear even by pure chance.
 *
 * </para><para>
 * Isolated vertices will have zero similarity to any other vertex.
 * Self-similarities are not calculated.
 *
 * </para><para>
 * Note that the presence of loop edges may yield counter-intuitive
 * results. A node with a loop edge is considered to be a neighbor of itself
 * \em twice (because there are two edge stems incident on the node). Adding a
 * loop edge to a node may decrease its similarity to other nodes, but it may
 * also \em increase it. For instance, if nodes A and B are connected but share
 * no common neighbors, their similarity is zero. However, if a loop edge is
 * added to B, then B itself becomes a common neighbor of A and B and thus the
 * similarity of A and B will be increased. Consider removing loop edges
 * explicitly before invoking this function using \ref igraph_simplify().
 *
 * </para><para>
 * See the following paper for more details: Lada A. Adamic and Eytan Adar:
 * Friends and neighbors on the Web. Social Networks, 25(3):211-230, 2003.
 * https://doi.org/10.1016/S0378-8733(03)00009-1
 *
 * \param graph The graph object to analyze.
 * \param res Pointer to a matrix, the result of the calculation will
 *        be stored here. The number of its rows is the same as the
 *        number of vertex IDs in \p vids, the number of
 *        columns is the number of vertices in the graph.
 * \param vids The vertex IDs of the vertices for which the
 *        calculation will be done.
 * \param mode The type of neighbors to be used for the calculation in
 *        directed graphs. Possible values:
 *        \clist
 *        \cli IGRAPH_OUT
 *          the outgoing edges will be considered for each node. Nodes
 *          will be weighted according to their in-degree.
 *        \cli IGRAPH_IN
 *          the incoming edges will be considered for each node. Nodes
 *          will be weighted according to their out-degree.
 *        \cli IGRAPH_ALL
 *          the directed graph is considered as an undirected one for the
 *          computation. Every node is weighted according to its undirected
 *          degree.
 *        \endclist
 * \return Error code:
 *         \c IGRAPH_EINVVID: invalid vertex ID.
 *
 * Time complexity: O(|V|d^2),
 * |V| is the number of vertices in
 * the graph, d is the (maximum)
 * degree of the vertices in the graph.
 *
 * \example examples/simple/igraph_similarity.c
 */

igraph_error_t igraph_similarity_inverse_log_weighted(const igraph_t *graph,
        igraph_matrix_t *res, const igraph_vs_t vids, igraph_neimode_t mode) {
    igraph_vector_t weights;
    igraph_vector_int_t degrees;
    igraph_neimode_t mode0 = IGRAPH_REVERSE_MODE(mode);
    igraph_integer_t no_of_nodes = igraph_vcount(graph);

    if (mode != IGRAPH_OUT && mode != IGRAPH_IN && mode != IGRAPH_ALL) {
        IGRAPH_ERROR("Invalid mode for inverse log weighted similarity.", IGRAPH_EINVMODE);
    }

    IGRAPH_VECTOR_INIT_FINALLY(&weights, no_of_nodes);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&degrees, no_of_nodes);
    IGRAPH_CHECK(igraph_degree(graph, &degrees, igraph_vss_all(), mode0, true));
    for (igraph_integer_t i = 0; i < no_of_nodes; i++) {
        VECTOR(weights)[i] = VECTOR(degrees)[i];
        if (VECTOR(weights)[i] > 1) {
            VECTOR(weights)[i] = 1.0 / log(VECTOR(weights)[i]);
        }
    }

    IGRAPH_CHECK(igraph_i_cocitation_real(graph, res, vids, mode0, &weights));
    igraph_vector_int_destroy(&degrees);
    igraph_vector_destroy(&weights);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_cocitation_real(const igraph_t *graph, igraph_matrix_t *res,
                           igraph_vs_t vids,
                           igraph_neimode_t mode,
                           igraph_vector_t *weights) {

    const igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_vids;
    igraph_integer_t i;
    igraph_vector_int_t neis;
    igraph_vector_int_t vid_reverse_index;
    igraph_vit_t vit;

    IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit);

    no_of_vids = IGRAPH_VIT_SIZE(vit);

    /* Create a mapping from vertex IDs to the row of the matrix where
     * the result for this vertex will appear */
    IGRAPH_VECTOR_INT_INIT_FINALLY(&vid_reverse_index, no_of_nodes);
    igraph_vector_int_fill(&vid_reverse_index, -1);
    for (IGRAPH_VIT_RESET(vit), i = 0; !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit), i++) {
        igraph_integer_t v = IGRAPH_VIT_GET(vit);
        if (v < 0 || v >= no_of_nodes) {
            IGRAPH_ERROR("Invalid vertex ID in vertex selector.", IGRAPH_EINVVID);
        }
        VECTOR(vid_reverse_index)[v] = i;
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&neis, 0);
    IGRAPH_CHECK(igraph_matrix_resize(res, no_of_vids, no_of_nodes));
    igraph_matrix_null(res);

    /* The result */

    for (igraph_integer_t from = 0; from < no_of_nodes; from++) {
        IGRAPH_ALLOW_INTERRUPTION();

        const igraph_real_t weight = weights ? VECTOR(*weights)[from] : 1;

        IGRAPH_CHECK(igraph_neighbors(graph, &neis, from, mode));
        const igraph_integer_t nei_count = igraph_vector_int_size(&neis);

        for (i = 0; i < nei_count - 1; i++) {
            igraph_integer_t u = VECTOR(neis)[i];
            igraph_integer_t k = VECTOR(vid_reverse_index)[u];
            for (igraph_integer_t j = i + 1; j < nei_count; j++) {
                igraph_integer_t v = VECTOR(neis)[j];
                igraph_integer_t l = VECTOR(vid_reverse_index)[v];
                if (k != -1) {
                    MATRIX(*res, k, v) += weight;
                }
                if (l != -1) {
                    MATRIX(*res, l, u) += weight;
                }
            }
        }
    }

    /* Clean up */
    igraph_vector_int_destroy(&neis);
    igraph_vector_int_destroy(&vid_reverse_index);
    igraph_vit_destroy(&vit);
    IGRAPH_FINALLY_CLEAN(3);

    return IGRAPH_SUCCESS;
}


static igraph_error_t igraph_i_neisets_intersect(
    const igraph_vector_int_t *v1, const igraph_vector_int_t *v2,
    igraph_integer_t *len_union, igraph_integer_t *len_intersection
) {
    /* ASSERT: v1 and v2 are sorted */
    igraph_integer_t i, j, i0, jj0;
    i0 = igraph_vector_int_size(v1); jj0 = igraph_vector_int_size(v2);
    *len_union = i0 + jj0; *len_intersection = 0;
    i = 0; j = 0;
    while (i < i0 && j < jj0) {
        if (VECTOR(*v1)[i] == VECTOR(*v2)[j]) {
            (*len_intersection)++; (*len_union)--;
            i++; j++;
        } else if (VECTOR(*v1)[i] < VECTOR(*v2)[j]) {
            i++;
        } else {
            j++;
        }
    }
    return IGRAPH_SUCCESS;
}

/**
 * \ingroup structural
 * \function igraph_similarity_jaccard
 * \brief Jaccard similarity coefficient for the given vertices.
 *
 * The Jaccard similarity coefficient of two vertices is the number of common
 * neighbors divided by the number of vertices that are neighbors of at
 * least one of the two vertices being considered. This function calculates
 * the pairwise Jaccard similarities for some (or all) of the vertices.
 *
 * \param graph The graph object to analyze
 * \param res Pointer to a matrix, the result of the calculation will
 *        be stored here. The number of its rows and columns is the same
 *        as the number of vertex IDs in \p vids.
 * \param vids The vertex IDs of the vertices for which the
 *        calculation will be done.
 * \param mode The type of neighbors to be used for the calculation in
 *        directed graphs. Possible values:
 *        \clist
 *        \cli IGRAPH_OUT
 *          the outgoing edges will be considered for each node.
 *        \cli IGRAPH_IN
 *          the incoming edges will be considered for each node.
 *        \cli IGRAPH_ALL
 *          the directed graph is considered as an undirected one for the
 *          computation.
 *        \endclist
 * \param loops Whether to include the vertices themselves in the neighbor
 *        sets.
 * \return Error code:
 *        \clist
 *        \cli IGRAPH_ENOMEM
 *           not enough memory for temporary data.
 *        \cli IGRAPH_EINVVID
 *           invalid vertex ID passed.
 *        \cli IGRAPH_EINVMODE
 *           invalid mode argument.
 *        \endclist
 *
 * Time complexity: O(|V|^2 d),
 * |V| is the number of vertices in the vertex iterator given, d is the
 * (maximum) degree of the vertices in the graph.
 *
 * \sa \ref igraph_similarity_dice(), a measure very similar to the Jaccard
 *   coefficient
 *
 * \example examples/simple/igraph_similarity.c
 */
igraph_error_t igraph_similarity_jaccard(const igraph_t *graph, igraph_matrix_t *res,
                              const igraph_vs_t vids, igraph_neimode_t mode, igraph_bool_t loops) {
    igraph_lazy_adjlist_t al;
    igraph_vit_t vit, vit2;
    igraph_integer_t i, j;
    igraph_integer_t len_union, len_intersection;
    igraph_vector_int_t *v1, *v2;
    igraph_integer_t k;

    IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit);
    IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit2));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit2);

    IGRAPH_CHECK(igraph_lazy_adjlist_init(graph, &al, mode, IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE));
    IGRAPH_FINALLY(igraph_lazy_adjlist_destroy, &al);

    IGRAPH_CHECK(igraph_matrix_resize(res, IGRAPH_VIT_SIZE(vit), IGRAPH_VIT_SIZE(vit)));

    if (loops) {
        for (IGRAPH_VIT_RESET(vit); !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit)) {
            i = IGRAPH_VIT_GET(vit);
            v1 = igraph_lazy_adjlist_get(&al, i);
            IGRAPH_CHECK_OOM(v1, "Failed to query neighbors.");
            if (!igraph_vector_int_binsearch(v1, i, &k)) {
                IGRAPH_CHECK(igraph_vector_int_insert(v1, k, i));
            }
        }
    }

    for (IGRAPH_VIT_RESET(vit), i = 0;
         !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit), i++) {
        MATRIX(*res, i, i) = 1.0;
        for (IGRAPH_VIT_RESET(vit2), j = 0;
             !IGRAPH_VIT_END(vit2); IGRAPH_VIT_NEXT(vit2), j++) {
            if (j <= i) {
                continue;
            }

            v1 = igraph_lazy_adjlist_get(&al, IGRAPH_VIT_GET(vit));
            IGRAPH_CHECK_OOM(v1, "Failed to query neighbors.");
            v2 = igraph_lazy_adjlist_get(&al, IGRAPH_VIT_GET(vit2));
            IGRAPH_CHECK_OOM(v2, "Failed to query neighbors.");

            IGRAPH_CHECK(igraph_i_neisets_intersect(v1, v2, &len_union, &len_intersection));
            if (len_union > 0) {
                MATRIX(*res, i, j) = ((igraph_real_t)len_intersection) / len_union;
            } else {
                MATRIX(*res, i, j) = 0.0;
            }
            MATRIX(*res, j, i) = MATRIX(*res, i, j);
        }
    }

    igraph_lazy_adjlist_destroy(&al);
    igraph_vit_destroy(&vit);
    igraph_vit_destroy(&vit2);
    IGRAPH_FINALLY_CLEAN(3);

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup structural
 * \function igraph_similarity_jaccard_pairs
 * \brief Jaccard similarity coefficient for given vertex pairs.
 *
 * The Jaccard similarity coefficient of two vertices is the number of common
 * neighbors divided by the number of vertices that are neighbors of at
 * least one of the two vertices being considered. This function calculates
 * the pairwise Jaccard similarities for a list of vertex pairs.
 *
 * \param graph The graph object to analyze
 * \param res Pointer to a vector, the result of the calculation will
 *        be stored here. The number of elements is the same as the number
 *        of pairs in \p pairs.
 * \param pairs A vector that contains the pairs for which the similarity
 *        will be calculated. Each pair is defined by two consecutive elements,
 *        i.e. the first and second element of the vector specifies the first
 *        pair, the third and fourth element specifies the second pair and so on.
 * \param mode The type of neighbors to be used for the calculation in
 *        directed graphs. Possible values:
 *        \clist
 *        \cli IGRAPH_OUT
 *          the outgoing edges will be considered for each node.
 *        \cli IGRAPH_IN
 *          the incoming edges will be considered for each node.
 *        \cli IGRAPH_ALL
 *          the directed graph is considered as an undirected one for the
 *          computation.
 *        \endclist
 * \param loops Whether to include the vertices themselves in the neighbor
 *        sets.
 * \return Error code:
 *        \clist
 *        \cli IGRAPH_ENOMEM
 *           not enough memory for temporary data.
 *        \cli IGRAPH_EINVVID
 *           invalid vertex ID passed.
 *        \cli IGRAPH_EINVMODE
 *           invalid mode argument.
 *        \endclist
 *
 * Time complexity: O(nd), n is the number of pairs in the given vector, d is
 * the (maximum) degree of the vertices in the graph.
 *
 * \sa \ref igraph_similarity_jaccard() to calculate the Jaccard similarity
 *   between all pairs of a vertex set, or \ref igraph_similarity_dice() and
 *   \ref igraph_similarity_dice_pairs() for a measure very similar to the
 *   Jaccard coefficient
 *
 * \example examples/simple/igraph_similarity.c
 */
igraph_error_t igraph_similarity_jaccard_pairs(const igraph_t *graph, igraph_vector_t *res,
                                    const igraph_vector_int_t *pairs, igraph_neimode_t mode, igraph_bool_t loops) {
    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_lazy_adjlist_t al;
    igraph_integer_t u, v;
    igraph_integer_t len_union, len_intersection;
    igraph_vector_int_t *v1, *v2;

    igraph_integer_t k = igraph_vector_int_size(pairs);
    if (k % 2 != 0) {
        IGRAPH_ERROR("Number of elements in `pairs' must be even.", IGRAPH_EINVAL);
    }
    if (!igraph_vector_int_isininterval(pairs, 0, no_of_nodes - 1)) {
        IGRAPH_ERROR("Invalid vertex ID in pairs.", IGRAPH_EINVVID);
    }
    IGRAPH_CHECK(igraph_vector_resize(res, k / 2));

    IGRAPH_CHECK(igraph_lazy_adjlist_init(graph, &al, mode, IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE));
    IGRAPH_FINALLY(igraph_lazy_adjlist_destroy, &al);

    if (loops) {
        /* Add the loop edges */

        igraph_vector_bool_t seen;
        IGRAPH_VECTOR_BOOL_INIT_FINALLY(&seen, no_of_nodes);

        for (igraph_integer_t i = 0; i < k; i++) {
            igraph_integer_t j = VECTOR(*pairs)[i];
            if (VECTOR(seen)[j]) {
                continue;
            }
            VECTOR(seen)[j] = true;
            v1 = igraph_lazy_adjlist_get(&al, j);
            IGRAPH_CHECK_OOM(v1, "Failed to query neighbors.");
            if (!igraph_vector_int_binsearch(v1, j, &u)) {
                IGRAPH_CHECK(igraph_vector_int_insert(v1, u, j));
            }
        }

        igraph_vector_bool_destroy(&seen);
        IGRAPH_FINALLY_CLEAN(1);
    }

    for (igraph_integer_t i = 0, j = 0; i < k; i += 2, j++) {
        u = VECTOR(*pairs)[i];
        v = VECTOR(*pairs)[i + 1];

        if (u == v) {
            VECTOR(*res)[j] = 1.0;
            continue;
        }

        v1 = igraph_lazy_adjlist_get(&al, u);
        IGRAPH_CHECK_OOM(v1, "Failed to query neighbors.");
        v2 = igraph_lazy_adjlist_get(&al, v);
        IGRAPH_CHECK_OOM(v2, "Failed to query neighbors.");

        IGRAPH_CHECK(igraph_i_neisets_intersect(v1, v2, &len_union, &len_intersection));
        if (len_union > 0) {
            VECTOR(*res)[j] = ((igraph_real_t)len_intersection) / len_union;
        } else {
            VECTOR(*res)[j] = 0.0;
        }
    }

    igraph_lazy_adjlist_destroy(&al);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup structural
 * \function igraph_similarity_jaccard_es
 * \brief Jaccard similarity coefficient for a given edge selector.
 *
 * The Jaccard similarity coefficient of two vertices is the number of common
 * neighbors divided by the number of vertices that are neighbors of at
 * least one of the two vertices being considered. This function calculates
 * the pairwise Jaccard similarities for the endpoints of edges in a given edge
 * selector.
 *
 * \param graph The graph object to analyze
 * \param res Pointer to a vector, the result of the calculation will
 *        be stored here. The number of elements is the same as the number
 *        of edges in \p es.
 * \param es An edge selector that specifies the edges to be included in the
 *        result.
 * \param mode The type of neighbors to be used for the calculation in
 *        directed graphs. Possible values:
 *        \clist
 *        \cli IGRAPH_OUT
 *          the outgoing edges will be considered for each node.
 *        \cli IGRAPH_IN
 *          the incoming edges will be considered for each node.
 *        \cli IGRAPH_ALL
 *          the directed graph is considered as an undirected one for the
 *          computation.
 *        \endclist
 * \param loops Whether to include the vertices themselves in the neighbor
 *        sets.
 * \return Error code:
 *        \clist
 *        \cli IGRAPH_ENOMEM
 *           not enough memory for temporary data.
 *        \cli IGRAPH_EINVVID
 *           invalid vertex ID passed.
 *        \cli IGRAPH_EINVMODE
 *           invalid mode argument.
 *        \endclist
 *
 * Time complexity: O(nd), n is the number of edges in the edge selector, d is
 * the (maximum) degree of the vertices in the graph.
 *
 * \sa \ref igraph_similarity_jaccard() and \ref igraph_similarity_jaccard_pairs()
 *   to calculate the Jaccard similarity between all pairs of a vertex set or
 *   some selected vertex pairs, or \ref igraph_similarity_dice(),
 *   \ref igraph_similarity_dice_pairs() and \ref igraph_similarity_dice_es() for a
 *   measure very similar to the Jaccard coefficient
 *
 * \example examples/simple/igraph_similarity.c
 */
igraph_error_t igraph_similarity_jaccard_es(const igraph_t *graph, igraph_vector_t *res,
                                 const igraph_es_t es, igraph_neimode_t mode, igraph_bool_t loops) {

    igraph_vector_int_t pairs;

    IGRAPH_VECTOR_INT_INIT_FINALLY(&pairs, 0);
    IGRAPH_CHECK(igraph_edges(graph, es, &pairs));
    IGRAPH_CHECK(igraph_similarity_jaccard_pairs(graph, res, &pairs, mode, loops));
    igraph_vector_int_destroy(&pairs);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup structural
 * \function igraph_similarity_dice
 * \brief Dice similarity coefficient.
 *
 * The Dice similarity coefficient of two vertices is twice the number of common
 * neighbors divided by the sum of the degrees of the vertices. This function
 * calculates the pairwise Dice similarities for some (or all) of the vertices.
 *
 * \param graph The graph object to analyze.
 * \param res Pointer to a matrix, the result of the calculation will
 *        be stored here. The number of its rows and columns is the same
 *        as the number of vertex IDs in \p vids.
 * \param vids The vertex IDs of the vertices for which the
 *        calculation will be done.
 * \param mode The type of neighbors to be used for the calculation in
 *        directed graphs. Possible values:
 *        \clist
 *        \cli IGRAPH_OUT
 *          the outgoing edges will be considered for each node.
 *        \cli IGRAPH_IN
 *          the incoming edges will be considered for each node.
 *        \cli IGRAPH_ALL
 *          the directed graph is considered as an undirected one for the
 *          computation.
 *        \endclist
 * \param loops Whether to include the vertices themselves as their own
 *        neighbors.
 * \return Error code:
 *        \clist
 *        \cli IGRAPH_ENOMEM
 *           not enough memory for temporary data.
 *        \cli IGRAPH_EINVVID
 *           invalid vertex ID passed.
 *        \cli IGRAPH_EINVMODE
 *           invalid mode argument.
 *        \endclist
 *
 * Time complexity: O(|V|^2 d),
 * where |V| is the number of vertices in the vertex iterator given, and
 * d is the (maximum) degree of the vertices in the graph.
 *
 * \sa \ref igraph_similarity_jaccard(), a measure very similar to the Dice
 *   coefficient
 *
 * \example examples/simple/igraph_similarity.c
 */
igraph_error_t igraph_similarity_dice(const igraph_t *graph, igraph_matrix_t *res,
                                      const igraph_vs_t vids,
                                      igraph_neimode_t mode, igraph_bool_t loops) {

    IGRAPH_CHECK(igraph_similarity_jaccard(graph, res, vids, mode, loops));

    igraph_integer_t nr = igraph_matrix_nrow(res);
    igraph_integer_t nc = igraph_matrix_ncol(res);
    for (igraph_integer_t i = 0; i < nr; i++) {
        for (igraph_integer_t j = 0; j < nc; j++) {
            igraph_real_t x = MATRIX(*res, i, j);
            MATRIX(*res, i, j) = 2 * x / (1 + x);
        }
    }

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup structural
 * \function igraph_similarity_dice_pairs
 * \brief Dice similarity coefficient for given vertex pairs.
 *
 * The Dice similarity coefficient of two vertices is twice the number of common
 * neighbors divided by the sum of the degrees of the vertices. This function
 * calculates the pairwise Dice similarities for a list of vertex pairs.
 *
 * \param graph The graph object to analyze
 * \param res Pointer to a vector, the result of the calculation will
 *        be stored here. The number of elements is the same as the number
 *        of pairs in \p pairs.
 * \param pairs A vector that contains the pairs for which the similarity
 *        will be calculated. Each pair is defined by two consecutive elements,
 *        i.e. the first and second element of the vector specifies the first
 *        pair, the third and fourth element specifies the second pair and so on.
 * \param mode The type of neighbors to be used for the calculation in
 *        directed graphs. Possible values:
 *        \clist
 *        \cli IGRAPH_OUT
 *          the outgoing edges will be considered for each node.
 *        \cli IGRAPH_IN
 *          the incoming edges will be considered for each node.
 *        \cli IGRAPH_ALL
 *          the directed graph is considered as an undirected one for the
 *          computation.
 *        \endclist
 * \param loops Whether to include the vertices themselves as their own
 *        neighbors.
 * \return Error code:
 *        \clist
 *        \cli IGRAPH_ENOMEM
 *           not enough memory for temporary data.
 *        \cli IGRAPH_EINVVID
 *           invalid vertex ID passed.
 *        \cli IGRAPH_EINVMODE
 *           invalid mode argument.
 *        \endclist
 *
 * Time complexity: O(nd), n is the number of pairs in the given vector, d is
 * the (maximum) degree of the vertices in the graph.
 *
 * \sa \ref igraph_similarity_dice() to calculate the Dice similarity
 *   between all pairs of a vertex set, or \ref igraph_similarity_jaccard(),
 *   \ref igraph_similarity_jaccard_pairs() and \ref igraph_similarity_jaccard_es()
 *   for a measure very similar to the Dice coefficient
 *
 * \example examples/simple/igraph_similarity.c
 */
igraph_error_t igraph_similarity_dice_pairs(const igraph_t *graph, igraph_vector_t *res,
                                 const igraph_vector_int_t *pairs, igraph_neimode_t mode, igraph_bool_t loops) {

    IGRAPH_CHECK(igraph_similarity_jaccard_pairs(graph, res, pairs, mode, loops));
    igraph_integer_t n = igraph_vector_size(res);
    for (igraph_integer_t i = 0; i < n; i++) {
        igraph_real_t x = VECTOR(*res)[i];
        VECTOR(*res)[i] = 2 * x / (1 + x);
    }

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup structural
 * \function igraph_similarity_dice_es
 * \brief Dice similarity coefficient for a given edge selector.
 *
 * The Dice similarity coefficient of two vertices is twice the number of common
 * neighbors divided by the sum of the degrees of the vertices. This function
 * calculates the pairwise Dice similarities for the endpoints of edges in a given
 * edge selector.
 *
 * \param graph The graph object to analyze
 * \param res Pointer to a vector, the result of the calculation will
 *        be stored here. The number of elements is the same as the number
 *        of edges in \p es.
 * \param es An edge selector that specifies the edges to be included in the
 *        result.
 * \param mode The type of neighbors to be used for the calculation in
 *        directed graphs. Possible values:
 *        \clist
 *        \cli IGRAPH_OUT
 *          the outgoing edges will be considered for each node.
 *        \cli IGRAPH_IN
 *          the incoming edges will be considered for each node.
 *        \cli IGRAPH_ALL
 *          the directed graph is considered as an undirected one for the
 *          computation.
 *        \endclist
 * \param loops Whether to include the vertices themselves as their own
 *        neighbors.
 * \return Error code:
 *        \clist
 *        \cli IGRAPH_ENOMEM
 *           not enough memory for temporary data.
 *        \cli IGRAPH_EINVVID
 *           invalid vertex ID passed.
 *        \cli IGRAPH_EINVMODE
 *           invalid mode argument.
 *        \endclist
 *
 * Time complexity: O(nd), n is the number of pairs in the given vector, d is
 * the (maximum) degree of the vertices in the graph.
 *
 * \sa \ref igraph_similarity_dice() and \ref igraph_similarity_dice_pairs()
 *   to calculate the Dice similarity between all pairs of a vertex set or
 *   some selected vertex pairs, or \ref igraph_similarity_jaccard(),
 *   \ref igraph_similarity_jaccard_pairs() and \ref igraph_similarity_jaccard_es()
 *   for a measure very similar to the Dice coefficient
 *
 * \example examples/simple/igraph_similarity.c
 */
igraph_error_t igraph_similarity_dice_es(const igraph_t *graph, igraph_vector_t *res,
                              const igraph_es_t es, igraph_neimode_t mode, igraph_bool_t loops) {

    IGRAPH_CHECK(igraph_similarity_jaccard_es(graph, res, es, mode, loops));
    igraph_integer_t n = igraph_vector_size(res);
    for (igraph_integer_t i = 0; i < n; i++) {
        igraph_real_t x = VECTOR(*res)[i];
        VECTOR(*res)[i] = 2 * x / (1 + x);
    }

    return IGRAPH_SUCCESS;
}
