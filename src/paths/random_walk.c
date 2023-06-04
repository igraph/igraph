/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2014  Gabor Csardi <csardi.gabor@gmail.com>
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

#include "igraph_paths.h"

#include "igraph_adjlist.h"
#include "igraph_interface.h"
#include "igraph_random.h"
#include "igraph_memory.h"
#include "igraph_vector_ptr.h"

#include "core/interruption.h"

/**
 * This function performs a random walk with a given length on a graph,
 * from the given start vertex.
 * It's used for igraph_random_walk when the given graph is unweighted,
 * and only vertex IDs of the vertices on the walk are needed (edge IDs are not needed).
 * \param vertices An allocated vector, the result is stored here as
 *   a list of vertex IDs. It will be resized as needed.
 *   It includes the starting vertex id as well.
 */
static igraph_error_t igraph_i_random_walk_adjlist(const igraph_t *graph,
                                    igraph_vector_int_t *vertices,
                                    igraph_integer_t start,
                                    igraph_neimode_t mode,
                                    igraph_integer_t steps,
                                    igraph_random_walk_stuck_t stuck) {
    igraph_integer_t i;
    igraph_lazy_adjlist_t adj;

    if (vertices == NULL) {
        /* Nothing to do */
        return IGRAPH_SUCCESS;
    }

    IGRAPH_CHECK(igraph_lazy_adjlist_init(graph, &adj, mode, IGRAPH_LOOPS, IGRAPH_MULTIPLE));
    IGRAPH_FINALLY(igraph_lazy_adjlist_destroy, &adj);

    IGRAPH_CHECK(igraph_vector_int_resize(vertices, steps + 1));

    RNG_BEGIN();

    VECTOR(*vertices)[0] = start;
    for (i = 1; i <= steps; i++) {
        igraph_vector_int_t *neis;
        igraph_integer_t nn;
        neis = igraph_lazy_adjlist_get(&adj, start);

        IGRAPH_CHECK_OOM(neis, "Failed to query neighbors.");

        nn = igraph_vector_int_size(neis);
        if (IGRAPH_UNLIKELY(nn == 0)) {
            igraph_vector_int_resize(vertices, i); /* shrinks */
            if (stuck == IGRAPH_RANDOM_WALK_STUCK_RETURN) {
                break;
            } else {
                IGRAPH_ERROR("Random walk got stuck.", IGRAPH_ERWSTUCK);
            }
        }
        start = VECTOR(*vertices)[i] = VECTOR(*neis)[RNG_INTEGER(0, nn - 1)];

        IGRAPH_ALLOW_INTERRUPTION();
    }

    RNG_END();

    igraph_lazy_adjlist_destroy(&adj);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}


/* Used as item destructor for 'cdfs' in igraph_i_random_walk_inclist(). */
static void vec_destr(igraph_vector_t *vec) {
    if (vec != NULL) {
        igraph_vector_destroy(vec);
    }
}


/**
 * This function performs a random walk with a given length on a graph,
 * from the given start vertex.
 * It's used for igraph_random_walk:
 *  - when weights are used or when edge IDs of the traversed edges
 *    and/or vertex IDs of the visited vertices are requested.
 * \param weights A vector of non-negative edge weights. It is assumed
 *   that at least one strictly positive weight is found among the
 *   outgoing edges of each vertex. Additionally, no edge weight may
 *   be NaN. If either case does not hold, an error is returned. If it
 *   is a NULL pointer, all edges are considered to have equal weight.
 * \param vertices An allocated vector, the result is stored here as
 *   a list of vertex IDs. It will be resized as needed.
 *   It includes the starting vertex id as well.
 * \param edges An initialized vector, the indices of traversed
 *   edges are stored here. It will be resized as needed.
 */
static igraph_error_t igraph_i_random_walk_inclist(
        const igraph_t *graph,
        const igraph_vector_t *weights,
        igraph_vector_int_t *vertices,
        igraph_vector_int_t *edges,
        igraph_integer_t start,
        igraph_neimode_t mode,
        igraph_integer_t steps,
        igraph_random_walk_stuck_t stuck) {

    igraph_integer_t vc = igraph_vcount(graph);
    igraph_integer_t i, next;
    igraph_vector_t weight_temp;
    igraph_lazy_inclist_t il;
    igraph_vector_ptr_t cdfs; /* cumulative distribution vectors for each node, used for weighted choice */

    if (vertices) {
        IGRAPH_CHECK(igraph_vector_int_resize(vertices, steps + 1)); /* size: steps + 1 because vertices includes start vertex */
    }
    if (edges) {
        IGRAPH_CHECK(igraph_vector_int_resize(edges, steps));
    }

    IGRAPH_CHECK(igraph_lazy_inclist_init(graph, &il, mode, IGRAPH_LOOPS));
    IGRAPH_FINALLY(igraph_lazy_inclist_destroy, &il);

    IGRAPH_VECTOR_INIT_FINALLY(&weight_temp, 0);

    /* cdf vectors will be computed lazily; that's why we are still using
     * igraph_vector_ptr_t as it does not require us to pre-initialize all
     * the vectors in the vector list */
    IGRAPH_CHECK(igraph_vector_ptr_init(&cdfs, vc));
    IGRAPH_FINALLY(igraph_vector_ptr_destroy_all, &cdfs);
    IGRAPH_VECTOR_PTR_SET_ITEM_DESTRUCTOR(&cdfs, vec_destr);
    for (i = 0; i < vc; ++i) {
        VECTOR(cdfs)[i] = NULL;
    }

    RNG_BEGIN();

    if (vertices) {
        VECTOR(*vertices)[0] = start;
    }
    for (i = 0; i < steps; ++i) {
        igraph_integer_t degree, edge, idx;
        igraph_vector_int_t *inc_edges = igraph_lazy_inclist_get(&il, start);

        IGRAPH_CHECK_OOM(inc_edges, "Failed to query incident edges.");
        degree = igraph_vector_int_size(inc_edges);

        /* are we stuck? */
        if (IGRAPH_UNLIKELY(degree == 0)) {
            /* can't fail since size is reduced, skip IGRAPH_CHECK */
            if (vertices) {
                igraph_vector_int_resize(vertices, i + 1); /* size: i + 1 because vertices includes start vertex */
            }
            if (edges) {
                igraph_vector_int_resize(edges, i);
            }
            if (stuck == IGRAPH_RANDOM_WALK_STUCK_RETURN) {
                break;
            } else {
                IGRAPH_ERROR("Random walk got stuck.", IGRAPH_ERWSTUCK);
            }
        }

        if (weights) { /* weighted: choose an out-edge with probability proportional to its weight */
            igraph_real_t r;
            igraph_vector_t **cd = (igraph_vector_t**) &(VECTOR(cdfs)[start]);

            /* compute out-edge cdf for this node if not already done */
            if (IGRAPH_UNLIKELY(! *cd)) {
                igraph_integer_t j;

                *cd = IGRAPH_CALLOC(1, igraph_vector_t);
                IGRAPH_CHECK_OOM(*cd, "Insufficient memory for random walk.");
                IGRAPH_CHECK(igraph_vector_init(*cd, degree));

                IGRAPH_CHECK(igraph_vector_resize(&weight_temp, degree));
                for (j = 0; j < degree; ++j) {
                    VECTOR(weight_temp)[j] = VECTOR(*weights)[VECTOR(*inc_edges)[j]];
                }

                IGRAPH_CHECK(igraph_vector_cumsum(*cd, &weight_temp));
            }

            r = RNG_UNIF(0, VECTOR(**cd)[degree - 1]);
            igraph_vector_binsearch(*cd, r, &idx);
        }
        else {
            idx = RNG_INTEGER(0, degree - 1);
        }

        edge = VECTOR(*inc_edges)[idx];
        if (edges) {
            VECTOR(*edges)[i] = edge;
        }

        /* travel along edge in a direction specified by 'mode' */
        /* note: 'mode' is always set to IGRAPH_ALL for undirected graphs */
        switch (mode) {
        case IGRAPH_OUT:
            next = IGRAPH_TO(graph, edge);
            break;
        case IGRAPH_IN:
            next = IGRAPH_FROM(graph, edge);
            break;
        case IGRAPH_ALL:
            next = IGRAPH_OTHER(graph, edge, start);
            break;
        }

        if (vertices) {
            VECTOR(*vertices)[i + 1] = next; /* index i + 1 because vertices includes start vertex at position 0 */
        }
        start = next;

        IGRAPH_ALLOW_INTERRUPTION();
    }

    RNG_END();

    igraph_vector_ptr_destroy_all(&cdfs);
    igraph_vector_destroy(&weight_temp);
    igraph_lazy_inclist_destroy(&il);
    IGRAPH_FINALLY_CLEAN(3);

    return IGRAPH_SUCCESS;
}


/**
 * \function igraph_random_walk
 * \brief Performs a random walk on a graph.
 *
 * Performs a random walk with a given length on a graph, from the given
 * start vertex. Edge directions are (potentially) considered, depending on
 * the \p mode argument.
 *
 * \param graph The input graph, it can be directed or undirected.
 *   Multiple edges are respected, so are loop edges.
 * \param weights A vector of non-negative edge weights. It is assumed
 *   that at least one strictly positive weight is found among the
 *   outgoing edges of each vertex. Additionally, no edge weight may
 *   be NaN. If either case does not hold, an error is returned. If it
 *   is \c NULL, all edges are considered to have equal weight.
 * \param vertices An allocated vector, the result is stored here as
 *   a list of vertex IDs. It will be resized as needed.
 *   It includes the vertex IDs of starting and ending vertices.
 *   Length of the vertices vector: \p steps + 1
 * \param edges An initialized vector, the indices of traversed
 *   edges are stored here. It will be resized as needed.
 *   Length of the edges vector: \p steps
 * \param start The start vertex for the walk.
 * \param steps The number of steps to take. If the random walk gets
 *   stuck, then the \p stuck argument specifies what happens.
 *   \p steps is the number of edges to traverse during the walk.
 * \param mode How to walk along the edges in directed graphs.
 *   \c IGRAPH_OUT means following edge directions, \c IGRAPH_IN means
 *   going opposite the edge directions, \c IGRAPH_ALL means ignoring
 *   edge directions. This argument is ignored for undirected graphs.
 * \param stuck What to do if the random walk gets stuck.
 *   \c IGRAPH_RANDOM_WALK_STUCK_RETURN means that the function returns
 *   with a shorter walk; \c IGRAPH_RANDOM_WALK_STUCK_ERROR means
 *   that an \c IGRAPH_ERWSTUCK error is reported.
 *   In both cases, \p vertices and \p edges are truncated to contain
 *   the actual interrupted walk.
 * \return Error code: \c IGRAPH_ERWSTUCK if the walk got stuck.
 *
 * Time complexity:
 *   O(l + d) for unweighted graphs and
 *   O(l * log(k) + d) for weighted graphs,
 *   where \c l is the length of the walk, \c d is the total degree of the visited nodes
 *   and \c k is the average degree of vertices of the given graph.
 */


igraph_error_t igraph_random_walk(const igraph_t *graph,
                       const igraph_vector_t *weights,
                       igraph_vector_int_t *vertices,
                       igraph_vector_int_t *edges,
                       igraph_integer_t start,
                       igraph_neimode_t mode,
                       igraph_integer_t steps,
                       igraph_random_walk_stuck_t stuck) {

    igraph_integer_t vc = igraph_vcount(graph);
    igraph_integer_t ec = igraph_ecount(graph);

    if (!(mode == IGRAPH_ALL || mode == IGRAPH_IN || mode == IGRAPH_OUT)) {
        IGRAPH_ERROR("Invalid mode parameter.", IGRAPH_EINVMODE);
    }

    if (start < 0 || start >= vc) {
        IGRAPH_ERRORF("Starting vertex must be between 0 and the "
                      "number of vertices in the graph (%" IGRAPH_PRId
                      "), got %" IGRAPH_PRId ".", IGRAPH_EINVAL,
                      vc, start);
    }
    if (steps < 0) {
        IGRAPH_ERRORF("Number of steps should be non-negative, got %"
                      IGRAPH_PRId ".", IGRAPH_EINVAL, steps);
    }

    if (weights) {
        if (igraph_vector_size(weights) != ec) {
            IGRAPH_ERROR("Invalid weight vector length.", IGRAPH_EINVAL);
        }
        if (ec > 0) {
            igraph_real_t min = igraph_vector_min(weights);
            if (min < 0) {
                IGRAPH_ERROR("Weights must be non-negative.", IGRAPH_EINVAL);
            } else if (isnan(min)) {
                IGRAPH_ERROR("Weights must not contain NaN values.", IGRAPH_EINVAL);
            }
        }
    }

    if (!igraph_is_directed(graph)) {
        mode = IGRAPH_ALL;
    }

    if (edges || weights) {
        return igraph_i_random_walk_inclist(graph, weights, vertices, edges,
                                            start, mode, steps, stuck);
    } else {
        return igraph_i_random_walk_adjlist(graph, vertices,
                                            start, mode, steps, stuck);
    }
}


/**
 * \function igraph_random_edge_walk
 * \brief Performs a random walk on a graph and returns the traversed edges.
 *
 * Performs a random walk with a given length on a graph, from the given
 * start vertex. Edge directions are (potentially) considered, depending on
 * the \p mode argument.
 *
 * \param graph The input graph, it can be directed or undirected.
 *   Multiple edges are respected, so are loop edges.
 * \param weights A vector of non-negative edge weights. It is assumed
 *   that at least one strictly positive weight is found among the
 *   outgoing edges of each vertex. Additionally, no edge weight may
 *   be NaN. If either case does not hold, an error is returned. If it
 *   is a NULL pointer, all edges are considered to have equal weight.
 * \param edgewalk An initialized vector; the indices of traversed
 *   edges are stored here. It will be resized as needed.
 * \param start The start vertex for the walk.
 * \param steps The number of steps to take. If the random walk gets
 *   stuck, then the \p stuck argument specifies what happens.
 * \param mode How to walk along the edges in directed graphs.
 *   \c IGRAPH_OUT means following edge directions, \c IGRAPH_IN means
 *   going opposite the edge directions, \c IGRAPH_ALL means ignoring
 *   edge directions. This argument is ignored for undirected graphs.
 * \param stuck What to do if the random walk gets stuck.
 *   \c IGRAPH_RANDOM_WALK_STUCK_RETURN means that the function returns
 *   with a shorter walk; \c IGRAPH_RANDOM_WALK_STUCK_ERROR means
 *   that an \c IGRAPH_ERWSTUCK error is reported. In both cases,
 *   \p edgewalk is truncated to contain the actual interrupted walk.
 *
 * \return Error code.
 *
 * \deprecated-by igraph_random_walk 0.10.0
 */
igraph_error_t igraph_random_edge_walk(
        const igraph_t *graph,
        const igraph_vector_t *weights,
        igraph_vector_int_t *edgewalk,
        igraph_integer_t start, igraph_neimode_t mode,
        igraph_integer_t steps,
        igraph_random_walk_stuck_t stuck) {

    return igraph_random_walk(graph, weights, NULL, edgewalk,
                              start, mode, steps, stuck);
}
