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
#include "core/interruption.h"

/**
 * \function igraph_random_walk
 * Perform a random walk on a graph
 *
 * Performs a random walk with a given length on a graph, from the given
 * start vertex. Edge directions are (potentially) considered, depending on
 * the \p mode argument.
 *
 * \param graph The input graph, it can be directed or undirected.
 *   Multiple edges are respected, so are loop edges.
 * \param walk An allocated vector, the result is stored here.
 *   It will be resized as needed.
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
 *   that an error is reported. In both cases \p walk is truncated
 *   to contain the actual interrupted walk.
 * \return Error code.
 *
 * Time complexity: O(l + d), where \c l is the length of the
 * walk, and \c d is the total degree of the visited nodes.
 */


int igraph_random_walk(const igraph_t *graph, igraph_vector_t *walk,
                       igraph_integer_t start, igraph_neimode_t mode,
                       igraph_integer_t steps,
                       igraph_random_walk_stuck_t stuck) {

    /* TODO:
       - multiple walks potentially from multiple start vertices
       - weights
    */

    igraph_lazy_adjlist_t adj;
    igraph_integer_t vc = igraph_vcount(graph);
    igraph_integer_t i;

    if (start < 0 || start >= vc) {
        IGRAPH_ERROR("Invalid start vertex", IGRAPH_EINVAL);
    }
    if (steps < 0) {
        IGRAPH_ERROR("Invalid number of steps", IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_lazy_adjlist_init(graph, &adj, mode, IGRAPH_LOOPS, IGRAPH_MULTIPLE));
    IGRAPH_FINALLY(igraph_lazy_adjlist_destroy, &adj);

    IGRAPH_CHECK(igraph_vector_resize(walk, steps));

    RNG_BEGIN();

    VECTOR(*walk)[0] = start;
    for (i = 1; i < steps; i++) {
        igraph_vector_int_t *neis;
        igraph_integer_t nn;
        neis = igraph_lazy_adjlist_get(&adj, start);
        nn = igraph_vector_int_size(neis);

        if (IGRAPH_UNLIKELY(nn == 0)) {
            igraph_vector_resize(walk, i);
            if (stuck == IGRAPH_RANDOM_WALK_STUCK_RETURN) {
                break;
            } else {
                IGRAPH_ERROR("Random walk got stuck", IGRAPH_ERWSTUCK);
            }
        }
        start = VECTOR(*walk)[i] = VECTOR(*neis)[ RNG_INTEGER(0, nn - 1) ];
    }

    RNG_END();

    igraph_lazy_adjlist_destroy(&adj);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}


/* Used as item destructor for 'cdfs' in igraph_random_edge_walk(). */
static void vec_destr(igraph_vector_t *vec) {
    if (vec != NULL) {
        igraph_vector_destroy(vec);
    }
}


/**
 * \function igraph_random_edge_walk
 * \brief Perform a random walk on a graph and return the traversed edges
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
 *   that an error is reported. In both cases, \p edgewalk is truncated
 *   to contain the actual interrupted walk.
 *
 * \return Error code.
 *
 */
int igraph_random_edge_walk(const igraph_t *graph,
                            const igraph_vector_t *weights,
                            igraph_vector_t *edgewalk,
                            igraph_integer_t start, igraph_neimode_t mode,
                            igraph_integer_t steps,
                            igraph_random_walk_stuck_t stuck) {
    igraph_integer_t vc = igraph_vcount(graph);
    igraph_integer_t ec = igraph_ecount(graph);
    igraph_integer_t i;
    igraph_inclist_t il;
    igraph_vector_t weight_temp;
    igraph_vector_ptr_t cdfs; /* cumulative distribution vectors for each node, used for weighted choice */

    /* the fourth igraph_neimode_t value, IGRAPH_TOTAL, is disallowed */
    if (! (mode == IGRAPH_ALL || mode == IGRAPH_IN || mode == IGRAPH_OUT)) {
        IGRAPH_ERROR("Invalid mode parameter", IGRAPH_EINVMODE);
    }

    /* ref switch statement at end of main loop */
    if (! igraph_is_directed(graph)) {
        mode = IGRAPH_ALL;
    }

    if (start < 0 || start >= vc) {
        IGRAPH_ERROR("Invalid start vertex", IGRAPH_EINVAL);
    }

    if (steps < 0) {
        IGRAPH_ERROR("Invalid number of steps", IGRAPH_EINVAL);
    }

    if (weights) {
        if (igraph_vector_size(weights) != ec) {
            IGRAPH_ERROR("Invalid weight vector length", IGRAPH_EINVAL);
        }
        if (ec > 0) {
            igraph_real_t min = igraph_vector_min(weights);
            if (min < 0) {
                IGRAPH_ERROR("Weights must be non-negative", IGRAPH_EINVAL);
            }
            else if (igraph_is_nan(min)) {
                IGRAPH_ERROR("Weights must not contain NaN values", IGRAPH_EINVAL);
            }
        }
    }

    IGRAPH_CHECK(igraph_vector_resize(edgewalk, steps));

    IGRAPH_CHECK(igraph_inclist_init(graph, &il, mode, IGRAPH_LOOPS));
    IGRAPH_FINALLY(igraph_inclist_destroy, &il);

    IGRAPH_VECTOR_INIT_FINALLY(&weight_temp, 0);

    /* cdf vectors will be computed lazily */
    IGRAPH_CHECK(igraph_vector_ptr_init(&cdfs, vc));
    IGRAPH_FINALLY(igraph_vector_ptr_destroy_all, &cdfs);
    IGRAPH_VECTOR_PTR_SET_ITEM_DESTRUCTOR(&cdfs, vec_destr);
    for (i = 0; i < vc; ++i) {
        VECTOR(cdfs)[i] = NULL;
    }

    RNG_BEGIN();

    for (i = 0; i < steps; ++i) {
        long degree, edge, idx;
        igraph_vector_int_t *edges = igraph_inclist_get(&il, start);

        degree = igraph_vector_int_size(edges);

        /* are we stuck? */
        if (IGRAPH_UNLIKELY(degree == 0)) {
            igraph_vector_resize(edgewalk, i); /* can't fail since size is reduced, skip IGRAPH_CHECK */
            if (stuck == IGRAPH_RANDOM_WALK_STUCK_RETURN) {
                break;
            } else {
                IGRAPH_ERROR("Random walk got stuck", IGRAPH_ERWSTUCK);
            }
        }

        if (weights) { /* weighted: choose an out-edge with probability proportional to its weight */
            igraph_real_t r;
            igraph_vector_t **cd = (igraph_vector_t **) & (VECTOR(cdfs)[start]);

            /* compute out-edge cdf for this node if not already done */
            if (IGRAPH_UNLIKELY(! *cd)) {
                long j;

                *cd = IGRAPH_CALLOC(1, igraph_vector_t);
                if (*cd == NULL) {
                    IGRAPH_ERROR("random edge walk failed", IGRAPH_ENOMEM);
                }
                IGRAPH_CHECK(igraph_vector_init(*cd, degree));

                IGRAPH_CHECK(igraph_vector_resize(&weight_temp, degree));
                for (j = 0; j < degree; ++j) {
                    VECTOR(weight_temp)[j] = VECTOR(*weights)[ VECTOR(*edges)[j] ];
                }

                IGRAPH_CHECK(igraph_vector_cumsum(*cd, &weight_temp));
            }

            r = RNG_UNIF(0, VECTOR( **cd )[degree - 1]);
            igraph_vector_binsearch(*cd, r, &idx);
        } else { /* unweighted: choose an out-edge at random */
            idx = RNG_INTEGER(0, degree - 1);
        }

        edge = VECTOR(*edges)[idx];
        VECTOR(*edgewalk)[i] = edge;

        /* travel along edge in a direction specified by 'mode' */
        /* note: 'mode' is always set to IGRAPH_ALL for undirected graphs */
        switch (mode) {
        case IGRAPH_OUT:
            start = IGRAPH_TO(graph, edge);
            break;
        case IGRAPH_IN:
            start = IGRAPH_FROM(graph, edge);
            break;
        case IGRAPH_ALL:
            start = IGRAPH_OTHER(graph, edge, start);
            break;
        }

        IGRAPH_ALLOW_INTERRUPTION();
    }

    RNG_END();

    igraph_vector_ptr_destroy_all(&cdfs);
    igraph_vector_destroy(&weight_temp);
    igraph_inclist_destroy(&il);
    IGRAPH_FINALLY_CLEAN(3);

    return IGRAPH_SUCCESS;
}
