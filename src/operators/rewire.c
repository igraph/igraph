/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
/*
   IGraph library.
   Copyright (C) 2005-2021 The igraph development team

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

#include "igraph_operators.h"

#include "igraph_adjlist.h"
#include "igraph_conversion.h"
#include "igraph_interface.h"
#include "igraph_iterators.h"
#include "igraph_progress.h"
#include "igraph_random.h"
#include "igraph_structural.h"

#include "core/interruption.h"
#include "operators/rewire_internal.h"

/* Threshold that defines when to switch over to using adjacency lists during
 * rewiring */
#define REWIRE_ADJLIST_THRESHOLD 10

/* Not declared static so that the testsuite can use it, but not part of the public API. */
int igraph_i_rewire(igraph_t *graph, igraph_integer_t n, igraph_rewiring_t mode, igraph_bool_t use_adjlist) {
    long int no_of_nodes = igraph_vcount(graph);
    long int no_of_edges = igraph_ecount(graph);
    char message[256];
    igraph_integer_t a, b, c, d, dummy, num_swaps, num_successful_swaps;
    igraph_vector_t eids, edgevec, alledges;
    igraph_bool_t directed, loops, ok;
    igraph_es_t es;
    igraph_adjlist_t al;

    if (no_of_nodes < 4) {
        IGRAPH_ERROR("graph unsuitable for rewiring", IGRAPH_EINVAL);
    }

    directed = igraph_is_directed(graph);
    loops = (mode & IGRAPH_REWIRING_SIMPLE_LOOPS);

    RNG_BEGIN();

    IGRAPH_VECTOR_INIT_FINALLY(&eids, 2);

    if (use_adjlist) {
        /* As well as the sorted adjacency list, we maintain an unordered
         * list of edges for picking a random edge in constant time.
         */
        IGRAPH_CHECK(igraph_adjlist_init(graph, &al, IGRAPH_OUT, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE));
        IGRAPH_FINALLY(igraph_adjlist_destroy, &al);
        IGRAPH_VECTOR_INIT_FINALLY(&alledges, no_of_edges * 2);
        igraph_get_edgelist(graph, &alledges, /*bycol=*/ 0);
    } else {
        IGRAPH_VECTOR_INIT_FINALLY(&edgevec, 4);
        es = igraph_ess_vector(&eids);
    }

    /* We don't want the algorithm to get stuck in an infinite loop when
     * it can't choose two edges satisfying the conditions. Instead of
     * this, we choose two arbitrary edges and if they have endpoints
     * in common, we just decrease the number of trials left and continue
     * (so unsuccessful rewirings still count as a trial)
     */

    num_swaps = num_successful_swaps = 0;
    while (num_swaps < n) {

        IGRAPH_ALLOW_INTERRUPTION();
        if (num_swaps % 1000 == 0) {
            snprintf(message, sizeof(message),
                     "Random rewiring (%.2f%% of the trials were successful)",
                     num_swaps > 0 ? ((100.0 * num_successful_swaps) / num_swaps) : 0.0);
            IGRAPH_PROGRESS(message, (100.0 * num_swaps) / n, 0);
        }

        switch (mode) {
        case IGRAPH_REWIRING_SIMPLE:
        case IGRAPH_REWIRING_SIMPLE_LOOPS:
            ok = 1;

            /* Choose two edges randomly */
            VECTOR(eids)[0] = RNG_INTEGER(0, no_of_edges - 1);
            do {
                VECTOR(eids)[1] = RNG_INTEGER(0, no_of_edges - 1);
            } while (VECTOR(eids)[0] == VECTOR(eids)[1]);

            /* Get the endpoints */
            if (use_adjlist) {
                a = VECTOR(alledges)[((igraph_integer_t)VECTOR(eids)[0]) * 2];
                b = VECTOR(alledges)[(((igraph_integer_t)VECTOR(eids)[0]) * 2) + 1];
                c = VECTOR(alledges)[((igraph_integer_t)VECTOR(eids)[1]) * 2];
                d = VECTOR(alledges)[(((igraph_integer_t)VECTOR(eids)[1]) * 2) + 1];
            } else {
                IGRAPH_CHECK(igraph_edge(graph, (igraph_integer_t) VECTOR(eids)[0],
                                         &a, &b));
                IGRAPH_CHECK(igraph_edge(graph, (igraph_integer_t) VECTOR(eids)[1],
                                         &c, &d));
            }

            /* For an undirected graph, we have two "variants" of each edge, i.e.
             * a -- b and b -- a. Since some rewirings can be performed only when we
             * "swap" the endpoints, we do it now with probability 0.5 */
            if (!directed && RNG_UNIF01() < 0.5) {
                dummy = c; c = d; d = dummy;
                if (use_adjlist) {
                    /* Flip the edge in the unordered edge-list, so the update later on
                     * hits the correct end. */
                    VECTOR(alledges)[((igraph_integer_t)VECTOR(eids)[1]) * 2] = c;
                    VECTOR(alledges)[(((igraph_integer_t)VECTOR(eids)[1]) * 2) + 1] = d;
                }
            }

            /* If we do not touch loops, check whether a == b or c == d and disallow
             * the swap if needed */
            if (!loops && (a == b || c == d)) {
                ok = 0;
            } else {
                /* Check whether they are suitable for rewiring */
                if (a == c || b == d) {
                    /* Swapping would have no effect */
                    ok = 0;
                } else {
                    /* a != c && b != d */
                    /* If a == d or b == c, the swap would generate at least one loop, so
                     * we disallow them unless we want to have loops */
                    ok = loops || (a != d && b != c);
                    /* Also, if a == b and c == d and we allow loops, doing the swap
                     * would result in a multiple edge if the graph is undirected */
                    ok = ok && (directed || a != b || c != d);
                }
            }

            /* All good so far. Now check for the existence of a --> d and c --> b to
             * disallow the creation of multiple edges */
            if (ok) {
                if (use_adjlist) {
                    if (igraph_adjlist_has_edge(&al, a, d, directed)) {
                        ok = 0;
                    }
                } else {
                    IGRAPH_CHECK(igraph_are_connected(graph, a, d, &ok));
                    ok = !ok;
                }
            }
            if (ok) {
                if (use_adjlist) {
                    if (igraph_adjlist_has_edge(&al, c, b, directed)) {
                        ok = 0;
                    }
                } else {
                    IGRAPH_CHECK(igraph_are_connected(graph, c, b, &ok));
                    ok = !ok;
                }
            }

            /* If we are still okay, we can perform the rewiring */
            if (ok) {
                /* printf("Deleting: %ld -> %ld, %ld -> %ld\n",
                              (long)a, (long)b, (long)c, (long)d); */
                if (use_adjlist) {
                    /* Replace entry in sorted adjlist: */
                    IGRAPH_CHECK(igraph_adjlist_replace_edge(&al, a, b, d, directed));
                    IGRAPH_CHECK(igraph_adjlist_replace_edge(&al, c, d, b, directed));
                    /* Also replace in unsorted edgelist: */
                    VECTOR(alledges)[(((igraph_integer_t)VECTOR(eids)[0]) * 2) + 1] = d;
                    VECTOR(alledges)[(((igraph_integer_t)VECTOR(eids)[1]) * 2) + 1] = b;
                } else {
                    IGRAPH_CHECK(igraph_delete_edges(graph, es));
                    VECTOR(edgevec)[0] = a; VECTOR(edgevec)[1] = d;
                    VECTOR(edgevec)[2] = c; VECTOR(edgevec)[3] = b;
                    /* printf("Adding: %ld -> %ld, %ld -> %ld\n",
                                (long)a, (long)d, (long)c, (long)b); */
                    igraph_add_edges(graph, &edgevec, 0);
                }
                num_successful_swaps++;
            }
            break;
        default:
            RNG_END();
            IGRAPH_ERROR("unknown rewiring mode", IGRAPH_EINVMODE);
        }
        num_swaps++;
    }

    if (use_adjlist) {
        /* Replace graph edges with the adjlist current state */
        IGRAPH_CHECK(igraph_delete_edges(graph, igraph_ess_all(IGRAPH_EDGEORDER_ID)));
        IGRAPH_CHECK(igraph_add_edges(graph, &alledges, 0));
    }

    IGRAPH_PROGRESS("Random rewiring: ", 100.0, 0);

    if (use_adjlist) {
        igraph_vector_destroy(&alledges);
        igraph_adjlist_destroy(&al);
    } else {
        igraph_vector_destroy(&edgevec);
    }

    igraph_vector_destroy(&eids);
    IGRAPH_FINALLY_CLEAN(use_adjlist ? 3 : 2);

    RNG_END();

    return 0;
}

/**
 * \ingroup structural
 * \function igraph_rewire
 * \brief Randomly rewires a graph while preserving the degree distribution.
 *
 * </para><para>
 * This function generates a new graph based on the original one by randomly
 * rewiring edges while preserving the original graph's degree distribution.
 * Please note that the rewiring is done "in place", so no new graph will
 * be allocated. If you would like to keep the original graph intact, use
 * \ref igraph_copy() beforehand.
 *
 * \param graph The graph object to be rewired.
 * \param n Number of rewiring trials to perform.
 * \param mode The rewiring algorithm to be used. It can be one of the following flags:
 *         \clist
 *           \cli IGRAPH_REWIRING_SIMPLE
 *                Simple rewiring algorithm which chooses two arbitrary edges
 *                in each step (namely (a,b) and (c,d)) and substitutes them
 *                with (a,d) and (c,b) if they don't exist.  The method will
 *                neither destroy nor create self-loops.
 *           \cli IGRAPH_REWIRING_SIMPLE_LOOPS
 *                Same as \c IGRAPH_REWIRING_SIMPLE but allows the creation or
 *                destruction of self-loops.
 *         \endclist
 *
 * \return Error code:
 *         \clist
 *           \cli IGRAPH_EINVMODE
 *                Invalid rewiring mode.
 *           \cli IGRAPH_EINVAL
 *                Graph unsuitable for rewiring (e.g. it has
 *                less than 4 nodes in case of \c IGRAPH_REWIRING_SIMPLE)
 *           \cli IGRAPH_ENOMEM
 *                Not enough memory for temporary data.
 *         \endclist
 *
 * Time complexity: TODO.
 */
int igraph_rewire(igraph_t *graph, igraph_integer_t n, igraph_rewiring_t mode) {
    igraph_bool_t use_adjlist = n >= REWIRE_ADJLIST_THRESHOLD;
    return igraph_i_rewire(graph, n, mode, use_adjlist);
}
