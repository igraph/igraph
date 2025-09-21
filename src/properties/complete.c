/*
  igraph library.
  Copyright (C) 2024 The igraph development team <igraph@igraph.org>

  This program is free software; you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free Software
  Foundation; either version 2 of the License, or (at your option) any later
  version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

  You should have received a copy of the GNU General Public License along with
  this program. If not, see <https://www.gnu.org/licenses/>.
*/


#include "igraph_interface.h"
#include "igraph_structural.h"

#include "core/interruption.h"

/**
 * \ingroup structural
 * \function igraph_is_complete
 * \brief Decides whether the graph is complete.
 *
 * A graph is considered complete if all pairs of different vertices are
 * adjacent.
 *
 * </para><para>
 * The null graph and the singleton graph are considered complete.
 *
 * \param graph The graph object to analyze.
 * \param res Pointer to a Boolean variable, the result will be stored here.
 *
 * \return Error code.
 *
 * Time complexity: O(|V| + |E|) at worst.
 */

igraph_error_t igraph_is_complete(const igraph_t *graph, igraph_bool_t *res) {

    const igraph_int_t vcount = igraph_vcount(graph);
    const igraph_int_t ecount = igraph_ecount(graph);
    igraph_int_t complete_ecount;
    igraph_bool_t simple, directed = igraph_is_directed(graph);
    igraph_vector_int_t neighbours;
    int iter = 0;

    /* If the graph is the null graph or the singleton graph, return early */
    if (vcount == 0 || vcount == 1) {
        *res = true;
        return IGRAPH_SUCCESS;
    }

    /* Compute the amount of edges a complete graph of vcount vertices would
       have */

    /* Depends on whether the graph is directed */

    /* We have to take care of integer overflowing */

#if IGRAPH_INTEGER_SIZE == 32
    if (directed) {
        /* Highest x s.t. x² - x < 2^31 - 1 */
        if (vcount > 46341) {
            *res = false;
            return IGRAPH_SUCCESS;
        } else {
            complete_ecount = vcount * (vcount - 1);
        }
    } else {
        /* Highest x s.t. (x² - x) / 2 < 2^31 - 1 */
        if (vcount > 65536) {
            *res = false;
            return IGRAPH_SUCCESS;
        } else {
            complete_ecount = vcount % 2 == 0 ?
                              (vcount / 2) * (vcount - 1) :
                              vcount * ((vcount - 1) / 2);
        }
    }
#elif IGRAPH_INTEGER_SIZE == 64
    if (directed) {
        /* Highest x s.t. x² - x < 2^63 - 1 */
        if (vcount > 3037000500) {
            *res = false;
            return IGRAPH_SUCCESS;
        } else {
            complete_ecount = vcount * (vcount - 1);
        }
    } else {
        /* Highest x s.t. (x² - x) / 2 < 2^63 - 1 */
        if (vcount > 4294967296) {
            *res = false;
            return IGRAPH_SUCCESS;
        } else {
            complete_ecount = vcount % 2 == 0 ?
                              (vcount / 2) * (vcount - 1) :
                              vcount * ((vcount - 1) / 2);
        }
    }
#else
    /* If values other than 32 or 64 become allowed,
     * this code will need to be updated. */
#  error "Unexpected IGRAPH_INTEGER_SIZE value."
#endif

    /* If the amount of edges is strictly lower than what it should be for a
       complete graph, return early */

    if (ecount < complete_ecount) {
        *res = false;
        return IGRAPH_SUCCESS;
    }

    /* If the graph is simple, compare and conclude */
    IGRAPH_CHECK(igraph_is_simple(graph, &simple, IGRAPH_DIRECTED));

    if (simple) {
        *res = (ecount == complete_ecount);
        return IGRAPH_SUCCESS;
    }

    /* Allocate memory for vector of size v */
    IGRAPH_VECTOR_INT_INIT_FINALLY(&neighbours, vcount);

    for (igraph_int_t i = 0; i < vcount; ++i) {
        IGRAPH_ALLOW_INTERRUPTION_LIMITED(iter, 1 << 8);

        IGRAPH_CHECK(igraph_neighbors(
            graph, &neighbours, i, IGRAPH_OUT, IGRAPH_NO_LOOPS,
            IGRAPH_NO_MULTIPLE
        ));

        if ((igraph_vector_int_size(&neighbours) < vcount - 1)) {
            *res = false;
            goto cleanup;
        }
    }

    /* If we arrive here, we have found no neighbour vector of size strictly
       less than vcount - 1. The graph is therefore complete */

    *res = true;

cleanup:

    igraph_vector_int_destroy(&neighbours);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}


/* Test for cliques or independent sets, depending on whether independent_set == true. */
static igraph_error_t is_clique(const igraph_t *graph, igraph_vs_t candidate,
                                igraph_bool_t directed, igraph_bool_t *res,
                                igraph_bool_t independent_set) {
    igraph_vector_int_t vids;
    igraph_int_t n; /* clique size */
    igraph_bool_t result = true; /* be optimistic */
    int iter = 0;

    /* The following implementation is optimized for testing for small cliques
     * in large graphs. */

    IGRAPH_VECTOR_INT_INIT_FINALLY(&vids, 0);
    IGRAPH_CHECK(igraph_vs_as_vector(graph, candidate, &vids));

    n = igraph_vector_int_size(&vids);

    for (igraph_int_t i = 0; i < n; i++) {
        igraph_int_t u = VECTOR(vids)[i];
        for (igraph_int_t j = directed ? 0 : i+1; j < n; j++) {
            igraph_int_t v = VECTOR(vids)[j];
            /* Compare u and v for equality instead of i and j in case
             * the vertex list contained duplicates. */
            if (u != v) {
                igraph_int_t eid;
                IGRAPH_CHECK(igraph_get_eid(graph, &eid, u, v, directed, false));
                if (independent_set) {
                    if (eid != -1) {
                        result = false;
                        goto done;
                    }
                } else {
                    if (eid == -1) {
                        result = false;
                        goto done;
                    }
                }
            }
        }
        IGRAPH_ALLOW_INTERRUPTION_LIMITED(iter, 1 << 8);
    }

done:

    *res = result;

    igraph_vector_int_destroy(&vids);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup structural
 * \function igraph_is_clique
 * \brief Does a set of vertices form a clique?
 *
 * Tests if all pairs within a set of vertices are adjacent, i.e. whether they
 * form a clique. An empty set and singleton set are considered to be a clique.
 *
 * \param graph The input graph.
 * \param candidate The vertex set to test for being a clique.
 * \param directed Whether to take edge directions into account in directed graphs.
 * \param res The result will be stored here.
 * \return Error code.
 *
 * \sa \ref igraph_is_complete() to test if a graph is complete;
 * \ref igraph_is_independent_vertex_set() to test for independent vertex sets;
 * \ref igraph_cliques(), \ref igraph_maximal_cliques() and
 * \ref igraph_largest_cliques() to find cliques.
 *
 * Time complexity: O(n^2 log(d)) where n is the number of vertices in the
 * candidate set and d is the typical vertex degree.
 */
igraph_error_t igraph_is_clique(const igraph_t *graph, igraph_vs_t candidate,
                                igraph_bool_t directed, igraph_bool_t *res) {

    if (! igraph_is_directed(graph)) {
        directed = false;
    }

    if (igraph_is_directed(graph) == directed && igraph_vs_is_all(&candidate)) {
        return igraph_is_complete(graph, res);
    }

    return is_clique(graph, candidate, directed, res, /* independent_set */ false);
}

/**
 * \ingroup structural
 * \function igraph_is_independent_vertex_set
 * \brief Does a set of vertices form an independent set?
 *
 * Tests if no pairs within a set of vertices are adjacenct, i.e. whether they
 * form an independent set. An empty set and singleton set are both considered
 * to be an independent set.
 *
 * \param graph The input graph.
 * \param candidate The vertex set to test for being an independent set.
 * \param res The result will be stored here.
 * \return Error code.
 *
 * \sa \ref igraph_is_clique() to test for cliques; \ref igraph_independent_vertex_sets(),
 * \ref igraph_maximal_independent_vertex_sets() and
 * \ref igraph_largest_independent_vertex_sets() to find independent vertex sets.
 *
 * Time complexity: O(n^2 log(d)) where n is the number of vertices in the
 * candidate set and d is the typical vertex degree.
 */
igraph_error_t igraph_is_independent_vertex_set(const igraph_t *graph, igraph_vs_t candidate,
                                         igraph_bool_t *res) {

    /* Note: igraph_count_loops() already makes use of the cache. */
    if (igraph_vs_is_all(&candidate)) {
        igraph_int_t loop_count;
        igraph_count_loops(graph, &loop_count);
        *res = (igraph_ecount(graph) - loop_count) == 0;
        return IGRAPH_SUCCESS;
    }

    return is_clique(graph, candidate, /* directed */ false, res, /* independent_set */ true);
}
