/*
   igraph library.
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

#include "igraph_structural.h"

#include "igraph_interface.h"

/**
 * \function igraph_has_loop
 * \brief Returns whether the graph has at least one loop edge.
 *
 * </para><para>
 * A loop edge is an edge from a vertex to itself.
 *
 * </para><para>
 * The return value of this function is cached in the graph itself; calling
 * the function multiple times with no modifications to the graph in between
 * will return a cached value in O(1) time.
 *
 * \param graph The input graph.
 * \param res Pointer to an initialized boolean vector for storing the result.
 * \return Error code.
 *
 * \sa \ref igraph_simplify() to get rid of loop edges.
 *
 * Time complexity: O(e), the number of edges to check.
 *
 * \example examples/simple/igraph_is_loop.c
 */
igraph_error_t igraph_has_loop(const igraph_t *graph, igraph_bool_t *res) {
    igraph_int_t i, m = igraph_ecount(graph);

    IGRAPH_RETURN_IF_CACHED_BOOL(graph, IGRAPH_PROP_HAS_LOOP, res);

    *res = false;

    for (i = 0; i < m; i++) {
        if (IGRAPH_FROM(graph, i) == IGRAPH_TO(graph, i)) {
            *res = true;
            break;
        }
    }

    igraph_i_property_cache_set_bool_checked(graph, IGRAPH_PROP_HAS_LOOP, *res);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_is_loop
 * \brief Find the loop edges in a graph.
 *
 * A loop edge, also called a self-loop, is an edge from a vertex to itself.
 *
 * \param graph The input graph.
 * \param res Pointer to an initialized boolean vector for storing the result,
 *         it will be resized as needed.
 * \param es The edges to check, for all edges supply \ref igraph_ess_all() here.
 * \return Error code.
 *
 * \sa \ref igraph_simplify() to get rid of loop edges.
 *
 * Time complexity: O(e), the number of edges to check.
 *
 * \example examples/simple/igraph_is_loop.c
 */
igraph_error_t igraph_is_loop(const igraph_t *graph, igraph_vector_bool_t *res,
                              igraph_es_t es) {
    igraph_eit_t eit;
    igraph_bool_t found_loop = false;

    IGRAPH_CHECK(igraph_eit_create(graph, es, &eit));
    IGRAPH_FINALLY(igraph_eit_destroy, &eit);

    IGRAPH_CHECK(igraph_vector_bool_resize(res, IGRAPH_EIT_SIZE(eit)));

    if (igraph_i_property_cache_has(graph, IGRAPH_PROP_HAS_LOOP) &&
        ! igraph_i_property_cache_get_bool(graph, IGRAPH_PROP_HAS_LOOP)) {
        igraph_vector_bool_null(res);
        goto done;
    }

    for (igraph_int_t i = 0; !IGRAPH_EIT_END(eit); i++, IGRAPH_EIT_NEXT(eit)) {
        igraph_int_t e = IGRAPH_EIT_GET(eit);
        igraph_bool_t is_loop = (IGRAPH_FROM(graph, e) == IGRAPH_TO(graph, e));
        VECTOR(*res)[i] = is_loop;
        if (is_loop) {
            found_loop = true;
        }
    }

    if (found_loop) {
        igraph_i_property_cache_set_bool_checked(graph, IGRAPH_PROP_HAS_LOOP, true);
    } else if (igraph_es_is_all(&es)) {
        igraph_i_property_cache_set_bool_checked(graph, IGRAPH_PROP_HAS_LOOP, false);
    }

done:
    igraph_eit_destroy(&eit);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_count_loops
 * \brief Counts the self-loops in the graph.
 *
 * Counts loop edges, i.e. edges whose two endpoints coincide.
 *
 * \param graph The input graph.
 * \param loop_count Pointer to an integer, the number of self-loops will be stored here.
 * \return Error code.
 *
 * Time complexity: O(|E|), linear in the number of edges.
 *
 * \example examples/simple/igraph_is_loop.c
 */
igraph_error_t igraph_count_loops(const igraph_t *graph, igraph_int_t *loop_count) {
    igraph_int_t no_of_edges = igraph_ecount(graph);
    igraph_int_t count;

    /* Nothing to do if we know that there are no loops. */
    if (igraph_i_property_cache_has(graph, IGRAPH_PROP_HAS_LOOP) &&
        !igraph_i_property_cache_get_bool(graph, IGRAPH_PROP_HAS_LOOP)) {
        *loop_count = 0;
        return IGRAPH_SUCCESS;
    }

    count = 0;
    for (igraph_int_t e=0; e < no_of_edges; e++) {
        if (IGRAPH_FROM(graph, e) == IGRAPH_TO(graph, e)) {
            count++;
        }
    }

    /* We already checked for loops, so take the opportunity to set the cache. */
    igraph_i_property_cache_set_bool_checked(graph, IGRAPH_PROP_HAS_LOOP, count > 0);

    *loop_count = count;

    return IGRAPH_SUCCESS;
}
