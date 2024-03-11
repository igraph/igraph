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
 *
 * \sa \ref igraph_simplify() to get rid of loop edges.
 *
 * Time complexity: O(e), the number of edges to check.
 *
 * \example examples/simple/igraph_has_loop.c
 */
igraph_error_t igraph_has_loop(const igraph_t *graph, igraph_bool_t *res) {
    igraph_integer_t i, m = igraph_ecount(graph);

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
 * </para><para>
 * A loop edge is an edge from a vertex to itself.
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
    igraph_integer_t i;

    IGRAPH_CHECK(igraph_eit_create(graph, es, &eit));
    IGRAPH_FINALLY(igraph_eit_destroy, &eit);

    IGRAPH_CHECK(igraph_vector_bool_resize(res, IGRAPH_EIT_SIZE(eit)));

    for (i = 0; !IGRAPH_EIT_END(eit); i++, IGRAPH_EIT_NEXT(eit)) {
        igraph_integer_t e = IGRAPH_EIT_GET(eit);
        VECTOR(*res)[i] = (IGRAPH_FROM(graph, e) == IGRAPH_TO(graph, e));
    }

    igraph_eit_destroy(&eit);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}
