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
 * \param graph The input graph.
 * \param res Pointer to an initialized boolean vector for storing the result.
 *
 * \sa \ref igraph_simplify() to get rid of loop edges.
 *
 * Time complexity: O(e), the number of edges to check.
 *
 * \example examples/simple/igraph_has_loop.c
 */
int igraph_has_loop(const igraph_t *graph, igraph_bool_t *res) {
    long int i, m = igraph_ecount(graph);

    *res = 0;

    for (i = 0; i < m; i++) {
        if (IGRAPH_FROM(graph, i) == IGRAPH_TO(graph, i)) {
            *res = 1;
            break;
        }
    }

    return 0;
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
int igraph_is_loop(const igraph_t *graph, igraph_vector_bool_t *res,
                   igraph_es_t es) {
    igraph_eit_t eit;
    long int i;

    IGRAPH_CHECK(igraph_eit_create(graph, es, &eit));
    IGRAPH_FINALLY(igraph_eit_destroy, &eit);

    IGRAPH_CHECK(igraph_vector_bool_resize(res, IGRAPH_EIT_SIZE(eit)));

    for (i = 0; !IGRAPH_EIT_END(eit); i++, IGRAPH_EIT_NEXT(eit)) {
        long int e = IGRAPH_EIT_GET(eit);
        VECTOR(*res)[i] = (IGRAPH_FROM(graph, e) == IGRAPH_TO(graph, e)) ? 1 : 0;
    }

    igraph_eit_destroy(&eit);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}
