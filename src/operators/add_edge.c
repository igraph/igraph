/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
/*
   IGraph library.
   Copyright (C) 2005-2020 The igraph development team

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

#include "igraph_interface.h"

/**
 * \function igraph_add_edge
 * \brief Adds a single edge to a graph.
 *
 * </para><para>
 * For directed graphs the edge points from \p from to \p to.
 *
 * </para><para>
 * Note that if you want to add many edges to a big graph, then it is
 * inefficient to add them one by one, it is better to collect them into
 * a vector and add all of them via a single \ref igraph_add_edges() call.
 * \param igraph The graph.
 * \param from The id of the first vertex of the edge.
 * \param to The id of the second vertex of the edge.
 * \return Error code.
 *
 * \sa \ref igraph_add_edges() to add many edges, \ref
 * igraph_delete_edges() to remove edges and \ref
 * igraph_add_vertices() to add vertices.
 *
 * Time complexity: O(|V|+|E|), the number of edges plus the number of
 * vertices.
 */
int igraph_add_edge(igraph_t *graph, igraph_integer_t from, igraph_integer_t to) {
    igraph_vector_t edges;
    int ret;

    IGRAPH_VECTOR_INIT_FINALLY(&edges, 2);

    VECTOR(edges)[0] = from;
    VECTOR(edges)[1] = to;
    IGRAPH_CHECK(ret = igraph_add_edges(graph, &edges, 0));

    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);
    return ret;
}
