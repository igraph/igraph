/* -*- mode: C -*-  */
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

#include "igraph_constructors.h"

#include "igraph_interface.h"

/**
 * \section about_generators
 *
 * <para>Graph generators create graphs.</para>
 *
 * <para>Almost all functions which create graph objects are documented
 * here. The exceptions are \ref igraph_induced_subgraph() and alike, these
 * create graphs based on another graph.</para>
 */


/**
 * \ingroup generators
 * \function igraph_create
 * \brief Creates a graph with the specified edges.
 *
 * \param graph An uninitialized graph object.
 * \param edges The edges to add, the first two elements are the first
 *        edge, etc.
 * \param n The number of vertices in the graph, if smaller or equal
 *        to the highest vertex id in the \p edges vector it
 *        will be increased automatically. So it is safe to give 0
 *        here.
 * \param directed Boolean, whether to create a directed graph or
 *        not. If yes, then the first edge points from the first
 *        vertex id in \p edges to the second, etc.
 * \return Error code:
 *         \c IGRAPH_EINVEVECTOR: invalid edges
 *         vector (odd number of vertices).
 *         \c IGRAPH_EINVVID: invalid (negative)
 *         vertex id.
 *
 * Time complexity: O(|V|+|E|),
 * |V| is the number of vertices,
 * |E| the number of edges in the
 * graph.
 *
 * \example examples/simple/igraph_create.c
 */
int igraph_create(igraph_t *graph, const igraph_vector_t *edges,
                  igraph_integer_t n, igraph_bool_t directed) {
    igraph_bool_t has_edges = igraph_vector_size(edges) > 0;
    igraph_real_t max = has_edges ? igraph_vector_max(edges) + 1 : 0;

    if (igraph_vector_size(edges) % 2 != 0) {
        IGRAPH_ERROR("Invalid (odd) edges vector", IGRAPH_EINVEVECTOR);
    }
    if (has_edges && !igraph_vector_isininterval(edges, 0, max - 1)) {
        IGRAPH_ERROR("Invalid (negative) vertex id", IGRAPH_EINVVID);
    }

    IGRAPH_CHECK(igraph_empty(graph, n, directed));
    IGRAPH_FINALLY(igraph_destroy, graph);
    if (has_edges) {
        igraph_integer_t vc = igraph_vcount(graph);
        if (vc < max) {
            IGRAPH_CHECK(igraph_add_vertices(graph, (igraph_integer_t) (max - vc), 0));
        }
        IGRAPH_CHECK(igraph_add_edges(graph, edges, 0));
    }

    IGRAPH_FINALLY_CLEAN(1);
    return 0;
}

/**
 * \function igraph_small
 * \brief Shorthand to create a small graph, giving the edges as arguments.
 *
 * </para><para>
 * This function is handy when a relatively small graph needs to be created.
 * Instead of giving the edges as a vector, they are given simply as
 * arguments and a '-1' needs to be given after the last meaningful
 * edge argument.
 *
 * </para><para>Note that only graphs which have vertices less than
 * the highest value of the 'int' type can be created this way. If you
 * give larger values then the result is undefined.
 *
 * \param graph Pointer to an uninitialized graph object. The result
 *        will be stored here.
 * \param n The number of vertices in the graph; a nonnegative integer.
 * \param directed Logical constant; gives whether the graph should be
 *        directed. Supported values are:
 *        \clist
 *        \cli IGRAPH_DIRECTED
 *          The graph to be created will be \em directed.
 *        \cli IGRAPH_UNDIRECTED
 *          The graph to be created will be \em undirected.
 *        \endclist
 * \param ... The additional arguments giving the edges of the
 *        graph. Don't forget to supply an additional '-1' after the last
 *        (meaningful) argument.
 * \return Error code.
 *
 * Time complexity: O(|V|+|E|), the number of vertices plus the number
 * of edges in the graph to create.
 *
 * \example examples/simple/igraph_small.c
 */

int igraph_small(igraph_t *graph, igraph_integer_t n, igraph_bool_t directed,
                 ...) {
    igraph_vector_t edges;
    va_list ap;

    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);

    va_start(ap, directed);
    while (1) {
        int num = va_arg(ap, int);
        if (num == -1) {
            break;
        }
        igraph_vector_push_back(&edges, num);
    }

    IGRAPH_CHECK(igraph_create(graph, &edges, n, directed));

    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);
    return 0;
}
