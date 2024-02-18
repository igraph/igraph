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
 *        to the highest vertex ID in the \p edges vector it
 *        will be increased automatically. So it is safe to give 0
 *        here.
 * \param directed Boolean, whether to create a directed graph or
 *        not. If yes, then the first edge points from the first
 *        vertex ID in \p edges to the second, etc.
 * \return Error code:
 *         \c IGRAPH_EINVEVECTOR: invalid edges
 *         vector (odd number of vertices).
 *         \c IGRAPH_EINVVID: invalid (negative)
 *         vertex ID.
 *
 * Time complexity: O(|V|+|E|),
 * |V| is the number of vertices,
 * |E| the number of edges in the
 * graph.
 *
 * \example examples/simple/igraph_create.c
 */
igraph_error_t igraph_create(igraph_t *graph, const igraph_vector_int_t *edges,
                  igraph_integer_t n, igraph_bool_t directed) {
    igraph_bool_t has_edges = igraph_vector_int_size(edges) > 0;
    igraph_integer_t max;

    if (igraph_vector_int_size(edges) % 2 != 0) {
        IGRAPH_ERROR("Invalid (odd) edges vector.", IGRAPH_EINVEVECTOR);
    }
    if (has_edges && !igraph_vector_int_isininterval(edges, 0, IGRAPH_VCOUNT_MAX-1)) {
        IGRAPH_ERROR("Invalid (negative or too large) vertex ID.", IGRAPH_EINVVID);
    }

    /* The + 1 here cannot overflow as above we have already
     * checked that vertex IDs are within range. */
    max = has_edges ? igraph_vector_int_max(edges) + 1 : 0;

    IGRAPH_CHECK(igraph_empty(graph, n, directed));
    IGRAPH_FINALLY(igraph_destroy, graph);
    if (has_edges) {
        n = igraph_vcount(graph);
        if (n < max) {
            IGRAPH_CHECK(igraph_add_vertices(graph, (max - n), 0));
        }
        IGRAPH_CHECK(igraph_add_edges(graph, edges, 0));
    }

    IGRAPH_FINALLY_CLEAN(1);
    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_small
 * \brief Shorthand to create a small graph, giving the edges as arguments.
 *
 * This function is handy when a relatively small graph needs to be created.
 * Instead of giving the edges as a vector, they are given simply as
 * arguments and a <code>-1</code> needs to be given after the last meaningful
 * edge argument.
 *
 * </para><para>
 * This function is intended to be used with vertex IDs that are entered as
 * literal integers. If you use a variable instead of a literal, make sure
 * that it is of type <type>int</type>, as this is the type that this function
 * assumes for all variadic arguments. Using a different integer type is
 * undefined behaviour and likely to cause platform-specific issues.
 *
 * \param graph Pointer to an uninitialized graph object. The result
 *        will be stored here.
 * \param n The number of vertices in the graph; a non-negative integer.
 * \param directed Logical constant; gives whether the graph should be
 *        directed. Supported values are:
 *        \clist
 *        \cli IGRAPH_DIRECTED
 *          The graph to be created will be \em directed.
 *        \cli IGRAPH_UNDIRECTED
 *          The graph to be created will be \em undirected.
 *        \endclist
 * \param ... The additional arguments giving the edges of the graph,
 *        and \em must be of type <type>int</type>. Don't forget to supply an
 *        additional <code>-1</code> after the last (meaningful) argument. The
 *        \p first parameter is present for technical reasons and represents
 *        the first variadic argument.
 * \return Error code.
 *
 * Time complexity: O(|V|+|E|), the number of vertices plus the number
 * of edges in the graph to create.
 *
 * \example examples/simple/igraph_small.c
 */

igraph_error_t igraph_small(igraph_t *graph, igraph_integer_t n, igraph_bool_t directed,
                            int first, ...) {
    igraph_vector_int_t edges;
    va_list ap;

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);

    va_start(ap, first);
    int num = first;
    while (num != -1) {
        igraph_vector_int_push_back(&edges, num);
        num = va_arg(ap, int);
    }
    va_end(ap);

    IGRAPH_CHECK(igraph_create(graph, &edges, n, directed));

    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);
    return IGRAPH_SUCCESS;
}
