/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2005-2021  The igraph development team

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

#include "igraph_datatype.h"
#include "igraph_interface.h"

/* Internal functions */

/* The functions in this file are sensible "default" implementations for some
 * of the core API functions that simply call other core API functions. If
 * you are implementing your own data type, chances are that you can use these
 * as is. */

/**
 * \ingroup interface
 * \function igraph_empty
 * \brief Creates an empty graph with some vertices and no edges.
 *
 * </para><para>
 * The most basic constructor, all the other constructors should call
 * this to create a minimal graph object. Our use of the term "empty graph"
 * in the above description should be distinguished from the mathematical
 * definition of the empty or null graph. Strictly speaking, the empty or null
 * graph in graph theory is the graph with no vertices and no edges. However
 * by "empty graph" as used in \c igraph we mean a graph having zero or more
 * vertices, but no edges.
 * \param graph Pointer to a not-yet initialized graph object.
 * \param n The number of vertices in the graph, a non-negative
 *          integer number is expected.
 * \param directed Boolean; whether the graph is directed or not. Supported
 *        values are:
 *        \clist
 *        \cli IGRAPH_DIRECTED
 *          The graph will be \em directed.
 *        \cli IGRAPH_UNDIRECTED
 *          The graph will be \em undirected.
 *        \endclist
 * \return Error code:
 *         \c IGRAPH_EINVAL: invalid number of vertices.
 *
 * Time complexity: O(|V|) for a graph with
 * |V| vertices (and no edges).
 *
 * \example examples/simple/creation.c
 */
igraph_error_t igraph_empty(igraph_t *graph, igraph_integer_t n, igraph_bool_t directed) {
    return igraph_empty_attrs(graph, n, directed, 0);
}

/**
 * \ingroup interface
 * \function igraph_delete_vertices
 * \brief Removes some vertices (with all their edges) from the graph.
 *
 * </para><para>
 * This function changes the IDs of the vertices (except in some very
 * special cases, but these should not be relied on anyway).
 *
 * </para><para>
 * This function invalidates all iterators.
 *
 * \param graph The graph to work on.
 * \param vertices The IDs of the vertices to remove, in a vector. The vector
 *     may contain the same ID more than once.
 * \return Error code:
 *         \c IGRAPH_EINVVID: invalid vertex ID.
 *
 * Time complexity: O(|V|+|E|), |V| and |E| are the number of vertices and
 * edges in the original graph.
 *
 * \example examples/simple/igraph_delete_vertices.c
 */
igraph_error_t igraph_delete_vertices(igraph_t *graph, const igraph_vs_t vertices) {
    return igraph_delete_vertices_idx(graph, vertices, /* idx= */ 0, /* invidx= */ 0);
}

/**
 * \function igraph_edge
 * \brief Returns the head and tail vertices of an edge.
 *
 * \param graph The graph object.
 * \param eid The edge ID.
 * \param from Pointer to an \type igraph_integer_t. The tail (source) of
 * the edge will be placed here.
 * \param to Pointer to an \type igraph_integer_t. The head (target) of the
 * edge will be placed here.
 * \return Error code.
 *
 * \sa \ref igraph_get_eid() for the opposite operation;
 *     \ref igraph_edges() to get the endpoints of several edges;
 *     \ref IGRAPH_TO(), \ref IGRAPH_FROM() and \ref IGRAPH_OTHER() for
 *     a faster but non-error-checked version.
 *
 * Added in version 0.2.</para><para>
 *
 * Time complexity: O(1).
 */
igraph_error_t igraph_edge(
    const igraph_t *graph, igraph_integer_t eid,
    igraph_integer_t *from, igraph_integer_t *to
) {

    if (eid < 0 || eid >= igraph_ecount(graph)) {
        IGRAPH_ERROR("Invalid edge ID when retrieving edge endpoints.", IGRAPH_EINVAL);
    }

    if (igraph_is_directed(graph)) {
        *from = IGRAPH_FROM(graph, eid);
        *to   = IGRAPH_TO(graph, eid);
    } else {
        *from = IGRAPH_TO(graph, eid);
        *to   = IGRAPH_FROM(graph, eid);
    }

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_edges
 * \brief Gives the head and tail vertices of a series of edges.
 *
 * \param graph The graph object.
 * \param eids  Edge selector, the series of edges.
 * \param edges Pointer to an initialized vector. The start and endpoints of
 *              each edge will be placed here.
 * \return Error code.
 * \sa \ref igraph_get_edgelist() to get the endpoints of all edges;
 *     \ref igraph_get_eids() for the opposite operation;
 *     \ref igraph_edge() for getting the endpoints of a single edge;
 *     \ref IGRAPH_TO(), \ref IGRAPH_FROM() and \ref IGRAPH_OTHER() for
 *     a faster but non-error-checked method.
 *
 * Time complexity: O(k) where k is the number of edges in the selector.
 */
igraph_error_t igraph_edges(const igraph_t *graph, igraph_es_t eids, igraph_vector_int_t *edges) {
    igraph_eit_t eit;
    igraph_integer_t n, ptr = 0;

    IGRAPH_CHECK(igraph_eit_create(graph, eids, &eit));
    IGRAPH_FINALLY(igraph_eit_destroy, &eit);
    n = IGRAPH_EIT_SIZE(eit);
    IGRAPH_CHECK(igraph_vector_int_resize(edges, n * 2));

    if (igraph_is_directed(graph)) {
        for (; !IGRAPH_EIT_END(eit); IGRAPH_EIT_NEXT(eit)) {
            igraph_integer_t e = IGRAPH_EIT_GET(eit);
            VECTOR(*edges)[ptr++] = IGRAPH_FROM(graph, e);
            VECTOR(*edges)[ptr++] = IGRAPH_TO(graph, e);
        }
    } else {
        for (; !IGRAPH_EIT_END(eit); IGRAPH_EIT_NEXT(eit)) {
            igraph_integer_t e = IGRAPH_EIT_GET(eit);
            VECTOR(*edges)[ptr++] = IGRAPH_TO(graph, e);
            VECTOR(*edges)[ptr++] = IGRAPH_FROM(graph, e);
        }
    }

    igraph_eit_destroy(&eit);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_invalidate_cache
 * \brief Invalidates the internal cache of an igraph graph.
 *
 * </para><para>
 * igraph graphs cache some basic properties about themselves in an internal
 * data structure. This function invalidates the contents of the cache and
 * forces a recalculation of the cached properties the next time they are
 * needed.
 *
 * </para><para>
 * You should not need to call this function during normal usage; however, we
 * might ask you to call this function explicitly if we suspect that you are
 * running into a bug in igraph's cache handling. A tell-tale sign of an invalid
 * cache entry is that the result of a cached igraph function (such as
 * \ref igraph_is_dag() or \ref igraph_is_simple()) is different before and
 * after a cache invalidation.
 *
 * \param graph The graph whose cache is to be invalidated.
 *
 * Time complexity: O(1).
 */
void igraph_invalidate_cache(const igraph_t* graph) {
    igraph_i_property_cache_invalidate_all(graph);
}
