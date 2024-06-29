/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2005-2021  The igraph development team
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

#include "igraph_datatype.h"
#include "igraph_interface.h"
#include "igraph_memory.h"

#include "graph/attributes.h"
#include "graph/caching.h"
#include "graph/internal.h"
#include "math/safe_intop.h"

/* Internal functions */

static igraph_error_t igraph_i_create_start_vectors(
        igraph_vector_int_t *res, igraph_vector_int_t *el,
        igraph_vector_int_t *index, igraph_integer_t nodes);

/**
 * \section about_basic_interface
 *
 * <para>This is the very minimal API in \a igraph. All the other
 * functions use this minimal set for creating and manipulating
 * graphs.</para>
 *
 * <para>This is a very important principle since it makes possible to
 * implement other data representations by implementing only this
 * minimal set.</para>
 *
 * <para>This section lists all the functions and macros that are considered
 * as part of the core API from the point of view of the \em users
 * of igraph. Some of these functions and macros have sensible default
 * implementations that simply call some other core function (e.g.,
 * \ref igraph_empty() calls \ref igraph_empty_attrs() with a null attribute
 * table pointer). If you wish to experiment with implementing an alternative
 * data type, the actual number of functions that you need to replace is lower
 * as you can rely on the same default implementations in most cases.</para>
 */

/**
 * \ingroup interface
 * \function igraph_empty_attrs
 * \brief Creates an empty graph with some vertices, no edges and some graph attributes.
 *
 * Use this instead of \ref igraph_empty() if you wish to add some graph
 * attributes right after initialization. This function is currently
 * not very interesting for the ordinary user. Just supply 0 here or
 * use \ref igraph_empty().
 *
 * </para><para>
 * This function does not set any vertex attributes. To create a graph which has
 * vertex attributes, call this function specifying 0 vertices, then use
 * \ref igraph_add_vertices() to add vertices and their attributes.
 *
 * \param graph Pointer to a not-yet initialized graph object.
 * \param n The number of vertices in the graph; a non-negative
 *          integer number is expected.
 * \param directed Boolean; whether the graph is directed or not. Supported
 *        values are:
 *        \clist
 *        \cli IGRAPH_DIRECTED
 *          Create a \em directed graph.
 *        \cli IGRAPH_UNDIRECTED
 *          Create an \em undirected graph.
 *        \endclist
 * \param attr The graph attributes. Supply \c NULL if not graph attributes
 *        are to be set.
 * \return Error code:
 *         \c IGRAPH_EINVAL: invalid number of vertices.
 *
 * \sa \ref igraph_empty() to create an empty graph without attributes;
 * \ref igraph_add_vertices() and \ref igraph_add_edges() to add vertices
 * and edges, possibly with associated attributes.
 *
 * Time complexity: O(|V|) for a graph with
 * |V| vertices (and no edges).
 */
igraph_error_t igraph_empty_attrs(igraph_t *graph, igraph_integer_t n, igraph_bool_t directed, void *attr) {

    if (n < 0) {
        IGRAPH_ERROR("Number of vertices must not be negative.", IGRAPH_EINVAL);
    }

    graph->n = 0;
    graph->directed = directed;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&graph->from, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&graph->to, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&graph->oi, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&graph->ii, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&graph->os, 1);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&graph->is, 1);

    /* init cache */
    graph->cache = IGRAPH_CALLOC(1, igraph_i_property_cache_t);
    IGRAPH_CHECK_OOM(graph->cache, "Cannot create graph.");
    IGRAPH_FINALLY(igraph_free, graph->cache);
    IGRAPH_CHECK(igraph_i_property_cache_init(graph->cache));
    IGRAPH_FINALLY(igraph_i_property_cache_destroy, graph->cache);

    VECTOR(graph->os)[0] = 0;
    VECTOR(graph->is)[0] = 0;

    /* init attributes */
    graph->attr = 0;
    IGRAPH_CHECK(igraph_i_attribute_init(graph, attr));

    /* add the vertices */
    IGRAPH_CHECK(igraph_add_vertices(graph, n, 0));

    IGRAPH_FINALLY_CLEAN(8);
    return IGRAPH_SUCCESS;
}

/**
 * \ingroup interface
 * \function igraph_destroy
 * \brief Frees the memory allocated for a graph object.
 *
 * </para><para>
 * This function should be called for every graph object exactly once.
 *
 * </para><para>
 * This function invalidates all iterators (of course), but the
 * iterators of a graph should be destroyed before the graph itself
 * anyway.
 * \param graph Pointer to the graph to free.
 *
 * Time complexity: operating system specific.
 */
void igraph_destroy(igraph_t *graph) {

    IGRAPH_I_ATTRIBUTE_DESTROY(graph);

    igraph_i_property_cache_destroy(graph->cache);
    IGRAPH_FREE(graph->cache);

    igraph_vector_int_destroy(&graph->from);
    igraph_vector_int_destroy(&graph->to);
    igraph_vector_int_destroy(&graph->oi);
    igraph_vector_int_destroy(&graph->ii);
    igraph_vector_int_destroy(&graph->os);
    igraph_vector_int_destroy(&graph->is);
}

/**
 * \ingroup interface
 * \function igraph_copy
 * \brief Creates an exact (deep) copy of a graph.
 *
 * </para><para>
 * This function deeply copies a graph object to create an exact
 * replica of it. The new replica should be destroyed by calling
 * \ref igraph_destroy() on it when not needed any more.
 *
 * </para><para>
 * You can also create a shallow copy of a graph by simply using the
 * standard assignment operator, but be careful and do \em not
 * destroy a shallow replica. To avoid this mistake, creating shallow
 * copies is not recommended.
 * \param to Pointer to an uninitialized graph object.
 * \param from Pointer to the graph object to copy.
 * \return Error code.
 *
 * Time complexity:  O(|V|+|E|) for a
 * graph with |V| vertices and
 * |E| edges.
 *
 * \example examples/simple/igraph_copy.c
 */

igraph_error_t igraph_copy(igraph_t *to, const igraph_t *from) {
    to->n = from->n;
    to->directed = from->directed;
    IGRAPH_CHECK(igraph_vector_int_init_copy(&to->from, &from->from));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &to->from);
    IGRAPH_CHECK(igraph_vector_int_init_copy(&to->to, &from->to));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &to->to);
    IGRAPH_CHECK(igraph_vector_int_init_copy(&to->oi, &from->oi));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &to->oi);
    IGRAPH_CHECK(igraph_vector_int_init_copy(&to->ii, &from->ii));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &to->ii);
    IGRAPH_CHECK(igraph_vector_int_init_copy(&to->os, &from->os));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &to->os);
    IGRAPH_CHECK(igraph_vector_int_init_copy(&to->is, &from->is));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &to->is);

    to->cache = IGRAPH_CALLOC(1, igraph_i_property_cache_t);
    IGRAPH_CHECK_OOM(to->cache, "Cannot copy graph.");
    IGRAPH_FINALLY(igraph_free, to->cache);
    IGRAPH_CHECK(igraph_i_property_cache_copy(to->cache, from->cache));
    IGRAPH_FINALLY(igraph_i_property_cache_destroy, to->cache);

    IGRAPH_I_ATTRIBUTE_COPY(to, from, true, true, true); /* does IGRAPH_CHECK */

    IGRAPH_FINALLY_CLEAN(8);
    return IGRAPH_SUCCESS;
}

/**
 * \ingroup interface
 * \function igraph_add_edges
 * \brief Adds edges to a graph object.
 *
 * </para><para>
 * The edges are given in a vector, the
 * first two elements define the first edge (the order is
 * <code>from</code>, <code>to</code> for directed
 * graphs). The vector
 * should contain even number of integer numbers between zero and the
 * number of vertices in the graph minus one (inclusive). If you also
 * want to add new vertices, call \ref igraph_add_vertices() first.
 * \param graph The graph to which the edges will be added.
 * \param edges The edges themselves.
 * \param attr The attributes of the new edges. You can supply a null pointer
 *        here if you do not need edge attributes.
 * \return Error code:
 *    \c IGRAPH_EINVEVECTOR: invalid (odd) edges vector length,
 *    \c IGRAPH_EINVVID: invalid vertex ID in edges vector.
 *
 * This function invalidates all iterators.
 *
 * </para><para>
 * Time complexity: O(|V|+|E|) where |V| is the number of vertices and
 * |E| is the number of edges in the \em new, extended graph.
 *
 * \example examples/simple/creation.c
 */
igraph_error_t igraph_add_edges(igraph_t *graph, const igraph_vector_int_t *edges,
                     void *attr) {
    igraph_integer_t no_of_edges = igraph_vector_int_size(&graph->from);
    igraph_integer_t edges_to_add = igraph_vector_int_size(edges) / 2;
    igraph_integer_t new_no_of_edges;
    igraph_integer_t i = 0;
    igraph_vector_int_t newoi, newii;
    igraph_bool_t directed = igraph_is_directed(graph);

    if (igraph_vector_int_size(edges) % 2 != 0) {
        IGRAPH_ERROR("Invalid (odd) length of edges vector.", IGRAPH_EINVEVECTOR);
    }
    if (!igraph_vector_int_isininterval(edges, 0, igraph_vcount(graph) - 1)) {
        IGRAPH_ERROR("Out-of-range vertex IDs when adding edges.", IGRAPH_EINVVID);
    }

    /* from & to */
    IGRAPH_SAFE_ADD(no_of_edges, edges_to_add, &new_no_of_edges);
    if (new_no_of_edges > IGRAPH_ECOUNT_MAX) {
        IGRAPH_ERRORF("Maximum edge count (%" IGRAPH_PRId ") exceeded.", IGRAPH_ERANGE,
                      IGRAPH_ECOUNT_MAX);
    }
    IGRAPH_CHECK(igraph_vector_int_reserve(&graph->from, no_of_edges + edges_to_add));
    IGRAPH_CHECK(igraph_vector_int_reserve(&graph->to, no_of_edges + edges_to_add));

    while (i < edges_to_add * 2) {
        if (directed || VECTOR(*edges)[i] > VECTOR(*edges)[i + 1]) {
            igraph_vector_int_push_back(&graph->from, VECTOR(*edges)[i++]); /* reserved */
            igraph_vector_int_push_back(&graph->to,   VECTOR(*edges)[i++]); /* reserved */
        } else {
            igraph_vector_int_push_back(&graph->to,   VECTOR(*edges)[i++]); /* reserved */
            igraph_vector_int_push_back(&graph->from, VECTOR(*edges)[i++]); /* reserved */
        }
    }

    /* If an error occurs while the edges are being added, we make the necessary fixup
     * to ensure that the graph is still in a consistent state when this function returns.
     * The graph may already be on the finally stack when calling this function. We use
     * a separate finally stack level to avoid its destructor from being called on error,
     * so that the fixup can succeed.
     */

#define CHECK_ERR(expr) \
    do { \
        igraph_error_t err = (expr); \
        if (err != IGRAPH_SUCCESS) { \
            igraph_vector_int_resize(&graph->from, no_of_edges); /* gets smaller, error safe */ \
            igraph_vector_int_resize(&graph->to, no_of_edges);   /* gets smaller, error safe */ \
            IGRAPH_FINALLY_EXIT(); \
            IGRAPH_ERROR("Cannot add edges.", err); \
        } \
    } while (0)

    /* oi & ii */
    IGRAPH_FINALLY_ENTER();
    {
        CHECK_ERR(igraph_vector_int_init(&newoi, no_of_edges));
        IGRAPH_FINALLY(igraph_vector_int_destroy, &newoi);
        CHECK_ERR(igraph_vector_int_init(&newii, no_of_edges));
        IGRAPH_FINALLY(igraph_vector_int_destroy, &newii);
        CHECK_ERR(igraph_vector_int_pair_order(&graph->from, &graph->to, &newoi, graph->n));
        CHECK_ERR(igraph_vector_int_pair_order(&graph->to, &graph->from, &newii, graph->n));

        /* Attributes */
        if (graph->attr) {
            /* TODO: Does this keep the attribute table in a consistent state upon failure? */
            CHECK_ERR(igraph_i_attribute_add_edges(graph, edges, attr));
        }

        /* os & is, its length does not change, error safe */
        igraph_i_create_start_vectors(&graph->os, &graph->from, &newoi, graph->n);
        igraph_i_create_start_vectors(&graph->is, &graph->to, &newii, graph->n);

        /* everything went fine */
        igraph_vector_int_destroy(&graph->oi);
        igraph_vector_int_destroy(&graph->ii);
        IGRAPH_FINALLY_CLEAN(2);

        graph->oi = newoi;
        graph->ii = newii;
    }
    IGRAPH_FINALLY_EXIT();

#undef CHECK_ERR

    /* modification successful, clear the cached properties of the graph.
     *
     * Adding one or more edges cannot make a strongly or weakly connected
     * graph disconnected, so we keep those flags if they are cached as true.
     *
     * Adding one or more edges may turn a DAG into a non-DAG or a forest into
     * a non-forest, so we can keep those flags only if they are cached as
     * false.
     *
     * Also, adding one or more edges does not change HAS_LOOP, HAS_MULTI and
     * HAS_MUTUAL if they were already true.
     */
    igraph_i_property_cache_invalidate_conditionally(
        graph,
        /* keep_always = */ 0,
        /* keep_when_false = */
        (1 << IGRAPH_PROP_IS_DAG) | (1 << IGRAPH_PROP_IS_FOREST),
        /* keep_when_true = */
        (1 << IGRAPH_PROP_IS_WEAKLY_CONNECTED) |
        (1 << IGRAPH_PROP_IS_STRONGLY_CONNECTED) |
        (1 << IGRAPH_PROP_HAS_LOOP) |
        (1 << IGRAPH_PROP_HAS_MULTI) |
        (1 << IGRAPH_PROP_HAS_MUTUAL)
    );

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup interface
 * \function igraph_add_vertices
 * \brief Adds vertices to a graph.
 *
 * </para><para>
 * This function invalidates all iterators.
 *
 * \param graph The graph object to extend.
 * \param nv Non-negative integer specifying the number of vertices to add.
 * \param attr The attributes of the new vertices. You can supply a null pointer
 *        here if you do not need vertex attributes.
 * \return Error code:
 *         \c IGRAPH_EINVAL: invalid number of new vertices.
 *
 * Time complexity: O(|V|) where |V| is the number of vertices in the \em new,
 * extended graph.
 *
 * \example examples/simple/creation.c
 */
igraph_error_t igraph_add_vertices(igraph_t *graph, igraph_integer_t nv, void *attr) {
    igraph_integer_t ec = igraph_ecount(graph);
    igraph_integer_t vc = igraph_vcount(graph);
    igraph_integer_t new_vc;
    igraph_integer_t i;

    if (nv < 0) {
        IGRAPH_ERROR("Cannot add negative number of vertices.", IGRAPH_EINVAL);
    }

    IGRAPH_SAFE_ADD(graph->n, nv, &new_vc);
    if (new_vc > IGRAPH_VCOUNT_MAX) {
        IGRAPH_ERRORF("Maximum vertex count (%" IGRAPH_PRId ") exceeded.", IGRAPH_ERANGE,
                      IGRAPH_VCOUNT_MAX);
    }
    IGRAPH_CHECK(igraph_vector_int_reserve(&graph->os, new_vc + 1));
    IGRAPH_CHECK(igraph_vector_int_reserve(&graph->is, new_vc + 1));

    igraph_vector_int_resize(&graph->os, new_vc + 1); /* reserved */
    igraph_vector_int_resize(&graph->is, new_vc + 1); /* reserved */
    for (i = graph->n + 1; i < new_vc + 1; i++) {
        VECTOR(graph->os)[i] = ec;
        VECTOR(graph->is)[i] = ec;
    }

    graph->n += nv;

    /* Add attributes if necessary. This section is protected with
     * FINALLY_ENTER/EXIT so that the graph would not be accidentally
     * free upon error until it could be restored to a consistant state. */

    if (graph->attr) {
        igraph_error_t err;
        IGRAPH_FINALLY_ENTER();
        err = igraph_i_attribute_add_vertices(graph, nv, attr);
        if (err != IGRAPH_SUCCESS) {
            /* Restore original vertex count on failure */
            graph->n = vc;
            igraph_vector_int_resize(&graph->os, vc + 1); /* shrinks */
            igraph_vector_int_resize(&graph->is, vc + 1); /* shrinks */
        }
        IGRAPH_FINALLY_EXIT();
        if (err != IGRAPH_SUCCESS) {
            IGRAPH_ERROR("Cannot add vertices.", err);
        }
    }

    /* modification successful, clear the cached properties of the graph.
     *
     * Adding one or more nodes does not change the following cached properties:
     *
     * - IGRAPH_PROP_HAS_LOOP
     * - IGRAPH_PROP_HAS_MULTI
     * - IGRAPH_PROP_HAS_MUTUAL
     * - IGRAPH_PROP_IS_DAG (adding a node does not create/destroy cycles)
     * - IGRAPH_PROP_IS_FOREST (same)
     *
     * Adding one or more nodes without any edges incident on them is sure to
     * make the graph disconnected (weakly or strongly), so we can keep the
     * connectivity-related properties if they are currently cached as false.
     * (Actually, even if they weren't cached as false, we could still set them
     * to false, but we don't have that functionality yet). The only exception
     * is when the graph had zero vertices and gained only one vertex, because
     * it then becomes connected. That's why we have the condition below in the
     * keep_when_false section.
     */
    igraph_i_property_cache_invalidate_conditionally(
        graph,
        /* keep_always = */
        (1 << IGRAPH_PROP_HAS_LOOP) |
        (1 << IGRAPH_PROP_HAS_MULTI) |
        (1 << IGRAPH_PROP_HAS_MUTUAL) |
        (1 << IGRAPH_PROP_IS_DAG) |
        (1 << IGRAPH_PROP_IS_FOREST),
        /* keep_when_false = */
        igraph_vcount(graph) >= 2 ? (
            (1 << IGRAPH_PROP_IS_STRONGLY_CONNECTED) |
            (1 << IGRAPH_PROP_IS_WEAKLY_CONNECTED)
        ) : 0,
        /* keep_when_true = */
        0
    );

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup interface
 * \function igraph_delete_edges
 * \brief Removes edges from a graph.
 *
 * </para><para>
 * The edges to remove are specified as an edge selector.
 *
 * </para><para>
 * This function cannot remove vertices; vertices will be kept even if they lose
 * all their edges.
 *
 * </para><para>
 * This function invalidates all iterators.
 * \param graph The graph to work on.
 * \param edges The edges to remove.
 * \return Error code.
 *
 * Time complexity: O(|V|+|E|) where |V| and |E| are the number of vertices
 * and edges in the \em original graph, respectively.
 *
 * \example examples/simple/igraph_delete_edges.c
 */
igraph_error_t igraph_delete_edges(igraph_t *graph, igraph_es_t edges) {
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t edges_to_remove = 0;
    igraph_integer_t remaining_edges;
    igraph_eit_t eit;

    igraph_vector_int_t newfrom, newto;
    igraph_vector_int_t newoi, newii;

    igraph_bool_t *mark;
    igraph_integer_t i, j;

    mark = IGRAPH_CALLOC(no_of_edges, igraph_bool_t);
    IGRAPH_CHECK_OOM(mark, "Cannot delete edges.");
    IGRAPH_FINALLY(igraph_free, mark);

    IGRAPH_CHECK(igraph_eit_create(graph, edges, &eit));
    IGRAPH_FINALLY(igraph_eit_destroy, &eit);

    for (IGRAPH_EIT_RESET(eit); !IGRAPH_EIT_END(eit); IGRAPH_EIT_NEXT(eit)) {
        igraph_integer_t e = IGRAPH_EIT_GET(eit);
        if (! mark[e]) {
            edges_to_remove++;
            mark[e] = true;
        }
    }
    remaining_edges = no_of_edges - edges_to_remove;

    /* We don't need the iterator any more */
    igraph_eit_destroy(&eit);
    IGRAPH_FINALLY_CLEAN(1);

    IGRAPH_VECTOR_INT_INIT_FINALLY(&newfrom, remaining_edges);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&newto, remaining_edges);

    /* Actually remove the edges, move from pos i to pos j in newfrom/newto */
    for (i = 0, j = 0; j < remaining_edges; i++) {
        if (! mark[i]) {
            VECTOR(newfrom)[j] = VECTOR(graph->from)[i];
            VECTOR(newto)[j] = VECTOR(graph->to)[i];
            j++;
        }
    }

    /* Create index, this might require additional memory */
    IGRAPH_VECTOR_INT_INIT_FINALLY(&newoi, remaining_edges);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&newii, remaining_edges);
    IGRAPH_CHECK(igraph_vector_int_pair_order(&newfrom, &newto, &newoi, no_of_nodes));
    IGRAPH_CHECK(igraph_vector_int_pair_order(&newto, &newfrom, &newii, no_of_nodes));

    /* Edge attributes, we need an index that gives the IDs of the
       original edges for every new edge.
    */
    if (graph->attr) {
        igraph_vector_int_t idx;
        IGRAPH_VECTOR_INT_INIT_FINALLY(&idx, remaining_edges);
        for (i = 0, j = 0; i < no_of_edges; i++) {
            if (! mark[i]) {
                VECTOR(idx)[j++] = i;
            }
        }
        IGRAPH_CHECK(igraph_i_attribute_permute_edges(graph, graph, &idx));
        igraph_vector_int_destroy(&idx);
        IGRAPH_FINALLY_CLEAN(1);
    }

    /* Ok, we've all memory needed, free the old structure  */
    igraph_vector_int_destroy(&graph->from);
    igraph_vector_int_destroy(&graph->to);
    igraph_vector_int_destroy(&graph->oi);
    igraph_vector_int_destroy(&graph->ii);
    graph->from = newfrom;
    graph->to = newto;
    graph->oi = newoi;
    graph->ii = newii;
    IGRAPH_FINALLY_CLEAN(4);

    IGRAPH_FREE(mark);
    IGRAPH_FINALLY_CLEAN(1);

    /* Create start vectors, no memory is needed for this */
    igraph_i_create_start_vectors(&graph->os, &graph->from, &graph->oi, no_of_nodes);
    igraph_i_create_start_vectors(&graph->is, &graph->to,   &graph->ii, no_of_nodes);

    /* modification successful, clear the cached properties of the graph.
     *
     * Deleting one or more edges cannot make a directed acyclic graph cyclic,
     * or an undirected forest into a cyclic graph, so we keep those flags if
     * they are cached as true.
     *
     * Similarly, deleting one or more edges cannot make a disconnected graph
     * connected, so we keep the connectivity flags if they are cached as false.
     *
     * Also, if the graph had no loop edges before the deletion, it will have
     * no loop edges after the deletion either. The same applies to reciprocal
     * edges or multiple edges as well.
     */
    igraph_i_property_cache_invalidate_conditionally(
        graph,
        /* keep_always = */ 0,
        /* keep_when_false = */
        (1 << IGRAPH_PROP_HAS_LOOP) |
        (1 << IGRAPH_PROP_HAS_MULTI) |
        (1 << IGRAPH_PROP_HAS_MUTUAL) |
        (1 << IGRAPH_PROP_IS_STRONGLY_CONNECTED) |
        (1 << IGRAPH_PROP_IS_WEAKLY_CONNECTED),
        /* keep_when_true = */
        (1 << IGRAPH_PROP_IS_DAG) |
        (1 << IGRAPH_PROP_IS_FOREST)
    );

    /* Nothing to deallocate... */
    return IGRAPH_SUCCESS;
}

/**
 * \ingroup interface
 * \function igraph_delete_vertices_idx
 * \brief Removes some vertices (with all their edges) from the graph.
 *
 * </para><para>
 * This function changes the IDs of the vertices (except in some very
 * special cases, but these should not be relied on anyway). You can use the
 * \c idx argument to obtain the mapping from old vertex IDs to the new ones,
 * and the \c newidx argument to obtain the reverse mapping.
 *
 * </para><para>
 * This function invalidates all iterators.
 *
 * \param graph The graph to work on.
 * \param vertices The IDs of the vertices to remove, in a vector. The vector
 *     may contain the same ID more than once.
 * \param idx An optional pointer to a vector that provides the mapping from
 *     the vertex IDs \em before the removal to the vertex IDs \em after
 *     the removal, \em plus one. Zero is used to represent vertices that were
 *     removed during the operation. You can supply \c NULL here if you are not
 *     interested.
 * \param invidx An optional pointer to a vector that provides the mapping from
 *     the vertex IDs \em after the removal to the vertex IDs \em before
 *     the removal. You can supply \c NULL here if you are not interested.
 * \return Error code:
 *         \c IGRAPH_EINVVID: invalid vertex ID.
 *
 * Time complexity: O(|V|+|E|), |V| and |E| are the number of vertices and
 * edges in the original graph.
 *
 * \example examples/simple/igraph_delete_vertices.c
 */
igraph_error_t igraph_delete_vertices_idx(
    igraph_t *graph, const igraph_vs_t vertices, igraph_vector_int_t *idx,
    igraph_vector_int_t *invidx
) {
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_vector_int_t edge_recoding, vertex_recoding;
    igraph_vector_int_t *my_vertex_recoding = &vertex_recoding;
    igraph_vit_t vit;
    igraph_t newgraph;
    igraph_integer_t i, j;
    igraph_integer_t remaining_vertices, remaining_edges;

    if (idx) {
        my_vertex_recoding = idx;
        IGRAPH_CHECK(igraph_vector_int_resize(idx, no_of_nodes));
        igraph_vector_int_null(idx);
    } else {
        IGRAPH_VECTOR_INT_INIT_FINALLY(&vertex_recoding, no_of_nodes);
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edge_recoding, no_of_edges);

    IGRAPH_CHECK(igraph_vit_create(graph, vertices, &vit));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit);

    /* mark the vertices to delete */
    for (; !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit) ) {
        igraph_integer_t vertex = IGRAPH_VIT_GET(vit);
        if (vertex < 0 || vertex >= no_of_nodes) {
            IGRAPH_ERROR("Cannot delete vertices", IGRAPH_EINVVID);
        }
        VECTOR(*my_vertex_recoding)[vertex] = 1;
    }
    /* create vertex recoding vector */
    for (remaining_vertices = 0, i = 0; i < no_of_nodes; i++) {
        if (VECTOR(*my_vertex_recoding)[i] == 0) {
            VECTOR(*my_vertex_recoding)[i] = remaining_vertices + 1;
            remaining_vertices++;
        } else {
            VECTOR(*my_vertex_recoding)[i] = 0;
        }
    }
    /* create edge recoding vector */
    for (remaining_edges = 0, i = 0; i < no_of_edges; i++) {
        igraph_integer_t from = VECTOR(graph->from)[i];
        igraph_integer_t to = VECTOR(graph->to)[i];
        if (VECTOR(*my_vertex_recoding)[from] != 0 &&
            VECTOR(*my_vertex_recoding)[to  ] != 0) {
            VECTOR(edge_recoding)[i] = remaining_edges + 1;
            remaining_edges++;
        }
    }

    /* start creating the graph */
    newgraph.n = remaining_vertices;
    newgraph.directed = graph->directed;

    /* allocate vectors */
    IGRAPH_VECTOR_INT_INIT_FINALLY(&newgraph.from, remaining_edges);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&newgraph.to, remaining_edges);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&newgraph.oi, remaining_edges);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&newgraph.ii, remaining_edges);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&newgraph.os, remaining_vertices + 1);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&newgraph.is, remaining_vertices + 1);

    /* Add the edges */
    for (i = 0, j = 0; j < remaining_edges; i++) {
        if (VECTOR(edge_recoding)[i] > 0) {
            igraph_integer_t from = VECTOR(graph->from)[i];
            igraph_integer_t to = VECTOR(graph->to  )[i];
            VECTOR(newgraph.from)[j] = VECTOR(*my_vertex_recoding)[from] - 1;
            VECTOR(newgraph.to  )[j] = VECTOR(*my_vertex_recoding)[to] - 1;
            j++;
        }
    }

    /* update oi & ii */
    IGRAPH_CHECK(igraph_vector_int_pair_order(&newgraph.from, &newgraph.to, &newgraph.oi,
                                         remaining_vertices));
    IGRAPH_CHECK(igraph_vector_int_pair_order(&newgraph.to, &newgraph.from, &newgraph.ii,
                                         remaining_vertices));

    IGRAPH_CHECK(igraph_i_create_start_vectors(&newgraph.os, &newgraph.from,
                                       &newgraph.oi, remaining_vertices));
    IGRAPH_CHECK(igraph_i_create_start_vectors(&newgraph.is, &newgraph.to,
                                       &newgraph.ii, remaining_vertices));

    newgraph.cache = IGRAPH_CALLOC(1, igraph_i_property_cache_t);
    IGRAPH_CHECK_OOM(newgraph.cache, "Cannot delete vertices.");
    IGRAPH_FINALLY(igraph_free, newgraph.cache);
    IGRAPH_CHECK(igraph_i_property_cache_init(newgraph.cache));
    IGRAPH_FINALLY(igraph_i_property_cache_destroy, newgraph.cache);

    /* attributes */
    IGRAPH_I_ATTRIBUTE_COPY(&newgraph, graph,
                            /*graph=*/ 1, /*vertex=*/0, /*edge=*/0);

    /* at this point igraph_destroy can take over the responsibility of
     * deallocating the graph */
    IGRAPH_FINALLY_CLEAN(8);    /* 2 for the property cache, 6 for the vectors */
    IGRAPH_FINALLY(igraph_destroy, &newgraph);

    if (newgraph.attr) {
        igraph_vector_int_t iidx;
        IGRAPH_VECTOR_INT_INIT_FINALLY(&iidx, remaining_vertices);
        for (i = 0; i < no_of_nodes; i++) {
            igraph_integer_t jj = VECTOR(*my_vertex_recoding)[i];
            if (jj != 0) {
                VECTOR(iidx)[ jj - 1 ] = i;
            }
        }
        IGRAPH_CHECK(igraph_i_attribute_permute_vertices(graph,
                     &newgraph,
                     &iidx));
        IGRAPH_CHECK(igraph_vector_int_resize(&iidx, remaining_edges));
        for (i = 0; i < no_of_edges; i++) {
            igraph_integer_t jj = VECTOR(edge_recoding)[i];
            if (jj != 0) {
                VECTOR(iidx)[ jj - 1 ] = i;
            }
        }
        IGRAPH_CHECK(igraph_i_attribute_permute_edges(graph, &newgraph, &iidx));
        igraph_vector_int_destroy(&iidx);
        IGRAPH_FINALLY_CLEAN(1);
    }

    igraph_vit_destroy(&vit);
    igraph_vector_int_destroy(&edge_recoding);
    igraph_destroy(graph);
    *graph = newgraph;

    IGRAPH_FINALLY_CLEAN(3);

    /* TODO: this is duplicate */
    if (invidx) {
        IGRAPH_CHECK(igraph_vector_int_resize(invidx, remaining_vertices));
        for (i = 0; i < no_of_nodes; i++) {
            igraph_integer_t newid = VECTOR(*my_vertex_recoding)[i];
            if (newid != 0) {
                VECTOR(*invidx)[newid - 1] = i;
            }
        }
    }

    if (!idx) {
        igraph_vector_int_destroy(my_vertex_recoding);
        IGRAPH_FINALLY_CLEAN(1);
    }

    /* modification successful, clear the cached properties of the graph.
     *
     * Deleting one or more vertices cannot make a directed acyclic graph cyclic,
     * or an undirected forest into a cyclic graph, so we keep those flags if
     * they are cached as true.
     *
     * Also, if the graph had no loop edges before the deletion, it will have
     * no loop edges after the deletion either. The same applies to reciprocal
     * edges or multiple edges as well.
     */
    igraph_i_property_cache_invalidate_conditionally(
        graph,
        /* keep_always = */ 0,
        /* keep_when_false = */
        (1 << IGRAPH_PROP_HAS_LOOP) |
        (1 << IGRAPH_PROP_HAS_MULTI) |
        (1 << IGRAPH_PROP_HAS_MUTUAL),
        /* keep_when_true = */
        (1 << IGRAPH_PROP_IS_DAG) |
        (1 << IGRAPH_PROP_IS_FOREST)
    );

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup interface
 * \function igraph_vcount
 * \brief The number of vertices in a graph.
 *
 * \param graph The graph.
 * \return Number of vertices.
 *
 * Time complexity: O(1)
 */
igraph_integer_t igraph_vcount(const igraph_t *graph) {
    return graph->n;
}

/**
 * \ingroup interface
 * \function igraph_ecount
 * \brief The number of edges in a graph.
 *
 * \param graph The graph.
 * \return Number of edges.
 *
 * Time complexity: O(1)
 */
igraph_integer_t igraph_ecount(const igraph_t *graph) {
    return igraph_vector_int_size(&graph->from);
}

/**
 * \ingroup interface
 * \function igraph_neighbors
 * \brief Adjacent vertices to a vertex.
 *
 * \param graph The graph to work on.
 * \param neis This vector will contain the result. The vector should
 *        be initialized beforehand and will be resized. Starting from igraph
 *        version 0.4 this vector is always sorted, the vertex IDs are
 *        in increasing order. If one neighbor is connected with multiple
 *        edges, the neighbor will be returned multiple times.
 * \param pnode The id of the node for which the adjacent vertices are
 *        to be searched.
 * \param mode Defines the way adjacent vertices are searched in
 *        directed graphs. It can have the following values:
 *        \c IGRAPH_OUT, vertices reachable by an
 *        edge from the specified vertex are searched;
 *        \c IGRAPH_IN, vertices from which the
 *        specified vertex is reachable are searched;
 *        \c IGRAPH_ALL, both kinds of vertices are
 *        searched.
 *        This parameter is ignored for undirected graphs.
 * \return Error code:
 *         \c IGRAPH_EINVVID: invalid vertex ID.
 *         \c IGRAPH_EINVMODE: invalid mode argument.
 *         \c IGRAPH_ENOMEM: not enough memory.
 *
 * Time complexity: O(d),
 * d is the number
 * of adjacent vertices to the queried vertex.
 *
 * \example examples/simple/igraph_neighbors.c
 */
igraph_error_t igraph_neighbors(const igraph_t *graph, igraph_vector_int_t *neis, igraph_integer_t pnode,
        igraph_neimode_t mode) {
    if (!igraph_is_directed(graph) || mode == IGRAPH_ALL) {
        return igraph_i_neighbors(graph, neis, pnode, mode, IGRAPH_LOOPS_TWICE, IGRAPH_MULTIPLE);
    } else {
        return igraph_i_neighbors(graph, neis, pnode, mode, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE);
    }
}

igraph_error_t igraph_i_neighbors(const igraph_t *graph, igraph_vector_int_t *neis, igraph_integer_t pnode,
        igraph_neimode_t mode, igraph_loops_t loops, igraph_multiple_t multiple) {
#define DEDUPLICATE_IF_NEEDED(vertex, n)                                                 \
    if (should_filter_duplicates) {                                                        \
        if (                                                                               \
            (loops == IGRAPH_NO_LOOPS && vertex == pnode) ||                               \
            (loops == IGRAPH_LOOPS_ONCE && vertex == pnode && last_added == pnode)         \
        ) {                                                                                \
            length -= n;                                                                   \
            if (loops == IGRAPH_LOOPS_ONCE) {                                              \
                last_added = -1;                                                           \
            }                                                                              \
            continue;                                                                      \
        } else if (multiple == IGRAPH_NO_MULTIPLE && vertex == last_added) {               \
            length -= n;                                                                   \
            continue;                                                                      \
        } else {                                                                           \
            last_added = vertex;                                                           \
        }                                                                                  \
    }

    igraph_integer_t length = 0, idx = 0;
    igraph_integer_t i, j;

    igraph_integer_t node = pnode;
    igraph_integer_t last_added = -1;
    igraph_bool_t should_filter_duplicates;

    if (node < 0 || node > igraph_vcount(graph) - 1) {
        IGRAPH_ERROR("Given vertex is not in the graph.", IGRAPH_EINVVID);
    }
    if (mode != IGRAPH_OUT && mode != IGRAPH_IN &&
            mode != IGRAPH_ALL) {
        IGRAPH_ERROR("Mode should be either IGRAPH_OUT, IGRAPH_IN or IGRAPH_ALL.", IGRAPH_EINVMODE);
    }

    if (!igraph_is_directed(graph)) {
        mode = IGRAPH_ALL;
    }

    if (mode != IGRAPH_ALL && loops == IGRAPH_LOOPS_TWICE) {
        IGRAPH_ERROR("For a directed graph (with directions not ignored), "
                     "IGRAPH_LOOPS_TWICE does not make sense.\n", IGRAPH_EINVAL);
    }
    /* Calculate needed space first & allocate it */
    /* Note that 'mode' is treated as a bit field here; it's okay because
     * IGRAPH_ALL = IGRAPH_IN | IGRAPH_OUT, bit-wise */
    if (mode & IGRAPH_OUT) {
        length += (VECTOR(graph->os)[node + 1] - VECTOR(graph->os)[node]);
    }
    if (mode & IGRAPH_IN) {
        length += (VECTOR(graph->is)[node + 1] - VECTOR(graph->is)[node]);
    }

    IGRAPH_CHECK(igraph_vector_int_resize(neis, length));

    /* The loops below produce an ordering what is consistent with the
     * ordering returned by igraph_neighbors(), and this should be preserved.
     * We are dealing with two sorted lists; one for the successors and one
     * for the predecessors. If we have requested only one of them, we have
     * an easy job. If we have requested both, we need to merge the two lists
     * to ensure that the output is sorted by the vertex IDs of the "other"
     * endpoint of the affected edges. We don't need to merge if the graph
     * is undirected, because in that case the data structure guarantees that
     * the "out-edges" contain only (u, v) pairs where u <= v and the
     * "in-edges" contains the rest, so the result is sorted even without
     * merging. */
    if (!igraph_is_directed(graph) || mode != IGRAPH_ALL) {
        /* graph is undirected or we did not ask for both directions in a
         * directed graph; this is the easy case */

        should_filter_duplicates = !(multiple == IGRAPH_MULTIPLE &&
                ((!igraph_is_directed(graph) && loops == IGRAPH_LOOPS_TWICE) ||
                 (igraph_is_directed(graph) && loops != IGRAPH_NO_LOOPS)));

        if (mode & IGRAPH_OUT) {
            j = VECTOR(graph->os)[node + 1];
            for (i = VECTOR(graph->os)[node]; i < j; i++) {
                igraph_integer_t to = VECTOR(graph->to)[ VECTOR(graph->oi)[i] ];
                DEDUPLICATE_IF_NEEDED(to, 1);
                VECTOR(*neis)[idx++] = to;
            }
        }

        if (mode & IGRAPH_IN) {
            j = VECTOR(graph->is)[node + 1];
            for (i = VECTOR(graph->is)[node]; i < j; i++) {
                igraph_integer_t from = VECTOR(graph->from)[ VECTOR(graph->ii)[i] ];
                DEDUPLICATE_IF_NEEDED(from, 1);
                VECTOR(*neis)[idx++] = from;
            }
        }
    } else {
        /* Both in- and out- neighbors in a directed graph,
           we need to merge the two 'vectors' so the result is
           correctly ordered. */
        igraph_integer_t j1 = VECTOR(graph->os)[node + 1];
        igraph_integer_t j2 = VECTOR(graph->is)[node + 1];
        igraph_integer_t i1 = VECTOR(graph->os)[node];
        igraph_integer_t i2 = VECTOR(graph->is)[node];
        igraph_integer_t eid1, eid2;
        igraph_integer_t n1, n2;

        should_filter_duplicates = !(multiple == IGRAPH_MULTIPLE &&
                loops == IGRAPH_LOOPS_TWICE);

        while (i1 < j1 && i2 < j2) {
            eid1 = VECTOR(graph->oi)[i1];
            eid2 = VECTOR(graph->ii)[i2];
            n1 = VECTOR(graph->to)[eid1];
            n2 = VECTOR(graph->from)[eid2];
            if (n1 < n2) {
                i1++;
                DEDUPLICATE_IF_NEEDED(n1, 1);
                VECTOR(*neis)[idx++] = n1;
            } else if (n1 > n2) {
                i2++;
                DEDUPLICATE_IF_NEEDED(n2, 1);
                VECTOR(*neis)[idx++] = n2;
            } else {
                i1++;
                i2++;
                DEDUPLICATE_IF_NEEDED(n1, 2);
                VECTOR(*neis)[idx++] = n1;
                if (should_filter_duplicates && ((loops == IGRAPH_LOOPS_ONCE && n1 == pnode && last_added == pnode) ||
                        (multiple == IGRAPH_NO_MULTIPLE))) {
                    length--;
                    if (loops == IGRAPH_LOOPS_ONCE) {
                        last_added = -1;
                    }
                    continue;
                }
                VECTOR(*neis)[idx++] = n2;
            }
        }

        while (i1 < j1) {
            eid1 = VECTOR(graph->oi)[i1++];
            igraph_integer_t to = VECTOR(graph->to)[eid1];
            DEDUPLICATE_IF_NEEDED(to, 1);
            VECTOR(*neis)[idx++] = to;
        }

        while (i2 < j2) {
            eid2 = VECTOR(graph->ii)[i2++];
            igraph_integer_t from = VECTOR(graph->from)[eid2];
            DEDUPLICATE_IF_NEEDED(from, 1);
            VECTOR(*neis)[idx++] = from;
        }

    }
    IGRAPH_CHECK(igraph_vector_int_resize(neis, length));

    return IGRAPH_SUCCESS;
#undef DEDUPLICATE_IF_NEEDED
}

/**
 * \ingroup internal
 */

static igraph_error_t igraph_i_create_start_vectors(
        igraph_vector_int_t *res, igraph_vector_int_t *el,
        igraph_vector_int_t *iindex, igraph_integer_t nodes) {

# define EDGE(i) (VECTOR(*el)[ VECTOR(*iindex)[(i)] ])

    igraph_integer_t no_of_nodes;
    igraph_integer_t no_of_edges;
    igraph_integer_t i, j, idx;

    no_of_nodes = nodes;
    no_of_edges = igraph_vector_int_size(el);

    /* result */

    IGRAPH_CHECK(igraph_vector_int_resize(res, nodes + 1));

    /* create the index */

    if (no_of_edges == 0) {
        /* empty graph */
        igraph_vector_int_null(res);
    } else {
        idx = -1;
        for (i = 0; i <= EDGE(0); i++) {
            idx++; VECTOR(*res)[idx] = 0;
        }
        for (i = 1; i < no_of_edges; i++) {
            igraph_integer_t n = EDGE(i) - EDGE(VECTOR(*res)[idx]);
            for (j = 0; j < n; j++) {
                idx++; VECTOR(*res)[idx] = i;
            }
        }
        j = EDGE(VECTOR(*res)[idx]);
        for (i = 0; i < no_of_nodes - j; i++) {
            idx++; VECTOR(*res)[idx] = no_of_edges;
        }
    }

    /* clean */

# undef EDGE
    return IGRAPH_SUCCESS;
}

/**
 * \ingroup interface
 * \function igraph_is_directed
 * \brief Is this a directed graph?
 *
 * \param graph The graph.
 * \return Logical value, \c true if the graph is directed,
 * \c false otherwise.
 *
 * Time complexity: O(1)
 *
 * \example examples/simple/igraph_is_directed.c
 */

igraph_bool_t igraph_is_directed(const igraph_t *graph) {
    return graph->directed;
}

/**
 * \ingroup interface
 * \function igraph_degree_1
 * \brief The degree of of a single vertex in the graph.
 *
 * This function calculates the in-, out- or total degree of a single vertex.
 * For a single vertex, it is more efficient than calling \ref igraph_degree().
 *
 * \param graph The graph.
 * \param deg Pointer to the integer where the computed degree will be stored.
 * \param vid The vertex for which the degree will be calculated.
 * \param mode Defines the type of the degree for directed graphs. Valid modes are:
 *        \c IGRAPH_OUT, out-degree;
 *        \c IGRAPH_IN, in-degree;
 *        \c IGRAPH_ALL, total degree (sum of the in- and out-degree).
 *        This parameter is ignored for undirected graphs.
 * \param loops Boolean, gives whether the self-loops should be
 *        counted.
 * \return Error code.
 *
 * \sa \ref igraph_degree() to compute the degree of several vertices at once.
 *
 * Time complexity: O(1) if \p loops is \c true, and
 * O(d) otherwise, where d is the degree.
 */
igraph_error_t igraph_degree_1(const igraph_t *graph, igraph_integer_t *deg,
                               igraph_integer_t vid, igraph_neimode_t mode, igraph_bool_t loops) {

    if (!igraph_is_directed(graph)) {
        mode = IGRAPH_ALL;
    }

    *deg = 0;
    if (mode & IGRAPH_OUT) {
        *deg += (VECTOR(graph->os)[vid + 1] - VECTOR(graph->os)[vid]);
    }
    if (mode & IGRAPH_IN) {
        *deg += (VECTOR(graph->is)[vid + 1] - VECTOR(graph->is)[vid]);
    }
    if (! loops) {
        /* When loops should not be counted, we remove their contribution from the
         * previously computed degree. */
        if (mode & IGRAPH_OUT) {
            for (igraph_integer_t i = VECTOR(graph->os)[vid]; i < VECTOR(graph->os)[vid + 1]; i++) {
                if (VECTOR(graph->to)[ VECTOR(graph->oi)[i] ] == vid) {
                    *deg -= 1;
                }
            }
        }
        if (mode & IGRAPH_IN) {
            for (igraph_integer_t i = VECTOR(graph->is)[vid]; i < VECTOR(graph->is)[vid + 1]; i++) {
                if (VECTOR(graph->from)[ VECTOR(graph->ii)[i] ] == vid) {
                    *deg -= 1;
                }
            }
        }
    }

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup interface
 * \function igraph_degree
 * \brief The degree of some vertices in a graph.
 *
 * </para><para>
 * This function calculates the in-, out- or total degree of the
 * specified vertices.
 *
 * </para><para>
 * This function returns the result as a vector of \c igraph_integer_t
 * values. In applications where \c igraph_real_t is desired, use
 * \ref igraph_strength() with \c NULL weights.
 *
 * \param graph The graph.
 * \param res Integer vector, this will contain the result. It should be
 *        initialized and will be resized to be the appropriate size.
 * \param vids Vertex selector, giving the vertex IDs of which the degree will
 *        be calculated.
 * \param mode Defines the type of the degree for directed graphs. Valid modes are:
 *        \c IGRAPH_OUT, out-degree;
 *        \c IGRAPH_IN, in-degree;
 *        \c IGRAPH_ALL, total degree (sum of the
 *        in- and out-degree).
 *        This parameter is ignored for undirected graphs.
 * \param loops Boolean, gives whether the self-loops should be
 *        counted.
 * \return Error code:
 *         \c IGRAPH_EINVVID: invalid vertex ID.
 *         \c IGRAPH_EINVMODE: invalid mode argument.
 *
 * Time complexity: O(v) if \p loops is \c true, and
 * O(v*d) otherwise. v is the number of
 * vertices for which the degree will be calculated, and
 * d is their (average) degree.
 *
 * \sa \ref igraph_strength() for the version that takes into account
 * edge weights; \ref igraph_degree_1() to efficiently compute the
 * degree of a single vertex; \ref igraph_maxdegree() if you only need
 * the largest degree.
 *
 * \example examples/simple/igraph_degree.c
 */
igraph_error_t igraph_degree(const igraph_t *graph, igraph_vector_int_t *res,
                  const igraph_vs_t vids,
                  igraph_neimode_t mode, igraph_bool_t loops) {

    igraph_integer_t nodes_to_calc;
    igraph_integer_t i, j;
    igraph_vit_t vit;

    IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit);

    if (mode != IGRAPH_OUT && mode != IGRAPH_IN && mode != IGRAPH_ALL) {
        IGRAPH_ERROR("Invalid mode for degree calculation.", IGRAPH_EINVMODE);
    }

    if (! loops) {
        /* If the graph is known not to have loops, we can use the faster
         * loops == true code path, which has O(1) complexity instead of of O(d). */
        if (igraph_i_property_cache_has(graph, IGRAPH_PROP_HAS_LOOP) &&
            !igraph_i_property_cache_get_bool(graph, IGRAPH_PROP_HAS_LOOP)) {
            loops = true;
        }
    }

    nodes_to_calc = IGRAPH_VIT_SIZE(vit);
    if (!igraph_is_directed(graph)) {
        mode = IGRAPH_ALL;
    }

    IGRAPH_CHECK(igraph_vector_int_resize(res, nodes_to_calc));
    igraph_vector_int_null(res);

    if (loops) {
        if (mode & IGRAPH_OUT) {
            for (IGRAPH_VIT_RESET(vit), i = 0;
                 !IGRAPH_VIT_END(vit);
                 IGRAPH_VIT_NEXT(vit), i++) {
                igraph_integer_t vid = IGRAPH_VIT_GET(vit);
                VECTOR(*res)[i] += (VECTOR(graph->os)[vid + 1] - VECTOR(graph->os)[vid]);
            }
        }
        if (mode & IGRAPH_IN) {
            for (IGRAPH_VIT_RESET(vit), i = 0;
                 !IGRAPH_VIT_END(vit);
                 IGRAPH_VIT_NEXT(vit), i++) {
                igraph_integer_t vid = IGRAPH_VIT_GET(vit);
                VECTOR(*res)[i] += (VECTOR(graph->is)[vid + 1] - VECTOR(graph->is)[vid]);
            }
        }
    } else if (igraph_vs_is_all(&vids)) { /* no loops, calculating degree for all vertices */
        // When calculating degree for all vertices, iterating over edges is faster
        igraph_integer_t no_of_edges = igraph_ecount(graph);

        if (mode & IGRAPH_OUT) {
            for (igraph_integer_t edge = 0; edge < no_of_edges; ++edge) {
                igraph_integer_t from = IGRAPH_FROM(graph, edge);
                if (from != IGRAPH_TO(graph, edge)) {
                    VECTOR(*res)[from]++;
                }
            }
        }
        if (mode & IGRAPH_IN) {
            for (igraph_integer_t edge = 0; edge < no_of_edges; ++edge) {
                igraph_integer_t to = IGRAPH_TO(graph, edge);
                if (IGRAPH_FROM(graph, edge) != to) {
                    VECTOR(*res)[to]++;
                }
            }
        }
    } else { /* no loops */
        if (mode & IGRAPH_OUT) {
            for (IGRAPH_VIT_RESET(vit), i = 0;
                 !IGRAPH_VIT_END(vit);
                 IGRAPH_VIT_NEXT(vit), i++) {
                igraph_integer_t vid = IGRAPH_VIT_GET(vit);
                VECTOR(*res)[i] += (VECTOR(graph->os)[vid + 1] - VECTOR(graph->os)[vid]);
                for (j = VECTOR(graph->os)[vid];
                     j < VECTOR(graph->os)[vid + 1]; j++) {
                    if (VECTOR(graph->to)[ VECTOR(graph->oi)[j] ] == vid) {
                        VECTOR(*res)[i] -= 1;
                    }
                }
            }
        }
        if (mode & IGRAPH_IN) {
            for (IGRAPH_VIT_RESET(vit), i = 0;
                 !IGRAPH_VIT_END(vit);
                 IGRAPH_VIT_NEXT(vit), i++) {
                igraph_integer_t vid = IGRAPH_VIT_GET(vit);
                VECTOR(*res)[i] += (VECTOR(graph->is)[vid + 1] - VECTOR(graph->is)[vid]);
                for (j = VECTOR(graph->is)[vid];
                     j < VECTOR(graph->is)[vid + 1]; j++) {
                    if (VECTOR(graph->from)[ VECTOR(graph->ii)[j] ] == vid) {
                        VECTOR(*res)[i] -= 1;
                    }
                }
            }
        }
    }  /* loops */

    igraph_vit_destroy(&vit);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

/* These are unsafe macros. Only supply variable names, i.e. no
   expressions as parameters, otherwise nasty things can happen.

   BINSEARCH is an inline binary search in the 'edgelist' vector, which is
   assumed to be sorted in the order of indices stored in the 'iindex' vector.
   (So, [edgelist[iindex[x]] for x in 0..] is assumed to be sorted). 'N' must
   be the same as 'end' when invoking the macro but it must be a separate
   variable as we want to modify 'end' independently of 'N'. Upon exiting the
   macro, 'result' is the index of the _leftmost_ item in the sorted 'edgelist'
   (i.e. indexed by 'iindex') where the value was found, if it was found;
   otherwise 'pos' is left intact.

   FIND_DIRECTED_EDGE looks for an edge from 'xfrom' to 'xto' in the graph, and
   stores the ID of the edge in 'eid' if it is found; otherwise 'eid' is left
   intact.

   FIND_UNDIRECTED_EDGE looks for an edge between 'xfrom' and 'xto' in an
   undirected graph, swapping them if necessary. It stores the ID of the edge
   in 'eid' if it is found; otherwise 'eid' is left intact.
   */

#define BINSEARCH(start, end, value, iindex, edgelist, N, result, result_pos) \
    do { \
        while ((start) < (end)) { \
            igraph_integer_t mid =(start)+((end)-(start))/2; \
            igraph_integer_t e = VECTOR((iindex))[mid]; \
            if (VECTOR((edgelist))[e] < (value)) { \
                (start) = mid+1; \
            } else { \
                (end) = mid; \
            } \
        } \
        if ((start) < (N)) { \
            igraph_integer_t e = VECTOR((iindex))[(start)]; \
            if (VECTOR((edgelist))[e] == (value)) { \
                *(result) = e; \
                if (result_pos != 0) { *(result_pos) = start; } \
            } \
        } \
    } while (0)

#define FIND_DIRECTED_EDGE(graph,xfrom,xto,eid) \
    do { \
        igraph_integer_t start = VECTOR(graph->os)[xfrom]; \
        igraph_integer_t end = VECTOR(graph->os)[xfrom+1]; \
        igraph_integer_t N = end; \
        igraph_integer_t start2 = VECTOR(graph->is)[xto]; \
        igraph_integer_t end2 = VECTOR(graph->is)[xto+1]; \
        igraph_integer_t N2 = end2; \
        igraph_integer_t *nullpointer = NULL; \
        if (end-start < end2-start2) { \
            BINSEARCH(start, end, xto, graph->oi, graph->to, N, eid, nullpointer); \
        } else { \
            BINSEARCH(start2, end2, xfrom, graph->ii, graph->from, N2, eid, nullpointer); \
        } \
    } while (0)

#define FIND_UNDIRECTED_EDGE(graph, from, to, eid) \
    do { \
        igraph_integer_t xfrom1 = from > to ? from : to; \
        igraph_integer_t xto1 = from > to ? to : from; \
        FIND_DIRECTED_EDGE(graph, xfrom1, xto1, eid); \
    } while (0)

/**
 * \function igraph_get_eid
 * \brief Get the edge ID from the endpoints of an edge.
 *
 * For undirected graphs \c from and \c to are exchangeable.
 *
 * \param graph The graph object.
 * \param eid Pointer to an integer, the edge ID will be stored here.
 *        If \p error is false and no edge was found, <code>-1</code>
 *        will be returned.
 * \param from The starting point of the edge.
 * \param to The end point of the edge.
 * \param directed Logical constant, whether to search for directed
 *        edges in a directed graph. Ignored for undirected graphs.
 * \param error Logical scalar, whether to report an error if the edge
 *        was not found. If it is false, then <code>-1</code> will be
 *        assigned to \p eid. Note that invalid vertex IDs in input
 *        arguments (\p from or \p to) always trigger an error,
 *        regardless of this setting.
 * \return Error code.
 * \sa \ref igraph_edge() for the opposite operation, \ref igraph_get_all_eids_between()
 *     to retrieve all edge IDs between a pair of vertices.
 *
 * Time complexity: O(log (d)), where d is smaller of the out-degree
 * of \c from and in-degree of \c to if \p directed is true. If \p directed
 * is false, then it is O(log(d)+log(d2)), where d is the same as before and
 * d2 is the minimum of the out-degree of \c to and the in-degree of \c from.
 *
 * \example examples/simple/igraph_get_eid.c
 *
 * Added in version 0.2.</para><para>
 */

igraph_error_t igraph_get_eid(const igraph_t *graph, igraph_integer_t *eid,
                   igraph_integer_t from, igraph_integer_t to,
                   igraph_bool_t directed, igraph_bool_t error) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);

    if (from < 0 || to < 0 || from >= no_of_nodes || to >= no_of_nodes) {
        IGRAPH_ERROR("Cannot get edge ID.", IGRAPH_EINVVID);
    }

    *eid = -1;
    if (igraph_is_directed(graph)) {

        /* Directed graph */
        FIND_DIRECTED_EDGE(graph, from, to, eid);
        if (!directed && *eid < 0) {
            FIND_DIRECTED_EDGE(graph, to, from, eid);
        }

    } else {

        /* Undirected graph, they only have one mode */
        FIND_UNDIRECTED_EDGE(graph, from, to, eid);

    }

    if (*eid < 0) {
        if (error) {
            IGRAPH_ERROR("Cannot get edge ID, no such edge", IGRAPH_EINVAL);
        }
    }

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_get_eids
 * Return edge IDs based on the adjacent vertices.
 *
 * The pairs of vertex IDs for which the edges are looked up are taken
 * consecutively from the \c pairs vector, i.e. <code>VECTOR(pairs)[0]</code>
 * and <code>VECTOR(pairs)[1]</code> specify the first pair,
 * <code>VECTOR(pairs)[2]</code> and <code>VECTOR(pairs)[3]</code> the second
 * pair, etc.
 *
 * </para><para>
 * If you have a sequence of vertex IDs that describe a \em path on the graph,
 * use \ref igraph_expand_path_to_pairs() to convert them to a list of vertex
 * pairs along the path.
 *
 * </para><para>
 * If the \c error argument is true, then it is an error to specify pairs
 * of vertices that are not connected. Otherwise -1 is reported for vertex pairs
 * without at least one edge between them.
 *
 * </para><para>
 * If there are multiple edges in the graph, then these are ignored;
 * i.e. for a given pair of vertex IDs, igraph always returns the same edge ID,
 * even if the pair appears multiple times in \c pairs.
 *
 * \param graph The input graph.
 * \param eids Pointer to an initialized vector, the result is stored
 *        here. It will be resized as needed.
 * \param pairs Vector giving pairs of vertices to fetch the edges for.
 * \param directed Logical scalar, whether to consider edge directions
 *        in directed graphs. This is ignored for undirected graphs.
 * \param error Logical scalar, whether it is an error to supply
 *        non-connected vertices. If false, then -1 is
 *        returned for non-connected pairs.
 * \return Error code.
 *
 * Time complexity: O(n log(d)), where n is the number of queried
 * edges and d is the average degree of the vertices.
 *
 * \sa \ref igraph_get_eid() for a single edge.
 *
 * \example examples/simple/igraph_get_eids.c
 */
igraph_error_t igraph_get_eids(const igraph_t *graph, igraph_vector_int_t *eids,
                    const igraph_vector_int_t *pairs,
                    igraph_bool_t directed, igraph_bool_t error) {

    igraph_integer_t n = pairs ? igraph_vector_int_size(pairs) : 0;
    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t i;
    igraph_integer_t eid = -1;

    if (n == 0) {
        igraph_vector_int_clear(eids);
        return IGRAPH_SUCCESS;
    }

    if (n % 2 != 0) {
        IGRAPH_ERROR("Cannot get edge IDs, invalid length of edge IDs",
                     IGRAPH_EINVAL);
    }

    if (!igraph_vector_int_isininterval(pairs, 0, no_of_nodes - 1)) {
        IGRAPH_ERROR("Cannot get edge IDs, invalid vertex ID", IGRAPH_EINVVID);
    }

    IGRAPH_CHECK(igraph_vector_int_resize(eids, n / 2));

    if (igraph_is_directed(graph)) {
        for (i = 0; i < n / 2; i++) {
            igraph_integer_t from = VECTOR(*pairs)[2 * i];
            igraph_integer_t to = VECTOR(*pairs)[2 * i + 1];

            eid = -1;
            FIND_DIRECTED_EDGE(graph, from, to, &eid);
            if (!directed && eid < 0) {
                FIND_DIRECTED_EDGE(graph, to, from, &eid);
            }

            VECTOR(*eids)[i] = eid;
            if (eid < 0 && error) {
                IGRAPH_ERROR("Cannot get edge ID, no such edge", IGRAPH_EINVAL);
            }
        }
    } else {
        for (i = 0; i < n / 2; i++) {
            igraph_integer_t from = VECTOR(*pairs)[2 * i];
            igraph_integer_t to = VECTOR(*pairs)[2 * i + 1];

            eid = -1;
            FIND_UNDIRECTED_EDGE(graph, from, to, &eid);
            VECTOR(*eids)[i] = eid;
            if (eid < 0 && error) {
                IGRAPH_ERROR("Cannot get edge ID, no such edge", IGRAPH_EINVAL);
            }
        }
    }

    return IGRAPH_SUCCESS;
}

#undef FIND_DIRECTED_EDGE
#undef FIND_UNDIRECTED_EDGE

#define FIND_ALL_DIRECTED_EDGES(graph, xfrom, xto, eidvec) \
    do { \
        igraph_integer_t start = VECTOR(graph->os)[xfrom]; \
        igraph_integer_t end = VECTOR(graph->os)[xfrom+1]; \
        igraph_integer_t N = end; \
        igraph_integer_t start2 = VECTOR(graph->is)[xto]; \
        igraph_integer_t end2 = VECTOR(graph->is)[xto+1]; \
        igraph_integer_t N2 = end2; \
        igraph_integer_t eid = -1; \
        igraph_integer_t pos = -1; \
        if (end-start < end2-start2) { \
            BINSEARCH(start, end, xto, graph->oi, graph->to, N, &eid, &pos); \
            while (pos >= 0 && pos < N) { \
                eid = VECTOR(graph->oi)[pos++]; \
                if (VECTOR(graph->to)[eid] != xto) { break; } \
                IGRAPH_CHECK(igraph_vector_int_push_back(eidvec, eid)); \
            } \
        } else { \
            BINSEARCH(start2, end2, xfrom, graph->ii, graph->from, N2, &eid, &pos); \
            while (pos >= 0 && pos < N2) { \
                eid = VECTOR(graph->ii)[pos++]; \
                if (VECTOR(graph->from)[eid] != xfrom) { break; } \
                IGRAPH_CHECK(igraph_vector_int_push_back(eidvec, eid)); \
            } \
        } \
    } while (0)

#define FIND_ALL_UNDIRECTED_EDGES(graph, from, to, eidvec) \
    do { \
        igraph_integer_t xfrom1 = from > to ? from : to; \
        igraph_integer_t xto1 = from > to ? to : from; \
        FIND_ALL_DIRECTED_EDGES(graph, xfrom1, xto1, eidvec); \
    } while (0)

/**
 * \function igraph_get_all_eids_between
 * \brief Returns all edge IDs between a pair of vertices.
 *
 * </para><para>
 * For undirected graphs \c source and \c target are exchangeable.
 *
 * \param graph The input graph.
 * \param eids Pointer to an initialized vector, the result is stored
 *        here. It will be resized as needed.
 * \param source The ID of the source vertex
 * \param target The ID of the target vertex
 * \param directed Logical scalar, whether to consider edge directions
 *        in directed graphs. This is ignored for undirected graphs.
 * \return Error code.
 *
 * Time complexity: TODO
 *
 * \sa \ref igraph_get_eid() for a single edge.
 */
igraph_error_t igraph_get_all_eids_between(
    const igraph_t *graph, igraph_vector_int_t *eids,
    igraph_integer_t source, igraph_integer_t target, igraph_bool_t directed
) {
    igraph_integer_t no_of_nodes = igraph_vcount(graph);

    if (source < 0 || source >= no_of_nodes) {
        IGRAPH_ERROR("Cannot get edge IDs, invalid source vertex ID", IGRAPH_EINVVID);
    }

    if (target < 0 || target >= no_of_nodes) {
        IGRAPH_ERROR("Cannot get edge IDs, invalid target vertex ID", IGRAPH_EINVVID);
    }

    igraph_vector_int_clear(eids);

    if (igraph_is_directed(graph)) {
        /* look in the specified direction first */
        FIND_ALL_DIRECTED_EDGES(graph, source, target, eids);
        if (!directed) {
            /* look in the reverse direction as well */
            FIND_ALL_DIRECTED_EDGES(graph, target, source, eids);
        }
    } else {
        FIND_ALL_UNDIRECTED_EDGES(graph, source, target, eids);
    }

    return IGRAPH_SUCCESS;
}

#undef FIND_DIRECTED_EDGE
#undef FIND_UNDIRECTED_EDGE
#undef BINSEARCH

/**
 * \function igraph_incident
 * \brief Gives the incident edges of a vertex.
 *
 * \param graph The graph object.
 * \param eids An initialized vector. It will be resized
 * to hold the result.
 * \param pnode A vertex ID.
 * \param mode Specifies what kind of edges to include for directed
 * graphs. \c IGRAPH_OUT means only outgoing edges, \c IGRAPH_IN only
 * incoming edges, \c IGRAPH_ALL both. This parameter is ignored for
 * undirected graphs.
 * \return Error code. \c IGRAPH_EINVVID: invalid \p pnode argument,
 *   \c IGRAPH_EINVMODE: invalid \p mode argument.
 *
 * Added in version 0.2.</para><para>
 *
 * Time complexity: O(d), the number of incident edges to \p pnode.
 */

igraph_error_t igraph_incident(const igraph_t *graph, igraph_vector_int_t *eids, igraph_integer_t pnode,
        igraph_neimode_t mode) {
    if (!igraph_is_directed(graph) || mode == IGRAPH_ALL) {
        return igraph_i_incident(graph, eids, pnode, mode, IGRAPH_LOOPS_TWICE);
    } else {
        return igraph_i_incident(graph, eids, pnode, mode, IGRAPH_LOOPS_ONCE);
    }
}

igraph_error_t igraph_i_incident(const igraph_t *graph, igraph_vector_int_t *eids, igraph_integer_t pnode,
        igraph_neimode_t mode, igraph_loops_t loops) {
    igraph_integer_t length = 0, idx = 0;
    igraph_integer_t i, j;
    igraph_integer_t node = pnode;
    igraph_bool_t directed = igraph_is_directed(graph);

    if (node < 0 || node > igraph_vcount(graph) - 1) {
        IGRAPH_ERROR("Given vertex is not in the graph.", IGRAPH_EINVVID);
    }
    if (mode != IGRAPH_OUT && mode != IGRAPH_IN &&
        mode != IGRAPH_ALL) {
        IGRAPH_ERROR("Mode should be either IGRAPH_OUT, IGRAPH_IN or IGRAPH_ALL.", IGRAPH_EINVMODE);
    }

    if (!directed) {
        mode = IGRAPH_ALL;
    }

    if (mode != IGRAPH_ALL && loops == IGRAPH_LOOPS_TWICE) {
        IGRAPH_ERROR("For a directed graph (with directions not ignored), "
                     "IGRAPH_LOOPS_TWICE does not make sense.\n", IGRAPH_EINVAL);
    }

    /* Calculate needed space first & allocate it */
    /* Note that 'mode' is treated as a bit field here; it's okay because
     * IGRAPH_ALL = IGRAPH_IN | IGRAPH_OUT, bit-wise */
    if (mode & IGRAPH_OUT) {
        length += (VECTOR(graph->os)[node + 1] - VECTOR(graph->os)[node]);
    }
    if (mode & IGRAPH_IN) {
        length += (VECTOR(graph->is)[node + 1] - VECTOR(graph->is)[node]);
    }

    IGRAPH_CHECK(igraph_vector_int_resize(eids, length));

    /* The loops below produce an ordering what is consistent with the
     * ordering returned by igraph_neighbors(), and this should be preserved.
     * We are dealing with two sorted lists; one for the successors and one
     * for the predecessors. If we have requested only one of them, we have
     * an easy job. If we have requested both, we need to merge the two lists
     * to ensure that the output is sorted by the vertex IDs of the "other"
     * endpoint of the affected edges */
    if (!directed || mode != IGRAPH_ALL) {
        /* We did not ask for both directions; this is the easy case */

        if (mode & IGRAPH_OUT) {
            j = VECTOR(graph->os)[node + 1];
            for (i = VECTOR(graph->os)[node]; i < j; i++) {
                igraph_integer_t edge = VECTOR(graph->oi)[i];
                igraph_integer_t other = VECTOR(graph->to)[edge];
                if (loops == IGRAPH_NO_LOOPS && other == pnode) {
                    length--;
                } else {
                    VECTOR(*eids)[idx++] = edge;
                }
            }
        }

        if (mode & IGRAPH_IN) {
            j = VECTOR(graph->is)[node + 1];
            for (i = VECTOR(graph->is)[node]; i < j; i++) {
                igraph_integer_t edge = VECTOR(graph->ii)[i];
                igraph_integer_t other = VECTOR(graph->from)[edge];
                if ((loops == IGRAPH_NO_LOOPS || (loops == IGRAPH_LOOPS_ONCE && !directed)) && other == pnode) {
                    length--;
                } else {
                    VECTOR(*eids)[idx++] = edge;
                }
            }
        }
    } else {
        /* both in- and out- neighbors in a directed graph,
           we need to merge the two 'vectors' */
        igraph_integer_t j1 = VECTOR(graph->os)[node + 1];
        igraph_integer_t j2 = VECTOR(graph->is)[node + 1];
        igraph_integer_t i1 = VECTOR(graph->os)[node];
        igraph_integer_t i2 = VECTOR(graph->is)[node];
        igraph_integer_t eid1, eid2;
        igraph_integer_t n1, n2;
        igraph_bool_t seen_loop_edge = false;

        while (i1 < j1 && i2 < j2) {
            eid1 = VECTOR(graph->oi)[i1];
            eid2 = VECTOR(graph->ii)[i2];
            n1 = VECTOR(graph->to)[eid1];
            n2 = VECTOR(graph->from)[eid2];
            if (n1 < n2) {
                i1++;
                VECTOR(*eids)[idx++] = eid1;
            } else if (n1 > n2) {
                i2++;
                VECTOR(*eids)[idx++] = eid2;
            } else if (n1 != pnode) {
                /* multiple edge */
                i1++;
                i2++;
                VECTOR(*eids)[idx++] = eid1;
                VECTOR(*eids)[idx++] = eid2;
            } else {
                /* loop edge */
                i1++;
                i2++;
                if (loops == IGRAPH_NO_LOOPS) {
                    length -= 2;
                } else if (loops == IGRAPH_LOOPS_ONCE) {
                    length--;
                    if (!seen_loop_edge) {
                        VECTOR(*eids)[idx++] = eid1;
                    } else {
                        VECTOR(*eids)[idx++] = eid2;
                    }
                    seen_loop_edge = !seen_loop_edge;
                } else {
                    VECTOR(*eids)[idx++] = eid1;
                    VECTOR(*eids)[idx++] = eid2;
                }
            }
        }

        while (i1 < j1) {
            eid1 = VECTOR(graph->oi)[i1++];
            VECTOR(*eids)[idx++] = eid1;
        }

        while (i2 < j2) {
            eid2 = VECTOR(graph->ii)[i2++];
            VECTOR(*eids)[idx++] = eid2;
        }
    }
    IGRAPH_CHECK(igraph_vector_int_resize(eids, length));
    return IGRAPH_SUCCESS;
#undef DEDUPLICATE_IF_NEEDED
}


/**
 * \function igraph_is_same_graph
 * \brief Are two graphs identical as labelled graphs?
 *
 * Two graphs are considered to be the same if they have the same vertex and edge sets.
 * Graphs which are the same may have multiple different representations in igraph,
 * hence the need for this function.
 *
 * </para><para>
 * This function verifies that the two graphs have the same directedness, the same
 * number of vertices, and that they contain precisely the same edges (regardless of their ordering)
 * when written in terms of vertex indices. Graph attributes are not taken into account.
 *
 * </para><para>
 * This concept is different from isomorphism. For example, the graphs
 * <code>0-1, 2-1</code> and <code>1-2, 0-1</code> are considered the same
 * because they only differ in the ordering of their edge lists and the ordering
 * of vertices in an undirected edge. However, they are not the same as
 * <code>0-2, 1-2</code>, even though they are isomorphic to it.
 * Note that this latter graph contains the edge <code>0-2</code>
 * while the former two do not — thus their edge sets differ.
 *
 * \param graph1 The first graph object.
 * \param graph2 The second graph object.
 * \param res The result will be stored here.
 * \return Error code.
 *
 * Time complexity: O(E), the number of edges in the graphs.
 *
 * \sa \ref igraph_isomorphic() to test if two graphs are isomorphic.
 */

igraph_error_t igraph_is_same_graph(const igraph_t *graph1, const igraph_t *graph2, igraph_bool_t *res) {
    igraph_integer_t nv1 = igraph_vcount(graph1);
    igraph_integer_t nv2 = igraph_vcount(graph2);
    igraph_integer_t ne1 = igraph_ecount(graph1);
    igraph_integer_t ne2 = igraph_ecount(graph2);
    igraph_integer_t i, eid1, eid2;

    *res = false; /* Assume that the graphs differ */

    /* Check for same number of vertices/edges */
    if ((nv1 != nv2) || (ne1 != ne2)) {
        return IGRAPH_SUCCESS;
    }

    /* Check for same directedness */
    if (igraph_is_directed(graph1) != igraph_is_directed(graph2)) {
        return IGRAPH_SUCCESS;
    }

    /* Vertices have no names, so they must be 0 to nv - 1 */

    /* Edges are double sorted in the current representations ii/oi of
     * igraph_t (ii: by incoming, then outgoing, oi: vice versa), so
     * we just need to check them one by one. If that representation
     * changes, this part will need to change too.
     *
     * Furthermore, in the current representation the "source" of undirected
     * edges always has a vertex index that is no larger than that of the
     * "target".
     */
    for (i = 0; i < ne1; i++) {
        eid1 = VECTOR(graph1->ii)[i];
        eid2 = VECTOR(graph2->ii)[i];

        /* Check they have the same source */
        if (IGRAPH_FROM(graph1, eid1) != IGRAPH_FROM(graph2, eid2)) {
            return IGRAPH_SUCCESS;
        }

        /* Check they have the same target */
        if (IGRAPH_TO(graph1, eid1) != IGRAPH_TO(graph2, eid2)) {
            return IGRAPH_SUCCESS;
        }
    }

    *res = true; /* No difference was found, graphs are the same */
    return IGRAPH_SUCCESS;
}


/* Reverses the direction of all edges in a directed graph.
 * The graph is modified in-place.
 * Attributes are preserved.
 */
igraph_error_t igraph_i_reverse(igraph_t *graph) {

    /* Nothing to do for undirected graphs. */
    if (! igraph_is_directed(graph)) {
        return IGRAPH_SUCCESS;
    }

    igraph_vector_int_swap(&graph->to, &graph->from);
    igraph_vector_int_swap(&graph->oi, &graph->ii);
    igraph_vector_int_swap(&graph->os, &graph->is);

    return IGRAPH_SUCCESS;
}
