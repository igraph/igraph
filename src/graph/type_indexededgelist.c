/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2005-2012  Gabor Csardi <csardi.gabor@gmail.com>
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
#include "graph/neighbors.h"

/* Internal functions */

static igraph_error_t igraph_i_create_start(
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
 */

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
 * \example examples/simple/igraph_empty.c
 */
igraph_error_t igraph_empty(igraph_t *graph, igraph_integer_t n, igraph_bool_t directed) {
    return igraph_empty_attrs(graph, n, directed, 0);
}


/**
 * \ingroup interface
 * \function igraph_empty_attrs
 * \brief Creates an empty graph with some vertices, no edges and some graph attributes.
 *
 * </para><para>
 * Use this instead of \ref igraph_empty() if you wish to add some graph
 * attributes right after initialization. This function is currently
 * not very interesting for the ordinary user. Just supply 0 here or
 * use \ref igraph_empty().
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
 * \param attr The attributes.
 * \return Error code:
 *         \c IGRAPH_EINVAL: invalid number of vertices.
 *
 * Time complexity: O(|V|) for a graph with
 * |V| vertices (and no edges).
 */
igraph_error_t igraph_empty_attrs(igraph_t *graph, igraph_integer_t n, igraph_bool_t directed, void* attr) {

    if (n < 0) {
        IGRAPH_ERROR("cannot create empty graph with negative number of vertices",
                     IGRAPH_EINVAL);
    }

    graph->n = 0;
    graph->directed = directed;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&graph->from, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&graph->to, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&graph->oi, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&graph->ii, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&graph->os, 1);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&graph->is, 1);

    VECTOR(graph->os)[0] = 0;
    VECTOR(graph->is)[0] = 0;

    /* init attributes */
    graph->attr = 0;
    IGRAPH_CHECK(igraph_i_attribute_init(graph, attr));

    /* add the vertices */
    IGRAPH_CHECK(igraph_add_vertices(graph, n, 0));

    IGRAPH_FINALLY_CLEAN(6);
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
    IGRAPH_CHECK(igraph_vector_int_copy(&to->from, &from->from));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &to->from);
    IGRAPH_CHECK(igraph_vector_int_copy(&to->to, &from->to));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &to->to);
    IGRAPH_CHECK(igraph_vector_int_copy(&to->oi, &from->oi));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &to->oi);
    IGRAPH_CHECK(igraph_vector_int_copy(&to->ii, &from->ii));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &to->ii);
    IGRAPH_CHECK(igraph_vector_int_copy(&to->os, &from->os));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &to->os);
    IGRAPH_CHECK(igraph_vector_int_copy(&to->is, &from->is));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &to->is);

    IGRAPH_I_ATTRIBUTE_COPY(to, from, 1, 1, 1); /* does IGRAPH_CHECK */

    IGRAPH_FINALLY_CLEAN(6);
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
 * want to add new vertices, call igraph_add_vertices() first.
 * \param graph The graph to which the edges will be added.
 * \param edges The edges themselves.
 * \param attr The attributes of the new edges, only used by high level
 *        interfaces currently, you can supply 0 here.
 * \return Error code:
 *    \c IGRAPH_EINVEVECTOR: invalid (odd)
 *    edges vector length, \c IGRAPH_EINVVID:
 *    invalid vertex ID in edges vector.
 *
 * This function invalidates all iterators.
 *
 * </para><para>
 * Time complexity: O(|V|+|E|) where
 * |V| is the number of vertices and
 * |E| is the number of
 * edges in the \em new, extended graph.
 *
 * \example examples/simple/igraph_add_edges.c
 */
igraph_error_t igraph_add_edges(igraph_t *graph, const igraph_vector_int_t *edges,
                     void *attr) {
    igraph_integer_t no_of_edges = igraph_vector_int_size(&graph->from);
    igraph_integer_t edges_to_add = igraph_vector_int_size(edges) / 2;
    igraph_integer_t i = 0;
    igraph_error_handler_t *oldhandler;
    igraph_error_t ret1, ret2;
    igraph_vector_int_t newoi, newii;
    igraph_bool_t directed = igraph_is_directed(graph);

    if (igraph_vector_int_size(edges) % 2 != 0) {
        IGRAPH_ERROR("invalid (odd) length of edges vector", IGRAPH_EINVEVECTOR);
    }
    if (!igraph_vector_int_isininterval(edges, 0, igraph_vcount(graph) - 1)) {
        IGRAPH_ERROR("cannot add edges", IGRAPH_EINVVID);
    }

    /* from & to */
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

    /* disable the error handler temporarily */
    oldhandler = igraph_set_error_handler(igraph_error_handler_ignore);

    /* oi & ii */
    ret1 = igraph_vector_int_init(&newoi, no_of_edges);
    ret2 = igraph_vector_int_init(&newii, no_of_edges);
    if (ret1 != 0 || ret2 != 0) {
        igraph_vector_int_resize(&graph->from, no_of_edges); /* gets smaller */
        igraph_vector_int_resize(&graph->to, no_of_edges);   /* gets smaller */
        igraph_set_error_handler(oldhandler);
        IGRAPH_ERROR("cannot add edges", IGRAPH_ERROR_SELECT_2(ret1, ret2));
    }
    ret1 = igraph_vector_int_pair_order(&graph->from, &graph->to, &newoi, graph->n);
    ret2 = igraph_vector_int_pair_order(&graph->to, &graph->from, &newii, graph->n);
    if (ret1 != 0 || ret2 != 0) {
        igraph_vector_int_resize(&graph->from, no_of_edges);
        igraph_vector_int_resize(&graph->to, no_of_edges);
        igraph_vector_int_destroy(&newoi);
        igraph_vector_int_destroy(&newii);
        igraph_set_error_handler(oldhandler);
        IGRAPH_ERROR("cannot add edges", IGRAPH_ERROR_SELECT_2(ret1, ret2));
    }

    /* Attributes */
    if (graph->attr) {
        igraph_set_error_handler(oldhandler);
        ret1 = igraph_i_attribute_add_edges(graph, edges, attr);
        igraph_set_error_handler(igraph_error_handler_ignore);
        if (ret1 != 0) {
            igraph_vector_int_resize(&graph->from, no_of_edges);
            igraph_vector_int_resize(&graph->to, no_of_edges);
            igraph_vector_int_destroy(&newoi);
            igraph_vector_int_destroy(&newii);
            igraph_set_error_handler(oldhandler);
            IGRAPH_ERROR("cannot add edges", ret1);
        }
    }

    /* os & is, its length does not change, error safe */
    igraph_i_create_start(&graph->os, &graph->from, &newoi, graph->n);
    igraph_i_create_start(&graph->is, &graph->to, &newii, graph->n);

    /* everything went fine  */
    igraph_vector_int_destroy(&graph->oi);
    igraph_vector_int_destroy(&graph->ii);
    graph->oi = newoi;
    graph->ii = newii;
    igraph_set_error_handler(oldhandler);

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
 * \param nv Non-negative integer giving the number of
 *           vertices to add.
 * \param attr The attributes of the new vertices, only used by
 *           high level interfaces, you can supply 0 here.
 * \return Error code:
 *         \c IGRAPH_EINVAL: invalid number of new
 *         vertices.
 *
 * Time complexity: O(|V|) where
 * |V| is
 * the number of vertices in the \em new, extended graph.
 *
 * \example examples/simple/igraph_add_vertices.c
 */
igraph_error_t igraph_add_vertices(igraph_t *graph, igraph_integer_t nv, void *attr) {
    igraph_integer_t ec = igraph_ecount(graph);
    igraph_integer_t i;

    if (nv < 0) {
        IGRAPH_ERROR("cannot add negative number of vertices", IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_vector_int_reserve(&graph->os, graph->n + nv + 1));
    IGRAPH_CHECK(igraph_vector_int_reserve(&graph->is, graph->n + nv + 1));

    igraph_vector_int_resize(&graph->os, graph->n + nv + 1); /* reserved */
    igraph_vector_int_resize(&graph->is, graph->n + nv + 1); /* reserved */
    for (i = graph->n + 1; i < graph->n + nv + 1; i++) {
        VECTOR(graph->os)[i] = ec;
        VECTOR(graph->is)[i] = ec;
    }

    graph->n += nv;

    if (graph->attr) {
        IGRAPH_CHECK(igraph_i_attribute_add_vertices(graph, nv, attr));
    }

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup interface
 * \function igraph_delete_edges
 * \brief Removes edges from a graph.
 *
 * </para><para>
 * The edges to remove are given as an edge selector.
 *
 * </para><para>
 * This function cannot remove vertices, they will be kept, even if
 * they lose all their edges.
 *
 * </para><para>
 * This function invalidates all iterators.
 * \param graph The graph to work on.
 * \param edges The edges to remove.
 * \return Error code.
 *
 * Time complexity: O(|V|+|E|) where
 * |V|
 * and |E| are the number of vertices
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
    igraph_vector_int_t newoi;

    igraph_bool_t *mark;
    igraph_integer_t i, j;

    mark = IGRAPH_CALLOC(no_of_edges, int);
    if (mark == 0) {
        IGRAPH_ERROR("Cannot delete edges", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, mark);

    IGRAPH_CHECK(igraph_eit_create(graph, edges, &eit));
    IGRAPH_FINALLY(igraph_eit_destroy, &eit);

    for (IGRAPH_EIT_RESET(eit); !IGRAPH_EIT_END(eit); IGRAPH_EIT_NEXT(eit)) {
        igraph_integer_t e = IGRAPH_EIT_GET(eit);
        if (mark[e] == 0) {
            edges_to_remove++;
            mark[e]++;
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
        if (mark[i] == 0) {
            VECTOR(newfrom)[j] = VECTOR(graph->from)[i];
            VECTOR(newto)[j] = VECTOR(graph->to)[i];
            j++;
        }
    }

    /* Create index, this might require additional memory */
    IGRAPH_VECTOR_INT_INIT_FINALLY(&newoi, remaining_edges);
    IGRAPH_CHECK(igraph_vector_int_pair_order(&newfrom, &newto, &newoi, no_of_nodes));
    IGRAPH_CHECK(igraph_vector_int_pair_order(&newto, &newfrom, &graph->ii, no_of_nodes));

    /* Edge attributes, we need an index that gives the IDs of the
       original edges for every new edge.
    */
    if (graph->attr) {
        igraph_vector_int_t idx;
        IGRAPH_VECTOR_INT_INIT_FINALLY(&idx, remaining_edges);
        for (i = 0, j = 0; i < no_of_edges; i++) {
            if (mark[i] == 0) {
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
    graph->from = newfrom;
    graph->to = newto;
    graph->oi = newoi;
    IGRAPH_FINALLY_CLEAN(3);

    IGRAPH_FREE(mark);
    IGRAPH_FINALLY_CLEAN(1);

    /* Create start vectors, no memory is needed for this */
    igraph_i_create_start(&graph->os, &graph->from, &graph->oi, no_of_nodes);
    igraph_i_create_start(&graph->is, &graph->to,   &graph->ii, no_of_nodes);

    /* Nothing to deallocate... */
    return IGRAPH_SUCCESS;
}

/**
 * \ingroup interface
 * \function igraph_delete_vertices
 * \brief Removes vertices (with all their edges) from the graph.
 *
 * </para><para>
 * This function changes the IDs of the vertices (except in some very
 * special cases, but these should not be relied on anyway).
 *
 * </para><para>
 * This function invalidates all iterators.
 *
 * \param graph The graph to work on.
 * \param vertices The IDs of the vertices to remove in a
 *                 vector. The vector may contain the same id more
 *                 than once.
 * \return Error code:
 *         \c IGRAPH_EINVVID: invalid vertex ID.
 *
 * Time complexity: O(|V|+|E|),
 * |V| and
 * |E| are the number of vertices and
 * edges in the original graph.
 *
 * \example examples/simple/igraph_delete_vertices.c
 */
igraph_error_t igraph_delete_vertices(igraph_t *graph, const igraph_vs_t vertices) {
    return igraph_delete_vertices_idx(graph, vertices, /* idx= */ 0,
                                      /* invidx= */ 0);
}

igraph_error_t igraph_delete_vertices_idx(igraph_t *graph, const igraph_vs_t vertices,
                               igraph_vector_int_t *idx,
                               igraph_vector_int_t *invidx) {

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

    IGRAPH_CHECK(igraph_i_create_start(&newgraph.os, &newgraph.from,
                                       &newgraph.oi, remaining_vertices));
    IGRAPH_CHECK(igraph_i_create_start(&newgraph.is, &newgraph.to,
                                       &newgraph.ii, remaining_vertices));

    /* attributes */
    IGRAPH_I_ATTRIBUTE_COPY(&newgraph, graph,
                            /*graph=*/ 1, /*vertex=*/0, /*edge=*/0);
    IGRAPH_FINALLY_CLEAN(6);
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
        if ((loops == IGRAPH_NO_LOOPS && vertex == pnode) ||                               \
                (loops == IGRAPH_LOOPS_ONCE && vertex == pnode && last_added == pnode) ||  \
                (multiple == IGRAPH_NO_MULTIPLE && vertex == last_added)) {                \
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
 *
 */

static igraph_error_t igraph_i_create_start(
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
 * \return Logical value, <code>TRUE</code> if the graph is directed,
 * <code>FALSE</code> otherwise.
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
 * \function igraph_degree
 * \brief The degree of some vertices in a graph.
 *
 * </para><para>
 * This function calculates the in-, out- or total degree of the
 * specified vertices.
 * \param graph The graph.
 * \param res Vector, this will contain the result. It should be
 *        initialized and will be resized to be the appropriate size.
 * \param vids Vector, giving the vertex IDs of which the degree will
 *        be calculated.
 * \param mode Defines the type of the degree. Valid modes are:
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
 * Time complexity: O(v) if
 * loops is
 * TRUE, and
 * O(v*d)
 * otherwise. v is the number of
 * vertices for which the degree will be calculated, and
 * d is their (average) degree.
 *
 * \sa \ref igraph_strength() for the version that takes into account
 * edge weights.
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
        IGRAPH_ERROR("degree calculation failed", IGRAPH_EINVMODE);
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

/**
 * \function igraph_edge
 * \brief Gives the head and tail vertices of an edge.
 *
 * \param graph The graph object.
 * \param eid The edge ID.
 * \param from Pointer to an \type igraph_integer_t. The tail (head) of
 * the edge will be placed here for undirected (directed) graphs.
 * \param to Pointer to an \type igraph_integer_t. The head (tail) of the
 * edge will be placed here for undirected (directed) graphs.
 * \return Error code. The current implementation always returns with
 * success.
 * \sa \ref igraph_get_eid() for the opposite operation;
 *     \ref igraph_edges() to get the endpoints of several edges;
 *     \ref IGRAPH_TO(), \ref IGRAPH_FROM() and \ref IGRAPH_OTHER() for
 *     a faster but non-error-checked version.
 *
 * Added in version 0.2.</para><para>
 *
 * Time complexity: O(1).
 */

igraph_error_t igraph_edge(const igraph_t *graph, igraph_integer_t eid,
                igraph_integer_t *from, igraph_integer_t *to) {

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
 *     \ref igraph_get_eids() and \ref igraph_get_eids_multi()
 *     for the opposite operation;
 *     \ref igraph_edge() for getting the endpoints of a single edge;
 *     \ref IGRAPH_TO(), \ref IGRAPH_FROM() and \ref IGRAPH_OTHER() for
 *     a faster but non-error-checked method.
 *
 * Time complexity: O(k) where k is the number of edges in the selector.
 */

igraph_error_t igraph_edges(const igraph_t *graph, igraph_es_t eids,
                 igraph_vector_int_t *edges) {

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

/* This is an unsafe macro. Only supply variable names, i.e. no
   expressions as parameters, otherwise nasty things can happen */

#define BINSEARCH(start,end,value,iindex,edgelist,N,pos)     \
    do {                                                      \
        while ((start) < (end)) {                                 \
            igraph_integer_t mid=(start)+((end)-(start))/2;                 \
            igraph_integer_t e= VECTOR((iindex))[mid];            \
            if (VECTOR((edgelist))[e] < (value)) {                  \
                (start)=mid+1;                                        \
            } else {                                                \
                (end)=mid;                                            \
            }                                                       \
        }                                                         \
        if ((start)<(N)) {                                        \
            igraph_integer_t e= VECTOR((iindex))[(start)];        \
            if (VECTOR((edgelist))[e] == (value)) {                 \
                *(pos)= e;                                         \
            }                                                       \
        } } while(0)

#define FIND_DIRECTED_EDGE(graph,xfrom,xto,eid)                     \
    do {                                                              \
        igraph_integer_t start= VECTOR(graph->os)[xfrom];         \
        igraph_integer_t end= VECTOR(graph->os)[xfrom+1];         \
        igraph_integer_t N=end;                                                 \
        igraph_integer_t start2= VECTOR(graph->is)[xto];          \
        igraph_integer_t end2= VECTOR(graph->is)[xto+1];          \
        igraph_integer_t N2=end2;                                               \
        if (end-start<end2-start2) {                                    \
            BINSEARCH(start,end,xto,graph->oi,graph->to,N,eid);           \
        } else {                                                        \
            BINSEARCH(start2,end2,xfrom,graph->ii,graph->from,N2,eid);    \
        }                                                               \
    } while (0)

#define FIND_UNDIRECTED_EDGE(graph,from,to,eid)                     \
    do {                                                              \
        igraph_integer_t xfrom1= from > to ? from : to;                         \
        igraph_integer_t xto1= from > to ? to : from;                           \
        FIND_DIRECTED_EDGE(graph,xfrom1,xto1,eid);                      \
    } while (0)

/**
 * \function igraph_get_eid
 * \brief Get the edge ID from the end points of an edge.
 *
 * For undirected graphs \c pfrom and \c pto are exchangeable.
 *
 * \param graph The graph object.
 * \param eid Pointer to an integer, the edge ID will be stored here.
 * \param pfrom The starting point of the edge.
 * \param pto The end point of the edge.
 * \param directed Logical constant, whether to search for directed
 *        edges in a directed graph. Ignored for undirected graphs.
 * \param error Logical scalar, whether to report an error if the edge
 *        was not found. If it is false, then -1 will be assigned to \p eid.
 * \return Error code.
 * \sa \ref igraph_edge() for the opposite operation.
 *
 * Time complexity: O(log (d)), where d is smaller of the out-degree
 * of \c pfrom and in-degree of \c pto if \p directed is true. If \p directed
 * is false, then it is O(log(d)+log(d2)), where d is the same as before and
 * d2 is the minimum of the out-degree of \c pto and the in-degree of \c pfrom.
 *
 * \example examples/simple/igraph_get_eid.c
 *
 * Added in version 0.2.</para><para>
 */

igraph_error_t igraph_get_eid(const igraph_t *graph, igraph_integer_t *eid,
                   igraph_integer_t pfrom, igraph_integer_t pto,
                   igraph_bool_t directed, igraph_bool_t error) {

    igraph_integer_t from = pfrom, to = pto;
    igraph_integer_t nov = igraph_vcount(graph);

    if (from < 0 || to < 0 || from > nov - 1 || to > nov - 1) {
        IGRAPH_ERROR("cannot get edge ID", IGRAPH_EINVVID);
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

static igraph_error_t igraph_i_get_eids_pairs(const igraph_t *graph, igraph_vector_int_t *eids,
                                     const igraph_vector_int_t *pairs,
                                     igraph_bool_t directed, igraph_bool_t error);

static igraph_error_t igraph_i_get_eids_path(const igraph_t *graph, igraph_vector_int_t *eids,
                                    const igraph_vector_int_t *path,
                                    igraph_bool_t directed, igraph_bool_t error);

static igraph_error_t igraph_i_get_eids_pairs(
    const igraph_t *graph, igraph_vector_int_t *eids,
    const igraph_vector_int_t *pairs, igraph_bool_t directed,
    igraph_bool_t error
) {
    igraph_integer_t n = igraph_vector_int_size(pairs);
    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t i;
    igraph_integer_t eid = -1;

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

static igraph_error_t igraph_i_get_eids_path(
    const igraph_t *graph, igraph_vector_int_t *eids,
    const igraph_vector_int_t *path, igraph_bool_t directed,
    igraph_bool_t error
) {

    igraph_integer_t n = igraph_vector_int_size(path);
    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t i;
    igraph_integer_t eid = -1;

    if (!igraph_vector_int_isininterval(path, 0, no_of_nodes - 1)) {
        IGRAPH_ERROR("Cannot get edge IDs, invalid vertex ID", IGRAPH_EINVVID);
    }

    IGRAPH_CHECK(igraph_vector_int_resize(eids, n == 0 ? 0 : n - 1));

    if (igraph_is_directed(graph)) {
        for (i = 0; i < n - 1; i++) {
            igraph_integer_t from = VECTOR(*path)[i];
            igraph_integer_t to = VECTOR(*path)[i + 1];

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
        for (i = 0; i < n - 1; i++) {
            igraph_integer_t from = VECTOR(*path)[i];
            igraph_integer_t to = VECTOR(*path)[i + 1];

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

/**
 * \function igraph_get_eids
 * Return edge IDs based on the adjacent vertices.
 *
 * This function operates in two modes. If the \c pairs argument is
 * not a null pointer, but the \c path argument is, then it searches
 * for the edge IDs of all pairs of vertices given in \c pairs. The
 * pairs of vertex IDs are taken consecutively from the vector,
 * i.e. <code>VECTOR(pairs)[0]</code> and
 * <code>VECTOR(pairs)[1]</code> give the first
 * pair, <code>VECTOR(pairs)[2]</code> and
 * <code>VECTOR(pairs)[3]</code> the second pair, etc.
 *
 * </para><para>
 * If the \c pairs argument is a null pointer, and \c path is not a
 * null pointer, then the \c path is interpreted as a path given by
 * vertex IDs and the edges along the path are returned.
 *
 * </para><para>
 * If neither \c pairs nor \c path are null pointers, then both are
 * considered (first \c pairs and then \c path), and the results are
 * concatenated.
 *
 * </para><para>
 * If the \c error argument is true, then it is an error to give pairs
 * of vertices that are not connected. Otherwise -1 is
 * reported for not connected vertices.
 *
 * </para><para>
 * If there are multiple edges in the graph, then these are ignored;
 * i.e. for a given pair of vertex IDs, always the same edge ID is
 * returned, even if the pair is given multiple time in \c pairs or in
 * \c path. See \ref igraph_get_eids_multi() for a similar function
 * that works differently in case of multiple edges.
 *
 * \param graph The input graph.
 * \param eids Pointer to an initialized vector, the result is stored
 *        here. It will be resized as needed.
 * \param pairs Vector giving pairs of vertices, or a null pointer.
 * \param path Vector giving vertex IDs along a path, or a null
 *        pointer.
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
 * \sa \ref igraph_get_eid() for a single edge, \ref
 * igraph_get_eids_multi() for a version that handles multiple edges
 * better (at a cost).
 *
 * \example examples/simple/igraph_get_eids.c
 */

igraph_error_t igraph_get_eids(const igraph_t *graph, igraph_vector_int_t *eids,
                    const igraph_vector_int_t *pairs,
                    const igraph_vector_int_t *path,
                    igraph_bool_t directed, igraph_bool_t error) {

    if (!pairs && !path) {
        igraph_vector_int_clear(eids);
        return IGRAPH_SUCCESS;
    } else if (pairs && !path) {
        return igraph_i_get_eids_pairs(graph, eids, pairs, directed, error);
    } else if (!pairs && path) {
        return igraph_i_get_eids_path(graph, eids, path, directed, error);
    } else {
        /* both */
        igraph_vector_int_t tmp;
        IGRAPH_VECTOR_INT_INIT_FINALLY(&tmp, 0);
        IGRAPH_CHECK(igraph_i_get_eids_pairs(graph, eids, pairs, directed, error));
        IGRAPH_CHECK(igraph_i_get_eids_path(graph, &tmp, path, directed, error));
        IGRAPH_CHECK(igraph_vector_int_append(eids, &tmp));
        igraph_vector_int_destroy(&tmp);
        IGRAPH_FINALLY_CLEAN(1);
        return IGRAPH_SUCCESS;
    }
}

#undef BINSEARCH
#undef FIND_DIRECTED_EDGE
#undef FIND_UNDIRECTED_EDGE

#define BINSEARCH(start,end,value,iindex,edgelist,N,pos,seen)    \
    do {                                                      \
        while ((start) < (end)) {                                 \
            igraph_integer_t mid=(start)+((end)-(start))/2;                 \
            igraph_integer_t e= VECTOR((iindex))[mid];        \
            if (VECTOR((edgelist))[e] < (value)) {                  \
                (start)=mid+1;                                        \
            } else {                                                \
                (end)=mid;                                            \
            }                                                       \
        }                                                         \
        if ((start)<(N)) {                                        \
            igraph_integer_t e = VECTOR((iindex))[(start)];       \
            while ((start)<(N) && seen[e] && VECTOR(edgelist)[e] == (value)) {  \
                (start)++;                        \
                e= VECTOR(iindex)[(start)];         \
            }                                           \
            if ((start)<(N) && !(seen[e]) && VECTOR(edgelist)[e] == (value)) {  \
                *(pos)=e;                  \
            }                                                       \
        } } while(0)

#define FIND_DIRECTED_EDGE(graph,xfrom,xto,eid,seen)            \
    do {                                                              \
        igraph_integer_t start = VECTOR(graph->os)[xfrom];         \
        igraph_integer_t end= VECTOR(graph->os)[xfrom+1];         \
        igraph_integer_t N=end;                                                 \
        igraph_integer_t start2= VECTOR(graph->is)[xto];          \
        igraph_integer_t end2= VECTOR(graph->is)[xto+1];          \
        igraph_integer_t N2=end2;                                               \
        if (end-start<end2-start2) {                                    \
            BINSEARCH(start,end,xto,graph->oi,graph->to,N,eid,seen);      \
        } else {                                                        \
            BINSEARCH(start2,end2,xfrom,graph->ii,graph->from,N2,eid,seen);   \
        }                                                               \
    } while (0)

#define FIND_UNDIRECTED_EDGE(graph,from,to,eid,seen)            \
    do {                                                              \
        igraph_integer_t xfrom1= from > to ? from : to;                         \
        igraph_integer_t xto1= from > to ? to : from;                           \
        FIND_DIRECTED_EDGE(graph,xfrom1,xto1,eid,seen);         \
    } while (0)


igraph_error_t igraph_i_get_eids_multipairs(const igraph_t *graph, igraph_vector_int_t *eids,
                                          const igraph_vector_int_t *pairs,
                                          igraph_bool_t directed, igraph_bool_t error);

igraph_error_t igraph_i_get_eids_multipath(const igraph_t *graph, igraph_vector_int_t *eids,
                              const igraph_vector_int_t *path,
                              igraph_bool_t directed, igraph_bool_t error);

igraph_error_t igraph_i_get_eids_multipairs(const igraph_t *graph, igraph_vector_int_t *eids,
                                          const igraph_vector_int_t *pairs,
                                          igraph_bool_t directed, igraph_bool_t error) {

    igraph_integer_t i, n = igraph_vector_int_size(pairs);
    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_bool_t *seen;
    igraph_integer_t eid = -1;

    if (n % 2 != 0) {
        IGRAPH_ERROR("Cannot get edge IDs, invalid length of edge IDs",
                     IGRAPH_EINVAL);
    }
    if (!igraph_vector_int_isininterval(pairs, 0, no_of_nodes - 1)) {
        IGRAPH_ERROR("Cannot get edge IDs, invalid vertex ID", IGRAPH_EINVVID);
    }

    seen = IGRAPH_CALLOC(no_of_edges, igraph_bool_t);
    if (seen == 0) {
        IGRAPH_ERROR("Cannot get edge IDs", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, seen);
    IGRAPH_CHECK(igraph_vector_int_resize(eids, n / 2));

    if (igraph_is_directed(graph)) {
        for (i = 0; i < n / 2; i++) {
            igraph_integer_t from = VECTOR(*pairs)[2 * i];
            igraph_integer_t to = VECTOR(*pairs)[2 * i + 1];

            eid = -1;
            FIND_DIRECTED_EDGE(graph, from, to, &eid, seen);
            if (!directed && eid < 0) {
                FIND_DIRECTED_EDGE(graph, to, from, &eid, seen);
            }

            VECTOR(*eids)[i] = eid;
            if (eid >= 0) {
                seen[(eid)] = 1;
            } else if (error) {
                IGRAPH_ERROR("Cannot get edge ID, no such edge", IGRAPH_EINVAL);
            }
        }
    } else {
        for (i = 0; i < n / 2; i++) {
            igraph_integer_t from = VECTOR(*pairs)[2 * i];
            igraph_integer_t to = VECTOR(*pairs)[2 * i + 1];

            eid = -1;
            FIND_UNDIRECTED_EDGE(graph, from, to, &eid, seen);
            VECTOR(*eids)[i] = eid;
            if (eid >= 0) {
                seen[(eid)] = 1;
            } else if (error) {
                IGRAPH_ERROR("Cannot get edge ID, no such edge", IGRAPH_EINVAL);
            }
        }
    }

    IGRAPH_FREE(seen);
    IGRAPH_FINALLY_CLEAN(1);
    return IGRAPH_SUCCESS;
}

igraph_error_t igraph_i_get_eids_multipath(const igraph_t *graph, igraph_vector_int_t *eids,
                                           const igraph_vector_int_t *path,
                                           igraph_bool_t directed, igraph_bool_t error) {

    igraph_integer_t n = igraph_vector_int_size(path);
    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_bool_t *seen;
    igraph_integer_t i;
    igraph_integer_t eid = -1;

    if (!igraph_vector_int_isininterval(path, 0, no_of_nodes - 1)) {
        IGRAPH_ERROR("Cannot get edge IDs, invalid vertex ID", IGRAPH_EINVVID);
    }

    seen = IGRAPH_CALLOC(no_of_edges, igraph_bool_t);
    if (!seen) {
        IGRAPH_ERROR("Cannot get edge IDs", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, seen);
    IGRAPH_CHECK(igraph_vector_int_resize(eids, n == 0 ? 0 : n - 1));

    if (igraph_is_directed(graph)) {
        for (i = 0; i < n - 1; i++) {
            igraph_integer_t from = VECTOR(*path)[i];
            igraph_integer_t to = VECTOR(*path)[i + 1];

            eid = -1;
            FIND_DIRECTED_EDGE(graph, from, to, &eid, seen);
            if (!directed && eid < 0) {
                FIND_DIRECTED_EDGE(graph, to, from, &eid, seen);
            }

            VECTOR(*eids)[i] = eid;
            if (eid >= 0) {
                seen[(eid)] = 1;
            } else if (error) {
                IGRAPH_ERROR("Cannot get edge ID, no such edge", IGRAPH_EINVAL);
            }
        }
    } else {
        for (i = 0; i < n - 1; i++) {
            igraph_integer_t from = VECTOR(*path)[i];
            igraph_integer_t to = VECTOR(*path)[i + 1];

            eid = -1;
            FIND_UNDIRECTED_EDGE(graph, from, to, &eid, seen);
            VECTOR(*eids)[i] = eid;
            if (eid >= 0) {
                seen[(eid)] = 1;
            } else if (error) {
                IGRAPH_ERROR("Cannot get edge ID, no such edge", IGRAPH_EINVAL);
            }
        }
    }

    IGRAPH_FREE(seen);
    IGRAPH_FINALLY_CLEAN(1);
    return IGRAPH_SUCCESS;
}

#undef BINSEARCH
#undef FIND_DIRECTED_EDGE
#undef FIND_UNDIRECTED_EDGE

/**
 * \function igraph_get_eids_multi
 * \brief Query edge IDs based on their adjacent vertices, handle multiple edges.
 *
 * This function operates in two modes. If the \c pairs argument is
 * not a null pointer, but the \c path argument is, then it searches
 * for the edge IDs of all pairs of vertices given in \c pairs. The
 * pairs of vertex IDs are taken consecutively from the vector,
 * i.e. <code>VECTOR(pairs)[0]</code> and
 * <code>VECTOR(pairs)[1]</code> give the first pair,
 * <code>VECTOR(pairs)[2]</code> and <code>VECTOR(pairs)[3]</code> the
 * second pair, etc.
 *
 * </para><para>
 * If the \c pairs argument is a null pointer, and \c path is not a
 * null pointer, then the \c path is interpreted as a path given by
 * vertex IDs and the edges along the path are returned.
 *
 * </para><para>
 * If the \c error argument is true, then it is an error to give pairs of
 * vertices that are not connected. Otherwise -1 is
 * returned for not connected vertex pairs.
 *
 * </para><para>
 * An error is triggered if both \c pairs and \c path are non-null
 * pointers.
 *
 * </para><para>
 * This function handles multiple edges properly, i.e. if the same
 * pair is given multiple times and they are indeed connected by
 * multiple edges, then each time a different edge ID is reported.
 *
 * \param graph The input graph.
 * \param eids Pointer to an initialized vector, the result is stored
 *        here. It will be resized as needed.
 * \param pairs Vector giving pairs of vertices, or a null pointer.
 * \param path Vector giving vertex IDs along a path, or a null
 *        pointer.
 * \param directed Logical scalar, whether to consider edge directions
 *        in directed graphs. This is ignored for undirected graphs.
 * \param error Logical scalar, whether to report an error if
 *        non-connected vertices are specified. If false, then -1
 *        is returned for non-connected vertex pairs.
 * \return Error code.
 *
 * Time complexity: O(|E|+n log(d)), where |E| is the number of edges
 * in the graph, n is the number of queried edges and d is the average
 * degree of the vertices.
 *
 * \sa \ref igraph_get_eid() for a single edge, \ref
 * igraph_get_eids() for a faster version that does not handle
 * multiple edges.
 */

igraph_error_t igraph_get_eids_multi(const igraph_t *graph, igraph_vector_int_t *eids,
                          const igraph_vector_int_t *pairs,
                          const igraph_vector_int_t *path,
                          igraph_bool_t directed, igraph_bool_t error) {

    if (!pairs && !path) {
        igraph_vector_int_clear(eids);
        return IGRAPH_SUCCESS;
    } else if (pairs && !path) {
        return igraph_i_get_eids_multipairs(graph, eids, pairs, directed, error);
    } else if (!pairs && path) {
        return igraph_i_get_eids_multipath(graph, eids, path, directed, error);
    } else { /* both */
        IGRAPH_ERROR("Give `pairs' or `path' but not both", IGRAPH_EINVAL);
    }
}

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
        return igraph_i_incident(graph, eids, pnode, mode, IGRAPH_LOOPS_TWICE, IGRAPH_MULTIPLE);
    } else {
        return igraph_i_incident(graph, eids, pnode, mode, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE);
    }
}

igraph_error_t igraph_i_incident(const igraph_t *graph, igraph_vector_int_t *eids, igraph_integer_t pnode,
        igraph_neimode_t mode, igraph_loops_t loops, igraph_multiple_t multiple) {
#define DEDUPLICATE_IF_NEEDED(vertex, n)                                                 \
    if (should_filter_duplicates) {                                                        \
        if ((loops == IGRAPH_NO_LOOPS && vertex == pnode) ||                               \
                (loops == IGRAPH_LOOPS_ONCE && vertex == pnode && last_added == pnode) ||  \
                (multiple == IGRAPH_NO_MULTIPLE && vertex == last_added)) {                \
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

    IGRAPH_CHECK(igraph_vector_int_resize(eids, length));

    /* The loops below produce an ordering what is consistent with the
     * ordering returned by igraph_neighbors(), and this should be preserved.
     * We are dealing with two sorted lists; one for the successors and one
     * for the predecessors. If we have requested only one of them, we have
     * an easy job. If we have requested both, we need to merge the two lists
     * to ensure that the output is sorted by the vertex IDs of the "other"
     * endpoint of the affected edges */
    if (!igraph_is_directed(graph) || mode != IGRAPH_ALL) {
        /* We did not ask for both directions; this is the easy case */

        should_filter_duplicates = !(multiple == IGRAPH_MULTIPLE &&
                ((!igraph_is_directed(graph) && loops == IGRAPH_LOOPS_TWICE) ||
                 (igraph_is_directed(graph) && loops != IGRAPH_NO_LOOPS)));

        if (mode & IGRAPH_OUT) {
            j = VECTOR(graph->os)[node + 1];
            for (i = VECTOR(graph->os)[node]; i < j; i++) {
                igraph_integer_t edge = VECTOR(graph->oi)[i];
                igraph_integer_t other = VECTOR(graph->to)[edge];
                DEDUPLICATE_IF_NEEDED(other, 1);
                VECTOR(*eids)[idx++] = edge;
            }
        }
        if (mode & IGRAPH_IN) {
            j = VECTOR(graph->is)[node + 1];
            for (i = VECTOR(graph->is)[node]; i < j; i++) {
                igraph_integer_t edge = VECTOR(graph->ii)[i];
                igraph_integer_t other = VECTOR(graph->from)[edge];
                DEDUPLICATE_IF_NEEDED(other, 1);
                VECTOR(*eids)[idx++] = edge;
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
                VECTOR(*eids)[idx++] = eid1;
            } else if (n1 > n2) {
                i2++;
                DEDUPLICATE_IF_NEEDED(n2, 1);
                VECTOR(*eids)[idx++] = eid2;
            } else {
                i1++;
                i2++;
                DEDUPLICATE_IF_NEEDED(n2, 2);
                VECTOR(*eids)[idx++] = eid1;
                if (should_filter_duplicates && ((loops == IGRAPH_LOOPS_ONCE && n1 == pnode && last_added == pnode) ||
                        (multiple == IGRAPH_NO_MULTIPLE))) {
                    length--;
                    continue;
                }
                VECTOR(*eids)[idx++] = eid2;
            }
        }

        while (i1 < j1) {
            eid1 = VECTOR(graph->oi)[i1++];
            igraph_integer_t to = VECTOR(graph->to)[eid1];
            DEDUPLICATE_IF_NEEDED(to, 1);
            VECTOR(*eids)[idx++] = eid1;
        }

        while (i2 < j2) {
            eid2 = VECTOR(graph->ii)[i2++];
            igraph_integer_t from = VECTOR(graph->from)[eid2];
            DEDUPLICATE_IF_NEEDED(from, 1);
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
 * \sa igraph_isomorphic() to test if two graphs are isomorphic.
 */

igraph_error_t igraph_is_same_graph(const igraph_t *graph1, const igraph_t *graph2, igraph_bool_t *res) {
    igraph_integer_t nv1 = igraph_vcount(graph1);
    igraph_integer_t nv2 = igraph_vcount(graph2);
    igraph_integer_t ne1 = igraph_ecount(graph1);
    igraph_integer_t ne2 = igraph_ecount(graph2);
    igraph_integer_t i, eid1, eid2;

    *res = 0; /* Assume that the graphs differ */

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

    *res = 1; /* No difference was found, graphs are the same */
    return IGRAPH_SUCCESS;
}
