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
#include "igraph_attributes.h"
#include "igraph_memory.h"
#include "config.h"

/* Internal functions */

static int igraph_i_create_start(
        igraph_vector_t *res, igraph_vector_t *el,
        igraph_vector_t *index, igraph_integer_t nodes);

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
int igraph_empty(igraph_t *graph, igraph_integer_t n, igraph_bool_t directed) {
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
int igraph_empty_attrs(igraph_t *graph, igraph_integer_t n, igraph_bool_t directed, void* attr) {

    if (n < 0) {
        IGRAPH_ERROR("cannot create empty graph with negative number of vertices",
                     IGRAPH_EINVAL);
    }

    if (!IGRAPH_FINITE(n)) {
        IGRAPH_ERROR("number of vertices is not finite (NA, NaN or Inf)", IGRAPH_EINVAL);
    }

    graph->n = 0;
    graph->directed = directed;
    IGRAPH_VECTOR_INIT_FINALLY(&graph->from, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&graph->to, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&graph->oi, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&graph->ii, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&graph->os, 1);
    IGRAPH_VECTOR_INIT_FINALLY(&graph->is, 1);

    VECTOR(graph->os)[0] = 0;
    VECTOR(graph->is)[0] = 0;

    /* init attributes */
    graph->attr = 0;
    IGRAPH_CHECK(igraph_i_attribute_init(graph, attr));

    /* add the vertices */
    IGRAPH_CHECK(igraph_add_vertices(graph, n, 0));

    IGRAPH_FINALLY_CLEAN(6);
    return 0;
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

    igraph_vector_destroy(&graph->from);
    igraph_vector_destroy(&graph->to);
    igraph_vector_destroy(&graph->oi);
    igraph_vector_destroy(&graph->ii);
    igraph_vector_destroy(&graph->os);
    igraph_vector_destroy(&graph->is);
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

int igraph_copy(igraph_t *to, const igraph_t *from) {
    to->n = from->n;
    to->directed = from->directed;
    IGRAPH_CHECK(igraph_vector_copy(&to->from, &from->from));
    IGRAPH_FINALLY(igraph_vector_destroy, &to->from);
    IGRAPH_CHECK(igraph_vector_copy(&to->to, &from->to));
    IGRAPH_FINALLY(igraph_vector_destroy, &to->to);
    IGRAPH_CHECK(igraph_vector_copy(&to->oi, &from->oi));
    IGRAPH_FINALLY(igraph_vector_destroy, &to->oi);
    IGRAPH_CHECK(igraph_vector_copy(&to->ii, &from->ii));
    IGRAPH_FINALLY(igraph_vector_destroy, &to->ii);
    IGRAPH_CHECK(igraph_vector_copy(&to->os, &from->os));
    IGRAPH_FINALLY(igraph_vector_destroy, &to->os);
    IGRAPH_CHECK(igraph_vector_copy(&to->is, &from->is));
    IGRAPH_FINALLY(igraph_vector_destroy, &to->is);

    IGRAPH_I_ATTRIBUTE_COPY(to, from, 1, 1, 1); /* does IGRAPH_CHECK */

    IGRAPH_FINALLY_CLEAN(6);
    return 0;
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
 *    invalid vertex id in edges vector.
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
int igraph_add_edges(igraph_t *graph, const igraph_vector_t *edges,
                     void *attr) {
    long int no_of_edges = igraph_vector_size(&graph->from);
    long int edges_to_add = igraph_vector_size(edges) / 2;
    long int i = 0;
    igraph_error_handler_t *oldhandler;
    int ret1, ret2;
    igraph_vector_t newoi, newii;
    igraph_bool_t directed = igraph_is_directed(graph);

    if (igraph_vector_size(edges) % 2 != 0) {
        IGRAPH_ERROR("invalid (odd) length of edges vector", IGRAPH_EINVEVECTOR);
    }
    if (!igraph_vector_isininterval(edges, 0, igraph_vcount(graph) - 1)) {
        IGRAPH_ERROR("cannot add edges", IGRAPH_EINVVID);
    }

    /* from & to */
    IGRAPH_CHECK(igraph_vector_reserve(&graph->from, no_of_edges + edges_to_add));
    IGRAPH_CHECK(igraph_vector_reserve(&graph->to, no_of_edges + edges_to_add));

    while (i < edges_to_add * 2) {
        if (directed || VECTOR(*edges)[i] > VECTOR(*edges)[i + 1]) {
            igraph_vector_push_back(&graph->from, VECTOR(*edges)[i++]); /* reserved */
            igraph_vector_push_back(&graph->to,   VECTOR(*edges)[i++]); /* reserved */
        } else {
            igraph_vector_push_back(&graph->to,   VECTOR(*edges)[i++]); /* reserved */
            igraph_vector_push_back(&graph->from, VECTOR(*edges)[i++]); /* reserved */
        }
    }

    /* disable the error handler temporarily */
    oldhandler = igraph_set_error_handler(igraph_error_handler_ignore);

    /* oi & ii */
    ret1 = igraph_vector_init(&newoi, no_of_edges);
    ret2 = igraph_vector_init(&newii, no_of_edges);
    if (ret1 != 0 || ret2 != 0) {
        igraph_vector_resize(&graph->from, no_of_edges); /* gets smaller */
        igraph_vector_resize(&graph->to, no_of_edges);   /* gets smaller */
        igraph_set_error_handler(oldhandler);
        IGRAPH_ERROR("cannot add edges", IGRAPH_ERROR_SELECT_2(ret1, ret2));
    }
    ret1 = igraph_vector_order(&graph->from, &graph->to, &newoi, graph->n);
    ret2 = igraph_vector_order(&graph->to, &graph->from, &newii, graph->n);
    if (ret1 != 0 || ret2 != 0) {
        igraph_vector_resize(&graph->from, no_of_edges);
        igraph_vector_resize(&graph->to, no_of_edges);
        igraph_vector_destroy(&newoi);
        igraph_vector_destroy(&newii);
        igraph_set_error_handler(oldhandler);
        IGRAPH_ERROR("cannot add edges", IGRAPH_ERROR_SELECT_2(ret1, ret2));
    }

    /* Attributes */
    if (graph->attr) {
        igraph_set_error_handler(oldhandler);
        ret1 = igraph_i_attribute_add_edges(graph, edges, attr);
        igraph_set_error_handler(igraph_error_handler_ignore);
        if (ret1 != 0) {
            igraph_vector_resize(&graph->from, no_of_edges);
            igraph_vector_resize(&graph->to, no_of_edges);
            igraph_vector_destroy(&newoi);
            igraph_vector_destroy(&newii);
            igraph_set_error_handler(oldhandler);
            IGRAPH_ERROR("cannot add edges", ret1);
        }
    }

    /* os & is, its length does not change, error safe */
    igraph_i_create_start(&graph->os, &graph->from, &newoi, graph->n);
    igraph_i_create_start(&graph->is, &graph->to, &newii, graph->n);

    /* everything went fine  */
    igraph_vector_destroy(&graph->oi);
    igraph_vector_destroy(&graph->ii);
    graph->oi = newoi;
    graph->ii = newii;
    igraph_set_error_handler(oldhandler);

    return 0;
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
int igraph_add_vertices(igraph_t *graph, igraph_integer_t nv, void *attr) {
    long int ec = igraph_ecount(graph);
    long int i;

    if (nv < 0) {
        IGRAPH_ERROR("cannot add negative number of vertices", IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_vector_reserve(&graph->os, graph->n + nv + 1));
    IGRAPH_CHECK(igraph_vector_reserve(&graph->is, graph->n + nv + 1));

    igraph_vector_resize(&graph->os, graph->n + nv + 1); /* reserved */
    igraph_vector_resize(&graph->is, graph->n + nv + 1); /* reserved */
    for (i = graph->n + 1; i < graph->n + nv + 1; i++) {
        VECTOR(graph->os)[i] = ec;
        VECTOR(graph->is)[i] = ec;
    }

    graph->n += nv;

    if (graph->attr) {
        IGRAPH_CHECK(igraph_i_attribute_add_vertices(graph, nv, attr));
    }

    return 0;
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
int igraph_delete_edges(igraph_t *graph, igraph_es_t edges) {
    long int no_of_edges = igraph_ecount(graph);
    long int no_of_nodes = igraph_vcount(graph);
    long int edges_to_remove = 0;
    long int remaining_edges;
    igraph_eit_t eit;

    igraph_vector_t newfrom, newto, newoi;

    int *mark;
    long int i, j;

    mark = igraph_Calloc(no_of_edges, int);
    if (mark == 0) {
        IGRAPH_ERROR("Cannot delete edges", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, mark);

    IGRAPH_CHECK(igraph_eit_create(graph, edges, &eit));
    IGRAPH_FINALLY(igraph_eit_destroy, &eit);

    for (IGRAPH_EIT_RESET(eit); !IGRAPH_EIT_END(eit); IGRAPH_EIT_NEXT(eit)) {
        long int e = IGRAPH_EIT_GET(eit);
        if (mark[e] == 0) {
            edges_to_remove++;
            mark[e]++;
        }
    }
    remaining_edges = no_of_edges - edges_to_remove;

    /* We don't need the iterator any more */
    igraph_eit_destroy(&eit);
    IGRAPH_FINALLY_CLEAN(1);

    IGRAPH_VECTOR_INIT_FINALLY(&newfrom, remaining_edges);
    IGRAPH_VECTOR_INIT_FINALLY(&newto, remaining_edges);

    /* Actually remove the edges, move from pos i to pos j in newfrom/newto */
    for (i = 0, j = 0; j < remaining_edges; i++) {
        if (mark[i] == 0) {
            VECTOR(newfrom)[j] = VECTOR(graph->from)[i];
            VECTOR(newto)[j] = VECTOR(graph->to)[i];
            j++;
        }
    }

    /* Create index, this might require additional memory */
    IGRAPH_VECTOR_INIT_FINALLY(&newoi, remaining_edges);
    IGRAPH_CHECK(igraph_vector_order(&newfrom, &newto, &newoi, no_of_nodes));
    IGRAPH_CHECK(igraph_vector_order(&newto, &newfrom, &graph->ii, no_of_nodes));

    /* Edge attributes, we need an index that gives the ids of the
       original edges for every new edge.
    */
    if (graph->attr) {
        igraph_vector_t idx;
        IGRAPH_VECTOR_INIT_FINALLY(&idx, remaining_edges);
        for (i = 0, j = 0; i < no_of_edges; i++) {
            if (mark[i] == 0) {
                VECTOR(idx)[j++] = i;
            }
        }
        IGRAPH_CHECK(igraph_i_attribute_permute_edges(graph, graph, &idx));
        igraph_vector_destroy(&idx);
        IGRAPH_FINALLY_CLEAN(1);
    }

    /* Ok, we've all memory needed, free the old structure  */
    igraph_vector_destroy(&graph->from);
    igraph_vector_destroy(&graph->to);
    igraph_vector_destroy(&graph->oi);
    graph->from = newfrom;
    graph->to = newto;
    graph->oi = newoi;
    IGRAPH_FINALLY_CLEAN(3);

    igraph_Free(mark);
    IGRAPH_FINALLY_CLEAN(1);

    /* Create start vectors, no memory is needed for this */
    igraph_i_create_start(&graph->os, &graph->from, &graph->oi,
                          (igraph_integer_t) no_of_nodes);
    igraph_i_create_start(&graph->is, &graph->to,   &graph->ii,
                          (igraph_integer_t) no_of_nodes);

    /* Nothing to deallocate... */
    return 0;
}

/**
 * \ingroup interface
 * \function igraph_delete_vertices
 * \brief Removes vertices (with all their edges) from the graph.
 *
 * </para><para>
 * This function changes the ids of the vertices (except in some very
 * special cases, but these should not be relied on anyway).
 *
 * </para><para>
 * This function invalidates all iterators.
 *
 * \param graph The graph to work on.
 * \param vertices The ids of the vertices to remove in a
 *                 vector. The vector may contain the same id more
 *                 than once.
 * \return Error code:
 *         \c IGRAPH_EINVVID: invalid vertex id.
 *
 * Time complexity: O(|V|+|E|),
 * |V| and
 * |E| are the number of vertices and
 * edges in the original graph.
 *
 * \example examples/simple/igraph_delete_vertices.c
 */
int igraph_delete_vertices(igraph_t *graph, const igraph_vs_t vertices) {
    return igraph_delete_vertices_idx(graph, vertices, /* idx= */ 0,
                                      /* invidx= */ 0);
}

int igraph_delete_vertices_idx(igraph_t *graph, const igraph_vs_t vertices,
                               igraph_vector_t *idx,
                               igraph_vector_t *invidx) {

    long int no_of_edges = igraph_ecount(graph);
    long int no_of_nodes = igraph_vcount(graph);
    igraph_vector_t edge_recoding, vertex_recoding;
    igraph_vector_t *my_vertex_recoding = &vertex_recoding;
    igraph_vit_t vit;
    igraph_t newgraph;
    long int i, j;
    long int remaining_vertices, remaining_edges;

    if (idx) {
        my_vertex_recoding = idx;
        IGRAPH_CHECK(igraph_vector_resize(idx, no_of_nodes));
        igraph_vector_null(idx);
    } else {
        IGRAPH_VECTOR_INIT_FINALLY(&vertex_recoding, no_of_nodes);
    }

    IGRAPH_VECTOR_INIT_FINALLY(&edge_recoding, no_of_edges);

    IGRAPH_CHECK(igraph_vit_create(graph, vertices, &vit));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit);

    /* mark the vertices to delete */
    for (; !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit) ) {
        long int vertex = IGRAPH_VIT_GET(vit);
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
        long int from = (long int) VECTOR(graph->from)[i];
        long int to = (long int) VECTOR(graph->to)[i];
        if (VECTOR(*my_vertex_recoding)[from] != 0 &&
            VECTOR(*my_vertex_recoding)[to  ] != 0) {
            VECTOR(edge_recoding)[i] = remaining_edges + 1;
            remaining_edges++;
        }
    }

    /* start creating the graph */
    newgraph.n = (igraph_integer_t) remaining_vertices;
    newgraph.directed = graph->directed;

    /* allocate vectors */
    IGRAPH_VECTOR_INIT_FINALLY(&newgraph.from, remaining_edges);
    IGRAPH_VECTOR_INIT_FINALLY(&newgraph.to, remaining_edges);
    IGRAPH_VECTOR_INIT_FINALLY(&newgraph.oi, remaining_edges);
    IGRAPH_VECTOR_INIT_FINALLY(&newgraph.ii, remaining_edges);
    IGRAPH_VECTOR_INIT_FINALLY(&newgraph.os, remaining_vertices + 1);
    IGRAPH_VECTOR_INIT_FINALLY(&newgraph.is, remaining_vertices + 1);

    /* Add the edges */
    for (i = 0, j = 0; j < remaining_edges; i++) {
        if (VECTOR(edge_recoding)[i] > 0) {
            long int from = (long int) VECTOR(graph->from)[i];
            long int to = (long int) VECTOR(graph->to  )[i];
            VECTOR(newgraph.from)[j] = VECTOR(*my_vertex_recoding)[from] - 1;
            VECTOR(newgraph.to  )[j] = VECTOR(*my_vertex_recoding)[to] - 1;
            j++;
        }
    }
    /* update oi & ii */
    IGRAPH_CHECK(igraph_vector_order(&newgraph.from, &newgraph.to, &newgraph.oi,
                                     remaining_vertices));
    IGRAPH_CHECK(igraph_vector_order(&newgraph.to, &newgraph.from, &newgraph.ii,
                                     remaining_vertices));

    IGRAPH_CHECK(igraph_i_create_start(&newgraph.os, &newgraph.from,
                                       &newgraph.oi, (igraph_integer_t)
                                       remaining_vertices));
    IGRAPH_CHECK(igraph_i_create_start(&newgraph.is, &newgraph.to,
                                       &newgraph.ii, (igraph_integer_t)
                                       remaining_vertices));

    /* attributes */
    IGRAPH_I_ATTRIBUTE_COPY(&newgraph, graph,
                            /*graph=*/ 1, /*vertex=*/0, /*edge=*/0);
    IGRAPH_FINALLY_CLEAN(6);
    IGRAPH_FINALLY(igraph_destroy, &newgraph);

    if (newgraph.attr) {
        igraph_vector_t iidx;
        IGRAPH_VECTOR_INIT_FINALLY(&iidx, remaining_vertices);
        for (i = 0; i < no_of_nodes; i++) {
            long int jj = (long int) VECTOR(*my_vertex_recoding)[i];
            if (jj != 0) {
                VECTOR(iidx)[ jj - 1 ] = i;
            }
        }
        IGRAPH_CHECK(igraph_i_attribute_permute_vertices(graph,
                     &newgraph,
                     &iidx));
        IGRAPH_CHECK(igraph_vector_resize(&iidx, remaining_edges));
        for (i = 0; i < no_of_edges; i++) {
            long int jj = (long int) VECTOR(edge_recoding)[i];
            if (jj != 0) {
                VECTOR(iidx)[ jj - 1 ] = i;
            }
        }
        IGRAPH_CHECK(igraph_i_attribute_permute_edges(graph, &newgraph, &iidx));
        igraph_vector_destroy(&iidx);
        IGRAPH_FINALLY_CLEAN(1);
    }

    igraph_vit_destroy(&vit);
    igraph_vector_destroy(&edge_recoding);
    igraph_destroy(graph);
    *graph = newgraph;

    IGRAPH_FINALLY_CLEAN(3);

    /* TODO: this is duplicate */
    if (invidx) {
        IGRAPH_CHECK(igraph_vector_resize(invidx, remaining_vertices));
        for (i = 0; i < no_of_nodes; i++) {
            long int newid = (long int) VECTOR(*my_vertex_recoding)[i];
            if (newid != 0) {
                VECTOR(*invidx)[newid - 1] = i;
            }
        }
    }

    if (!idx) {
        igraph_vector_destroy(my_vertex_recoding);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return 0;
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
    return (igraph_integer_t) igraph_vector_size(&graph->from);
}

/**
 * \ingroup interface
 * \function igraph_neighbors
 * \brief Adjacent vertices to a vertex.
 *
 * \param graph The graph to work on.
 * \param neis This vector will contain the result. The vector should
 *        be initialized beforehand and will be resized. Starting from igraph
 *        version 0.4 this vector is always sorted, the vertex ids are
 *        in increasing order.
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
 *         \c IGRAPH_EINVVID: invalid vertex id.
 *         \c IGRAPH_EINVMODE: invalid mode argument.
 *         \c IGRAPH_ENOMEM: not enough memory.
 *
 * Time complexity: O(d),
 * d is the number
 * of adjacent vertices to the queried vertex.
 *
 * \example examples/simple/igraph_neighbors.c
 */
int igraph_neighbors(const igraph_t *graph, igraph_vector_t *neis, igraph_integer_t pnode,
                     igraph_neimode_t mode) {

    long int length = 0, idx = 0;
    long int i, j;

    long int node = pnode;

    if (node < 0 || node > igraph_vcount(graph) - 1) {
        IGRAPH_ERROR("cannot get neighbors", IGRAPH_EINVVID);
    }
    if (mode != IGRAPH_OUT && mode != IGRAPH_IN &&
        mode != IGRAPH_ALL) {
        IGRAPH_ERROR("cannot get neighbors", IGRAPH_EINVMODE);
    }

    if (! graph->directed) {
        mode = IGRAPH_ALL;
    }

    /* Calculate needed space first & allocate it*/

    if (mode & IGRAPH_OUT) {
        length += (VECTOR(graph->os)[node + 1] - VECTOR(graph->os)[node]);
    }
    if (mode & IGRAPH_IN) {
        length += (VECTOR(graph->is)[node + 1] - VECTOR(graph->is)[node]);
    }

    IGRAPH_CHECK(igraph_vector_resize(neis, length));

    if (!igraph_is_directed(graph) || mode != IGRAPH_ALL) {

        if (mode & IGRAPH_OUT) {
            j = (long int) VECTOR(graph->os)[node + 1];
            for (i = (long int) VECTOR(graph->os)[node]; i < j; i++) {
                VECTOR(*neis)[idx++] =
                    VECTOR(graph->to)[ (long int)VECTOR(graph->oi)[i] ];
            }
        }
        if (mode & IGRAPH_IN) {
            j = (long int) VECTOR(graph->is)[node + 1];
            for (i = (long int) VECTOR(graph->is)[node]; i < j; i++) {
                VECTOR(*neis)[idx++] =
                    VECTOR(graph->from)[ (long int)VECTOR(graph->ii)[i] ];
            }
        }
    } else {
        /* both in- and out- neighbors in a directed graph,
           we need to merge the two 'vectors' */
        long int jj1 = (long int) VECTOR(graph->os)[node + 1];
        long int j2 = (long int) VECTOR(graph->is)[node + 1];
        long int i1 = (long int) VECTOR(graph->os)[node];
        long int i2 = (long int) VECTOR(graph->is)[node];
        while (i1 < jj1 && i2 < j2) {
            long int n1 = (long int) VECTOR(graph->to)[
                   (long int)VECTOR(graph->oi)[i1] ];
            long int n2 = (long int) VECTOR(graph->from)[
                   (long int)VECTOR(graph->ii)[i2] ];
            if (n1 < n2) {
                VECTOR(*neis)[idx++] = n1;
                i1++;
            } else if (n1 > n2) {
                VECTOR(*neis)[idx++] = n2;
                i2++;
            } else {
                VECTOR(*neis)[idx++] = n1;
                VECTOR(*neis)[idx++] = n2;
                i1++;
                i2++;
            }
        }
        while (i1 < jj1) {
            long int n1 = (long int) VECTOR(graph->to)[
                   (long int)VECTOR(graph->oi)[i1] ];
            VECTOR(*neis)[idx++] = n1;
            i1++;
        }
        while (i2 < j2) {
            long int n2 = (long int) VECTOR(graph->from)[
                   (long int)VECTOR(graph->ii)[i2] ];
            VECTOR(*neis)[idx++] = n2;
            i2++;
        }
    }

    return 0;
}

/**
 * \ingroup internal
 *
 */

static int igraph_i_create_start(
        igraph_vector_t *res, igraph_vector_t *el,
        igraph_vector_t *iindex, igraph_integer_t nodes) {

# define EDGE(i) (VECTOR(*el)[ (long int) VECTOR(*iindex)[(i)] ])

    long int no_of_nodes;
    long int no_of_edges;
    long int i, j, idx;

    no_of_nodes = nodes;
    no_of_edges = igraph_vector_size(el);

    /* result */

    IGRAPH_CHECK(igraph_vector_resize(res, nodes + 1));

    /* create the index */

    if (igraph_vector_size(el) == 0) {
        /* empty graph */
        igraph_vector_null(res);
    } else {
        idx = -1;
        for (i = 0; i <= EDGE(0); i++) {
            idx++; VECTOR(*res)[idx] = 0;
        }
        for (i = 1; i < no_of_edges; i++) {
            long int n = (long int) (EDGE(i) - EDGE((long int)VECTOR(*res)[idx]));
            for (j = 0; j < n; j++) {
                idx++; VECTOR(*res)[idx] = i;
            }
        }
        j = (long int) EDGE((long int)VECTOR(*res)[idx]);
        for (i = 0; i < no_of_nodes - j; i++) {
            idx++; VECTOR(*res)[idx] = no_of_edges;
        }
    }

    /* clean */

# undef EDGE
    return 0;
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
 * \param vids Vector, giving the vertex ids of which the degree will
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
 *         \c IGRAPH_EINVVID: invalid vertex id.
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
int igraph_degree(const igraph_t *graph, igraph_vector_t *res,
                  const igraph_vs_t vids,
                  igraph_neimode_t mode, igraph_bool_t loops) {

    long int nodes_to_calc;
    long int i, j;
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

    IGRAPH_CHECK(igraph_vector_resize(res, nodes_to_calc));
    igraph_vector_null(res);

    if (loops) {
        if (mode & IGRAPH_OUT) {
            for (IGRAPH_VIT_RESET(vit), i = 0;
                 !IGRAPH_VIT_END(vit);
                 IGRAPH_VIT_NEXT(vit), i++) {
                long int vid = IGRAPH_VIT_GET(vit);
                VECTOR(*res)[i] += (VECTOR(graph->os)[vid + 1] - VECTOR(graph->os)[vid]);
            }
        }
        if (mode & IGRAPH_IN) {
            for (IGRAPH_VIT_RESET(vit), i = 0;
                 !IGRAPH_VIT_END(vit);
                 IGRAPH_VIT_NEXT(vit), i++) {
                long int vid = IGRAPH_VIT_GET(vit);
                VECTOR(*res)[i] += (VECTOR(graph->is)[vid + 1] - VECTOR(graph->is)[vid]);
            }
        }
    } else { /* no loops */
        if (mode & IGRAPH_OUT) {
            for (IGRAPH_VIT_RESET(vit), i = 0;
                 !IGRAPH_VIT_END(vit);
                 IGRAPH_VIT_NEXT(vit), i++) {
                long int vid = IGRAPH_VIT_GET(vit);
                VECTOR(*res)[i] += (VECTOR(graph->os)[vid + 1] - VECTOR(graph->os)[vid]);
                for (j = (long int) VECTOR(graph->os)[vid];
                     j < VECTOR(graph->os)[vid + 1]; j++) {
                    if (VECTOR(graph->to)[ (long int)VECTOR(graph->oi)[j] ] == vid) {
                        VECTOR(*res)[i] -= 1;
                    }
                }
            }
        }
        if (mode & IGRAPH_IN) {
            for (IGRAPH_VIT_RESET(vit), i = 0;
                 !IGRAPH_VIT_END(vit);
                 IGRAPH_VIT_NEXT(vit), i++) {
                long int vid = IGRAPH_VIT_GET(vit);
                VECTOR(*res)[i] += (VECTOR(graph->is)[vid + 1] - VECTOR(graph->is)[vid]);
                for (j = (long int) VECTOR(graph->is)[vid];
                     j < VECTOR(graph->is)[vid + 1]; j++) {
                    if (VECTOR(graph->from)[ (long int)VECTOR(graph->ii)[j] ] == vid) {
                        VECTOR(*res)[i] -= 1;
                    }
                }
            }
        }
    }  /* loops */

    igraph_vit_destroy(&vit);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

/**
 * \function igraph_edge
 * \brief Gives the head and tail vertices of an edge.
 *
 * \param graph The graph object.
 * \param eid The edge id.
 * \param from Pointer to an \type igraph_integer_t. The tail of the edge
 * will be placed here.
 * \param to Pointer to an \type igraph_integer_t. The head of the edge
 * will be placed here.
 * \return Error code. The current implementation always returns with
 * success.
 * \sa \ref igraph_get_eid() for the opposite operation.
 *
 * Added in version 0.2.</para><para>
 *
 * Time complexity: O(1).
 */

int igraph_edge(const igraph_t *graph, igraph_integer_t eid,
                igraph_integer_t *from, igraph_integer_t *to) {

    if (igraph_is_directed(graph)) {
        *from = (igraph_integer_t) VECTOR(graph->from)[(long int)eid];
        *to   = (igraph_integer_t) VECTOR(graph->to  )[(long int)eid];
    } else {
        *from = (igraph_integer_t) VECTOR(graph->to  )[(long int)eid];
        *to   = (igraph_integer_t) VECTOR(graph->from)[(long int)eid];
    }

    return 0;
}

int igraph_edges(const igraph_t *graph, igraph_es_t eids,
                 igraph_vector_t *edges) {

    igraph_eit_t eit;
    long int n, ptr = 0;

    IGRAPH_CHECK(igraph_eit_create(graph, eids, &eit));
    IGRAPH_FINALLY(igraph_eit_destroy, &eit);
    n = IGRAPH_EIT_SIZE(eit);
    IGRAPH_CHECK(igraph_vector_resize(edges, n * 2));
    if (igraph_is_directed(graph)) {
        for (; !IGRAPH_EIT_END(eit); IGRAPH_EIT_NEXT(eit)) {
            long int e = IGRAPH_EIT_GET(eit);
            VECTOR(*edges)[ptr++] = IGRAPH_FROM(graph, e);
            VECTOR(*edges)[ptr++] = IGRAPH_TO(graph, e);
        }
    } else {
        for (; !IGRAPH_EIT_END(eit); IGRAPH_EIT_NEXT(eit)) {
            long int e = IGRAPH_EIT_GET(eit);
            VECTOR(*edges)[ptr++] = IGRAPH_TO(graph, e);
            VECTOR(*edges)[ptr++] = IGRAPH_FROM(graph, e);
        }
    }

    igraph_eit_destroy(&eit);
    IGRAPH_FINALLY_CLEAN(1);
    return 0;
}

/* This is an unsafe macro. Only supply variable names, i.e. no
   expressions as parameters, otherwise nasty things can happen */

#define BINSEARCH(start,end,value,iindex,edgelist,N,pos)     \
    do {                                                      \
        while ((start) < (end)) {                                 \
            long int mid=(start)+((end)-(start))/2;                 \
            long int e=(long int) VECTOR((iindex))[mid];            \
            if (VECTOR((edgelist))[e] < (value)) {                  \
                (start)=mid+1;                                        \
            } else {                                                \
                (end)=mid;                                            \
            }                                                       \
        }                                                         \
        if ((start)<(N)) {                                        \
            long int e=(long int) VECTOR((iindex))[(start)];        \
            if (VECTOR((edgelist))[e] == (value)) {                 \
                *(pos)=(igraph_integer_t) e;              \
            }                                                       \
        } } while(0)

#define FIND_DIRECTED_EDGE(graph,xfrom,xto,eid)                     \
    do {                                                              \
        long int start=(long int) VECTOR(graph->os)[xfrom];         \
        long int end=(long int) VECTOR(graph->os)[xfrom+1];         \
        long int N=end;                                                 \
        long int start2=(long int) VECTOR(graph->is)[xto];          \
        long int end2=(long int) VECTOR(graph->is)[xto+1];          \
        long int N2=end2;                                               \
        if (end-start<end2-start2) {                                    \
            BINSEARCH(start,end,xto,graph->oi,graph->to,N,eid);           \
        } else {                                                        \
            BINSEARCH(start2,end2,xfrom,graph->ii,graph->from,N2,eid);    \
        }                                                               \
    } while (0)

#define FIND_UNDIRECTED_EDGE(graph,from,to,eid)                     \
    do {                                                              \
        long int xfrom1= from > to ? from : to;                         \
        long int xto1= from > to ? to : from;                           \
        FIND_DIRECTED_EDGE(graph,xfrom1,xto1,eid);                      \
    } while (0)

/**
 * \function igraph_get_eid
 * \brief Get the edge id from the end points of an edge.
 *
 * For undirected graphs \c pfrom and \c pto are exchangeable.
 *
 * \param graph The graph object.
 * \param eid Pointer to an integer, the edge id will be stored here.
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

int igraph_get_eid(const igraph_t *graph, igraph_integer_t *eid,
                   igraph_integer_t pfrom, igraph_integer_t pto,
                   igraph_bool_t directed, igraph_bool_t error) {

    long int from = pfrom, to = pto;
    long int nov = igraph_vcount(graph);

    if (from < 0 || to < 0 || from > nov - 1 || to > nov - 1) {
        IGRAPH_ERROR("cannot get edge id", IGRAPH_EINVVID);
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
            IGRAPH_ERROR("Cannot get edge id, no such edge", IGRAPH_EINVAL);
        }
    }

    return IGRAPH_SUCCESS;
}

int igraph_get_eids_pairs(const igraph_t *graph, igraph_vector_t *eids,
                          const igraph_vector_t *pairs,
                          igraph_bool_t directed, igraph_bool_t error);

int igraph_get_eids_path(const igraph_t *graph, igraph_vector_t *eids,
                         const igraph_vector_t *path,
                         igraph_bool_t directed, igraph_bool_t error);

int igraph_get_eids_pairs(const igraph_t *graph, igraph_vector_t *eids,
                          const igraph_vector_t *pairs,
                          igraph_bool_t directed, igraph_bool_t error) {
    long int n = igraph_vector_size(pairs);
    long int no_of_nodes = igraph_vcount(graph);
    long int i;
    igraph_integer_t eid = -1;

    if (n % 2 != 0) {
        IGRAPH_ERROR("Cannot get edge ids, invalid length of edge ids",
                     IGRAPH_EINVAL);
    }
    if (!igraph_vector_isininterval(pairs, 0, no_of_nodes - 1)) {
        IGRAPH_ERROR("Cannot get edge ids, invalid vertex id", IGRAPH_EINVVID);
    }

    IGRAPH_CHECK(igraph_vector_resize(eids, n / 2));

    if (igraph_is_directed(graph)) {
        for (i = 0; i < n / 2; i++) {
            long int from = (long int) VECTOR(*pairs)[2 * i];
            long int to = (long int) VECTOR(*pairs)[2 * i + 1];

            eid = -1;
            FIND_DIRECTED_EDGE(graph, from, to, &eid);
            if (!directed && eid < 0) {
                FIND_DIRECTED_EDGE(graph, to, from, &eid);
            }

            VECTOR(*eids)[i] = eid;
            if (eid < 0 && error) {
                IGRAPH_ERROR("Cannot get edge id, no such edge", IGRAPH_EINVAL);
            }
        }
    } else {
        for (i = 0; i < n / 2; i++) {
            long int from = (long int) VECTOR(*pairs)[2 * i];
            long int to = (long int) VECTOR(*pairs)[2 * i + 1];

            eid = -1;
            FIND_UNDIRECTED_EDGE(graph, from, to, &eid);
            VECTOR(*eids)[i] = eid;
            if (eid < 0 && error) {
                IGRAPH_ERROR("Cannot get edge id, no such edge", IGRAPH_EINVAL);
            }
        }
    }

    return 0;
}

int igraph_get_eids_path(const igraph_t *graph, igraph_vector_t *eids,
                         const igraph_vector_t *path,
                         igraph_bool_t directed, igraph_bool_t error) {

    long int n = igraph_vector_size(path);
    long int no_of_nodes = igraph_vcount(graph);
    long int i;
    igraph_integer_t eid = -1;

    if (!igraph_vector_isininterval(path, 0, no_of_nodes - 1)) {
        IGRAPH_ERROR("Cannot get edge ids, invalid vertex id", IGRAPH_EINVVID);
    }

    IGRAPH_CHECK(igraph_vector_resize(eids, n == 0 ? 0 : n - 1));

    if (igraph_is_directed(graph)) {
        for (i = 0; i < n - 1; i++) {
            long int from = (long int) VECTOR(*path)[i];
            long int to = (long int) VECTOR(*path)[i + 1];

            eid = -1;
            FIND_DIRECTED_EDGE(graph, from, to, &eid);
            if (!directed && eid < 0) {
                FIND_DIRECTED_EDGE(graph, to, from, &eid);
            }

            VECTOR(*eids)[i] = eid;
            if (eid < 0 && error) {
                IGRAPH_ERROR("Cannot get edge id, no such edge", IGRAPH_EINVAL);
            }
        }
    } else {
        for (i = 0; i < n - 1; i++) {
            long int from = (long int) VECTOR(*path)[i];
            long int to = (long int) VECTOR(*path)[i + 1];

            eid = -1;
            FIND_UNDIRECTED_EDGE(graph, from, to, &eid);
            VECTOR(*eids)[i] = eid;
            if (eid < 0 && error) {
                IGRAPH_ERROR("Cannot get edge id, no such edge", IGRAPH_EINVAL);
            }
        }
    }

    return 0;
}

/**
 * \function igraph_get_eids
 * Return edge ids based on the adjacent vertices.
 *
 * This function operates in two modes. If the \c pairs argument is
 * not a null pointer, but the \c path argument is, then it searches
 * for the edge ids of all pairs of vertices given in \c pairs. The
 * pairs of vertex ids are taken consecutively from the vector,
 * i.e. <code>VECTOR(pairs)[0]</code> and
 * <code>VECTOR(pairs)[1]</code> give the first
 * pair, <code>VECTOR(pairs)[2]</code> and
 * <code>VECTOR(pairs)[3]</code> the second pair, etc.
 *
 * </para><para>
 * If the \c pairs argument is a null pointer, and \c path is not a
 * null pointer, then the \c path is interpreted as a path given by
 * vertex ids and the edges along the path are returned.
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
 * i.e. for a given pair of vertex ids, always the same edge id is
 * returned, even if the pair is given multiple time in \c pairs or in
 * \c path. See \ref igraph_get_eids_multi() for a similar function
 * that works differently in case of multiple edges.
 *
 * \param graph The input graph.
 * \param eids Pointer to an initialized vector, the result is stored
 *        here. It will be resized as needed.
 * \param pairs Vector giving pairs of vertices, or a null pointer.
 * \param path Vector giving vertex ids along a path, or a null
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

int igraph_get_eids(const igraph_t *graph, igraph_vector_t *eids,
                    const igraph_vector_t *pairs,
                    const igraph_vector_t *path,
                    igraph_bool_t directed, igraph_bool_t error) {

    if (!pairs && !path) {
        igraph_vector_clear(eids);
        return 0;
    } else if (pairs && !path) {
        return igraph_get_eids_pairs(graph, eids, pairs, directed, error);
    } else if (!pairs && path) {
        return igraph_get_eids_path(graph, eids, path, directed, error);
    } else {
        /* both */
        igraph_vector_t tmp;
        IGRAPH_VECTOR_INIT_FINALLY(&tmp, 0);
        IGRAPH_CHECK(igraph_get_eids_pairs(graph, eids, pairs, directed, error));
        IGRAPH_CHECK(igraph_get_eids_path(graph, &tmp, path, directed, error));
        IGRAPH_CHECK(igraph_vector_append(eids, &tmp));
        igraph_vector_destroy(&tmp);
        IGRAPH_FINALLY_CLEAN(1);
        return 0;
    }
}

#undef BINSEARCH
#undef FIND_DIRECTED_EDGE
#undef FIND_UNDIRECTED_EDGE

#define BINSEARCH(start,end,value,iindex,edgelist,N,pos,seen)    \
    do {                                                      \
        while ((start) < (end)) {                                 \
            long int mid=(start)+((end)-(start))/2;                 \
            long int e=(long int) VECTOR((iindex))[mid];        \
            if (VECTOR((edgelist))[e] < (value)) {                  \
                (start)=mid+1;                                        \
            } else {                                                \
                (end)=mid;                                            \
            }                                                       \
        }                                                         \
        if ((start)<(N)) {                                        \
            long int e=(long int) VECTOR((iindex))[(start)];        \
            while ((start)<(N) && seen[e] && VECTOR(edgelist)[e] == (value)) {  \
                (start)++;                        \
                e=(long int) VECTOR(iindex)[(start)];         \
            }                                           \
            if ((start)<(N) && !(seen[e]) && VECTOR(edgelist)[e] == (value)) {  \
                *(pos)=(igraph_integer_t) e;                  \
            }                                                       \
        } } while(0)

#define FIND_DIRECTED_EDGE(graph,xfrom,xto,eid,seen)            \
    do {                                                              \
        long int start=(long int) VECTOR(graph->os)[xfrom];         \
        long int end=(long int) VECTOR(graph->os)[xfrom+1];         \
        long int N=end;                                                 \
        long int start2=(long int) VECTOR(graph->is)[xto];          \
        long int end2=(long int) VECTOR(graph->is)[xto+1];          \
        long int N2=end2;                                               \
        if (end-start<end2-start2) {                                    \
            BINSEARCH(start,end,xto,graph->oi,graph->to,N,eid,seen);      \
        } else {                                                        \
            BINSEARCH(start2,end2,xfrom,graph->ii,graph->from,N2,eid,seen);   \
        }                                                               \
    } while (0)

#define FIND_UNDIRECTED_EDGE(graph,from,to,eid,seen)            \
    do {                                                              \
        long int xfrom1= from > to ? from : to;                         \
        long int xto1= from > to ? to : from;                           \
        FIND_DIRECTED_EDGE(graph,xfrom1,xto1,eid,seen);         \
    } while (0)


int igraph_get_eids_multipairs(const igraph_t *graph, igraph_vector_t *eids,
                               const igraph_vector_t *pairs,
                               igraph_bool_t directed, igraph_bool_t error);

int igraph_get_eids_multipath(const igraph_t *graph, igraph_vector_t *eids,
                              const igraph_vector_t *path,
                              igraph_bool_t directed, igraph_bool_t error);

int igraph_get_eids_multipairs(const igraph_t *graph, igraph_vector_t *eids,
                               const igraph_vector_t *pairs,
                               igraph_bool_t directed, igraph_bool_t error) {

    long int n = igraph_vector_size(pairs);
    long int no_of_nodes = igraph_vcount(graph);
    long int no_of_edges = igraph_ecount(graph);
    igraph_bool_t *seen;
    long int i;
    igraph_integer_t eid = -1;

    if (n % 2 != 0) {
        IGRAPH_ERROR("Cannot get edge ids, invalid length of edge ids",
                     IGRAPH_EINVAL);
    }
    if (!igraph_vector_isininterval(pairs, 0, no_of_nodes - 1)) {
        IGRAPH_ERROR("Cannot get edge ids, invalid vertex id", IGRAPH_EINVVID);
    }

    seen = igraph_Calloc(no_of_edges, igraph_bool_t);
    if (seen == 0) {
        IGRAPH_ERROR("Cannot get edge ids", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, seen);
    IGRAPH_CHECK(igraph_vector_resize(eids, n / 2));

    if (igraph_is_directed(graph)) {
        for (i = 0; i < n / 2; i++) {
            long int from = (long int) VECTOR(*pairs)[2 * i];
            long int to = (long int) VECTOR(*pairs)[2 * i + 1];

            eid = -1;
            FIND_DIRECTED_EDGE(graph, from, to, &eid, seen);
            if (!directed && eid < 0) {
                FIND_DIRECTED_EDGE(graph, to, from, &eid, seen);
            }

            VECTOR(*eids)[i] = eid;
            if (eid >= 0) {
                seen[(long int)(eid)] = 1;
            } else if (error) {
                IGRAPH_ERROR("Cannot get edge id, no such edge", IGRAPH_EINVAL);
            }
        }
    } else {
        for (i = 0; i < n / 2; i++) {
            long int from = (long int) VECTOR(*pairs)[2 * i];
            long int to = (long int) VECTOR(*pairs)[2 * i + 1];

            eid = -1;
            FIND_UNDIRECTED_EDGE(graph, from, to, &eid, seen);
            VECTOR(*eids)[i] = eid;
            if (eid >= 0) {
                seen[(long int)(eid)] = 1;
            } else if (error) {
                IGRAPH_ERROR("Cannot get edge id, no such edge", IGRAPH_EINVAL);
            }
        }
    }

    igraph_Free(seen);
    IGRAPH_FINALLY_CLEAN(1);
    return 0;
}

int igraph_get_eids_multipath(const igraph_t *graph, igraph_vector_t *eids,
                              const igraph_vector_t *path,
                              igraph_bool_t directed, igraph_bool_t error) {

    long int n = igraph_vector_size(path);
    long int no_of_nodes = igraph_vcount(graph);
    long int no_of_edges = igraph_ecount(graph);
    igraph_bool_t *seen;
    long int i;
    igraph_integer_t eid = -1;

    if (!igraph_vector_isininterval(path, 0, no_of_nodes - 1)) {
        IGRAPH_ERROR("Cannot get edge ids, invalid vertex id", IGRAPH_EINVVID);
    }

    seen = igraph_Calloc(no_of_edges, igraph_bool_t);
    if (!seen) {
        IGRAPH_ERROR("Cannot get edge ids", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, seen);
    IGRAPH_CHECK(igraph_vector_resize(eids, n == 0 ? 0 : n - 1));

    if (igraph_is_directed(graph)) {
        for (i = 0; i < n - 1; i++) {
            long int from = (long int) VECTOR(*path)[i];
            long int to = (long int) VECTOR(*path)[i + 1];

            eid = -1;
            FIND_DIRECTED_EDGE(graph, from, to, &eid, seen);
            if (!directed && eid < 0) {
                FIND_DIRECTED_EDGE(graph, to, from, &eid, seen);
            }

            VECTOR(*eids)[i] = eid;
            if (eid >= 0) {
                seen[(long int)(eid)] = 1;
            } else if (error) {
                IGRAPH_ERROR("Cannot get edge id, no such edge", IGRAPH_EINVAL);
            }
        }
    } else {
        for (i = 0; i < n - 1; i++) {
            long int from = (long int) VECTOR(*path)[i];
            long int to = (long int) VECTOR(*path)[i + 1];

            eid = -1;
            FIND_UNDIRECTED_EDGE(graph, from, to, &eid, seen);
            VECTOR(*eids)[i] = eid;
            if (eid >= 0) {
                seen[(long int)(eid)] = 1;
            } else if (error) {
                IGRAPH_ERROR("Cannot get edge id, no such edge", IGRAPH_EINVAL);
            }
        }
    }

    igraph_Free(seen);
    IGRAPH_FINALLY_CLEAN(1);
    return 0;
}

#undef BINSEARCH
#undef FIND_DIRECTED_EDGE
#undef FIND_UNDIRECTED_EDGE

/**
 * \function igraph_get_eids_multi
 * \brief Query edge ids based on their adjacent vertices, handle multiple edges.
 *
 * This function operates in two modes. If the \c pairs argument is
 * not a null pointer, but the \c path argument is, then it searches
 * for the edge ids of all pairs of vertices given in \c pairs. The
 * pairs of vertex ids are taken consecutively from the vector,
 * i.e. <code>VECTOR(pairs)[0]</code> and
 * <code>VECTOR(pairs)[1]</code> give the first pair,
 * <code>VECTOR(pairs)[2]</code> and <code>VECTOR(pairs)[3]</code> the
 * second pair, etc.
 *
 * </para><para>
 * If the \c pairs argument is a null pointer, and \c path is not a
 * null pointer, then the \c path is interpreted as a path given by
 * vertex ids and the edges along the path are returned.
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
 * multiple edges, then each time a different edge id is reported.
 *
 * \param graph The input graph.
 * \param eids Pointer to an initialized vector, the result is stored
 *        here. It will be resized as needed.
 * \param pairs Vector giving pairs of vertices, or a null pointer.
 * \param path Vector giving vertex ids along a path, or a null
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

int igraph_get_eids_multi(const igraph_t *graph, igraph_vector_t *eids,
                          const igraph_vector_t *pairs,
                          const igraph_vector_t *path,
                          igraph_bool_t directed, igraph_bool_t error) {

    if (!pairs && !path) {
        igraph_vector_clear(eids);
        return 0;
    } else if (pairs && !path) {
        return igraph_get_eids_multipairs(graph, eids, pairs, directed, error);
    } else if (!pairs && path) {
        return igraph_get_eids_multipath(graph, eids, path, directed, error);
    } else { /* both */
        IGRAPH_ERROR("Give `pairs' or `path' but not both", IGRAPH_EINVAL);
    }
}

/**
 * \function igraph_adjacent
 * \brief Gives the incident edges of a vertex.
 *
 * This function was superseded by \ref igraph_incident() in igraph 0.6.
 * Please use \ref igraph_incident() instead of this function.
 *
 * </para><para>
 * Added in version 0.2, deprecated in version 0.6.
 */
int igraph_adjacent(const igraph_t *graph, igraph_vector_t *eids,
                    igraph_integer_t pnode, igraph_neimode_t mode) {
    IGRAPH_WARNING("igraph_adjacent is deprecated, use igraph_incident");
    return igraph_incident(graph, eids, pnode, mode);
}

/**
 * \function igraph_incident
 * \brief Gives the incident edges of a vertex.
 *
 * \param graph The graph object.
 * \param eids An initialized \type vector_t object. It will be resized
 * to hold the result.
 * \param pnode A vertex id.
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

int igraph_incident(const igraph_t *graph, igraph_vector_t *eids,
                    igraph_integer_t pnode, igraph_neimode_t mode) {

    long int length = 0, idx = 0;
    long int i, j;

    long int node = pnode;

    if (node < 0 || node > igraph_vcount(graph) - 1) {
        IGRAPH_ERROR("cannot get neighbors", IGRAPH_EINVVID);
    }
    if (mode != IGRAPH_OUT && mode != IGRAPH_IN &&
        mode != IGRAPH_ALL) {
        IGRAPH_ERROR("cannot get neighbors", IGRAPH_EINVMODE);
    }

    if (! graph->directed) {
        mode = IGRAPH_ALL;
    }

    /* Calculate needed space first & allocate it*/

    if (mode & IGRAPH_OUT) {
        length += (VECTOR(graph->os)[node + 1] - VECTOR(graph->os)[node]);
    }
    if (mode & IGRAPH_IN) {
        length += (VECTOR(graph->is)[node + 1] - VECTOR(graph->is)[node]);
    }

    IGRAPH_CHECK(igraph_vector_resize(eids, length));

    if (mode & IGRAPH_OUT) {
        j = (long int) VECTOR(graph->os)[node + 1];
        for (i = (long int) VECTOR(graph->os)[node]; i < j; i++) {
            VECTOR(*eids)[idx++] = VECTOR(graph->oi)[i];
        }
    }
    if (mode & IGRAPH_IN) {
        j = (long int) VECTOR(graph->is)[node + 1];
        for (i = (long int) VECTOR(graph->is)[node]; i < j; i++) {
            VECTOR(*eids)[idx++] = VECTOR(graph->ii)[i];
        }
    }

    return 0;
}
