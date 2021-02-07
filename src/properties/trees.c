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

#include "igraph_adjlist.h"
#include "igraph_constructors.h"
#include "igraph_dqueue.h"
#include "igraph_interface.h"
#include "igraph_stack.h"

/**
 * \function igraph_unfold_tree
 * Unfolding a graph into a tree, by possibly multiplicating its vertices.
 *
 * A graph is converted into a tree (or forest, if it is unconnected),
 * by performing a breadth-first search on it, and replicating
 * vertices that were found a second, third, etc. time.
 * \param graph The input graph, it can be either directed or
 *   undirected.
 * \param tree Pointer to an uninitialized graph object, the result is
 *   stored here.
 * \param mode For directed graphs; whether to follow paths along edge
 *    directions (\c IGRAPH_OUT), or the opposite (\c IGRAPH_IN), or
 *    ignore edge directions completely (\c IGRAPH_ALL). It is ignored
 *    for undirected graphs.
 * \param roots A numeric vector giving the root vertex, or vertices
 *   (if the graph is not connected), to start from.
 * \param vertex_index Pointer to an initialized vector, or a null
 *   pointer. If not a null pointer, then a mapping from the vertices
 *   in the new graph to the ones in the original is created here.
 * \return Error code.
 *
 * Time complexity: O(n+m), linear in the number vertices and edges.
 *
 */
int igraph_unfold_tree(const igraph_t *graph, igraph_t *tree,
                       igraph_neimode_t mode, const igraph_vector_t *roots,
                       igraph_vector_t *vertex_index) {

    long int no_of_nodes = igraph_vcount(graph);
    long int no_of_edges = igraph_ecount(graph);
    long int no_of_roots = igraph_vector_size(roots);
    long int tree_vertex_count = no_of_nodes;

    igraph_vector_t edges;
    igraph_vector_bool_t seen_vertices;
    igraph_vector_bool_t seen_edges;

    igraph_dqueue_t Q;
    igraph_vector_t neis;

    long int i, n, r, v_ptr = no_of_nodes;

    /* TODO: handle not-connected graphs, multiple root vertices */

    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
    igraph_vector_reserve(&edges, no_of_edges * 2);
    IGRAPH_DQUEUE_INIT_FINALLY(&Q, 100);
    IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
    IGRAPH_VECTOR_BOOL_INIT_FINALLY(&seen_vertices, no_of_nodes);
    IGRAPH_VECTOR_BOOL_INIT_FINALLY(&seen_edges, no_of_edges);

    if (vertex_index) {
        IGRAPH_CHECK(igraph_vector_resize(vertex_index, no_of_nodes));
        for (i = 0; i < no_of_nodes; i++) {
            VECTOR(*vertex_index)[i] = i;
        }
    }

    for (r = 0; r < no_of_roots; r++) {

        long int root = (long int) VECTOR(*roots)[r];
        VECTOR(seen_vertices)[root] = 1;
        igraph_dqueue_push(&Q, root);

        while (!igraph_dqueue_empty(&Q)) {
            long int actnode = (long int) igraph_dqueue_pop(&Q);

            IGRAPH_CHECK(igraph_incident(graph, &neis, (igraph_integer_t) actnode, mode));
            n = igraph_vector_size(&neis);
            for (i = 0; i < n; i++) {

                long int edge = (long int) VECTOR(neis)[i];
                long int from = IGRAPH_FROM(graph, edge);
                long int to = IGRAPH_TO(graph, edge);
                long int nei = IGRAPH_OTHER(graph, edge, actnode);

                if (! VECTOR(seen_edges)[edge]) {

                    VECTOR(seen_edges)[edge] = 1;

                    if (! VECTOR(seen_vertices)[nei]) {

                        igraph_vector_push_back(&edges, from);
                        igraph_vector_push_back(&edges, to);

                        VECTOR(seen_vertices)[nei] = 1;
                        IGRAPH_CHECK(igraph_dqueue_push(&Q, nei));

                    } else {

                        tree_vertex_count++;
                        if (vertex_index) {
                            IGRAPH_CHECK(igraph_vector_push_back(vertex_index, nei));
                        }

                        if (from == nei) {
                            igraph_vector_push_back(&edges, v_ptr++);
                            igraph_vector_push_back(&edges, to);
                        } else {
                            igraph_vector_push_back(&edges, from);
                            igraph_vector_push_back(&edges, v_ptr++);
                        }
                    }
                }

            } /* for i<n */

        } /* ! igraph_dqueue_empty(&Q) */

    } /* r < igraph_vector_size(roots) */

    igraph_vector_bool_destroy(&seen_edges);
    igraph_vector_bool_destroy(&seen_vertices);
    igraph_vector_destroy(&neis);
    igraph_dqueue_destroy(&Q);
    IGRAPH_FINALLY_CLEAN(4);

    IGRAPH_CHECK(igraph_create(tree, &edges, tree_vertex_count, igraph_is_directed(graph)));
    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

/* igraph_is_tree -- check if a graph is a tree */

/* count the number of vertices reachable from the root */
static int igraph_i_is_tree_visitor(igraph_integer_t root, const igraph_adjlist_t *al, igraph_integer_t *visited_count) {
    igraph_stack_int_t stack;
    igraph_vector_bool_t visited;
    long i;

    IGRAPH_CHECK(igraph_vector_bool_init(&visited, igraph_adjlist_size(al)));
    IGRAPH_FINALLY(igraph_vector_bool_destroy, &visited);

    IGRAPH_CHECK(igraph_stack_int_init(&stack, 0));
    IGRAPH_FINALLY(igraph_stack_int_destroy, &stack);

    *visited_count = 0;

    /* push the root into the stack */
    IGRAPH_CHECK(igraph_stack_int_push(&stack, root));

    while (! igraph_stack_int_empty(&stack)) {
        igraph_integer_t u;
        igraph_vector_int_t *neighbors;
        long ncount;

        /* take a vertex from the stack, mark it as visited */
        u = igraph_stack_int_pop(&stack);
        if (IGRAPH_LIKELY(! VECTOR(visited)[u])) {
            VECTOR(visited)[u] = 1;
            *visited_count += 1;
        }

        /* register all its yet-unvisited neighbours for future processing */
        neighbors = igraph_adjlist_get(al, u);
        ncount = igraph_vector_int_size(neighbors);
        for (i = 0; i < ncount; ++i) {
            igraph_integer_t v = VECTOR(*neighbors)[i];
            if (! VECTOR(visited)[v]) {
                IGRAPH_CHECK(igraph_stack_int_push(&stack, v));
            }
        }
    }

    igraph_stack_int_destroy(&stack);
    igraph_vector_bool_destroy(&visited);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}


/**
 * \ingroup structural
 * \function igraph_is_tree
 * \brief Decides whether the graph is a tree.
 *
 * An undirected graph is a tree if it is connected and has no cycles.
 * </para><para>
 *
 * In the directed case, a possible additional requirement is that all
 * edges are oriented away from a root (out-tree or arborescence) or all edges
 * are oriented towards a root (in-tree or anti-arborescence).
 * This test can be controlled using the \p mode parameter.
 * </para><para>
 *
 * By convention, the null graph (i.e. the graph with no vertices) is considered not to be a tree.
 *
 * \param graph The graph object to analyze.
 * \param res Pointer to a logical variable, the result will be stored
 *        here.
 * \param root If not \c NULL, the root node will be stored here. When \p mode
 *        is \c IGRAPH_ALL or the graph is undirected, any vertex can be the root
 *        and \p root is set to 0 (the first vertex). When \p mode is \c IGRAPH_OUT
 *        or \c IGRAPH_IN, the root is set to the vertex with zero in- or out-degree,
 *        respectively.
 * \param mode For a directed graph this specifies whether to test for an
 *        out-tree, an in-tree or ignore edge directions. The respective
 *        possible values are:
 *        \c IGRAPH_OUT, \c IGRAPH_IN, \c IGRAPH_ALL. This argument is
 *        ignored for undirected graphs.
 * \return Error code:
 *        \c IGRAPH_EINVAL: invalid mode argument.
 *
 * Time complexity: At most O(|V|+|E|), the
 * number of vertices plus the number of edges in the graph.
 *
 * \sa igraph_is_weakly_connected()
 *
 * \example examples/simple/igraph_tree.c
 */
int igraph_is_tree(const igraph_t *graph, igraph_bool_t *res, igraph_integer_t *root, igraph_neimode_t mode) {
    igraph_adjlist_t al;
    igraph_integer_t iroot = 0;
    igraph_integer_t visited_count;
    igraph_integer_t vcount, ecount;

    vcount = igraph_vcount(graph);
    ecount = igraph_ecount(graph);

    /* A tree must have precisely vcount-1 edges. */
    /* By convention, the zero-vertex graph will not be considered a tree. */
    if (ecount != vcount - 1) {
        *res = 0;
        return IGRAPH_SUCCESS;
    }

    /* The single-vertex graph is a tree, provided it has no edges (checked in the previous if (..)) */
    if (vcount == 1) {
        *res = 1;
        if (root) {
            *root = 0;
        }
        return IGRAPH_SUCCESS;
    }

    /* For higher vertex counts we cannot short-circuit due to the possibility
     * of loops or multi-edges even when the edge count is correct. */

    /* Ignore mode for undirected graphs. */
    if (! igraph_is_directed(graph)) {
        mode = IGRAPH_ALL;
    }

    IGRAPH_CHECK(igraph_adjlist_init(graph, &al, mode, IGRAPH_LOOPS, IGRAPH_MULTIPLE));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &al);

    /* The main algorithm:
     * We find a root and check that all other vertices are reachable from it.
     * We have already checked the number of edges, so with the additional
     * reachability condition we can verify if the graph is a tree.
     *
     * For directed graphs, the root is the node with no incoming/outgoing
     * connections, depending on 'mode'. For undirected, it is arbitrary, so
     * we choose 0.
     */

    *res = 1; /* assume success */

    switch (mode) {
    case IGRAPH_ALL:
        iroot = 0;
        break;

    case IGRAPH_IN:
    case IGRAPH_OUT: {
        igraph_vector_t degree;
        igraph_integer_t i;

        IGRAPH_CHECK(igraph_vector_init(&degree, 0));
        IGRAPH_FINALLY(igraph_vector_destroy, &degree);

        IGRAPH_CHECK(igraph_degree(graph, &degree, igraph_vss_all(), mode == IGRAPH_IN ? IGRAPH_OUT : IGRAPH_IN, /* loops = */ 1));

        for (i = 0; i < vcount; ++i)
            if (VECTOR(degree)[i] == 0) {
                break;
            }

        /* if no suitable root is found, the graph is not a tree */
        if (i == vcount) {
            *res = 0;
        } else {
            iroot = i;
        }

        igraph_vector_destroy(&degree);
        IGRAPH_FINALLY_CLEAN(1);
    }

    break;
    default:
        IGRAPH_ERROR("Invalid mode", IGRAPH_EINVMODE);
    }

    /* if no suitable root was found, skip visiting vertices */
    if (*res) {
        IGRAPH_CHECK(igraph_i_is_tree_visitor(iroot, &al, &visited_count));
        *res = visited_count == vcount;
    }

    if (root) {
        *root = iroot;
    }

    igraph_adjlist_destroy(&al);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}
