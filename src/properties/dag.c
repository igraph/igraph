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

#include "igraph_topology.h"

#include "igraph_constructors.h"
#include "igraph_dqueue.h"
#include "igraph_interface.h"
#include "igraph_stack.h"

/**
 * \function igraph_topological_sorting
 * \brief Calculate a possible topological sorting of the graph.
 *
 * </para><para>
 * A topological sorting of a directed acyclic graph (DAG) is a linear ordering
 * of its vertices where each vertex comes before all nodes to which it has
 * edges. Every DAG has at least one topological sort, and may have many.
 * This function returns one possible topological sort among them. If the
 * graph contains any cycles that are not self-loops, an error is raised.
 *
 * \param graph The input graph.
 * \param res Pointer to a vector, the result will be stored here.
 *   It will be resized if needed.
 * \param mode Specifies how to use the direction of the edges.
 *   For \c IGRAPH_OUT, the sorting order ensures that each vertex comes
 *   before all vertices to which it has edges, so vertices with no incoming
 *   edges go first. For \c IGRAPH_IN, it is quite the opposite: each
 *   vertex comes before all vertices from which it receives edges. Vertices
 *   with no outgoing edges go first.
 * \return Error code.
 *
 * Time complexity: O(|V|+|E|), where |V| and |E| are the number of
 * vertices and edges in the original input graph.
 *
 * \sa \ref igraph_is_dag() if you are only interested in whether a given
 *     graph is a DAG or not, or \ref igraph_feedback_arc_set() to find a
 *     set of edges whose removal makes the graph acyclic.
 *
 * \example examples/simple/igraph_topological_sorting.c
 */
igraph_error_t igraph_topological_sorting(
        const igraph_t* graph, igraph_vector_int_t *res, igraph_neimode_t mode) {

    /* Note: This function ignores self-loops, there it cannot
     * use the IGRAPH_PROP_IS_DAG property cache entry. */

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_vector_int_t degrees;
    igraph_vector_int_t neis;
    igraph_dqueue_int_t sources;
    igraph_neimode_t deg_mode;
    igraph_integer_t node, i, j;

    if (mode == IGRAPH_ALL || !igraph_is_directed(graph)) {
        IGRAPH_ERROR("Topological sorting does not make sense for undirected graphs.",
                     IGRAPH_EINVAL);
    } else if (mode == IGRAPH_OUT) {
        deg_mode = IGRAPH_IN;
    } else if (mode == IGRAPH_IN) {
        deg_mode = IGRAPH_OUT;
    } else {
        IGRAPH_ERROR("Invalid mode for topological sorting.", IGRAPH_EINVMODE);
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&degrees, no_of_nodes);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&neis, 0);
    IGRAPH_CHECK(igraph_dqueue_int_init(&sources, 0));
    IGRAPH_FINALLY(igraph_dqueue_int_destroy, &sources);
    IGRAPH_CHECK(igraph_degree(graph, &degrees, igraph_vss_all(), deg_mode, 0));

    igraph_vector_int_clear(res);

    /* Do we have nodes with no incoming vertices? */
    for (i = 0; i < no_of_nodes; i++) {
        if (VECTOR(degrees)[i] == 0) {
            IGRAPH_CHECK(igraph_dqueue_int_push(&sources, i));
        }
    }

    /* Take all nodes with no incoming vertices and remove them */
    while (!igraph_dqueue_int_empty(&sources)) {
        node = igraph_dqueue_int_pop(&sources);
        /* Add the node to the result vector */
        IGRAPH_CHECK(igraph_vector_int_push_back(res, node));
        /* Exclude the node from further source searches */
        VECTOR(degrees)[node] = -1;
        /* Get the neighbors and decrease their degrees by one */
        IGRAPH_CHECK(igraph_neighbors(graph, &neis, node, mode));
        j = igraph_vector_int_size(&neis);
        for (i = 0; i < j; i++) {
            VECTOR(degrees)[ VECTOR(neis)[i] ]--;
            if (VECTOR(degrees)[ VECTOR(neis)[i] ] == 0) {
                IGRAPH_CHECK(igraph_dqueue_int_push(&sources, VECTOR(neis)[i]));
            }
        }
    }

    igraph_vector_int_destroy(&degrees);
    igraph_vector_int_destroy(&neis);
    igraph_dqueue_int_destroy(&sources);
    IGRAPH_FINALLY_CLEAN(3);

    if (igraph_vector_int_size(res) < no_of_nodes) {
        IGRAPH_ERROR("The graph has cycles; "
                     "topological sorting is only possible in acyclic graphs.",
                     IGRAPH_EINVAL);
    }

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_is_dag
 * \brief Checks whether a graph is a directed acyclic graph (DAG).
 *
 * </para><para>
 * A directed acyclic graph (DAG) is a directed graph with no cycles.
 *
 * </para><para>
 * This function returns false for undirected graphs.
 *
 * </para><para>
 * The return value of this function is cached in the graph itself; calling
 * the function multiple times with no modifications to the graph in between
 * will return a cached value in O(1) time.
 *
 * \param graph The input graph.
 * \param res Pointer to a boolean constant, the result
 *     is stored here.
 * \return Error code.
 *
 * Time complexity: O(|V|+|E|), where |V| and |E| are the number of
 * vertices and edges in the original input graph.
 *
 * \sa \ref igraph_topological_sorting() to get a possible topological
 *     sorting of a DAG.
 */
igraph_error_t igraph_is_dag(const igraph_t* graph, igraph_bool_t *res) {
    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_vector_int_t degrees;
    igraph_vector_int_t neis;
    igraph_dqueue_int_t sources;

    if (!igraph_is_directed(graph)) {
        *res = false;
        return IGRAPH_SUCCESS;
    }

    IGRAPH_RETURN_IF_CACHED_BOOL(graph, IGRAPH_PROP_IS_DAG, res);

    IGRAPH_VECTOR_INT_INIT_FINALLY(&degrees, no_of_nodes);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&neis, 0);
    IGRAPH_DQUEUE_INT_INIT_FINALLY(&sources, 0);

    IGRAPH_CHECK(igraph_degree(graph, &degrees, igraph_vss_all(), IGRAPH_IN, /* loops */ true));

    igraph_integer_t vertices_left = no_of_nodes;

    /* Do we have nodes with no incoming edges? */
    for (igraph_integer_t i = 0; i < no_of_nodes; i++) {
        if (VECTOR(degrees)[i] == 0) {
            IGRAPH_CHECK(igraph_dqueue_int_push(&sources, i));
        }
    }

    /* Take all nodes with no incoming edges and remove them */
    while (!igraph_dqueue_int_empty(&sources)) {
        igraph_integer_t node = igraph_dqueue_int_pop(&sources);
        /* Exclude the node from further source searches */
        VECTOR(degrees)[node] = -1;
        vertices_left--;
        /* Get the neighbors and decrease their degrees by one */
        IGRAPH_CHECK(igraph_neighbors(graph, &neis, node, IGRAPH_OUT));
        igraph_integer_t n = igraph_vector_int_size(&neis);
        for (igraph_integer_t i = 0; i < n; i++) {
            igraph_integer_t nei = VECTOR(neis)[i];
            if (nei == node) {
                /* Found a self-loop, graph is not a DAG */
                *res = false;
                goto finalize;
            }
            VECTOR(degrees)[nei]--;
            if (VECTOR(degrees)[nei] == 0) {
                IGRAPH_CHECK(igraph_dqueue_int_push(&sources, nei));
            }
        }
    }

    IGRAPH_ASSERT(vertices_left >= 0);
    *res = (vertices_left == 0);

finalize:
    igraph_vector_int_destroy(&degrees);
    igraph_vector_int_destroy(&neis);
    igraph_dqueue_int_destroy(&sources);
    IGRAPH_FINALLY_CLEAN(3);

    igraph_i_property_cache_set_bool_checked(graph, IGRAPH_PROP_IS_DAG, *res);

    return IGRAPH_SUCCESS;
}

/* Create the transitive closure of a tree graph.
   This is fairly simple, we just collect all ancestors of a vertex
   using a depth-first search.
 */
igraph_error_t igraph_transitive_closure_dag(const igraph_t *graph, igraph_t *closure) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_vector_int_t deg;
    igraph_vector_int_t new_edges;
    igraph_vector_int_t ancestors;
    igraph_integer_t root;
    igraph_vector_int_t neighbors;
    igraph_stack_int_t path;
    igraph_vector_bool_t done;

    if (!igraph_is_directed(graph)) {
        IGRAPH_ERROR("Tree transitive closure of a directed graph",
                     IGRAPH_EINVAL);
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&new_edges, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&deg, no_of_nodes);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&ancestors, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&neighbors, 0);
    IGRAPH_CHECK(igraph_stack_int_init(&path, 0));
    IGRAPH_FINALLY(igraph_stack_int_destroy, &path);
    IGRAPH_CHECK(igraph_vector_bool_init(&done, no_of_nodes));
    IGRAPH_FINALLY(igraph_vector_bool_destroy, &done);

    IGRAPH_CHECK(igraph_degree(graph, &deg, igraph_vss_all(),
                               IGRAPH_OUT, IGRAPH_LOOPS));

#define STAR (-1)

    for (root = 0; root < no_of_nodes; root++) {
        if (VECTOR(deg)[root] != 0) {
            continue;
        }
        IGRAPH_CHECK(igraph_stack_int_push(&path, root));

        while (!igraph_stack_int_empty(&path)) {
            igraph_integer_t node = igraph_stack_int_top(&path);
            if (node == STAR) {
                /* Leaving a node */
                igraph_integer_t j, n;
                igraph_stack_int_pop(&path);
                node = igraph_stack_int_pop(&path);
                if (!VECTOR(done)[node]) {
                    igraph_vector_int_pop_back(&ancestors);
                    VECTOR(done)[node] = true;
                }
                n = igraph_vector_int_size(&ancestors);
                for (j = 0; j < n; j++) {
                    IGRAPH_CHECK(igraph_vector_int_push_back(&new_edges, node));
                    IGRAPH_CHECK(igraph_vector_int_push_back(&new_edges,
                                                         VECTOR(ancestors)[j]));
                }
            } else {
                /* Getting into a node */
                igraph_integer_t n, j;
                if (!VECTOR(done)[node]) {
                    IGRAPH_CHECK(igraph_vector_int_push_back(&ancestors, node));
                }
                IGRAPH_CHECK(igraph_neighbors(graph, &neighbors,
                                              node, IGRAPH_IN));
                n = igraph_vector_int_size(&neighbors);
                IGRAPH_CHECK(igraph_stack_int_push(&path, STAR));
                for (j = 0; j < n; j++) {
                    igraph_integer_t nei = VECTOR(neighbors)[j];
                    IGRAPH_CHECK(igraph_stack_int_push(&path, nei));
                }
            }
        }
    }

#undef STAR

    igraph_vector_bool_destroy(&done);
    igraph_stack_int_destroy(&path);
    igraph_vector_int_destroy(&neighbors);
    igraph_vector_int_destroy(&ancestors);
    igraph_vector_int_destroy(&deg);
    IGRAPH_FINALLY_CLEAN(5);

    IGRAPH_CHECK(igraph_create(closure, &new_edges, no_of_nodes,
                               IGRAPH_DIRECTED));

    igraph_vector_int_destroy(&new_edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}
