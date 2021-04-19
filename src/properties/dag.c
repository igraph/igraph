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
 * graph is not acyclic (it has at least one cycle), an error is raised.
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
int igraph_topological_sorting(const igraph_t* graph, igraph_vector_t *res,
                               igraph_neimode_t mode) {
    long int no_of_nodes = igraph_vcount(graph);
    igraph_vector_t degrees, neis;
    igraph_dqueue_t sources;
    igraph_neimode_t deg_mode;
    long int node, i, j;

    if (mode == IGRAPH_ALL || !igraph_is_directed(graph)) {
        IGRAPH_ERROR("Topological sorting does not make sense for undirected graphs", IGRAPH_EINVAL);
    } else if (mode == IGRAPH_OUT) {
        deg_mode = IGRAPH_IN;
    } else if (mode == IGRAPH_IN) {
        deg_mode = IGRAPH_OUT;
    } else {
        IGRAPH_ERROR("Invalid mode", IGRAPH_EINVAL);
    }

    IGRAPH_VECTOR_INIT_FINALLY(&degrees, no_of_nodes);
    IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
    IGRAPH_CHECK(igraph_dqueue_init(&sources, 0));
    IGRAPH_FINALLY(igraph_dqueue_destroy, &sources);
    IGRAPH_CHECK(igraph_degree(graph, &degrees, igraph_vss_all(), deg_mode, 0));

    igraph_vector_clear(res);

    /* Do we have nodes with no incoming vertices? */
    for (i = 0; i < no_of_nodes; i++) {
        if (VECTOR(degrees)[i] == 0) {
            IGRAPH_CHECK(igraph_dqueue_push(&sources, i));
        }
    }

    /* Take all nodes with no incoming vertices and remove them */
    while (!igraph_dqueue_empty(&sources)) {
        igraph_real_t tmp = igraph_dqueue_pop(&sources); node = (long) tmp;
        /* Add the node to the result vector */
        igraph_vector_push_back(res, node);
        /* Exclude the node from further source searches */
        VECTOR(degrees)[node] = -1;
        /* Get the neighbors and decrease their degrees by one */
        IGRAPH_CHECK(igraph_neighbors(graph, &neis, (igraph_integer_t) node, mode));
        j = igraph_vector_size(&neis);
        for (i = 0; i < j; i++) {
            VECTOR(degrees)[(long)VECTOR(neis)[i]]--;
            if (VECTOR(degrees)[(long)VECTOR(neis)[i]] == 0) {
                IGRAPH_CHECK(igraph_dqueue_push(&sources, VECTOR(neis)[i]));
            }
        }
    }

    if (igraph_vector_size(res) < no_of_nodes) {
        IGRAPH_ERROR("The graph has cycles; topological sorting is only possible in acyclic graphs", IGRAPH_EINVAL);
    }

    igraph_vector_destroy(&degrees);
    igraph_vector_destroy(&neis);
    igraph_dqueue_destroy(&sources);
    IGRAPH_FINALLY_CLEAN(3);

    return 0;
}

/**
 * \function igraph_is_dag
 * Checks whether a graph is a directed acyclic graph (DAG) or not.
 *
 * </para><para>
 * A directed acyclic graph (DAG) is a directed graph with no cycles.
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
int igraph_is_dag(const igraph_t* graph, igraph_bool_t *res) {
    long int no_of_nodes = igraph_vcount(graph);
    igraph_vector_t degrees, neis;
    igraph_dqueue_t sources;
    long int node, i, j, nei, vertices_left;

    if (!igraph_is_directed(graph)) {
        *res = 0;
        return IGRAPH_SUCCESS;
    }

    IGRAPH_VECTOR_INIT_FINALLY(&degrees, no_of_nodes);
    IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
    IGRAPH_CHECK(igraph_dqueue_init(&sources, 0));
    IGRAPH_FINALLY(igraph_dqueue_destroy, &sources);
    IGRAPH_CHECK(igraph_degree(graph, &degrees, igraph_vss_all(), IGRAPH_OUT, 1));

    vertices_left = no_of_nodes;

    /* Do we have nodes with no incoming edges? */
    for (i = 0; i < no_of_nodes; i++) {
        if (VECTOR(degrees)[i] == 0) {
            IGRAPH_CHECK(igraph_dqueue_push(&sources, i));
        }
    }

    /* Take all nodes with no incoming edges and remove them */
    while (!igraph_dqueue_empty(&sources)) {
        igraph_real_t tmp = igraph_dqueue_pop(&sources); node = (long) tmp;
        /* Exclude the node from further source searches */
        VECTOR(degrees)[node] = -1;
        vertices_left--;
        /* Get the neighbors and decrease their degrees by one */
        IGRAPH_CHECK(igraph_neighbors(graph, &neis, (igraph_integer_t) node,
                                      IGRAPH_IN));
        j = igraph_vector_size(&neis);
        for (i = 0; i < j; i++) {
            nei = (long)VECTOR(neis)[i];
            if (nei == node) {
                continue;
            }
            VECTOR(degrees)[nei]--;
            if (VECTOR(degrees)[nei] == 0) {
                IGRAPH_CHECK(igraph_dqueue_push(&sources, nei));
            }
        }
    }

    *res = (vertices_left == 0);
    if (vertices_left < 0) {
        IGRAPH_WARNING("vertices_left < 0 in igraph_is_dag, possible bug");
    }

    igraph_vector_destroy(&degrees);
    igraph_vector_destroy(&neis);
    igraph_dqueue_destroy(&sources);
    IGRAPH_FINALLY_CLEAN(3);

    return IGRAPH_SUCCESS;
}

/* Create the transitive closure of a tree graph.
   This is fairly simple, we just collect all ancestors of a vertex
   using a depth-first search.
 */
int igraph_transitive_closure_dag(const igraph_t *graph,
                                  igraph_t *closure) {

    long int no_of_nodes = igraph_vcount(graph);
    igraph_vector_t deg;
    igraph_vector_t new_edges;
    igraph_vector_t ancestors;
    long int root;
    igraph_vector_t neighbors;
    igraph_stack_t path;
    igraph_vector_bool_t done;

    if (!igraph_is_directed(graph)) {
        IGRAPH_ERROR("Tree transitive closure of a directed graph",
                     IGRAPH_EINVAL);
    }

    IGRAPH_VECTOR_INIT_FINALLY(&new_edges, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&deg, no_of_nodes);
    IGRAPH_VECTOR_INIT_FINALLY(&ancestors, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&neighbors, 0);
    IGRAPH_CHECK(igraph_stack_init(&path, 0));
    IGRAPH_FINALLY(igraph_stack_destroy, &path);
    IGRAPH_CHECK(igraph_vector_bool_init(&done, no_of_nodes));
    IGRAPH_FINALLY(igraph_vector_bool_destroy, &done);

    IGRAPH_CHECK(igraph_degree(graph, &deg, igraph_vss_all(),
                               IGRAPH_OUT, IGRAPH_LOOPS));

#define STAR (-1)

    for (root = 0; root < no_of_nodes; root++) {
        if (VECTOR(deg)[root] != 0) {
            continue;
        }
        IGRAPH_CHECK(igraph_stack_push(&path, root));

        while (!igraph_stack_empty(&path)) {
            long int node = (long int) igraph_stack_top(&path);
            if (node == STAR) {
                /* Leaving a node */
                long int j, n;
                igraph_stack_pop(&path);
                node = (long int) igraph_stack_pop(&path);
                if (!VECTOR(done)[node]) {
                    igraph_vector_pop_back(&ancestors);
                    VECTOR(done)[node] = 1;
                }
                n = igraph_vector_size(&ancestors);
                for (j = 0; j < n; j++) {
                    IGRAPH_CHECK(igraph_vector_push_back(&new_edges, node));
                    IGRAPH_CHECK(igraph_vector_push_back(&new_edges,
                                                         VECTOR(ancestors)[j]));
                }
            } else {
                /* Getting into a node */
                long int n, j;
                if (!VECTOR(done)[node]) {
                    IGRAPH_CHECK(igraph_vector_push_back(&ancestors, node));
                }
                IGRAPH_CHECK(igraph_neighbors(graph, &neighbors,
                                              (igraph_integer_t) node, IGRAPH_IN));
                n = igraph_vector_size(&neighbors);
                IGRAPH_CHECK(igraph_stack_push(&path, STAR));
                for (j = 0; j < n; j++) {
                    long int nei = (long int) VECTOR(neighbors)[j];
                    IGRAPH_CHECK(igraph_stack_push(&path, nei));
                }
            }
        }
    }

#undef STAR

    igraph_vector_bool_destroy(&done);
    igraph_stack_destroy(&path);
    igraph_vector_destroy(&neighbors);
    igraph_vector_destroy(&ancestors);
    igraph_vector_destroy(&deg);
    IGRAPH_FINALLY_CLEAN(5);

    IGRAPH_CHECK(igraph_create(closure, &new_edges, (igraph_integer_t)no_of_nodes,
                               IGRAPH_DIRECTED));

    igraph_vector_destroy(&new_edges);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}
