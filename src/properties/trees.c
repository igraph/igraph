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
#include "igraph_topology.h"

#include "igraph_constructors.h"
#include "igraph_dqueue.h"
#include "igraph_interface.h"
#include "igraph_stack.h"

/**
 * \function igraph_unfold_tree
 * \brief Unfolding a graph into a tree, by possibly multiplicating its vertices.
 *
 * A graph is converted into a tree (or forest, if it is unconnected),
 * by performing a breadth-first search on it, and replicating
 * vertices that were found a second, third, etc. time.
 *
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
igraph_error_t igraph_unfold_tree(const igraph_t *graph, igraph_t *tree,
                       igraph_neimode_t mode, const igraph_vector_int_t *roots,
                       igraph_vector_int_t *vertex_index) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_integer_t no_of_roots = igraph_vector_int_size(roots);
    igraph_integer_t tree_vertex_count = no_of_nodes;

    igraph_vector_int_t edges;
    igraph_vector_bool_t seen_vertices;
    igraph_vector_bool_t seen_edges;

    igraph_dqueue_int_t Q;
    igraph_vector_int_t neis;

    igraph_integer_t v_ptr = no_of_nodes;

    if (! igraph_vector_int_isininterval(roots, 0, no_of_nodes-1)) {
        IGRAPH_ERROR("All roots should be vertices of the graph.", IGRAPH_EINVVID);
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);
    IGRAPH_CHECK(igraph_vector_int_reserve(&edges, no_of_edges * 2));
    IGRAPH_DQUEUE_INT_INIT_FINALLY(&Q, 100);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&neis, 0);
    IGRAPH_VECTOR_BOOL_INIT_FINALLY(&seen_vertices, no_of_nodes);
    IGRAPH_VECTOR_BOOL_INIT_FINALLY(&seen_edges, no_of_edges);

    if (vertex_index) {
        IGRAPH_CHECK(igraph_vector_int_range(vertex_index, 0, no_of_nodes));
    }

    for (igraph_integer_t r = 0; r < no_of_roots; r++) {

        igraph_integer_t root = VECTOR(*roots)[r];
        VECTOR(seen_vertices)[root] = true;
        IGRAPH_CHECK(igraph_dqueue_int_push(&Q, root));

        while (!igraph_dqueue_int_empty(&Q)) {
            igraph_integer_t actnode = igraph_dqueue_int_pop(&Q);

            IGRAPH_CHECK(igraph_incident(graph, &neis, actnode, mode));

            igraph_integer_t n = igraph_vector_int_size(&neis);
            for (igraph_integer_t i = 0; i < n; i++) {

                igraph_integer_t edge = VECTOR(neis)[i];
                igraph_integer_t from = IGRAPH_FROM(graph, edge);
                igraph_integer_t to = IGRAPH_TO(graph, edge);
                igraph_integer_t nei = IGRAPH_OTHER(graph, edge, actnode);

                if (! VECTOR(seen_edges)[edge]) {

                    VECTOR(seen_edges)[edge] = true;

                    if (! VECTOR(seen_vertices)[nei]) {

                        IGRAPH_CHECK(igraph_vector_int_push_back(&edges, from));
                        IGRAPH_CHECK(igraph_vector_int_push_back(&edges, to));

                        VECTOR(seen_vertices)[nei] = true;
                        IGRAPH_CHECK(igraph_dqueue_int_push(&Q, nei));

                    } else {

                        tree_vertex_count++;
                        if (vertex_index) {
                            IGRAPH_CHECK(igraph_vector_int_push_back(vertex_index, nei));
                        }

                        if (from == nei) {
                            IGRAPH_CHECK(igraph_vector_int_push_back(&edges, v_ptr++));
                            IGRAPH_CHECK(igraph_vector_int_push_back(&edges, to));
                        } else {
                            IGRAPH_CHECK(igraph_vector_int_push_back(&edges, from));
                            IGRAPH_CHECK(igraph_vector_int_push_back(&edges, v_ptr++));
                        }
                    }
                }

            } /* for i<n */

        } /* ! igraph_dqueue_int_empty(&Q) */

    } /* r < igraph_vector_int_size(roots) */

    igraph_vector_bool_destroy(&seen_edges);
    igraph_vector_bool_destroy(&seen_vertices);
    igraph_vector_int_destroy(&neis);
    igraph_dqueue_int_destroy(&Q);
    IGRAPH_FINALLY_CLEAN(4);

    IGRAPH_CHECK(igraph_create(tree, &edges, tree_vertex_count, igraph_is_directed(graph)));
    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

/* igraph_is_tree -- check if a graph is a tree */

/* count the number of vertices reachable from the root */
static igraph_error_t igraph_i_is_tree_visitor(const igraph_t *graph, igraph_integer_t root, igraph_neimode_t mode, igraph_integer_t *visited_count) {
    igraph_stack_int_t stack;
    igraph_vector_bool_t visited;
    igraph_vector_int_t neighbors;
    igraph_integer_t i;

    IGRAPH_VECTOR_INT_INIT_FINALLY(&neighbors, 0);

    IGRAPH_VECTOR_BOOL_INIT_FINALLY(&visited, igraph_vcount(graph));

    IGRAPH_CHECK(igraph_stack_int_init(&stack, 0));
    IGRAPH_FINALLY(igraph_stack_int_destroy, &stack);

    *visited_count = 0;

    /* push the root into the stack */
    IGRAPH_CHECK(igraph_stack_int_push(&stack, root));

    while (! igraph_stack_int_empty(&stack)) {
        igraph_integer_t u;
        igraph_integer_t ncount;

        /* take a vertex from the stack, mark it as visited */
        u = igraph_stack_int_pop(&stack);
        if (IGRAPH_LIKELY(! VECTOR(visited)[u])) {
            VECTOR(visited)[u] = true;
            *visited_count += 1;
        }

        /* register all its yet-unvisited neighbours for future processing */
        IGRAPH_CHECK(igraph_neighbors(graph, &neighbors, u, mode));
        ncount = igraph_vector_int_size(&neighbors);
        for (i = 0; i < ncount; ++i) {
            igraph_integer_t v = VECTOR(neighbors)[i];
            if (! VECTOR(visited)[v]) {
                IGRAPH_CHECK(igraph_stack_int_push(&stack, v));
            }
        }
    }

    igraph_vector_int_destroy(&neighbors);
    igraph_stack_int_destroy(&stack);
    igraph_vector_bool_destroy(&visited);
    IGRAPH_FINALLY_CLEAN(3);

    return IGRAPH_SUCCESS;
}


/**
 * \ingroup structural
 * \function igraph_is_tree
 * \brief Decides whether the graph is a tree.
 *
 * An undirected graph is a tree if it is connected and has no cycles.
 *
 * </para><para>
 * In the directed case, an additional requirement is that all edges
 * are oriented away from a root (out-tree or arborescence) or all edges
 * are oriented towards a root (in-tree or anti-arborescence).
 * This test can be controlled using the \p mode parameter.
 *
 * </para><para>
 * By convention, the null graph (i.e. the graph with no vertices) is considered
 * not to be connected, and therefore not a tree.
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
 * \sa \ref igraph_is_connected()
 *
 * \example examples/simple/igraph_kary_tree.c
 */
igraph_error_t igraph_is_tree(const igraph_t *graph, igraph_bool_t *res, igraph_integer_t *root, igraph_neimode_t mode) {
    igraph_bool_t is_tree = false;
    igraph_bool_t treat_as_undirected = !igraph_is_directed(graph) || mode == IGRAPH_ALL;
    igraph_integer_t iroot = 0;
    igraph_integer_t visited_count;
    igraph_integer_t vcount, ecount;

    vcount = igraph_vcount(graph);
    ecount = igraph_ecount(graph);

    if (igraph_i_property_cache_has(graph, IGRAPH_PROP_IS_WEAKLY_CONNECTED)) {
        igraph_bool_t weakly_connected = igraph_i_property_cache_get_bool(graph, IGRAPH_PROP_IS_WEAKLY_CONNECTED);

        if (weakly_connected) {
            /* For undirected graphs and for directed graphs with mode == IGRAPH_ALL,
             * we can return early if we know from the cache that the graph is weakly
             * connected and is a forest. We can do this even if the user wants the
             * root vertex because we always return zero as the root vertex for
             * undirected graphs */
            if (treat_as_undirected &&
                igraph_i_property_cache_has(graph, IGRAPH_PROP_IS_FOREST) &&
                igraph_i_property_cache_get_bool(graph, IGRAPH_PROP_IS_FOREST)
            ) {
                is_tree = true;
                iroot = 0;
                goto success;
            }
        } else /* ! weakly_connected */ {
            /* If the graph is not weakly connected, then it is neither an undirected
             * not a directed tree. There is no root, so we can return early. */
            is_tree = false;
            goto success;
        }
    }

    /* A tree must have precisely vcount-1 edges. */
    /* By convention, the zero-vertex graph will not be considered a tree. */
    if (ecount != vcount - 1) {
        is_tree = false;
        goto success;
    }

    /* The single-vertex graph is a tree, provided it has no edges (checked in the previous if (..)) */
    if (vcount == 1) {
        is_tree = true;
        iroot = 0;
        goto success;
    }

    /* For higher vertex counts we cannot short-circuit due to the possibility
     * of loops or multi-edges even when the edge count is correct. */

    /* Ignore mode for undirected graphs. */
    if (! igraph_is_directed(graph)) {
        mode = IGRAPH_ALL;
    }

    /* The main algorithm:
     * We find a root and check that all other vertices are reachable from it.
     * We have already checked the number of edges, so with the additional
     * reachability condition we can verify if the graph is a tree.
     *
     * For directed graphs, the root is the node with no incoming/outgoing
     * connections, depending on 'mode'. For undirected, it is arbitrary, so
     * we choose 0.
     */

    is_tree = true; /* assume success */

    switch (mode) {
    case IGRAPH_ALL:
        iroot = 0;
        break;

    case IGRAPH_IN:
    case IGRAPH_OUT: {
        igraph_vector_int_t degree;
        igraph_integer_t i;

        IGRAPH_CHECK(igraph_vector_int_init(&degree, 0));
        IGRAPH_FINALLY(igraph_vector_int_destroy, &degree);

        IGRAPH_CHECK(igraph_degree(graph, &degree, igraph_vss_all(), mode == IGRAPH_IN ? IGRAPH_OUT : IGRAPH_IN, IGRAPH_LOOPS));

        for (i = 0; i < vcount; ++i) {
            if (VECTOR(degree)[i] == 0) {
                break;
            }
            if (VECTOR(degree)[i] > 1) {
                /* In an out-tree, all vertices have in-degree 1, except for the root,
                 * which has in-degree 0. Thus, if we encounter a larger in-degree,
                 * the graph cannot be an out-tree.
                 * We could perform this check for all degrees, but that would not
                 * improve performance when the graph is indeed a tree, persumably
                 * the most common case. Thus we only check until finding the root.
                 */
                is_tree = false;
                break;
            }
        }

        /* If no suitable root is found, the graph is not a tree. */
        if (is_tree && i == vcount) {
            is_tree = false;
        } else {
            iroot = i;
        }

        igraph_vector_int_destroy(&degree);
        IGRAPH_FINALLY_CLEAN(1);
    }

    break;
    default:
        IGRAPH_ERROR("Invalid mode.", IGRAPH_EINVMODE);
    }

    /* if no suitable root was found, skip visiting vertices */
    if (is_tree) {
        IGRAPH_CHECK(igraph_i_is_tree_visitor(graph, iroot, mode, &visited_count));
        is_tree = visited_count == vcount;
    }

success:
    if (res) {
        *res = is_tree;
    }

    if (root) {
        *root = iroot;
    }

    if (is_tree) {
        /* A graph that is a directed tree is also an undirected tree.
         * An undirected tree is weakly connected and is a forest,
         * so we can cache this. */
        igraph_i_property_cache_set_bool_checked(graph, IGRAPH_PROP_IS_FOREST, true);
        igraph_i_property_cache_set_bool_checked(graph, IGRAPH_PROP_IS_WEAKLY_CONNECTED, true);
    }

    return IGRAPH_SUCCESS;
}

/* igraph_is_forest() -- check if a graph is a forest */

/* Verify that the graph has no cycles and count the number of reachable vertices.
 * This function performs a DFS starting from 'root'.
 * If it finds a cycle, it sets *res to false, otherwise it does not change it.
 * *visited_count will be incremented by the number of vertices reachable from 'root',
 * including 'root' itself.
 */
static igraph_error_t igraph_i_is_forest_visitor(
        const igraph_t *graph, igraph_integer_t root, igraph_neimode_t mode,
        igraph_vector_bool_t *visited, igraph_stack_int_t *stack, igraph_vector_int_t *neis,
        igraph_integer_t *visited_count, igraph_bool_t *res)
{
    igraph_integer_t i;

    igraph_stack_int_clear(stack);

    /* push the root onto the stack */
    IGRAPH_CHECK(igraph_stack_int_push(stack, root));

    while (! igraph_stack_int_empty(stack)) {
        igraph_integer_t u;
        igraph_integer_t ncount;

        /* Take a vertex from stack and check if it is already visited.
         * If yes, then we found a cycle: the graph is not a forest.
         * Otherwise mark it as visited and continue.
         */
        u = igraph_stack_int_pop(stack);
        if (IGRAPH_LIKELY(! VECTOR(*visited)[u])) {
            VECTOR(*visited)[u] = true;
            *visited_count += 1;
        }
        else {
            *res = false;
            break;
        }

        /* Vertex discovery: Register all its neighbours for future processing */
        IGRAPH_CHECK(igraph_neighbors(graph, neis, u, mode));
        ncount = igraph_vector_int_size(neis);

        for (i = 0; i < ncount; ++i) {
            igraph_integer_t v = VECTOR(*neis)[i];

            if (mode == IGRAPH_ALL) {
                /* In the undirected case, we avoid returning to the predecessor
                 * vertex of 'v' in the DFS tree by skipping visited vertices.
                 *
                 * Note that in order to succcessfully detect a cycle, a vertex
                 * within that cycle must end up on the stack more than once.
                 * Does skipping visited vertices preclude this sometimes?
                 * No, because any visited vertex can only be accessed through
                 * an already discovered vertex (i.e. one that has already been
                 * pushed onto the stack).
                 */
                if (IGRAPH_LIKELY(! VECTOR(*visited)[v])) {
                    IGRAPH_CHECK(igraph_stack_int_push(stack, v));
                }
                /* To check for a self-loop in undirected graph */
                else if (v == u) {
                    *res = false;
                    break;
                }
            }
            else {
                IGRAPH_CHECK(igraph_stack_int_push(stack, v));
            }
        }
    }

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_is_forest(
    const igraph_t *graph, igraph_bool_t *res,
    igraph_vector_int_t *roots, igraph_neimode_t mode
);

/**
 * \ingroup structural
 * \function igraph_is_forest
 * \brief Decides whether the graph is a forest.
 *
 * An undirected graph is a forest if it has no cycles. Equivalently,
 * a graph is a forest if all connected components are trees.
 *
 * </para><para>
 * In the directed case, an additional requirement is that edges in each
 * tree are oriented away from the root (out-trees or arborescences) or all edges
 * are oriented towards the root (in-trees or anti-arborescences).
 * This test can be controlled using the \p mode parameter.
 *
 * </para><para>
 * By convention, the null graph (i.e. the graph with no vertices) is considered to be a forest.
 *
 * </para><para>
 * The \p res return value of this function is cached in the graph itself if
 * \p mode is set to \c IGRAPH_ALL or if the graph is undirected. Calling the
 * function multiple times with no modifications to the graph in between
 * will return a cached value in O(1) time if the roots are not requested.
 *
 * \param graph The graph object to analyze.
 * \param res Pointer to a logical variable. If not \c NULL, then the result will be stored
 *        here.
 * \param roots If not \c NULL, the root nodes will be stored here. When \p mode
 *        is \c IGRAPH_ALL or the graph is undirected, any one vertex from each
 *        component can be the root. When \p mode is \c IGRAPH_OUT
 *        or \c IGRAPH_IN, all the vertices with zero in- or out-degree,
 *        respectively are considered as root nodes.
 * \param mode For a directed graph this specifies whether to test for an
 *        out-forest, an in-forest or ignore edge directions. The respective
 *        possible values are:
 *        \c IGRAPH_OUT, \c IGRAPH_IN, \c IGRAPH_ALL. This argument is
 *        ignored for undirected graphs.
 * \return Error code:
 *        \c IGRAPH_EINVMODE: invalid mode argument.
 *
 * Time complexity: At most O(|V|+|E|), the
 * number of vertices plus the number of edges in the graph.
 */
igraph_error_t igraph_is_forest(const igraph_t *graph, igraph_bool_t *res,
                                igraph_vector_int_t *roots, igraph_neimode_t mode) {

    const igraph_bool_t treat_as_undirected = !igraph_is_directed(graph) || mode == IGRAPH_ALL;

    if (!roots && !res) {
        return IGRAPH_SUCCESS;
    }

    /* Note on cache use:
     *
     * The IGRAPH_PROP_IS_FOREST cached property is equivalent to this function's
     * result ONLY in the undirected case. Keep in mind that a graph that is not
     * a directed forest may still be an undirected forest, i.e. may still be free
     * of undirected cycles. Example: 1->2<-3->4.
     */

    if (igraph_i_property_cache_has(graph, IGRAPH_PROP_IS_FOREST)) {
        const igraph_bool_t no_undirected_cycles = igraph_i_property_cache_get_bool(graph, IGRAPH_PROP_IS_FOREST);
        if (treat_as_undirected && res && ! roots) {
            /* If the graph is treated as undirected and no roots are requested,
             * we can directly use the cached IGRAPH_PROP_IS_FOREST value. */
            *res = no_undirected_cycles;
            return IGRAPH_SUCCESS;
        } else {
            /* Otherwise we can use negative cached values (i.e. "false"):
             *  - A graph with undirected cycles cannot be a directed forest.
             *  - If the graph is not a forest, we don't need to look for roots.
             */
            if (! no_undirected_cycles) {
                if (res) { res = false; }
                if (roots) { igraph_vector_int_clear(roots); }
                return IGRAPH_SUCCESS;
            }
        }
    }

    IGRAPH_CHECK(igraph_i_is_forest(graph, res, roots, mode));

    /* At this point we know whether the graph is an (undirected or directed) forest
     * as we have at least one of 'res' or 'roots'. The case when both are NULL was
     * caught above. */
    igraph_bool_t is_forest;
    if (res != NULL) {
        is_forest = *res;
    } else /* roots != NULL */ {
        is_forest = igraph_vcount(graph) == 0 || !igraph_vector_int_empty(roots);
    }
    if (is_forest) {
        /* If the graph is a directed forest, then it has no undirected cycles.
         * We can enter positive results in the cache unconditionally. */
        igraph_i_property_cache_set_bool_checked(graph, IGRAPH_PROP_IS_FOREST, true);
    } else if (treat_as_undirected) {
        /* However, if the graph is not a directed forest, it might still be
         * an undirected forest. We can only enter negative results in the cache
         * when edge directions were ignored, but NOT in the directed case. */
        igraph_i_property_cache_set_bool_checked(graph, IGRAPH_PROP_IS_FOREST, false);
    }

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_is_forest(
    const igraph_t *graph, igraph_bool_t *res,
    igraph_vector_int_t *roots, igraph_neimode_t mode
) {
    igraph_vector_bool_t visited;
    igraph_vector_int_t neis;
    igraph_stack_int_t stack;
    igraph_integer_t visited_count = 0;
    igraph_integer_t vcount, ecount;
    igraph_integer_t v;
    igraph_bool_t result;

    vcount = igraph_vcount(graph);
    ecount = igraph_ecount(graph);

    if (roots) {
        igraph_vector_int_clear(roots);
    }

    /* Any graph with 0 edges is a forest. */
    if (ecount == 0) {
        if (res) {
            *res = true;
        }
        if (roots) {
            for (v = 0; v < vcount; v++) {
                IGRAPH_CHECK(igraph_vector_int_push_back(roots, v));
            }
        }
        return IGRAPH_SUCCESS;
    }

    /* A forest can have at most vcount-1 edges. */
    if (ecount > vcount - 1) {
        if (res) {
            *res = false;
        }
        return IGRAPH_SUCCESS;
    }

    /* Ignore mode for undirected graphs. */
    if (! igraph_is_directed(graph)) {
        mode = IGRAPH_ALL;
    }

    result = true; /* assume success */

    IGRAPH_VECTOR_BOOL_INIT_FINALLY(&visited, vcount);

    IGRAPH_CHECK(igraph_stack_int_init(&stack, 0));
    IGRAPH_FINALLY(igraph_stack_int_destroy, &stack);

    IGRAPH_VECTOR_INT_INIT_FINALLY(&neis, 0);

    /* The main algorithm:
     *
     * Undirected Graph:- We add each unvisited vertex to the roots vector, and
     * mark all other vertices that are reachable from it as visited.
     *
     * Directed Graph:- For each tree, the root is the node with no
     * incoming/outgoing connections, depending on 'mode'. We add each vertex
     * with zero degree to the roots vector and mark all other vertices that are
     * reachable from it as visited.
     *
     * If all the vertices are visited exactly once, then the graph is a forest.
     */

    switch (mode) {
        case IGRAPH_ALL:
        {
            for (v = 0; v < vcount; ++v) {
                if (!result) {
                    break;
                }
                if (! VECTOR(visited)[v]) {
                    if (roots) {
                        IGRAPH_CHECK(igraph_vector_int_push_back(roots, v));
                    }
                    IGRAPH_CHECK(igraph_i_is_forest_visitor(
                                     graph, v, mode,
                                     &visited, &stack, &neis,
                                     &visited_count, &result));
                }
            }
            break;
        }

        case IGRAPH_IN:
        case IGRAPH_OUT:
        {
            igraph_vector_int_t degree;

            IGRAPH_VECTOR_INT_INIT_FINALLY(&degree, 0);
            IGRAPH_CHECK(igraph_degree(graph, &degree, igraph_vss_all(),
                            IGRAPH_REVERSE_MODE(mode), /* loops = */ 1));

            for (v = 0; v < vcount; ++v) {
                /* In an out-tree, roots have in-degree 0,
                 * and all other vertices have in-degree 1. */
                if (VECTOR(degree)[v] > 1 || !result) {
                    result = false;
                    break;
                }
                if (VECTOR(degree)[v] == 0) {
                    if (roots) {
                        IGRAPH_CHECK(igraph_vector_int_push_back(roots, v));
                    }
                    IGRAPH_CHECK(igraph_i_is_forest_visitor(
                                     graph, v, mode,
                                     &visited, &stack, &neis,
                                     &visited_count, &result));
                }
            }

            igraph_vector_int_destroy(&degree);
            IGRAPH_FINALLY_CLEAN(1);
            break;
        }

        default:
            IGRAPH_ERROR("Invalid mode.", IGRAPH_EINVMODE);
    }

    if (result) {
        /* In a forest, all vertices are reachable from the roots. */
        result = (visited_count == vcount);
    }

    if (res) {
        *res = result;
    }

    /* If the graph is not a forest then the root vector will be empty. */
    if (!result && roots) {
        igraph_vector_int_clear(roots);
    }

    igraph_vector_int_destroy(&neis);
    igraph_stack_int_destroy(&stack);
    igraph_vector_bool_destroy(&visited);
    IGRAPH_FINALLY_CLEAN(3);

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup structural
 * \function igraph_is_acyclic
 * \brief Checks whether a graph is acyclic or not.
 *
 * This function checks whether a graph is acyclic or not.
 *
 * \param graph The input graph.
 * \param res Pointer to a boolean constant, the result
        is stored here.
 * \return Error code.
 *
 * Time complexity: O(|V|+|E|), where |V| and |E| are the number of
 * vertices and edges in the original input graph.
 */
igraph_error_t igraph_is_acyclic(const igraph_t *graph, igraph_bool_t *res) {
    if (igraph_is_directed(graph)) {
        /* igraph_is_dag is cached */
        return igraph_is_dag(graph, res);
    } else {
        /* igraph_is_forest is cached if mode == IGRAPH_ALL and we don't need
         * the roots */
        return igraph_is_forest(graph, res, NULL, IGRAPH_ALL);
    }
}
