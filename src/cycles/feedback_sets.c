/*
   igraph library.
   Copyright (C) 2011-2024  The igraph development team <igraph@igraph.org>

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include "igraph_cycles.h"
#include "cycles/feedback_sets.h"

#include "igraph_bitset.h"
#include "igraph_components.h"
#include "igraph_dqueue.h"
#include "igraph_interface.h"
#include "igraph_memory.h"
#include "igraph_stack.h"
#include "igraph_structural.h"
#include "igraph_vector.h"
#include "igraph_vector_list.h"
#include "igraph_visitor.h"

#include "internal/glpk_support.h"
#include "math/safe_intop.h"

#include <limits.h>

/**
 * \param found If not NULL, will indicate if any cycles were found.
 * \param  A bitset the same length as the edge count, indicating which edges
 *    to consider non-existent.
 *
 * See igraph_find_cycle() for the other parameters.
 *
 * \return Error code.
 */
static igraph_error_t igraph_i_find_cycle(const igraph_t *graph,
                                          igraph_vector_int_t *vertices,
                                          igraph_vector_int_t *edges,
                                          igraph_bool_t *found,
                                          igraph_neimode_t mode,
                                          const igraph_bitset_t *removed) {

    const igraph_int_t vcount = igraph_vcount(graph);
    const igraph_int_t ecount = igraph_ecount(graph);
    igraph_stack_int_t stack;
    igraph_vector_int_t inc;
    igraph_vector_int_t vpath, epath;
    igraph_vector_char_t seen; /* 0 = unseen, 1 = acestor of current, 2 = seen, non-ancestor */
    igraph_int_t ea, va;
    igraph_int_t depth;

    if (vertices) {
        igraph_vector_int_clear(vertices);
    }
    if (edges) {
        igraph_vector_int_clear(edges);
    }
    if (found) {
        *found = false;
    }

    if (ecount == 0) {
        return IGRAPH_SUCCESS;
    }

#define PATH_PUSH(v, e) \
    do { \
        IGRAPH_CHECK(igraph_vector_int_push_back(&epath, e)); \
        IGRAPH_CHECK(igraph_vector_int_push_back(&vpath, v)); \
        VECTOR(seen)[v] = 1; \
    } while (0)

#define PATH_POP() \
    do { \
        igraph_vector_int_pop_back(&epath); \
        VECTOR(seen)[igraph_vector_int_pop_back(&vpath)] = 2; \
    } while (0)

    IGRAPH_VECTOR_INT_INIT_FINALLY(&vpath, 100);
    igraph_vector_int_clear(&vpath);

    IGRAPH_VECTOR_INT_INIT_FINALLY(&epath, 100);
    igraph_vector_int_clear(&epath);

    IGRAPH_VECTOR_INT_INIT_FINALLY(&inc, 10);
    IGRAPH_VECTOR_CHAR_INIT_FINALLY(&seen, vcount);
    IGRAPH_STACK_INT_INIT_FINALLY(&stack, 200);

    for (igraph_int_t v=0; v < vcount; v++) {
        if (VECTOR(seen)[v]) {
            continue;
        }

        IGRAPH_CHECK(igraph_stack_int_push(&stack, -1));
        IGRAPH_CHECK(igraph_stack_int_push(&stack, v));

        while (! igraph_stack_int_empty(&stack)) {
            igraph_int_t x = igraph_stack_int_pop(&stack);
            if (x == -1) {
                PATH_POP();
                continue;
            } else {
                va = x;
                ea = igraph_stack_int_pop(&stack);

                if (VECTOR(seen)[va] == 1) {
                    goto finish;
                } else if (VECTOR(seen)[va] == 2) {
                    continue;
                }
            }

            PATH_PUSH(va, ea);

            IGRAPH_CHECK(igraph_stack_int_push(&stack, -1));
            IGRAPH_CHECK(igraph_incident(graph, &inc, va, mode, IGRAPH_LOOPS));
            igraph_int_t n = igraph_vector_int_size(&inc);
            for (igraph_int_t i=0; i < n; i++) {
                igraph_int_t eb = VECTOR(inc)[i];
                igraph_int_t vb = IGRAPH_OTHER(graph, eb, va);
                if (eb == ea) continue;
                if (VECTOR(seen)[vb] == 2) continue;
                if (removed && IGRAPH_BIT_TEST(*removed, eb)) continue;
                IGRAPH_CHECK(igraph_stack_int_push(&stack, eb));
                IGRAPH_CHECK(igraph_stack_int_push(&stack, vb));
            }
        }
    }


finish:

    igraph_stack_int_destroy(&stack);
    igraph_vector_char_destroy(&seen);
    igraph_vector_int_destroy(&inc);
    IGRAPH_FINALLY_CLEAN(3);

    depth = igraph_vector_int_size(&vpath);

    if (depth > 0) {
        igraph_int_t i = depth;
        while (VECTOR(vpath)[i-1] != va) i--;
        for (; i < depth; i++) {
            if (vertices) {
                IGRAPH_CHECK(igraph_vector_int_push_back(vertices, VECTOR(vpath)[i]));
            }
            if (edges) {
                IGRAPH_CHECK(igraph_vector_int_push_back(edges, VECTOR(epath)[i]));
            }
        }
        if (vertices) {
            IGRAPH_CHECK(igraph_vector_int_push_back(vertices, va));
        }
        if (edges) {
            IGRAPH_CHECK(igraph_vector_int_push_back(edges, ea));
        }
        if (found) {
            *found = true;
        }
    }

    igraph_vector_int_destroy(&epath);
    igraph_vector_int_destroy(&vpath);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;

#undef PATH_PUSH
#undef PATH_POP
}


/**
 * \function igraph_find_cycle
 * \brief Finds a single cycle in the graph.
 *
 * This function returns a cycle of the graph, if there is one. If the graph
 * is acyclic, it returns empty vectors.
 *
 * \param graph The input graph.
 * \param vertices Pointer to an integer vector. If a cycle is found, its
 *    vertices will be stored here. Otherwise the vector will be empty.
 * \param edges Pointer to an integer vector. If a cycle is found, its
 *    edges will be stored here. Otherwise the vector will be empty.
 * \param mode A constant specifying how edge directions are
 *        considered in directed graphs. Valid modes are:
 *        \c IGRAPH_OUT, follows edge directions;
 *        \c IGRAPH_IN, follows the opposite directions; and
 *        \c IGRAPH_ALL, ignores edge directions. This argument is
 *        ignored for undirected graphs.
 * \return Error code.
 *
 * Time complexity: O(|V|+|E|), where |V| and |E| are the number of
 * vertices and edges in the original input graph.
 *
 * \sa \ref igraph_is_acyclic() to determine if a graph is acyclic,
 * without returning a specific cycle; \ref igraph_simple_cycles()
 * to list all cycles in a graph.
 */
igraph_error_t igraph_find_cycle(const igraph_t *graph,
                                 igraph_vector_int_t *vertices,
                                 igraph_vector_int_t *edges,
                                 igraph_neimode_t mode) {

    /* If the graph is cached to be acyclic, we don't need to run the algorithm. */

    igraph_bool_t known_acyclic = false;
    igraph_bool_t found;

    if (mode != IGRAPH_OUT && mode != IGRAPH_IN && mode != IGRAPH_ALL) {
        IGRAPH_ERROR("Invalid mode for finding cycles.", IGRAPH_EINVAL);
    }

    if (! igraph_is_directed(graph)) {
        mode = IGRAPH_ALL;
    }

    if (mode == IGRAPH_ALL) /* undirected */ {
        if (igraph_i_property_cache_has(graph, IGRAPH_PROP_IS_FOREST) &&
            igraph_i_property_cache_get_bool(graph, IGRAPH_PROP_IS_FOREST)) {
            known_acyclic = true;
        }
    } else /* directed */ {
        if (igraph_i_property_cache_has(graph, IGRAPH_PROP_IS_DAG) &&
            igraph_i_property_cache_get_bool(graph, IGRAPH_PROP_IS_DAG)) {
            known_acyclic = true;
        }
    }

    if (known_acyclic) {
        if (vertices) {
            igraph_vector_int_clear(vertices);
        }
        if (edges) {
            igraph_vector_int_clear(edges);
        }
        return IGRAPH_SUCCESS;
    }

    IGRAPH_CHECK(igraph_i_find_cycle(graph, vertices, edges, &found, mode, NULL));

    if (! found) {
        if (mode == IGRAPH_ALL) /* undirected */ {
            igraph_i_property_cache_set_bool_checked(graph, IGRAPH_PROP_IS_FOREST, true);
        } else /* directed */ {
            igraph_i_property_cache_set_bool_checked(graph, IGRAPH_PROP_IS_DAG, true);
        }
    }

    return IGRAPH_SUCCESS;
}


/**
 * \ingroup structural
 * \function igraph_feedback_arc_set
 * \brief Feedback arc set of a graph using exact or heuristic methods.
 *
 * A feedback arc set is a set of edges whose removal makes the graph acyclic.
 * We are usually interested in \em minimum feedback arc sets, i.e. sets of edges
 * whose total weight is the smallest among all the feedback arc sets.
 *
 * </para><para>
 * For undirected graphs, the solution is simple: one has to find a maximum weight
 * spanning tree and then remove all the edges not in the spanning tree. For directed
 * graphs, this is an NP-complete problem, and various heuristics are usually used to
 * find an approximate solution to the problem. This function implements both exact
 * methods and heuristics, selectable with the \p algo parameter.
 *
 * </para><para>
 * References:
 *
 * </para><para>
 * Eades P, Lin X and Smyth WF:
 * A fast and effective heuristic for the feedback arc set problem.
 * Information Processing Letters 47(6), pp 319-323 (1993).
 * https://doi.org/10.1016/0020-0190(93)90079-O
 *
 * </para><para>
 * Baharev A, Hermann S, Arnold N and Tobias A:
 * An Exact Method for the Minimum Feedback Arc Set Problem.
 * ACM Journal of Experimental Algorithmics 26, 1–28 (2021).
 * https://doi.org/10.1145/3446429.
 *
 * \param graph The graph object.
 * \param result An initialized vector, the result will be written here.
 * \param weights Weight vector or \c NULL if no weights are specified.
 * \param algo The algorithm to use to solve the problem if the graph is directed.
 *        Possible values:
 *        \clist
 *        \cli IGRAPH_FAS_EXACT_IP
 *          Finds a \em minimum feedback arc set using integer programming (IP),
 *          automatically selecting the best method of this type (currently
 *          always \c IGRAPH_FAS_EXACT_IP_CG). The complexity is of course
 *          at least exponential.
 *        \cli IGRAPH_FAS_EXACT_IP_CG
 *          This is an integer programming approach based on a minimum set cover
 *          formulation and using incremental constraint generation (CG), added
 *          in igraph 0.10.14. We minimize <code>sum_e w_e b_e</code> subject to
 *          the constraints <code>sum_e c_e b_e &gt;= 1</code> for all cycles \c c.
 *          Here \c w_e is the weight of edge \c e, \c b_e is a binary variable
 *          (0 or 1) indicating whether edge \c e is in the feedback set,
 *          and \c c_e is a binary coefficient indicating whether edge \c e
 *          is in cycle \c c. The constraint expresses the requirement that all
 *          cycles must intersect with (be broken by) the edge set represented
 *          by \c b. Since there are a very large number of cycles in the graph,
 *          constraints are generated incrementally, iteratively adding some cycles
 *          that do not intersect with the current edge set \c b, then solving for
 *          \c b again, until finally no unbroken cycles remain. This approach is
 *          similar to that described by Baharev et al (though with a simpler
 *          cycle generation scheme), and to what is implemented by SageMath's.
 *          \c feedback_edge_set function.
 *        \cli IGRAPH_FAS_EXACT_IP_TI
 *          This is another integer programming approach based on finding a
 *          maximum (largest weight) edge set that adheres to a topological
 *          order. It uses the common formulation through triangle inequalities
 *          (TI), see Section 3.1 of Baharev et al (2021) for an overview. This
 *          method was used before igraph 0.10.14, and is typically much slower
 *          than \c IGRAPH_FAS_EXACT_IP_CG.
 *        \cli IGRAPH_FAS_APPROX_EADES
 *          Finds a feedback arc set using the heuristic of Eades, Lin and
 *          Smyth (1993). This is guaranteed to be smaller than |E|/2 - |V|/6,
 *          and it is linear in the number of edges (i.e. O(|E|)) to compute.
 *        \endclist
 *
 * \return Error code:
 *         \c IGRAPH_EINVAL if an unknown method was specified or the weight vector
 *            is invalid.
 *
 * \example examples/simple/igraph_feedback_arc_set.c
 * \example examples/simple/igraph_feedback_arc_set_ip.c
 *
 * Time complexity: depends on \p algo, see the time complexities there.
 */
igraph_error_t igraph_feedback_arc_set(
        const igraph_t *graph,
        igraph_vector_int_t *result,
        const igraph_vector_t *weights,
        igraph_fas_algorithm_t algo) {

    if (weights) {
        if (igraph_vector_size(weights) != igraph_ecount(graph)) {
            IGRAPH_ERROR("Weight vector length must match the number of edges.", IGRAPH_EINVAL);
        }
        if (! igraph_vector_is_all_finite(weights)) {
            IGRAPH_ERROR("Weights must not be infinite or NaN.", IGRAPH_EINVAL);
        }
    }

    if (!igraph_is_directed(graph)) {
        return igraph_i_feedback_arc_set_undirected(graph, result, weights, NULL);
    }

    switch (algo) {
    case IGRAPH_FAS_EXACT_IP:
    case IGRAPH_FAS_EXACT_IP_CG:
        return igraph_i_feedback_arc_set_ip_cg(graph, result, weights);

    case IGRAPH_FAS_EXACT_IP_TI:
        return igraph_i_feedback_arc_set_ip_ti(graph, result, weights);

    case IGRAPH_FAS_APPROX_EADES:
        return igraph_i_feedback_arc_set_eades(graph, result, weights, NULL);

    default:
        IGRAPH_ERROR("Invalid feedback arc set algorithm.", IGRAPH_EINVAL);
    }
}


/**
 * \function igraph_feedback_vertex_set
 * \brief Feedback vertex set of a graph.
 *
 * A feedback vertex set is a set of vertices whose removal makes the graph
 * acyclic. Finding a \em minimum feedback vertex set is an NP-complete
 * problem, both on directed and undirected graphs.
 *
 * \param graph The graph.
 * \param result An initialized vector, the result will be written here.
 * \param vertex_weights Vertex weight vector or \c NULL if no weights are specified.
 * \param algo Algorithm to use. Possible values:
 *        \clist
 *        \cli IGRAPH_FVS_EXACT_IP
 *         Finds a \em miniumum feedback vertex set using integer programming
 *         (IP). The complexity is of course at least exponential. Currently
 *         this method uses an approach analogous to that of the
 *         \c IGRAPH_FAS_EXACT_IP_CG algorithm of \ref  igraph_feedback_arc_set().
 *        \endclist
 *
 * \return Error code.
 *
 * Time complexity: depends on \p algo, see the time complexities there.
 */
igraph_error_t igraph_feedback_vertex_set(
    const igraph_t *graph, igraph_vector_int_t *result,
    const igraph_vector_t *vertex_weights, igraph_fvs_algorithm_t algo) {

    if (vertex_weights) {
        if (igraph_vector_size(vertex_weights) != igraph_vcount(graph)) {
            IGRAPH_ERROR("Vertex weight vector length must match the number of vertices.", IGRAPH_EINVAL);
        }
        if (! igraph_vector_is_all_finite(vertex_weights)) {
            IGRAPH_ERROR("Vertex weights must not be infinite or NaN.", IGRAPH_EINVAL);
        }
    }

    switch (algo) {
    case IGRAPH_FVS_EXACT_IP:
        return igraph_i_feedback_vertex_set_ip_cg(graph, result, vertex_weights);

    default:
        IGRAPH_ERROR("Invalid feedback vertex set algorithm.", IGRAPH_EINVAL);
    }
}


/**
 * Solves the feedback arc set problem for undirected graphs.
 */
igraph_error_t igraph_i_feedback_arc_set_undirected(const igraph_t *graph, igraph_vector_int_t *result,
        const igraph_vector_t *weights, igraph_vector_int_t *layering) {

    const igraph_int_t no_of_nodes = igraph_vcount(graph);
    const igraph_int_t no_of_edges = igraph_ecount(graph);
    igraph_vector_int_t edges;

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, no_of_nodes > 0 ? no_of_nodes - 1 : 0);
    if (weights) {
        /* Find a maximum weight spanning tree. igraph has a routine for minimum
         * spanning trees, so we negate the weights */
        igraph_vector_t vcopy;
        IGRAPH_CHECK(igraph_vector_init_copy(&vcopy, weights));
        IGRAPH_FINALLY(igraph_vector_destroy, &vcopy);
        igraph_vector_scale(&vcopy, -1);
        IGRAPH_CHECK(igraph_minimum_spanning_tree(graph, &edges, &vcopy, IGRAPH_MST_AUTOMATIC));
        igraph_vector_destroy(&vcopy);
        IGRAPH_FINALLY_CLEAN(1);
    } else {
        /* Any spanning tree will do */
        IGRAPH_CHECK(igraph_minimum_spanning_tree(graph, &edges, NULL, IGRAPH_MST_AUTOMATIC));
    }

    /* Now we have a bunch of edges that constitute a spanning forest. We have
     * to come up with a layering, and return those edges that are not in the
     * spanning forest */
    igraph_vector_int_sort(&edges);
    IGRAPH_CHECK(igraph_vector_int_push_back(&edges, -1));  /* guard element */

    if (result) {
        igraph_vector_int_clear(result);
        for (igraph_int_t i = 0, j = 0; i < no_of_edges; i++) {
            if (i == VECTOR(edges)[j]) {
                j++;
                continue;
            }
            IGRAPH_CHECK(igraph_vector_int_push_back(result, i));
        }
    }

    if (layering) {
        igraph_vector_t degrees;
        igraph_vector_int_t roots;

        IGRAPH_VECTOR_INIT_FINALLY(&degrees, no_of_nodes);
        IGRAPH_VECTOR_INT_INIT_FINALLY(&roots, no_of_nodes);
        IGRAPH_CHECK(igraph_strength(graph, &degrees, igraph_vss_all(),
                                     IGRAPH_ALL, IGRAPH_NO_LOOPS, weights));
        IGRAPH_CHECK(igraph_vector_sort_ind(&degrees, &roots, IGRAPH_DESCENDING));

        IGRAPH_CHECK(igraph_bfs(graph,
                                /* root = */ 0,
                                /* roots = */ &roots,
                                /* mode = */ IGRAPH_OUT,
                                /* unreachable = */ false,
                                /* restricted = */ NULL,
                                /* order = */ NULL,
                                /* rank = */ NULL,
                                /* parents = */ NULL,
                                /* pred = */ NULL,
                                /* succ = */ NULL,
                                /* dist = */ layering,
                                /* callback = */ NULL,
                                /* extra = */ NULL));

        igraph_vector_destroy(&degrees);
        igraph_vector_int_destroy(&roots);
        IGRAPH_FINALLY_CLEAN(2);
    }

    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}


/**
 * Solves the feedback arc set problem using the heuristics of Eades et al.
 */
igraph_error_t igraph_i_feedback_arc_set_eades(const igraph_t *graph, igraph_vector_int_t *result,
                                    const igraph_vector_t *weights, igraph_vector_int_t *layers) {
    const igraph_int_t no_of_nodes = igraph_vcount(graph);
    const igraph_int_t no_of_edges = igraph_ecount(graph);
    igraph_int_t nodes_left;
    igraph_int_t neis_size;
    igraph_dqueue_int_t sources, sinks;
    igraph_vector_int_t neis;
    igraph_vector_int_t indegrees, outdegrees;
    igraph_vector_t instrengths, outstrengths;
    igraph_vector_int_t ordering;
    igraph_int_t order_next_pos = 0, order_next_neg = -1;

    IGRAPH_VECTOR_INT_INIT_FINALLY(&ordering, no_of_nodes);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&neis, 0);

    IGRAPH_DQUEUE_INT_INIT_FINALLY(&sources, 0);
    IGRAPH_DQUEUE_INT_INIT_FINALLY(&sinks, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&indegrees, no_of_nodes);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&outdegrees, no_of_nodes);
    IGRAPH_VECTOR_INIT_FINALLY(&instrengths, no_of_nodes);
    IGRAPH_VECTOR_INIT_FINALLY(&outstrengths, no_of_nodes);

    IGRAPH_CHECK(igraph_degree(graph, &indegrees, igraph_vss_all(), IGRAPH_IN, IGRAPH_NO_LOOPS));
    IGRAPH_CHECK(igraph_degree(graph, &outdegrees, igraph_vss_all(), IGRAPH_OUT, IGRAPH_NO_LOOPS));

    if (weights) {
        IGRAPH_CHECK(igraph_strength(graph, &instrengths, igraph_vss_all(), IGRAPH_IN, IGRAPH_NO_LOOPS, weights));
        IGRAPH_CHECK(igraph_strength(graph, &outstrengths, igraph_vss_all(), IGRAPH_OUT, IGRAPH_NO_LOOPS, weights));
    } else {
        for (igraph_int_t u = 0; u < no_of_nodes; u++) {
            VECTOR(instrengths)[u] = VECTOR(indegrees)[u];
            VECTOR(outstrengths)[u] = VECTOR(outdegrees)[u];
        }
    }

    /* Find initial sources and sinks */
    nodes_left = no_of_nodes;
    for (igraph_int_t u = 0; u < no_of_nodes; u++) {
        if (VECTOR(indegrees)[u] == 0) {
            if (VECTOR(outdegrees)[u] == 0) {
                /* Isolated vertex, we simply ignore it */
                nodes_left--;
                VECTOR(ordering)[u] = order_next_pos++;
                VECTOR(indegrees)[u] = VECTOR(outdegrees)[u] = -1;
            } else {
                /* This is a source */
                IGRAPH_CHECK(igraph_dqueue_int_push(&sources, u));
            }
        } else if (VECTOR(outdegrees)[u] == 0) {
            /* This is a sink */
            IGRAPH_CHECK(igraph_dqueue_int_push(&sinks, u));
        }
    }

    /* While we have any nodes left... */
    while (nodes_left > 0) {

        /* (1) Remove the sources one by one */
        while (!igraph_dqueue_int_empty(&sources)) {
            const igraph_int_t u = igraph_dqueue_int_pop(&sources);
            /* Add the node to the ordering */
            VECTOR(ordering)[u] = order_next_pos++;
            /* Exclude the node from further searches */
            VECTOR(indegrees)[u] = VECTOR(outdegrees)[u] = -1;
            /* Get the neighbors and decrease their degrees */
            IGRAPH_CHECK(igraph_incident(graph, &neis, u, IGRAPH_OUT, IGRAPH_LOOPS));
            neis_size = igraph_vector_int_size(&neis);
            for (igraph_int_t i = 0; i < neis_size; i++) {
                const igraph_int_t eid = VECTOR(neis)[i];
                const igraph_int_t w = IGRAPH_TO(graph, eid);
                if (VECTOR(indegrees)[w] <= 0) {
                    /* Already removed, continue */
                    continue;
                }
                VECTOR(indegrees)[w]--;
                VECTOR(instrengths)[w] -= (weights ? VECTOR(*weights)[eid] : 1.0);
                if (VECTOR(indegrees)[w] == 0) {
                    IGRAPH_CHECK(igraph_dqueue_int_push(&sources, w));
                }
            }
            nodes_left--;
        }

        /* (2) Remove the sinks one by one */
        while (!igraph_dqueue_int_empty(&sinks)) {
            const igraph_int_t u = igraph_dqueue_int_pop(&sinks);
            /* Maybe the vertex became sink and source at the same time, hence it
             * was already removed in the previous iteration. Check it. */
            if (VECTOR(indegrees)[u] < 0) {
                continue;
            }
            /* Add the node to the ordering */
            VECTOR(ordering)[u] = order_next_neg--;
            /* Exclude the node from further searches */
            VECTOR(indegrees)[u] = VECTOR(outdegrees)[u] = -1;
            /* Get the neighbors and decrease their degrees */
            IGRAPH_CHECK(igraph_incident(graph, &neis, u, IGRAPH_IN, IGRAPH_LOOPS));
            neis_size = igraph_vector_int_size(&neis);
            for (igraph_int_t i = 0; i < neis_size; i++) {
                const igraph_int_t eid = VECTOR(neis)[i];
                const igraph_int_t w = IGRAPH_FROM(graph, eid);
                if (VECTOR(outdegrees)[w] <= 0) {
                    /* Already removed, continue */
                    continue;
                }
                VECTOR(outdegrees)[w]--;
                VECTOR(outstrengths)[w] -= (weights ? VECTOR(*weights)[eid] : 1.0);
                if (VECTOR(outdegrees)[w] == 0) {
                    IGRAPH_CHECK(igraph_dqueue_int_push(&sinks, w));
                }
            }
            nodes_left--;
        }

        /* (3) No more sources or sinks. Find the node with the largest
         * difference between its out-strength and in-strength */
        igraph_int_t v = -1;
        igraph_real_t maxdiff = -IGRAPH_INFINITY;
        for (igraph_int_t u = 0; u < no_of_nodes; u++) {
            if (VECTOR(outdegrees)[u] < 0) {
                continue;
            }
            igraph_real_t diff = VECTOR(outstrengths)[u] - VECTOR(instrengths)[u];
            if (diff > maxdiff) {
                maxdiff = diff;
                v = u;
            }
        }
        if (v >= 0) {
            /* Remove vertex v */
            VECTOR(ordering)[v] = order_next_pos++;
            /* Remove outgoing edges */
            IGRAPH_CHECK(igraph_incident(graph, &neis, v, IGRAPH_OUT, IGRAPH_LOOPS));
            neis_size = igraph_vector_int_size(&neis);
            for (igraph_int_t i = 0; i < neis_size; i++) {
                const igraph_int_t eid = VECTOR(neis)[i];
                const igraph_int_t w = IGRAPH_TO(graph, eid);
                if (VECTOR(indegrees)[w] <= 0) {
                    /* Already removed, continue */
                    continue;
                }
                VECTOR(indegrees)[w]--;
                VECTOR(instrengths)[w] -= (weights ? VECTOR(*weights)[eid] : 1.0);
                if (VECTOR(indegrees)[w] == 0) {
                    IGRAPH_CHECK(igraph_dqueue_int_push(&sources, w));
                }
            }
            /* Remove incoming edges */
            IGRAPH_CHECK(igraph_incident(graph, &neis, v, IGRAPH_IN, IGRAPH_LOOPS));
            neis_size = igraph_vector_int_size(&neis);
            for (igraph_int_t i = 0; i < neis_size; i++) {
                const igraph_int_t eid = VECTOR(neis)[i];
                const igraph_int_t w = IGRAPH_FROM(graph, eid);
                if (VECTOR(outdegrees)[w] <= 0) {
                    /* Already removed, continue */
                    continue;
                }
                VECTOR(outdegrees)[w]--;
                VECTOR(outstrengths)[w] -= (weights ? VECTOR(*weights)[eid] : 1.0);
                if (VECTOR(outdegrees)[w] == 0 && VECTOR(indegrees)[w] > 0) {
                    IGRAPH_CHECK(igraph_dqueue_int_push(&sinks, w));
                }
            }

            VECTOR(outdegrees)[v] = -1;
            VECTOR(indegrees)[v] = -1;
            nodes_left--;
        }
    }

    igraph_vector_destroy(&outstrengths);
    igraph_vector_destroy(&instrengths);
    igraph_vector_int_destroy(&outdegrees);
    igraph_vector_int_destroy(&indegrees);
    igraph_dqueue_int_destroy(&sinks);
    igraph_dqueue_int_destroy(&sources);
    IGRAPH_FINALLY_CLEAN(6);

    /* Tidy up the ordering */
    for (igraph_int_t i = 0; i < no_of_nodes; i++) {
        if (VECTOR(ordering)[i] < 0) {
            VECTOR(ordering)[i] += no_of_nodes;
        }
    }

    /* Find the feedback edges based on the ordering */
    if (result) {
        igraph_vector_int_clear(result);
        for (igraph_int_t eid = 0; eid < no_of_edges; eid++) {
            igraph_int_t from = IGRAPH_FROM(graph, eid);
            igraph_int_t to = IGRAPH_TO(graph, eid);
            if (from == to || VECTOR(ordering)[from] > VECTOR(ordering)[to]) {
                IGRAPH_CHECK(igraph_vector_int_push_back(result, eid));
            }
        }
    }

    /* If we have also requested a layering, return that as well */
    if (layers) {
        igraph_vector_int_t ranks;

        IGRAPH_CHECK(igraph_vector_int_resize(layers, no_of_nodes));
        igraph_vector_int_null(layers);

        IGRAPH_VECTOR_INT_INIT_FINALLY(&ranks, 0);

        IGRAPH_CHECK(igraph_vector_int_sort_ind(&ordering, &ranks, IGRAPH_ASCENDING));

        for (igraph_int_t i = 0; i < no_of_nodes; i++) {
            igraph_int_t from = VECTOR(ranks)[i];
            IGRAPH_CHECK(igraph_neighbors(
                graph, &neis, from, IGRAPH_OUT, IGRAPH_LOOPS, IGRAPH_MULTIPLE
            ));
            neis_size = igraph_vector_int_size(&neis);
            for (igraph_int_t j = 0; j < neis_size; j++) {
                igraph_int_t to = VECTOR(neis)[j];
                if (from == to) {
                    continue;
                }
                if (VECTOR(ordering)[from] > VECTOR(ordering)[to]) {
                    continue;
                }
                if (VECTOR(*layers)[to] < VECTOR(*layers)[from] + 1) {
                    VECTOR(*layers)[to] = VECTOR(*layers)[from] + 1;
                }
            }
        }

        igraph_vector_int_destroy(&ranks);
        IGRAPH_FINALLY_CLEAN(1);
    }

    /* Free the ordering vector */
    igraph_vector_int_destroy(&neis);
    igraph_vector_int_destroy(&ordering);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}


/**
 * Solves the feedback arc set problem with integer programming,
 * using the triangle inequalities formulation.
 */
igraph_error_t igraph_i_feedback_arc_set_ip_ti(
        const igraph_t *graph, igraph_vector_int_t *result,
        const igraph_vector_t *weights) {
#ifndef HAVE_GLPK
    IGRAPH_ERROR("GLPK is not available.", IGRAPH_UNIMPLEMENTED);
#else

    const igraph_int_t no_of_vertices = igraph_vcount(graph);
    const igraph_int_t no_of_edges = igraph_ecount(graph);
    igraph_int_t no_of_components;
    igraph_vector_int_t membership, *vec;
    igraph_vector_int_t ordering, vertex_remapping;
    igraph_vector_int_list_t vertices_by_components, edges_by_components;
    igraph_int_t i, j, k, l, m, n, from, to, no_of_rows, n_choose_2;
    igraph_real_t weight;
    glp_prob *ip;
    glp_iocp parm;

    IGRAPH_VECTOR_INT_INIT_FINALLY(&membership, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&ordering, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&vertex_remapping, no_of_vertices);

    igraph_vector_int_clear(result);

    /* Decompose the graph into connected components */
    IGRAPH_CHECK(igraph_connected_components(graph, &membership, NULL, &no_of_components, IGRAPH_WEAK));

    /* Construct vertex and edge lists for each of the components */
    IGRAPH_VECTOR_INT_LIST_INIT_FINALLY(&vertices_by_components, no_of_components);
    IGRAPH_VECTOR_INT_LIST_INIT_FINALLY(&edges_by_components, no_of_components);
    for (i = 0; i < no_of_vertices; i++) {
        j = VECTOR(membership)[i];
        vec = igraph_vector_int_list_get_ptr(&vertices_by_components, j);
        IGRAPH_CHECK(igraph_vector_int_push_back(vec, i));
    }
    for (i = 0; i < no_of_edges; i++) {
        j = VECTOR(membership)[IGRAPH_FROM(graph, i)];
        vec = igraph_vector_int_list_get_ptr(&edges_by_components, j);
        IGRAPH_CHECK(igraph_vector_int_push_back(vec, i));
    }

#define VAR2IDX(i, j) (i*(n-1)+j-(i+1)*i/2)

    /* Configure GLPK */
    IGRAPH_GLPK_SETUP();
    glp_init_iocp(&parm);
    parm.br_tech = GLP_BR_DTH;
    parm.bt_tech = GLP_BT_BLB;
    parm.pp_tech = GLP_PP_ALL;
    parm.presolve = GLP_ON;
    parm.binarize = GLP_OFF;
    parm.cb_func = igraph_i_glpk_interruption_hook;

    /* Solve an IP for feedback arc sets in each of the components */
    for (i = 0; i < no_of_components; i++) {
        igraph_vector_int_t *vertices_in_comp = igraph_vector_int_list_get_ptr(&vertices_by_components, i);
        igraph_vector_int_t *edges_in_comp = igraph_vector_int_list_get_ptr(&edges_by_components, i);

        /*
         * Let x_ij denote whether layer(i) < layer(j).
         *
         * The standard formulation of the problem is as follows:
         *
         * max sum_{i,j} w_ij x_ij
         *
         * subject to
         *
         * (1) x_ij + x_ji = 1   (i.e. either layer(i) < layer(j) or layer(i) > layer(j))
         *     for all i < j
         * (2) x_ij + x_jk + x_ki <= 2 for all i < j, i < k, j != k
         *
         * Note that x_ij = 1 implies that x_ji = 0 and vice versa; in other words,
         * x_ij = 1 - x_ji. Thus, we can get rid of the (1) constraints and half of the
         * x_ij variables (where j < i) if we rewrite constraints of type (2) as follows:
         *
         * (2a) x_ij + x_jk - x_ik <= 1 for all i < j, i < k, j < k
         * (2b) x_ij - x_kj - x_ik <= 0 for all i < j, i < k, j > k
         *
         * The goal function then becomes:
         *
         * max sum_{i<j} (w_ij-w_ji) x_ij
         */
        n = igraph_vector_int_size(vertices_in_comp);
        ip = glp_create_prob();
        IGRAPH_FINALLY(igraph_i_glp_delete_prob, ip);
        glp_set_obj_dir(ip, GLP_MAX);

        /* Construct a mapping from vertex IDs to the [0; n-1] range */
        for (j = 0; j < n; j++) {
            VECTOR(vertex_remapping)[VECTOR(*vertices_in_comp)[j]] = j;
        }

        /* Set up variables */
        IGRAPH_SAFE_N_CHOOSE_2(n, &n_choose_2);
        if (n_choose_2 > INT_MAX) {
            IGRAPH_ERROR("Feedback arc set problem too large for GLPK.", IGRAPH_EOVERFLOW);
        }

        if (n_choose_2 > 0) {
            glp_add_cols(ip, (int) n_choose_2);
            for (j = 1; j <= n_choose_2; j++) {
                glp_set_col_kind(ip, (int) j, GLP_BV);
            }
        }

        /* Set up coefficients in the goal function */
        k = igraph_vector_int_size(edges_in_comp);
        for (j = 0; j < k; j++) {
            l = VECTOR(*edges_in_comp)[j];
            from = VECTOR(vertex_remapping)[IGRAPH_FROM(graph, l)];
            to = VECTOR(vertex_remapping)[IGRAPH_TO(graph, l)];
            if (from == to) {
                continue;
            }

            weight = weights ? VECTOR(*weights)[l] : 1;

            if (from < to) {
                l = VAR2IDX(from, to);
                glp_set_obj_coef(ip, (int) l, glp_get_obj_coef(ip, (int) l) + weight);
            } else {
                l = VAR2IDX(to, from);
                glp_set_obj_coef(ip, (int) l, glp_get_obj_coef(ip, (int) l) - weight);
            }
        }

        /* Add constraints */
        if (n > 1) {
            {
                /* Overflow-safe block for:
                 *   no_of_rows = n * (n - 1) / 2 + n * (n - 1) * (n - 2) / 3
                 */

                /* res = n * (n - 1) * (n - 2) / 3 */
                igraph_int_t mod = n % 3;
                igraph_int_t res = n / 3; /* same as (n - mod) / 3 */

                mod = (mod + 1) % 3;
                IGRAPH_SAFE_MULT(res, n - mod, &res);
                mod = (mod + 1) % 3;
                IGRAPH_SAFE_MULT(res, n - mod, &res);

                /* no_of_rows = n * (n - 1) / 2 + res */
                IGRAPH_SAFE_ADD(n_choose_2, res, &no_of_rows);
            }
            if (no_of_rows > INT_MAX) {
                IGRAPH_ERROR("Feedback arc set problem too large for GLPK.", IGRAPH_EOVERFLOW);
            }
            glp_add_rows(ip, (int) no_of_rows);
            m = 1;
            for (j = 0; j < n; j++) {
                int ind[4];
                double val[4] = {0, 1, 1, -1};
                for (k = j + 1; k < n; k++) {
                    ind[1] = (int) VAR2IDX(j, k);
                    /* Type (2a) */
                    val[2] = 1;
                    for (l = k + 1; l < n; l++, m++) {
                        ind[2] = (int) VAR2IDX(k, l);
                        ind[3] = (int) VAR2IDX(j, l);
                        glp_set_row_bnds(ip, (int) m, GLP_UP, 1, 1);
                        glp_set_mat_row(ip, (int) m, 3, ind, val);
                    }
                    /* Type (2b) */
                    val[2] = -1;
                    for (l = j + 1; l < k; l++, m++) {
                        ind[2] = (int) VAR2IDX(l, k);
                        ind[3] = (int) VAR2IDX(j, l);
                        glp_set_row_bnds(ip, (int) m, GLP_UP, 0, 0);
                        glp_set_mat_row(ip, (int) m, 3, ind, val);
                    }
                }
            }
        }

        /* Solve the problem */
        IGRAPH_GLPK_CHECK(glp_intopt(ip, &parm),
                          "Feedback arc set using IP with triangle inequalities failed");

        /* Find the ordering of the vertices */
        IGRAPH_CHECK(igraph_vector_int_resize(&ordering, n));
        igraph_vector_int_null(&ordering);
        j = 0; k = 1;
        for (l = 1; l <= n_choose_2; l++) {
            /* variable l always corresponds to the (j, k) vertex pair */
            /* printf("(%ld, %ld) = %g\n", i, j, glp_mip_col_val(ip, l)); */
            if (glp_mip_col_val(ip, (int) l) > 0) {
                /* j comes earlier in the ordering than k */
                VECTOR(ordering)[j]++;
            } else {
                /* k comes earlier in the ordering than j */
                VECTOR(ordering)[k]++;
            }
            k++;
            if (k == n) {
                j++; k = j + 1;
            }
        }

        /* Find the feedback edges */
        k = igraph_vector_int_size(edges_in_comp);
        for (j = 0; j < k; j++) {
            l = VECTOR(*edges_in_comp)[j];
            from = VECTOR(vertex_remapping)[IGRAPH_FROM(graph, l)];
            to = VECTOR(vertex_remapping)[IGRAPH_TO(graph, l)];
            if (from == to || VECTOR(ordering)[from] < VECTOR(ordering)[to]) {
                IGRAPH_CHECK(igraph_vector_int_push_back(result, l));
            }
        }

        /* Clean up */
        glp_delete_prob(ip);
        IGRAPH_FINALLY_CLEAN(1);
    }

    igraph_vector_int_list_destroy(&vertices_by_components);
    igraph_vector_int_list_destroy(&edges_by_components);
    igraph_vector_int_destroy(&vertex_remapping);
    igraph_vector_int_destroy(&ordering);
    igraph_vector_int_destroy(&membership);
    IGRAPH_FINALLY_CLEAN(5);

    return IGRAPH_SUCCESS;
#endif
}


/**
 * Incremental constraint generation based integer programming implementation
 * for feedback arc set (FAS) and feedback vertex set (FVS).
 *
 * b_i are binary variables indicating the presence of edge/vertex i in the
 * FAS/FVS. w_i is the weight of edge/vertex i.
 *
 * We minimize
 *
 * sum_i w_i b_i
 *
 * subject to the constraints
 *
 * sum_i c^k_i b_i >= 1
 *
 * where c^k_i is a binary coefficient indicating if edge/vertex i is present
 * in cycle k.
 *
 * While this must hold for all cycles (all cycles must be broken),
 * we generate cycles incrementally, re-solving the problem after
 * each step. New cycles are generated in such a way as to avoid
 * the feedback set from the previous solution step.
 */

#define VAR_TO_ID(j) ((j) - 1)

/* Helper data structure for adding rows to GLPK problems.
 * ind[] and val[] use one-based indexing, in line with GLPK's convention.
 * Storing the zero-based ind0/val0 and the offset ind/val is necessary
 * to avoid GCC warnings. */
typedef struct {
    int alloc_size;
    int *ind0, *ind;
    double *val0, *val;
} rowdata_t;

static igraph_error_t rowdata_init(rowdata_t *rd, int size) {
    int *ind0 = IGRAPH_CALLOC(size, int);
    IGRAPH_CHECK_OOM(ind0, "Insufficient memory for feedback arc set.");
    IGRAPH_FINALLY(igraph_free, ind0);
    double *val0 = IGRAPH_CALLOC(size, double);
    IGRAPH_CHECK_OOM(val0, "Insufficient memory for feedback arc set.");
    for (int i=0; i < size; i++) {
        val0[i] = 1.0;
    }
    rd->alloc_size = size;
    rd->ind0 = ind0;
    rd->ind = ind0 - 1;
    rd->val0 = val0;
    rd->val = val0 - 1;
    IGRAPH_FINALLY_CLEAN(1);
    return IGRAPH_SUCCESS;
}

static igraph_error_t rowdata_set(rowdata_t *rd, const igraph_vector_int_t *idx) {
    int size = igraph_vector_int_size(idx);

    /* Expand size if needed */
    if (size > rd->alloc_size) {
        int new_alloc_size = 2 * rd->alloc_size;
        if (size > new_alloc_size) {
            new_alloc_size = size;
        }

        int *ind0 = rd->ind0;
        double *val0 = rd->val0;

        ind0 = IGRAPH_REALLOC(ind0, new_alloc_size, int);
        IGRAPH_CHECK_OOM(ind0, "Insufficient memory for feedback arc set.");
        rd->ind0 = ind0;
        rd->ind = ind0 - 1;

        val0 = IGRAPH_REALLOC(val0, new_alloc_size, double);
        IGRAPH_CHECK_OOM(val0, "Insufficient memory for feedback arc set.");
        for (int i = rd->alloc_size; i < new_alloc_size; i++) {
            val0[i] = 1.0;
        }
        rd->val0 = val0;
        rd->val = val0 - 1;

        rd->alloc_size = new_alloc_size;
    }

    for (int i = 0; i < size; i++) {
        rd->ind0[i] = VECTOR(*idx)[i] + 1;
    }

    return IGRAPH_SUCCESS;
}

static void rowdata_destroy(rowdata_t *rd) {
    igraph_free(rd->ind0);
    igraph_free(rd->val0);
}

igraph_error_t igraph_i_feedback_arc_set_ip_cg(
        const igraph_t *graph, igraph_vector_int_t *result,
        const igraph_vector_t *weights) {

#ifndef HAVE_GLPK
    IGRAPH_ERROR("GLPK is not available.", IGRAPH_UNIMPLEMENTED);
#else
    const igraph_int_t ecount = igraph_ecount(graph);
    igraph_bool_t is_dag;
    igraph_bitset_t removed;
    igraph_vector_int_t cycle;
    glp_prob *ip;
    glp_iocp parm;
    rowdata_t rd;
    int var_count;

    /* Avoid starting up the IP machinery for DAGs. */
    IGRAPH_CHECK(igraph_is_dag(graph, &is_dag));
    if (is_dag) {
        igraph_vector_int_clear(result);
        return IGRAPH_SUCCESS;
    }

    if (ecount > INT_MAX) {
        IGRAPH_ERROR("Feedback arc set problem too large for GLPK.", IGRAPH_EOVERFLOW);
    }

    var_count = (int) ecount;

    /* TODO: In-depth investigation of whether decomposing to SCCs helps performance.
     * Basic benchmarking on sparse random graphs with mean degrees between 1-2
     * indicate no benefit from avoiding creating GLPK variables for non-cycle edges. */

    IGRAPH_BITSET_INIT_FINALLY(&removed, ecount);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&cycle, 0);

    IGRAPH_CHECK(rowdata_init(&rd, 20));
    IGRAPH_FINALLY(rowdata_destroy, &rd);

    /* Configure GLPK */
    IGRAPH_GLPK_SETUP();
    glp_init_iocp(&parm);
    parm.br_tech = GLP_BR_MFV;
    parm.bt_tech = GLP_BT_BLB;
    parm.pp_tech = GLP_PP_ALL;
    parm.presolve = GLP_ON;
    parm.cb_func = igraph_i_glpk_interruption_hook;

    ip = glp_create_prob();
    IGRAPH_FINALLY(igraph_i_glp_delete_prob, ip);

    glp_set_obj_dir(ip, GLP_MIN);

    glp_add_cols(ip, var_count);
    for (int j = 1; j <= var_count; j++) {
        glp_set_obj_coef(ip, j, weights ? VECTOR(*weights)[ VAR_TO_ID(j) ] : 1);
        glp_set_col_kind(ip, j, GLP_BV);
    }

    while (true) {
        int cycle_size, row;

        IGRAPH_CHECK(igraph_i_find_cycle(graph, NULL, &cycle, NULL, IGRAPH_OUT, &removed));

        cycle_size = (int) igraph_vector_int_size(&cycle);

        if (cycle_size == 0) break; /* no more cycles, we're done */

        IGRAPH_CHECK(rowdata_set(&rd, &cycle));

        row = glp_add_rows(ip, 1);
        glp_set_row_bnds(ip, row, GLP_LO, 1, 0);
        glp_set_mat_row(ip, row, cycle_size, rd.ind, rd.val);

        /* Add as many edge-disjoint cycles at once as possible. */
        while (true) {
            for (int i=0; i < cycle_size; i++) {
                IGRAPH_BIT_SET(removed, VECTOR(cycle)[i]);
            }
            IGRAPH_CHECK(igraph_i_find_cycle(graph, NULL, &cycle, NULL, IGRAPH_OUT, &removed));

            cycle_size = (int) igraph_vector_int_size(&cycle);
            if (cycle_size == 0) break; /* no more edge disjoint cycles */

            IGRAPH_CHECK(rowdata_set(&rd, &cycle));

            row = glp_add_rows(ip, 1);
            glp_set_row_bnds(ip, row, GLP_LO, 1, 0);
            glp_set_mat_row(ip, row, cycle_size, rd.ind, rd.val);
        }

        IGRAPH_GLPK_CHECK(glp_intopt(ip, &parm),
                          "Feedback arc set using IP with incremental cycle generation failed");

        igraph_vector_int_clear(result);
        igraph_bitset_null(&removed);
        for (int j=1; j <= var_count; j++) {
            if (glp_mip_col_val(ip, j) > 0) {
                igraph_int_t i = VAR_TO_ID(j);
                IGRAPH_CHECK(igraph_vector_int_push_back(result, i));
                IGRAPH_BIT_SET(removed, i);
            }
        }
    }

    /* Clean up */
    glp_delete_prob(ip);
    rowdata_destroy(&rd);
    igraph_vector_int_destroy(&cycle);
    igraph_bitset_destroy(&removed);
    IGRAPH_FINALLY_CLEAN(4);

    return IGRAPH_SUCCESS;
#endif
}


igraph_error_t igraph_i_feedback_vertex_set_ip_cg(
        const igraph_t *graph, igraph_vector_int_t *result,
        const igraph_vector_t *vertex_weights) {
#ifndef HAVE_GLPK
    IGRAPH_ERROR("GLPK is not available.", IGRAPH_UNIMPLEMENTED);
#else

    const igraph_int_t vcount = igraph_vcount(graph);
    const igraph_int_t ecount = igraph_ecount(graph);
    igraph_bool_t is_acyclic;
    igraph_bitset_t removed;
    igraph_vector_int_t cycle;
    igraph_vector_int_t incident;
    glp_prob *ip;
    glp_iocp parm;
    rowdata_t rd;
    int var_count;

    /* Avoid starting up the IP machinery for acyclic graphs. */
    IGRAPH_CHECK(igraph_is_acyclic(graph, &is_acyclic));

    if (is_acyclic) {
        igraph_vector_int_clear(result);
        return IGRAPH_SUCCESS;
    }

    if (vcount > INT_MAX) {
        IGRAPH_ERROR("Feedback vertex set problem too large for GLPK.", IGRAPH_EOVERFLOW);
    }

    var_count = (int) vcount;

    IGRAPH_BITSET_INIT_FINALLY(&removed, ecount);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&cycle, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&incident, 0);
    IGRAPH_CHECK(rowdata_init(&rd, 20));
    IGRAPH_FINALLY(rowdata_destroy, &rd);

    /* Configure GLPK */
    IGRAPH_GLPK_SETUP();
    glp_init_iocp(&parm);
    parm.br_tech = GLP_BR_MFV;
    parm.bt_tech = GLP_BT_BLB;
    parm.pp_tech = GLP_PP_ALL;
    parm.presolve = GLP_ON;
    parm.cb_func = igraph_i_glpk_interruption_hook;

    ip = glp_create_prob();
    IGRAPH_FINALLY(igraph_i_glp_delete_prob, ip);

    glp_set_obj_dir(ip, GLP_MIN);

    glp_add_cols(ip, var_count);
    for (int j = 1; j <= var_count; j++) {
        glp_set_obj_coef(ip, j, vertex_weights ? VECTOR(*vertex_weights)[ VAR_TO_ID(j) ] : 1);
        glp_set_col_kind(ip, j, GLP_BV);
    }

    while (true) {
        int cycle_size, row;

        IGRAPH_CHECK(igraph_i_find_cycle(graph, &cycle, NULL, NULL, IGRAPH_OUT, &removed));

        cycle_size = (int) igraph_vector_int_size(&cycle);

        if (cycle_size == 0) break; /* no more cycles, we're done */

        IGRAPH_CHECK(rowdata_set(&rd, &cycle));

        row = glp_add_rows(ip, 1);
        glp_set_row_bnds(ip, row, GLP_LO, 1, 0);
        glp_set_mat_row(ip, row, cycle_size, rd.ind, rd.val);

        /* Add as many vertex-disjoint cycles at once as possible. */
        while (true) {
            for (int i=0; i < cycle_size; i++) {
                IGRAPH_CHECK(igraph_incident(graph, &incident, VECTOR(cycle)[i], IGRAPH_ALL, IGRAPH_LOOPS));
                const igraph_int_t incident_size = igraph_vector_int_size(&incident);
                for (igraph_int_t j = 0; j < incident_size; j++) {
                    igraph_int_t eid = VECTOR(incident)[j];
                    IGRAPH_BIT_SET(removed, eid);
                }
            }
            IGRAPH_CHECK(igraph_i_find_cycle(graph, &cycle, NULL, NULL, IGRAPH_OUT, &removed));

            cycle_size = (int) igraph_vector_int_size(&cycle);
            if (cycle_size == 0) break; /* no more vertex disjoint cycles */

            IGRAPH_CHECK(rowdata_set(&rd, &cycle));

            row = glp_add_rows(ip, 1);
            glp_set_row_bnds(ip, row, GLP_LO, 1, 0);
            glp_set_mat_row(ip, row, cycle_size, rd.ind, rd.val);
        }

        IGRAPH_GLPK_CHECK(glp_intopt(ip, &parm),
                          "Feedback vertex set using IP with incremental cycle generation failed");

        igraph_vector_int_clear(result);
        igraph_bitset_null(&removed);

        for (int j=1; j <= var_count; j++) {
            if (glp_mip_col_val(ip, j) > 0) {
                igraph_int_t i = VAR_TO_ID(j);
                IGRAPH_CHECK(igraph_vector_int_push_back(result, i));

                IGRAPH_CHECK(igraph_incident(graph, &incident, i, IGRAPH_ALL, IGRAPH_LOOPS));

                const igraph_int_t incident_size = igraph_vector_int_size(&incident);
                for (igraph_int_t k = 0; k < incident_size; k++) {
                    igraph_int_t eid = VECTOR(incident)[k];
                    IGRAPH_BIT_SET(removed, eid);
                }
            }
        }
    }

    /* Clean up */
    glp_delete_prob(ip);
    rowdata_destroy(&rd);
    igraph_vector_int_destroy(&cycle);
    igraph_vector_int_destroy(&incident);
    igraph_bitset_destroy(&removed);
    IGRAPH_FINALLY_CLEAN(5);

    return IGRAPH_SUCCESS;
#endif
}

#undef VAR_TO_ID
