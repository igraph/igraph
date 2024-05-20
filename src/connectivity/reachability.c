/*
   IGraph library.
   Copyright (C) 2024  The igraph development team <igraph@igraph.org>

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

#include "igraph_reachability.h"

#include "igraph_adjlist.h"
#include "igraph_bitset_list.h"
#include "igraph_components.h"
#include "igraph_constructors.h"
#include "igraph_interface.h"


/**
 * \ingroup structural
 * \function igraph_reachability
 * \brief Calculates which vertices are reachable from each vertex in the graph.
 *
 * \experimental
 *
 * The resulting list will contain one bitset for each strongly connected component.
 * The bitset for component i will have its j-th bit set, if vertex j is reachable
 * from some vertex in component i in 0 or more steps.
 * In particular, a vertex is always reachable from itself.
 *
 * \param graph The graph object to analyze.
 * \param membership Pointer to an integer vector. For every vertex,
 *    the ID of its component is given. The vector will be resized as needed.
 *    This parameter must not be \c NULL.
 * \param csize Pointer to an integer vector or \c NULL. For every component, it
 *    gives its size (vertex count), the order being defined by the component
 *    IDs. The vector will be resized as needed.
 * \param no_of_components Pointer to an integer or \c NULL. The number of
 *    components will be stored here.
 * \param reach A list of bitsets representing the result. It will be resized
 *    as needed. <code>reach[membership[u]][v]</code> is set to \c true if
 *    vertex \c v is reachable from vertex \c u.
 * \param mode In directed graphs, controls the treatment of edge directions.
 *    Ignored in undirected graphs. With \c IGRAPH_OUT, reachability is computed
 *    by traversing edges along their direction. With \c IGRAPH_IN, edges are
 *    traversed opposite to their direction. With \c IGRAPH_ALL, edge directions
 *    are ignored and the graph is treated as undirected.
 * \return Error code:
 *         \c IGRAPH_ENOMEM if there is not enough memory
 *         to perform the operation.
 *
 * \sa \ref igraph_connected_components() to find the connnected components
 * of a graph; \ref igraph_count_reachable() to count how many vertices
 * are reachable from each vertex; \ref igraph_subcomponent() to find
 * which vertices are rechable from a single vertex.
 *
 * Time complexity: O(|C||V|/w + |V| + |E|), where
 * |C| is the number of strongly connected components (at most |V|),
 * |V| is the number of vertices, and
 * |E| is the number of edges respectively,
 * and w is the bit width of \type igraph_integer_t, typically the
 * word size of the machine (32 or 64).
 */

igraph_error_t igraph_reachability(
        const igraph_t *graph,
        igraph_vector_int_t *membership,
        igraph_vector_int_t *csize,
        igraph_integer_t *no_of_components,
        igraph_bitset_list_t *reach,
        igraph_neimode_t mode) {

    const igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_comps;
    igraph_adjlist_t adjlist, dag;

    if (mode != IGRAPH_ALL && mode != IGRAPH_OUT && mode != IGRAPH_IN) {
        IGRAPH_ERROR("Invalid mode for reachability.", IGRAPH_EINVMODE);
    }

    if (! igraph_is_directed(graph)) {
        mode = IGRAPH_ALL;
    }

    IGRAPH_CHECK(igraph_connected_components(graph,
                                             membership, csize, &no_of_comps,
                                             mode == IGRAPH_ALL ? IGRAPH_WEAK : IGRAPH_STRONG));

    if (no_of_components) {
        *no_of_components = no_of_comps;
    }

    IGRAPH_CHECK(igraph_bitset_list_resize(reach, no_of_comps));

    for (igraph_integer_t comp = 0; comp < no_of_comps; comp++) {
        IGRAPH_CHECK(igraph_bitset_resize(igraph_bitset_list_get_ptr(reach, comp), no_of_nodes));
    }
    for (igraph_integer_t v = 0; v < no_of_nodes; v++) {
        IGRAPH_BIT_SET(*igraph_bitset_list_get_ptr(reach, VECTOR(*membership)[v]), v);
    }

    if (mode == IGRAPH_ALL) {
        return IGRAPH_SUCCESS;
    }

    IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist, mode, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist);

    IGRAPH_CHECK(igraph_adjlist_init_empty(&dag, no_of_comps));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &dag);

    for (igraph_integer_t v = 0; v < no_of_nodes; v++) {
        const igraph_vector_int_t *neighbours = igraph_adjlist_get(&adjlist, v);
        igraph_vector_int_t *dag_neighbours = igraph_adjlist_get(&dag, VECTOR(*membership)[v]);
        const igraph_integer_t n = igraph_vector_int_size(neighbours);
        for (igraph_integer_t i = 0; i < n; i++) {
            igraph_integer_t w = VECTOR(*neighbours)[i];
            if (VECTOR(*membership)[v] != VECTOR(*membership)[w]) {
                IGRAPH_CHECK(igraph_vector_int_push_back(dag_neighbours, VECTOR(*membership)[w]));
            }
        }
    }

    /* Iterate through strongly connected components in reverser topological order,
     * exploiting the fact that they are indexed in topological order. */
    for (igraph_integer_t i = 0; i < no_of_comps; i++) {
        const igraph_integer_t comp = mode == IGRAPH_IN ? i : no_of_comps - i - 1;
        const igraph_vector_int_t *dag_neighbours = igraph_adjlist_get(&dag, comp);
        igraph_bitset_t *from_bitset = igraph_bitset_list_get_ptr(reach, comp);
        const igraph_integer_t n = igraph_vector_int_size(dag_neighbours);
        for (igraph_integer_t j = 0; j < n; j++) {
            const igraph_bitset_t *to_bitset = igraph_bitset_list_get_ptr(reach, VECTOR(*dag_neighbours)[j]);
            igraph_bitset_or(from_bitset, from_bitset, to_bitset);
        }
    }

    igraph_adjlist_destroy(&adjlist);
    igraph_adjlist_destroy(&dag);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}


/**
 * \ingroup structural
 * \function igraph_count_reachable
 * \brief The number of vertices reachable from each vertex in the graph.
 *
 * \experimental
 *
 * \param graph The graph object to analyze.
 * \param counts Integer vector. <code>counts[v]</code> will store the number
 *    of vertices reachable from vertex \c v, including \c v itself.
 * \param mode In directed graphs, controls the treatment of edge directions.
 *    Ignored in undirected graphs. With \c IGRAPH_OUT, reachability is computed
 *    by traversing edges along their direction. With \c IGRAPH_IN, edges are
 *    traversed opposite to their direction. With \c IGRAPH_ALL, edge directions
 *    are ignored and the graph is treated as undirected.
 * \return Error code:
 *         \c IGRAPH_ENOMEM if there is not enough memory
 *         to perform the operation.
 *
 * \sa \ref igraph_connected_components(), \ref igraph_transitive_closure()
 *
 * Time complexity: O(|C||V|/w + |V| + |E|), where
 * |C| is the number of strongly connected components (at most |V|),
 * |V| is the number of vertices, and
 * |E| is the number of edges respectively,
 * and w is the bit width of \type igraph_integer_t, typically the
 * word size of the machine (32 or 64).
 */

igraph_error_t igraph_count_reachable(const igraph_t *graph,
                                      igraph_vector_int_t *counts,
                                      igraph_neimode_t mode) {

    igraph_vector_int_t membership;
    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_bitset_list_t reach;

    IGRAPH_VECTOR_INT_INIT_FINALLY(&membership, 0);
    IGRAPH_BITSET_LIST_INIT_FINALLY(&reach, 0);

    IGRAPH_CHECK(igraph_reachability(graph, &membership, NULL, NULL, &reach, mode));

    IGRAPH_CHECK(igraph_vector_int_resize(counts, igraph_vcount(graph)));
    for (igraph_integer_t i = 0; i < no_of_nodes; i++) {
        VECTOR(*counts)[i] = igraph_bitset_popcount(igraph_bitset_list_get_ptr(&reach, VECTOR(membership)[i]));
    }

    igraph_bitset_list_destroy(&reach);
    igraph_vector_int_destroy(&membership);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}


/**
 * \ingroup structural
 * \function igraph_transitive_closure
 * \brief Computes the transitive closure of a graph.
 *
 * \experimental
 *
 * The resulting graph will have an edge from vertex \c i to vertex \c j
 * if \c j is reachable from \c i.
 *
 * \param graph The graph object to analyze.
 * \param closure The resulting graph representing the transitive closure.
 * \return Error code:
 *         \c IGRAPH_ENOMEM if there is not enough memory
 *         to perform the operation.
 *
 * \sa \ref igraph_connected_components(), \ref igraph_count_reachable()
 *
 * Time complexity: O(|V|^2 + |E|), where
 * |V| is the number of vertices, and
 * |E| is the number of edges, respectively.
 */
igraph_error_t igraph_transitive_closure(const igraph_t *graph, igraph_t *closure) {
    const igraph_integer_t no_of_nodes = igraph_vcount(graph);
    const igraph_bool_t directed = igraph_is_directed(graph);
    igraph_vector_int_t membership, edges;
    igraph_bitset_list_t reach;

    IGRAPH_VECTOR_INT_INIT_FINALLY(&membership, 0);
    IGRAPH_BITSET_LIST_INIT_FINALLY(&reach, 0);

    IGRAPH_CHECK(igraph_reachability(graph, &membership, NULL, NULL, &reach, IGRAPH_OUT));

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);
    for (igraph_integer_t u = 0; u < no_of_nodes; u++) {
        const igraph_bitset_t *row = igraph_bitset_list_get_ptr(&reach, VECTOR(membership)[u]);
        for (igraph_integer_t v = directed ? 0 : u + 1; v < no_of_nodes; v++) {
            if (u != v && IGRAPH_BIT_TEST(*row, v)) {
                IGRAPH_CHECK(igraph_vector_int_push_back(&edges, u));
                IGRAPH_CHECK(igraph_vector_int_push_back(&edges, v));
            }
        }
    }

    igraph_bitset_list_destroy(&reach);
    igraph_vector_int_destroy(&membership);
    IGRAPH_FINALLY_CLEAN(2);

    IGRAPH_CHECK(igraph_create(closure, &edges, no_of_nodes, directed));

    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}
