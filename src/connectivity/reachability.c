/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
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
 * </para><para>
 * The resulting list will contain one bitset for each strongly connected component.
 * The bitset for component i will have its j-th bit set, if vertex j is reachable from some vertex in component i in 0 or more steps.
 * In particular, a vertex is always reachable from itself.
 * Use the membership vector to get the reachability for a particular vertex.
 *
 * \param graph The graph object to analyze.
 * \param membership See \ref igraph_connected_components .
 * \param csize See \ref igraph_connected_components .
 * \param no_of_components Pointer to an integer. The number of
 *        components will be stored here.
 * \param reach A list of bitsets representing the result.
 * \param directed For a directed graph, determines whether edges are treated as directed or undirected.
 * \return Error code:
 *         \c IGRAPH_ENOMEM if there is not enough memory
 *         to perform the operation.
 *
 * Time complexity: O(|C||V|/w + |V| + |E|),
 * |C| is the number of strongly connected components (at most |V|),
 * |V| is the number of vertices,
 * |E| is the number of edges and
 * edges in the graph, respectively and
 * w is the word size of the machine (32 or 64).
 */

igraph_error_t igraph_reachability(
    const igraph_t *graph,
    igraph_vector_int_t *membership,
    igraph_vector_int_t *csize,
    igraph_integer_t *no_of_components,
    igraph_bitset_list_t *reach,
    igraph_bool_t directed)
{
    igraph_integer_t no_of_nodes;
    igraph_integer_t i, j, n;
    igraph_adjlist_t adjlist, dag;
    igraph_vector_int_t *dag_neighbours, *neighbours;
    igraph_bitset_t *from_bitset, *to_bitset;

    no_of_nodes = igraph_vcount(graph);

    if (directed) {
        IGRAPH_CHECK(igraph_connected_components(graph, membership, csize, no_of_components, IGRAPH_STRONG));
    }
    else {
        IGRAPH_CHECK(igraph_connected_components(graph, membership, csize, no_of_components, IGRAPH_WEAK));
    }

    IGRAPH_CHECK(igraph_bitset_list_resize(reach, *no_of_components));

    for (i = 0; i < *no_of_components; ++i)
    {
        IGRAPH_CHECK(igraph_bitset_resize(igraph_bitset_list_get_ptr(reach, i), no_of_nodes));
    }
    for (i = 0; i < no_of_nodes; ++i)
    {
        IGRAPH_BIT_SET(*igraph_bitset_list_get_ptr(reach, VECTOR(*membership)[i]), i);
    }

    if (!directed) {
        return IGRAPH_SUCCESS;
    }

    IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist, IGRAPH_OUT, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist);
    IGRAPH_CHECK(igraph_adjlist_init_empty(&dag, *no_of_components));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &dag);

    for (i = 0; i < no_of_nodes; ++i)
    {
        neighbours = igraph_adjlist_get(&adjlist, i);
        dag_neighbours = igraph_adjlist_get(&dag, VECTOR(*membership)[i]);
        n = igraph_vector_int_size(neighbours);
        for (j = 0; j < n; j++)
        {
            if (VECTOR(*membership)[i] != VECTOR(*membership)[VECTOR(*neighbours)[j]])
            {
                igraph_vector_int_push_back(dag_neighbours, VECTOR(*membership)[VECTOR(*neighbours)[j]]);
            }
        }
    }

    for (i = *no_of_components - 1; i >= 0; --i)
    {
        dag_neighbours = igraph_adjlist_get(&dag, i);
        from_bitset = igraph_bitset_list_get_ptr(reach, i);
        n = igraph_vector_int_size(dag_neighbours);
        for (j = 0; j < n; j++)
        {
            to_bitset = igraph_bitset_list_get_ptr(reach, VECTOR(*dag_neighbours)[j]);
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
 * \brief Calculates the number of vertices reachable from each vertex in the graph.
 *
 * \experimental
 * </para><para>
 * The resulting vector will store how many vertices are reachable from vertex i at index i.
 *
 * \param graph The graph object to analyze.
 * \param counts A vector of integers representing the result.
 * \param directed For a directed graph, determines whether edges are treated as directed or undirected.
 * \return Error code:
 *         \c IGRAPH_ENOMEM if there is not enough memory
 *         to perform the operation.
 *
 * Time complexity: O(|C||V|/w + |V| + |E|),
 * |C| is the number of strongly connected components (at most |V|),
 * |V| is the number of vertices,
 * |E| is the number of edges and
 * edges in the graph, respectively and
 * w is the word size of the machine (32 or 64).
 */

igraph_error_t igraph_count_reachable(
    const igraph_t *graph,
    igraph_vector_int_t *counts,
    igraph_bool_t directed)
{
    igraph_vector_int_t membership, csize;
    igraph_integer_t no_of_components, no_of_nodes = igraph_vcount(graph);
    igraph_bitset_list_t reach;

    IGRAPH_CHECK(igraph_vector_int_init(&membership, 0));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &membership);

    IGRAPH_CHECK(igraph_vector_int_init(&csize, 0));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &csize);

    IGRAPH_CHECK(igraph_bitset_list_init(&reach, 0));
    IGRAPH_FINALLY(igraph_bitset_list_destroy, &reach);

    IGRAPH_CHECK(igraph_reachability(graph, &membership, &csize, &no_of_components, &reach, directed));

    IGRAPH_CHECK(igraph_vector_int_resize(counts, igraph_vcount(graph)));
    for (igraph_integer_t i = 0; i < no_of_nodes; i++)
    {
        VECTOR(*counts)[i] = igraph_bitset_popcount(igraph_bitset_list_get_ptr(&reach, VECTOR(membership)[i]));
    }

    igraph_bitset_list_destroy(&reach);
    igraph_vector_int_destroy(&csize);
    igraph_vector_int_destroy(&membership);
    IGRAPH_FINALLY_CLEAN(3);

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
 * Time complexity: O(|V|^2 + |E|), where
 * |V| is the number of vertices, and
 * |E| is the number of edges, respectively.
 */
igraph_error_t igraph_transitive_closure(const igraph_t *graph, igraph_t *closure) {
    const igraph_integer_t no_of_nodes = igraph_vcount(graph);
    const igraph_bool_t directed = igraph_is_directed(graph);
    igraph_vector_int_t membership, csize, edges;
    igraph_integer_t no_of_components;
    igraph_bitset_list_t reach;

    IGRAPH_VECTOR_INT_INIT_FINALLY(&membership, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&csize, 0);
    IGRAPH_BITSET_LIST_INIT_FINALLY(&reach, 0);

    IGRAPH_CHECK(igraph_reachability(graph, &membership, &csize, &no_of_components, &reach, /* directed */ true));

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

    IGRAPH_CHECK(igraph_create(closure, &edges, no_of_nodes, directed));

    igraph_bitset_list_destroy(&reach);
    igraph_vector_int_destroy(&csize);
    igraph_vector_int_destroy(&membership);
    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(4);

    return IGRAPH_SUCCESS;
}
