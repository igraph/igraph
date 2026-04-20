/*
   IGraph library.
   Copyright (C) 2025  The igraph development team

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

#include "igraph_vector_list.h"

static igraph_integer_t igraph_i_number_of_subsets(igraph_integer_t n, igraph_integer_t k) {
    igraph_integer_t r = (k < n - k ? k : n - k), number_of_subsets = 1;

    for(igraph_integer_t i = 1 ; i <= r ; i++) {
        number_of_subsets *= n - i + 1;
        number_of_subsets /= i;
    }

    return number_of_subsets;
}

static igraph_error_t igraph_i_generate_subsets(igraph_integer_t no_of_subsets, igraph_integer_t k,
                                                igraph_bitset_list_t *vertices) {
    igraph_integer_t i;
    igraph_bitset_t *vertex = igraph_bitset_list_get_ptr(vertices, 0);
    igraph_vector_int_t subset;
    IGRAPH_CHECK(igraph_vector_int_init(&subset, k));

    for(igraph_integer_t i = 0; i < k; i++) {
        VECTOR(subset)[i] = i;
        IGRAPH_BIT_SET(*vertex, i);
    }

    i = 1;
    while(i < no_of_subsets) {
        igraph_integer_t j = k - 1;
        vertex = igraph_bitset_list_get_ptr(vertices, i);

        while(j >= 0 && VECTOR(subset)[j] == n - (k - 1 - j)) {
            j--;
        }

        VECTOR(subset)[j]++;
        while(++j < k) {
            VECTOR(subset)[j] = VECTOR(subset)[j - 1] + 1;
        }
        for(igraph_integer_t pos = 0; pos < k; pos++) {
            IGRAPH_BIT_SET(*vertex, VECTOR(subset)[pos]);
        }

        i++;
    }

    igraph_vector_int_destroy(&subset);
    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_kneser
 * \brief Generate a Kneser graph.
 *
 * A Kneser graph is a graph in which each node is a \c k element subset of
 * a set of n elements. Two vertices representing two distinct subsets are
 * connected by an edge only if the two subsets have no element common between
 * them.
 *
 * </para><para>
 * Please note that the graph will have \c nCk vertices and edges of the same
 * order, so probably you don't want to supply too big numbers for
 * \c m and \c n.
 *
 * </para><para>
 * Kneser graphs are a well-studied graph family in graph theory and have some
 * interesting properties, see e.g. Wikipedia for details.
 *
 * \param graph Pointer to an uninitialized graph object, the result will be
 *        stored here.
 * \param n Integer, the total number of elements
 * \param k Integer, the size of each subset
 * \return Error code.
 *
 * Time complexity: O((nCk)^2), the number of vertices plus the number of edges.
 */
igraph_error_t igraph_kneser(igraph_t *graph, igraph_int_t n, igraph_int_t k) {
    igraph_int_t no_of_nodes;
    igraph_bitset_list_t nodes;
    igraph_vector_int_t edges;

    if (n < 0 || k < 0) {
        IGRAPH_ERROR("`n' and `k' should be positive in a Kneser graph",
                     IGRAPH_EINVAL);
    }

    if (n < k) {
        IGRAPH_ERROR("`n' must be at least as large as `k' in a Kneser graph",
                     IGRAPH_EINVAL);
    }

    if (n != 2 && n == k + 1) {
        return igraph_empty(graph, n, IGRAPH_UNDIRECTED);
    }

    no_of_nodes = igraph_i_number_of_subsets(n, k);
    IGRAPH_BITSET_LIST_INIT_FINALLY(&nodes, no_of_nodes);
    for (igraph_int_t i = 0; i < no_of_nodes; ++i) {
        IGRAPH_CHECK(igraph_bitset_resize(igraph_bitset_list_get_ptr(&nodes, i), n));
    }

    igraph_i_generate_subsets(no_of_nodes, k, &nodes);

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);

    for (igraph_int_t i = 0; i < no_of_nodes; ++i) {
        for (igraph_int_t j = i + 1; j < no_of_nodes; ++j) {
            igraph_bitset_t result;
            igraph_bitset_and(&result, igraph_bitset_list_get_ptr(&nodes, i),
                              igraph_bitset_list_get_ptr(&nodes, j));
            if (igraph_bitset_is_all_zero(&result)) {
                IGRAPH_CHECK(igraph_vector_int_push_back(&edges, i));
                IGRAPH_CHECK(igraph_vector_int_push_back(&edges, j));
            }
        }
    }

    IGRAPH_CHECK(igraph_create(graph, &edges, no_of_nodes, IGRAPH_UNDIRECTED));

    igraph_vector_int_destroy(&edges);
    igraph_bitset_list_destroy(&nodes);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}
