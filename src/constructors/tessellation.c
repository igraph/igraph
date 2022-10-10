/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2005-2022 The igraph development team

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

#include "igraph_constructors.h"

/**
 * \ingroup generators
 * \function igraph_triangulated_mesh
 * \brief A triangulated mesh with the given shape.
 *
 * Creates a triangulated mesh whose vertices have the form (i, j) for non-negative integers i and j
 * and (i, j) is generally connected with (i + 1, j), (i, j + 1), and (i - 1, j + 1).
 *
 *
 * </para><para>
 * In the zero-dimensional case, the singleton graph is returned.
 *
 * </para><para>
 * The vertices of the resulting graph are ordered lexicographically with the 2nd coordinate being
 *  more significant, e.g., (i, j) < (i + 1, j) and (i + 1, j) < (i, j + 1)
 *
 * \param graph An uninitialized graph object.
 * \param directed Boolean, whether to create a directed graph.
 *        If the \c mutual argument is  not set to true,
 *        edges will be directed from lower-index vertices towards
 *        higher-index ones.
 * \param mutual Boolean, if the graph is directed this gives whether
 *        to create all connections as mutual.
 * \param row_lengths_vector Integer vector, defines the number of vertices with
 *        the second coordinate equal to the index. The length of this vector must match
 *        the length of \p row_start_vector.
 * \param row_start_vector Integer vector, defines the leftmost coordinate of
 *        the vertex with the second coordinate equal to the index.
 * \param is_vertex_in_graph A funtion, a predicate that takes coordinates of a vertex (i, j) and some additional \p condition_params.
 *        The predicate is evaluated to true if and only if the vertex corresponding to (i, j) belongs to the graph being constructed.
 *        Typically, composed of linear inqualities.
 *        Can be deduced from \p row_lengths_vector and \p row_start_vector.
 * \param condition_params Integer vector, that contains parameters supplied to \p is_vertex_in_graph.
 *
 * \return Error code:
 *         \c IGRAPH_EINVAL: invalid (negative) length of row_lengths_vector does not match the length of the
            row_start_vector.
 *
 * Time complexity:  O(|V|), where |V| is the number of vertices in the generated graph.
 */
igraph_error_t igraph_triangulated_mesh(igraph_t *graph, igraph_bool_t directed, igraph_bool_t mutual, igraph_integer_t row_count,
                                        igraph_vector_int_t row_lengths_vector, igraph_vector_int_t row_start_vector,
                                        igraph_bool_t(is_vertex_in_graph)(igraph_integer_t, igraph_integer_t, igraph_vector_int_t),
                                        igraph_vector_int_t condition_params)
{
    igraph_vector_int_t edges = IGRAPH_VECTOR_NULL;
    igraph_integer_t no_of_nodes;
    igraph_vector_int_t row_lengths_prefix_sum_vector;
    igraph_integer_t i, j;

    if (igraph_vector_int_size(&row_lengths_vector) != igraph_vector_int_size(&row_start_vector))
    {
        IGRAPH_ERRORF(
            "Length of row_lengths_vector vector must match the length of the "
            "row_start_vector (%" IGRAPH_PRId ").",
            IGRAPH_EINVAL, igraph_vector_int_size(&row_start_vector));
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);
    igraph_vector_int_init(&row_lengths_prefix_sum_vector, row_count + 1);
    VECTOR(row_lengths_prefix_sum_vector)[0] = 0;
    for (i = 1; i < row_count + 1; i++)
    {
        VECTOR(row_lengths_prefix_sum_vector)[i] = VECTOR(row_lengths_prefix_sum_vector)[i - 1] + VECTOR(row_lengths_vector)[i - 1];
    }
    no_of_nodes = VECTOR(row_lengths_prefix_sum_vector)[row_count];

#define VERTEX_INDEX(i, j) (VECTOR(row_lengths_prefix_sum_vector)[j] + i - VECTOR(row_start_vector)[j])
#define ADD_EDGE_IJ_KL_IF_EXISTS(i, j, k, l)                             \
    if (is_vertex_in_graph(k, l, condition_params))                      \
    {                                                                    \
        igraph_vector_int_push_back(&edges, VERTEX_INDEX((i), (j)));     \
        igraph_vector_int_push_back(&edges, VERTEX_INDEX((k), (l)));     \
        if (directed && mutual)                                          \
        {                                                                \
            igraph_vector_int_push_back(&edges, VERTEX_INDEX((k), (l))); \
            igraph_vector_int_push_back(&edges, VERTEX_INDEX((i), (j))); \
        }                                                                \
    }

    igraph_integer_t k;
    for (j = 0; j < row_count; j++)
    {
        for (i = 0; i < VECTOR(row_lengths_vector)[j]; i++)
        {
            k = VECTOR(row_start_vector)[j] + i;
            ADD_EDGE_IJ_KL_IF_EXISTS(k, j, (k + 1), j);
            ADD_EDGE_IJ_KL_IF_EXISTS(k, j, k, (j + 1));
            ADD_EDGE_IJ_KL_IF_EXISTS(k, j, (k - 1), (j + 1));
        }
    }

    IGRAPH_CHECK(igraph_create(graph, &edges, no_of_nodes, directed));
    igraph_vector_int_destroy(&row_lengths_prefix_sum_vector);

    return IGRAPH_SUCCESS;
}
