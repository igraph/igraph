/*
   IGraph library.
   Copyright (C) 2022  The igraph development team <igraph@igraph.org>

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

#include "igraph_constructors.h"

#include "core/interruption.h"
#include "math/safe_intop.h"

/**
 * \ingroup generators
 * \function igraph_triangle_lattice
 * \brief A triangle lattice with the given shape.
 *
 * Creates a triangle lattice whose vertices have the form (i, j) for non-negative integers i and j
 * and (i, j) is connected with (i + 1, j), (i, j + 1), and (i - 1, j + 1) provided a vertex
 * exists. Thus, all vertices have degree at most 6.
 *
 * </para><para>
 * The vertices of the resulting graph are ordered lexicographically with the 2nd coordinate being
 * more significant, e.g., (i, j) < (i + 1, j) and (i + 1, j) < (i, j + 1).
 *
 * \param graph An uninitialized graph object.
 * \param directed Boolean, whether to create a directed graph.
 *        If the \c mutual argument is not set to true,
 *        edges will be directed from lower-index vertices towards
 *        higher-index ones.
 * \param mutual Boolean, if the graph is directed this gives whether
 *        to create all connections as mutual.
 * \param row_lengths_vector Integer vector, defines the number of vertices with
 *        the second coordinate equal to the index. The length of this vector must match
 *        the length of \p row_start_vector. All coordinates must be non-negative.
 * \param row_start_vector Integer vector, defines the leftmost coordinate of
 *        the vertex with the second coordinate equal to the index.
 *
 * \return Error code:
 *         \c IGRAPH_EINVAL: invalid (negative) length of row_lengths_vector does not match the length of the
            row_start_vector.
 *
 * Time complexity:  O(|V|), where |V| is the number of vertices in the generated graph.
 */
static igraph_error_t triangle_lattice(
    igraph_t *graph, igraph_bool_t directed, igraph_bool_t mutual,
    const igraph_vector_int_t *row_lengths_vector, const igraph_vector_int_t *row_start_vector
) {
    igraph_vector_int_t edges = IGRAPH_VECTOR_NULL;
    igraph_integer_t row_count = igraph_vector_int_size(row_lengths_vector);
    igraph_integer_t no_of_nodes;
    igraph_vector_int_t row_lengths_prefix_sum_vector;
    igraph_integer_t i, j;

    if (igraph_vector_int_size(row_lengths_vector) != igraph_vector_int_size(row_start_vector)) {
        IGRAPH_ERRORF(
            "Length of row_lengths_vector vector (%" IGRAPH_PRId ") must match the length of the "
            "row_start_vector (%" IGRAPH_PRId ").",
            IGRAPH_EINVAL,
            igraph_vector_int_size(row_lengths_vector),
            igraph_vector_int_size(row_start_vector)
        );
    }

    for (i = 0; i < row_count; i++) {
        if (VECTOR(*row_lengths_vector)[i] < 0) {
            IGRAPH_ERRORF(
                "row_lengths_vector vector must have non-negative coordinates, "
                "was (%" IGRAPH_PRId ") for the (%" IGRAPH_PRId ")-th row.",
                IGRAPH_EINVAL, VECTOR(*row_lengths_vector)[i], i);
        }
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);

    IGRAPH_VECTOR_INT_INIT_FINALLY(&row_lengths_prefix_sum_vector, row_count + 1);
    VECTOR(row_lengths_prefix_sum_vector)[0] = 0;
    for (i = 1; i < row_count + 1; i++)
    {
        IGRAPH_SAFE_ADD(VECTOR(row_lengths_prefix_sum_vector)[i - 1], VECTOR(*row_lengths_vector)[i - 1], &(VECTOR(row_lengths_prefix_sum_vector)[i]));
    }
    no_of_nodes = VECTOR(row_lengths_prefix_sum_vector)[row_count];

#define VERTEX_INDEX(i, j) (VECTOR(row_lengths_prefix_sum_vector)[j] + i - VECTOR(*row_start_vector)[j])
#define ROW_END(j) (VECTOR(*row_start_vector)[j] + VECTOR(*row_lengths_vector)[j] - 1)
#define ADD_EDGE_IJ_KL_IF_EXISTS(i, j, k, l)                                                                                                      \
    if (VECTOR(*row_start_vector)[l] <= k && k <= ROW_END(l) && 0 <= l && l <= row_count - 1)                                                     \
    {                                                                                                                                             \
        igraph_vector_int_push_back(&edges, VERTEX_INDEX((i), (j)));      /* reserved */                                                          \
        igraph_vector_int_push_back(&edges, VERTEX_INDEX((k), (l)));      /* reserved */                                                          \
        if (directed && mutual)                                                                                                                   \
        {                                                                                                                                         \
            igraph_vector_int_push_back(&edges, VERTEX_INDEX((k), (l)));    /* reserved */                                                        \
            igraph_vector_int_push_back(&edges, VERTEX_INDEX((i), (j)));    /* reserved */                                                        \
        }                                                                                                                                         \
    }

     /* computing the number of edges in the constructed triangle lattice */
    igraph_integer_t no_of_edges2 = VECTOR(*row_lengths_vector)[row_count - 1] - 1;
    igraph_integer_t multiplier = mutual && directed ? 4 : 2;
    for (j = 0; j < row_count - 1; j++) {
        IGRAPH_SAFE_ADD(no_of_edges2, VECTOR(*row_lengths_vector)[j], &no_of_edges2);
        IGRAPH_SAFE_ADD(no_of_edges2, fmin(ROW_END(j), ROW_END((j + 1))) - fmax(VECTOR(*row_start_vector)[j], VECTOR(*row_start_vector)[j + 1]) + 1,
                        &no_of_edges2);
        IGRAPH_SAFE_ADD(no_of_edges2, fmin(ROW_END(j), ROW_END((j + 1)) + 1) - fmax(VECTOR(*row_start_vector)[j], VECTOR(*row_start_vector)[j + 1] + 1) + 1,
                        &no_of_edges2);
    }
    IGRAPH_SAFE_MULT(no_of_edges2, multiplier, &no_of_edges2);
    IGRAPH_CHECK(igraph_vector_int_reserve(&edges, no_of_edges2));

    /* constructing the edge array */
    igraph_integer_t k;
    for (j = 0; j < row_count; j++) {
        IGRAPH_ALLOW_INTERRUPTION();
        for (i = 0; i < VECTOR(*row_lengths_vector)[j]; i++) {
            k = VECTOR(*row_start_vector)[j] + i;
            ADD_EDGE_IJ_KL_IF_EXISTS(k, j, (k + 1), j);
            if (j < row_count - 1) {
                ADD_EDGE_IJ_KL_IF_EXISTS(k, j, k, (j + 1));
                ADD_EDGE_IJ_KL_IF_EXISTS(k, j, (k - 1), (j + 1));
            }
        }
    }

    IGRAPH_CHECK(igraph_create(graph, &edges, no_of_nodes, directed));

    igraph_vector_int_destroy(&row_lengths_prefix_sum_vector);
    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}

static igraph_error_t triangle_lattice_triangle_shape(igraph_t *graph, igraph_integer_t size, igraph_bool_t directed, igraph_bool_t mutual) {
    igraph_integer_t row_count = size;
    igraph_vector_int_t row_lengths_vector;
    igraph_vector_int_t row_start_vector;
    igraph_integer_t i;

    IGRAPH_VECTOR_INT_INIT_FINALLY(&row_lengths_vector, row_count);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&row_start_vector, row_count);

    for (i = 0; i < row_count; i++) {
        VECTOR(row_lengths_vector)[i] = size - i;
        VECTOR(row_start_vector)[i] = 0;
    }

    IGRAPH_CHECK(triangle_lattice(graph, directed, mutual, &row_lengths_vector, &row_start_vector));

    igraph_vector_int_destroy(&row_lengths_vector);
    igraph_vector_int_destroy(&row_start_vector);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}

static igraph_error_t triangle_lattice_rectangle_shape(
    igraph_t *graph, igraph_integer_t size_x, igraph_integer_t size_y,
    igraph_bool_t directed, igraph_bool_t mutual
) {
    igraph_integer_t row_count = size_y;
    igraph_vector_int_t row_lengths_vector;
    igraph_vector_int_t row_start_vector;
    igraph_integer_t i;

    IGRAPH_VECTOR_INT_INIT_FINALLY(&row_lengths_vector, row_count);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&row_start_vector, row_count);

    for (i = 0; i < row_count; i++) {
        VECTOR(row_lengths_vector)[i] = size_x;
        VECTOR(row_start_vector)[i] = (row_count - i) / 2;
    }

    IGRAPH_CHECK(triangle_lattice(graph, directed, mutual, &row_lengths_vector, &row_start_vector));

    igraph_vector_int_destroy(&row_lengths_vector);
    igraph_vector_int_destroy(&row_start_vector);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}

static igraph_error_t triangle_lattice_hex_shape(
    igraph_t *graph, igraph_integer_t size_x, igraph_integer_t size_y,
    igraph_integer_t size_z, igraph_bool_t directed, igraph_bool_t mutual
) {
    igraph_integer_t row_count = size_y + size_z - 1;
    igraph_vector_int_t row_lengths_vector;
    igraph_vector_int_t row_start_vector;
    igraph_integer_t i;

    IGRAPH_VECTOR_INT_INIT_FINALLY(&row_lengths_vector, row_count);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&row_start_vector, row_count);

    igraph_integer_t row_length = size_x;
    igraph_integer_t row_start = size_y - 1;
    igraph_integer_t first_threshold = fmin(size_y - 1, size_z - 1);
    igraph_integer_t second_threshold = fmax(size_y - 1, size_z - 1);
    igraph_integer_t sgn_flag = size_y < size_z ? 0 : -1;

    for (i = 0; i < row_count; i++) {
        VECTOR(row_lengths_vector)[i] = row_length;
        VECTOR(row_start_vector)[i] = row_start;

        if (i < first_threshold) {
            row_length++;
            row_start--;
        } else if (i < second_threshold) {
            row_start += sgn_flag;
        } else {
            row_length--;
        }
    }

    IGRAPH_CHECK(triangle_lattice(graph, directed, mutual, &row_lengths_vector, &row_start_vector));

    igraph_vector_int_destroy(&row_lengths_vector);
    igraph_vector_int_destroy(&row_start_vector);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup generators
 * \function igraph_triangle_lattice
 * \brief A triangle lattice with the given shape.
 *
 * Creates a triangle lattice whose vertices have the form (i, j) for non-negative integers i and j
 * and (i, j) is generally connected with (i + 1, j), (i, j + 1), and (i - 1, j + 1).
 *
 * </para><para>
 * The vertices of the resulting graph are ordered lexicographically with the 2nd coordinate being
 * more significant, e.g., (i, j) < (i + 1, j) and (i + 1, j) < (i, j + 1)
 *
 * \param graph An uninitialized graph object.
 * \param dims Integer vector, defines the shape of the lattice. (Below the "edge length"s are in terms of graph theoretical path lengths.)
 *        If \p dims is of length 1, the resulting lattice has a triangular shape
 *        where each side of the triangle contains \c "dims[0]" vertices.
 *        If \p dims is of length 2, the resulting lattice has a
 *        "quasi rectangular" shape with the sides containing \c "dims[0]" and
 *        \c "dims[1]" vertices, respectively.
 *        If \p dims is of length 3, the resulting lattice has a hexagonal shape
 *        where the sides of the hexagon contain \c "dims[0]", \c "dims[1]" and
 *        \c "dims[2]" vertices. See https://github.com/igraph/igraph/issues/1842 for
 *        lattice visualizations.
 * \param directed Boolean, whether to create a directed graph.
 *        If the \c mutual argument is not set to true,
 *        edges will be directed from lower-index vertices towards
 *        higher-index ones.
 * \param mutual Boolean, if the graph is directed this gives whether
 *        to create all connections as mutual.
 * \return Error code:
 *         \c IGRAPH_EINVAL: The size of \p dims must be either 1, 2, or 3 with all the components
 *         at least 1.
 *
 * Time complexity:  O(|V|), where |V| is the number of vertices in the generated graph.
 *
 */
igraph_error_t igraph_triangle_lattice(
    igraph_t *graph, const igraph_vector_int_t *dims, igraph_bool_t directed,
    igraph_bool_t mutual
) {
    igraph_integer_t num_dims = igraph_vector_int_size(dims);
    if (igraph_vector_int_any_smaller(dims, 1)) {
        IGRAPH_ERROR("Invalid dimension vector.", IGRAPH_EINVAL);
    }

    switch (num_dims) {
    case 1:
        IGRAPH_CHECK(triangle_lattice_triangle_shape(graph, VECTOR(*dims)[0], directed, mutual));
        break;
    case 2:
        IGRAPH_CHECK(triangle_lattice_rectangle_shape(graph, VECTOR(*dims)[0], VECTOR(*dims)[1], directed, mutual));
        break;
    case 3:
        IGRAPH_CHECK(triangle_lattice_hex_shape(graph, VECTOR(*dims)[0], VECTOR(*dims)[1], VECTOR(*dims)[2], directed, mutual));
        break;
    default:
        IGRAPH_ERRORF(
            "The size of the dimension vector must be 1, 2 or 3, got %" IGRAPH_PRId ".",
            IGRAPH_EINVAL, num_dims
        );
    }

    return IGRAPH_SUCCESS;
}
