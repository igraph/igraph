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
#include "igraph_interface.h"

#include "core/interruption.h"
#include "math/safe_intop.h"

#define MIN(n, m) (n < m ? n : m)
#define MAX(n, m) (n < m ? m : n)

#define VERTEX_INDEX(i, j) \
    lex_ordering ? row_count * (i - VECTOR(*row_start_vector)[j]) + j : (VECTOR(row_lengths_prefix_sum_vector)[j] + i - VECTOR(*row_start_vector)[j])

#define ROW_END(j) (VECTOR(*row_start_vector)[j] + VECTOR(*row_lengths_vector)[j] - 1)
#define ADD_EDGE_IJ_KL_IF_EXISTS(i, j, k, l)                                                  \
    if (VECTOR(*row_start_vector)[l] <= k && k <= ROW_END(l) && 0 <= l && l <= row_count - 1) \
    {                                                                                         \
        igraph_vector_int_push_back(&edges, VERTEX_INDEX((i), (j))); /* reserved */           \
        igraph_vector_int_push_back(&edges, VERTEX_INDEX((k), (l))); /* reserved */           \
        if (directed && mutual)                                                               \
        {                                                                                     \
            igraph_vector_int_push_back(&edges, VERTEX_INDEX((k), (l))); /* reserved */       \
            igraph_vector_int_push_back(&edges, VERTEX_INDEX((i), (j))); /* reserved */       \
        }                                                                                     \
    }

#define COMPUTE_NUMBER_OF_VERTICES()                                                                                                                        \
    do                                                                                                                                                      \
    {                                                                                                                                                       \
        IGRAPH_VECTOR_INT_INIT_FINALLY(&row_lengths_prefix_sum_vector, row_count + 1);                                                                      \
        VECTOR(row_lengths_prefix_sum_vector)[0] = 0;                                                                                                       \
        for (i = 1; i < row_count + 1; i++)                                                                                                                 \
        {                                                                                                                                                   \
            IGRAPH_SAFE_ADD(VECTOR(row_lengths_prefix_sum_vector)[i - 1], VECTOR(*row_lengths_vector)[i - 1], &(VECTOR(row_lengths_prefix_sum_vector)[i])); \
        }                                                                                                                                                   \
        no_of_nodes = VECTOR(row_lengths_prefix_sum_vector)[row_count];                                                                                     \
    } while (0)

/**
 * Creates a triangular lattice whose vertices have the form (i, j) for non-negative integers i and j
 * and (i, j) is connected with (i + 1, j), (i, j + 1), and (i - 1, j + 1) provided a vertex
 * exists. Thus, all vertices have degree at most 6.
 *
 * </para><para>
 * The vertices of the resulting graph are ordered lexicographically with the 2nd coordinate being
 * more significant, e.g., (i, j) &lt; (i + 1, j) and (i + 1, j) &lt; (i, j + 1) unless
 * \c lex_ordering is set to true in which case the roles of the coordinates are reversed.
 *
 * \param graph An uninitialized graph object.
 * \param directed Boolean, whether to create a directed graph.
 *        If the \c mutual argument is not set to true,
 *        edges will be directed from lower-index vertices towards
 *        higher-index ones.
 * \param mutual Boolean, if the graph is directed this gives whether
 *        to create all connections as mutual.
 * \param lex_ordering Boolean, set to true if the vertices of the resulting graph are ordered
 *        lexicographically with the 1st coordinate being more significant. Use only when all the
 *        rows have the number of vertices.
 * \param row_lengths_vector Integer vector, defines the number of vertices with
 *        the second coordinate equal to the index. The length of this vector must match
 *        the length of \p row_start_vector. All coordinates must be non-negative.
 * \param row_start_vector Integer vector, defines the leftmost coordinate of
 *        the vertex with the second coordinate equal to the index.
 *
 * \return Error code:
 *         \c IGRAPH_EINVAL: invalid (negative) length of row_lengths_vector does not match the length of the
 *         row_start_vector.
 *
 * Time complexity:  O(|V|), where |V| is the number of vertices in the generated graph.
 */
static igraph_error_t triangular_lattice(
    igraph_t *graph, igraph_bool_t directed, igraph_bool_t mutual, igraph_bool_t lex_ordering,
    const igraph_vector_int_t *row_lengths_vector, const igraph_vector_int_t *row_start_vector) {
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
            igraph_vector_int_size(row_start_vector));
    }


    if (row_count > 0 && lex_ordering && !igraph_vector_int_isininterval(row_lengths_vector, VECTOR(*row_lengths_vector)[0], VECTOR(*row_lengths_vector)[0])) {
        IGRAPH_ERROR(
            "row_lengths_vector must have all the coordinates the same", IGRAPH_EINVAL);
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

    COMPUTE_NUMBER_OF_VERTICES();

    /* computing the number of edges in the constructed triangular lattice */
    igraph_integer_t no_of_edges2 = VECTOR(*row_lengths_vector)[row_count - 1] - 1;
    igraph_integer_t multiplier = mutual && directed ? 4 : 2;
    for (j = 0; j < row_count - 1; j++) {
        IGRAPH_SAFE_ADD(no_of_edges2, VECTOR(*row_lengths_vector)[j] - 1, &no_of_edges2);
        IGRAPH_SAFE_ADD(no_of_edges2, MIN(ROW_END(j), ROW_END((j + 1))) - MAX(VECTOR(*row_start_vector)[j], VECTOR(*row_start_vector)[j + 1]) + 1,
                        &no_of_edges2);
        IGRAPH_SAFE_ADD(no_of_edges2, MIN(ROW_END(j), ROW_END((j + 1)) + 1) - MAX(VECTOR(*row_start_vector)[j], VECTOR(*row_start_vector)[j + 1] + 1) + 1,
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

static igraph_error_t triangular_lattice_triangle_shape(igraph_t *graph, igraph_integer_t size, igraph_bool_t directed, igraph_bool_t mutual) {
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

    IGRAPH_CHECK(triangular_lattice(graph, directed, mutual, false, &row_lengths_vector, &row_start_vector));

    igraph_vector_int_destroy(&row_lengths_vector);
    igraph_vector_int_destroy(&row_start_vector);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}

static igraph_error_t triangular_lattice_rectangle_shape(
    igraph_t *graph, igraph_integer_t size_x, igraph_integer_t size_y,
    igraph_bool_t directed, igraph_bool_t mutual) {
    igraph_integer_t row_count = size_x;
    igraph_vector_int_t row_lengths_vector;
    igraph_vector_int_t row_start_vector;
    igraph_integer_t i;

    IGRAPH_VECTOR_INT_INIT_FINALLY(&row_lengths_vector, row_count);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&row_start_vector, row_count);

    for (i = 0; i < row_count; i++) {
        VECTOR(row_lengths_vector)[i] = size_y;
        VECTOR(row_start_vector)[i] = (row_count - i) / 2;
    }

    IGRAPH_CHECK(triangular_lattice(graph, directed, mutual, false, &row_lengths_vector, &row_start_vector));

    igraph_vector_int_destroy(&row_lengths_vector);
    igraph_vector_int_destroy(&row_start_vector);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}

static igraph_error_t triangular_lattice_hex_shape(
    igraph_t *graph, igraph_integer_t size_x, igraph_integer_t size_y,
    igraph_integer_t size_z, igraph_bool_t directed, igraph_bool_t mutual) {
    igraph_integer_t row_count = size_y + size_z - 1;
    igraph_vector_int_t row_lengths_vector;
    igraph_vector_int_t row_start_vector;
    igraph_integer_t i;

    IGRAPH_VECTOR_INT_INIT_FINALLY(&row_lengths_vector, row_count);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&row_start_vector, row_count);

    igraph_integer_t row_length = size_x;
    igraph_integer_t row_start = size_y - 1;
    igraph_integer_t first_threshold = MIN(size_y - 1, size_z - 1);
    igraph_integer_t second_threshold = MAX(size_y - 1, size_z - 1);
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

    IGRAPH_CHECK(triangular_lattice(graph, directed, mutual, false, &row_lengths_vector, &row_start_vector));

    igraph_vector_int_destroy(&row_lengths_vector);
    igraph_vector_int_destroy(&row_start_vector);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_triangular_lattice
 * \brief A triangular lattice with the given shape.
 *
 * \experimental
 *
 * Creates a triangular lattice whose vertices have the form (i, j) for non-negative integers i and j
 * and (i, j) is generally connected with (i + 1, j), (i, j + 1), and (i - 1, j + 1).
 * The function constructs a planar dual of the graph constructed by \ref igraph_hexagonal_lattice().
 * In particular, there a one-to-one correspondence between the vertices in the constructed graph
 * and the cycles of length 6 in the graph constructed by \ref igraph_hexagonal_lattice()
 * with the same \p dims parameter.
 *
 * </para><para>
 * The vertices of the resulting graph are ordered lexicographically with the 2nd coordinate being
 * more significant, e.g., (i, j) &lt; (i + 1, j) and (i + 1, j) &lt; (i, j + 1)
 *
 * \param graph An uninitialized graph object.
 * \param dims Integer vector, defines the shape of the lattice. (Below the "edge length"s are in terms of graph theoretical path lengths.)
 *        If \p dims is of length 1, the resulting lattice has a triangular shape
 *        where each side of the triangle contains <code>dims[0]</code> vertices.
 *        If \p dims is of length 2, the resulting lattice has a
 *        "quasi rectangular" shape with the sides containing <code>dims[0]</code> and
 *        <code>dims[1]</code> vertices, respectively.
 *        If \p dims is of length 3, the resulting lattice has a hexagonal shape
 *        where the sides of the hexagon contain <code>dims[0]</code>, <code>dims[1]</code> and
 *        <code>dims[2]</code> vertices.
 *        All coordinates must be non-negative.
 * \param directed Boolean, whether to create a directed graph.
 *        If the \c mutual argument is not set to true,
 *        edges will be directed from lower-index vertices towards
 *        higher-index ones.
 * \param mutual Boolean, if the graph is directed this gives whether
 *        to create all connections as mutual.
 * \return Error code:
 *         \c IGRAPH_EINVAL: The size of \p dims must be either 1, 2, or 3 with all the components
 *         at least 1.
 * \sa \ref igraph_hexagonal_lattice() for creating a triangular lattice.
 *
 * Time complexity:  O(|V|), where |V| is the number of vertices in the generated graph.
 *
 */
igraph_error_t igraph_triangular_lattice(
    igraph_t *graph, const igraph_vector_int_t *dims, igraph_bool_t directed,
    igraph_bool_t mutual) {
    igraph_integer_t num_dims = igraph_vector_int_size(dims);
    if (igraph_vector_int_any_smaller(dims, 0)) {
        IGRAPH_ERROR("Invalid dimension vector.", IGRAPH_EINVAL);
    }
    /* If a coordinate of dims is 0 the result is an empty graph. */
    if (igraph_vector_int_contains(dims, 0)) {
        return igraph_empty(graph, 0, directed);
    }

    switch (num_dims) {
    case 1:
        IGRAPH_CHECK(triangular_lattice_triangle_shape(graph, VECTOR(*dims)[0], directed, mutual));
        break;
    case 2:
        IGRAPH_CHECK(triangular_lattice_rectangle_shape(graph, VECTOR(*dims)[0], VECTOR(*dims)[1], directed, mutual));
        break;
    case 3:
        IGRAPH_CHECK(triangular_lattice_hex_shape(graph, VECTOR(*dims)[0], VECTOR(*dims)[1], VECTOR(*dims)[2], directed, mutual));
        break;
    default:
        IGRAPH_ERRORF(
            "The size of the dimension vector must be 1, 2 or 3, got %" IGRAPH_PRId ".",
            IGRAPH_EINVAL, num_dims);
    }

    return IGRAPH_SUCCESS;
}


/**
 * Creates a hexagonal lattice whose vertices have the form (i, j) for non-negative integers i and j
 * and (i, j) is connected with (i + 1, j), and if i is odd also with (i - 1, j + 1) provided a vertex
 * exists. Thus, all vertices have degree at most 3.
 *
 * </para><para>
 * The vertices of the resulting graph are ordered lexicographically with the 2nd coordinate being
 * more significant, e.g., (i, j) &lt; (i + 1, j) and (i + 1, j) &lt; (i, j + 1).
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
 *         row_start_vector.
 *
 * Time complexity:  O(|V|), where |V| is the number of vertices in the generated graph.
 */
static igraph_error_t hexagonal_lattice(
    igraph_t *graph, igraph_bool_t directed, igraph_bool_t mutual,
    const igraph_vector_int_t *row_lengths_vector, const igraph_vector_int_t *row_start_vector
) {
    igraph_vector_int_t edges = IGRAPH_VECTOR_NULL;
    igraph_integer_t row_count = igraph_vector_int_size(row_lengths_vector);
    igraph_integer_t no_of_nodes;
    igraph_vector_int_t row_lengths_prefix_sum_vector;
    igraph_integer_t i, j;
    igraph_bool_t lex_ordering = false;

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

    COMPUTE_NUMBER_OF_VERTICES();

    /* computing the number of edges in the constructed hex lattice */
    igraph_integer_t no_of_edges2 = VECTOR(*row_lengths_vector)[row_count - 1] - 1;
    igraph_integer_t multiplier = mutual && directed ? 4 : 2, low, high;
    for (j = 0; j < row_count - 1; j++) {
        IGRAPH_SAFE_ADD(no_of_edges2, VECTOR(*row_lengths_vector)[j] - 1, &no_of_edges2);
        low = MAX((VECTOR(*row_start_vector)[j] - 1), (VECTOR(*row_start_vector)[j + 1]));
        low = low % 2 ? low + 1 : low;
        high = MIN((ROW_END(j) - 1), (ROW_END(j + 1)));
        high = high % 2 ? high - 1 : high;
        IGRAPH_SAFE_ADD(no_of_edges2, (high - low) / 2 + 1, &no_of_edges2);
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
            if (j < row_count - 1 && k % 2 == 1) {
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

static igraph_error_t hexagonal_lattice_triangle_shape(igraph_t *graph, igraph_integer_t size, igraph_bool_t directed, igraph_bool_t mutual) {
    igraph_integer_t row_count;
    IGRAPH_SAFE_ADD(size, 2, &row_count);
    igraph_vector_int_t row_lengths_vector;
    igraph_vector_int_t row_start_vector;
    igraph_integer_t i;

    IGRAPH_VECTOR_INT_INIT_FINALLY(&row_lengths_vector, row_count - 1);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&row_start_vector, row_count - 1);

    for (i = 0; i < row_count - 1; i++) {
        VECTOR(row_lengths_vector)[i] = 2 * (row_count - i) - (i ? 1 : 3);
        VECTOR(row_start_vector)[i] = (i ? 0 : 1);
    }

    IGRAPH_CHECK(hexagonal_lattice(graph, directed, mutual, &row_lengths_vector, &row_start_vector));

    igraph_vector_int_destroy(&row_lengths_vector);
    igraph_vector_int_destroy(&row_start_vector);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}

static igraph_error_t hexagonal_lattice_rectangle_shape(
    igraph_t *graph, igraph_integer_t size_x, igraph_integer_t size_y,
    igraph_bool_t directed, igraph_bool_t mutual
) {
    igraph_integer_t row_count;
    IGRAPH_SAFE_ADD(size_x, 1, &row_count);
    igraph_vector_int_t row_lengths_vector;
    igraph_vector_int_t row_start_vector;
    igraph_integer_t actual_size_y;
    IGRAPH_SAFE_ADD(size_y, 1, &actual_size_y);
    IGRAPH_SAFE_MULT(actual_size_y, 2, &actual_size_y);
    igraph_integer_t i;

    IGRAPH_VECTOR_INT_INIT_FINALLY(&row_lengths_vector, row_count);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&row_start_vector, row_count);

    igraph_bool_t is_first_row, is_last_row, is_start_odd;

    for (i = 0; i < row_count; i++) {
        is_first_row = (i == 0);
        is_last_row = i == row_count - 1;
        is_start_odd = (row_count - i - 1) % 2;
        VECTOR(row_lengths_vector)[i] = actual_size_y - (is_first_row || is_last_row ? 1 : 0);
        VECTOR(row_start_vector)[i] = row_count - i - 1 + (is_first_row && !is_start_odd ? 1 : 0);
    }

    IGRAPH_CHECK(hexagonal_lattice(graph, directed, mutual, &row_lengths_vector, &row_start_vector));

    igraph_vector_int_destroy(&row_lengths_vector);
    igraph_vector_int_destroy(&row_start_vector);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}

static igraph_error_t hexagonal_lattice_hex_shape(
    igraph_t *graph, igraph_integer_t size_x, igraph_integer_t size_y,
    igraph_integer_t size_z, igraph_bool_t directed, igraph_bool_t mutual
) {
    igraph_integer_t row_count = size_y + size_z;
    igraph_vector_int_t row_lengths_vector;
    igraph_vector_int_t row_start_vector;
    igraph_integer_t i;

    IGRAPH_VECTOR_INT_INIT_FINALLY(&row_lengths_vector, row_count);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&row_start_vector, row_count);

    igraph_integer_t row_length;
    IGRAPH_SAFE_MULT(size_x, 2, &row_length);
    IGRAPH_SAFE_ADD(row_length, 1, &row_length);
    igraph_integer_t row_start;
    IGRAPH_SAFE_MULT(size_y, 2, &row_start);
    IGRAPH_SAFE_ADD(row_start, -1, &row_start);
    igraph_integer_t first_threshold = MIN(size_y - 1, size_z - 1);
    igraph_integer_t second_threshold = MAX(size_y - 1, size_z - 1);
    igraph_integer_t sgn_flag = size_y < size_z ? 0 : -2;

    for (i = 0; i < row_count; i++) {
        VECTOR(row_lengths_vector)[i] = row_length;
        VECTOR(row_start_vector)[i] = row_start;

        if (i < first_threshold) {
            row_length += 2;
            row_start -= 2;
        } else if (i < second_threshold) {
            row_start += sgn_flag;
        } else {
            row_length -= 2;
        }
        if (i == size_y - 1) {
            row_start--;
            row_length++;
        }
        if (i == size_z - 1) {
            row_length++;
        }
    }

    IGRAPH_CHECK(hexagonal_lattice(graph, directed, mutual, &row_lengths_vector, &row_start_vector));

    igraph_vector_int_destroy(&row_lengths_vector);
    igraph_vector_int_destroy(&row_start_vector);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_hexagonal_lattice
 * \brief A hexagonal lattice with the given shape.
 *
 * \experimental
 *
 * Creates a hexagonal lattice whose vertices have the form (i, j) for non-negative integers i and j
 * and (i, j) is generally connected with (i + 1, j), and if i is odd also with (i - 1, j + 1).
 * The function constructs a planar dual of the graph constructed by \ref igraph_triangular_lattice().
 * In particular, there a one-to-one correspondence between the cycles of length 6 in the constructed graph
 * and the vertices of the graph constructed by \ref igraph_triangular_lattice() function
 * with the same \p dims parameter.
 *
 * </para><para>
 * The vertices of the resulting graph are ordered lexicographically with the 2nd coordinate being
 * more significant, e.g., (i, j) &lt; (i + 1, j) and (i + 1, j) &lt; (i, j + 1)
 *
 * \param graph An uninitialized graph object.
 * \param dims Integer vector, defines the shape of the lattice. (Below the "edge length"s are in terms of graph theoretical path lengths.)
 *        If \p dims is of length 1, the resulting lattice has a triangular shape
 *        where each side of the triangle contains <code>dims[0]</code> vertices.
 *        If \p dims is of length 2, the resulting lattice has a
 *        "quasi rectangular" shape with the sides containing <code>dims[0]</code> and
 *        <code>dims[1]</code> vertices, respectively.
 *        If \p dims is of length 3, the resulting lattice has a hexagonal shape
 *        where the sides of the hexagon contain <code>dims[0]</code>, <code>dims[1]</code> and
 *        <code>dims[2]</code> vertices.
 *        All coordinates must be non-negative.
 * \param directed Boolean, whether to create a directed graph.
 *        If the \c mutual argument is not set to true,
 *        edges will be directed from lower-index vertices towards
 *        higher-index ones.
 * \param mutual Boolean, if the graph is directed this gives whether
 *        to create all connections as mutual.
 * \return Error code:
 *         \c IGRAPH_EINVAL: The size of \p dims must be either 1, 2, or 3 with all the components
 *         at least 1.
 * \sa \ref igraph_triangular_lattice() for creating a triangular lattice.
 *
 * Time complexity:  O(|V|), where |V| is the number of vertices in the generated graph.
 *
 */
igraph_error_t igraph_hexagonal_lattice(
    igraph_t *graph, const igraph_vector_int_t *dims, igraph_bool_t directed,
    igraph_bool_t mutual
) {
    igraph_integer_t num_dims = igraph_vector_int_size(dims);
    if (igraph_vector_int_any_smaller(dims, 0)) {
        IGRAPH_ERROR("Invalid dimension vector.", IGRAPH_EINVAL);
    }
    /* If a coordinate of dims is 0 the result is an empty graph. */
    if (igraph_vector_int_any_smaller(dims, 1)) {
        return igraph_empty(graph, 0, directed);
    }

    switch (num_dims) {
    case 1:
        IGRAPH_CHECK(hexagonal_lattice_triangle_shape(graph, VECTOR(*dims)[0], directed, mutual));
        break;
    case 2:
        IGRAPH_CHECK(hexagonal_lattice_rectangle_shape(graph, VECTOR(*dims)[0], VECTOR(*dims)[1], directed, mutual));
        break;
    case 3:
        IGRAPH_CHECK(hexagonal_lattice_hex_shape(graph, VECTOR(*dims)[0], VECTOR(*dims)[1], VECTOR(*dims)[2], directed, mutual));
        break;
    default:
        IGRAPH_ERRORF(
            "The size of the dimension vector must be 1, 2 or 3, got %" IGRAPH_PRId ".",
            IGRAPH_EINVAL, num_dims
        );
    }
    return IGRAPH_SUCCESS;
}
