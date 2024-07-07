/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2005-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge, MA 02139 USA

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

#include "igraph_interface.h"
#include "igraph_nongraph.h"
#include "igraph_paths.h"

#include "core/interruption.h"

/**
 * \ingroup nongraph
 * \function igraph_running_mean
 * \brief Calculates the running mean of a vector.
 *
 * </para><para>
 * The running mean is defined by the mean of the
 * previous \p binwidth values.
 * \param data The vector containing the data.
 * \param res The vector containing the result. This should be
 *        initialized before calling this function and will be
 *        resized.
 * \param binwidth Integer giving the width of the bin for the running
 *        mean calculation.
 * \return Error code.
 *
 * Time complexity: O(n),
 * n is the length of
 * the data vector.
 */

igraph_error_t igraph_running_mean(const igraph_vector_t *data, igraph_vector_t *res,
                        igraph_integer_t binwidth) {

    double sum = 0;
    igraph_integer_t i;

    /* Check */
    if (igraph_vector_size(data) < binwidth) {
        IGRAPH_ERRORF("Data vector length (%" IGRAPH_PRId ") smaller than bin width (%" IGRAPH_PRId ").", IGRAPH_EINVAL, igraph_vector_size(data), binwidth);
    }
    if (binwidth < 1) {
        IGRAPH_ERRORF("Bin width for running mean should be at least 1, got %" IGRAPH_PRId ".", IGRAPH_EINVAL, binwidth);
    }

    /* Memory for result */

    IGRAPH_CHECK(igraph_vector_resize(res, (igraph_vector_size(data) - binwidth + 1)));

    /* Initial bin */
    for (i = 0; i < binwidth; i++) {
        sum += VECTOR(*data)[i];
    }

    VECTOR(*res)[0] = sum / binwidth;

    for (i = 1; i < igraph_vector_size(data) - binwidth + 1; i++) {
        IGRAPH_ALLOW_INTERRUPTION();
        sum -= VECTOR(*data)[i - 1];
        sum += VECTOR(*data)[ (i + binwidth - 1)];
        VECTOR(*res)[i] = sum / binwidth;
    }

    return IGRAPH_SUCCESS;
}


/**
 * \ingroup nongraph
 * \function igraph_convex_hull
 * \brief Determines the convex hull of a given set of points in the 2D plane.
 *
 * </para><para>
 * The convex hull is determined by the Graham scan algorithm.
 * See the following reference for details:
 *
 * </para><para>
 * Thomas H. Cormen, Charles E. Leiserson, Ronald L. Rivest, and Clifford
 * Stein. Introduction to Algorithms, Second Edition. MIT Press and
 * McGraw-Hill, 2001. ISBN 0262032937. Pages 949-955 of section 33.3:
 * Finding the convex hull.
 *
 * \param data vector containing the coordinates. The length of the
 *        vector must be even, since it contains X-Y coordinate pairs.
 * \param resverts the vector containing the result, e.g. the vector of
 *        vertex indices used as the corners of the convex hull. Supply
 *        \c NULL here if you are only interested in the coordinates of
 *        the convex hull corners.
 * \param rescoords the matrix containing the coordinates of the selected
 *        corner vertices. Supply \c NULL here if you are only interested in
 *        the vertex indices.
 * \return Error code:
 *         \c IGRAPH_ENOMEM: not enough memory
 *
 * Time complexity: O(n log(n)) where n is the number of vertices.
 */
igraph_error_t igraph_convex_hull(
    const igraph_matrix_t *data, igraph_vector_int_t *resverts,
    igraph_matrix_t *rescoords
) {
    igraph_integer_t no_of_nodes;
    igraph_integer_t i, pivot_idx = 0, last_idx, before_last_idx, next_idx, j;
    igraph_vector_t angles;
    igraph_vector_int_t order, stack;
    igraph_real_t px, py, cp;

    no_of_nodes = igraph_matrix_nrow(data);
    if (igraph_matrix_ncol(data) != 2) {
        IGRAPH_ERROR("Only two-dimensional point sets are supports, matrix must have two columns.", IGRAPH_EINVAL);
    }
    if (no_of_nodes == 0) {
        if (resverts) {
            igraph_vector_int_clear(resverts);
        }
        if (rescoords) {
            IGRAPH_CHECK(igraph_matrix_resize(rescoords, 0, 2));
        }
        /**************************** this is an exit here *********/
        return IGRAPH_SUCCESS;
    }

    IGRAPH_VECTOR_INIT_FINALLY(&angles, no_of_nodes);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&stack, 0);

    /* Search for the pivot vertex */
    for (i = 1; i < no_of_nodes; i++) {
        if (MATRIX(*data, i, 1) < MATRIX(*data, pivot_idx, 1)) {
            pivot_idx = i;
        } else if (MATRIX(*data, i, 1) == MATRIX(*data, pivot_idx, 1) &&
                   MATRIX(*data, i, 0) < MATRIX(*data, pivot_idx, 0)) {
            pivot_idx = i;
        }
    }
    px = MATRIX(*data, pivot_idx, 0);
    py = MATRIX(*data, pivot_idx, 1);

    /* Create angle array */
    for (i = 0; i < no_of_nodes; i++) {
        if (i == pivot_idx) {
            /* We can't calculate the angle of the pivot point with itself,
             * so we use 10 here. This way, after sorting the angle vector,
             * the pivot point will always be the first one, since the range
             * of atan2 is -3.14..3.14 */
            VECTOR(angles)[i] = 10;
        } else {
            VECTOR(angles)[i] = atan2(MATRIX(*data, i, 1) - py, MATRIX(*data, i, 0) - px);
        }
    }

    /* Sort points by angles */
    IGRAPH_VECTOR_INT_INIT_FINALLY(&order, no_of_nodes);
    IGRAPH_CHECK(igraph_vector_qsort_ind(&angles, &order, IGRAPH_ASCENDING));

    /* Check if two points have the same angle. If so, keep only the point that
     * is farthest from the pivot */
    j = 0;
    last_idx = VECTOR(order)[0];
    pivot_idx = VECTOR(order)[no_of_nodes - 1];
    for (i = 1; i < no_of_nodes; i++) {
        next_idx = VECTOR(order)[i];
        if (VECTOR(angles)[last_idx] == VECTOR(angles)[next_idx]) {
            /* Keep the vertex that is farther from the pivot, drop the one that is
             * closer */
            px = pow(MATRIX(*data, last_idx, 0) - MATRIX(*data, pivot_idx, 0), 2) +
                 pow(MATRIX(*data, last_idx, 1) - MATRIX(*data, pivot_idx, 1), 2);
            py = pow(MATRIX(*data, next_idx, 0) - MATRIX(*data, pivot_idx, 0), 2) +
                 pow(MATRIX(*data, next_idx, 1) - MATRIX(*data, pivot_idx, 1), 2);
            if (px > py) {
                VECTOR(order)[i] = -1;
            } else {
                VECTOR(order)[j] = -1;
                last_idx = next_idx;
                j = i;
            }
        } else {
            last_idx = next_idx;
            j = i;
        }
    }

    j = 0;
    last_idx = -1;
    before_last_idx = -1;
    while (!igraph_vector_int_empty(&order)) {
        next_idx = igraph_vector_int_tail(&order);
        if (next_idx < 0) {
            /* This vertex should be skipped; was excluded in an earlier step */
            igraph_vector_int_pop_back(&order);
            continue;
        }
        /* Determine whether we are at a left or right turn */
        if (j < 2) {
            /* Pretend that we are turning into the right direction if we have less
             * than two items in the stack */
            cp = -1;
        } else {
            cp = (MATRIX(*data, last_idx, 0) - MATRIX(*data, before_last_idx, 0)) *
                 (MATRIX(*data, next_idx, 1) - MATRIX(*data, before_last_idx, 1)) -
                 (MATRIX(*data, next_idx, 0) - MATRIX(*data, before_last_idx, 0)) *
                 (MATRIX(*data, last_idx, 1) - MATRIX(*data, before_last_idx, 1));
        }

        if (cp < 0) {
            /* We are turning into the right direction */
            igraph_vector_int_pop_back(&order);
            IGRAPH_CHECK(igraph_vector_int_push_back(&stack, next_idx));
            before_last_idx = last_idx;
            last_idx = next_idx;
            j++;
        } else {
            /* No, skip back and try again in the next iteration */
            igraph_vector_int_pop_back(&stack);
            j--;
            last_idx = before_last_idx;
            before_last_idx = (j >= 2) ? VECTOR(stack)[j - 2] : -1;
        }
    }

    /* Create result vector */
    if (resverts != 0) {
        igraph_vector_int_clear(resverts);
        IGRAPH_CHECK(igraph_vector_int_append(resverts, &stack));
    }
    if (rescoords != 0) {
        igraph_matrix_select_rows(data, rescoords, &stack);
    }

    /* Free everything */
    igraph_vector_int_destroy(&order);
    igraph_vector_int_destroy(&stack);
    igraph_vector_destroy(&angles);
    IGRAPH_FINALLY_CLEAN(3);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_expand_path_to_pairs
 * \brief Helper function to convert a sequence of vertex IDs describing a path into a "pairs" vector.
 *
 * </para><para>
 * This function is useful when you have a sequence of vertex IDs in a graph and
 * you would like to retrieve the IDs of the edges between them. The function
 * duplicates all but the first and the last elements in the vector, effectively
 * converting the path into a vector of vertex IDs that can be passed to
 * \ref igraph_get_eids().
 *
 * \param  path  the input vector. It will be modified in-place and it will be
 *         resized as needed. When the vector contains less than two vertex IDs,
 *         it will be cleared.
 * \return Error code: \c IGRAPH_ENOMEM if there is not enough memory to expand
 *         the vector.
 */
igraph_error_t igraph_expand_path_to_pairs(igraph_vector_int_t* path) {
    igraph_integer_t no_of_vertices = igraph_vector_int_size(path);
    igraph_integer_t i, j, no_of_items = (no_of_vertices - 1) * 2;

    if (no_of_vertices <= 1) {
        igraph_vector_int_clear(path);
    } else {
        IGRAPH_CHECK(igraph_vector_int_resize(path, no_of_items));

        i = no_of_vertices - 1;
        j = no_of_items - 1;
        VECTOR(*path)[j] = VECTOR(*path)[i];
        while (i > 1) {
            i--; j--;
            VECTOR(*path)[j] = VECTOR(*path)[i];
            j--;
            VECTOR(*path)[j] = VECTOR(*path)[i];
        }
    }

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_vertex_path_from_edge_path
 * \brief Converts a path of edge IDs to the traversed vertex IDs.
 *
 * </para><para>
 * This function is useful when you have a sequence of edge IDs representing a
 * continuous path in a graph and you would like to obtain the vertex IDs that
 * the path traverses. The function is used implicitly by several shortest path
 * related functions to convert a path of edge IDs to the corresponding
 * representation that describes the path in terms of vertex IDs instead.
 *
 * \param  graph  the graph that the edge IDs refer to
 * \param  start  the start vertex of the path
 * \param  edge_path  the sequence of edge IDs that describe the path
 * \param  vertex_path  the sequence of vertex IDs traversed will be returned here
 * \return Error code: \c IGRAPH_ENOMEM if there is not enough memory,
 *         \c IGRAPH_EINVAL if the edge path does not start at the given vertex
 *         or if there is at least one edge whose start vertex does not match
 *         the end vertex of the previous edge
 */
igraph_error_t igraph_vertex_path_from_edge_path(
   const igraph_t *graph, igraph_integer_t start,
   const igraph_vector_int_t *edge_path, igraph_vector_int_t *vertex_path,
   igraph_neimode_t mode
) {
    igraph_integer_t i, no_of_edges;
    igraph_bool_t directed = igraph_is_directed(graph);
    igraph_bool_t next_edge_ok;
    igraph_integer_t next_start;

    igraph_vector_int_clear(vertex_path);

    no_of_edges = igraph_vector_int_size(edge_path);
    IGRAPH_CHECK(igraph_vector_int_reserve(vertex_path, no_of_edges + 1));

    if (!directed) {
        mode = IGRAPH_ALL;
    }

    for (i = 0; i < no_of_edges; i++) {
        igraph_integer_t from = IGRAPH_FROM(graph, VECTOR(*edge_path)[i]);
        igraph_integer_t to = IGRAPH_TO(graph, VECTOR(*edge_path)[i]);

        igraph_vector_int_push_back(vertex_path, start);  /* reserved */

        switch (mode) {
            case IGRAPH_OUT:
                next_edge_ok = from == start;
                next_start = to;
                break;

            case IGRAPH_IN:
                next_edge_ok = to == start;
                next_start = from;
                break;

            case IGRAPH_ALL:
                if (from == start) {
                    next_edge_ok = true;
                    next_start = to;
                } else if (to == start) {
                    next_edge_ok = true;
                    next_start = from;
                } else {
                    next_edge_ok = false;
                }
                break;

            default:
                IGRAPH_ERROR("Invalid neighborhood mode.", IGRAPH_EINVMODE);
        }

        if (!next_edge_ok) {
            IGRAPH_ERROR("Edge IDs do not form a continuous path.", IGRAPH_EINVAL);
        }

        start = next_start;
    }

    igraph_vector_int_push_back(vertex_path, start);  /* reserved */

    return IGRAPH_SUCCESS;
}
