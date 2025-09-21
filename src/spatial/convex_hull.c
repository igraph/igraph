/*
   igraph library.
   Copyright (C) 2025  The igraph development team <igraph@igraph.org>

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

#include "igraph_spatial.h"

/**
 * \function igraph_convex_hull_2d
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
igraph_error_t igraph_convex_hull_2d(
        const igraph_matrix_t *data,
        igraph_vector_int_t *resverts,
        igraph_matrix_t *rescoords)
{
    igraph_int_t no_of_nodes;
    igraph_int_t i, pivot_idx = 0, last_idx, before_last_idx, next_idx, j;
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
    IGRAPH_CHECK(igraph_vector_sort_ind(&angles, &order, IGRAPH_ASCENDING));

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
