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

#include "igraph_layout.h"

#include "igraph_blas.h"
#include "igraph_interface.h"
#include "igraph_lapack.h"

#include "core/interruption.h"


/**
 * \ingroup layout
 * \function igraph_layout_align
 * \brief Aligns a graph layout with the coordinate axes.
 *
 * This function centers a vertex layout on the coordinate system origin and
 * rotates the layout to achieve a visually pleasing alignment with the coordinate
 * axes. Doing this is particularly useful with force-directed layouts such as
 * \ref igraph_layout_fruchterman_reingold(). Layouts in arbitrary dimensional
 * spaces are supported.
 *
 * \param graph The graph whose layout is to be aligned.
 * \param layout A matrix whose rows are the coordinates of vertices. It will
 *    be modified in-place.
 * \return Error code.
 *
 * Time complexity: O(|E| + |V|), linear in the number of edges and vertices.
 */

igraph_error_t igraph_layout_align(const igraph_t *graph, igraph_matrix_t *layout) {
    const igraph_int_t vcount = igraph_vcount(graph);
    const igraph_int_t ecount = igraph_ecount(graph);
    const igraph_int_t dim = igraph_matrix_ncol(layout);
    igraph_matrix_t M, Q; /* nematic tensor */
    igraph_real_t norm2_sum; /* sum of squared norms of alignment vectors */
    igraph_matrix_t R; /* rotation matrix consisting of Q's eigenvectors */
    igraph_vector_t lambda; /* Q's eigenvalues */
    igraph_matrix_t temp_layout;

    /* The correction terms are used to remove one vector/edge from the computation
     * of the nematic tensor when it is necessary to break symmetries. */
    igraph_bool_t correction_saved = false;
    igraph_matrix_t M_correction;
    igraph_real_t norm2_sum_correction;

    igraph_bool_t retried = false;
    int iter = 0;

    /* General approach:
     *
     * - Shift the center of mass of the vertices to the origin of our coordinate
     *   system.
     *
     * - Compute the nematic tensor Q_ij = <v_i v_j> - delta_ij / d, where v
     *   are d-dimensional unit vectors used for the alignment, usually along
     *   the edges of the graph, <...> denotes averaging over all such vectors
     *   (i.e. all edges), and delta_ij is the Kronecker symbol.
     *
     *   We use a weighted average, using the squared edge lengths as weights.
     *
     * - Find the principal axes of the nematic tensor and align them to our
     *   coordinate axes.
     *
     * - If there are no edges of non-zero length then we use the vector
     *   representing vertex positions (after centering) to compute the nematic
     *   tensor.
     *
     * - If the nematic tensor is close to zero (i.e. the layout is rotationally
     *   symmetric), we remove a vector v to break the symmetry and re-try.
     *   (See the "correction" terms.)
     */

    if (igraph_matrix_nrow(layout) != vcount) {
        IGRAPH_ERROR("Number of points in layout does not match vertex count.", IGRAPH_EINVAL);
    }

    if (vcount == 0) {
        /* Null graph, nothing to do. */
        return IGRAPH_SUCCESS;
    }

    /* This check is not done for the null graph, which was handled above.
     * Skipping the check is necessary due to the different handling of
     * zero-by-n matrices i various systems. */
    if (dim == 0) {
        IGRAPH_ERROR("Vertex coordinates must be at least one dimensional, "
                     "but received a zero-dimensional input.", IGRAPH_EINVAL);
    }

    /* Shift layout to origin. */
    {
        igraph_vector_t center;
        IGRAPH_VECTOR_INIT_FINALLY(&center, dim);

        IGRAPH_CHECK(igraph_matrix_colsum(layout, &center));
        igraph_vector_scale(&center, 1.0 / vcount);

        for (igraph_int_t j=0; j < dim; j++) {
            for (igraph_int_t i=0; i < vcount; i++) {
                MATRIX(*layout, i, j) -= VECTOR(center)[j];
            }
        }

        igraph_vector_destroy(&center);
        IGRAPH_FINALLY_CLEAN(1);
    }

    /* Nothing more to do for a 1D layout. */
    if (dim == 1) {
        return IGRAPH_SUCCESS;
    }

    IGRAPH_MATRIX_INIT_FINALLY(&M, dim, dim);
    IGRAPH_MATRIX_INIT_FINALLY(&M_correction, dim, dim);

    /* Compute M_ij = sum v_i v_j based on edge vectors. */
    {
        igraph_vector_t edge_vec;
        IGRAPH_VECTOR_INIT_FINALLY(&edge_vec, dim);

        norm2_sum = 0;
        for (igraph_int_t eid=0; eid < ecount; eid++) {
            const igraph_int_t from = IGRAPH_FROM(graph, eid);
            const igraph_int_t to   = IGRAPH_TO(graph, eid);

            if (from == to) continue; /* skip self-loops */

            for (igraph_int_t i=0; i < dim; i++) {
                VECTOR(edge_vec)[i] = MATRIX(*layout, from, i) - MATRIX(*layout, to, i);
            }

            for (igraph_int_t i=0; i < dim; i++) {
                for (igraph_int_t j=0; j < dim; j++) {
                    igraph_real_t m = VECTOR(edge_vec)[i] * VECTOR(edge_vec)[j];
                    MATRIX(M, i, j) += m;
                    if (i == j) {
                        norm2_sum += m;
                    }
                }
            }

            /* Save first non-zero term as a "correction" for later potential removal. */
            if (!correction_saved && norm2_sum > 0) {
                correction_saved = true;
                norm2_sum_correction = norm2_sum;
                IGRAPH_CHECK(igraph_matrix_update(&M_correction, &M));
            }

            IGRAPH_ALLOW_INTERRUPTION_LIMITED(iter, 1 << 12);
        }

        igraph_vector_destroy(&edge_vec);
        IGRAPH_FINALLY_CLEAN(1);
    }

    /* If there are no edges, or all edges are of length zero, the sum of
     * squared norms is zero and cannot be used for normalizing M. We resort
     * to using the position vectors of vertices relative to the center of mass
     * to compute M_ij. Note that the layout has already been centered. */
    if (norm2_sum == 0) {
        /* If norm2_sum == 0 then M is also all-zero, no need to null it explicitly. */
        for (igraph_int_t vid=0; vid < vcount; vid++) {
            for (igraph_int_t i=0; i < dim; i++) {
                for (igraph_int_t j=0; j < dim; j++) {
                    igraph_real_t m = MATRIX(*layout, vid, i) * MATRIX(*layout, vid, j);
                    MATRIX(M, i, j) += m;
                    if (i == j) {
                        norm2_sum += m;
                    }
                }
            }

            /* Save first non-zero term as a "correction" for later potential removal. */
            if (!correction_saved && norm2_sum > 0) {
                correction_saved = true;
                norm2_sum_correction = norm2_sum;
                IGRAPH_CHECK(igraph_matrix_update(&M_correction, &M));
            }

            IGRAPH_ALLOW_INTERRUPTION_LIMITED(iter, 1 << 12);
        }
    }

    /* If all vertices are at the origin, no alignment is possible or necessary.
     * We are done. */
    if (norm2_sum == 0) {
        goto done;
    }

    IGRAPH_ASSERT(correction_saved);

    IGRAPH_MATRIX_INIT_FINALLY(&Q, dim, dim);
    IGRAPH_MATRIX_INIT_FINALLY(&R, dim, dim);
    IGRAPH_VECTOR_INIT_FINALLY(&lambda, dim);

    while (true) {

        /* Finish computing the nematic tensor as
         * Q_ij = M_ij / norm2 - delta_ij / d. */
        IGRAPH_CHECK(igraph_matrix_update(&Q, &M));
        igraph_matrix_scale(&Q, 1.0 / norm2_sum);
        for (igraph_int_t i=0; i < dim; i++) {
            MATRIX(Q, i, i) -= 1.0 / dim;
        }

        /* Compute the eigenvectors of the nematic tensor, which together form
         * the appropriate rotation matrix for alignment. Note that the nematic
         * tensor is symmetric, so its eigenvectors form an orthogonal basis.
         * LAPACK's DSYEVR guarantees that eigenvectors are normalized,
         * "the first M columns of Z contain the orthonormal eigenvectors of the
         * matrix A".
         */
        IGRAPH_CHECK(igraph_lapack_dsyevr(
                &Q, IGRAPH_LAPACK_DSYEV_ALL,
                0, 0, 0, 0, 0, /* ignored when computing all eigenvectors */
                /* abstol */ 1e-6,
                /* eigenvalues */ &lambda, /* eigenvectors */ &R,
                NULL));

        /* Compute the matrix norm, i.e. the largest eigenvalue magnitude,
         * to determine if the nematic tensor Q is close to zero. */
        igraph_real_t matrix_norm = 0;
        for (igraph_int_t i=0; i < dim; i++) {
            igraph_real_t magnitude = fabs(VECTOR(lambda)[i]);
            if (magnitude > matrix_norm) {
                matrix_norm = magnitude;
            }
        }

        /* If the nematic tensor is close to zero, we try to remove a vector once
         * and recompute. The 'retried' variable prevents more than a single retry. */
        if (matrix_norm > 1e-3 || retried) {
            break;
        }

        IGRAPH_CHECK(igraph_matrix_sub(&M, &M_correction));
        norm2_sum -= norm2_sum_correction;

        retried = true;
    }

    /* Rotate the layout. */
    IGRAPH_MATRIX_INIT_FINALLY(&temp_layout, vcount, dim);
    IGRAPH_CHECK(igraph_blas_dgemm(/* transpose 'layout' */ false,
                                   /* transpose 'R' */ false,
                                   /* alpha */ 1.0,
                                   layout, &R,
                                   /* beta */ 0.0,
                                   &temp_layout));

    /* We compute the extent of the coordinate values along each axis,
     * and re-order axes from greatest to smallest extent. This ensures
     * that with a standard plotting setup, plots will be wide rather
     * than tall.
     *
     * At the same time, we copy back 'temp_layout' to 'layout'.
     *
     * Note that the extents are not necessarily ordered the same way
     * as the eigenvalues of the nematic tensor, as there may be isolated
     * vertices which do not contribute to the nematic tensor.
     */
    {
        igraph_vector_t extent;
        igraph_vector_int_t permutation;

        IGRAPH_VECTOR_INIT_FINALLY(&extent, dim);
        IGRAPH_VECTOR_INT_INIT_FINALLY(&permutation, dim);

        for (igraph_int_t j=0; j < dim; j++) {
            igraph_real_t min = IGRAPH_INFINITY, max = -IGRAPH_INFINITY;
            for (igraph_int_t i=0; i < vcount; i++) {
                igraph_real_t c = MATRIX(temp_layout, i, j);
                if (c < min) min = c;
                if (c > max) max = c;
            }
            VECTOR(extent)[j] = max - min;
        }
        IGRAPH_CHECK(igraph_vector_sort_ind(&extent, &permutation, IGRAPH_DESCENDING));

        for (igraph_int_t j=0; j < dim; j++) {
            for (igraph_int_t i=0; i < vcount; i++) {
                MATRIX(*layout, i, j) = MATRIX(temp_layout, i, VECTOR(permutation)[j]);
            }
        }

        igraph_vector_int_destroy(&permutation);
        igraph_vector_destroy(&extent);
        IGRAPH_FINALLY_CLEAN(2);
    }

    igraph_matrix_destroy(&temp_layout);
    igraph_vector_destroy(&lambda);
    igraph_matrix_destroy(&R);
    igraph_matrix_destroy(&Q);
    IGRAPH_FINALLY_CLEAN(4);

done:
    igraph_matrix_destroy(&M_correction);
    igraph_matrix_destroy(&M);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}
