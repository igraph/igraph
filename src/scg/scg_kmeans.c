/*
 *  SCGlib : A C library for the spectral coarse graining of matrices
 *  as described in the paper: Shrinking Matrices while preserving their
 *  eigenpairs with Application to the Spectral Coarse Graining of Graphs.
 *  Preprint available at <http://people.epfl.ch/david.morton>
 *
 *  Copyright (C) 2008 David Morton de Lachapelle <david.morton@a3.epfl.ch>
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 *  02110-1301 USA
 *
 *  DESCRIPTION
 *  -----------
 *    The kmeans_Lloyd function is adapted from the R-stats package.
 *    It perfoms Lloyd's k-means clustering on a p x n data matrix
 *    stored row-wise in a vector 'x'. 'cen' contains k initial centers.
 *    The group label to which each object belongs is stored in 'cl'.
 *    Labels are positive consecutive integers starting from 0.
 *    See also Section 5.3.3 of the above reference.
 */

#include "scg_headers.h"

int igraph_i_kmeans_Lloyd(const igraph_vector_t *x, int n, int p,
                          igraph_vector_t *cen, int k, int *cl, int maxiter) {

    int iter, i, j, c, it, inew = 0;
    igraph_real_t best, dd, tmp;
    int updated;
    igraph_vector_int_t nc;

    IGRAPH_CHECK(igraph_vector_int_init(&nc, k));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &nc);

    for (i = 0; i < n; i++) {
        cl[i] = -1;
    }
    for (iter = 0; iter < maxiter; iter++) {
        updated = 0;
        for (i = 0; i < n; i++) {
            /* find nearest centre for each point */
            best = IGRAPH_INFINITY;
            for (j = 0; j < k; j++) {
                dd = 0.0;
                for (c = 0; c < p; c++) {
                    tmp = VECTOR(*x)[i + n * c] - VECTOR(*cen)[j + k * c];
                    dd += tmp * tmp;
                }
                if (dd < best) {
                    best = dd;
                    inew = j + 1;
                }
            }
            if (cl[i] != inew) {
                updated = 1;
                cl[i] = inew;
            }
        }
        if (!updated) {
            break;
        }

        /* update each centre */
        for (j = 0; j < k * p; j++) {
            VECTOR(*cen)[j] = 0.0;
        }
        for (j = 0; j < k; j++) {
            VECTOR(nc)[j] = 0;
        }
        for (i = 0; i < n; i++) {
            it = cl[i] - 1;
            VECTOR(nc)[it]++;
            for (c = 0; c < p; c++) {
                VECTOR(*cen)[it + c * k] += VECTOR(*x)[i + c * n];
            }
        }
        for (j = 0; j < k * p; j++) {
            VECTOR(*cen)[j] /= VECTOR(nc)[j % k];
        }
    }
    igraph_vector_int_destroy(&nc);
    IGRAPH_FINALLY_CLEAN(1);

    /* convervenge check */
    if (iter >= maxiter - 1) {
        IGRAPH_ERROR("Lloyd k-means did not converge", IGRAPH_FAILURE);
    }

    return 0;
}
