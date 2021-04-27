/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2011-12  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard st, Cambridge, MA, 02138 USA

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
 *    The intervals_method and intervals_plus_kmeans implements the
 *    methods of sec. 5.3.2 and sec. 5.3.3 of the above reference.
 *    They take an eigenvector 'v' as parameter and a vector 'breaks'
 *    of length 'nb', which provide the intervals used to cut 'v'.
 *    Then all components of 'v' that fall into the same interval are
 *    assigned the same group label in 'gr'. The group labels are
 *    positive consecutive integers starting from 0.
 *    The intervals_method function is adapted from bincode of the R
 *    base package.
 *    The intervals_plus_kmeans is initialized with regularly-spaced
 *    breaks, which rougly corresponds to the intervals_method. Then
 *    kmeans minimizes iteratively the objective function until it gets
 *    stuck in a (usually) local minimum, or until 'itermax' is reached.
 *    So far, the breaks_computation function allows computation of
 *    constant bins, as used in intervals_method, and of equidistant
 *    centers as used in intervals_plus_kmeans.
 */

#include "scg_headers.h"

#include "igraph_error.h"
#include "igraph_types.h"
#include "igraph_vector.h"

int igraph_i_intervals_plus_kmeans(const igraph_vector_t *v, int *gr,
                                   int n, int n_interv,
                                   int maxiter) {
    int i;
    igraph_vector_t centers;

    IGRAPH_VECTOR_INIT_FINALLY(&centers, n_interv);

    igraph_i_breaks_computation(v, &centers, n_interv, 2);
    IGRAPH_CHECK(igraph_i_kmeans_Lloyd(v, n, 1, &centers, n_interv, gr,
                                       maxiter));

    /*renumber the groups*/
    for (i = 0; i < n; i++) {
        gr[i] = gr[i] - 1;
    }

    igraph_vector_destroy(&centers);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

int igraph_i_intervals_method(const igraph_vector_t *v, int *gr, int n,
                              int n_interv) {
    int i, lo, hi, new;
    const int lft = 1;
    const int include_border = 1;
    igraph_vector_t breaks;

    IGRAPH_VECTOR_INIT_FINALLY(&breaks, n_interv + 1);

    IGRAPH_CHECK(igraph_i_breaks_computation(v, &breaks, n_interv + 1, 1));

    for (i = 0; i < n; i++) {
        lo = 0;
        hi = n_interv;
        if (VECTOR(*v)[i] <  VECTOR(breaks)[lo] ||
            VECTOR(breaks)[hi] < VECTOR(*v)[i] ||
            (VECTOR(*v)[i] == VECTOR(breaks)[lft ? hi : lo] && !include_border)) {
            /* Do nothing */
        } else {
            while (hi - lo >= 2) {
                new = (hi + lo) / 2;
                if (VECTOR(*v)[i] > VECTOR(breaks)[new] ||
                    (lft && VECTOR(*v)[i] == VECTOR(breaks)[new])) {
                    lo = new;
                } else {
                    hi = new;
                }
            }
            gr[i] = lo;
        }
    }
    igraph_vector_destroy(&breaks);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

int igraph_i_breaks_computation(const igraph_vector_t *v,
                                igraph_vector_t *breaks,
                                int nb, int method) {
    int i;
    igraph_real_t eps, vmin, vmax;
    igraph_vector_minmax(v, &vmin, &vmax);

    if (vmax == vmin) {
        IGRAPH_ERROR("There is only one (repeated) value in argument 'v' "
                     "of bin_size_computation()", IGRAPH_EINVAL);
    }

    if (nb < 2) {
        IGRAPH_ERROR("'nb' in bin_size_computation() must be >= 2",
                     IGRAPH_EINVAL);
    }

    switch (method) {
    case 1: /* constant bins for fixed-size intervals method */
        eps = (vmax - vmin) / (igraph_real_t)(nb - 1);
        VECTOR(*breaks)[0] = vmin;
        for (i = 1; i < nb - 1; i++) {
            VECTOR(*breaks)[i] = VECTOR(*breaks)[i - 1] + eps;
        }
        VECTOR(*breaks)[nb - 1] = vmax;
        break;
    case 2: /* equidistant centers for kmeans */
        eps = (vmax - vmin) / (igraph_real_t)nb;
        VECTOR(*breaks)[0] = vmin + eps / 2.;
        for (i = 1; i < nb; i++) {
            VECTOR(*breaks)[i] = VECTOR(*breaks)[i - 1] + eps;
        }
        break;
    /* TODO: implement logarithmic binning for power-law-like distributions */
    default:
        IGRAPH_ERROR("Internal SCG error, this should ot happen",
                     IGRAPH_FAILURE);
    }

    return 0;
}
