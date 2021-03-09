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
 *    This file implements algorithm 5.8 of the above reference.
 *    The optimal_partition function returns the minimizing partition
 *    with size 'nt' of the objective function ||v-Pv||, where P is
 *    a problem-specific projector. So far, Symmetric (matrix=1),
 *    Laplacian (matrix=2) and Stochastic (matrix=3) projectors
 *    have been implemented (the cost_matrix function below).
 *    In the stochastic case, 'p' is expected to be a valid propability
 *    vector. In all other cases, 'p' is ignored and can be set to NULL.
 *    The group labels are given in 'gr' as positive consecutive integers
 *    starting from 0.
 */

#include "scg_headers.h"

#include "igraph_error.h"
#include "igraph_memory.h"
#include "igraph_matrix.h"
#include "igraph_vector.h"
#include "igraph_qsort.h"

int igraph_i_optimal_partition(const igraph_real_t *v, int *gr, int n,
                               int nt, int matrix, const igraph_real_t *p,
                               igraph_real_t *value) {

    int i, non_ties, q, j, l, part_ind, col;
    igraph_i_scg_indval_t *vs = IGRAPH_CALLOC(n, igraph_i_scg_indval_t);
    igraph_real_t *Cv, temp, sumOfSquares;
    igraph_vector_t ps;
    igraph_matrix_t F;
    igraph_matrix_int_t Q;

    /*-----------------------------------------------
      -----Sorts v and counts non-ties-----------------
      -----------------------------------------------*/

    if (!vs) {
        IGRAPH_ERROR("SCG error", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, vs);

    for (i = 0; i < n; i++) {
        vs[i].val = v[i];
        vs[i].ind = i;
    }

    igraph_qsort(vs, (size_t) n, sizeof(igraph_i_scg_indval_t),
          igraph_i_compare_ind_val);

    non_ties = 1;
    for (i = 1; i < n; i++) {
        if (vs[i].val < vs[i - 1].val - 1e-14 ||
            vs[i].val > vs[i - 1].val + 1e-14) {
            non_ties++;
        }
    }

    if (nt >= non_ties) {
        IGRAPH_ERROR("`Invalid number of intervals, should be smaller than "
                     "number of unique values in V", IGRAPH_EINVAL);
    }

    /*------------------------------------------------
      ------Computes Cv, the matrix of costs------------
      ------------------------------------------------*/
    Cv = igraph_i_real_sym_matrix(n);
    if (!Cv) {
        IGRAPH_ERROR("SCG error", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, Cv);

    /* if stochastic SCG orders p */
    if (matrix == 3) {
        IGRAPH_VECTOR_INIT_FINALLY(&ps, n);
        for (i = 0; i < n; i++) {
            VECTOR(ps)[i] = p[vs[i].ind];
        }
    }

    IGRAPH_CHECK(igraph_i_cost_matrix(Cv, vs, n, matrix, &ps));
    if (matrix == 3) {
        igraph_vector_destroy(&ps);
        IGRAPH_FINALLY_CLEAN(1);
    }
    /*-------------------------------------------------
      -------Fills up matrices F and Q-------------------
      -------------------------------------------------*/
    /*here j also is a counter but the use of unsigned variables
      is to be proscribed in "for (unsigned int j=...;j>=0;j--)",
      for such loops never ends!*/

    IGRAPH_MATRIX_INIT_FINALLY(&F, nt, n);
    IGRAPH_CHECK(igraph_matrix_int_init(&Q, nt, n));
    IGRAPH_FINALLY(igraph_matrix_int_destroy, &Q);

    for (i = 0; i < n; i++) {
        MATRIX(Q, 0, i)++;
    }
    for (i = 0; i < nt; i++) {
        MATRIX(Q, i, i) = i + 1;
    }

    for (i = 0; i < n; i++) {
        MATRIX(F, 0, i) = igraph_i_real_sym_mat_get(Cv, 0, i);
    }

    for (i = 1; i < nt; i++)
        for (j = i + 1; j < n; j++) {
            MATRIX(F, i, j) = MATRIX(F, i - 1, i - 1) + igraph_i_real_sym_mat_get(Cv, i, j);
            MATRIX(Q, i, j) = 2;

            for (q = i - 1; q <= j - 1; q++) {
                temp = MATRIX(F, i - 1, q) + igraph_i_real_sym_mat_get(Cv, q + 1, j);
                if (temp < MATRIX(F, i, j)) {
                    MATRIX(F, i, j) = temp;
                    MATRIX(Q, i, j) = q + 2;
                }
            }
        }
    igraph_i_free_real_sym_matrix(Cv);
    IGRAPH_FINALLY_CLEAN(1);

    /*--------------------------------------------------
      -------Back-tracks through Q to work out the groups-
      --------------------------------------------------*/
    part_ind = nt;
    col = n - 1;

    for (j = nt - 1; j >= 0; j--) {
        for (i = MATRIX(Q, j, col) - 1; i <= col; i++) {
            gr[vs[i].ind] = part_ind - 1;
        }
        if (MATRIX(Q, j, col) != 2) {
            col = MATRIX(Q, j, col) - 2;
            part_ind -= 1;
        } else {
            if (j > 1) {
                for (l = 0; l <= (j - 1); l++) {
                    gr[vs[l].ind] = l;
                }
                break;
            } else {
                col = MATRIX(Q, j, col) - 2;
                part_ind -= 1;
            }
        }
    }

    sumOfSquares = MATRIX(F, nt - 1, n - 1);

    igraph_matrix_destroy(&F);
    igraph_matrix_int_destroy(&Q);
    IGRAPH_FREE(vs);
    IGRAPH_FINALLY_CLEAN(3);

    if (value) {
        *value = sumOfSquares;
    }
    return 0;
}

int igraph_i_cost_matrix(igraph_real_t*Cv, const igraph_i_scg_indval_t *vs,
                         int n,  int matrix, const igraph_vector_t *ps) {

    /* if symmetric of Laplacian SCG -> same Cv */
    if (matrix == 1 || matrix == 2) {
        int i, j;
        igraph_vector_t w, w2;

        IGRAPH_VECTOR_INIT_FINALLY(&w, n + 1);
        IGRAPH_VECTOR_INIT_FINALLY(&w2, n + 1);

        VECTOR(w)[1] = vs[0].val;
        VECTOR(w2)[1] = vs[0].val * vs[0].val;

        for (i = 2; i <= n; i++) {
            VECTOR(w)[i] = VECTOR(w)[i - 1] + vs[i - 1].val;
            VECTOR(w2)[i] = VECTOR(w2)[i - 1] + vs[i - 1].val * vs[i - 1].val;
        }

        for (i = 0; i < n; i++) {
            for (j = i + 1; j < n; j++) {
                igraph_real_t v = (VECTOR(w2)[j + 1] - VECTOR(w2)[i]) -
                                  (VECTOR(w)[j + 1] - VECTOR(w)[i]) * (VECTOR(w)[j + 1] - VECTOR(w)[i]) /
                                  (j - i + 1);
                igraph_i_real_sym_mat_set(Cv, i, j, v);
            }
        }

        igraph_vector_destroy(&w);
        igraph_vector_destroy(&w2);
        IGRAPH_FINALLY_CLEAN(2);
    }
    /* if stochastic */
    /* TODO: optimize it to O(n^2) instead of O(n^3) (as above) */
    if (matrix == 3) {
        int i, j, k;
        igraph_real_t t1, t2;
        for (i = 0; i < n; i++) {
            for (j = i + 1; j < n; j++) {
                t1 = t2 = 0;
                for (k = i; k < j; k++) {
                    t1 += VECTOR(*ps)[k];
                    t2 += VECTOR(*ps)[k] * vs[k].val;
                }
                t1 = t2 / t1;
                t2 = 0;
                for (k = i; k < j; k++) {
                    t2 += (vs[k].val - t1) * (vs[k].val - t1);
                }
                igraph_i_real_sym_mat_set(Cv, i, j, t2);
            }
        }
    }

    return 0;
}
