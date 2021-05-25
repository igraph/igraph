/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2007-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#include "igraph_nongraph.h"
#include "core/interruption.h"
#include "igraph_statusbar.h"

#include <math.h>

/* This is from GNU R's optim.c, slightly adapted to igraph */

#define stepredn    0.2
#define acctol      0.0001
#define reltest     10.0
#define FALSE           0
#define TRUE            1

/*  BFGS variable-metric method, based on Pascal code
in J.C. Nash, `Compact Numerical Methods for Computers', 2nd edition,
converted by p2c then re-crafted by B.D. Ripley */

int
igraph_bfgs(igraph_vector_t *b, igraph_real_t *Fmin,
            igraph_scalar_function_t fminfn, igraph_vector_function_t fmingr,
            int maxit, int trace,
            igraph_real_t abstol, igraph_real_t reltol, int nREPORT, void *ex,
            igraph_integer_t *fncount, igraph_integer_t *grcount) {
    int n = (int) igraph_vector_size(b);
    igraph_bool_t accpoint, enough;
    igraph_vector_t g, t, X, c;
    igraph_matrix_t B;        /* Lmatrix really */
    int   count, funcount, gradcount;
    igraph_real_t f, gradproj;
    int   i, j, ilast, iter = 0;
    igraph_real_t s, steplength;
    igraph_real_t D1, D2;

    if (maxit <= 0) {
        *Fmin = fminfn(b, 0, ex);
        *fncount = 1;
        *grcount = 0;
        return 0;
    }

    if (nREPORT <= 0) {
        IGRAPH_ERROR("REPORT must be > 0 (method = \"BFGS\")", IGRAPH_EINVAL);
    }
    IGRAPH_VECTOR_INIT_FINALLY(&g, n);
    IGRAPH_VECTOR_INIT_FINALLY(&t, n);
    IGRAPH_VECTOR_INIT_FINALLY(&X, n);
    IGRAPH_VECTOR_INIT_FINALLY(&c, n);
    IGRAPH_MATRIX_INIT_FINALLY(&B, n, n);
    f = fminfn(b, 0, ex);
    if (!IGRAPH_FINITE(f)) {
        IGRAPH_ERROR("initial value in 'BFGS' is not finite", IGRAPH_DIVERGED);
    }
    if (trace) {
        igraph_statusf("initial  value %f ", 0, f);
    }
    *Fmin = f;
    funcount = gradcount = 1;
    fmingr(b, 0, &g, ex);
    iter++;
    ilast = gradcount;

    do {

        IGRAPH_ALLOW_INTERRUPTION();

        if (ilast == gradcount) {
            for (i = 0; i < n; i++) {
                for (j = 0; j < i; j++) {
                    MATRIX(B, i, j) = 0.0;
                }
                MATRIX(B, i, i) = 1.0;
            }
        }
        for (i = 0; i < n; i++) {
            VECTOR(X)[i] = VECTOR(*b)[i];
            VECTOR(c)[i] = VECTOR(g)[i];
        }
        gradproj = 0.0;
        for (i = 0; i < n; i++) {
            s = 0.0;
            for (j = 0; j <= i; j++) {
                s -= MATRIX(B, i, j) * VECTOR(g)[j];
            }
            for (j = i + 1; j < n; j++) {
                s -= MATRIX(B, j, i) * VECTOR(g)[j];
            }
            VECTOR(t)[i] = s;
            gradproj += s * VECTOR(g)[i];
        }

        if (gradproj < 0.0) {   /* search direction is downhill */
            steplength = 1.0;
            accpoint = FALSE;
            do {
                count = 0;
                for (i = 0; i < n; i++) {
                    VECTOR(*b)[i] = VECTOR(X)[i] + steplength * VECTOR(t)[i];
                    if (reltest + VECTOR(X)[i] == reltest + VECTOR(*b)[i]) { /* no change */
                        count++;
                    }
                }
                if (count < n) {
                    f = fminfn(b, 0, ex);
                    funcount++;
                    accpoint = IGRAPH_FINITE(f) &&
                               (f <= *Fmin + gradproj * steplength * acctol);
                    if (!accpoint) {
                        steplength *= stepredn;
                    }
                }
            } while (!(count == n || accpoint));
            enough = (f > abstol) &&
                     fabs(f - *Fmin) > reltol * (fabs(*Fmin) + reltol);
            /* stop if value if small or if relative change is low */
            if (!enough) {
                count = n;
                *Fmin = f;
            }
            if (count < n) {/* making progress */
                *Fmin = f;
                fmingr(b, 0, &g, ex);
                gradcount++;
                iter++;
                D1 = 0.0;
                for (i = 0; i < n; i++) {
                    VECTOR(t)[i] = steplength * VECTOR(t)[i];
                    VECTOR(c)[i] = VECTOR(g)[i] - VECTOR(c)[i];
                    D1 += VECTOR(t)[i] * VECTOR(c)[i];
                }
                if (D1 > 0) {
                    D2 = 0.0;
                    for (i = 0; i < n; i++) {
                        s = 0.0;
                        for (j = 0; j <= i; j++) {
                            s += MATRIX(B, i, j) * VECTOR(c)[j];
                        }
                        for (j = i + 1; j < n; j++) {
                            s += MATRIX(B, j, i) * VECTOR(c)[j];
                        }
                        VECTOR(X)[i] = s;
                        D2 += s * VECTOR(c)[i];
                    }
                    D2 = 1.0 + D2 / D1;
                    for (i = 0; i < n; i++) {
                        for (j = 0; j <= i; j++)
                            MATRIX(B, i, j) += (D2 * VECTOR(t)[i] * VECTOR(t)[j]
                                                - VECTOR(X)[i] * VECTOR(t)[j]
                                                - VECTOR(t)[i] * VECTOR(X)[j]) / D1;
                    }
                } else {    /* D1 < 0 */
                    ilast = gradcount;
                }
            } else {  /* no progress */
                if (ilast < gradcount) {
                    count = 0;
                    ilast = gradcount;
                }
            }
        } else {        /* uphill search */
            count = 0;
            if (ilast == gradcount) {
                count = n;
            } else {
                ilast = gradcount;
            }
            /* Resets unless has just been reset */
        }
        if (trace && (iter % nREPORT == 0)) {
            igraph_statusf("iter%4d value %f", 0, iter, f);
        }
        if (iter >= maxit) {
            break;
        }
        if (gradcount - ilast > 2 * n) {
            ilast = gradcount;    /* periodic restart */
        }
    } while (count != n || ilast != gradcount);
    if (trace) {
        igraph_statusf("final  value %f ", 0, *Fmin);
        if (iter < maxit) {
            igraph_status("converged", 0);
        } else {
            igraph_statusf("stopped after %i iterations", 0, iter);
        }
    }
    *fncount = funcount;
    *grcount = gradcount;

    igraph_matrix_destroy(&B);
    igraph_vector_destroy(&c);
    igraph_vector_destroy(&X);
    igraph_vector_destroy(&t);
    igraph_vector_destroy(&g);
    IGRAPH_FINALLY_CLEAN(5);

    return (iter < maxit) ? 0 : IGRAPH_DIVERGED;
}
