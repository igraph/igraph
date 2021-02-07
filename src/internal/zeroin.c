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

/* from GNU R's zeroin.c, minor modifications by Gabor Csardi */

/* from NETLIB c/brent.shar with max.iter, add'l info and convergence
   details hacked in by Peter Dalgaard */

/*************************************************************************
 *              C math library
 * function ZEROIN - obtain a function zero within the given range
 *
 * Input
 *  double zeroin(ax,bx,f,info,Tol,Maxit)
 *  double ax;          Root will be seeked for within
 *  double bx;          a range [ax,bx]
 *  double (*f)(double x, void *info); Name of the function whose zero
 *                  will be seeked for
 *  void *info;         Add'l info passed to f
 *  double *Tol;            Acceptable tolerance for the root
 *                  value.
 *                  May be specified as 0.0 to cause
 *                  the program to find the root as
 *                  accurate as possible
 *
 *  int *Maxit;         Max. iterations
 *
 *
 * Output
 *  Zeroin returns an estimate for the root with accuracy
 *  4*EPSILON*abs(x) + tol
 *  *Tol returns estimated precision
 *  *Maxit returns actual # of iterations, or -1 if maxit was
 *      reached without convergence.
 *
 * Algorithm
 *  G.Forsythe, M.Malcolm, C.Moler, Computer methods for mathematical
 *  computations. M., Mir, 1980, p.180 of the Russian edition
 *
 *  The function makes use of the bisection procedure combined with
 *  the linear or quadric inverse interpolation.
 *  At every step program operates on three abscissae - a, b, and c.
 *  b - the last and the best approximation to the root
 *  a - the last but one approximation
 *  c - the last but one or even earlier approximation than a that
 *      1) |f(b)| <= |f(c)|
 *      2) f(b) and f(c) have opposite signs, i.e. b and c confine
 *         the root
 *  At every step Zeroin selects one of the two new approximations, the
 *  former being obtained by the bisection procedure and the latter
 *  resulting in the interpolation (if a,b, and c are all different
 *  the quadric interpolation is utilized, otherwise the linear one).
 *  If the latter (i.e. obtained by the interpolation) point is
 *  reasonable (i.e. lies within the current interval [b,c] not being
 *  too close to the boundaries) it is accepted. The bisection result
 *  is used in the other case. Therefore, the range of uncertainty is
 *  ensured to be reduced at least by the factor 1.6
 *
 ************************************************************************
 */

#include "igraph_nongraph.h"
#include "igraph_types.h"

#include "core/interruption.h"

#include <float.h>
#include <math.h>

#define EPSILON DBL_EPSILON

int igraph_zeroin(              /* An estimate of the root */
    igraph_real_t *ax,          /* Left border | of the range   */
    igraph_real_t *bx,          /* Right border| the root is seeked*/
    igraph_real_t (*f)(igraph_real_t x, void *info),    /* Function under investigation */
    void *info,             /* Add'l info passed on to f    */
    igraph_real_t *Tol,         /* Acceptable tolerance     */
    int *Maxit,             /* Max # of iterations */
    igraph_real_t *res) {               /* Result is stored here */
    igraph_real_t a, b, c,      /* Abscissae, descr. see above  */
                  fa, fb, fc;         /* f(a), f(b), f(c) */
    igraph_real_t tol;
    int maxit;

    a = *ax;  b = *bx;  fa = (*f)(a, info);  fb = (*f)(b, info);
    c = a;   fc = fa;
    maxit = *Maxit + 1; tol = * Tol;

    /* First test if we have found a root at an endpoint */
    if (fa == 0.0) {
        *Tol = 0.0;
        *Maxit = 0;
        *res = a;
        return 0;
    }
    if (fb ==  0.0) {
        *Tol = 0.0;
        *Maxit = 0;
        *res = b;
        return 0;
    }

    while (maxit--) {   /* Main iteration loop  */
        igraph_real_t prev_step = b - a;  /* Distance from the last but one
                       to the last approximation    */
        igraph_real_t tol_act;      /* Actual tolerance     */
        igraph_real_t p;        /* Interpolation step is calcu- */
        igraph_real_t q;        /* lated in the form p/q; divi-
                     * sion operations is delayed
                     * until the last moment    */
        igraph_real_t new_step;     /* Step at this iteration   */

        IGRAPH_ALLOW_INTERRUPTION();

        if ( fabs(fc) < fabs(fb) ) {
            /* Swap data for b to be the    */
            a = b;  b = c;  c = a;  /* best approximation       */
            fa = fb;  fb = fc;  fc = fa;
        }
        tol_act = 2 * EPSILON * fabs(b) + tol / 2;
        new_step = (c - b) / 2;

        if ( fabs(new_step) <= tol_act || fb == (igraph_real_t)0 ) {
            *Maxit -= maxit;
            *Tol = fabs(c - b);
            *res = b;
            return 0;           /* Acceptable approx. is found  */
        }

        /* Decide if the interpolation can be tried */
        if ( fabs(prev_step) >= tol_act /* If prev_step was large enough*/
             && fabs(fa) > fabs(fb) ) {
            /* and was in true direction,
                         * Interpolation may be tried   */
            register igraph_real_t t1, cb, t2;
            cb = c - b;
            if ( a == c ) {     /* If we have only two distinct */
                /* points linear interpolation  */
                t1 = fb / fa;   /* can only be applied      */
                p = cb * t1;
                q = 1.0 - t1;
            } else {        /* Quadric inverse interpolation*/

                q = fa / fc;  t1 = fb / fc;  t2 = fb / fa;
                p = t2 * ( cb * q * (q - t1) - (b - a) * (t1 - 1.0) );
                q = (q - 1.0) * (t1 - 1.0) * (t2 - 1.0);
            }
            if ( p > (igraph_real_t)0 ) { /* p was calculated with the */
                q = -q;    /* opposite sign; make p positive */
            } else {        /* and assign possible minus to */
                p = -p;    /* q              */
            }

            if ( p < (0.75 * cb * q - fabs(tol_act * q) / 2) /* If b+p/q falls in [b,c]*/
                 && p < fabs(prev_step * q / 2) ) { /* and isn't too large  */
                new_step = p / q;
            }         /* it is accepted
                         * If p/q is too large then the
                         * bisection procedure can
                         * reduce [b,c] range to more
                         * extent */
        }

        if ( fabs(new_step) < tol_act) { /* Adjust the step to be not less*/
            if ( new_step > (igraph_real_t)0 ) { /* than tolerance       */
                new_step = tol_act;
            } else {
                new_step = -tol_act;
            }
        }
        a = b;  fa = fb;            /* Save the previous approx. */
        b += new_step;  fb = (*f)(b, info); /* Do step to a new approxim. */
        if ( (fb > 0 && fc > 0) || (fb < 0 && fc < 0) ) {
            /* Adjust c for it to have a sign opposite to that of b */
            c = a;  fc = fa;
        }

    }
    /* failed! */
    *Tol = fabs(c - b);
    *Maxit = -1;
    *res = b;
    return IGRAPH_DIVERGED;
}

