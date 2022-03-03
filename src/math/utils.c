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

#include "igraph_types.h"
#include "igraph_nongraph.h"

#include "core/math.h"

#include "config.h"

#include <math.h>
#include <float.h>

int igraph_finite(double x) {
    return isfinite(x);
}

int igraph_chebyshev_init(const double *dos, int nos, double eta) {
    int i, ii;
    double err;

    if (nos < 1) {
        return 0;
    }

    err = 0.0;
    i = 0;          /* just to avoid compiler warnings */
    for (ii = 1; ii <= nos; ii++) {
        i = nos - ii;
        err += fabs(dos[i]);
        if (err > eta) {
            return i;
        }
    }
    return i;
}

double igraph_chebyshev_eval(double x, const double *a, const int n) {
    double b0, b1, b2, twox;
    int i;

    if (n < 1 || n > 1000) {
        IGRAPH_WARNING("chebyshev_eval: argument out of domain");
        return IGRAPH_NAN;
    }

    if (x < -1.1 || x > 1.1) {
        IGRAPH_WARNING("chebyshev_eval: argument out of domain");
        return IGRAPH_NAN;
    }

    twox = x * 2;
    b2 = b1 = 0;
    b0 = 0;
    for (i = 1; i <= n; i++) {
        b2 = b1;
        b1 = b0;
        b0 = twox * b1 - b2 + a[n - i];
    }
    return (b0 - b2) * 0.5;
}

int igraph_is_nan(double x) {
    return isnan(x);
}

int igraph_is_inf(double x) {
    return isinf(x) != 0;
}

int igraph_is_posinf(double x) {
    return isinf(x) && x > 0;
}

int igraph_is_neginf(double x) {
    return isinf(x) && x < 0;
}

/**
 * \function igraph_almost_equals
 * \brief Compare two double-precision floats with a tolerance.
 *
 * Determines whether two double-precision floats are "almost equal"
 * to each other with a given level of tolerance on the relative error.
 *
 * \param  a  The first float.
 * \param  b  The second float.
 * \param  eps  The level of tolerance on the relative error. The relative
 *         error is defined as <code>abs(a-b) / (abs(a) + abs(b))</code>. The
 *         two numbers are considered equal if this is less than \c eps.
 *
 * \return True if the two floats are nearly equal to each other within
 *         the given level of tolerance, false otherwise.
 */
igraph_bool_t igraph_almost_equals(double a, double b, double eps) {
    return igraph_cmp_epsilon(a, b, eps) == 0 ? 1 : 0;
}

/* Use value-safe floating point math for igraph_cmp_epsilon() with
 * the Intel compiler.
 *
 * The Intel compiler rewrites arithmetic expressions for faster
 * evaluation by default. In the below function, it will evaluate
 * (eps * fabs(a) + eps * fabs(b)) as eps*(fabs(a) + fabs(b)).
 * However, this code path is taken precisely when fabs(a) + fabs(b)
 * overflows, thus this rearrangement of the expression causes
 * the function to return incorrect results, and some test failures.
 * To avoid this, we switch the Intel compiler to "precise" mode.
 */
#ifdef __INTEL_COMPILER
#pragma float_control(push)
#pragma float_control (precise, on)
#endif

/**
 * \function igraph_cmp_epsilon
 * \brief Compare two double-precision floats with a tolerance.
 *
 * Determines whether two double-precision floats are "almost equal"
 * to each other with a given level of tolerance on the relative error.
 *
 * \param  a  The first float.
 * \param  b  The second float.
 * \param  eps  The level of tolerance on the relative error. The relative
 *         error is defined as <code>abs(a-b) / (abs(a) + abs(b))</code>. The
 *         two numbers are considered equal if this is less than \c eps.
 *
 * \return Zero if the two floats are nearly equal to each other within
 *         the given level of tolerance, positive number if the first float is
 *         larger, negative number if the second float is larger.
 */
int igraph_cmp_epsilon(double a, double b, double eps) {
    double diff;
    double abs_diff;
    double sum;

    if (a == b) {
        /* shortcut, handles infinities */
        return 0;
    }

    diff = a - b;
    abs_diff = fabs(diff);
    sum = fabs(a) + fabs(b);

    if (a == 0 || b == 0 || sum < DBL_MIN) {
        /* a or b is zero or both are extremely close to it; relative
         * error is less meaningful here so just compare it with
         * epsilon */
        return abs_diff < (eps * DBL_MIN) ? 0 : (diff < 0 ? -1 : 1);
    } else if (!isfinite(sum)) {
        /* addition overflow, so presumably |a| and |b| are both large; use a
         * different formulation */
        return (abs_diff < (eps * fabs(a) + eps * fabs(b))) ? 0 : (diff < 0 ? -1 : 1);
    } else {
        return (abs_diff / sum < eps) ? 0 : (diff < 0 ? -1 : 1);
    }
}

#ifdef __INTEL_COMPILER
#pragma float_control(pop)
#endif
