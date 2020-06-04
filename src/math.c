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

#include <math.h>
#include <float.h>
#include <stdarg.h>
#include "config.h"
#include "igraph_math.h"
#include "igraph_types.h"

#ifdef _MSC_VER
    #define isinf(x) (!_finite(x) && !_isnan(x))
#endif

int igraph_finite(double x) {
#if HAVE_DECL_ISFINITE
    return isfinite(x);
#elif HAVE_FINITE == 1
    return finite(x);
#else
    /* neither finite nor isfinite work. Do we really need the AIX exception? */
# ifdef _AIX
#  include <fp.h>
    return FINITE(x);
# else
    return (!isnan(x) & (x != IGRAPH_POSINFINITY) & (x != IGRAPH_NEGINFINITY));
# endif
#endif
}

double igraph_log2(const double a) {
    return log(a) / log(2.0);
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
        IGRAPH_NAN;
    }

    if (x < -1.1 || x > 1.1) {
        IGRAPH_NAN;
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

double igraph_log1p(double x) {
    /* series for log1p on the interval -.375 to .375
     *                   with weighted error   6.35e-32
     *                    log weighted error  31.20
     *              significant figures required  30.93
     *               decimal places required  32.01
     */
    static const double alnrcs[43] = {
        +.10378693562743769800686267719098e+1,
            -.13364301504908918098766041553133e+0,
            +.19408249135520563357926199374750e-1,
            -.30107551127535777690376537776592e-2,
            +.48694614797154850090456366509137e-3,
            -.81054881893175356066809943008622e-4,
            +.13778847799559524782938251496059e-4,
            -.23802210894358970251369992914935e-5,
            +.41640416213865183476391859901989e-6,
            -.73595828378075994984266837031998e-7,
            +.13117611876241674949152294345011e-7,
            -.23546709317742425136696092330175e-8,
            +.42522773276034997775638052962567e-9,
            -.77190894134840796826108107493300e-10,
            +.14075746481359069909215356472191e-10,
            -.25769072058024680627537078627584e-11,
            +.47342406666294421849154395005938e-12,
            -.87249012674742641745301263292675e-13,
            +.16124614902740551465739833119115e-13,
            -.29875652015665773006710792416815e-14,
            +.55480701209082887983041321697279e-15,
            -.10324619158271569595141333961932e-15,
            +.19250239203049851177878503244868e-16,
            -.35955073465265150011189707844266e-17,
            +.67264542537876857892194574226773e-18,
            -.12602624168735219252082425637546e-18,
            +.23644884408606210044916158955519e-19,
            -.44419377050807936898878389179733e-20,
            +.83546594464034259016241293994666e-21,
            -.15731559416479562574899253521066e-21,
            +.29653128740247422686154369706666e-22,
            -.55949583481815947292156013226666e-23,
            +.10566354268835681048187284138666e-23,
            -.19972483680670204548314999466666e-24,
            +.37782977818839361421049855999999e-25,
            -.71531586889081740345038165333333e-26,
            +.13552488463674213646502024533333e-26,
            -.25694673048487567430079829333333e-27,
            +.48747756066216949076459519999999e-28,
            -.92542112530849715321132373333333e-29,
            +.17578597841760239233269760000000e-29,
            -.33410026677731010351377066666666e-30,
            +.63533936180236187354180266666666e-31,
        };

    static IGRAPH_THREAD_LOCAL int nlnrel = 0;
    static IGRAPH_THREAD_LOCAL double xmin = 0.0;

    if (xmin == 0.0) {
        xmin = -1 + sqrt(DBL_EPSILON);    /*was sqrt(d1mach(4)); */
    }
    if (nlnrel == 0) { /* initialize chebychev coefficients */
        nlnrel = igraph_chebyshev_init(alnrcs, 43, DBL_EPSILON / 20);    /*was .1*d1mach(3)*/
    }

    if (x == 0.) {
        return 0.;    /* speed */
    }
    if (x == -1) {
        return (IGRAPH_NEGINFINITY);
    }
    if (x  < -1) {
        return (IGRAPH_NAN);
    }

    if (fabs(x) <= .375) {
        /* Improve on speed (only);
        again give result accurate to IEEE double precision: */
        if (fabs(x) < .5 * DBL_EPSILON) {
            return x;
        }

        if ( (0 < x && x < 1e-8) || (-1e-9 < x && x < 0)) {
            return x * (1 - .5 * x);
        }
        /* else */
        return x * (1 - x * igraph_chebyshev_eval(x / .375, alnrcs, nlnrel));
    }
    /* else */
    /*     if (x < xmin) { */
    /*  /\* answer less than half precision because x too near -1 *\/ */
    /*         ML_ERROR(ME_PRECISION, "log1p"); */
    /*     } */
    return log(1 + x);
}

long double igraph_fabsl(long double a) {
    if (a < 0) {
        return -a;
    } else {
        return a;
    }
}

double igraph_fmin(double a, double b) {
    if (b < a) {
        return b;
    } else {
        return a;
    }
}

double igraph_i_round(double X) {

    /* NaN */
    if (X != X) {
        return X;
    }

    if (X < 0.0) {
        return floor(X);
    }

    return ceil(X);
}

#ifdef _MSC_VER
/**
 * Internal function, replacement for snprintf
 * Used only in case of the Microsoft Visual C compiler which does not
 * provide a proper sprintf implementation.
 *
 * This implementation differs from the standard in the value returned
 * when the number of characters needed by the output, excluding the
 * terminating '\0' is larger than count
 */
int igraph_i_snprintf(char *buffer, size_t count, const char *format, ...) {
    int n;
    va_list args;
    if (count > 0) {
        va_start(args, format);
        n = _vsnprintf(buffer, count, format, args);
        buffer[count - 1] = 0;
        va_end(args);
    } else {
        n = 0;
    }
    return n;
}

#endif

int igraph_is_nan(double x) {
    return isnan(x);
}

int igraph_is_inf(double x) {
    return isinf(x) != 0;
}

int igraph_is_posinf(double x) {
    return isinf(x) == 1;
}

int igraph_is_neginf(double x) {
    return isinf(x) == -1;
}

/**
 * \function igraph_almost_equals
 * Compare two double-precision floats with a tolerance
 *
 * Determines whether two double-precision floats are "almost equal"
 * to each other with a given level of tolerance on the relative error.
 *
 * \param  a  the first float
 * \param  b  the second float
 * \param  eps  the level of tolerance on the relative error. The relative
 *         error is defined as \c "abs(a-b) / (abs(a) + abs(b))". The
 *         two numbers are considered equal if this is less than \c eps.
 *
 * \return nonzero if the two floats are nearly equal to each other within
 *         the given level of tolerance, zero otherwise
 */
int igraph_almost_equals(double a, double b, double eps) {
    return igraph_cmp_epsilon(a, b, eps) == 0 ? 1 : 0;
}


/**
 * \function igraph_cmp_epsilon
 * Compare two double-precision floats with a tolerance
 *
 * Determines whether two double-precision floats are "almost equal"
 * to each other with a given level of tolerance on the relative error.
 *
 * \param  a  the first float
 * \param  b  the second float
 * \param  eps  the level of tolerance on the relative error. The relative
 *         error is defined as \c "abs(a-b) / (abs(a) + abs(b))". The
 *         two numbers are considered equal if this is less than \c eps.
 *
 * \return zero if the two floats are nearly equal to each other within
 *         the given level of tolerance, positive number if the first float is
 *         larger, negative number if the second float is larger
 */
int igraph_cmp_epsilon(double a, double b, double eps) {
    double diff;
    double abs_diff;

    if (a == b) {
        /* shortcut, handles infinities */
        return 0;
    }

    diff = a - b;
    abs_diff = fabs(diff);

    if (a == 0 || b == 0 || diff < DBL_MIN) {
        /* a or b is zero or both are extremely close to it; relative
         * error is less meaningful here so just compare it with
         * epsilon */
        return abs_diff < (eps * DBL_MIN) ? 0 : (diff < 0 ? -1 : 1);
    } else {
        /* use relative error */
        return (abs_diff / (fabs(a) + fabs(b)) < eps) ? 0 : (diff < 0 ? -1 : 1);
    }
}

