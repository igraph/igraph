/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2010-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#include <igraph.h>
/* This definition ensures math symbols are also available when
 * compiling with MSVC.
 */
#define _USE_MATH_DEFINES
#include <math.h>

#include "test_utilities.inc"

#define ARE 4
#define AIM 5
#define BRE 6
#define BIM 2

int main() {

    igraph_complex_t a = igraph_complex(ARE, AIM);
    igraph_complex_t b = igraph_complex(BRE, BIM);
    igraph_complex_t c, d, e;

    /* polar, mod, arg */
    c = igraph_complex_polar(igraph_complex_mod(a), igraph_complex_arg(a));
    IGRAPH_ASSERT(igraph_complex_eq_tol(a, c, 1e-14));

    /* add */
    c = igraph_complex_add(a, b);
    IGRAPH_ASSERT(IGRAPH_REAL(c) == ARE + BRE && IGRAPH_IMAG(c) == AIM + BIM);

    /* sub */
    c = igraph_complex_sub(a, b);
    IGRAPH_ASSERT(IGRAPH_REAL(c) == ARE - BRE && IGRAPH_IMAG(c) == AIM - BIM);

    /* mul */
    c = igraph_complex_mul(a, b);
    IGRAPH_ASSERT(IGRAPH_REAL(c) == ARE * BRE - AIM * BIM);
    IGRAPH_ASSERT(IGRAPH_IMAG(c) == ARE * BIM + AIM * BRE);

    /* div */
    c = igraph_complex_div(a, b);
    c = igraph_complex_mul(c, b);
    IGRAPH_ASSERT(igraph_complex_eq_tol(a, c, 1e-14));

    /* add_real */
    c = igraph_complex_add_real(a, IGRAPH_REAL(b));
    IGRAPH_ASSERT(IGRAPH_REAL(c) == IGRAPH_REAL(a) + IGRAPH_REAL(b));
    IGRAPH_ASSERT(IGRAPH_IMAG(c) == IGRAPH_IMAG(a));

    /* add_imag */
    c = igraph_complex_add_imag(a, IGRAPH_IMAG(b));
    IGRAPH_ASSERT(IGRAPH_REAL(c) == IGRAPH_REAL(a));
    IGRAPH_ASSERT(IGRAPH_IMAG(c) == IGRAPH_IMAG(a) + IGRAPH_IMAG(b));

    /* sub_real */
    c = igraph_complex_sub_real(a, IGRAPH_REAL(b));
    IGRAPH_ASSERT(IGRAPH_REAL(c) == IGRAPH_REAL(a) - IGRAPH_REAL(b));
    IGRAPH_ASSERT(IGRAPH_IMAG(c) == IGRAPH_IMAG(a));

    /* sub_imag */
    c = igraph_complex_sub_imag(a, IGRAPH_IMAG(b));
    IGRAPH_ASSERT(IGRAPH_REAL(c) == IGRAPH_REAL(a));
    IGRAPH_ASSERT(IGRAPH_IMAG(c) == IGRAPH_IMAG(a) - IGRAPH_IMAG(b));

    /* mul_real */
    c = igraph_complex_mul_real(a, IGRAPH_REAL(b));
    IGRAPH_ASSERT(IGRAPH_REAL(c) == IGRAPH_REAL(a) * IGRAPH_REAL(b));
    IGRAPH_ASSERT(IGRAPH_IMAG(c) == IGRAPH_IMAG(a) * IGRAPH_REAL(b));

    /* mul_imag */
    c = igraph_complex_mul_imag(a, IGRAPH_REAL(b));
    IGRAPH_ASSERT(IGRAPH_REAL(c) == - IGRAPH_IMAG(a) * IGRAPH_REAL(b));
    IGRAPH_ASSERT(IGRAPH_IMAG(c) == IGRAPH_REAL(a) * IGRAPH_REAL(b));

    /* div_real */
    c = igraph_complex_div_real(a, IGRAPH_REAL(b));
    IGRAPH_ASSERT(fabs(IGRAPH_REAL(c) - IGRAPH_REAL(a) / IGRAPH_REAL(b)) < 1e-15);
    IGRAPH_ASSERT(fabs(IGRAPH_IMAG(c) - IGRAPH_IMAG(a) / IGRAPH_REAL(b)) < 1e-15);

    /* div_imag */
    c = igraph_complex_div_imag(a, IGRAPH_IMAG(b));
    IGRAPH_ASSERT(IGRAPH_REAL(c) == IGRAPH_IMAG(a) / IGRAPH_IMAG(b));
    IGRAPH_ASSERT(IGRAPH_IMAG(c) == - IGRAPH_REAL(a) / IGRAPH_IMAG(b));

    /* conj */
    c = igraph_complex_conj(a);
    IGRAPH_ASSERT(IGRAPH_REAL(c) == ARE && IGRAPH_IMAG(c) == -AIM);

    /* neg */
    c = igraph_complex_neg(a);
    IGRAPH_ASSERT(IGRAPH_REAL(c) == - IGRAPH_REAL(a));
    IGRAPH_ASSERT(IGRAPH_IMAG(c) == - IGRAPH_IMAG(a));

    /* inv */
    c = igraph_complex_inv(a);
    d = igraph_complex(1.0, 0.0);
    e = igraph_complex_div(d, a);
    IGRAPH_ASSERT(igraph_complex_eq_tol(c, e, 1e-14));

    /* abs */
    IGRAPH_ASSERT(igraph_complex_abs(a) == igraph_complex_mod(a));

    /* logabs */

    /* sqrt */
    c = igraph_complex_sqrt(a);
    d = igraph_complex_mul(c, c);
    IGRAPH_ASSERT(igraph_complex_eq_tol(a, d, 1e-14));

    /* sqrt_real */
    c = igraph_complex_sqrt(igraph_complex(-1.0, 0.0));
    d = igraph_complex_sqrt_real(-1.0);
    IGRAPH_ASSERT(igraph_complex_eq_tol(c, d, 1e-14));

    /* exp */
    c = igraph_complex_exp(igraph_complex(0.0, M_PI));
    IGRAPH_ASSERT(igraph_complex_eq_tol(c, igraph_complex(-1.0, 0.0), 1e-14));

    /* pow */
    c = igraph_complex_pow(igraph_complex(M_E, 0.0), igraph_complex(0.0, M_PI));
    IGRAPH_ASSERT(igraph_complex_eq_tol(c, igraph_complex(-1.0, 0.0), 1e-14));

    /* pow_real */
    c = igraph_complex_pow_real(a, 2.0);
    d = igraph_complex_mul(a, a);
    IGRAPH_ASSERT(igraph_complex_eq_tol(c, d, 1e-12));

    /* log */
    c = igraph_complex_exp(igraph_complex_log(a));
    IGRAPH_ASSERT(igraph_complex_eq_tol(a, c, 1e-14));

    /* log10 */
    c = igraph_complex_pow(igraph_complex(10.0, 0), igraph_complex_log10(a));
    IGRAPH_ASSERT(igraph_complex_eq_tol(a, c, 1e-14));

    /* log_b */
    c = igraph_complex_pow(b, igraph_complex_log_b(a, b));
    IGRAPH_ASSERT(igraph_complex_eq_tol(a, c, 1e-14));

    /* sin, cos */
    c = igraph_complex_sin(a);
    d = igraph_complex_cos(a);
    e = igraph_complex_add(igraph_complex_mul(c, c), igraph_complex_mul(d, d));
    IGRAPH_ASSERT(igraph_complex_eq_tol(e, igraph_complex(1.0, 0.0), 1e-11));

    /* tan */
    c = igraph_complex_tan(a);
    d = igraph_complex_div(igraph_complex_sin(a), igraph_complex_cos(a));
    IGRAPH_ASSERT(igraph_complex_eq_tol(c, d, 1e-14));

    /* sec */
    c = igraph_complex_sec(a);
    d = igraph_complex_inv(igraph_complex_cos(a));
    IGRAPH_ASSERT(igraph_complex_eq_tol(c, d, 1e-14));

    /* csc */
    c = igraph_complex_csc(a);
    d = igraph_complex_inv(igraph_complex_sin(a));
    IGRAPH_ASSERT(igraph_complex_eq_tol(c, d, 1e-14));

    /* cot */
    c = igraph_complex_tan(a);
    d = igraph_complex_div(igraph_complex_sin(a), igraph_complex_cos(a));
    IGRAPH_ASSERT(igraph_complex_eq_tol(d, c, 1e-14));

    VERIFY_FINALLY_STACK();
    return 0;
}
