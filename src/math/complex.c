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

#include "igraph_complex.h"
#include "core/math.h"
#include <math.h>

/**
 * \example igraph_complex.c
 */

igraph_complex_t igraph_complex(igraph_real_t x, igraph_real_t y) {
    igraph_complex_t res;
    IGRAPH_REAL(res) = x;
    IGRAPH_IMAG(res) = y;
    return res;
}

igraph_complex_t igraph_complex_polar(igraph_real_t r, igraph_real_t theta) {
    igraph_complex_t res;
    IGRAPH_REAL(res) = r * cos(theta);
    IGRAPH_IMAG(res) = r * sin(theta);
    return res;
}

igraph_bool_t igraph_complex_eq_tol(igraph_complex_t z1,
                                    igraph_complex_t z2,
                                    igraph_real_t tol) {
    if (fabs(IGRAPH_REAL(z1) - IGRAPH_REAL(z2)) > tol ||
        fabs(IGRAPH_IMAG(z1) - IGRAPH_IMAG(z2)) > tol) {
        return 0;
    }
    return 1;
}

igraph_real_t igraph_complex_mod(igraph_complex_t z) {
    igraph_real_t x = IGRAPH_REAL(z);
    igraph_real_t y = IGRAPH_IMAG(z);
    return hypot(x, y);
}

igraph_real_t igraph_complex_arg(igraph_complex_t z) {
    igraph_real_t x = IGRAPH_REAL(z);
    igraph_real_t y = IGRAPH_IMAG(z);
    if (x == 0.0 && y == 0.0) {
        return 0.0;
    }
    return atan2(y, x);
}

igraph_complex_t igraph_complex_add(igraph_complex_t z1,
                                    igraph_complex_t z2) {
    igraph_complex_t res;
    IGRAPH_REAL(res) = IGRAPH_REAL(z1) + IGRAPH_REAL(z2);
    IGRAPH_IMAG(res) = IGRAPH_IMAG(z1) + IGRAPH_IMAG(z2);
    return res;
}

igraph_complex_t igraph_complex_sub(igraph_complex_t z1,
                                    igraph_complex_t z2) {
    igraph_complex_t res;
    IGRAPH_REAL(res) = IGRAPH_REAL(z1) - IGRAPH_REAL(z2);
    IGRAPH_IMAG(res) = IGRAPH_IMAG(z1) - IGRAPH_IMAG(z2);
    return res;
}

igraph_complex_t igraph_complex_mul(igraph_complex_t z1,
                                    igraph_complex_t z2) {
    igraph_complex_t res;
    IGRAPH_REAL(res) = IGRAPH_REAL(z1) * IGRAPH_REAL(z2) -
                       IGRAPH_IMAG(z1) * IGRAPH_IMAG(z2);
    IGRAPH_IMAG(res) = IGRAPH_REAL(z1) * IGRAPH_IMAG(z2) +
                       IGRAPH_IMAG(z1) * IGRAPH_REAL(z2);
    return res;
}

igraph_complex_t igraph_complex_div(igraph_complex_t z1,
                                    igraph_complex_t z2) {
    igraph_complex_t res;
    igraph_real_t z1r = IGRAPH_REAL(z1), z1i = IGRAPH_IMAG(z1);
    igraph_real_t z2r = IGRAPH_REAL(z2), z2i = IGRAPH_IMAG(z2);
    igraph_real_t s = 1.0 / igraph_complex_abs(z2);
    igraph_real_t sz2r = s * z2r;
    igraph_real_t sz2i = s * z2i;
    IGRAPH_REAL(res) = (z1r * sz2r + z1i * sz2i) * s;
    IGRAPH_IMAG(res) = (z1i * sz2r - z1r * sz2i) * s;
    return res;
}

igraph_complex_t igraph_complex_add_real(igraph_complex_t z,
        igraph_real_t x) {
    igraph_complex_t res;
    IGRAPH_REAL(res) = IGRAPH_REAL(z) + x;
    IGRAPH_IMAG(res) = IGRAPH_IMAG(z);
    return res;
}

igraph_complex_t igraph_complex_add_imag(igraph_complex_t z,
        igraph_real_t y) {
    igraph_complex_t res;
    IGRAPH_REAL(res) = IGRAPH_REAL(z);
    IGRAPH_IMAG(res) = IGRAPH_IMAG(z) + y;
    return res;
}

igraph_complex_t igraph_complex_sub_real(igraph_complex_t z,
        igraph_real_t x) {
    igraph_complex_t res;
    IGRAPH_REAL(res) = IGRAPH_REAL(z) - x;
    IGRAPH_IMAG(res) = IGRAPH_IMAG(z);
    return res;
}

igraph_complex_t igraph_complex_sub_imag(igraph_complex_t z,
        igraph_real_t y) {
    igraph_complex_t res;
    IGRAPH_REAL(res) = IGRAPH_REAL(z);
    IGRAPH_IMAG(res) = IGRAPH_IMAG(z) - y;
    return res;
}

igraph_complex_t igraph_complex_mul_real(igraph_complex_t z,
        igraph_real_t x) {
    igraph_complex_t res;
    IGRAPH_REAL(res) = IGRAPH_REAL(z) * x;
    IGRAPH_IMAG(res) = IGRAPH_IMAG(z) * x;
    return res;
}

igraph_complex_t igraph_complex_mul_imag(igraph_complex_t z,
        igraph_real_t y) {
    igraph_complex_t res;
    IGRAPH_REAL(res) = - IGRAPH_IMAG(z) * y;
    IGRAPH_IMAG(res) =   IGRAPH_REAL(z) * y;
    return res;
}

igraph_complex_t igraph_complex_div_real(igraph_complex_t z,
        igraph_real_t x) {
    igraph_complex_t res;
    IGRAPH_REAL(res) = IGRAPH_REAL(z) / x;
    IGRAPH_IMAG(res) = IGRAPH_IMAG(z) / x;
    return res;
}

igraph_complex_t igraph_complex_div_imag(igraph_complex_t z,
        igraph_real_t y) {
    igraph_complex_t res;
    IGRAPH_REAL(res) =   IGRAPH_IMAG(z) / y;
    IGRAPH_IMAG(res) = - IGRAPH_REAL(z) / y;
    return res;
}

igraph_complex_t igraph_complex_conj(igraph_complex_t z) {
    igraph_complex_t res;
    IGRAPH_REAL(res) =   IGRAPH_REAL(z);
    IGRAPH_IMAG(res) = - IGRAPH_IMAG(z);
    return res;
}

igraph_complex_t igraph_complex_neg(igraph_complex_t z) {
    igraph_complex_t res;
    IGRAPH_REAL(res) = - IGRAPH_REAL(z);
    IGRAPH_IMAG(res) = - IGRAPH_IMAG(z);
    return res;
}

igraph_complex_t igraph_complex_inv(igraph_complex_t z) {
    igraph_complex_t res;
    igraph_real_t s = 1.0 / igraph_complex_abs(z);
    IGRAPH_REAL(res) =   (IGRAPH_REAL(z) * s) * s;
    IGRAPH_IMAG(res) = - (IGRAPH_IMAG(z) * s) * s;
    return res;
}

igraph_real_t igraph_complex_abs(igraph_complex_t z) {
    return hypot(IGRAPH_REAL(z), IGRAPH_IMAG(z));
}

igraph_real_t igraph_complex_logabs(igraph_complex_t z) {
    igraph_real_t xabs = fabs(IGRAPH_REAL(z));
    igraph_real_t yabs = fabs(IGRAPH_IMAG(z));
    igraph_real_t max, u;
    if (xabs >= yabs) {
        max = xabs;
        u = yabs / xabs;
    } else {
        max = yabs;
        u = xabs / yabs;
    }
    return log (max) + 0.5 * log1p (u * u);
}

igraph_complex_t igraph_complex_sqrt(igraph_complex_t z) {
    igraph_complex_t res;

    if (IGRAPH_REAL(z) == 0.0 && IGRAPH_IMAG(z) == 0.0) {
        IGRAPH_REAL(res) = IGRAPH_IMAG(res) = 0.0;
    } else {
        igraph_real_t x = fabs (IGRAPH_REAL(z));
        igraph_real_t y = fabs (IGRAPH_IMAG(z));
        igraph_real_t w;
        if (x >= y)  {
            igraph_real_t t = y / x;
            w = sqrt (x) * sqrt (0.5 * (1.0 + sqrt (1.0 + t * t)));
        } else {
            igraph_real_t t = x / y;
            w = sqrt (y) * sqrt (0.5 * (t + sqrt (1.0 + t * t)));
        }

        if (IGRAPH_REAL(z) >= 0.0) {
            igraph_real_t ai = IGRAPH_IMAG(z);
            IGRAPH_REAL(res) = w;
            IGRAPH_IMAG(res) = ai / (2.0 * w);
        } else {
            igraph_real_t ai = IGRAPH_IMAG(z);
            igraph_real_t vi = (ai >= 0) ? w : -w;
            IGRAPH_REAL(res) = ai / (2.0 * vi);
            IGRAPH_IMAG(res) = vi;
        }
    }

    return res;
}

igraph_complex_t igraph_complex_sqrt_real(igraph_real_t x) {
    igraph_complex_t res;
    if (x >= 0) {
        IGRAPH_REAL(res) = sqrt(x);
        IGRAPH_IMAG(res) = 0.0;
    } else {
        IGRAPH_REAL(res) = 0.0;
        IGRAPH_IMAG(res) = sqrt(-x);
    }
    return res;
}

igraph_complex_t igraph_complex_exp(igraph_complex_t z) {
    igraph_real_t rho   = exp(IGRAPH_REAL(z));
    igraph_real_t theta = IGRAPH_IMAG(z);
    igraph_complex_t res;
    IGRAPH_REAL(res) = rho * cos(theta);
    IGRAPH_IMAG(res) = rho * sin(theta);
    return res;
}

igraph_complex_t igraph_complex_pow(igraph_complex_t z1,
                                    igraph_complex_t z2) {
    igraph_complex_t res;

    if (IGRAPH_REAL(z1) == 0 && IGRAPH_IMAG(z1) == 0.0) {
        if (IGRAPH_REAL(z2) == 0 && IGRAPH_IMAG(z2) == 0.0) {
            IGRAPH_REAL(res) = 1.0;
            IGRAPH_IMAG(res) = 0.0;
        } else {
            IGRAPH_REAL(res) = IGRAPH_IMAG(res) = 0.0;
        }
    } else if (IGRAPH_REAL(z2) == 1.0 && IGRAPH_IMAG(z2) == 0.0) {
        IGRAPH_REAL(res) = IGRAPH_REAL(z1);
        IGRAPH_IMAG(res) = IGRAPH_IMAG(z1);
    } else if (IGRAPH_REAL(z2) == -1.0 && IGRAPH_IMAG(z2) == 0.0) {
        res = igraph_complex_inv(z1);
    } else {
        igraph_real_t logr = igraph_complex_logabs (z1);
        igraph_real_t theta = igraph_complex_arg (z1);
        igraph_real_t z2r = IGRAPH_REAL(z2), z2i = IGRAPH_IMAG(z2);
        igraph_real_t rho = exp (logr * z2r - z2i * theta);
        igraph_real_t beta = theta * z2r + z2i * logr;
        IGRAPH_REAL(res) = rho * cos(beta);
        IGRAPH_IMAG(res) = rho * sin(beta);
    }

    return res;
}

igraph_complex_t igraph_complex_pow_real(igraph_complex_t z,
        igraph_real_t x) {
    igraph_complex_t res;
    if (IGRAPH_REAL(z) == 0.0 && IGRAPH_IMAG(z) == 0.0) {
        if (x == 0) {
            IGRAPH_REAL(res) = 1.0;
            IGRAPH_IMAG(res) = 0.0;
        } else {
            IGRAPH_REAL(res) = IGRAPH_IMAG(res) = 0.0;
        }
    } else {
        igraph_real_t logr = igraph_complex_logabs(z);
        igraph_real_t theta = igraph_complex_arg(z);
        igraph_real_t rho = exp (logr * x);
        igraph_real_t beta = theta * x;
        IGRAPH_REAL(res) = rho * cos(beta);
        IGRAPH_IMAG(res) = rho * sin(beta);
    }
    return res;
}

igraph_complex_t igraph_complex_log(igraph_complex_t z) {
    igraph_complex_t res;
    IGRAPH_REAL(res) = igraph_complex_logabs(z);
    IGRAPH_IMAG(res) = igraph_complex_arg(z);
    return res;
}

igraph_complex_t igraph_complex_log10(igraph_complex_t z) {
    return igraph_complex_mul_real(igraph_complex_log(z), 1 / log(10.0));
}

igraph_complex_t igraph_complex_log_b(igraph_complex_t z,
                                      igraph_complex_t b) {
    return igraph_complex_div (igraph_complex_log(z), igraph_complex_log(b));
}

igraph_complex_t igraph_complex_sin(igraph_complex_t z) {
    igraph_real_t zr = IGRAPH_REAL(z);
    igraph_real_t zi = IGRAPH_IMAG(z);
    igraph_complex_t res;
    if (zi == 0.0) {
        IGRAPH_REAL(res) = sin(zr);
        IGRAPH_IMAG(res) = 0.0;
    } else {
        IGRAPH_REAL(res) = sin(zr) * cosh(zi);
        IGRAPH_IMAG(res) = cos(zr) * sinh(zi);
    }
    return res;
}

igraph_complex_t igraph_complex_cos(igraph_complex_t z) {
    igraph_real_t zr = IGRAPH_REAL(z);
    igraph_real_t zi = IGRAPH_IMAG(z);
    igraph_complex_t res;
    if (zi == 0.0) {
        IGRAPH_REAL(res) = cos(zr);
        IGRAPH_IMAG(res) = 0.0;
    } else {
        IGRAPH_REAL(res) = cos(zr) * cosh(zi);
        IGRAPH_IMAG(res) = sin(zr) * sinh(-zi);
    }
    return res;
}

igraph_complex_t igraph_complex_tan(igraph_complex_t z) {
    igraph_real_t zr = IGRAPH_REAL(z);
    igraph_real_t zi = IGRAPH_IMAG(z);
    igraph_complex_t res;
    if (fabs (zi) < 1) {
        igraph_real_t D = pow (cos (zr), 2.0) + pow (sinh (zi), 2.0);
        IGRAPH_REAL(res) = 0.5 * sin (2 * zr) / D;
        IGRAPH_IMAG(res) = 0.5 * sinh (2 * zi) / D;
    } else {
        igraph_real_t u = exp (-zi);
        igraph_real_t C = 2 * u / (1 - pow (u, 2.0));
        igraph_real_t D = 1 + pow (cos (zr), 2.0) * pow (C, 2.0);
        igraph_real_t S = pow (C, 2.0);
        igraph_real_t T = 1.0 / tanh (zi);
        IGRAPH_REAL(res) = 0.5 * sin (2 * zr) * S / D;
        IGRAPH_IMAG(res) = T / D;
    }
    return res;
}

igraph_complex_t igraph_complex_sec(igraph_complex_t z) {
    return igraph_complex_inv(igraph_complex_cos(z));
}

igraph_complex_t igraph_complex_csc(igraph_complex_t z) {
    return igraph_complex_inv(igraph_complex_sin(z));
}

igraph_complex_t igraph_complex_cot(igraph_complex_t z) {
    return igraph_complex_inv(igraph_complex_tan(z));
}

