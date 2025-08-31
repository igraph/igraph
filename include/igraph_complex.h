/*
   igraph library.
   Copyright (C) 2010-2025  The igraph development team <igraph@igraph.org>

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

#ifndef IGRAPH_COMPLEX_H
#define IGRAPH_COMPLEX_H

#include "igraph_decls.h"
#include "igraph_types.h"

IGRAPH_BEGIN_C_DECLS

typedef struct igraph_complex_t {
    igraph_real_t dat[2];
} igraph_complex_t;

#define IGRAPH_REAL(x) ((x).dat[0])
#define IGRAPH_IMAG(x) ((x).dat[1])
#define IGRAPH_COMPLEX_EQ(x,y) ((x).dat[0]==(y).dat[0] && (x).dat[1]==(y).dat[1])

IGRAPH_EXPORT IGRAPH_FUNCATTR_CONST igraph_complex_t igraph_complex(igraph_real_t x, igraph_real_t y);
IGRAPH_EXPORT IGRAPH_FUNCATTR_CONST igraph_complex_t igraph_complex_polar(igraph_real_t r, igraph_real_t theta);

IGRAPH_EXPORT IGRAPH_FUNCATTR_CONST igraph_bool_t igraph_complex_almost_equals(igraph_complex_t z1,
                                                         igraph_complex_t z2,
                                                         igraph_real_t eps);

IGRAPH_EXPORT IGRAPH_FUNCATTR_CONST igraph_real_t igraph_complex_arg(igraph_complex_t z);

IGRAPH_EXPORT IGRAPH_FUNCATTR_CONST igraph_real_t igraph_complex_abs(igraph_complex_t z);
IGRAPH_EXPORT IGRAPH_FUNCATTR_CONST igraph_real_t igraph_complex_logabs(igraph_complex_t z);

IGRAPH_EXPORT IGRAPH_FUNCATTR_CONST igraph_complex_t igraph_complex_add(igraph_complex_t z1,
                                                                        igraph_complex_t z2);
IGRAPH_EXPORT IGRAPH_FUNCATTR_CONST igraph_complex_t igraph_complex_sub(igraph_complex_t z1,
                                                                        igraph_complex_t z2);
IGRAPH_EXPORT IGRAPH_FUNCATTR_CONST igraph_complex_t igraph_complex_mul(igraph_complex_t z1,
                                                                        igraph_complex_t z2);
IGRAPH_EXPORT IGRAPH_FUNCATTR_CONST igraph_complex_t igraph_complex_div(igraph_complex_t z1,
                                                                        igraph_complex_t z2);

IGRAPH_EXPORT IGRAPH_FUNCATTR_CONST igraph_complex_t igraph_complex_add_real(igraph_complex_t z,
                                                                             igraph_real_t x);
IGRAPH_EXPORT IGRAPH_FUNCATTR_CONST igraph_complex_t igraph_complex_add_imag(igraph_complex_t z,
                                                                             igraph_real_t y);
IGRAPH_EXPORT IGRAPH_FUNCATTR_CONST igraph_complex_t igraph_complex_sub_real(igraph_complex_t z,
                                                                             igraph_real_t x);
IGRAPH_EXPORT IGRAPH_FUNCATTR_CONST igraph_complex_t igraph_complex_sub_imag(igraph_complex_t z,
                                                                             igraph_real_t y);
IGRAPH_EXPORT IGRAPH_FUNCATTR_CONST igraph_complex_t igraph_complex_mul_real(igraph_complex_t z,
                                                                             igraph_real_t x);
IGRAPH_EXPORT IGRAPH_FUNCATTR_CONST igraph_complex_t igraph_complex_mul_imag(igraph_complex_t z,
                                                                             igraph_real_t y);
IGRAPH_EXPORT IGRAPH_FUNCATTR_CONST igraph_complex_t igraph_complex_div_real(igraph_complex_t z,
                                                                             igraph_real_t x);
IGRAPH_EXPORT IGRAPH_FUNCATTR_CONST igraph_complex_t igraph_complex_div_imag(igraph_complex_t z,
                                                                             igraph_real_t y);

IGRAPH_EXPORT IGRAPH_FUNCATTR_CONST igraph_complex_t igraph_complex_conj(igraph_complex_t z);
IGRAPH_EXPORT IGRAPH_FUNCATTR_CONST igraph_complex_t igraph_complex_neg(igraph_complex_t z);
IGRAPH_EXPORT IGRAPH_FUNCATTR_CONST igraph_complex_t igraph_complex_inv(igraph_complex_t z);

IGRAPH_EXPORT IGRAPH_FUNCATTR_CONST igraph_complex_t igraph_complex_sqrt(igraph_complex_t z);
IGRAPH_EXPORT IGRAPH_FUNCATTR_CONST igraph_complex_t igraph_complex_sqrt_real(igraph_real_t x);
IGRAPH_EXPORT IGRAPH_FUNCATTR_CONST igraph_complex_t igraph_complex_exp(igraph_complex_t z);
IGRAPH_EXPORT IGRAPH_FUNCATTR_CONST igraph_complex_t igraph_complex_pow(igraph_complex_t z1,
                                                                        igraph_complex_t z2);
IGRAPH_EXPORT IGRAPH_FUNCATTR_CONST igraph_complex_t igraph_complex_pow_real(igraph_complex_t z,
                                                                             igraph_real_t x);
IGRAPH_EXPORT IGRAPH_FUNCATTR_CONST igraph_complex_t igraph_complex_log(igraph_complex_t z);
IGRAPH_EXPORT IGRAPH_FUNCATTR_CONST igraph_complex_t igraph_complex_log10(igraph_complex_t z);
IGRAPH_EXPORT IGRAPH_FUNCATTR_CONST igraph_complex_t igraph_complex_log_b(igraph_complex_t z,
                                                                          igraph_complex_t b);

IGRAPH_EXPORT IGRAPH_FUNCATTR_CONST igraph_complex_t igraph_complex_sin(igraph_complex_t z);
IGRAPH_EXPORT IGRAPH_FUNCATTR_CONST igraph_complex_t igraph_complex_cos(igraph_complex_t z);
IGRAPH_EXPORT IGRAPH_FUNCATTR_CONST igraph_complex_t igraph_complex_tan(igraph_complex_t z);
IGRAPH_EXPORT IGRAPH_FUNCATTR_CONST igraph_complex_t igraph_complex_sec(igraph_complex_t z);
IGRAPH_EXPORT IGRAPH_FUNCATTR_CONST igraph_complex_t igraph_complex_csc(igraph_complex_t z);
IGRAPH_EXPORT IGRAPH_FUNCATTR_CONST igraph_complex_t igraph_complex_cot(igraph_complex_t z);

IGRAPH_EXPORT int igraph_complex_printf(igraph_complex_t val);
IGRAPH_EXPORT int igraph_complex_fprintf(FILE *file, igraph_complex_t val);
IGRAPH_EXPORT int igraph_complex_printf_aligned(int width, igraph_complex_t val);
IGRAPH_EXPORT int igraph_complex_fprintf_aligned(FILE *file, int width, igraph_complex_t val);
IGRAPH_EXPORT int igraph_complex_snprintf(char *str, size_t size, igraph_complex_t val);

IGRAPH_END_C_DECLS

#endif
