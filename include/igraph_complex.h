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

#ifndef IGRAPH_COMPLEX_H
#define IGRAPH_COMPLEX_H

#include "igraph_decls.h"
#include "igraph_types.h"

__BEGIN_DECLS

typedef struct igraph_complex_t {
    igraph_real_t dat[2];
} igraph_complex_t;

#define IGRAPH_REAL(x) ((x).dat[0])
#define IGRAPH_IMAG(x) ((x).dat[1])
#define IGRAPH_COMPLEX_EQ(x,y) ((x).dat[0]==(y).dat[0] && (x).dat[1]==(y).dat[1])

DECLDIR igraph_complex_t igraph_complex(igraph_real_t x, igraph_real_t y);
DECLDIR igraph_complex_t igraph_complex_polar(igraph_real_t r, igraph_real_t theta);

DECLDIR igraph_bool_t igraph_complex_eq_tol(igraph_complex_t z1,
        igraph_complex_t z2,
        igraph_real_t tol);

DECLDIR igraph_real_t igraph_complex_mod(igraph_complex_t z);
DECLDIR igraph_real_t igraph_complex_arg(igraph_complex_t z);

DECLDIR igraph_real_t igraph_complex_abs(igraph_complex_t z);
DECLDIR igraph_real_t igraph_complex_logabs(igraph_complex_t z);

DECLDIR igraph_complex_t igraph_complex_add(igraph_complex_t z1,
        igraph_complex_t z2);
DECLDIR igraph_complex_t igraph_complex_sub(igraph_complex_t z1,
        igraph_complex_t z2);
DECLDIR igraph_complex_t igraph_complex_mul(igraph_complex_t z1,
        igraph_complex_t z2);
DECLDIR igraph_complex_t igraph_complex_div(igraph_complex_t z1,
        igraph_complex_t z2);

DECLDIR igraph_complex_t igraph_complex_add_real(igraph_complex_t z,
        igraph_real_t x);
DECLDIR igraph_complex_t igraph_complex_add_imag(igraph_complex_t z,
        igraph_real_t y);
DECLDIR igraph_complex_t igraph_complex_sub_real(igraph_complex_t z,
        igraph_real_t x);
DECLDIR igraph_complex_t igraph_complex_sub_imag(igraph_complex_t z,
        igraph_real_t y);
DECLDIR igraph_complex_t igraph_complex_mul_real(igraph_complex_t z,
        igraph_real_t x);
DECLDIR igraph_complex_t igraph_complex_mul_imag(igraph_complex_t z,
        igraph_real_t y);
DECLDIR igraph_complex_t igraph_complex_div_real(igraph_complex_t z,
        igraph_real_t x);
DECLDIR igraph_complex_t igraph_complex_div_imag(igraph_complex_t z,
        igraph_real_t y);

DECLDIR igraph_complex_t igraph_complex_conj(igraph_complex_t z);
DECLDIR igraph_complex_t igraph_complex_neg(igraph_complex_t z);
DECLDIR igraph_complex_t igraph_complex_inv(igraph_complex_t z);

DECLDIR igraph_complex_t igraph_complex_sqrt(igraph_complex_t z);
DECLDIR igraph_complex_t igraph_complex_sqrt_real(igraph_real_t x);
DECLDIR igraph_complex_t igraph_complex_exp(igraph_complex_t z);
DECLDIR igraph_complex_t igraph_complex_pow(igraph_complex_t z1,
        igraph_complex_t z2);
DECLDIR igraph_complex_t igraph_complex_pow_real(igraph_complex_t z,
        igraph_real_t x);
DECLDIR igraph_complex_t igraph_complex_log(igraph_complex_t z);
DECLDIR igraph_complex_t igraph_complex_log10(igraph_complex_t z);
DECLDIR igraph_complex_t igraph_complex_log_b(igraph_complex_t z,
        igraph_complex_t b);

DECLDIR igraph_complex_t igraph_complex_sin(igraph_complex_t z);
DECLDIR igraph_complex_t igraph_complex_cos(igraph_complex_t z);
DECLDIR igraph_complex_t igraph_complex_tan(igraph_complex_t z);
DECLDIR igraph_complex_t igraph_complex_sec(igraph_complex_t z);
DECLDIR igraph_complex_t igraph_complex_csc(igraph_complex_t z);
DECLDIR igraph_complex_t igraph_complex_cot(igraph_complex_t z);

__END_DECLS

#endif
