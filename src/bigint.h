/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2009-2012  Gabor Csardi <csardi.gabor@gmail.com>
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
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#ifndef IGRAPH_BIGINT_H
#define IGRAPH_BIGINT_H

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
    #define __BEGIN_DECLS extern "C" {
    #define __END_DECLS }
#else
    #define __BEGIN_DECLS /* empty */
    #define __END_DECLS /* empty */
#endif

#include "igraph_types.h"
#include "igraph_vector.h"
#include "bignum.h"

#include <stdio.h>

/* Arbitrary precision integer */

#define BASE_LIMB
#include "igraph_pmt.h"
#include "igraph_vector_type.h"
#include "igraph_vector_pmt.h"
#include "igraph_pmt_off.h"
#undef BASE_LIMB

__BEGIN_DECLS

typedef struct igraph_biguint_t {
    igraph_vector_limb_t v;
} igraph_biguint_t;

#define IGRAPH_BIGUINT_DEFAULT_SIZE 5

int igraph_biguint_init(igraph_biguint_t *b);
void igraph_biguint_destroy(igraph_biguint_t *b);
int igraph_biguint_copy(igraph_biguint_t *to, igraph_biguint_t *from);

int igraph_biguint_extend(igraph_biguint_t *b, limb_t l);

int igraph_biguint_size(igraph_biguint_t *b);
int igraph_biguint_resize(igraph_biguint_t *b, int newlength);
int igraph_biguint_reserve(igraph_biguint_t *b, int length);

int igraph_biguint_zero(igraph_biguint_t *b);
int igraph_biguint_set_limb(igraph_biguint_t *b, int value);

igraph_real_t igraph_biguint_get(igraph_biguint_t *b);

int igraph_biguint_compare_limb(igraph_biguint_t *b, limb_t l);
int igraph_biguint_compare(igraph_biguint_t *left, igraph_biguint_t *right);
igraph_bool_t igraph_biguint_equal(igraph_biguint_t *left, igraph_biguint_t *right);
igraph_bool_t igraph_biguint_bigger(igraph_biguint_t *left,
                                    igraph_biguint_t *right);
igraph_bool_t igraph_biguint_biggerorequal(igraph_biguint_t *left,
        igraph_biguint_t *right);

int igraph_biguint_inc(igraph_biguint_t *res, igraph_biguint_t *b);
int igraph_biguint_dec(igraph_biguint_t *res, igraph_biguint_t *b);

int igraph_biguint_add_limb(igraph_biguint_t *res, igraph_biguint_t *b,
                            limb_t l);
int igraph_biguint_sub_limb(igraph_biguint_t *res, igraph_biguint_t *b,
                            limb_t l);
int igraph_biguint_mul_limb(igraph_biguint_t *res, igraph_biguint_t *b,
                            limb_t l);

int igraph_biguint_add(igraph_biguint_t *res, igraph_biguint_t *left,
                       igraph_biguint_t *right);
int igraph_biguint_sub(igraph_biguint_t *res, igraph_biguint_t *left,
                       igraph_biguint_t *right);
int igraph_biguint_mul(igraph_biguint_t *res, igraph_biguint_t *left,
                       igraph_biguint_t *right);
int igraph_biguint_div(igraph_biguint_t *q, igraph_biguint_t *r,
                       igraph_biguint_t *u, igraph_biguint_t *v);

int igraph_biguint_print(igraph_biguint_t *b);
int igraph_biguint_fprint(igraph_biguint_t *b, FILE *file);

__END_DECLS

#endif
