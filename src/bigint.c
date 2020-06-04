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
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include "bigint.h"
#include "igraph_error.h"
#include "igraph_memory.h"

int igraph_biguint_init(igraph_biguint_t *b) {
    IGRAPH_CHECK(igraph_vector_limb_init(&b->v, IGRAPH_BIGUINT_DEFAULT_SIZE));
    igraph_vector_limb_clear(&b->v);
    return 0;
}

void igraph_biguint_destroy(igraph_biguint_t *b) {
    igraph_vector_limb_destroy(&b->v);
}

int igraph_biguint_copy(igraph_biguint_t *to, igraph_biguint_t *from) {
    return igraph_vector_limb_copy(&to->v, &from->v);
}

int igraph_biguint_extend(igraph_biguint_t *b, limb_t l) {
    return igraph_vector_limb_push_back(&b->v, l);
}

int igraph_biguint_size(igraph_biguint_t *b) {
    return (int) igraph_vector_limb_size(&b->v);
}

int igraph_biguint_resize(igraph_biguint_t *b, int newlength) {
    int origlen = igraph_biguint_size(b);
    IGRAPH_CHECK(igraph_vector_limb_resize(&b->v, newlength));
    if (newlength > origlen) {
        memset(VECTOR(b->v) + origlen, 0,
               (size_t) (newlength - origlen) * sizeof(limb_t));
    }
    return 0;
}

int igraph_biguint_reserve(igraph_biguint_t *b, int length) {
    return igraph_vector_limb_reserve(&b->v, length);
}

int igraph_biguint_zero(igraph_biguint_t *b) {
    igraph_vector_limb_clear(&b->v);
    return 0;
}

int igraph_biguint_set_limb(igraph_biguint_t *b, int value) {
    IGRAPH_CHECK(igraph_vector_limb_resize(&b->v, 1));
    VECTOR(b->v)[0] = (limb_t) value;
    return 0;
}

igraph_real_t igraph_biguint_get(igraph_biguint_t *b) {
    int size = igraph_biguint_size(b);
    int i;
    double val = VECTOR(b->v)[size - 1];
    if (size == 0) {
        return 0.0;
    }
    for (i = size - 2; i >= 0; i--) {
        val = val * LIMBMASK + VECTOR(b->v)[i];
        if (!IGRAPH_FINITE(val)) {
            break;
        }
    }
    return val;
}

int igraph_biguint_compare_limb(igraph_biguint_t *b, limb_t l) {
    int n = igraph_biguint_size(b);
    return bn_cmp_limb(VECTOR(b->v), l, (count_t) n);
}

int igraph_biguint_compare(igraph_biguint_t *left, igraph_biguint_t *right) {
    /* bn_cmp requires the two numbers to have the same number of limbs,
       so we do this partially by hand here */
    int size_left = igraph_biguint_size(left);
    int size_right = igraph_biguint_size(right);
    while (size_left > size_right) {
        if (VECTOR(left->v)[--size_left] > 0) {
            return +1;
        }
    }
    while (size_right > size_left) {
        if (VECTOR(right->v)[--size_right] > 0) {
            return -1;
        }
    }
    return bn_cmp( VECTOR(left->v), VECTOR(right->v), (count_t) size_right );
}


igraph_bool_t igraph_biguint_equal(igraph_biguint_t *left, igraph_biguint_t *right) {
    return 0 == igraph_biguint_compare(left, right);
}


igraph_bool_t igraph_biguint_bigger(igraph_biguint_t *left,
                                    igraph_biguint_t *right) {
    return 0 < igraph_biguint_compare(left, right);
}


igraph_bool_t igraph_biguint_biggerorequal(igraph_biguint_t *left,
        igraph_biguint_t *right) {
    return 0 <= igraph_biguint_compare(left, right);
}

int igraph_biguint_inc(igraph_biguint_t *res, igraph_biguint_t *b) {
    return igraph_biguint_add_limb(res, b, 1);
}

int igraph_biguint_dec(igraph_biguint_t *res, igraph_biguint_t *b) {
    return igraph_biguint_sub_limb(res, b, 1);
}


int igraph_biguint_add_limb(igraph_biguint_t *res, igraph_biguint_t *b,
                            limb_t l) {
    int nlimb = igraph_biguint_size(b);
    limb_t carry;

    if (res != b) {
        IGRAPH_CHECK(igraph_biguint_resize(res, nlimb));
    }

    carry = bn_add_limb( VECTOR(res->v), VECTOR(b->v), l, (count_t) nlimb);
    if (carry) {
        IGRAPH_CHECK(igraph_biguint_extend(res, carry));
    }
    return 0;
}

int igraph_biguint_sub_limb(igraph_biguint_t *res, igraph_biguint_t *b,
                            limb_t l) {
    int nlimb = igraph_biguint_size(b);

    if (res != b) {
        IGRAPH_CHECK(igraph_biguint_resize(res, nlimb));
    }

    /* We don't check the return value here */
    bn_sub_limb( VECTOR(res->v), VECTOR(b->v), l, (count_t) nlimb);

    return 0;
}

int igraph_biguint_mul_limb(igraph_biguint_t *res, igraph_biguint_t *b,
                            limb_t l) {
    int nlimb = igraph_biguint_size(b);
    limb_t carry;

    if (res != b) {
        IGRAPH_CHECK(igraph_biguint_resize(res, nlimb));
    }

    carry = bn_mul_limb( VECTOR(res->v), VECTOR(b->v), l, (count_t) nlimb);
    if (carry) {
        IGRAPH_CHECK(igraph_biguint_extend(res, carry));
    }
    return 0;
}

int igraph_biguint_add(igraph_biguint_t *res, igraph_biguint_t *left,
                       igraph_biguint_t *right) {

    int size_left = igraph_biguint_size(left);
    int size_right = igraph_biguint_size(right);
    limb_t carry;

    if (size_left > size_right) {
        IGRAPH_CHECK(igraph_biguint_resize(right, size_left));
        size_right = size_left;
    } else if (size_left < size_right) {
        IGRAPH_CHECK(igraph_biguint_resize(left, size_right));
        size_left = size_right;
    }
    IGRAPH_CHECK(igraph_biguint_resize(res, size_left));

    carry = bn_add( VECTOR(res->v), VECTOR(left->v), VECTOR(right->v),
                    (count_t) size_left);
    if (carry) {
        IGRAPH_CHECK(igraph_biguint_extend(res, carry));
    }
    return 0;
}

int igraph_biguint_sub(igraph_biguint_t *res, igraph_biguint_t *left,
                       igraph_biguint_t *right) {

    int size_left = igraph_biguint_size(left);
    int size_right = igraph_biguint_size(right);

    if (size_left > size_right) {
        IGRAPH_CHECK(igraph_biguint_resize(right, size_left));
        size_right = size_left;
    } else if (size_left < size_right) {
        IGRAPH_CHECK(igraph_biguint_resize(left, size_right));
        size_left = size_right;
    }
    IGRAPH_CHECK(igraph_biguint_resize(res, size_left));

    /* We don't check return value, left should not be smaller than right! */
    bn_sub( VECTOR(res->v), VECTOR(left->v), VECTOR(right->v),
            (count_t) size_left);

    return 0;
}

int igraph_biguint_mul(igraph_biguint_t *res, igraph_biguint_t *left,
                       igraph_biguint_t *right) {

    int size_left = igraph_biguint_size(left);
    int size_right = igraph_biguint_size(right);

    if (size_left > size_right) {
        IGRAPH_CHECK(igraph_biguint_resize(right, size_left));
        size_right = size_left;
    } else if (size_left < size_right) {
        IGRAPH_CHECK(igraph_biguint_resize(left, size_right));
        size_left = size_right;
    }
    IGRAPH_CHECK(igraph_biguint_resize(res, 2 * size_left));

    bn_mul( VECTOR(res->v), VECTOR(left->v), VECTOR(right->v),
            (count_t) size_left );
    return 0;
}

int igraph_biguint_div(igraph_biguint_t *q, igraph_biguint_t *r,
                       igraph_biguint_t *u, igraph_biguint_t *v) {

    int ret;
    int size_q = igraph_biguint_size(q);
    int size_r = igraph_biguint_size(r);
    int size_u = igraph_biguint_size(u);
    int size_v = igraph_biguint_size(v);
    int size_qru = size_q > size_r ? size_q : size_r;
    size_qru = size_u > size_qru ? size_u : size_qru;

    if (size_q < size_qru) {
        IGRAPH_CHECK(igraph_biguint_resize(q, size_qru));
    }
    if (size_r < size_qru) {
        IGRAPH_CHECK(igraph_biguint_resize(r, size_qru));
    }
    if (size_u < size_qru) {
        IGRAPH_CHECK(igraph_biguint_resize(u, size_qru));
    }

    ret = bn_div( VECTOR(q->v), VECTOR(r->v), VECTOR(u->v), VECTOR(v->v),
                  (count_t) size_qru, (count_t) size_v );

    if (ret) {
        IGRAPH_ERROR("Bigint division by zero", IGRAPH_EDIVZERO);
    }

    return 0;
}

#ifndef USING_R
int igraph_biguint_print(igraph_biguint_t *b) {
    return igraph_biguint_fprint(b, stdout);
}
#endif

int igraph_biguint_fprint(igraph_biguint_t *b, FILE *file) {

    /* It is hard to control memory allocation for the bn2d function,
       so we do our own version */

    int n = igraph_biguint_size(b);
    long int size = 12 * n + 1;
    igraph_biguint_t tmp;
    char *dst;
    limb_t r;

    /* Zero? */
    if (!bn_cmp_limb(VECTOR(b->v), 0, (count_t) n)) {
        fputs("0", file);
        return 0;
    }

    IGRAPH_CHECK(igraph_biguint_copy(&tmp, b));
    IGRAPH_FINALLY(igraph_biguint_destroy, &tmp);
    dst = igraph_Calloc(size, char);
    if (!dst) {
        IGRAPH_ERROR("Cannot print big number", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, dst);

    size--;
    dst[size] = '\0';
    while (0 != bn_cmp_limb(VECTOR(tmp.v), 0, (count_t) n)) {
        r = bn_div_limb(VECTOR(tmp.v), VECTOR(tmp.v), 10, (count_t) n);
        dst[--size] = '0' + (char) r;
    }

    fputs(&dst[size], file);

    igraph_Free(dst);
    igraph_biguint_destroy(&tmp);
    IGRAPH_FINALLY_CLEAN(2);

    return 0;
}

