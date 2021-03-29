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
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include "igraph_error.h"
#include "igraph_types.h"
#include "igraph_vector.h"
#include "igraph_complex.h"

#include <float.h>

#define BASE_IGRAPH_REAL
#include "igraph_pmt.h"
#include "vector.pmt"
#include "igraph_pmt_off.h"
#undef BASE_IGRAPH_REAL

#define BASE_FLOAT
#include "igraph_pmt.h"
#include "vector.pmt"
#include "igraph_pmt_off.h"
#undef BASE_FLOAT

#define BASE_LONG
#include "igraph_pmt.h"
#include "vector.pmt"
#include "igraph_pmt_off.h"
#undef BASE_LONG

#define BASE_CHAR
#include "igraph_pmt.h"
#include "vector.pmt"
#include "igraph_pmt_off.h"
#undef BASE_CHAR

#define BASE_BOOL
#include "igraph_pmt.h"
#include "vector.pmt"
#include "igraph_pmt_off.h"
#undef BASE_BOOL

#define BASE_INT
#include "igraph_pmt.h"
#include "vector.pmt"
#include "igraph_pmt_off.h"
#undef BASE_INT

#define BASE_COMPLEX
#include "igraph_pmt.h"
#include "vector.pmt"
#include "igraph_pmt_off.h"
#undef BASE_COMPLEX

#include "core/math.h"

#include "core/indheap.h"

int igraph_vector_floor(const igraph_vector_t *from, igraph_vector_long_t *to) {
    long int i, n = igraph_vector_size(from);

    IGRAPH_CHECK(igraph_vector_long_resize(to, n));
    for (i = 0; i < n; i++) {
        VECTOR(*to)[i] = (long int) floor(VECTOR(*from)[i]);
    }
    return 0;
}

int igraph_vector_round(const igraph_vector_t *from, igraph_vector_long_t *to) {
    long int i, n = igraph_vector_size(from);

    IGRAPH_CHECK(igraph_vector_long_resize(to, n));
    for (i = 0; i < n; i++) {
        VECTOR(*to)[i] = (long int) round(VECTOR(*from)[i]);
    }
    return 0;
}

int igraph_vector_order2(igraph_vector_t *v) {

    igraph_indheap_t heap;

    igraph_indheap_init_array(&heap, VECTOR(*v), igraph_vector_size(v));
    IGRAPH_FINALLY(igraph_indheap_destroy, &heap);

    igraph_vector_clear(v);
    while (!igraph_indheap_empty(&heap)) {
        IGRAPH_CHECK(igraph_vector_push_back(v, igraph_indheap_max_index(&heap) - 1));
        igraph_indheap_delete_max(&heap);
    }

    igraph_indheap_destroy(&heap);
    IGRAPH_FINALLY_CLEAN(1);
    return 0;
}

/**
 * \ingroup vector
 * \function igraph_vector_order
 * \brief Calculate the order of the elements in a vector.
 *
 * </para><para>
 * The smallest element will have order zero, the second smallest
 * order one, etc.
 * \param v The original \type igraph_vector_t object.
 * \param v2 A secondary key, another \type igraph_vector_t object.
 * \param res An initialized \type igraph_vector_t object, it will be
 *    resized to match the size of \p v. The
 *    result of the computation will be stored here.
 * \param nodes Hint, the largest element in \p v.
 * \return Error code:
 *         \c IGRAPH_ENOMEM: out of memory
 *
 * Time complexity: O()
 */

int igraph_vector_order(const igraph_vector_t* v,
                        const igraph_vector_t *v2,
                        igraph_vector_t* res, igraph_real_t nodes) {
    long int edges = igraph_vector_size(v);
    igraph_vector_t ptr;
    igraph_vector_t rad;
    long int i, j;

    IGRAPH_ASSERT(v != NULL);
    IGRAPH_ASSERT(v->stor_begin != NULL);

    IGRAPH_VECTOR_INIT_FINALLY(&ptr, (long int) nodes + 1);
    IGRAPH_VECTOR_INIT_FINALLY(&rad, edges);
    IGRAPH_CHECK(igraph_vector_resize(res, edges));

    for (i = 0; i < edges; i++) {
        long int radix = (long int) v2->stor_begin[i];
        if (VECTOR(ptr)[radix] != 0) {
            VECTOR(rad)[i] = VECTOR(ptr)[radix];
        }
        VECTOR(ptr)[radix] = i + 1;
    }

    j = 0;
    for (i = 0; i < nodes + 1; i++) {
        if (VECTOR(ptr)[i] != 0) {
            long int next = (long int) VECTOR(ptr)[i] - 1;
            res->stor_begin[j++] = next;
            while (VECTOR(rad)[next] != 0) {
                next = (long int) VECTOR(rad)[next] - 1;
                res->stor_begin[j++] = next;
            }
        }
    }

    igraph_vector_null(&ptr);
    igraph_vector_null(&rad);

    for (i = 0; i < edges; i++) {
        long int edge = (long int) VECTOR(*res)[edges - i - 1];
        long int radix = (long int) VECTOR(*v)[edge];
        if (VECTOR(ptr)[radix] != 0) {
            VECTOR(rad)[edge] = VECTOR(ptr)[radix];
        }
        VECTOR(ptr)[radix] = edge + 1;
    }

    j = 0;
    for (i = 0; i < nodes + 1; i++) {
        if (VECTOR(ptr)[i] != 0) {
            long int next = (long int) VECTOR(ptr)[i] - 1;
            res->stor_begin[j++] = next;
            while (VECTOR(rad)[next] != 0) {
                next = (long int) VECTOR(rad)[next] - 1;
                res->stor_begin[j++] = next;
            }
        }
    }

    igraph_vector_destroy(&ptr);
    igraph_vector_destroy(&rad);
    IGRAPH_FINALLY_CLEAN(2);

    return 0;
}

int igraph_vector_order1(const igraph_vector_t* v,
                         igraph_vector_t* res, igraph_real_t nodes) {
    long int edges = igraph_vector_size(v);
    igraph_vector_t ptr;
    igraph_vector_t rad;
    long int i, j;

    IGRAPH_ASSERT(v != NULL);
    IGRAPH_ASSERT(v->stor_begin != NULL);

    IGRAPH_VECTOR_INIT_FINALLY(&ptr, (long int) nodes + 1);
    IGRAPH_VECTOR_INIT_FINALLY(&rad, edges);
    IGRAPH_CHECK(igraph_vector_resize(res, edges));

    for (i = 0; i < edges; i++) {
        long int radix = (long int) v->stor_begin[i];
        if (VECTOR(ptr)[radix] != 0) {
            VECTOR(rad)[i] = VECTOR(ptr)[radix];
        }
        VECTOR(ptr)[radix] = i + 1;
    }

    j = 0;
    for (i = 0; i < nodes + 1; i++) {
        if (VECTOR(ptr)[i] != 0) {
            long int next = (long int) VECTOR(ptr)[i] - 1;
            res->stor_begin[j++] = next;
            while (VECTOR(rad)[next] != 0) {
                next = (long int) VECTOR(rad)[next] - 1;
                res->stor_begin[j++] = next;
            }
        }
    }

    igraph_vector_destroy(&ptr);
    igraph_vector_destroy(&rad);
    IGRAPH_FINALLY_CLEAN(2);

    return 0;
}

int igraph_vector_order1_int(const igraph_vector_t* v,
                             igraph_vector_int_t* res,
                             igraph_real_t nodes) {
    long int edges = igraph_vector_size(v);
    igraph_vector_t ptr;
    igraph_vector_t rad;
    long int i, j;

    IGRAPH_ASSERT(v != NULL);
    IGRAPH_ASSERT(v->stor_begin != NULL);

    IGRAPH_VECTOR_INIT_FINALLY(&ptr, (long int) nodes + 1);
    IGRAPH_VECTOR_INIT_FINALLY(&rad, edges);
    IGRAPH_CHECK(igraph_vector_int_resize(res, edges));

    for (i = 0; i < edges; i++) {
        long int radix = (long int) v->stor_begin[i];
        if (VECTOR(ptr)[radix] != 0) {
            VECTOR(rad)[i] = VECTOR(ptr)[radix];
        }
        VECTOR(ptr)[radix] = i + 1;
    }

    j = 0;
    for (i = 0; i < nodes + 1; i++) {
        if (VECTOR(ptr)[i] != 0) {
            long int next = (long int) VECTOR(ptr)[i] - 1;
            res->stor_begin[j++] = next;
            while (VECTOR(rad)[next] != 0) {
                next = (long int) VECTOR(rad)[next] - 1;
                res->stor_begin[j++] = next;
            }
        }
    }

    igraph_vector_destroy(&ptr);
    igraph_vector_destroy(&rad);
    IGRAPH_FINALLY_CLEAN(2);

    return 0;
}

int igraph_vector_rank(const igraph_vector_t *v, igraph_vector_t *res,
                       long int nodes) {

    igraph_vector_t rad;
    igraph_vector_t ptr;
    long int edges = igraph_vector_size(v);
    long int i, c = 0;

    IGRAPH_VECTOR_INIT_FINALLY(&rad, nodes);
    IGRAPH_VECTOR_INIT_FINALLY(&ptr, edges);
    IGRAPH_CHECK(igraph_vector_resize(res, edges));

    for (i = 0; i < edges; i++) {
        long int elem = (long int) VECTOR(*v)[i];
        VECTOR(ptr)[i] = VECTOR(rad)[elem];
        VECTOR(rad)[elem] = i + 1;
    }

    for (i = 0; i < nodes; i++) {
        long int p = (long int) VECTOR(rad)[i];
        while (p != 0) {
            VECTOR(*res)[p - 1] = c++;
            p = (long int) VECTOR(ptr)[p - 1];
        }
    }

    igraph_vector_destroy(&ptr);
    igraph_vector_destroy(&rad);
    IGRAPH_FINALLY_CLEAN(2);
    return 0;
}

#ifndef USING_R
int igraph_vector_complex_print(const igraph_vector_complex_t *v) {
    long int i, n = igraph_vector_complex_size(v);
    if (n != 0) {
        igraph_complex_t z = VECTOR(*v)[0];
        printf("%g%+gi", IGRAPH_REAL(z), IGRAPH_IMAG(z));
    }
    for (i = 1; i < n; i++) {
        igraph_complex_t z = VECTOR(*v)[i];
        printf(" %g%+gi", IGRAPH_REAL(z), IGRAPH_IMAG(z));
    }
    printf("\n");
    return 0;
}
#endif

int igraph_vector_complex_fprint(const igraph_vector_complex_t *v,
                                 FILE *file) {
    long int i, n = igraph_vector_complex_size(v);
    if (n != 0) {
        igraph_complex_t z = VECTOR(*v)[0];
        fprintf(file, "%g%+g", IGRAPH_REAL(z), IGRAPH_IMAG(z));
    }
    for (i = 1; i < n; i++) {
        igraph_complex_t z = VECTOR(*v)[i];
        fprintf(file, " %g%+g", IGRAPH_REAL(z), IGRAPH_IMAG(z));
    }
    fprintf(file, "\n");
    return 0;
}

int igraph_vector_complex_real(const igraph_vector_complex_t *v,
                               igraph_vector_t *real) {
    long int i, n = igraph_vector_complex_size(v);
    IGRAPH_CHECK(igraph_vector_resize(real, n));
    for (i = 0; i < n; i++) {
        VECTOR(*real)[i] = IGRAPH_REAL(VECTOR(*v)[i]);
    }

    return 0;
}

int igraph_vector_complex_imag(const igraph_vector_complex_t *v,
                               igraph_vector_t *imag) {
    long int i, n = igraph_vector_complex_size(v);
    IGRAPH_CHECK(igraph_vector_resize(imag, n));
    for (i = 0; i < n; i++) {
        VECTOR(*imag)[i] = IGRAPH_IMAG(VECTOR(*v)[i]);
    }

    return 0;
}

int igraph_vector_complex_realimag(const igraph_vector_complex_t *v,
                                   igraph_vector_t *real,
                                   igraph_vector_t *imag) {
    long int i, n = igraph_vector_complex_size(v);
    IGRAPH_CHECK(igraph_vector_resize(real, n));
    IGRAPH_CHECK(igraph_vector_resize(imag, n));
    for (i = 0; i < n; i++) {
        igraph_complex_t z = VECTOR(*v)[i];
        VECTOR(*real)[i] = IGRAPH_REAL(z);
        VECTOR(*imag)[i] = IGRAPH_IMAG(z);
    }

    return 0;
}

int igraph_vector_complex_create(igraph_vector_complex_t *v,
                                 const igraph_vector_t *real,
                                 const igraph_vector_t *imag) {
    long int i, n = igraph_vector_size(real);
    if (n != igraph_vector_size(imag)) {
        IGRAPH_ERROR("Real and imag vector sizes don't match", IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_vector_complex_init(v, n));
    /* FINALLY not needed */

    for (i = 0; i < n; i++) {
        VECTOR(*v)[i] = igraph_complex(VECTOR(*real)[i], VECTOR(*imag)[i]);
    }

    return 0;
}

int igraph_vector_complex_create_polar(igraph_vector_complex_t *v,
                                       const igraph_vector_t *r,
                                       const igraph_vector_t *theta) {
    long int i, n = igraph_vector_size(r);
    if (n != igraph_vector_size(theta)) {
        IGRAPH_ERROR("'r' and 'theta' vector sizes don't match", IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_vector_complex_init(v, n));
    /* FINALLY not needed */

    for (i = 0; i < n; i++) {
        VECTOR(*v)[i] = igraph_complex_polar(VECTOR(*r)[i], VECTOR(*theta)[i]);
    }

    return 0;
}

igraph_bool_t igraph_vector_e_tol(const igraph_vector_t *lhs,
                                  const igraph_vector_t *rhs,
                                  igraph_real_t tol) {
    long int i, s;
    IGRAPH_ASSERT(lhs != 0);
    IGRAPH_ASSERT(rhs != 0);
    IGRAPH_ASSERT(lhs->stor_begin != 0);
    IGRAPH_ASSERT(rhs->stor_begin != 0);

    s = igraph_vector_size(lhs);
    if (s != igraph_vector_size(rhs)) {
        return 0;
    } else {
        if (tol == 0) {
            tol = DBL_EPSILON;
        }
        for (i = 0; i < s; i++) {
            igraph_real_t l = VECTOR(*lhs)[i];
            igraph_real_t r = VECTOR(*rhs)[i];
            if (l < r - tol || l > r + tol) {
                return 0;
            }
        }
        return 1;
    }
}

int igraph_vector_zapsmall(igraph_vector_t *v, igraph_real_t tol) {
    long int i, n = igraph_vector_size(v);
    if (tol < 0.0) {
        IGRAPH_ERROR("`tol' tolerance must be non-negative", IGRAPH_EINVAL);
    }
    if (tol == 0.0) {
        tol = sqrt(DBL_EPSILON);
    }
    for (i = 0; i < n; i++) {
        igraph_real_t val = VECTOR(*v)[i];
        if (val < tol && val > -tol) {
            VECTOR(*v)[i] = 0.0;
        }
    }
    return 0;
}

/**
 * \ingroup vector
 * \function igraph_vector_is_nan
 * \brief Check for each element if it is NaN.
 *
 * </para><para>
 * \param v The \type igraph_vector_t object to check.
 * \param is_nan The resulting boolean vector indicating for each element
 *               whether it is NaN or not.
 * \return Error code,
 *         \c IGRAPH_ENOMEM if there is not enough
 *         memory. Note that this function \em never returns an error
 *         if the vector \p is_nan will already be large enough.
 * Time complexity: O(n), the number of elements.
 */
int igraph_vector_is_nan(const igraph_vector_t *v, igraph_vector_bool_t *is_nan)
{
    igraph_real_t *ptr;
    igraph_bool_t *ptr_nan;
    IGRAPH_ASSERT(v != NULL);
    IGRAPH_ASSERT(v->stor_begin != NULL);
    IGRAPH_ASSERT(is_nan != NULL);
    IGRAPH_ASSERT(is_nan->stor_begin != NULL);
    IGRAPH_CHECK(igraph_vector_bool_resize(is_nan, igraph_vector_size(v)));
    for (ptr = v->stor_begin, ptr_nan = is_nan->stor_begin; ptr < v->end; ptr++, ptr_nan++) {
        *ptr_nan = igraph_is_nan(*ptr);
    }
    return IGRAPH_SUCCESS;
}

/**
 * \ingroup vector
 * \function igraph_vector_is_any_nan
 * \brief Check if any element is NaN.
 *
 * </para><para>
 * \param v The \type igraph_vector_t object to check.
 * \return 1 if any element is NaN, 0 otherwise.
 *
 * Time complexity: O(n), the number of elements.
 */
igraph_bool_t igraph_vector_is_any_nan(const igraph_vector_t *v)
{
    igraph_real_t *ptr;
    IGRAPH_ASSERT(v != NULL);
    IGRAPH_ASSERT(v->stor_begin != NULL);
    ptr = v->stor_begin;
    while (ptr < v->end) {
        if (igraph_is_nan(*ptr)) {
            return 1;
        }
        ptr++;
    }
    return 0;
}
