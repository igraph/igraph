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

#include "igraph_vector.h"

#include "igraph_complex.h"
#include "igraph_types.h"
#include "igraph_nongraph.h"

#include <float.h>

#define BASE_IGRAPH_REAL
#include "igraph_pmt.h"
#include "vector.pmt"
#include "igraph_pmt_off.h"
#undef BASE_IGRAPH_REAL

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

#include "core/indheap.h"

/**
 * \ingroup vector
 * \function igraph_vector_floor
 * \brief Transform a real vector to an integer vector by flooring each element.
 *
 * Flooring means rounding down to the nearest integer.
 *
 * \param from The original real vector object.
 * \param to Pointer to an initialized integer vector. The result will be stored here.
 * \return Error code:
 *         \c IGRAPH_ENOMEM: out of memory
 *
 * Time complexity: O(n), where n is the number of elements in the vector.
 */
igraph_error_t igraph_vector_floor(const igraph_vector_t *from, igraph_vector_int_t *to) {
    const igraph_integer_t n = igraph_vector_size(from);

    IGRAPH_CHECK(igraph_vector_int_resize(to, n));
    for (igraph_integer_t i = 0; i < n; i++) {
        VECTOR(*to)[i] = floor(VECTOR(*from)[i]);
    }

    return IGRAPH_SUCCESS;
}

igraph_error_t igraph_vector_round(const igraph_vector_t *from, igraph_vector_int_t *to) {
    const igraph_integer_t n = igraph_vector_size(from);

    IGRAPH_CHECK(igraph_vector_int_resize(to, n));
    for (igraph_integer_t i = 0; i < n; i++) {
        VECTOR(*to)[i] = round(VECTOR(*from)[i]);
    }

    return IGRAPH_SUCCESS;
}

igraph_error_t igraph_vector_order2(igraph_vector_t *v) {
    igraph_indheap_t heap;

    IGRAPH_CHECK(igraph_indheap_init_array(&heap, VECTOR(*v), igraph_vector_size(v)));
    IGRAPH_FINALLY(igraph_indheap_destroy, &heap);

    igraph_vector_clear(v);
    while (!igraph_indheap_empty(&heap)) {
        IGRAPH_CHECK(igraph_vector_push_back(v, igraph_indheap_max_index(&heap) - 1));
        igraph_indheap_delete_max(&heap);
    }

    igraph_indheap_destroy(&heap);
    IGRAPH_FINALLY_CLEAN(1);
    return IGRAPH_SUCCESS;
}

/**
 * \ingroup vector
 * \function igraph_vector_int_pair_order
 * \brief Calculates the order of the elements in a pair of integer vectors of equal length.
 *
 * The smallest element will have order zero, the second smallest
 * order one, etc.
 *
 * \param v The original \ref igraph_vector_int_t object.
 * \param v2 A secondary key, another \ref igraph_vector_int_t object.
 * \param res An initialized \ref igraph_vector_int_t object, it will be
 *    resized to match the size of \p v. The result of the computation will
 *    be stored here.
 * \param nodes Hint, the largest element in \p v.
 * \return Error code:
 *         \c IGRAPH_ENOMEM: out of memory
 *
 * Time complexity: O()
 */

igraph_error_t igraph_vector_int_pair_order(const igraph_vector_int_t* v,
                                       const igraph_vector_int_t* v2,
                                       igraph_vector_int_t* res, igraph_integer_t nodes) {
    igraph_integer_t edges = igraph_vector_int_size(v);
    igraph_vector_int_t ptr;
    igraph_vector_int_t rad;
    igraph_integer_t i, j;

    IGRAPH_ASSERT(v != NULL);
    IGRAPH_ASSERT(v->stor_begin != NULL);

    IGRAPH_VECTOR_INT_INIT_FINALLY(&ptr, nodes + 1);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&rad, edges);
    IGRAPH_CHECK(igraph_vector_int_resize(res, edges));

    for (i = 0; i < edges; i++) {
        igraph_integer_t radix = VECTOR(*v2)[i];
        if (VECTOR(ptr)[radix] != 0) {
            VECTOR(rad)[i] = VECTOR(ptr)[radix];
        }
        VECTOR(ptr)[radix] = i + 1;
    }

    j = 0;
    for (i = 0; i < nodes + 1; i++) {
        if (VECTOR(ptr)[i] != 0) {
            igraph_integer_t next = VECTOR(ptr)[i] - 1;
            VECTOR(*res)[j++] = next;
            while (VECTOR(rad)[next] != 0) {
                next = VECTOR(rad)[next] - 1;
                VECTOR(*res)[j++] = next;
            }
        }
    }

    igraph_vector_int_null(&ptr);
    igraph_vector_int_null(&rad);

    for (i = 0; i < edges; i++) {
        igraph_integer_t edge = VECTOR(*res)[edges - i - 1];
        igraph_integer_t radix = VECTOR(*v)[edge];
        if (VECTOR(ptr)[radix] != 0) {
            VECTOR(rad)[edge] = VECTOR(ptr)[radix];
        }
        VECTOR(ptr)[radix] = edge + 1;
    }

    j = 0;
    for (i = 0; i < nodes + 1; i++) {
        if (VECTOR(ptr)[i] != 0) {
            igraph_integer_t next = VECTOR(ptr)[i] - 1;
            VECTOR(*res)[j++] = next;
            while (VECTOR(rad)[next] != 0) {
                next = VECTOR(rad)[next] - 1;
                VECTOR(*res)[j++] = next;
            }
        }
    }

    igraph_vector_int_destroy(&ptr);
    igraph_vector_int_destroy(&rad);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}

igraph_error_t igraph_vector_int_order1(const igraph_vector_int_t* v,
                             igraph_vector_int_t* res,
                             igraph_integer_t nodes) {
    igraph_integer_t edges = igraph_vector_int_size(v);
    igraph_vector_int_t ptr;
    igraph_vector_int_t rad;
    igraph_integer_t i, j;

    IGRAPH_ASSERT(v != NULL);
    IGRAPH_ASSERT(v->stor_begin != NULL);

    IGRAPH_VECTOR_INT_INIT_FINALLY(&ptr, nodes + 1);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&rad, edges);
    IGRAPH_CHECK(igraph_vector_int_resize(res, edges));

    for (i = 0; i < edges; i++) {
        igraph_integer_t radix = v->stor_begin[i];
        if (VECTOR(ptr)[radix] != 0) {
            VECTOR(rad)[i] = VECTOR(ptr)[radix];
        }
        VECTOR(ptr)[radix] = i + 1;
    }

    j = 0;
    for (i = 0; i < nodes + 1; i++) {
        if (VECTOR(ptr)[i] != 0) {
            igraph_integer_t next = VECTOR(ptr)[i] - 1;
            res->stor_begin[j++] = next;
            while (VECTOR(rad)[next] != 0) {
                next = VECTOR(rad)[next] - 1;
                res->stor_begin[j++] = next;
            }
        }
    }

    igraph_vector_int_destroy(&ptr);
    igraph_vector_int_destroy(&rad);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}

igraph_error_t igraph_vector_rank(
        const igraph_vector_t *v, igraph_vector_int_t *res, igraph_integer_t nodes) {

    igraph_vector_int_t rad;
    igraph_vector_int_t ptr;
    igraph_integer_t edges = igraph_vector_size(v);
    igraph_integer_t i, c = 0;

    IGRAPH_VECTOR_INT_INIT_FINALLY(&rad, nodes);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&ptr, edges);
    IGRAPH_CHECK(igraph_vector_int_resize(res, edges));

    for (i = 0; i < edges; i++) {
        igraph_integer_t elem = VECTOR(*v)[i];
        VECTOR(ptr)[i] = VECTOR(rad)[elem];
        VECTOR(rad)[elem] = i + 1;
    }

    for (i = 0; i < nodes; i++) {
        igraph_integer_t p = VECTOR(rad)[i];
        while (p != 0) {
            VECTOR(*res)[p - 1] = c++;
            p = VECTOR(ptr)[p - 1];
        }
    }

    igraph_vector_int_destroy(&ptr);
    igraph_vector_int_destroy(&rad);
    IGRAPH_FINALLY_CLEAN(2);
    return IGRAPH_SUCCESS;
}

igraph_error_t igraph_vector_int_rank(
        const igraph_vector_int_t *v, igraph_vector_int_t *res, igraph_integer_t nodes) {

    igraph_vector_int_t rad;
    igraph_vector_int_t ptr;
    igraph_integer_t edges = igraph_vector_int_size(v);
    igraph_integer_t i, c = 0;

    IGRAPH_VECTOR_INT_INIT_FINALLY(&rad, nodes);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&ptr, edges);
    IGRAPH_CHECK(igraph_vector_int_resize(res, edges));

    for (i = 0; i < edges; i++) {
        igraph_integer_t elem = VECTOR(*v)[i];
        VECTOR(ptr)[i] = VECTOR(rad)[elem];
        VECTOR(rad)[elem] = i + 1;
    }

    for (i = 0; i < nodes; i++) {
        igraph_integer_t p = VECTOR(rad)[i];
        while (p != 0) {
            VECTOR(*res)[p - 1] = c++;
            p = VECTOR(ptr)[p - 1];
        }
    }

    igraph_vector_int_destroy(&ptr);
    igraph_vector_int_destroy(&rad);
    IGRAPH_FINALLY_CLEAN(2);
    return IGRAPH_SUCCESS;
}

/**
 * \ingroup vector
 * \function igraph_vector_complex_real
 * \brief Gives the real part of a complex vector.
 *
 * \param v Pointer to a complex vector.
 * \param real Pointer to an initialized vector. The result will be stored here.
 * \return Error code.
 *
 * Time complexity: O(n), n is the number of elements in the vector.
 */

igraph_error_t igraph_vector_complex_real(const igraph_vector_complex_t *v,
                               igraph_vector_t *real) {
    igraph_integer_t i, n = igraph_vector_complex_size(v);
    IGRAPH_CHECK(igraph_vector_resize(real, n));
    for (i = 0; i < n; i++) {
        VECTOR(*real)[i] = IGRAPH_REAL(VECTOR(*v)[i]);
    }

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup vector
 * \function igraph_vector_complex_imag
 * \brief Gives the imaginary part of a complex vector.
 *
 * \param v Pointer to a complex vector.
 * \param real Pointer to an initialized vector. The result will be stored here.
 * \return Error code.
 *
 * Time complexity: O(n), n is the number of elements in the vector.
 */

igraph_error_t igraph_vector_complex_imag(const igraph_vector_complex_t *v,
                               igraph_vector_t *imag) {
    igraph_integer_t i, n = igraph_vector_complex_size(v);
    IGRAPH_CHECK(igraph_vector_resize(imag, n));
    for (i = 0; i < n; i++) {
        VECTOR(*imag)[i] = IGRAPH_IMAG(VECTOR(*v)[i]);
    }

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup vector
 * \function igraph_vector_complex_realimag
 * \brief Gives the real and imaginary parts of a complex vector.
 *
 * \param v Pointer to a complex vector.
 * \param real Pointer to an initialized vector. The real part will be stored here.
 * \param imag Pointer to an initialized vector. The imaginary part will be stored here.
 * \return Error code.
 *
 * Time complexity: O(n), n is the number of elements in the vector.
 */

igraph_error_t igraph_vector_complex_realimag(const igraph_vector_complex_t *v,
                                   igraph_vector_t *real,
                                   igraph_vector_t *imag) {
    igraph_integer_t i, n = igraph_vector_complex_size(v);
    IGRAPH_CHECK(igraph_vector_resize(real, n));
    IGRAPH_CHECK(igraph_vector_resize(imag, n));
    for (i = 0; i < n; i++) {
        igraph_complex_t z = VECTOR(*v)[i];
        VECTOR(*real)[i] = IGRAPH_REAL(z);
        VECTOR(*imag)[i] = IGRAPH_IMAG(z);
    }

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup vector
 * \function igraph_vector_complex_create
 * \brief Creates a complex vector from a real and imaginary part.
 *
 * \param v Pointer to an uninitialized complex vector.
 * \param real Pointer to the real part of the complex vector.
 * \param imag Pointer to the imaginary part of the complex vector.
 * \return Error code.
 *
 * Time complexity: O(n), n is the number of elements in the vector.
 */

igraph_error_t igraph_vector_complex_create(igraph_vector_complex_t *v,
                                 const igraph_vector_t *real,
                                 const igraph_vector_t *imag) {
    igraph_integer_t i, n = igraph_vector_size(real);
    if (n != igraph_vector_size(imag)) {
        IGRAPH_ERROR("Real and imag vector sizes don't match", IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_vector_complex_init(v, n));
    /* FINALLY not needed */

    for (i = 0; i < n; i++) {
        VECTOR(*v)[i] = igraph_complex(VECTOR(*real)[i], VECTOR(*imag)[i]);
    }

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup vector
 * \function igraph_vector_complex_create_polar
 * \brief Creates a complex matrix from a magnitude and an angle.
 *
 * \param m Pointer to an uninitialized complex vector.
 * \param r Pointer to a real vector containing magnitudes.
 * \param theta Pointer to a real vector containing arguments (phase angles).
 * \return Error code.
 *
 * Time complexity: O(n), n is the number of elements in the vector.
 */

igraph_error_t igraph_vector_complex_create_polar(igraph_vector_complex_t *v,
                                       const igraph_vector_t *r,
                                       const igraph_vector_t *theta) {
    igraph_integer_t i, n = igraph_vector_size(r);
    if (n != igraph_vector_size(theta)) {
        IGRAPH_ERROR("'r' and 'theta' vector sizes don't match", IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_vector_complex_init(v, n));
    /* FINALLY not needed */

    for (i = 0; i < n; i++) {
        VECTOR(*v)[i] = igraph_complex_polar(VECTOR(*r)[i], VECTOR(*theta)[i]);
    }

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_vector_complex_all_almost_e
 * \brief Are all elements almost equal?
 *
 * Checks if the elements of two complex vectors are equal within a relative tolerance.
 *
 * \param lhs The first vector.
 * \param rhs The second vector.
 * \param eps Relative tolerance, see \ref igraph_complex_almost_equals() for details.
 * \return True if the two vectors are almost equal, false if there is at least
 *     one differing element or if the vectors are not of the same size.
 */
igraph_bool_t igraph_vector_complex_all_almost_e(const igraph_vector_complex_t *lhs,
                                                 const igraph_vector_complex_t *rhs,
                                                 igraph_real_t eps) {

    igraph_integer_t n = igraph_vector_complex_size(lhs);

    if (lhs == rhs) {
        return true;
    }

    if (igraph_vector_complex_size(rhs) != n) {
        return false;
    }

    for (igraph_integer_t i=0; i < n; i++) {
        if (! igraph_complex_almost_equals(VECTOR(*lhs)[i], VECTOR(*rhs)[i], eps))
            return false;
    }

    return true;
}

/**
 * Deprecated in favour of \ref igraph_vector_all_almost_e() which uses
 * relative tolerances. Will be removed in 0.11.
 *
 * Checks if two vectors are equal within an absolute tolerance.
 */
igraph_bool_t igraph_vector_e_tol(const igraph_vector_t *lhs,
                                  const igraph_vector_t *rhs,
                                  igraph_real_t tol) {
    igraph_integer_t i, s;
    IGRAPH_ASSERT(lhs != 0);
    IGRAPH_ASSERT(rhs != 0);
    IGRAPH_ASSERT(lhs->stor_begin != 0);
    IGRAPH_ASSERT(rhs->stor_begin != 0);

    s = igraph_vector_size(lhs);
    if (s != igraph_vector_size(rhs)) {
        return false;
    } else {
        if (tol == 0) {
            tol = DBL_EPSILON;
        }
        for (i = 0; i < s; i++) {
            igraph_real_t l = VECTOR(*lhs)[i];
            igraph_real_t r = VECTOR(*rhs)[i];
            if (l < r - tol || l > r + tol) {
                return false;
            }
        }
        return true;
    }
}

/**
 * \function igraph_vector_all_almost_e
 * \brief Are all elements almost equal?
 *
 * Checks if the elements of two vectors are equal within a relative tolerance.
 *
 * \param lhs The first vector.
 * \param rhs The second vector.
 * \param eps Relative tolerance, see \ref igraph_almost_equals() for details.
 * \return True if the two vectors are almost equal, false if there is at least
 *     one differing element or if the vectors are not of the same size.
 */
igraph_bool_t igraph_vector_all_almost_e(const igraph_vector_t *lhs,
                                         const igraph_vector_t *rhs,
                                         igraph_real_t eps) {

    igraph_integer_t n = igraph_vector_size(lhs);

    if (lhs == rhs) {
        return true;
    }

    if (igraph_vector_size(rhs) != n) {
        return false;
    }

    for (igraph_integer_t i=0; i < n; i++) {
        if (! igraph_almost_equals(VECTOR(*lhs)[i], VECTOR(*rhs)[i], eps))
            return false;
    }

    return true;
}

/**
 * \function igraph_vector_zapsmall
 * \brief Replaces small elements of a vector by exact zeros.
 *
 * Vector elements which are smaller in magnitude than the given absolute
 * tolerance will be replaced by exact zeros. The default tolerance
 * corresponds to two-thirds of the representable digits of \type igraph_real_t,
 * i.e. <code>DBL_EPSILON^(2/3)</code> which is approximately <code>10^-10</code>.
 *
 * \param v   The vector to process, it will be changed in-place.
 * \param tol Tolerance value. Numbers smaller than this in magnitude will
 *            be replaced by zeros. Pass in zero to use the default tolerance.
 *            Must not be negative.
 * \return Error code.
 *
 * \sa \ref igraph_vector_all_almost_e() and \ref igraph_almost_equals() to
 * perform comparisons with relative tolerances.
 */
igraph_error_t igraph_vector_zapsmall(igraph_vector_t *v, igraph_real_t tol) {
    igraph_integer_t i, n = igraph_vector_size(v);
    if (tol < 0.0) {
        IGRAPH_ERROR("Tolerance must be positive or zero.", IGRAPH_EINVAL);
    }
    if (tol == 0.0) {
        tol = pow(DBL_EPSILON, 2.0/3);
    }
    for (i = 0; i < n; i++) {
        igraph_real_t val = VECTOR(*v)[i];
        if (val < tol && val > -tol) {
            VECTOR(*v)[i] = 0.0;
        }
    }
    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_vector_complex_zapsmall
 * \brief Replaces small elements of a complex vector by exact zeros.
 *
 * Similarly to \ref igraph_vector_zapsmall(), small elements will be replaced
 * by zeros. The operation is performed separately on the real and imaginary
 * parts of the numbers. This way, complex numbers with a large real part and
 * tiny imaginary part will effectively be transformed to real numbers.
 * The default tolerance
 * corresponds to two-thirds of the representable digits of \type igraph_real_t,
 * i.e. <code>DBL_EPSILON^(2/3)</code> which is approximately <code>10^-10</code>.
 *
 * \param v   The vector to process, it will be changed in-place.
 * \param tol Tolerance value. Real and imaginary parts smaller than this in
 *            magnitude will be replaced by zeros. Pass in zero to use the default
 *            tolerance. Must not be negative.
 * \return Error code.
 *
 * \sa \ref igraph_vector_complex_all_almost_e() and
 * \ref igraph_complex_almost_equals() to perform comparisons with relative
 * tolerances.
 */
igraph_error_t igraph_vector_complex_zapsmall(igraph_vector_complex_t *v, igraph_real_t tol) {
    igraph_integer_t i, n = igraph_vector_complex_size(v);
    if (tol < 0.0) {
        IGRAPH_ERROR("Tolerance must be positive or zero.", IGRAPH_EINVAL);
    }
    if (tol == 0.0) {
        tol = pow(DBL_EPSILON, 2.0/3);
    }
    for (i = 0; i < n; i++) {
        igraph_complex_t val = VECTOR(*v)[i];
        igraph_bool_t zapped = false;
        if (IGRAPH_REAL(val) < tol && IGRAPH_REAL(val) > -tol) {
            IGRAPH_REAL(val) = 0.0;
            zapped = true;
        }
        if (IGRAPH_IMAG(val) < tol && IGRAPH_IMAG(val) > -tol) {
            IGRAPH_IMAG(val) = 0.0;
            zapped = true;
        }
        if (zapped) {
            VECTOR(*v)[i] = val;
        }
    }
    return IGRAPH_SUCCESS;
}

/**
 * \ingroup vector
 * \function igraph_vector_is_nan
 * \brief Check for each element if it is NaN.
 *
 * \param v The \type igraph_vector_t object to check.
 * \param is_nan The resulting boolean vector indicating for each element
 *               whether it is NaN or not.
 * \return Error code,
 *         \c IGRAPH_ENOMEM if there is not enough memory.
 *         Note that this function \em never returns an error
 *         if the vector \p is_nan will already be large enough.
 *
 * Time complexity: O(n), the number of elements.
 */
igraph_error_t igraph_vector_is_nan(const igraph_vector_t *v, igraph_vector_bool_t *is_nan)
{
    igraph_real_t *ptr;
    igraph_bool_t *ptr_nan;
    IGRAPH_ASSERT(v != NULL);
    IGRAPH_ASSERT(v->stor_begin != NULL);
    IGRAPH_ASSERT(is_nan != NULL);
    IGRAPH_ASSERT(is_nan->stor_begin != NULL);
    IGRAPH_CHECK(igraph_vector_bool_resize(is_nan, igraph_vector_size(v)));
    for (ptr = v->stor_begin, ptr_nan = is_nan->stor_begin; ptr < v->end; ptr++, ptr_nan++) {
        *ptr_nan = isnan(*ptr);
    }
    return IGRAPH_SUCCESS;
}

/**
 * \ingroup vector
 * \function igraph_vector_is_any_nan
 * \brief Check if any element is NaN.
 *
 * \param v The \type igraph_vector_t object to check.
 * \return True if any element is NaN, false otherwise.
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
        if (isnan(*ptr)) {
            return true;
        }
        ptr++;
    }
    return false;
}
