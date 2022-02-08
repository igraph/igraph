/*
   IGraph library.
   Copyright (C) 2020 The igraph development team

   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received _safe_a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA
*/

#ifndef IGRAPH_MATH_SAFE_INTOP_H
#define IGRAPH_MATH_SAFE_INTOP_H

#include "igraph_decls.h"
#include "igraph_error.h"
#include "igraph_types.h"
#include "igraph_vector.h"

__BEGIN_DECLS

/* These macros raise an error if the operation would result in an overflow.
 * They must only be used in functions that return an igraph_error_t.
 *
 * This code is based on the recommendation of
 * https://wiki.sei.cmu.edu/confluence/display/c/SEI+CERT+C+Coding+Standard
 */

/* TODO: re-enable implementations in terms of GCC intrinsics */
#if /* defined(__GNUC__) */ 0

#define IGRAPH_SAFE_ADD(a, b, res) \
    do { \
        igraph_integer_t _safe_a = (a), _safe_b = (b); \
        igraph_integer_t _safe_sum; \
        if (__builtin_add_overflow(_safe_a, _safe_b, &_safe_sum)) { \
            IGRAPH_ERRORF("Overflow when adding %" IGRAPH_PRId " and %" IGRAPH_PRId ".", IGRAPH_EOVERFLOW, _safe_a, _safe_b); \
        } \
        *(res) = _safe_sum; \
    } while (0)

#define IGRAPH_SAFE_MULT(a, b, res) \
    do { \
        igraph_integer_t _safe_a = (a), _safe_b = (b); \
        igraph_integer_t _safe_prod; \
        if (__builtin_mul_overflow(_safe_a, _safe_b, &_safe_prod)) { \
            IGRAPH_ERRORF("Overflow when multiplying %" IGRAPH_PRId " and %" IGRAPH_PRId ".", IGRAPH_EOVERFLOW, _safe_a, _safe_b); \
        } \
        *(res) = _safe_prod; \
    } while (0)

#else

#define IGRAPH_SAFE_ADD(a, b, res) \
    do { \
        igraph_integer_t _safe_a = (a), _safe_b = (b); \
        igraph_integer_t _safe_sum; \
        if (((_safe_b > 0) && (_safe_a > (IGRAPH_INTEGER_MAX - _safe_b))) || \
            ((_safe_b < 0) && (_safe_a < (IGRAPH_INTEGER_MIN - _safe_b)))) { \
            IGRAPH_ERRORF("Overflow when adding %" IGRAPH_PRId " and %" IGRAPH_PRId ".", IGRAPH_EOVERFLOW, _safe_a, _safe_b); \
        } \
        _safe_sum = _safe_a+_safe_b; \
        *(res) = _safe_sum; \
    } while (0)

#define IGRAPH_SAFE_MULT(a, b, res) \
    do { \
        igraph_integer_t _safe_a = (a), _safe_b = (b); \
        igraph_integer_t _safe_prod; \
        int err=0; \
        if (_safe_a > 0) {  /* _safe_a is positive */ \
            if (_safe_b > 0) {  /* _safe_a and _safe_b are positive */ \
                if (_safe_a > (IGRAPH_INTEGER_MAX / _safe_b)) { \
                    err=1; \
                } \
            } else { /* _safe_a positive, _safe_b nonpositive */ \
                if (_safe_b < (IGRAPH_INTEGER_MIN / _safe_a)) { \
                    err=1; \
                } \
            } /* _safe_a positive, _safe_b nonpositive */ \
        } else { /* _safe_a is nonpositive */ \
            if (_safe_b > 0) { /* _safe_a is nonpositive, _safe_b is positive */ \
                if (_safe_a < (IGRAPH_INTEGER_MIN / _safe_b)) { \
                    err=1; \
                } \
            } else { /* _safe_a and _safe_b are nonpositive */ \
                if ( (_safe_a != 0) && (_safe_b < (IGRAPH_INTEGER_MAX / _safe_a))) { \
                    err=1; \
                } \
            } /* End if _safe_a and _safe_b are nonpositive */ \
        } /* End if _safe_a is nonpositive */ \
        if (err) { \
            IGRAPH_ERRORF("Overflow when multiplying %" IGRAPH_PRId " and %" IGRAPH_PRId ".", IGRAPH_EOVERFLOW, _safe_a, _safe_b); \
        } \
        _safe_prod = _safe_a*_safe_b; \
        *(res) = _safe_prod; \
    } while (0)

#endif

IGRAPH_PRIVATE_EXPORT igraph_error_t igraph_i_safe_add(igraph_integer_t a, igraph_integer_t b, igraph_integer_t *res);
IGRAPH_PRIVATE_EXPORT igraph_error_t igraph_i_safe_mult(igraph_integer_t a, igraph_integer_t b, igraph_integer_t *res);
igraph_error_t igraph_i_safe_vector_int_sum(const igraph_vector_int_t *vec, igraph_integer_t *res);
igraph_error_t igraph_i_safe_vector_int_prod(const igraph_vector_int_t *vec, igraph_integer_t *res);

__END_DECLS

#endif /* IGRAPH_MATH_SAFE_INTOP_H */
