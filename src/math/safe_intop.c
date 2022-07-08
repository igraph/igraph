/*
   IGraph library.
   Copyright (C) 2022 The igraph development team

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

#include "math/safe_intop.h"

/* Use IGRAPH_SAFE_ADD() instead unless there is a need to intercept errors. */
igraph_error_t igraph_i_safe_add(igraph_integer_t a, igraph_integer_t b, igraph_integer_t *res) {
    IGRAPH_SAFE_ADD(a, b, res);
    return IGRAPH_SUCCESS;
}

/* Use IGRAPH_SAFE_MULT() instead unless there is a need to intercept errors. */
igraph_error_t igraph_i_safe_mult(igraph_integer_t a, igraph_integer_t b, igraph_integer_t *res) {
    IGRAPH_SAFE_MULT(a, b, res);
    return IGRAPH_SUCCESS;
}

/* Overflow-safe sum of integer vector elements. */
igraph_error_t igraph_i_safe_vector_int_sum(const igraph_vector_int_t *vec, igraph_integer_t *res) {
    igraph_integer_t i, n = igraph_vector_int_size(vec);
    igraph_integer_t sum = 0;
    for (i=0; i < n; ++i) {
        IGRAPH_SAFE_ADD(sum, VECTOR(*vec)[i], &sum);
    }
    *res = sum;
    return IGRAPH_SUCCESS;
}

/* Overflow-safe product of integer vector elements. */
igraph_error_t igraph_i_safe_vector_int_prod(const igraph_vector_int_t *vec, igraph_integer_t *res) {
    igraph_integer_t i, n = igraph_vector_int_size(vec);
    igraph_integer_t prod = 1;
    for (i=0; i < n; ++i) {
        IGRAPH_SAFE_MULT(prod, VECTOR(*vec)[i], &prod);
    }
    *res = prod;
    return IGRAPH_SUCCESS;
}

/**
 *  Rounds up an integer to the next power of 2, with overflow check.
 *  The result for 2, 3 and 4, respectively, would be 2, 4, and 4.
 *  This function must not be called with negative input.
 *  Based on https://graphics.stanford.edu/~seander/bithacks.html#RoundUpPowerOf2
 */
igraph_error_t igraph_i_safe_next_pow_2(igraph_integer_t k, igraph_integer_t *res) {
    IGRAPH_ASSERT(k >= 0);
    if (k == 0) {
        *res = 0;
        return IGRAPH_SUCCESS;
    }
    k--;
    k |= k >> 1;
    k |= k >> 2;
    k |= k >> 4;
    k |= k >> 8;
    k |= k >> 16;
#if IGRAPH_INTEGER_SIZE == 32
    /* Nothing else to do. */
#elif IGRAPH_INTEGER_SIZE == 64
    k |= k >> 32;
#else
    /* If values other than 32 or 64 become allowed,
     * this code will need to be updated. */
#  error "Unexpected IGRAPH_INTEGER_SIZE value."
#endif
    if (k < IGRAPH_INTEGER_MAX) {
        *res = k+1;
        return IGRAPH_SUCCESS;
    } else {
        IGRAPH_ERRORF("Overflow when computing next power of 2 for %" IGRAPH_PRId ".",
                      IGRAPH_EOVERFLOW, k);
    }
}

/**
 * Computes 2^k as an integer, with overflow check.
 * This function must not be called with negative input.
 */
igraph_error_t igraph_i_safe_exp2(igraph_integer_t k, igraph_integer_t *res) {
    IGRAPH_ASSERT(k >= 0);
    if (k > IGRAPH_INTEGER_SIZE-2) {
        IGRAPH_ERRORF("Overflow when raising 2 to power %" IGRAPH_PRId ".",
                      IGRAPH_EOVERFLOW, k);
    }
    *res = (igraph_integer_t) 1 << k;
    return IGRAPH_SUCCESS;
}
