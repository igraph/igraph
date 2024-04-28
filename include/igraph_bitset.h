/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2024  The igraph development team

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

#ifndef IGRAPH_BITSET_H
#define IGRAPH_BITSET_H

#include "igraph_decls.h"
#include "igraph_error.h"
#include "igraph_vector.h"

/* Required for MSVC intrinsics such as __popcnt and __popcnt64 */
#ifdef _MSC_VER
    #include "intrin.h"
#endif

__BEGIN_DECLS

/**
 * \ingroup bitset
 * \section igraph_bitset_accessing_elements Accessing elements
 *
 * <para>The simplest way to access an element of a bitset is
 * to use the \ref IGRAPH_BIT_TEST(), \ref IGRAPH_BIT_SET() and \ref IGRAPH_BIT_CLEAR() macros.
 * </para>
 *
 * <para>There are a few other macros which allow manual manipulation of bitsets.
 * Those are \ref VECTOR(), \ref IGRAPH_BIT_SLOT(), \ref IGRAPH_BIT_MASK() and
 * \ref IGRAPH_BIT_NSLOTS().</para>
 */

/**
 * \ingroup bitset
 * \define IGRAPH_BIT_MASK
 * \brief Computes mask used to access a specific bit of an integer.
 *
 * \experimental
 *
 * Used in combination with \ref IGRAPH_BIT_SLOT() to access an element of a bitset.
 *
 * Usage:
 * \verbatim IGRAPH_BIT_MASK(10) \endverbatim
 * to obtain an integer where only the 11th least significant bit is set.
 *
 * Note that passing negative values here results in undefined behaviour.
 *
 * \param b The only bit index that should have its bit set.
 *
 * Time complexity: O(1).
 */
#define IGRAPH_BIT_MASK(b) ((igraph_uint_t)(1) << ((b) % IGRAPH_INTEGER_SIZE))

/**
 * \ingroup bitset
 * \define IGRAPH_BIT_SLOT
 * \brief Computes index used to access a specific slot of a bitset.
 *
 * \experimental
 *
 * Used in combination with \ref IGRAPH_BIT_MASK to access an element of a bitset.
 *
 * Usage:
 * \verbatim IGRAPH_BIT_SLOT(70) \endverbatim
 * will return 1 if using 64-bit words or 2 if using 32-bit words.
 *
 * \param b The bit index whose slot should be determined.
 *
 * Time complexity: O(1).
 */
#define IGRAPH_BIT_SLOT(b) ((b) / IGRAPH_INTEGER_SIZE)

/**
 * \ingroup bitset
 * \define IGRAPH_BIT_SET
 * \brief Sets a specific bit in a bitset to 1 without altering other bits.
 *
 * \experimental
 *
 * Usage:
 * \verbatim IGRAPH_BIT_SET(bitset, 3) \endverbatim
 * will set the fourth least significant bit in the bitset to 1.
 *
 * \param a The bitset
 * \param b The bit index that should have its bit set to 1 after the operation.
 *
 * Time complexity: O(1).
 */
#define IGRAPH_BIT_SET(a, b) (VECTOR((a))[IGRAPH_BIT_SLOT(b)] |= IGRAPH_BIT_MASK(b))

/**
 * \ingroup bitset
 * \define IGRAPH_BIT_CLEAR
 * \brief Sets a specific bit in a bitset to 0 without altering other bits.
 *
 * \experimental
 *
 * Usage:
 * \verbatim IGRAPH_BIT_CLEAR(bitset, 4) \endverbatim
 * will set the fifth least significant bit in the bitset to 0.
 *
 * \param a The bitset
 * \param b The bit index that should have its bit set to 0 after the operation.
 *
 * Time complexity: O(1).
 */
#define IGRAPH_BIT_CLEAR(a, b) (VECTOR((a))[IGRAPH_BIT_SLOT(b)] &= ~IGRAPH_BIT_MASK(b))

/**
 * \ingroup bitset
 * \define IGRAPH_BIT_TEST
 * \brief Tests whether a bit is set in a bitset.
 *
 * \experimental
 *
 * Returns 0 if the bit at the specified bit index is not set,
 * otherwise returns a non-zero value.
 *
 * Usage:
 * \verbatim IGRAPH_BIT_TEST(bitset, 7) \endverbatim
 * will test the eighth least significant bit in the bitset.
 *
 * \param a The bitset
 * \param b The bit index that should have its bit tested.
 *
 * Time complexity: O(1).
 */
#define IGRAPH_BIT_TEST(a, b) (VECTOR((a))[IGRAPH_BIT_SLOT(b)] & IGRAPH_BIT_MASK(b))

/**
 * \ingroup bitset
 * \define IGRAPH_BIT_NSLOTS
 * \brief Computes the number of slots required to store a specified number of bits.
 *
 * \experimental
 *
 * Usage:
 * \verbatim IGRAPH_BIT_NSLOTS(70) \endverbatim
 * will return 2 if using 64-bit words and 3 if using 32-bit words.
 * \verbatim IGRAPH_BIT_NSLOTS(128) \endverbatim
 * will return 2 if using 64-bit words and 4 if using 32-bit words.
 *
 * \param nb The specified number of bits.
 *
 * Time complexity: O(1).
 */
#define IGRAPH_BIT_NSLOTS(nb) ((nb + IGRAPH_INTEGER_SIZE - (igraph_integer_t)(1)) / IGRAPH_INTEGER_SIZE)


#if defined(__GNUC__)
    /* GCC and Clang support these six builtins from very early versions. */

    #define IGRAPH_I_POPCOUNT32(x) __builtin_popcount(x)
    #define IGRAPH_I_POPCOUNT64(x) __builtin_popcountll(x)

    /* The result of the following four builtins is undefined for zero input,
     * therefore we handle this specially.
     * See https://gcc.gnu.org/onlinedocs/gcc/Other-Builtins.html */
    #define IGRAPH_I_CTZ32(x) (x ? __builtin_ctz(x) : 32)
    #define IGRAPH_I_CTZ64(x) (x ? __builtin_ctzll(x) : 64)
    #define IGRAPH_I_CLZ32(x) (x ? __builtin_clz(x) : 32)
    #define IGRAPH_I_CLZ64(x) (x ? __builtin_clzll(x) : 64)

#else
    /* Non-GNU compilers, i.e. MSVC and others. Attempt to use MSVC intrinsics
     * when available, otherwise use fallback implementation. */

    #if IGRAPH_INTEGER_SIZE == 32
        #ifdef HAVE__POPCNT
            #define IGRAPH_I_POPCOUNT32(x) __popcnt(x)
        #else
            igraph_integer_t igraph_i_popcnt(igraph_uint_t x);
            #define IGRAPH_I_POPCOUNT32(x) igraph_i_popcnt(x)
        #endif
    #elif IGRAPH_INTEGER_SIZE == 64
        #ifdef HAVE__POPCNT64
            #define IGRAPH_I_POPCOUNT64(x) __popcnt64(x)
        #else
            igraph_integer_t igraph_i_popcnt(igraph_uint_t x);
            #define IGRAPH_I_POPCOUNT64(x) igraph_i_popcnt(x)
        #endif
    #else
        #error "Unexpected IGRAPH_INTEGER_SIZE value."
    #endif

    igraph_integer_t igraph_i_ctz32(igraph_uint_t x);
    igraph_integer_t igraph_i_ctz64(igraph_uint_t x);
    igraph_integer_t igraph_i_clz32(igraph_uint_t x);
    igraph_integer_t igraph_i_clz64(igraph_uint_t x);
    #define IGRAPH_I_CTZ32(x) igraph_i_ctz32(x)
    #define IGRAPH_I_CTZ64(x) igraph_i_ctz64(x)
    #define IGRAPH_I_CLZ32(x) igraph_i_clz32(x)
    #define IGRAPH_I_CLZ64(x) igraph_i_clz64(x)

#endif


#if IGRAPH_INTEGER_SIZE == 32
    #define IGRAPH_POPCOUNT IGRAPH_I_POPCOUNT32
    #define IGRAPH_CLZ IGRAPH_I_CLZ32
    #define IGRAPH_CTZ IGRAPH_I_CTZ32
#elif IGRAPH_INTEGER_SIZE == 64
    #define IGRAPH_POPCOUNT IGRAPH_I_POPCOUNT64
    #define IGRAPH_CLZ IGRAPH_I_CLZ64
    #define IGRAPH_CTZ IGRAPH_I_CTZ64
#else
    #error "Unexpected IGRAPH_INTEGER_SIZE value."
#endif

#define IGRAPH_CLO(x) IGRAPH_CLZ(~(x))
#define IGRAPH_CTO(x) IGRAPH_CTZ(~(x))

typedef struct {
    igraph_integer_t size;
    igraph_uint_t *stor_begin;
    igraph_uint_t *stor_end;
} igraph_bitset_t;

IGRAPH_EXPORT igraph_error_t igraph_bitset_init(igraph_bitset_t *bitset, igraph_integer_t size);
IGRAPH_EXPORT void igraph_bitset_destroy(igraph_bitset_t *bitset);
IGRAPH_EXPORT igraph_error_t igraph_bitset_init_copy(igraph_bitset_t *dest, const igraph_bitset_t *src);
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE igraph_integer_t igraph_bitset_capacity(const igraph_bitset_t *bitset);
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE igraph_integer_t igraph_bitset_size(const igraph_bitset_t *bitset);
IGRAPH_EXPORT igraph_error_t igraph_bitset_reserve(igraph_bitset_t *bitset, igraph_integer_t capacity);
IGRAPH_EXPORT igraph_error_t igraph_bitset_resize(igraph_bitset_t *bitset, igraph_integer_t new_size);
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE igraph_integer_t igraph_bitset_popcount(const igraph_bitset_t *bitset);
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE igraph_integer_t igraph_bitset_countl_zero(const igraph_bitset_t *bitset);
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE igraph_integer_t igraph_bitset_countl_one(const igraph_bitset_t *bitset);
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE igraph_integer_t igraph_bitset_countr_zero(const igraph_bitset_t *bitset);
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE igraph_integer_t igraph_bitset_countr_one(const igraph_bitset_t *bitset);
IGRAPH_EXPORT void igraph_bitset_or(igraph_bitset_t *dest, const igraph_bitset_t *src1, const igraph_bitset_t *src2);
IGRAPH_EXPORT void igraph_bitset_and(igraph_bitset_t *dest, const igraph_bitset_t *src1, const igraph_bitset_t *src2);
IGRAPH_EXPORT void igraph_bitset_xor(igraph_bitset_t *dest, const igraph_bitset_t *src1, const igraph_bitset_t *src2);
IGRAPH_EXPORT void igraph_bitset_not(igraph_bitset_t *dest, const igraph_bitset_t *src);
IGRAPH_EXPORT igraph_error_t igraph_bitset_fprint(const igraph_bitset_t *bitset, FILE *file);
IGRAPH_EXPORT igraph_error_t igraph_bitset_print(const igraph_bitset_t *bitset);
__END_DECLS

#endif /* IGRAPH_BITSET_H */
