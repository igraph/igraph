/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
/*
   IGraph library.
   Copyright (C) 2024  The igraph development team <igraph@igraph.org>

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

#include "string.h"

#include "igraph_bitset.h"
#include "igraph_memory.h"

igraph_integer_t igraph_i_ctz32(igraph_uint_t x) {
#ifdef HAVE__BITSCANFORWARD
    unsigned long index;
    return _BitScanForward(&index, x) ? index : 32;
#else
    for (igraph_integer_t i = 0; i < 32; ++i) {
        if (IGRAPH_BIT_MASK(i) & x) {
            return i;
        }
    }
    return 32;
#endif
}

igraph_integer_t igraph_i_clz32(igraph_uint_t x) {
#ifdef HAVE_BITSCANREVERSE
    unsigned long index;
    return _BitScanReverse(&index, x) ? 31 - index : 32;
#else
    for (igraph_integer_t i = 0; i >= 0; --i) {
        if (IGRAPH_BIT_MASK(i) & x) {
            return 31 - i;
        }
    }
    return 32;
#endif
}

igraph_integer_t igraph_i_popcnt(igraph_uint_t x) {
    igraph_integer_t result = 0;
    while (x) {
        result++;
        x = x & (x - 1);
    }
    return result;
}

/* Fallbacks for 64-bit word (and igraph_integer_t/igraph_uint_t) size */
#if IGRAPH_INTEGER_SIZE == 64
igraph_integer_t igraph_i_ctz64(igraph_uint_t x) {
#ifdef HAVE_BITSCANFORWARD64
    unsigned long index;
    return _BitScanForward64(&index, x) ? index : 64;
#else
    for (igraph_integer_t i = 0; i < 64; ++i) {
        if (IGRAPH_BIT_MASK(i) & x) {
            return i;
        }
    }
    return 64;
#endif
}

igraph_integer_t igraph_i_clz64(igraph_uint_t x) {
#ifdef HAVE_BITSCANREVERSE64
    unsigned long index;
    return _BitScanReverse64(&index, x) ? 63 - index : 64;
#else
    for (igraph_integer_t i = 63; i >= 0; --i) {
        if (IGRAPH_BIT_MASK(i) & x) {
            return 63 - i;
        }
    }
    return 64;
#endif
}
#endif /* IGRAPH_INTEGER_SIZE == 64 */

/**
 * \ingroup bitset
 * \section about_igraph_bitset_t_objects About \type igraph_bitset_t objects
 *
 * <para>The \type igraph_bitset_t data type is a simple and efficient
 * interface to arrays containing boolean values. It is similar to
 * the \type bitset template in the C++ standard library, although the main
 * difference being the C++ version's size is initialized at compile time.</para>
 *
 * <para>The \type igraph_bitset_t type and use O(n/w) space
 * to store n elements, where w is the bit width of \type igraph_integer_t,
 * the integer type used throughout the library (either 32 or 64).
 * Sometimes they use more, this is because bitsets can
 * shrink, but even if they shrink, the current implementation does not free a
 * single bit of memory.</para>
 *
 * <para>The elements in an \type igraph_bitset_t object and its variants are
 * indexed from zero, we follow the usual C convention here. Bitsets are indexed
 * from right to left, meaning index 0 is the least significant bit and index n-1
 * is the most significant bit.</para>
 *
 * <para>The elements of a bitset always occupy a single block of
 * memory, the starting address of this memory block can be queried
 * with the \ref VECTOR() macro. This way, bitset objects can be used
 * with standard mathematical libraries, like the GNU Scientific
 * Library.</para>
 */

/**
 * \ingroup bitset
 * \section igraph_bitset_constructors_and_destructors Constructors and
 * destructors
 *
 * <para>\type igraph_bitset_t objects have to be initialized before using
 * them, this is analogous to calling a constructor on them. There are two
 * \type igraph_bitset_t constructors, for your convenience.
 * \ref igraph_bitset_init() is the basic constructor, it
 * creates a bitset of the given length, filled with zeros.
 * \ref igraph_bitset_init_copy() creates a new identical copy
 * of an already existing and initialized bitset.</para>
 *
 * <para>If a \type igraph_bitset_t object is not needed any more, it
 * should be destroyed to free its allocated memory by calling the
 * \type igraph_bitset_t destructor, \ref igraph_bitset_destroy().</para>
 */

/**
 * \ingroup bitset
 * \function igraph_bitset_init
 * \brief Initializes a bitset object (constructor).
 *
 * \experimental
 *
 * </para><para>
 * Every bitset needs to be initialized before it can be used, and
 * there are a number of initialization functions or otherwise called
 * constructors. This function constructs a bitset of the given size and
 * initializes each entry to 0.
 *
 * </para><para>
 * Every bitset object initialized by this function should be
 * destroyed (ie. the memory allocated for it should be freed) when it
 * is not needed anymore, the \ref igraph_bitset_destroy() function is
 * responsible for this.
 *
 * \param bitset Pointer to a not yet initialized bitset object.
 * \param size The size of the bitset.
 * \return error code:
 *       \c IGRAPH_ENOMEM if there is not enough memory.
 *
 * Time complexity: operating system dependent, the amount of
 * \quote time \endquote required to allocate
 * O(n/w) elements,
 * n is the number of elements.
 * w is the word size of the machine (32 or 64).
 */

igraph_error_t igraph_bitset_init(igraph_bitset_t *bitset, igraph_integer_t size) {
    igraph_integer_t alloc_size = IGRAPH_BIT_NSLOTS(size);
    bitset->stor_begin = IGRAPH_CALLOC(alloc_size, igraph_uint_t);
    IGRAPH_CHECK_OOM(bitset->stor_begin, "Cannot initialize bitset");
    bitset->size = size;
    bitset->stor_end = bitset->stor_begin + alloc_size;
    return IGRAPH_SUCCESS;
}

/**
 * \ingroup bitset
 * \function igraph_bitset_destroy
 * \brief Destroys a bitset object.
 *
 * \experimental
 *
 * </para><para>
 * All bitsets initialized by \ref igraph_bitset_init() should be properly
 * destroyed by this function. A destroyed bitset needs to be
 * reinitialized by \ref igraph_bitset_init() or
 * another constructor.
 *
 * \param bitset Pointer to the (previously initialized) bitset object to
 *        destroy.
 *
 * Time complexity: operating system dependent.
 */

void igraph_bitset_destroy(igraph_bitset_t *bitset) {
    IGRAPH_ASSERT(bitset != NULL);
    IGRAPH_FREE(bitset->stor_begin);
    bitset->size = 0;
}

/**
 * \ingroup bitset
 * \function igraph_bitset_init_copy
 * \brief Initializes a bitset from another bitset object (constructor).
 *
 * \experimental
 *
 * </para><para>
 *
 * The contents of the existing bitset object will be copied to
 * the new one.
 * \param dest Pointer to a not yet initialized bitset object.
 * \param src The original bitset object to copy.
 * \return Error code:
 *         \c IGRAPH_ENOMEM if there is not enough memory.
 *
 * Time complexity: operating system dependent, usually
 * O(n/w),
 * n is the size of the bitset,
 * w is the word size of the machine (32 or 64).
 */

igraph_error_t igraph_bitset_init_copy(igraph_bitset_t *dest, const igraph_bitset_t *src) {
    IGRAPH_ASSERT(src != NULL);
    IGRAPH_ASSERT(src->stor_begin != NULL);
    IGRAPH_CHECK(igraph_bitset_init(dest, src->size));
    for (igraph_integer_t i = 0; i < IGRAPH_BIT_NSLOTS(dest->size); ++i) {
        VECTOR(*dest)[i] = VECTOR(*src)[i];
    }
    return IGRAPH_SUCCESS;
}

/**
 * \ingroup bitset
 * \function igraph_bitset_capacity
 * \brief Returns the allocated capacity of the bitset.
 *
 * \experimental
 *
 * Note that this might be different from the size of the bitset (as
 * queried by \ref igraph_bitset_size()), and specifies how many elements
 * the bitset can hold, without reallocation.
 *
 * \param bitset Pointer to the (previously initialized) bitset object
 *          to query.
 * \return The allocated capacity.
 *
 * \sa \ref igraph_bitset_size().
 *
 * Time complexity: O(1).
 */

igraph_integer_t igraph_bitset_capacity(const igraph_bitset_t *bitset) {
    return IGRAPH_INTEGER_SIZE * (bitset->stor_end - bitset->stor_begin);
}

/**
 * \ingroup bitset
 * \function igraph_bitset_size
 * \brief Returns the length of the bitset.
 *
 * \experimental
 *
 * \param bitset The bitset object
 * \return The size of the bitset.
 *
 * Time complexity: O(1).
 */

igraph_integer_t igraph_bitset_size(const igraph_bitset_t *bitset) {
    return bitset->size;
}

/**
 * \ingroup bitset
 * \function igraph_bitset_reserve
 * \brief Reserves memory for a bitset.
 *
 * \experimental
 *
 * </para><para>
 * \a igraph bitsets are flexible, they can grow and
 * shrink. Growing
 * however occasionally needs the data in the bitset to be copied.
 * In order to avoid this, you can call this function to reserve space for
 * future growth of the bitset.
 *
 * </para><para>
 * Note that this function does \em not change the size of the
 * bitset. Let us see a small example to clarify things: if you
 * reserve space for 100 elements and the size of your
 * bitset was (and still is) 60, then you can surely add additional 40
 * elements to your bitset before it will be copied.
 * \param bitset The bitset object.
 * \param capacity The new \em allocated size of the bitset.
 * \return Error code:
 *         \c IGRAPH_ENOMEM if there is not enough memory.
 *
 * Time complexity: operating system dependent, should be around
 * O(n/w),
 * n is the new allocated size of the bitset,
 * w is the word size of the machine (32 or 64).
 */

igraph_error_t igraph_bitset_reserve(igraph_bitset_t *bitset, igraph_integer_t capacity) {
    igraph_integer_t current_capacity;
    igraph_uint_t *tmp;

    IGRAPH_ASSERT(bitset != NULL);
    IGRAPH_ASSERT(bitset->stor_begin != NULL);
    IGRAPH_ASSERT(capacity >= 0);

    current_capacity = igraph_bitset_capacity(bitset);

    if (IGRAPH_BIT_NSLOTS(capacity) <= IGRAPH_BIT_NSLOTS(current_capacity)) {
        return IGRAPH_SUCCESS;
    }

    tmp = IGRAPH_REALLOC(bitset->stor_begin, IGRAPH_BIT_NSLOTS(capacity), igraph_uint_t);
    IGRAPH_CHECK_OOM(tmp, "Cannot reserve space for bitset.");

    bitset->stor_begin = tmp;
    bitset->stor_end = bitset->stor_begin + IGRAPH_BIT_NSLOTS(capacity);

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup bitset
 * \function igraph_bitset_resize
 * \brief Resizes the bitset.
 *
 * \experimental
 *
 * </para><para>
 * Note that this function does not free any memory, just sets the
 * size of the bitset to the given one. It may, on the other hand,
 * allocate more memory if the new size is larger than the previous
 * one. In this case the newly appeared elements in the bitset are
 * set to zero.
 *
 * \param bitset The bitset object
 * \param new_size The new size of the bitset.
 * \return Error code,
 *         \c IGRAPH_ENOMEM if there is not enough
 *         memory. Note that this function \em never returns an error
 *         if the bitset is made smaller.
 * \sa \ref igraph_bitset_reserve() for allocating memory for future
 * extensions of a bitset.
 *
 * Time complexity: O(1) if the new
 * size is smaller, operating system dependent if it is larger. In the
 * latter case it is usually around
 * O(n/w),
 * n is the new size of the bitset,
 * w is the word size of the machine (32 or 64).
 */

igraph_error_t igraph_bitset_resize(igraph_bitset_t *bitset, igraph_integer_t new_size) {
    IGRAPH_ASSERT(bitset != NULL);
    IGRAPH_ASSERT(bitset->stor_begin != NULL);
    IGRAPH_CHECK(igraph_bitset_reserve(bitset, new_size));

    if (new_size > bitset->size) {
        for (igraph_integer_t i = bitset->size; i % IGRAPH_INTEGER_SIZE != 0; ++i) {
            IGRAPH_BIT_CLEAR(*bitset, i);
        }
        memset(bitset->stor_begin + IGRAPH_BIT_NSLOTS(bitset->size), 0,
               sizeof(igraph_integer_t) * (IGRAPH_BIT_NSLOTS(new_size) - IGRAPH_BIT_NSLOTS(bitset->size)));
    }
    bitset->size = new_size;

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup bitset
 * \function igraph_bitset_popcount
 * \brief The population count of the bitset.
 *
 * \experimental
 *
 * Returns the number of set bits, also called the population count,
 * of the bitset.
 *
 * \param bitset The bitset object
 * \return The population count of the bitset.
 *
 * Time complexity: O(n/w).
 */

igraph_integer_t igraph_bitset_popcount(const igraph_bitset_t *bitset) {
    const igraph_integer_t final_block_size = bitset->size % IGRAPH_INTEGER_SIZE ? bitset->size % IGRAPH_INTEGER_SIZE : IGRAPH_INTEGER_SIZE;
    const igraph_integer_t slots = IGRAPH_BIT_NSLOTS(bitset->size);
    const igraph_uint_t one = 1, zero = 0; /* to avoid the need to cast 1 and 0 to igraph_uint_t below */
    const igraph_uint_t mask = final_block_size == IGRAPH_INTEGER_SIZE ? ~zero : ((one << final_block_size) - 1);
    igraph_integer_t count = 0;

    for (igraph_integer_t i = 0; i + 1 < slots; ++i) {
        count += IGRAPH_POPCOUNT(VECTOR(*bitset)[i]);
    }
    if (bitset->size) {
        count += IGRAPH_POPCOUNT(mask & VECTOR(*bitset)[slots - 1]);
    }

    return count;
}

/**
 * \ingroup bitset
 * \function igraph_bitset_countl_zero
 * \brief The number of leading zeros in the bitset.
 *
 * \experimental
 *
 * Returns the number of leading (starting at the most significant bit)
 * zeros in the bitset before the first one is encountered. If the bitset
 * is all zeros, then its size is returned.
 *
 * \param bitset The bitset object
 * \return The number of leading zeros in the bitset.
 *
 * Time complexity: O(n/w).
 */

igraph_integer_t igraph_bitset_countl_zero(const igraph_bitset_t *bitset) {
    const igraph_integer_t final_block_size = bitset->size % IGRAPH_INTEGER_SIZE ? bitset->size % IGRAPH_INTEGER_SIZE : IGRAPH_INTEGER_SIZE;
    const igraph_integer_t padding = IGRAPH_INTEGER_SIZE - final_block_size;
    const igraph_integer_t slots = IGRAPH_BIT_NSLOTS(bitset->size);
    const igraph_uint_t one = 1, zero = 0;
    const igraph_uint_t mask = final_block_size == IGRAPH_INTEGER_SIZE ? ~zero : ((one << final_block_size) - one);

    if (bitset->size && (mask & VECTOR(*bitset)[slots - 1]) != 0) {
        return IGRAPH_CLZ(mask & VECTOR(*bitset)[slots - 1]) - padding;
    }
    for (igraph_integer_t i = 1; i < slots; ++i) {
        if (VECTOR(*bitset)[slots - i - 1] != 0) {
            return IGRAPH_INTEGER_SIZE * i + IGRAPH_CLZ(VECTOR(*bitset)[slots - i - 1]) - padding;
        }
    }

    return bitset->size;
}

/**
 * \ingroup bitset
 * \function igraph_bitset_countl_one
 * \brief The number of leading ones in the bitset.
 *
 * \experimental
 *
 * Returns the number of leading ones (starting at the most significant bit)
 * in the bitset before the first zero is encountered.
 * If the bitset is all ones, then its size is returned.
 *
 * \param bitset The bitset object
 * \return The number of leading ones in the bitset.
 *
 * Time complexity: O(n/w).
 */

igraph_integer_t igraph_bitset_countl_one(const igraph_bitset_t *bitset) {
    const igraph_integer_t final_block_size = bitset->size % IGRAPH_INTEGER_SIZE ? bitset->size % IGRAPH_INTEGER_SIZE : IGRAPH_INTEGER_SIZE;
    const igraph_integer_t padding = IGRAPH_INTEGER_SIZE - final_block_size;
    const igraph_integer_t slots = IGRAPH_BIT_NSLOTS(bitset->size);
    const igraph_uint_t one = 1, zero = 0; /* to avoid the need to cast 1 and 0 to igraph_uint_t below */
    const igraph_uint_t mask = final_block_size == IGRAPH_INTEGER_SIZE ? zero : ~((one << final_block_size) - one);

    if (bitset->size && (mask | VECTOR(*bitset)[slots - 1]) != ~zero) {
        return IGRAPH_CLO(mask | VECTOR(*bitset)[slots - 1]) - padding;
    }
    for (igraph_integer_t i = 1; i < slots; ++i) {
        if (VECTOR(*bitset)[slots - i - 1] != ~zero) {
            return IGRAPH_INTEGER_SIZE * i + IGRAPH_CLO(VECTOR(*bitset)[slots - i - 1]) - padding;
        }
    }

    return bitset->size;
}

/**
 * \ingroup bitset
 * \function igraph_bitset_countr_zero
 * \brief The number of trailing zeros in the bitset.
 *
 * \experimental
 *
 * Returns the number of trailing (starting at the least significant bit)
 * zeros in the bitset before the first one is encountered.
 * If the bitset is all zeros, then its size is returned.
 *
 * \param bitset The bitset object
 * \return The number of trailing zeros in the bitset.
 *
 * Time complexity: O(n/w).
 */

igraph_integer_t igraph_bitset_countr_zero(const igraph_bitset_t *bitset) {
    const igraph_integer_t final_block_size = bitset->size % IGRAPH_INTEGER_SIZE ? bitset->size % IGRAPH_INTEGER_SIZE : IGRAPH_INTEGER_SIZE;
    const igraph_integer_t slots = IGRAPH_BIT_NSLOTS(bitset->size);
    const igraph_uint_t one = 1, zero = 0; /* to avoid the need to cast 1 and 0 to igraph_uint_t below */
    const igraph_uint_t mask = final_block_size == IGRAPH_INTEGER_SIZE ? ~zero : ((one << final_block_size) - one);

    for (igraph_integer_t i = 0; i + 1 < slots; ++i) {
        if (VECTOR(*bitset)[i] != zero) {
            return IGRAPH_INTEGER_SIZE * i + IGRAPH_CTZ(VECTOR(*bitset)[i]);
        }
    }
    if (bitset->size && (mask & VECTOR(*bitset)[slots - 1]) != zero) {
        return IGRAPH_INTEGER_SIZE * (slots - 1) + IGRAPH_CTZ(mask & VECTOR(*bitset)[slots - 1]);
    }

    return bitset->size;
}

/**
 * \ingroup bitset
 * \function igraph_bitset_countr_one
 * \brief The number of trailing ones in the bitset.
 *
 * \experimental
 *
 * Returns the number of trailing ones (starting at the least significant bit)
 * in the bitset before the first zero is encountered.
 * If the bitset is all ones, then its size is returned.
 *
 * \param bitset The bitset object
 * \return The number of trailing ones in the bitset.
 *
 * Time complexity: O(n/w).
 */

igraph_integer_t igraph_bitset_countr_one(const igraph_bitset_t *bitset) {
    const igraph_integer_t final_block_size = bitset->size % IGRAPH_INTEGER_SIZE ? bitset->size % IGRAPH_INTEGER_SIZE : IGRAPH_INTEGER_SIZE;
    const igraph_integer_t slots = IGRAPH_BIT_NSLOTS(bitset->size);
    const igraph_uint_t one = 1, zero = 0; /* to avoid the need to cast 1 and 0 to igraph_uint_t below */
    const igraph_uint_t mask = final_block_size == IGRAPH_INTEGER_SIZE ? zero : ~((one << final_block_size) - one);

    for (igraph_integer_t i = 0; i + 1 < slots; ++i) {
        if (VECTOR(*bitset)[i] != ~zero) {
            return IGRAPH_INTEGER_SIZE * i + IGRAPH_CTO(VECTOR(*bitset)[i]);
        }
    }
    if (bitset->size && (mask | VECTOR(*bitset)[slots - 1]) != ~zero) {
        return IGRAPH_INTEGER_SIZE * (slots - 1) + IGRAPH_CTO(mask | VECTOR(*bitset)[slots - 1]);
    }

    return bitset->size;
}

/**
 * \ingroup bitset
 * \function igraph_bitset_or
 * \brief Bitwise OR of two bitsets.
 *
 * \experimental
 *
 * Applies a bitwise or to the contents of two bitsets and stores it in an
 * already initialized bitset. The destination bitset may be equal to one
 * (or even both) of the sources. When working with bitsets, it is common
 * that those created are of the same size fixed size. Therefore, this
 * function does not check the sizes of the bitsets passed to it, the caller
 * must do so if necessary.
 *
 * \param dest The bitset object where the result is stored
 * \param src1 A bitset
 * \param src2 A bitset
 *
 * Time complexity: O(n/w).
 */

void igraph_bitset_or(igraph_bitset_t *dest,
                      const igraph_bitset_t *src1, const igraph_bitset_t *src2) {
    for (igraph_integer_t i = 0; i < IGRAPH_BIT_NSLOTS(dest->size); ++i) {
        VECTOR(*dest)[i] = VECTOR(*src1)[i] | VECTOR(*src2)[i];
    }
}

/**
 * \ingroup bitset
 * \function igraph_bitset_and
 * \brief Bitwise AND of two bitsets.
 *
 * \experimental
 *
 * Applies a bitwise and to the contents of two bitsets and stores it in an
 * already initialized bitset. The destination bitset may be equal to one
 * (or even both) of the sources. When working with bitsets, it is common
 * that those created are of the same size fixed size. Therefore, this
 * function does not check the sizes of the bitsets passed to it, the caller
 * must do so if necessary.
 *
 * \param dest The bitset object where the result is stored
 * \param src1 A bitset
 * \param src2 A bitset
 *
 * Time complexity: O(n/w).
 */

void igraph_bitset_and(igraph_bitset_t *dest, const igraph_bitset_t *src1, const igraph_bitset_t *src2) {
    for (igraph_integer_t i = 0; i < IGRAPH_BIT_NSLOTS(dest->size); ++i) {
        VECTOR(*dest)[i] = VECTOR(*src1)[i] & VECTOR(*src2)[i];
    }
}

/**
 * \ingroup bitset
 * \function igraph_bitset_xor
 * \brief Bitwise XOR of two bitsets.
 *
 * \experimental
 *
 * Applies a bitwise xor to the contents of two bitsets and stores it in
 * an already initialized bitset. The destination bitset may be equal to
 * one (or even both) of the sources. When working with bitsets, it is common
 * that those created are of the same size fixed size. Therefore, this
 * function does not check the sizes of the bitsets passed to it, the caller
 * must do so if necessary.
 *
 * \param dest The bitset object where the result is stored
 * \param src1 A bitset
 * \param src2 A bitset
 *
 * Time complexity: O(n/w).
 */

void igraph_bitset_xor(igraph_bitset_t *dest,
                       const igraph_bitset_t *src1, const igraph_bitset_t *src2) {
    for (igraph_integer_t i = 0; i < IGRAPH_BIT_NSLOTS(dest->size); ++i) {
        VECTOR(*dest)[i] = VECTOR(*src1)[i] ^ VECTOR(*src2)[i];
    }
}

/**
 * \ingroup bitset
 * \function igraph_bitset_not
 * \brief Bitwise negation of a bitset.
 *
 * \experimental
 *
 * Applies a bitwise not to the contents of a bitset and stores it in an
 * already initialized bitset. The destination bitset may be equal to the
 * source. When working with bitsets, it is common that those created are
 * of the same size fixed size. Therefore, this function does not check the
 * sizes of the bitsets passed to it, the caller must do so if necessary.
 *
 * \param dest The bitset object where the result is stored
 * \param src A bitset
 *
 * Time complexity: O(n/w).
 */

void igraph_bitset_not(igraph_bitset_t *dest, const igraph_bitset_t *src) {
    for (igraph_integer_t i = 0; i < IGRAPH_BIT_NSLOTS(dest->size); ++i) {
        VECTOR(*dest)[i] = ~VECTOR(*src)[i];
    }
}

/**
 * \ingroup bitset
 * \function igraph_bitset_fprint
 * \brief Prints the bits of a bitset.
 *
 * \experimental
 *
 * Outputs the contents of a bitset to a file.
 * The bitset is written from index n-1 to index 0, left to right,
 * such that index 0 is the least significant bit and index n-1 is
 * the most significant bit, where n is the size of the bitset.
 * This is the reverse of how sequential structures are usually written,
 * such as vectors, but consistent with the bitset being a binary
 * representation of an integer and how they are usually written.
 *
 * </para><para>
 * No newline is printed at the end.
 *
 * \param bitset The bitset to be printed.
 * \param file The file to be written.
 * \return Error code.
 *
 * Time complexity: O(n).
 */
igraph_error_t igraph_bitset_fprint(const igraph_bitset_t *bitset, FILE *file) {
    for (igraph_integer_t i = bitset->size - 1; i >= 0; i--) {
        fputc(IGRAPH_BIT_TEST(*bitset, i) ? '1' : '0', file);
    }
    return IGRAPH_SUCCESS;
}

#ifndef USING_R
igraph_error_t igraph_bitset_print(const igraph_bitset_t *bitset) {
    return igraph_bitset_fprint(bitset, stdout);
}
#endif
