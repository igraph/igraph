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

#include "igraph_bitset.h"
#include "igraph_memory.h"

#ifdef _MSC_VER
igraph_integer_t igraph_i_ctz32(igraph_integer_t x) {
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

igraph_integer_t igraph_i_ctz64(igraph_integer_t x) {
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

igraph_integer_t igraph_i_clz32(igraph_integer_t x) {
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

igraph_integer_t igraph_i_clz64(igraph_integer_t x) {
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

#if (!defined(HAVE__POPCNT64) && IGRAPH_INTEGER_SIZE==64) || (!defined(HAVE_POPCNT) && IGRAPH_INTEGER_SIZE==32)
igraph_integer_t igraph_i_popcnt(igraph_integer_t x) {
    igraph_integer_t result = 0;
    while (x) {
        result++;
        x = x & (x-1);
    }
    return result;
}
#endif
#endif

igraph_error_t igraph_bitset_init(igraph_bitset_t *bitset, igraph_integer_t size) {
    igraph_integer_t alloc_size = IGRAPH_BIT_NSLOTS(size);
    bitset->stor_begin = IGRAPH_CALLOC(alloc_size, igraph_integer_t);
    IGRAPH_CHECK_OOM(bitset->stor_begin, "Cannot initialize bitset");
    bitset->size = size;
    bitset->stor_end = bitset->stor_begin + alloc_size;
    return IGRAPH_SUCCESS;
}

void igraph_bitset_destroy(igraph_bitset_t *bitset) {
    IGRAPH_ASSERT(bitset != NULL);
    if (bitset->stor_begin != NULL) {
        IGRAPH_FREE(bitset->stor_begin);
        bitset->size = 0;
        bitset->stor_begin = NULL;
    }
}

igraph_error_t igraph_bitset_init_copy(igraph_bitset_t *dest, const igraph_bitset_t* src) {
    IGRAPH_ASSERT(src != NULL);
    IGRAPH_ASSERT(src->stor_begin != NULL);
    IGRAPH_CHECK(igraph_bitset_init(dest, src->size));
    for (igraph_integer_t i = 0; i < IGRAPH_BIT_NSLOTS(dest->size); ++i)
    {
        VECTOR(*dest)[i] = VECTOR(*src)[i];
    }
    return IGRAPH_SUCCESS;
}

igraph_integer_t igraph_bitset_capacity(igraph_bitset_t *bitset) {
    return IGRAPH_INTEGER_SIZE * (bitset->stor_end - bitset->stor_begin);
}

igraph_integer_t igraph_bitset_size(igraph_bitset_t *bitset) {
    return bitset->size;
}

igraph_error_t igraph_bitset_reserve(igraph_bitset_t *bitset, igraph_integer_t capacity) {
    igraph_integer_t current_capacity;
    igraph_integer_t *tmp;

    IGRAPH_ASSERT(bitset != NULL);
    IGRAPH_ASSERT(bitset->stor_begin != NULL);
    IGRAPH_ASSERT(capacity >= 0);

    current_capacity = igraph_bitset_capacity(bitset);

    if (IGRAPH_BIT_NSLOTS(capacity) <= IGRAPH_BIT_NSLOTS(current_capacity)) {
        return IGRAPH_SUCCESS;
    }

    tmp = IGRAPH_REALLOC(bitset->stor_begin, IGRAPH_BIT_NSLOTS(capacity), igraph_integer_t);
    IGRAPH_CHECK_OOM(tmp, "Cannot reserve space for bitset.");

    bitset->stor_begin = tmp;
    bitset->stor_end = bitset->stor_begin + IGRAPH_BIT_NSLOTS(capacity);

    return IGRAPH_SUCCESS;
}

igraph_error_t igraph_bitset_resize(igraph_bitset_t *bitset, igraph_integer_t new_size) {
    IGRAPH_ASSERT(bitset != NULL);
    IGRAPH_ASSERT(bitset->stor_begin != NULL);
    IGRAPH_CHECK(igraph_bitset_reserve(bitset, new_size));
    bitset->size = new_size;
    return IGRAPH_SUCCESS;
}

igraph_integer_t igraph_bitset_popcount(igraph_bitset_t *bitset)
{
    const igraph_integer_t final_block_size = bitset->size % IGRAPH_INTEGER_SIZE ? bitset->size % IGRAPH_INTEGER_SIZE : IGRAPH_INTEGER_SIZE;
    const igraph_integer_t slots = IGRAPH_BIT_NSLOTS(bitset->size);
    const igraph_integer_t mask = final_block_size == IGRAPH_INTEGER_SIZE ? ~0 : ((1 << final_block_size) - 1);
    igraph_integer_t count = 0;
    for (igraph_integer_t i = 0; i + 1 < slots; ++i)
    {
        count += IGRAPH_POPCOUNT(VECTOR(*bitset)[i]);
    }
    if (bitset->size) {
        count += IGRAPH_POPCOUNT(mask & VECTOR(*bitset)[slots - 1]);
    }
    return count;
}

igraph_integer_t igraph_bitset_countl_zero(igraph_bitset_t *bitset)
{
    const igraph_integer_t final_block_size = bitset->size % IGRAPH_INTEGER_SIZE ? bitset->size % IGRAPH_INTEGER_SIZE : IGRAPH_INTEGER_SIZE;
    const igraph_integer_t padding = IGRAPH_INTEGER_SIZE - final_block_size;
    const igraph_integer_t slots = IGRAPH_BIT_NSLOTS(bitset->size);
    const igraph_integer_t one = 1, zero = 0;
    const igraph_integer_t mask = final_block_size == IGRAPH_INTEGER_SIZE ? ~zero : ((one << final_block_size) - one);
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

igraph_integer_t igraph_bitset_countl_one(igraph_bitset_t *bitset)
{
    const igraph_integer_t final_block_size = bitset->size % IGRAPH_INTEGER_SIZE ? bitset->size % IGRAPH_INTEGER_SIZE : IGRAPH_INTEGER_SIZE;
    const igraph_integer_t padding = IGRAPH_INTEGER_SIZE - final_block_size;
    const igraph_integer_t slots = IGRAPH_BIT_NSLOTS(bitset->size);
    const igraph_integer_t one = 1, zero = 0;
    const igraph_integer_t mask = final_block_size == IGRAPH_INTEGER_SIZE ? zero : ~((one << final_block_size) - one);
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

igraph_integer_t igraph_bitset_countr_zero(igraph_bitset_t *bitset)
{
    const igraph_integer_t final_block_size = bitset->size % IGRAPH_INTEGER_SIZE ? bitset->size % IGRAPH_INTEGER_SIZE : IGRAPH_INTEGER_SIZE;
    const igraph_integer_t slots = IGRAPH_BIT_NSLOTS(bitset->size);
    const igraph_integer_t one = 1, zero = 0;
    const igraph_integer_t mask = final_block_size == IGRAPH_INTEGER_SIZE ? ~zero : ((one << final_block_size) - one);
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

igraph_integer_t igraph_bitset_countr_one(igraph_bitset_t *bitset)
{
    const igraph_integer_t final_block_size = bitset->size % IGRAPH_INTEGER_SIZE ? bitset->size % IGRAPH_INTEGER_SIZE : IGRAPH_INTEGER_SIZE;
    const igraph_integer_t slots = IGRAPH_BIT_NSLOTS(bitset->size);
    const igraph_integer_t one = 1, zero = 0;
    const igraph_integer_t mask = final_block_size == IGRAPH_INTEGER_SIZE ? zero : ~((one << final_block_size) - one);
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

void igraph_bitset_or(igraph_bitset_t *dest, igraph_bitset_t *src1, igraph_bitset_t *src2) {
    for (igraph_integer_t i = 0; i < IGRAPH_BIT_NSLOTS(dest->size); ++i)
    {
        VECTOR(*dest)[i] = VECTOR(*src1)[i] | VECTOR(*src2)[i];
    }
}

void igraph_bitset_and(igraph_bitset_t *dest, igraph_bitset_t *src1, igraph_bitset_t *src2) {
    for (igraph_integer_t i = 0; i < IGRAPH_BIT_NSLOTS(dest->size); ++i)
    {
        VECTOR(*dest)[i] = VECTOR(*src1)[i] & VECTOR(*src2)[i];
    }
}

void igraph_bitset_xor(igraph_bitset_t *dest, igraph_bitset_t *src1, igraph_bitset_t *src2) {
    for (igraph_integer_t i = 0; i < IGRAPH_BIT_NSLOTS(dest->size); ++i)
    {
        VECTOR(*dest)[i] = VECTOR(*src1)[i] ^ VECTOR(*src2)[i];
    }
}

void igraph_bitset_not(igraph_bitset_t *dest, igraph_bitset_t *src) {
    for (igraph_integer_t i = 0; i < IGRAPH_BIT_NSLOTS(dest->size); ++i)
    {
        VECTOR(*dest)[i] = ~VECTOR(*src)[i];
    }
}
