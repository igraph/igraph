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

igraph_error_t igraph_bitset_init(igraph_vector_int_t *bitset, igraph_integer_t n) {
    IGRAPH_CHECK(igraph_vector_int_init(bitset, IGRAPH_BITNSLOTS(n)));
}

void igraph_bitset_destroy(igraph_vector_int_t *bitset) {
    igraph_vector_int_destroy(bitset);
}

igraph_integer_t igraph_bitset_popcount(igraph_vector_int_t *bitset, igraph_integer_t n)
{
    const igraph_integer_t final_block_size = n % IGRAPH_INTEGER_SIZE ? n % IGRAPH_INTEGER_SIZE : IGRAPH_INTEGER_SIZE;
    const igraph_integer_t padding = IGRAPH_INTEGER_SIZE - final_block_size;
    const igraph_integer_t slots = IGRAPH_BITNSLOTS(n);
    const igraph_integer_t mask = final_block_size == IGRAPH_INTEGER_SIZE ? ~0 : ((1 << final_block_size) - 1);
    igraph_integer_t count = 0;
    for (igraph_integer_t i = 0; i + 1 < slots; ++i)
    {
        count += IGRAPH_POPCOUNT(VECTOR(*bitset)[i]);
    }
    if (n) {
        count += IGRAPH_POPCOUNT(mask & VECTOR(*bitset)[slots - 1]);
    }
    return count;
}

igraph_integer_t igraph_bitset_countl_zero(igraph_vector_int_t *bitset, igraph_integer_t n)
{
    const igraph_integer_t final_block_size = n % IGRAPH_INTEGER_SIZE ? n % IGRAPH_INTEGER_SIZE : IGRAPH_INTEGER_SIZE;
    const igraph_integer_t padding = IGRAPH_INTEGER_SIZE - final_block_size;
    const igraph_integer_t slots = IGRAPH_BITNSLOTS(n);
    const igraph_integer_t one = 1, zero = 0;
    const igraph_integer_t mask = final_block_size == IGRAPH_INTEGER_SIZE ? ~zero : ((one << final_block_size) - one);
    if (n && (mask & VECTOR(*bitset)[slots - 1]) != 0) {
        return IGRAPH_CLZ(mask & VECTOR(*bitset)[slots - 1]) - padding;
    }
    for (igraph_integer_t i = 1; i < slots; ++i) {
        if (VECTOR(*bitset)[slots - i - 1] != 0) {
            const igraph_integer_t result = IGRAPH_INTEGER_SIZE * i + IGRAPH_CLZ(VECTOR(*bitset)[slots - i - 1]);
            return result - padding;
        }
    }
    return n;
}

igraph_integer_t igraph_bitset_countl_one(igraph_vector_int_t *bitset, igraph_integer_t n)
{
    const igraph_integer_t final_block_size = n % IGRAPH_INTEGER_SIZE ? n % IGRAPH_INTEGER_SIZE : IGRAPH_INTEGER_SIZE;
    const igraph_integer_t padding = IGRAPH_INTEGER_SIZE - final_block_size;
    const igraph_integer_t slots = IGRAPH_BITNSLOTS(n);
    const igraph_integer_t one = 1, zero = 0;
    const igraph_integer_t mask = final_block_size == IGRAPH_INTEGER_SIZE ? zero : ~((one << final_block_size) - one);
    if (n && (mask | VECTOR(*bitset)[slots - 1]) != ~zero) {
        const igraph_integer_t result = IGRAPH_CLO(mask | VECTOR(*bitset)[slots - 1]);
        return IGRAPH_CLO(mask | VECTOR(*bitset)[slots - 1]) - padding;
    }
    for (igraph_integer_t i = 1; i < slots; ++i) {
        if (VECTOR(*bitset)[slots - i - 1] != ~zero) {
            const igraph_integer_t result = IGRAPH_INTEGER_SIZE * i + IGRAPH_CLO(VECTOR(*bitset)[slots - i - 1]);
            return result - padding;
        }
    }
    return n;
}

igraph_integer_t igraph_bitset_countr_zero(igraph_vector_int_t *bitset, igraph_integer_t n)
{
    const igraph_integer_t final_block_size = n % IGRAPH_INTEGER_SIZE ? n % IGRAPH_INTEGER_SIZE : IGRAPH_INTEGER_SIZE;
    const igraph_integer_t slots = IGRAPH_BITNSLOTS(n);
    const igraph_integer_t one = 1, zero = 0;
    const igraph_integer_t mask = final_block_size == IGRAPH_INTEGER_SIZE ? ~zero : ((one << final_block_size) - one);
    for (igraph_integer_t i = 0; i + 1 < IGRAPH_BITNSLOTS(n); ++i) {
        if (VECTOR(*bitset)[i] != zero) {
            const igraph_integer_t result = IGRAPH_INTEGER_SIZE * i + IGRAPH_CTZ(VECTOR(*bitset)[i]);
            return result;
        }
    }
    if (n && (mask & VECTOR(*bitset)[slots - 1]) != zero) {
        return IGRAPH_INTEGER_SIZE * (slots - 1) + IGRAPH_CTZ(mask & VECTOR(*bitset)[slots - 1]);
    }
    return n;
}

igraph_integer_t igraph_bitset_countr_one(igraph_vector_int_t *bitset, igraph_integer_t n)
{
    const igraph_integer_t final_block_size = n % IGRAPH_INTEGER_SIZE ? n % IGRAPH_INTEGER_SIZE : IGRAPH_INTEGER_SIZE;
    const igraph_integer_t slots = IGRAPH_BITNSLOTS(n);
    const igraph_integer_t one = 1, zero = 0;
    const igraph_integer_t mask = final_block_size == IGRAPH_INTEGER_SIZE ? zero : ~((one << final_block_size) - one);
    for (igraph_integer_t i = 0; i + 1 < slots; ++i) {
        if (VECTOR(*bitset)[i] != ~zero) {
            const igraph_integer_t result = IGRAPH_INTEGER_SIZE * i + IGRAPH_CTO(VECTOR(*bitset)[i]);
            return result;
        }
    }
    if (n && (mask | VECTOR(*bitset)[slots - 1]) != ~zero) {
        return IGRAPH_INTEGER_SIZE * (slots - 1) + IGRAPH_CTO(mask | VECTOR(*bitset)[slots - 1]);
    }
    return n;
}

void igraph_bitset_or(igraph_vector_int_t *dest, igraph_vector_int_t *src1, igraph_vector_int_t *src2, igraph_integer_t n) {
    for (igraph_integer_t i = 0; i < IGRAPH_BITNSLOTS(n); ++i)
    {
        VECTOR(*dest)[i] = VECTOR(*src1)[i] | VECTOR(*src2)[i];
    }
}

void igraph_bitset_and(igraph_vector_int_t *dest, igraph_vector_int_t *src1, igraph_vector_int_t *src2, igraph_integer_t n) {
    for (igraph_integer_t i = 0; i < IGRAPH_BITNSLOTS(n); ++i)
    {
        VECTOR(*dest)[i] = VECTOR(*src1)[i] & VECTOR(*src2)[i];
    }
}

void igraph_bitset_xor(igraph_vector_int_t *dest, igraph_vector_int_t *src1, igraph_vector_int_t *src2, igraph_integer_t n) {
    for (igraph_integer_t i = 0; i < IGRAPH_BITNSLOTS(n); ++i)
    {
        VECTOR(*dest)[i] = VECTOR(*src1)[i] ^ VECTOR(*src2)[i];
    }
}

void igraph_bitset_not(igraph_vector_int_t *dest, igraph_vector_int_t *src, igraph_integer_t n) {
    for (igraph_integer_t i = 0; i < IGRAPH_BITNSLOTS(n); ++i)
    {
        VECTOR(*dest)[i] = ~VECTOR(*src)[i];
    }
}
