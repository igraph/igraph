/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
/*
   IGraph library.
   Copyright (C) 2009-2020  Gabor Csardi <csardi.gabor@gmail.com>

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

#ifndef IGRAPH_BITSET_H
#define IGRAPH_BITSET_H

#include "igraph_decls.h"
#include "igraph_error.h"
#include "igraph_vector.h"

#ifdef _MSC_VER
#include "intrin.h"
#endif

#include "limits.h"

__BEGIN_DECLS

#define IGRAPH_BITMASK(b) ((igraph_integer_t)(1) << ((b) % IGRAPH_INTEGER_SIZE))
#define IGRAPH_BITSLOT(b) ((b) / IGRAPH_INTEGER_SIZE)
#define IGRAPH_BITSET(a, b) (VECTOR((a))[IGRAPH_BITSLOT(b)] |= IGRAPH_BITMASK(b))
#define IGRAPH_BITCLEAR(a, b) (VECTOR((a))[IGRAPH_BITSLOT(b)] &= ~IGRAPH_BITMASK(b))
#define IGRAPH_BITTEST(a, b) (VECTOR((a))[IGRAPH_BITSLOT(b)] & IGRAPH_BITMASK(b))
#define IGRAPH_BITNSLOTS(nb) ((nb + IGRAPH_INTEGER_SIZE - (igraph_integer_t)(1)) / IGRAPH_INTEGER_SIZE)

#ifdef _MSC_VER
igraph_integer_t igraph_msvc_ctz32(igraph_integer_t x) {
    unsigned long index;
    return _BitScanForward(&index, x) ? index : 32;
}
igraph_integer_t igraph_msvc_ctz64(igraph_integer_t x) {
    unsigned long index;
    return _BitScanForward(&index, x) ? index : 64;
}
// TODO: Check if __cpuid claims the operation is supported by CPU
// https://learn.microsoft.com/en-us/cpp/intrinsics/popcnt16-popcnt-popcnt64?view=msvc-170
#define IGRAPH_POPCOUNT32(x) __popcnt(x)
#define IGRAPH_POPCOUNT64(x) __popcnt64(x)
#define IGRAPH_CTZ32(x) igraph_msvc_ctz32(x)
#define IGRAPH_CTZ64(x) igraph_msvc_ctz64(x)
#else
#define IGRAPH_POPCOUNT32(x) __builtin_popcount(x)
#define IGRAPH_POPCOUNT64(x) __builtin_popcountll(x)
#define IGRAPH_CTZ32(x) (x ? __builtin_ctz(x) : 32)
#define IGRAPH_CTZ64(x) (x ? __builtin_ctzll(x) : 64)
#define IGRAPH_CLZ32(x) (x ? __builtin_clz(x) : 32)
#define IGRAPH_CLZ64(x) (x ? __builtin_clzll(x) : 64)
#endif

#if IGRAPH_INTEGER_SIZE==32
#define IGRAPH_POPCOUNT IGRAPH_POPCOUNT32
#define IGRAPH_CLZ IGRAPH_CLZ32
#define IGRAPH_CTZ IGRAPH_CTZ32
#else
#define IGRAPH_POPCOUNT IGRAPH_POPCOUNT64
#define IGRAPH_CLZ IGRAPH_CLZ64
#define IGRAPH_CTZ IGRAPH_CTZ64
#endif

#define IGRAPH_CLO(x) IGRAPH_CLZ(~(x))
#define IGRAPH_CTO(x) IGRAPH_CTZ(~(x))

typedef struct s_bitset {
    igraph_integer_t size;
    igraph_integer_t* stor_begin;
    igraph_integer_t* stor_end;
} igraph_bitset_t;

IGRAPH_EXPORT igraph_error_t igraph_bitset_init(igraph_bitset_t *bitset, igraph_integer_t size);
IGRAPH_EXPORT void igraph_bitset_destroy(igraph_bitset_t *bitset);
IGRAPH_EXPORT igraph_integer_t igraph_bitset_popcount(igraph_bitset_t* bitset);
IGRAPH_EXPORT igraph_integer_t igraph_bitset_countl_zero(igraph_bitset_t* bitset);
IGRAPH_EXPORT igraph_integer_t igraph_bitset_countl_one(igraph_bitset_t* bitset);
IGRAPH_EXPORT igraph_integer_t igraph_bitset_countr_zero(igraph_bitset_t* bitset);
IGRAPH_EXPORT igraph_integer_t igraph_bitset_countr_one(igraph_bitset_t* bitset);
IGRAPH_EXPORT void igraph_bitset_or(igraph_bitset_t* dest, igraph_bitset_t* src1, igraph_bitset_t* src2);
IGRAPH_EXPORT void igraph_bitset_and(igraph_bitset_t* dest, igraph_bitset_t* src1, igraph_bitset_t* src2);
IGRAPH_EXPORT void igraph_bitset_xor(igraph_bitset_t* dest, igraph_bitset_t* src1, igraph_bitset_t* src2);
IGRAPH_EXPORT void igraph_bitset_not(igraph_bitset_t* dest, igraph_bitset_t* src);

__END_DECLS

#endif // IGRAPH_BITSET_H
