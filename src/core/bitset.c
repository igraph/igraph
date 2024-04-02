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

igraph_integer_t igraph_bitset_popcount(igraph_vector_int_t *bitset, igraph_integer_t n)
{
    igraph_integer_t count = 0;
    for (igraph_integer_t i = 0; i < IGRAPH_BITNSLOTS(n); ++i)
    {
        count += IGRAPH_POPCOUNT(VECTOR(*bitset)[i]);
    }
    return count;
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
