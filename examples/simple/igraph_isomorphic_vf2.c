/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2009-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard st, Cambridge MA, 02139 USA

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
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include <igraph.h>
#include <stdio.h>
#include <stdlib.h>

int main(void) {

    igraph_t ring1, ring2;
    igraph_vector_int_t color1, color2;
    igraph_vector_int_t perm;
    igraph_bool_t iso;
    igraph_integer_t count;
    igraph_integer_t i;

    igraph_rng_seed(igraph_rng_default(), 12345);

    igraph_ring(&ring1, 100, /*directed=*/ 0, /*mutual=*/ 0, /*circular=*/1);
    igraph_vector_int_init_range(&perm, 0, igraph_vcount(&ring1));
    igraph_vector_int_shuffle(&perm);
    igraph_permute_vertices(&ring1, &ring2, &perm);

    /* Everything has the same colors */
    igraph_vector_int_init(&color1, igraph_vcount(&ring1));
    igraph_vector_int_init(&color2, igraph_vcount(&ring2));
    igraph_isomorphic_vf2(&ring1, &ring2, &color1, &color2, 0, 0, &iso, 0, 0, 0, 0, 0);
    if (!iso) {
        fprintf(stderr, "Single color failed.\n");
        return 1;
    }

    /* Two colors, just counting */
    for (i = 0; i < igraph_vector_int_size(&color1); i += 2) {
        VECTOR(color1)[i] = VECTOR(color2)[VECTOR(perm)[i]] = 1;
    }
    igraph_count_isomorphisms_vf2(&ring1, &ring2, &color1, &color2, 0, 0, &count, 0, 0, 0);
    if (count != 100) {
        fprintf(stderr, "Count with two colors failed, expected 100, got %" IGRAPH_PRId ".\n", count);
        return 2;
    }

    igraph_destroy(&ring1);
    igraph_destroy(&ring2);
    igraph_vector_int_destroy(&color2);
    igraph_vector_int_destroy(&perm);

    /* Two colors, count subisomorphisms */
    igraph_ring(&ring1, 100, /*directed=*/ 0, /*mutual=*/ 0, /*circular=*/0);
    igraph_ring(&ring2, 80, /*directed=*/ 0, /*mutual=*/ 0, /*circular=*/0);

    igraph_vector_int_init(&color2, igraph_vcount(&ring2));
    for (i = 0; i < igraph_vector_int_size(&color1); i += 2) {
        VECTOR(color1)[i]   = 0;
        VECTOR(color1)[i + 1] = 1;
    }
    for (i = 0; i < igraph_vector_int_size(&color2); i += 2) {
        VECTOR(color2)[i]   = 0;
        VECTOR(color2)[i + 1] = 1;
    }
    igraph_count_subisomorphisms_vf2(&ring1, &ring2, &color1, &color2, 0, 0,
                                     &count, 0, 0, 0);
    if (count != 21) {
        fprintf(stderr, "Count with two colors failed, expected 21, got %" IGRAPH_PRId ".\n", count);
        return 3;
    }

    igraph_vector_int_destroy(&color1);
    igraph_vector_int_destroy(&color2);

    igraph_destroy(&ring1);
    igraph_destroy(&ring2);

    return 0;
}
