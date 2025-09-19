/*
   igraph library.
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

#include <igraph.h>

#include "bench.h"

/* This program benchmarks igraph_qsort() indirectly through vector_sort() */

int main(void) {
    igraph_vector_int_t vec;

    igraph_rng_seed(igraph_rng_default(), 137);
    BENCH_INIT();


    igraph_vector_int_init(&vec, 0);

#define N 10000000

    igraph_vector_int_resize(&vec, N);
    for (igraph_int_t i=0; i < N; i++) {
        VECTOR(vec)[i] = RNG_INTEGER(0, N-1);
    }
    BENCH("Sort vector of length " IGRAPH_I_STRINGIFY(N), igraph_vector_int_sort(&vec));

#undef N
#define N 1000000

    igraph_vector_int_resize(&vec, N);
    for (igraph_int_t i=0; i < N; i++) {
        VECTOR(vec)[i] = RNG_INTEGER(0, N-1);
    }
    BENCH("Sort vector of length " IGRAPH_I_STRINGIFY(N), igraph_vector_int_sort(&vec));

#undef N
#define N 100000

    igraph_vector_int_resize(&vec, N);
    for (igraph_int_t i=0; i < N; i++) {
        VECTOR(vec)[i] = RNG_INTEGER(0, N-1);
    }
    BENCH("Sort vector of length " IGRAPH_I_STRINGIFY(N), igraph_vector_int_sort(&vec));

    igraph_vector_int_destroy(&vec);

    return 0;
}
