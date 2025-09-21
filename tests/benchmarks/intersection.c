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

void rand_vec(igraph_vector_int_t *v, igraph_int_t n, igraph_int_t k) {
    igraph_vector_int_resize(v, n);
    for (igraph_int_t i=0; i < n; i++) {
        VECTOR(*v)[i] = RNG_INTEGER(0, k);
    }
}

void run_bench(int i, int n, int r) {
    igraph_vector_int_t a, b;
    igraph_int_t na = n, nb = r*n;
    int rep = 300000000 / nb;
    char msg[255];

    igraph_vector_int_init(&a, na);
    igraph_vector_int_init(&b, nb);

    rand_vec(&a, na, nb);
    igraph_vector_int_sort(&a);
    rand_vec(&b, nb, nb);
    igraph_vector_int_sort(&b);

    snprintf(msg, sizeof(msg) / sizeof(msg[0]),
             "%2d n = %5d, r = %3d, %dx", i, n, r, rep);

    /* 'volatile' is needed to prevent the compiler from optimizing
     * away multiple calls to the function with the same parameters within REPEAT.
     * This would normally happen due to the use of IGRAPH_FUNCATTR_PURE.
     * ATTENTION! 'volatile', when used this way, may not prevent this optimization
     * with future compiler versions. */
    volatile igraph_int_t res;
    BENCH(msg, REPEAT(res = igraph_vector_int_intersection_size_sorted(&a, &b), rep));
    (void) res;  /* silence unused-but-set-variable warning */

    igraph_vector_int_destroy(&a);
    igraph_vector_int_destroy(&b);
}

int main(void) {
    int i = 0;

    igraph_rng_seed(igraph_rng_default(), 137);
    BENCH_INIT();

#define BENCHSET(n) \
    run_bench(++i, n, 1); \
    run_bench(++i, n, 3); \
    run_bench(++i, n, 10); \
    run_bench(++i, n, 30); \
    run_bench(++i, n, 100); \
    printf("\n");

    BENCHSET(1);
    BENCHSET(3);
    BENCHSET(10);
    BENCHSET(30);
    BENCHSET(100);
    BENCHSET(300);
    BENCHSET(1000);
    BENCHSET(3000);
    BENCHSET(10000);
    BENCHSET(30000);
    BENCHSET(100000);

    return 0;
}
