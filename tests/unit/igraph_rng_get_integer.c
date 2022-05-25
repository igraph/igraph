/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2011-2021  The igraph development team

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

void simple_tests() {
    int i;

    /* Seed the RNG, generate 10 random integers */
    igraph_rng_seed(igraph_rng_default(), 42);
    for (i = 0; i < 10; i++) {
        printf("%" IGRAPH_PRId "\n", igraph_rng_get_integer(igraph_rng_default(), 10, 100));
    }

    printf("========\n");

    /* Seed the RNG again with the same seed, verify that we get the same
     * numbers */
    igraph_rng_seed(igraph_rng_default(), 42);
    for (i = 0; i < 10; i++) {
        printf("%" IGRAPH_PRId "\n", igraph_rng_get_integer(igraph_rng_default(), 10, 100));
    }

    printf("========\n");

    /* Seed the RNG again with a different seed, verify that we get different
     * numbers */
    igraph_rng_seed(igraph_rng_default(), 84);
    for (i = 0; i < 10; i++) {
        printf("%" IGRAPH_PRId "\n", igraph_rng_get_integer(igraph_rng_default(), 10, 100));
    }

    printf("========\n");
}

void generate_random_vector(
    igraph_rng_t* rng, igraph_vector_int_t* numbers, igraph_integer_t lo, igraph_integer_t hi
) {
    igraph_integer_t i, n = igraph_vector_int_size(numbers);

    for (i = 0; i < n; i++) {
        VECTOR(*numbers)[i] = igraph_rng_get_integer(rng, lo, hi);
    }
    IGRAPH_ASSERT(igraph_vector_int_min(numbers) >= lo);
    IGRAPH_ASSERT(igraph_vector_int_max(numbers) <= hi);
}

void check_occurrences(const igraph_vector_int_t* numbers, igraph_integer_t lo, igraph_integer_t hi) {
    igraph_integer_t i;

    for (i = lo; i <= hi; i++) {
        IGRAPH_ASSERT(igraph_vector_int_contains(numbers, i));
    }
}

void stress_tests() {
    igraph_rng_t rng;
    const igraph_rng_type_t* rng_types[] = {
        &igraph_rngtype_mt19937,
        &igraph_rngtype_glibc2,
    };
    igraph_integer_t i;
    igraph_vector_int_t numbers;
    const igraph_integer_t N = 1000;

    igraph_vector_int_init(&numbers, N);

    for (i = 0; i < sizeof(rng_types) / sizeof(rng_types[0]); i++) {
        igraph_rng_init(&rng, rng_types[i]);
        igraph_rng_seed(&rng, 42);

        /* We are going to test multiple ranges. In each range, we generate
         * 1000 random numbers and test whether all the values are in the
         * specified range. For the small ranges, we also test whether each
         * value occurred at least once as the chances of this not happening is
         * astronomically small */

        /* Test integer generation in a small 0-based range */
        generate_random_vector(&rng, &numbers, 0, 5);
        check_occurrences(&numbers, 0, 5);

        /* Test integer generation in a small non-0-based range */
        generate_random_vector(&rng, &numbers, 8, 13);
        check_occurrences(&numbers, 8, 13);

        /* Test integer generation in [0; IGRAPH_INTEGER_MAX] */
        generate_random_vector(&rng, &numbers, 0, IGRAPH_INTEGER_MAX);

        /* Test integer generation in [-5; IGRAPH_INTEGER_MAX-5] */
        generate_random_vector(&rng, &numbers, -5, IGRAPH_INTEGER_MAX-5);

        /* Test integer generation in [-5; IGRAPH_INTEGER_MAX] */
        generate_random_vector(&rng, &numbers, -5, IGRAPH_INTEGER_MAX);

        /* Test integer generation in [IGRAPH_INTEGER_MIN; -5] */
        generate_random_vector(&rng, &numbers, IGRAPH_INTEGER_MIN, -5);

        /* Test integer generation in [IGRAPH_INTEGER_MIN; IGRAPH_INTEGER_MAX] */
        generate_random_vector(&rng, &numbers, IGRAPH_INTEGER_MIN, IGRAPH_INTEGER_MAX);

        igraph_rng_destroy(&rng);
    }

    igraph_vector_int_destroy(&numbers);
}

int main() {
    simple_tests();
    stress_tests();

    return 0;
}
