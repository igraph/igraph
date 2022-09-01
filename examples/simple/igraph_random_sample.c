/* -*- mode: C -*-  */
/*
  Test suite for random sampling.
  Copyright (C) 2011 Minh Van Nguyen <nguyenminh2@gmail.com>

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

#include <igraph.h>
#include <math.h>
#include <stdio.h>

#define R_INTEGER(a,b) (igraph_rng_get_integer(igraph_rng_default(), (a), (b)))

/* test parameters */
typedef struct {
    igraph_integer_t low;
    igraph_integer_t high;
    igraph_integer_t length;
    igraph_error_t retval;
} sampling_test_t;

/* Get a few random samples and test their properties. */
int random_sample_test() {
    const igraph_integer_t min = -1000;
    const igraph_integer_t max = 1000;
    igraph_integer_t low;       /* lower limit */
    igraph_integer_t high;      /* upper limit */
    igraph_integer_t length;    /* sample size */
    igraph_integer_t poolsize;  /* size of candidate pool */
    igraph_real_t sP;           /* population total sum */
    igraph_real_t ss;           /* sample total sum */
    igraph_vector_int_t V;
    igraph_integer_t i;

    igraph_rng_seed(igraph_rng_default(), 57); /* make tests deterministic */

    /* The generated sequence of numbers must be in increasing order. */
    igraph_vector_int_init(&V, /*size*/ 0);
    do {
        high = (igraph_integer_t)R_INTEGER(min, max);
    } while (high == min);
    do {
        low = (igraph_integer_t)R_INTEGER(min, max);
    } while (low >= high);
    poolsize = (igraph_integer_t)fabs((double)high - (double)low);
    do {
        length = (igraph_integer_t)R_INTEGER(1, max);
    } while (length > poolsize);
    igraph_random_sample(&V, low, high, length);
    if (length != igraph_vector_int_size(&V)) {
        printf("Requested vector length and resulting length mismatch.\n");
        return 1;
    }
    for (i = 0; i < length - 1; i++) {
        if (VECTOR(V)[i] >= VECTOR(V)[i + 1]) {
            printf("Sample not in increasing order.\n");
            return 1;
        }
    }
    igraph_vector_int_destroy(&V);

    /* Let P be a candidate pool of positive integers with total sum s_P. */
    /* Let S be a random sample from P and having total sum s_S. Then we */
    /* have the bound s_s <= s_P. */
    igraph_vector_int_init(&V, /*size*/ 0);
    low = 1;
    do {
        high = (igraph_integer_t)R_INTEGER(low, max);
    } while (high == low);
    poolsize = (igraph_integer_t)fabs((double)high - (double)low);
    do {
        length = (igraph_integer_t)R_INTEGER(low, max);
    } while (length > poolsize);
    igraph_random_sample(&V, low, high, length);
    /* Use Gauss' formula to sum all consecutive positive integers from 1 */
    /* up to and including an upper limit. In LaTeX, Gauss' formula is */
    /* \sum_{i=1}^n i = \frac{n(n+1)}{2} where n is the upper limit. */
    sP = (high * (high + 1)) / 2;
    ss = igraph_vector_int_sum(&V);
    if (ss > sP) {
        printf("Sum of sampled sequence exceeds sum of whole population.\n");
        return 1;
    }
    igraph_vector_int_destroy(&V);

    return 0;
}

int equal_test() {
    igraph_vector_int_t V;

    igraph_vector_int_init(&V, 0);

    igraph_random_sample(&V, 0, 0, 1);
    if (igraph_vector_int_size(&V) != 1) {
        return 1;
    }
    if (VECTOR(V)[0] != 0) {
        return 2;
    }

    igraph_random_sample(&V, 10, 10, 1);
    if (igraph_vector_int_size(&V) != 1) {
        return 3;
    }
    if (VECTOR(V)[0] != 10) {
        return 4;
    }

    igraph_random_sample(&V, 2, 12, 11);
    if (igraph_vector_int_size(&V) != 11) {
        return 5;
    }
    for (igraph_integer_t i = 0; i < 11; i++)
        if (VECTOR(V)[i] != i + 2) {
            return 6;
        }

    igraph_vector_int_destroy(&V);
    return 0;
}

int rare_test() {
    igraph_vector_int_t V;

    igraph_vector_int_init(&V, 0);

    igraph_random_sample(&V, 0, 0, 1);
    if (igraph_vector_int_size(&V) != 1) {
        return 1;
    }
    if (VECTOR(V)[0] != 0) {
        return 2;
    }

    igraph_random_sample(&V, 10, 10, 1);
    if (igraph_vector_int_size(&V) != 1) {
        return 3;
    }
    if (VECTOR(V)[0] != 10) {
        return 4;
    }

    igraph_vector_int_destroy(&V);
    return 0;
}

int main() {
    int ret;

    ret = random_sample_test();
    if (ret) {
        return 2;
    }
    ret = equal_test();
    if (ret) {
        return 3;
    }
    ret = rare_test();
    if (ret) {
        return 4;
    }

    return 0;
}
