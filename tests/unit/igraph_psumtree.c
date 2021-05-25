
/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
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
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include <igraph.h>
#include <stdlib.h>

#include "test_utilities.inc"

int main() {
    igraph_psumtree_t tree;
    igraph_vector_t vec;
    long int i;
    igraph_real_t sum;

    /* Uniform random numbers */
    igraph_vector_init(&vec, 16);
    igraph_psumtree_init(&tree, 16);
    sum = igraph_psumtree_sum(&tree);
    if (sum != 0) {
        printf("Sum: %f instead of 0.\n", sum);
        return 1;
    }

    for (i = 0; i < 16; i++) {
        igraph_psumtree_update(&tree, i, 1);
    }
    if ((sum = igraph_psumtree_sum(&tree)) != 16) {
        printf("Sum: %f instead of 16.\n", sum);
        return 2;
    }

    for (i = 0; i < 16000; i++) {
        igraph_real_t r = ((double)rand()) / RAND_MAX * sum;
        long int idx;
        igraph_psumtree_search(&tree, &idx, r);
        VECTOR(vec)[idx] += 1;
    }
    for (i = 0; i < 16; i++) {
        if (VECTOR(vec)[i] < 800 || VECTOR(vec)[i] > 1200) {
            return 3;
        }
    }

    /* Nonuniform, even indices have twice as much chance */
    for (i = 0; i < 16; i += 2) {
        igraph_psumtree_update(&tree, i, 2);
    }
    if ((sum = igraph_psumtree_sum(&tree)) != 24) {
        printf("Sum: %f instead of 24.\n", sum);
        return 4;
    }

    igraph_vector_null(&vec);
    for (i = 0; i < 24000; i++) {
        igraph_real_t r = ((double)rand()) / RAND_MAX * sum;
        long int idx;
        igraph_psumtree_search(&tree, &idx, r);
        VECTOR(vec)[idx] += 1;
    }
    for (i = 0; i < 16; i++) {
        if (i % 2 == 0 && (VECTOR(vec)[i] < 1800 || VECTOR(vec)[i] > 2200)) {
            return 5;
        }
        if (i % 2 != 0 && (VECTOR(vec)[i] < 800 || VECTOR(vec)[i] > 1200)) {
            return 6;
        }
    }

    /* Test zero probabilities */
    igraph_psumtree_update(&tree, 0, 0);
    igraph_psumtree_update(&tree, 5, 0);
    igraph_psumtree_update(&tree, 15, 0);
    sum = igraph_psumtree_sum(&tree);

    igraph_vector_null(&vec);
    for (i = 0; i < 20000; i++) {
        igraph_real_t r = ((double)rand()) / RAND_MAX * sum;
        long int idx;
        igraph_psumtree_search(&tree, &idx, r);
        VECTOR(vec)[idx] += 1;
    }
    if (VECTOR(vec)[0] != 0 || VECTOR(vec)[5] != 0 || VECTOR(vec)[15] != 0) {
        return 7;
    }

    igraph_vector_destroy(&vec);
    igraph_psumtree_destroy(&tree);

    /****************************************************/
    /* Non power-of-two vector size                     */
    /****************************************************/

    igraph_vector_init(&vec, 9);
    igraph_psumtree_init(&tree, 9);

    for (i = 0; i < 9; i++) {
        igraph_psumtree_update(&tree, i, 1);
    }
    sum = igraph_psumtree_sum(&tree);

    for (i = 0; i < 9000; i++) {
        igraph_real_t r = ((double)rand()) / RAND_MAX * sum;
        long int idx;
        igraph_psumtree_search(&tree, &idx, r);
        VECTOR(vec)[idx] += 1;
    }
    for (i = 0; i < 9; i++) {
        if (VECTOR(vec)[i] < 800 || VECTOR(vec)[i] > 1200) {
            return 8;
        }
    }

    /* Nonuniform, even indices have twice as much chance */
    for (i = 0; i < 9; i += 2) {
        igraph_psumtree_update(&tree, i, 2);
    }
    sum = igraph_psumtree_sum(&tree);

    igraph_vector_null(&vec);
    for (i = 0; i < 14000; i++) {
        igraph_real_t r = ((double)rand()) / RAND_MAX * sum;
        long int idx;
        igraph_psumtree_search(&tree, &idx, r);
        VECTOR(vec)[idx] += 1;
    }
    for (i = 0; i < 9; i++) {
        if (i % 2 == 0 && (VECTOR(vec)[i] < 1800 || VECTOR(vec)[i] > 2200)) {
            return 9;
        }
        if (i % 2 != 0 && (VECTOR(vec)[i] < 800 || VECTOR(vec)[i] > 1200)) {
            return 10;
        }
    }

    /* Test query */
    for (i = 0; i < igraph_psumtree_size(&tree); i++) {
        if (i % 2 == 0 && igraph_psumtree_get(&tree, i) != 2) {
            return 11;
        }
        if (i % 2 != 0 && igraph_psumtree_get(&tree, i) != 1) {
            return 12;
        }
    }

    /* Test zero probabilities */
    igraph_psumtree_update(&tree, 0, 0);
    igraph_psumtree_update(&tree, 5, 0);
    igraph_psumtree_update(&tree, 8, 0);
    sum = igraph_psumtree_sum(&tree);

    igraph_vector_null(&vec);
    for (i = 0; i < 9000; i++) {
        igraph_real_t r = ((double)rand()) / RAND_MAX * sum;
        long int idx;
        igraph_psumtree_search(&tree, &idx, r);
        VECTOR(vec)[idx] += 1;
    }
    if (VECTOR(vec)[0] != 0 || VECTOR(vec)[5] != 0 || VECTOR(vec)[8] != 0) {
        return 11;
    }

    igraph_vector_destroy(&vec);
    igraph_psumtree_destroy(&tree);

    /****************************************************/
    /* Error handling                                   */
    /****************************************************/

    igraph_psumtree_init(&tree, 9);
    if (igraph_psumtree_update(&tree, 2, -2) == IGRAPH_SUCCESS) {
        return 12;
    }
    if (igraph_psumtree_update(&tree, 2, -INFINITY) == IGRAPH_SUCCESS) {
        return 13;
    }
    if (igraph_psumtree_update(&tree, 2, IGRAPH_NAN) == IGRAPH_SUCCESS) {
        return 14;
    }
    igraph_psumtree_destroy(&tree);

    VERIFY_FINALLY_STACK();

    return 0;
}
