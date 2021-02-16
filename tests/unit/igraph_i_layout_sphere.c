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
#include <math.h>

#include "layout/layout_internal.h"

#include "test_utilities.inc"

int main () {
    long int i;
    igraph_matrix_t m;
    igraph_real_t x, y, z, r;

    srand(42); /* make tests deterministic */

    /* 2D */
    igraph_matrix_init(&m, 1000, 2);
    for (i = 0; i < igraph_matrix_nrow(&m); i++) {
        MATRIX(m, i, 0) = rand() / (double)RAND_MAX;
        MATRIX(m, i, 1) = rand() / (double)RAND_MAX;
    }
    igraph_i_layout_sphere_2d(&m, &x, &y, &r);

    for (i = 0; i < igraph_matrix_nrow(&m); i++) {
        igraph_real_t dist = sqrt((MATRIX(m, i, 0) - x) * (MATRIX(m, i, 0) - x) +
                                  (MATRIX(m, i, 1) - y) * (MATRIX(m, i, 1) - y));
        if (dist > r) {
            printf("x: %f y: %f r: %f\n", x, y, r);
            printf("x: %f y: %f dist: %f (%li)\n",
                   MATRIX(m, i, 0), MATRIX(m, i, 1), dist, i);
            return 1;
        }
    }
    igraph_matrix_destroy(&m);

    /* 3D */
    igraph_matrix_init(&m, 1000, 3);
    for (i = 0; i < igraph_matrix_nrow(&m); i++) {
        MATRIX(m, i, 0) = rand() / (double)RAND_MAX;
        MATRIX(m, i, 1) = rand() / (double)RAND_MAX;
        MATRIX(m, i, 2) = rand() / (double)RAND_MAX;
    }
    igraph_i_layout_sphere_3d(&m, &x, &y, &z, &r);

    for (i = 0; i < igraph_matrix_nrow(&m); i++) {
        igraph_real_t dist = sqrt((MATRIX(m, i, 0) - x) * (MATRIX(m, i, 0) - x) +
                                  (MATRIX(m, i, 1) - y) * (MATRIX(m, i, 1) - y) +
                                  (MATRIX(m, i, 2) - z) * (MATRIX(m, i, 2) - z));
        if (dist > r) {
            printf("x: %f y: %f z: %f r: %f\n", x, y, z, r);
            printf("x: %f y: %f z: %f dist: %f (%li)\n",
                   MATRIX(m, i, 0), MATRIX(m, i, 1), MATRIX(m, i, 2), dist, i);
            return 1;
        }
    }
    igraph_matrix_destroy(&m);

    VERIFY_FINALLY_STACK();

    return 0;
}
