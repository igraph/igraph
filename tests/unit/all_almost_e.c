/*
   IGraph library.
   Copyright (C) 2022  The igraph development team <igraph@igraph.org>

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

#include <math.h>
#include <float.h>

#include "test_utilities.h"

int main(void) {
    igraph_matrix_t rm1, rm2;
    igraph_matrix_complex_t cm1, cm2;
    const igraph_integer_t nrow = 50, ncol = 60;

    igraph_rng_seed(igraph_rng_default(), 97);

    /* Real matrices */

    igraph_matrix_init(&rm1, nrow, ncol);
    for (igraph_integer_t i=0; i < nrow; i++) {
        for (igraph_integer_t j=0; j < ncol; j++) {
            MATRIX(rm1, i, j) = RNG_UNIF(0.0, 3.0);
        }
    }

    igraph_matrix_init_copy(&rm2, &rm1);
    for (igraph_integer_t i=0; i < nrow; i++) {
        for (igraph_integer_t j=0; j < ncol; j++) {
            MATRIX(rm2, i, j) = pow(MATRIX(rm2, i, j), 1 / 7.0);
            MATRIX(rm2, i, j) = pow(MATRIX(rm2, i, j), 7.0);
        }
    }

    IGRAPH_ASSERT(igraph_matrix_all_almost_e(&rm1, &rm2, 4*DBL_EPSILON));
    MATRIX(rm2, 0, 0) *= 2;
    IGRAPH_ASSERT(! igraph_matrix_all_almost_e(&rm1, &rm2, 4*DBL_EPSILON));

    igraph_matrix_destroy(&rm2);
    igraph_matrix_destroy(&rm1);

    /* Complex matrices */

    igraph_matrix_complex_init(&cm1, nrow, ncol);

    for (igraph_integer_t i=0; i < nrow; i++) {
        for (igraph_integer_t j=0; j < ncol; j++) {
            IGRAPH_REAL(MATRIX(cm1, i, j)) = RNG_NORMAL(0,1);
            IGRAPH_IMAG(MATRIX(cm1, i, j)) = RNG_NORMAL(0,1);
        }
    }


    igraph_matrix_complex_init_copy(&cm2, &cm1);
    for (igraph_integer_t i=0; i < nrow; i++) {
        for (igraph_integer_t j=0; j < ncol; j++) {
            MATRIX(cm2, i, j) = igraph_complex_pow_real(MATRIX(cm2, i, j), 1 / 7.0);
            MATRIX(cm2, i, j) = igraph_complex_pow_real(MATRIX(cm2, i, j), 7.0);
        }
    }

    IGRAPH_ASSERT(igraph_matrix_complex_all_almost_e(&cm1, &cm2, 8*DBL_EPSILON));
    MATRIX(cm2, 0, 0) = igraph_complex_mul(MATRIX(cm2, 0, 0), MATRIX(cm2, 1, 1));
    IGRAPH_ASSERT(! igraph_matrix_complex_all_almost_e(&cm1, &cm2, 8*DBL_EPSILON));

    igraph_matrix_complex_destroy(&cm2);
    igraph_matrix_complex_destroy(&cm1);

    VERIFY_FINALLY_STACK();

    return 0;
}
