/*
   igraph library.
   Copyright (C) 2021  The igraph development team <igraph@igraph.org>

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
#include "test_utilities.h"

void check(igraph_bool_t volume, igraph_int_t dim, igraph_int_t n, igraph_real_t radius, igraph_bool_t positive) {
    igraph_matrix_t samples;

    igraph_matrix_init(&samples, 0, 0);
    if (volume) {
        igraph_rng_sample_sphere_volume(igraph_rng_default(), dim, n, radius, positive, &samples);
    } else {
        igraph_rng_sample_sphere_surface(igraph_rng_default(), dim, n, radius, positive, &samples);
    }
    IGRAPH_ASSERT(igraph_matrix_ncol(&samples) == n);
    IGRAPH_ASSERT(igraph_matrix_nrow(&samples) == dim);
    for (igraph_int_t col = 0; col < n; col++) {
        igraph_real_t sum = 0;
        for (igraph_int_t row = 0; row < dim; row++) {
            if (positive) {
                IGRAPH_ASSERT(MATRIX(samples, row, col) >= 0);
            }
            sum += MATRIX(samples, row, col) * MATRIX(samples, row, col);
        }
        if (volume) {
            IGRAPH_ASSERT(sum <= radius * radius);
        } else {
            IGRAPH_ASSERT(igraph_almost_equals(sum, radius * radius, 0.00001));
        }
    }
    igraph_matrix_destroy(&samples);
}

int main(void) {

    igraph_rng_seed(igraph_rng_default(), 42);

    //No samples
    check(0, 2, 0, 1, 0);
    check(1, 2, 0, 1, 0);
    //Five samples, four-dimensions, radius 2
    check(0, 4, 5, 2, 0);
    check(1, 4, 5, 2, 0);
    //Same, positive orthant
    check(0, 4, 5, 2, 1);
    check(1, 4, 5, 2, 1);

    VERIFY_FINALLY_STACK();
    return 0;
}
