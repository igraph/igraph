/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2010-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#define DIM1 10
#define DIM2 5

#define INT(a) (igraph_rng_get_integer(igraph_rng_default(), 0, (a)))

int main() {
    igraph_matrix_t mat;
    igraph_sparsemat_t spmat, spmat2;
    int i;
    igraph_real_t m1, m2;

    igraph_rng_seed(igraph_rng_default(), 42);

    igraph_sparsemat_init(&spmat, DIM1, DIM2, 20);
    igraph_sparsemat_entry(&spmat, 1, 2, -1.0);
    igraph_sparsemat_entry(&spmat, 3, 2, 10.0);
    for (i = 0; i < 10; i++) {
        igraph_sparsemat_entry(&spmat, INT(DIM1 - 1), INT(DIM2 - 1), 1.0);
    }
    igraph_sparsemat_entry(&spmat, 1, 2, -1.0);
    igraph_sparsemat_entry(&spmat, 3, 2, 10.0);

    igraph_sparsemat_compress(&spmat, &spmat2);
    igraph_matrix_init(&mat, 0, 0);
    igraph_sparsemat_as_matrix(&mat, &spmat2);
    m1 = igraph_sparsemat_min(&spmat2);
    m2 = igraph_matrix_min(&mat);
    if (m1 != m2) {
        printf("%f %f\n", m1, m2);
        return 1;
    }
    m1 = igraph_sparsemat_max(&spmat2);
    m2 = igraph_matrix_max(&mat);
    if (m1 != m2) {
        printf("%f %f\n", m1, m2);
        return 2;
    }

    igraph_sparsemat_minmax(&spmat2, &m1, &m2);
    if (m1 != igraph_matrix_min(&mat)) {
        return 3;
    }
    if (m2 != igraph_matrix_max(&mat)) {
        return 4;
    }

    igraph_matrix_destroy(&mat);
    igraph_sparsemat_destroy(&spmat);
    igraph_sparsemat_destroy(&spmat2);

    return 0;
}




