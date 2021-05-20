/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2011-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#include "test_utilities.inc"

#define DUMP() do {             \
        igraph_vector_complex_print(&values);   \
        igraph_vector_complex_print(&values2);  \
    } while(0)

int main() {

    const int nodes = 10, skip = 3;
    igraph_matrix_t mat2;
    igraph_vector_complex_t values, values2;
    igraph_matrix_complex_t vectors, vectors2;
    igraph_eigen_which_t which;
    int i;

    igraph_rng_seed(igraph_rng_default(), 42);
    igraph_matrix_init(&mat2, nodes, nodes);
    for (i = 0; i < nodes; i++) {
        int j;
        for (j = 0; j < nodes; j++) {
            MATRIX(mat2, i, j) = igraph_rng_get_integer(igraph_rng_default(), 1, 10);
        }
    }

    which.pos = IGRAPH_EIGEN_SELECT;
    which.il = skip;
    which.iu = nodes - skip;

    igraph_vector_complex_init(&values, 0);
    igraph_matrix_complex_init(&vectors, 0, 0);
    igraph_eigen_matrix(&mat2, /*sparsemat=*/ 0, /*fun=*/ 0, nodes,
                        /*extra=*/ 0, IGRAPH_EIGEN_LAPACK, &which,
                        /*options=*/ 0, /*storage=*/ 0, &values, &vectors);

    which.pos = IGRAPH_EIGEN_ALL;

    igraph_vector_complex_init(&values2, 0);
    igraph_matrix_complex_init(&vectors2, 0, 0);
    igraph_eigen_matrix(&mat2, /*sparsemat=*/ 0, /*fun=*/ 0, nodes,
                        /*extra=*/ 0, IGRAPH_EIGEN_LAPACK, &which,
                        /*options=*/ 0, /*storage=*/ 0, &values2, &vectors2);

    for (i = 0; i < nodes - skip * 2 + 1; i++) {
        int j;
        igraph_real_t d =
            igraph_complex_abs(igraph_complex_sub(VECTOR(values)[i],
                               VECTOR(values2)[i + skip - 1]));
        if (d > 1e-15) {
            DUMP();
            return 2;
        }
        for (j = 0; j < nodes; j++) {
            igraph_real_t d =
                igraph_complex_abs(igraph_complex_sub(MATRIX(vectors, j, i),
                                   MATRIX(vectors2, j, i + skip - 1)));
            if (d > 1e-15) {
                DUMP();
                return 3;
            }
        }
    }

    igraph_vector_complex_destroy(&values);
    igraph_matrix_complex_destroy(&vectors);
    igraph_vector_complex_destroy(&values2);
    igraph_matrix_complex_destroy(&vectors2);
    igraph_matrix_destroy(&mat2);

    VERIFY_FINALLY_STACK();

    return 0;
}
