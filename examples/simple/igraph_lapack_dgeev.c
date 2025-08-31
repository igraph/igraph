/*
   igraph library.
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
#include <stdio.h>

int main(void) {
    igraph_matrix_t A;
    igraph_matrix_t vectors_left, vectors_right;
    igraph_vector_t values_real, values_imag;
    igraph_vector_complex_t values;
    int info = 1;

    /* Initialize the library. */
    igraph_setup();

    igraph_matrix_init(&A, 2, 2);
    igraph_matrix_init(&vectors_left, 0, 0);
    igraph_matrix_init(&vectors_right, 0, 0);
    igraph_vector_init(&values_real, 0);
    igraph_vector_init(&values_imag, 0);
    MATRIX(A, 0, 0) = 1.0;
    MATRIX(A, 0, 1) = 1.0;
    MATRIX(A, 1, 0) = -1.0;
    MATRIX(A, 1, 1) = 1.0;

    igraph_lapack_dgeev(&A, &values_real, &values_imag,
                        &vectors_left, &vectors_right, &info);
    igraph_vector_complex_create(&values, &values_real, &values_imag);
    printf("eigenvalues:\n");
    igraph_vector_complex_print(&values);
    printf("left eigenvectors:\n");
    igraph_matrix_print(&vectors_left);
    printf("right eigenvectors:\n");
    igraph_matrix_print(&vectors_right);

    igraph_vector_destroy(&values_imag);
    igraph_vector_destroy(&values_real);
    igraph_vector_complex_destroy(&values);
    igraph_matrix_destroy(&vectors_right);
    igraph_matrix_destroy(&vectors_left);
    igraph_matrix_destroy(&A);

    return 0;
}
