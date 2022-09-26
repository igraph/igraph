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
#include <stdio.h>

void matrix_to_complex_vectors(igraph_vector_complex_t *c1, igraph_vector_complex_t *c2, igraph_matrix_t *m) {
    IGRAPH_REAL(VECTOR(*c1)[0]) = MATRIX(*m, 0, 0);
    IGRAPH_REAL(VECTOR(*c1)[1]) = MATRIX(*m, 1, 0);
    IGRAPH_IMAG(VECTOR(*c1)[0]) = MATRIX(*m, 0, 1);
    IGRAPH_IMAG(VECTOR(*c1)[1]) = MATRIX(*m, 1, 1);

    IGRAPH_REAL(VECTOR(*c2)[0]) = MATRIX(*m, 0, 0);
    IGRAPH_REAL(VECTOR(*c2)[1]) = MATRIX(*m, 1, 0);
    IGRAPH_IMAG(VECTOR(*c2)[0]) = -MATRIX(*m, 0, 1);
    IGRAPH_IMAG(VECTOR(*c2)[1]) = -MATRIX(*m, 1, 1);
}

int main(void) {

    igraph_matrix_t A;
    igraph_matrix_t vectors_left, vectors_right;
    igraph_vector_t values_real, values_imag;
    igraph_vector_complex_t values;
    igraph_vector_complex_t eigenvector1;
    igraph_vector_complex_t eigenvector2;
    int info = 1;
    igraph_real_t abnrm;

    igraph_matrix_init(&A, 2, 2);
    igraph_matrix_init(&vectors_left, 0, 0);
    igraph_matrix_init(&vectors_right, 0, 0);
    igraph_vector_init(&values_real, 0);
    igraph_vector_init(&values_imag, 0);
    igraph_vector_complex_init(&eigenvector1, 2);
    igraph_vector_complex_init(&eigenvector2, 2);
    MATRIX(A, 0, 0) = 1.0;
    MATRIX(A, 0, 1) = 1.0;
    MATRIX(A, 1, 0) = -1.0;
    MATRIX(A, 1, 1) = 1.0;

    igraph_lapack_dgeevx(IGRAPH_LAPACK_DGEEVX_BALANCE_BOTH,
                         &A, &values_real, &values_imag,
                         &vectors_left, &vectors_right, NULL, NULL,
                         /*scale=*/ NULL, &abnrm, /*rconde=*/ NULL,
                         /*rcondv=*/ NULL, &info);

    igraph_vector_complex_create(&values, &values_real, &values_imag);
    printf("eigenvalues:\n");
    igraph_vector_complex_print(&values);

    printf("\nleft eigenvectors:\n");
    /*matrix_to_complex_vectors only works because we have two complex
      conjugate eigenvalues */
    matrix_to_complex_vectors(&eigenvector1, &eigenvector2, &vectors_left);
    igraph_vector_complex_print(&eigenvector1);
    igraph_vector_complex_print(&eigenvector2);

    printf("\nright eigenvectors:\n");
    matrix_to_complex_vectors(&eigenvector1, &eigenvector2, &vectors_right);
    igraph_vector_complex_print(&eigenvector1);
    igraph_vector_complex_print(&eigenvector2);
    printf("\nOne-norm of the balanced matrix:\n%g\n", abnrm);

    igraph_vector_destroy(&values_imag);
    igraph_vector_destroy(&values_real);
    igraph_vector_complex_destroy(&values);
    igraph_vector_complex_destroy(&eigenvector1);
    igraph_vector_complex_destroy(&eigenvector2);
    igraph_matrix_destroy(&vectors_right);
    igraph_matrix_destroy(&vectors_left);
    igraph_matrix_destroy(&A);

    return 0;
}
