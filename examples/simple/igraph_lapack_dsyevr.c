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

int main(void) {

    igraph_matrix_t A;
    igraph_matrix_t vectors;
    igraph_vector_t values;

    igraph_matrix_init(&A, 2, 2);
    igraph_matrix_init(&vectors, 0, 0);
    igraph_vector_init(&values, 0);

    MATRIX(A, 0, 0) = 2.0;
    MATRIX(A, 0, 1) = -1.0;
    MATRIX(A, 1, 0) = -1.0;
    MATRIX(A, 1, 1) = 3.0;

    printf("Take a subset:\n");

    igraph_lapack_dsyevr(&A, IGRAPH_LAPACK_DSYEV_SELECT, /*vl=*/ 0, /*vu=*/ 0,
                         /*vestimate=*/ 0, /*il=*/ 1, /*iu=*/ 1,
                         /*abstol=*/ 1e-10, &values, &vectors,
                         /*support=*/ 0);
    printf("eigenvalues:\n");
    igraph_vector_print(&values);
    printf("eigenvectors:\n");
    igraph_matrix_print(&vectors);

    printf("\nTake a subset based on an interval:\n");

    igraph_lapack_dsyevr(&A, IGRAPH_LAPACK_DSYEV_INTERVAL, /*vl*/ 3, /*vu*/ 4,
                         /*vestimate=*/ 1, /*il=*/ 0, /*iu=*/ 0,
                         /*abstol=*/ 1e-10, &values, &vectors,
                         /*support=*/ 0);

    printf("eigenvalues:\n");
    igraph_vector_print(&values);
    printf("eigenvectors:\n");
    igraph_matrix_print(&vectors);

    igraph_vector_destroy(&values);
    igraph_matrix_destroy(&vectors);
    igraph_matrix_destroy(&A);

    return 0;
}
