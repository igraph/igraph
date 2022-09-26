/*
   IGraph library.
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

int main(void) {
    igraph_matrix_t a, b, c, d, e;

    igraph_matrix_init(&c, 0, 0);

    igraph_real_t elemA[] = {1, 2, 3, 4};
    igraph_real_t elemB[] = {5, 6, 7, 8};
    igraph_real_t elemD[] = {5, 6, 7, 8, 9, 10};
    matrix_init_real_row_major(&a, 2, 2, elemA);
    matrix_init_real_row_major(&b, 2, 2, elemB);
    matrix_init_real_row_major(&d, 2, 3, elemD);
    matrix_init_real_row_major(&e, 3, 2, elemD);

    printf("matrix multiplication, A={{1,2},{3,4}} B ={{5,6},{7,8}}\n");
    igraph_blas_dgemm(0, 0, 1, &a, &b, 0, &c);
    igraph_matrix_print(&c);

    printf("transpose a first\n");
    igraph_blas_dgemm(1, 0, 1, &a, &b, 0, &c);
    igraph_matrix_print(&c);

    printf("transpose b first\n");
    igraph_blas_dgemm(0, 1, 1, &a, &b, 0, &c);
    igraph_matrix_print(&c);

    printf("transpose both matrices first\n");
    igraph_blas_dgemm(1, 1, 1, &a, &b, 0, &c);
    igraph_matrix_print(&c);

    printf("multiply by 0.5\n");
    igraph_blas_dgemm(0, 0, 0.5, &a, &b, 0, &c);
    igraph_matrix_print(&c);

    printf("multiply by 1.5 by using previous result\n");
    igraph_blas_dgemm(0, 0, 1, &a, &b, 1, &c);
    igraph_matrix_print(&c);

    printf("matrix multiplication, A={{1,2},{3,4}} B={{5,6,7},{8,9,10}}\n");
    igraph_blas_dgemm(0, 0, 1, &a, &d, 0, &c);
    igraph_matrix_print(&c);

    printf("matrix multiplication, A={{5,8},{6,9},{7,10}} B={{1,2},{3,4}}\n");
    igraph_blas_dgemm(1, 0, 1, &d, &a, 0, &c);
    igraph_matrix_print(&c);

    printf("check error when matrix sizes don't match\n");
    CHECK_ERROR(igraph_blas_dgemm(0, 0, 1, &d, &a, 0, &c), IGRAPH_EINVAL);
    CHECK_ERROR(igraph_blas_dgemm(0, 1, 1, &a, &d, 0, &c), IGRAPH_EINVAL);
    CHECK_ERROR(igraph_blas_dgemm(0, 0, 1, &a, &b, 1, &d), IGRAPH_EINVAL);
    CHECK_ERROR(igraph_blas_dgemm(0, 0, 1, &a, &b, 1, &e), IGRAPH_EINVAL);

    igraph_matrix_destroy(&a);
    igraph_matrix_destroy(&b);
    igraph_matrix_destroy(&c);
    igraph_matrix_destroy(&d);
    igraph_matrix_destroy(&e);

    VERIFY_FINALLY_STACK();
    return 0;
}
