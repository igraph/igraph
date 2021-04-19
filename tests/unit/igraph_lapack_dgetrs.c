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
#include "test_utilities.inc"

void check_and_destroy(igraph_matrix_t *a, igraph_matrix_t *b, igraph_vector_int_t *ipiv,
                       igraph_bool_t transpose, igraph_error_type_t error) {
    printf("LU matrix:\n");
    igraph_matrix_print(a);
    printf("Pivot vector:\n");
    igraph_vector_int_print(ipiv);
    printf("B matrix:\n");
    igraph_matrix_print(b);
    IGRAPH_ASSERT(igraph_lapack_dgetrs(transpose, a, ipiv, b) == error);
    if (error == IGRAPH_SUCCESS) {
        printf("Returned matrix:\n");
        igraph_matrix_print(b);
    }
    igraph_vector_int_destroy(ipiv);
    igraph_matrix_destroy(a);
    igraph_matrix_destroy(b);
    printf("\n");
}

int main() {
    igraph_matrix_t a, b;
    igraph_vector_int_t ipiv;

    igraph_set_error_handler(igraph_error_handler_ignore);

    printf("Checking empty matrices:\n");
    igraph_matrix_init(&a, 0, 0);
    igraph_matrix_init(&b, 0, 0);
    igraph_vector_int_init_int(&ipiv, 0);
    check_and_destroy(&a, &b, &ipiv, 0, IGRAPH_SUCCESS);

    {
        printf("Checking 3x3 matrix:\n");
        double a_elements[] = {7, 8, 9, 2./7., -2./7., 3./7., 1./7., 1./2., -1./2.};
        int b_elements[] = {1, 1, 1};
        matrix_init_real_row_major(&a, 3, 3, a_elements);
        matrix_init_int_row_major(&b, 3, 1, b_elements);
        igraph_vector_int_init_int(&ipiv, 3, 1, 2, 3);
        check_and_destroy(&a, &b, &ipiv, 0, IGRAPH_SUCCESS);
    }
    {
        printf("Checking transpose and pivot:\n");
        double a_elements[] = {9, 3, 1, 8./9., -2./3., 1./9., 7./9., 1./2., 1./6.};
        int b_elements[] = {1, 1, 1};
        matrix_init_real_row_major(&a, 3, 3, a_elements);
        matrix_init_int_row_major(&b, 3, 1, b_elements);
        igraph_vector_int_init_int(&ipiv, 3, 3, 2, 3);
        check_and_destroy(&a, &b, &ipiv, 1, IGRAPH_SUCCESS);
    }
    {
        printf("Checking 2x3 matrix, expected to fail:\n");
        double a_elements[] = {4, 5, 6, 1./4., 3./4., 3./2.};
        int b_elements[] = {1, 1};
        matrix_init_real_row_major(&a, 2, 3, a_elements);
        matrix_init_int_row_major(&b, 2, 1, b_elements);
        igraph_vector_int_init_int(&ipiv, 2, 2, 2);
        check_and_destroy(&a, &b, &ipiv, 0, IGRAPH_NONSQUARE);
    }
    {
        printf("Checking singular matrix, this gives random output, so we just check for memory problems.\n\n");
        double a_elements[] = {6, 8, 7, 0, 1, 2, 0.5, 0.5, 0};
        int b_elements[] = {1, 1, 1};
        matrix_init_real_row_major(&a, 3, 3, a_elements);
        matrix_init_int_row_major(&b, 3, 1, b_elements);
        igraph_vector_int_init_int(&ipiv, 3, 1, 2, 3);
        IGRAPH_ASSERT(igraph_lapack_dgetrs(0, &a, &ipiv, &b) == IGRAPH_SUCCESS);
        igraph_vector_int_destroy(&ipiv);
        igraph_matrix_destroy(&a);
        igraph_matrix_destroy(&b);
    }
    {
        printf("Checking wrong size of B matrix, should fail:\n");
        double a_elements[] = {7, 8, 9, 2./7., -2./7., 3./7., 1./7., 1./2., -1./2.};
        int b_elements[] = {1, 1};
        matrix_init_real_row_major(&a, 3, 3, a_elements);
        matrix_init_int_row_major(&b, 2, 1, b_elements);
        igraph_vector_int_init_int(&ipiv, 3, 1, 2, 3);
        check_and_destroy(&a, &b, &ipiv, 0, IGRAPH_EINVAL);
    }
    {
        printf("Checking nonexisting pivots, should fail:\n");
        double a_elements[] = {7, 8, 9, 2./7., -2./7., 3./7., 1./7., 1./2., -1./2.};
        int b_elements[] = {1, 1, 1};
        matrix_init_real_row_major(&a, 3, 3, a_elements);
        matrix_init_int_row_major(&b, 3, 1, b_elements);
        igraph_vector_int_init_int(&ipiv, 3, 5, 6, 7);
        check_and_destroy(&a, &b, &ipiv, 0, IGRAPH_EINVAL);
    }
    {
        printf("Checking too few pivots, should fail:\n");
        double a_elements[] = {7, 8, 9, 2./7., -2./7., 3./7., 1./7., 1./2., -1./2.};
        int b_elements[] = {1, 1, 1};
        matrix_init_real_row_major(&a, 3, 3, a_elements);
        matrix_init_int_row_major(&b, 3, 1, b_elements);
        igraph_vector_int_init_int(&ipiv, 0);
        check_and_destroy(&a, &b, &ipiv, 0, IGRAPH_EINVAL);
    }

    VERIFY_FINALLY_STACK();
    return 0;
}
