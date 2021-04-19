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

void print_and_destroy(igraph_matrix_t *vectors, igraph_matrix_t *values, int nev, igraph_error_type_t error) {
    printf("vectors in:\n");
    print_matrix_format(vectors, stdout, "%6.2f");
    printf("values in:\n");
    print_matrix_format(values, stdout, "%6.2f");
    IGRAPH_ASSERT(igraph_arpack_unpack_complex(vectors, values, nev) == (int)error);
    printf("vectors out:\n");
    print_matrix_format(vectors, stdout, "%6.2f");
    printf("values out:\n");
    print_matrix_format(values, stdout, "%6.2f");
    igraph_matrix_destroy(vectors);
    igraph_matrix_destroy(values);
    printf("\n");
}

int main() {
    igraph_matrix_t vectors;
    igraph_matrix_t values;

    printf("Empty vectors and values:\n");
    matrix_init_int_row_major(&vectors, 0, 0, NULL);
    matrix_init_int_row_major(&values, 0, 0, NULL);
    print_and_destroy(&vectors, &values, 0, IGRAPH_SUCCESS);

    {
        printf("Real vectors and values:\n");
        int vectors_elem[4] = {-1, 0, 9, 10};
        int values_elem[4] = {-6, 0, 3, 0};
        matrix_init_int_row_major(&vectors, 2, 2, vectors_elem);
        matrix_init_int_row_major(&values, 2, 2, values_elem);
        print_and_destroy(&vectors, &values, 2, IGRAPH_SUCCESS);
    }
    {
        printf("Complex vectors and values:\n");
        igraph_real_t vectors_elem[36] = {
            0.123938, 0.3411, 0.114301, -0.134822, -0.421672, -0.484969,
            -0.268889, 0.00766665, 0.413844, 0.200565, -0.0336028, -0.133362,
            -0.192782, -0.140134, 0.579782, -0.0853149, -0.0684855, 0.117105,
            0.175547, 0.1833, 0.156218, 0.0623488, 0.422265, -0.257261,
            -0.266691, 0.404647, -0.462498, -0.0885737, 0.203893, -0.135195,
            0.662813, -0.022972, -0.193704, 0.355354, -0.0405741, 0.493652};
        igraph_real_t values_elem[12] = {
            -2.58338, 9.66092, -2.58338, -9.66092, 7.07998, 6.51033,
            7.07998, -6.51033, -7.9966, 2.74368, -7.9966, -2.74368};
        matrix_init_real_row_major(&vectors, 6, 6, vectors_elem);
        matrix_init_real_row_major(&values, 6, 2, values_elem);
        print_and_destroy(&vectors, &values, 6, IGRAPH_SUCCESS);
    }
    {
        printf("Both complex and real vectors and values:\n");
        int vectors_elem[16] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
        int values_elem[8] = {1, 0, 2, 1, 2, -1, 3, 0};
        matrix_init_int_row_major(&vectors, 4, 4, vectors_elem);
        matrix_init_int_row_major(&values, 4, 2, values_elem);
        print_and_destroy(&vectors, &values, 4, IGRAPH_SUCCESS);
    }
    {
        printf("Both complex and real vectors and values, but nev = 2:\n");
        int vectors_elem[16] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
        int values_elem[8] = {1, 0, 2, 1, 2, -1, 3, 0};
        matrix_init_int_row_major(&vectors, 4, 4, vectors_elem);
        matrix_init_int_row_major(&values, 4, 2, values_elem);
        print_and_destroy(&vectors, &values, 2, IGRAPH_SUCCESS);
    }
    {
        printf("No vectors but there are values:\n");
        int values_elem[8] = {1, 0, 2, 1, 2, -1, 3, 0};
        matrix_init_int_row_major(&vectors, 0, 0, NULL);
        matrix_init_int_row_major(&values, 4, 2, values_elem);
        print_and_destroy(&vectors, &values, 4, IGRAPH_SUCCESS);
    }
    {
        printf("No values, but there are vectors:\n");
        int vectors_elem[16] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
        matrix_init_int_row_major(&vectors, 4, 4, vectors_elem);
        matrix_init_int_row_major(&values, 0, 0, NULL);
        print_and_destroy(&vectors, &values, 0, IGRAPH_SUCCESS);
    }

    VERIFY_FINALLY_STACK();
    return 0;
}
