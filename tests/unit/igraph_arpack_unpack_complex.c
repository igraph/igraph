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
    print_matrix(vectors);
    printf("values in:\n");
    print_matrix(values);
    IGRAPH_ASSERT(igraph_arpack_unpack_complex(vectors, values, nev) == error);
    printf("vectors out:\n");
    print_matrix(vectors);
    printf("values out:\n");
    print_matrix(values);
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
        int vectors_elem[4] = {0, 3, -3, 0};
        int values_elem[4] = {2, 3, 2, -3};
        matrix_init_int_row_major(&vectors, 2, 2, vectors_elem);
        matrix_init_int_row_major(&values, 2, 2, values_elem);
        print_and_destroy(&vectors, &values, 2, IGRAPH_SUCCESS);
    }

    VERIFY_FINALLY_STACK();
    return 0;
}
