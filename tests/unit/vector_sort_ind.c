/*
   igraph library.
   Copyright (C) 2006-2021  The igraph development team

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
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include <igraph.h>
#include <stdlib.h>

#include "test_utilities.h"

int main(void) {
    igraph_vector_t vector;
    igraph_vector_int_t indices;
    igraph_real_t values[] = { 87, 23, 8, 82, 94, 56, 36, 33, 76, 66 };
    igraph_real_t values2[] = { 87, 23, 8, 82, 94, 56, 36, 33, 76, 66 };

    /* Special case: empty vector */
    igraph_vector_init(&vector, 0);
    igraph_vector_int_init(&indices, 0);
    igraph_vector_sort_ind(&vector, &indices, IGRAPH_ASCENDING);
    print_vector_int(&indices);
    igraph_vector_int_destroy(&indices);
    igraph_vector_destroy(&vector);

    /* Non-empty vector, descending */
    vector = igraph_vector_view(values, sizeof(values) / sizeof(values[0]));
    igraph_vector_int_init(&indices, 0);
    igraph_vector_sort_ind(&vector, &indices, IGRAPH_DESCENDING);
    print_vector_int(&indices);

    /* Non-empty vector, ascending */
    vector = igraph_vector_view(values2, sizeof(values2) / sizeof(values2[0]));
    igraph_vector_sort_ind(&vector, &indices, IGRAPH_ASCENDING);
    print_vector_int(&indices);

    /* Permute the vector by the index vector */
    igraph_vector_permute(&vector, &indices);
    print_vector(&vector);

    /* Print and clean up */
    igraph_vector_int_destroy(&indices);

    /* Check finalizer stack */
    VERIFY_FINALLY_STACK();

    return 0;
}
