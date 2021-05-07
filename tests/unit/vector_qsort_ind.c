/* -*- mode: C -*-  */
/*
   IGraph library.
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

#include "test_utilities.inc"

static int compare_first_items(const void* a, const void* b) {
    igraph_vector_t *vec1 = (igraph_vector_t*) a;
    igraph_vector_t *vec2 = (igraph_vector_t*) b;

    return VECTOR(*vec1)[0] - VECTOR(*vec2)[0];
}

int main() {
    igraph_vector_t vector;
    igraph_vector_t indices;
    igraph_real_t values[] = { 87, 23, 8, 82, 94, 56, 36, 33, 76, 66 };
    igraph_real_t values2[] = { 87, 23, 8, 82, 94, 56, 36, 33, 76, 66 };
    int i, *ptr;

    /* Special case: empty vector */
    igraph_vector_init(&vector, 0);
    igraph_vector_init(&indices, 0);
    igraph_vector_qsort_ind(&vector, &indices, /* descending = */ 0);
    print_vector(&indices);
    igraph_vector_destroy(&indices);
    igraph_vector_destroy(&vector);

    /* Non-empty vector, descending */
    igraph_vector_view(&vector, values, sizeof(values) / sizeof(igraph_real_t));
    igraph_vector_init(&indices, 0);
    igraph_vector_qsort_ind(&vector, &indices, /* descending = */ 1);
    print_vector(&indices);

    /* Non-empty vector, ascending */
    igraph_vector_view(&vector, values2, sizeof(values2) / sizeof(igraph_real_t));
    igraph_vector_qsort_ind(&vector, &indices, /* descending = */ 0);
    print_vector(&indices);

    /* Permute the vector by the index vector */
    igraph_vector_permute(&vector, &indices);
    print_vector(&vector);

    /* Print and clean up */
    igraph_vector_destroy(&indices);

    /* Check finalizer stack */
    VERIFY_FINALLY_STACK();

    return 0;
}
