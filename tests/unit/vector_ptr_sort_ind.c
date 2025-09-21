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

static int compare_first_items(const void* a, const void* b) {
    igraph_vector_t *vec1 = (igraph_vector_t*) a;
    igraph_vector_t *vec2 = (igraph_vector_t*) b;

    return VECTOR(*vec1)[0] - VECTOR(*vec2)[0];
}

int main(void) {
    igraph_vector_ptr_t vectors;
    igraph_vector_t* vec;
    igraph_vector_int_t indices;
    int values[] = { 3, 5, 0, 7, 9, 6, 2, 8, 1, 4, -1 };
    int i, *ptr;

    /* Special case: empty vector */
    igraph_vector_ptr_init(&vectors, 0);
    igraph_vector_int_init(&indices, 0);
    igraph_vector_ptr_sort_ind(&vectors, &indices, compare_first_items);
    print_vector_int(&indices);
    igraph_vector_int_destroy(&indices);
    igraph_vector_ptr_destroy_all(&vectors);

    /* Create vectors of length 1, each containing a value from 'values', and
     * put them in a vector of pointers */
    igraph_vector_ptr_init(&vectors, 0);
    for (ptr = values; *ptr >= 0; ptr++) {
        vec = IGRAPH_CALLOC(1, igraph_vector_t);
        igraph_vector_init(vec, 1);
        VECTOR(*vec)[0] = *ptr;
        igraph_vector_ptr_push_back(&vectors, vec);
    }

    /* Sort the vector of vectors by the first item of each vector, and get
     * the index vector */
    igraph_vector_int_init(&indices, 0);
    igraph_vector_ptr_sort_ind(&vectors, &indices, compare_first_items);
    print_vector_int(&indices);

    /* Permute the vector of vectors by the index vector */
    igraph_vector_ptr_permute(&vectors, &indices);

    /* Print and clean up */
    for (i = 0; i < igraph_vector_ptr_size(&vectors); i++) {
        print_vector(VECTOR(vectors)[i]);
        igraph_vector_destroy(VECTOR(vectors)[i]);
    }
    igraph_vector_ptr_destroy_all(&vectors);
    igraph_vector_int_destroy(&indices);

    /* Check finalizer stack */
    VERIFY_FINALLY_STACK();

    return 0;
}
