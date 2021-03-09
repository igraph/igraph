/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge MA, 02139 USA

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

igraph_vector_ptr_t custom_destructor_stack;

void custom_destructor(void* ptr) {
    igraph_vector_ptr_push_back(&custom_destructor_stack, ptr);
}

int main() {

    igraph_vector_ptr_t v1, v2;
    igraph_vector_ptr_t v3 = IGRAPH_VECTOR_PTR_NULL;
    int i;
    void ** ptr;
    int d1 = 1, d2 = 2, d3 = 3, d4 = 4, d5 = 5;
    char *block1 = 0, *block2 = 0;

    /* igraph_vector_ptr_init, igraph_vector_ptr_destroy */
    igraph_vector_ptr_init(&v1, 10);
    igraph_vector_ptr_destroy(&v1);
    igraph_vector_ptr_init(&v1, 0);
    igraph_vector_ptr_destroy(&v1);

    /* igraph_vector_ptr_free_all, igraph_vector_ptr_destroy_all */
    igraph_vector_ptr_init(&v1, 5);
    for (i = 0; i < igraph_vector_ptr_size(&v1); i++) {
        VECTOR(v1)[i] = (void*)malloc(i * 10);
    }
    igraph_vector_ptr_free_all(&v1);
    for (i = 0; i < igraph_vector_ptr_size(&v1); i++) {
        VECTOR(v1)[i] = (void*)malloc(i * 10);
    }
    igraph_vector_ptr_destroy_all(&v1);

    /* igraph_vector_ptr_reserve */
    igraph_vector_ptr_init(&v1, 0);
    igraph_vector_ptr_reserve(&v1, 5);
    igraph_vector_ptr_reserve(&v1, 15);
    igraph_vector_ptr_reserve(&v1, 1);
    igraph_vector_ptr_reserve(&v1, 0);
    igraph_vector_ptr_destroy(&v1);

    /* igraph_vector_ptr_empty, igraph_vector_ptr_clear */
    igraph_vector_ptr_init(&v1, 10);
    if (igraph_vector_ptr_empty(&v1)) {
        return 1;
    }
    igraph_vector_ptr_clear(&v1);
    if (!igraph_vector_ptr_empty(&v1)) {
        return 2;
    }

    /* igraph_vector_ptr_size */
    if (igraph_vector_ptr_size(&v1) != 0) {
        return 3;
    }
    igraph_vector_ptr_resize(&v1, 10);
    if (igraph_vector_ptr_size(&v1) != 10) {
        return 4;
    }
    igraph_vector_ptr_destroy(&v1);

    /* igraph_vector_ptr_push_back */
    igraph_vector_ptr_init(&v1, 0);
    for (i = 0; i < 10; i++) {
        igraph_vector_ptr_push_back(&v1, (void*)malloc(i * 10));
    }
    igraph_vector_ptr_destroy_all(&v1);

    /* igraph_vector_ptr_e */
    igraph_vector_ptr_init(&v1, 5);
    VECTOR(v1)[0] = &d1;
    VECTOR(v1)[1] = &d2;
    VECTOR(v1)[2] = &d3;
    VECTOR(v1)[3] = &d4;
    VECTOR(v1)[4] = &d5;
    if (igraph_vector_ptr_e(&v1, 0) != &d1) {
        return 5;
    }
    if (igraph_vector_ptr_e(&v1, 1) != &d2) {
        return 6;
    }
    if (igraph_vector_ptr_e(&v1, 2) != &d3) {
        return 7;
    }
    if (igraph_vector_ptr_e(&v1, 3) != &d4) {
        return 8;
    }
    if (igraph_vector_ptr_e(&v1, 4) != &d5) {
        return 9;
    }
    igraph_vector_ptr_destroy(&v1);

    /* igraph_vector_ptr_set */
    igraph_vector_ptr_init(&v1, 5);
    igraph_vector_ptr_set(&v1, 0, &d1);
    igraph_vector_ptr_set(&v1, 1, &d2);
    igraph_vector_ptr_set(&v1, 2, &d3);
    igraph_vector_ptr_set(&v1, 3, &d4);
    igraph_vector_ptr_set(&v1, 4, &d5);
    if (igraph_vector_ptr_e(&v1, 0) != &d1) {
        return 5;
    }
    if (igraph_vector_ptr_e(&v1, 1) != &d2) {
        return 6;
    }
    if (igraph_vector_ptr_e(&v1, 2) != &d3) {
        return 7;
    }
    if (igraph_vector_ptr_e(&v1, 3) != &d4) {
        return 8;
    }
    if (igraph_vector_ptr_e(&v1, 4) != &d5) {
        return 9;
    }
    igraph_vector_ptr_destroy(&v1);

    /* igraph_vector_ptr_null */
    igraph_vector_ptr_init(&v1, 5);
    igraph_vector_ptr_set(&v1, 0, &d1);
    igraph_vector_ptr_set(&v1, 1, &d2);
    igraph_vector_ptr_set(&v1, 2, &d3);
    igraph_vector_ptr_set(&v1, 3, &d4);
    igraph_vector_ptr_set(&v1, 4, &d5);
    igraph_vector_ptr_null(&v1);
    for (i = 0; i < igraph_vector_ptr_size(&v1); i++) {
        if (VECTOR(v1)[i] != 0) {
            return 10;
        }
    }
    igraph_vector_ptr_destroy(&v1);

    /* igraph_vector_ptr_resize */
    igraph_vector_ptr_init(&v1, 10);
    igraph_vector_ptr_set(&v1, 0, &d1);
    igraph_vector_ptr_set(&v1, 1, &d2);
    igraph_vector_ptr_set(&v1, 2, &d3);
    igraph_vector_ptr_set(&v1, 3, &d4);
    igraph_vector_ptr_set(&v1, 4, &d5);
    igraph_vector_ptr_resize(&v1, 10);
    igraph_vector_ptr_resize(&v1, 15);
    igraph_vector_ptr_resize(&v1, 5);
    if (igraph_vector_ptr_size(&v1) != 5) {
        return 11;
    }
    if (igraph_vector_ptr_e(&v1, 0) != &d1) {
        return 12;
    }
    if (igraph_vector_ptr_e(&v1, 1) != &d2) {
        return 13;
    }
    if (igraph_vector_ptr_e(&v1, 2) != &d3) {
        return 14;
    }
    if (igraph_vector_ptr_e(&v1, 3) != &d4) {
        return 15;
    }
    if (igraph_vector_ptr_e(&v1, 4) != &d5) {
        return 16;
    }
    igraph_vector_ptr_destroy(&v1);

    /* igraph_vector_ptr_view */
    ptr = (void**) malloc(5 * sizeof(void*));
    igraph_vector_ptr_view(&v3, ptr, 5);
    ptr[0] = &d1;
    ptr[1] = &d2;
    ptr[2] = &d3;
    ptr[3] = &d4;
    ptr[4] = &d5;
    for (i = 0; i < igraph_vector_ptr_size(&v3); i++) {
        if ( *((int*)VECTOR(v3)[i]) != i + 1) {
            return 17;
        }
    }

    /* igraph_vector_ptr_init_copy */
    igraph_vector_ptr_init_copy(&v1, ptr, 5);
    for (i = 0; i < igraph_vector_ptr_size(&v1); i++) {
        if ( *((int*)VECTOR(v1)[i]) != i + 1) {
            return 18;
        }
    }

    /* igraph_vector_ptr_copy_to */
    igraph_vector_ptr_copy_to(&v1, ptr);
    for (i = 0; i < igraph_vector_ptr_size(&v1); i++) {
        if ( *((int*)ptr[i]) != i + 1) {
            return 19;
        }
    }
    free(ptr);
    igraph_vector_ptr_destroy(&v1);

    /* igraph_vector_ptr_copy */
    igraph_vector_ptr_init(&v1, 5);
    igraph_vector_ptr_set(&v1, 0, &d1);
    igraph_vector_ptr_set(&v1, 1, &d2);
    igraph_vector_ptr_set(&v1, 2, &d3);
    igraph_vector_ptr_set(&v1, 3, &d4);
    igraph_vector_ptr_set(&v1, 4, &d5);
    igraph_vector_ptr_copy(&v2, &v1);
    igraph_vector_ptr_destroy(&v1);
    for (i = 0; i < igraph_vector_ptr_size(&v2); i++) {
        if ( *((int*)VECTOR(v2)[i]) != i + 1) {
            return 20;
        }
    }

    /* igraph_vector_ptr_remove */
    igraph_vector_ptr_remove(&v2, 0);
    igraph_vector_ptr_remove(&v2, 3);
    if ( *((int*)VECTOR(v2)[0]) != 2) {
        return 21;
    }
    if ( *((int*)VECTOR(v2)[1]) != 3) {
        return 22;
    }
    if ( *((int*)VECTOR(v2)[2]) != 4) {
        return 23;
    }

    igraph_vector_ptr_destroy(&v2);

    /* Testing destructor */
    igraph_vector_ptr_init(&custom_destructor_stack, 0);
    igraph_vector_ptr_init(&v1, 2);
    block1 = IGRAPH_CALLOC(32, char);
    block2 = IGRAPH_CALLOC(64, char);
    VECTOR(v1)[0] = block1;
    VECTOR(v1)[1] = block2;
    if (igraph_vector_ptr_get_item_destructor(&v1) != 0) {
        return 24;
    }
    if (igraph_vector_ptr_set_item_destructor(&v1, &custom_destructor) != 0) {
        return 25;
    }
    /* Okay, let's clear the vector. This should push the blocks in the
     * custom destructor stack */
    igraph_vector_ptr_clear(&v1);
    /* Put the blocks back and destroy the vector */
    igraph_vector_ptr_push_back(&v1, block1);
    igraph_vector_ptr_push_back(&v1, block2);
    igraph_vector_ptr_destroy_all(&v1);

    if (VECTOR(custom_destructor_stack)[0] != block1 ||
        VECTOR(custom_destructor_stack)[1] != block2 ||
        VECTOR(custom_destructor_stack)[2] != block1 ||
        VECTOR(custom_destructor_stack)[3] != block2
       ) {
        return 26;
    }

    igraph_vector_ptr_destroy(&custom_destructor_stack);

    VERIFY_FINALLY_STACK();

    return 0;
}
