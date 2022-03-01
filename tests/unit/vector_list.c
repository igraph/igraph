/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2022 The igraph development team

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

int main() {
    igraph_vector_int_list_t list;
    igraph_vector_int_t v;
    igraph_integer_t i;

    printf("Initialise empty vector list\n");
    igraph_vector_int_list_init(&list, 0);
    print_vector_int_list(&list);
    igraph_vector_int_list_destroy(&list);

    printf("Initialise vector list of length 10\n");
    igraph_vector_int_list_init(&list, 10);
    print_vector_int_list(&list);
    igraph_vector_int_list_destroy(&list);

    printf("Test igraph_vector_int_list_get_ptr and igraph_vector_int_list_size\n");
    igraph_vector_int_list_init(&list, 10);
    for (i = 0; i < igraph_vector_int_list_size(&list); i++) {
        igraph_vector_int_push_back(igraph_vector_int_list_get_ptr(&list, i), 10 - i);
    }
    print_vector_int_list(&list);
    igraph_vector_int_list_destroy(&list);

    printf("Test igraph_vector_int_list_reserve and igraph_vector_int_list_push_back_new\n");
    igraph_vector_int_list_init(&list, 0);
    igraph_vector_int_list_reserve(&list, 10);
    for (i = 0; i < 10; i++) {
        igraph_vector_int_list_push_back_new(&list, /* item = */ 0);
    }

    printf("Test igraph_vector_int_list_empty and igraph_vector_int_list_clear\n");
    IGRAPH_ASSERT(!igraph_vector_int_list_empty(&list));
    igraph_vector_int_list_clear(&list);
    IGRAPH_ASSERT(igraph_vector_int_list_empty(&list));
    igraph_vector_int_list_destroy(&list);

    printf("Test igraph_vector_list_set\n");
    igraph_vector_int_list_init(&list, 5);
    for (i = 0; i < igraph_vector_int_list_size(&list); i++) {
        igraph_vector_int_init(&v, 1);
        VECTOR(v)[0] = 20 * i;
        igraph_vector_int_list_set(&list, i, &v);  /* ownership taken */
    }
    print_vector_int_list(&list);
    igraph_vector_int_list_destroy(&list);

    printf("Test igraph_vector_int_list_tail_ptr, igraph_vector_int_list_pop_back\n");
    igraph_vector_int_list_init(&list, 10);
    for (i = 0; i < igraph_vector_int_list_size(&list); i++) {
        igraph_vector_int_push_back(igraph_vector_int_list_get_ptr(&list, i), i + 1);
    }
    while (!igraph_vector_int_list_empty(&list)) {
        print_vector_int(igraph_vector_int_list_tail_ptr(&list));
        IGRAPH_ASSERT(igraph_vector_int_list_pop_back(&list, &v) == IGRAPH_SUCCESS);
        /* v is now owned by us, not the vector_int_list */
        print_vector_int(&v);
        igraph_vector_int_destroy(&v);
    }
    printf("\n");
    igraph_vector_int_list_destroy(&list);

    printf("Test igraph_vector_int_list_resize, igraph_vector_int_list_sort\n");
    igraph_vector_int_list_init(&list, 20);
    for (i = 0; i < 10; i++) {
        igraph_vector_int_push_back(igraph_vector_int_list_get_ptr(&list, i), 10 - i);
    }
    igraph_vector_int_list_resize(&list, 10);
    igraph_vector_int_list_sort(&list, igraph_vector_int_lex_cmp);
    print_vector_int_list(&list);
    igraph_vector_int_list_destroy(&list);

    printf("Test igraph_vector_int_list_remove\n");
    igraph_vector_int_list_init(&list, 10);
    for (i = 0; i < 10; i++) {
        igraph_vector_int_push_back(igraph_vector_int_list_get_ptr(&list, i), i + 1);
    }
    igraph_vector_int_list_remove(&list, 9, /* item = */ &v);
    print_vector_int(&v);  /* v is now owned by us */
    igraph_vector_int_destroy(&v);
    igraph_vector_int_list_remove(&list, 0, /* item = */ &v);
    print_vector_int(&v);  /* v is now owned by us */
    igraph_vector_int_destroy(&v);
    igraph_vector_int_list_remove(&list, 4, /* item = */ &v);
    print_vector_int(&v);  /* v is now owned by us */
    igraph_vector_int_destroy(&v);
    print_vector_int_list(&list);
    igraph_vector_int_list_destroy(&list);

    printf("Test igraph_vector_int_list_remove_fast\n");
    igraph_vector_int_list_init(&list, 10);
    for (i = 0; i < 10; i++) {
        igraph_vector_int_push_back(igraph_vector_int_list_get_ptr(&list, i), i + 1);
    }
    igraph_vector_int_list_remove_fast(&list, 9, /* item = */ &v);
    print_vector_int(&v);  /* v is now owned by us */
    igraph_vector_int_destroy(&v);
    igraph_vector_int_list_remove_fast(&list, 0, /* item = */ &v);
    print_vector_int(&v);  /* v is now owned by us */
    igraph_vector_int_destroy(&v);
    igraph_vector_int_list_remove_fast(&list, 4, /* item = */ &v);
    print_vector_int(&v);  /* v is now owned by us */
    igraph_vector_int_destroy(&v);
    print_vector_int_list(&list);
    igraph_vector_int_list_destroy(&list);

    printf("Test igraph_vector_int_list_discard\n");
    igraph_vector_int_list_init(&list, 10);
    for (i = 0; i < 10; i++) {
        igraph_vector_int_push_back(igraph_vector_int_list_get_ptr(&list, i), i + 1);
    }
    igraph_vector_int_list_discard(&list, 9);
    igraph_vector_int_list_discard(&list, 0);
    igraph_vector_int_list_discard(&list, 4);
    print_vector_int_list(&list);
    igraph_vector_int_list_destroy(&list);

    printf("Test igraph_vector_int_list_discard_fast\n");
    igraph_vector_int_list_init(&list, 10);
    for (i = 0; i < 10; i++) {
        igraph_vector_int_push_back(igraph_vector_int_list_get_ptr(&list, i), i + 1);
    }
    igraph_vector_int_list_discard_fast(&list, 9);
    igraph_vector_int_list_discard_fast(&list, 0);
    igraph_vector_int_list_discard_fast(&list, 4);
    print_vector_int_list(&list);
    igraph_vector_int_list_destroy(&list);

    VERIFY_FINALLY_STACK();

    return 0;
}
