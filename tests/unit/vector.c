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
#include <assert.h>

int main() {

    igraph_vector_t v, v2, v3;
    int i;
    igraph_real_t *ptr;
    long int pos;

    printf("Initialise empty vector\n");
    igraph_vector_init(&v, 0);
    igraph_vector_destroy(&v);

    printf("Initialise vector of length 10\n");
    igraph_vector_init(&v, 10);
    print_vector_format(&v, stdout, "%g");
    igraph_vector_destroy(&v);

    printf("Test VECTOR() and igraph_vector_size\n");
    igraph_vector_init(&v, 10);
    for (i = 0; i < igraph_vector_size(&v); i++) {
        VECTOR(v)[i] = 10 - i;
    }
    print_vector_format(&v, stdout, "%g");
    igraph_vector_destroy(&v);

    printf("Test igraph_vector_reserve and igraph_vector_push_back\n");
    igraph_vector_init(&v, 0);
    igraph_vector_reserve(&v, 10);
    for (i = 0; i < 10; i++) {
        igraph_vector_push_back(&v, i);
    }

    printf("Test igraph_vector_empty and igraph_vector_clear\n");
    assert(!igraph_vector_empty(&v));
    igraph_vector_clear(&v);
    assert(igraph_vector_empty(&v));
    igraph_vector_destroy(&v);

    printf("Test igraph_vector_e and igraph_vector_e_ptr\n");
    igraph_vector_init(&v, 5);
    for (i = 0; i < igraph_vector_size(&v); i++) {
        *igraph_vector_e_ptr(&v, i) = 100 * i;
    }
    for (i = 0; i < igraph_vector_size(&v); i++) {
        printf(" %li", (long int)igraph_vector_e(&v, i));
    }
    printf("\n");
    igraph_vector_destroy(&v);

    printf("Test igraph_vector_set\n");
    igraph_vector_init(&v, 5);
    for (i = 0; i < igraph_vector_size(&v); i++) {
        igraph_vector_set(&v, i, 20 * i);
    }
    print_vector_format(&v, stdout, "%g");
    igraph_vector_destroy(&v);

    printf("Test igraph_vector_null\n");
    igraph_vector_init(&v, 0);
    igraph_vector_null(&v);
    igraph_vector_destroy(&v);
    igraph_vector_init(&v, 10);
    for (i = 0; i < igraph_vector_size(&v); i++) {
        VECTOR(v)[i] = i + 1;
    }
    igraph_vector_null(&v);
    print_vector_format(&v, stdout, "%g");
    igraph_vector_destroy(&v);

    printf("Test igraph_vector_tail, igraph_vector_pop_back\n");
    igraph_vector_init(&v, 10);
    for (i = 0; i < igraph_vector_size(&v); i++) {
        VECTOR(v)[i] = i + 1;
    }
    while (!igraph_vector_empty(&v)) {
        printf(" %li", (long int)igraph_vector_tail(&v));
        printf(" %li", (long int)igraph_vector_pop_back(&v));
    }
    printf("\n");
    igraph_vector_destroy(&v);

    printf("Test igraph_vector_init_seq, igraph_vector_order\n");
    igraph_vector_init_seq(&v, 1, 10);
    igraph_vector_init(&v2, 0);
    igraph_vector_order1(&v, &v2, 10);
    print_vector_format(&v2, stdout, "%g");
    igraph_vector_destroy(&v2);
    igraph_vector_destroy(&v);

    printf("Test igraph_vector_resize, igraph_vector_sort\n");
    igraph_vector_init(&v, 20);
    for (i = 0; i < 10; i++) {
        VECTOR(v)[i] = 10 - i;
    }
    igraph_vector_resize(&v, 10);
    igraph_vector_sort(&v);
    print_vector_format(&v, stdout, "%g");
    igraph_vector_destroy(&v);

    printf("Test igraph_vector_{which}_{min, max}\n");
    igraph_vector_init(&v, 10);
    for (i = 0; i < igraph_vector_size(&v); i++) {
        VECTOR(v)[i] = 100 - i;
    }
    for (i = 0; i < 10; i++) {
        printf(" %li", (long int)VECTOR(v)[i]);
    }
    printf("\n");
    assert(igraph_vector_max(&v) == 100);
    assert(igraph_vector_which_max(&v) == 0);
    assert(igraph_vector_min(&v) == 91);
    assert(igraph_vector_which_min(&v) == 9);

    printf("Test NaN values\n");
    igraph_vector_push_back(&v, IGRAPH_NAN);
    igraph_vector_push_back(&v, IGRAPH_NAN);
    igraph_vector_push_back(&v, 1);
    assert(igraph_is_nan(igraph_vector_max(&v)));
    /* Index should be to first NaN value */
    assert(igraph_vector_which_max(&v) == 10);
    assert(igraph_is_nan(igraph_vector_min(&v)));
    /* Index should be to first NaN value */
    assert(igraph_vector_which_min(&v) == 10);

    printf("Test igraph_vector_init_copy\n");
    igraph_vector_destroy(&v);
    ptr = (igraph_real_t*) malloc(10 * sizeof(igraph_real_t));
    igraph_vector_init_copy(&v, ptr, 10);
    free(ptr);
    for (i = 0; i < 10; i++) {
        VECTOR(v)[i] = 100 - i;
    }
    print_vector_format(&v, stdout, "%g");
    igraph_vector_destroy(&v);

    printf("Test igraph_vector_copy_to\n");
    ptr = (igraph_real_t*) malloc(10 * sizeof(igraph_real_t));
    igraph_vector_init_seq(&v, 11, 20);
    igraph_vector_copy_to(&v, ptr);
    for (i = 0; i < 10; i++) {
        printf(" %li", (long int)ptr[i]);
    }
    printf("\n");
    free(ptr);
    igraph_vector_destroy(&v);

    printf("Test igraph_vector_init_seq, igraph_vector_sum, igraph_vector_prod\n");
    igraph_vector_init_seq(&v, 1, 5);
    printf(" %li", (long int)igraph_vector_sum(&v));
    printf(" %li\n", (long int)igraph_vector_prod(&v));

    printf("Test igraph_vector_remove_section\n");
    igraph_vector_remove_section(&v, 2, 4);
    printf(" %li", (long int)igraph_vector_sum(&v));
    printf(" %li\n", (long int)igraph_vector_prod(&v));
    igraph_vector_destroy(&v);

    printf("Test igraph_vector_remove\n");
    igraph_vector_init_seq(&v, 1, 10);
    igraph_vector_remove(&v, 9);
    igraph_vector_remove(&v, 0);
    igraph_vector_remove(&v, 4);
    printf(" %li\n", (long int)igraph_vector_sum(&v));
    igraph_vector_destroy(&v);

    printf("Test igraph_vector_move_interval\n");
    igraph_vector_init_seq(&v, 0, 9);
    igraph_vector_move_interval(&v, 5, 10, 0);
    assert(igraph_vector_sum(&v) == 70);
    igraph_vector_destroy(&v);

    printf("Test igraph_vector_isininterval\n");
    igraph_vector_init_seq(&v, 1, 10);
    assert(igraph_vector_isininterval(&v, 1, 10));
    assert(!igraph_vector_isininterval(&v, 2, 10));
    assert(!igraph_vector_isininterval(&v, 1, 9));

    printf("Test igraph_vector_any_smaller\n");
    assert(!igraph_vector_any_smaller(&v, 1));
    assert(igraph_vector_any_smaller(&v, 2));
    igraph_vector_destroy(&v);

    printf("Test igraph_vector_all_e\n");

    printf("Test igraph_vector_binsearch\n");
    igraph_vector_init_seq(&v, 0, 9);
    for (i = 0; i < igraph_vector_size(&v); i++) {
        assert(igraph_vector_binsearch(&v, 0, 0));
    }
    assert(!igraph_vector_binsearch(&v, 10, 0));
    assert(!igraph_vector_binsearch(&v, -1, 0));
    for (i = 0; i < igraph_vector_size(&v); i++) {
        VECTOR(v)[i] = 2 * i;
    }
    for (i = 0; i < igraph_vector_size(&v); i++) {
        long int pos;
        assert(igraph_vector_binsearch(&v, VECTOR(v)[i], &pos));
        assert(pos == i);
        assert(!igraph_vector_binsearch(&v, VECTOR(v)[i] + 1, &pos));
    }
    igraph_vector_destroy(&v);

    printf("Test Binsearch in empty vector\n");
    igraph_vector_init(&v, 0);
    assert(!igraph_vector_binsearch2(&v, 0));
    assert(!igraph_vector_binsearch(&v, 1, &pos));
    assert(pos == 0);
    igraph_vector_destroy(&v);

    printf("Test igraph_vector_init_real\n");
    igraph_vector_init_real(&v, 10, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0);
    print_vector_format(&v, stdout, "%g");
    igraph_vector_destroy(&v);

    printf("Test igraph_vector_init_int\n");
    igraph_vector_init_int(&v, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10);
    print_vector_format(&v, stdout, "%g");
    igraph_vector_destroy(&v);

    printf("Test igraph_vector_init_real\n");
    igraph_vector_init_real_end(&v, -1, 1.0, 2.0, 3.0, 4.0, 5.0,
                                6.0, 7.0, 8.0, 9.0, 10.0, -1.0);
    print_vector_format(&v, stdout, "%g");
    igraph_vector_destroy(&v);

    printf("Test igraph_vector_init_int\n");
    igraph_vector_init_int_end(&v, -1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, -1);
    print_vector_format(&v, stdout, "%g");
    igraph_vector_destroy(&v);

    printf("Test igraph_vector_permdelete\n");
    printf("Test igraph_vector_remove_negidx\n");

    printf("Test order2\n");
    igraph_vector_init_int_end(&v, -1, 10, 9, 8, 7, 6, 7, 8, 9, 10, -1);
    igraph_vector_order2(&v);
    print_vector_format(&v, stdout, "%g");
    igraph_vector_destroy(&v);

    printf("Test filter_smaller, quite special....\n");
    igraph_vector_init_int_end(&v, -1, 0, 1, 2, 3, 4, 4, 4, 4, 5, 6, 7, 8, -1);
    igraph_vector_filter_smaller(&v, 4);
    print_vector_format(&v, stdout, "%g");
    igraph_vector_destroy(&v);
    igraph_vector_init_int_end(&v, -1, 1, 2, 3, 4, 4, 4, 4, 5, 6, 7, 8, -1);
    igraph_vector_filter_smaller(&v, 0);
    print_vector_format(&v, stdout, "%g");
    igraph_vector_destroy(&v);
    igraph_vector_init_int_end(&v, -1, 0, 0, 1, 2, 3, 4, 4, 4, 4, 5, 6, 7, 8, -1);
    igraph_vector_filter_smaller(&v, 0);
    print_vector_format(&v, stdout, "%g");
    igraph_vector_destroy(&v);

    printf("Test rank\n");
    igraph_vector_init_int_end(&v, -1, 0, 1, 2, 6, 5, 2, 1, 0, -1);
    igraph_vector_init(&v2, 0);
    igraph_vector_rank(&v, &v2, 7);
    print_vector_format(&v, stdout, "%g");
    print_vector_format(&v2, stdout, "%g");
    igraph_vector_destroy(&v);
    igraph_vector_destroy(&v2);

    printf("Test order\n");
    igraph_vector_init_int_end(&v,  -1, 1, 1, 2, 2, -1);
    igraph_vector_init_int_end(&v2, -1, 2, 3, 1, 3, -1);
    igraph_vector_init(&v3, 0);
    igraph_vector_order(&v, &v2, &v3, 3);
    print_vector_format(&v3, stdout, "%g");
    igraph_vector_destroy(&v);
    igraph_vector_destroy(&v2);
    igraph_vector_destroy(&v3);

    printf("Test fill\n");

    igraph_vector_init(&v, 100);
    igraph_vector_fill(&v, 1.234567);
    for (i = 0; i < igraph_vector_size(&v); i++) {
        assert(VECTOR(v)[i] == 1.234567);
    }
    igraph_vector_destroy(&v);

    VERIFY_FINALLY_STACK();

    return 0;
}

