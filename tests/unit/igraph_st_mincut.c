/*
   igraph library.
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

void sort_and_print(igraph_vector_int_t *vec) {
    igraph_vector_int_sort(vec);
    print_vector_int(vec);
}

int main(void) {

    igraph_t g;
    igraph_vector_int_t cut, partition, partition2;
    igraph_vector_t capacity;
    igraph_real_t value;
    int source = 0;
    int target = 4;

    igraph_vector_int_init(&partition, 0);
    igraph_vector_int_init(&partition2, 0);
    igraph_vector_int_init(&cut, 0);

    igraph_small(&g, 5, IGRAPH_DIRECTED, 0, 1, 1, 2, 1, 3, 2, 4, 3, 4, -1);
    igraph_vector_init_int_end(&capacity, -1, 8, 2, 3, 3, 2, -1);

    /*    test without capacity    */

    igraph_st_mincut(&g, &value, &cut, &partition, &partition2, source, target, /*capacity*/ NULL);

    /* cut and partition should have only one element */
    print_vector_int(&cut);
    print_vector_int(&partition);
    sort_and_print(&partition2);

    /*    test with capacity    */

    igraph_st_mincut(&g, &value, &cut, &partition, &partition2, source, target, &capacity);

    sort_and_print(&cut);
    sort_and_print(&partition);
    sort_and_print(&partition2);

    /*     cleanup     */

    igraph_vector_int_destroy(&cut);
    igraph_vector_int_destroy(&partition);
    igraph_vector_int_destroy(&partition2);
    igraph_vector_destroy(&capacity);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();
    return 0;
}
