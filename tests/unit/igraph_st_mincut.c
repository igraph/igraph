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

void sort_and_print(igraph_vector_t *vec) {
    igraph_vector_sort(vec);
    print_vector_round(vec);
}

int main() {

    igraph_t g;
    igraph_vector_t cut, partition, partition2, capacity;
    igraph_real_t value;
    int source = 0;
    int target = 4;
    int cut_edge;
    int partition_vertex;

    igraph_vector_init(&partition, 0);
    igraph_vector_init(&partition2, 0);
    igraph_vector_init(&cut, 0);

    igraph_small(&g, 5, IGRAPH_DIRECTED, 0, 1, 1, 2, 1, 3, 2, 4, 3, 4, -1);
    igraph_vector_init_int_end(&capacity, -1, 8, 2, 3, 3, 2, -1);

    /*    test without capacity    */

    igraph_st_mincut(&g, &value, &cut, &partition, &partition2, source, target, /*capacity*/ NULL);

    /* cut and partition should have only one element */
    print_vector_round(&cut);
    print_vector_round(&partition);
    sort_and_print(&partition2);

    /*    test with capacity    */

    igraph_st_mincut(&g, &value, &cut, &partition, &partition2, source, target, &capacity);

    sort_and_print(&cut);
    sort_and_print(&partition);
    sort_and_print(&partition2);

    /*     cleanup     */

    igraph_vector_destroy(&cut);
    igraph_vector_destroy(&partition);
    igraph_vector_destroy(&partition2);
    igraph_vector_destroy(&capacity);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();
    return 0;
}
