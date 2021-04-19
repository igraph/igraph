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

int main() {

    igraph_heap_t h_max;
    igraph_heap_min_t h_min;
    igraph_integer_t i;
    igraph_real_t list[] = {-2, -9.999, 0, 6, 235, -2, -1000, -1, 4, 2000, 6, 0.5, 1, -9, 10};
    const igraph_integer_t l_size = sizeof(list) / sizeof(igraph_real_t);

    /* max heap init & destroy*/
    printf("Create empty max heap & destroy\n");
    igraph_heap_init(&h_max, 0);
    IGRAPH_ASSERT(igraph_heap_empty(&h_max));
    igraph_heap_destroy(&h_max);
    printf("Create empty max heap but allocate size for some elements\n");
    igraph_heap_init(&h_max, 10);
    IGRAPH_ASSERT(igraph_heap_empty(&h_max));
    igraph_heap_destroy(&h_max);

    /* min heap init & destroy*/
    printf("Create empty min heap & destroy\n");
    igraph_heap_min_init(&h_min, 0);
    IGRAPH_ASSERT(igraph_heap_min_empty(&h_min));
    igraph_heap_min_destroy(&h_min);
    printf("Create empty min heap but allocate size for some elements\n");
    igraph_heap_min_init(&h_min, 10);
    IGRAPH_ASSERT(igraph_heap_min_empty(&h_min));
    igraph_heap_min_destroy(&h_min);

    /*  max heap_reserve, heap_size and heap_empty*/
    printf("Test max heap_reserve, heap_size and heap_empty\n");
    igraph_heap_init(&h_max, 5);
    IGRAPH_ASSERT(igraph_heap_empty(&h_max));
    IGRAPH_ASSERT(igraph_heap_size(&h_max) == 0);
    igraph_heap_reserve(&h_max, 10);
    IGRAPH_ASSERT(igraph_heap_empty(&h_max));
    IGRAPH_ASSERT(igraph_heap_size(&h_max) == 0);
    for (i=0; i < 15; i++){
        igraph_heap_push(&h_max,i);
    }
    IGRAPH_ASSERT(igraph_heap_size(&h_max) == 15);
    IGRAPH_ASSERT(!igraph_heap_empty(&h_max));
    igraph_heap_reserve(&h_max, 5);
    IGRAPH_ASSERT(igraph_heap_size(&h_max) == 15);
    IGRAPH_ASSERT(!igraph_heap_empty(&h_max));
    igraph_heap_destroy(&h_max);

    /*  min heap reserve, heap_size and heap_empty*/
    printf("Test min heap_reserve, heap_size and heap_empty\n");
    igraph_heap_min_init(&h_min, 5);
    IGRAPH_ASSERT(igraph_heap_min_empty(&h_min));
    IGRAPH_ASSERT(igraph_heap_min_size(&h_min) == 0);
    igraph_heap_min_reserve(&h_min, 10);
    IGRAPH_ASSERT(igraph_heap_min_empty(&h_min));
    IGRAPH_ASSERT(igraph_heap_min_size(&h_min) == 0);
    for (i=0; i < 15; i++){
        igraph_heap_min_push(&h_min, i);
    }
    IGRAPH_ASSERT(igraph_heap_min_size(&h_min) == 15);
    IGRAPH_ASSERT(!igraph_heap_min_empty(&h_min));
    igraph_heap_min_reserve(&h_min, 5);
    IGRAPH_ASSERT(igraph_heap_min_size(&h_min) == 15);
    IGRAPH_ASSERT(!igraph_heap_min_empty(&h_min));
    igraph_heap_min_destroy(&h_min);

    /* max heap init_array and delete_top */
    printf("Test max heap_init array and delete_top\n");
    igraph_heap_init_array(&h_max,list,l_size);
    while (igraph_heap_size(&h_max) > 0){
        printf("%g ", igraph_heap_delete_top(&h_max));
    }
    printf("\n");
    igraph_heap_destroy(&h_max);

    /* min heap init_array and delete_top */
    printf("Test min heap init_array and delete_top\n");
    igraph_heap_min_init_array(&h_min, list, l_size);
    while (igraph_heap_min_size(&h_min) > 0) {
        printf("%g ", igraph_heap_min_delete_top(&h_min));
    }
    printf("\n");
    igraph_heap_min_destroy(&h_min);

    /* max heap top and push */
    printf("Test max heap top and push\n");
    igraph_heap_init(&h_max, 0);
    for (i=0; i < l_size; i++){
        igraph_heap_push(&h_max, list[i]);
        printf("%g ", igraph_heap_top(&h_max));
    }
    printf("\n");
    while (igraph_heap_size(&h_max)>0){
        printf("%g ", igraph_heap_delete_top(&h_max));
    }
    printf("\n");

    /* min heap top and push */
    printf("Test min heap top and push\n");
    igraph_heap_min_init(&h_min, 0);
    for (i=0; i < l_size; i++){
        igraph_heap_min_push(&h_min, list[i]);
        printf("%g ", igraph_heap_min_top(&h_min));
    }
    printf("\n");
    while (igraph_heap_min_size(&h_min) > 0){
        printf("%g ", igraph_heap_min_delete_top(&h_min));
    }
    printf("\n");

    igraph_heap_destroy(&h_max);
    igraph_heap_min_destroy(&h_min);

    VERIFY_FINALLY_STACK();

    return 0;
}
