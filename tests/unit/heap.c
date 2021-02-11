/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge, MA 02139 USA

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

#include "test_utilities.inc"

int main() {

    igraph_heap_t h;
    igraph_integer_t i;
    igraph_real_t *arr =(igraph_real_t*) calloc(5, sizeof(igraph_real_t));

    /* igraph_heap_init, igraph_heap_destroy  */
    printf("Test igraph_heap_init, igraph_heap_destroy\n");
    igraph_heap_init(&h, 0);
    igraph_heap_destroy(&h);
    igraph_heap_init(&h, 10);
    igraph_heap_destroy(&h);

    /* igraph_heap_reserve */
    printf("Test igraph_heap_reserve\n");
    igraph_heap_init(&h, 0);
    igraph_heap_reserve(&h, 10);
    igraph_heap_reserve(&h, 5);
    igraph_heap_destroy(&h);

    /* igraph_heap_init_array*/
    printf("Test igraph_heap_init_array\n");
    igraph_heap_init_array(&h,arr,5);

    /* igraph_heap_empty igraph_heap_size */
    printf("Test igraph_heap_size, igraph_heap_empty\n");
    IGRAPH_ASSERT((igraph_heap_size(&h)==5) && (!igraph_heap_empty(&h)));
    igraph_heap_push(&h,100);
    IGRAPH_ASSERT((igraph_heap_size(&h)==6) && (!igraph_heap_empty(&h)));
    igraph_heap_delete_top(&h);
    IGRAPH_ASSERT((igraph_heap_size(&h)==5) && (!igraph_heap_empty(&h)));
    igraph_heap_destroy(&h);
    igraph_heap_init(&h, 0);
    IGRAPH_ASSERT((igraph_heap_size(&h)==0) && (igraph_heap_empty(&h)));

    /* igraph_heap_push igraph_heap_top  */
    printf("Test igraph_heap_push, igraph_heap_top\n");
    igraph_real_t random_list[13]={-9.999,0,6,235,-2,-1000,-1,4,2000,0.1,1,-9,10};
    igraph_vector_t top_sequence;
    igraph_vector_init(&top_sequence,13);
    for (i=0;i<13;i++)
    {
        igraph_heap_push(&h,random_list[i]);
        VECTOR(top_sequence)[i]=igraph_heap_top(&h);
    }
    print_vector_format(&top_sequence,stdout,"%g");

    /* igraph_heap_delete_top  */
    printf("Test igraph_heap_delete_top\n");
    for (i=0;i<13;i++)  VECTOR(top_sequence)[i]=(&top_sequence,igraph_heap_delete_top(&h));
    print_vector_format(&top_sequence,stdout,"%g");
    IGRAPH_ASSERT(igraph_heap_empty(&h));

    free(arr);
    igraph_heap_destroy(&h);
    igraph_vector_destroy(&top_sequence);

    VERIFY_FINALLY_STACK();

    return 0;
}
