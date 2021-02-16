/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard st, Cambridge MA, 02139 USA

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

int main() {

    igraph_t g;
    const igraph_vector_t v = IGRAPH_VECTOR_NULL;
    igraph_real_t edges1[] = { 0, 1, 1, 2, 2, 2, 2, 3, 2, 4, 3, 4 };
    igraph_vector_t from, to;
    igraph_es_t it;
    igraph_integer_t size;
    long int i;

    igraph_vector_view(&v, edges1, sizeof(edges1) / sizeof(igraph_real_t));

    /******************************************/
    /* Directed graph                         */
    /******************************************/

    igraph_create(&g, &v, 0, IGRAPH_DIRECTED);

    /* {0,1} -> {2,3}, result should be { 1->2 } */
    igraph_vector_init(&from, 2);
    VECTOR(from)[0] = 0;
    VECTOR(from)[1] = 1;
    igraph_vector_init(&to, 2);
    VECTOR(to)  [0] = 2;
    VECTOR(to)  [1] = 3;
    igraph_es_fromto(&g, &it, IGRAPH_VS_VECTOR(&g, &from),
                     IGRAPH_VS_VECTOR(&g, &to), IGRAPH_DIRECTED);
    igraph_vector_clear(&from);
    igraph_vector_clear(&to);
    igraph_es_size(&g, &it, &size);
    printf("%ld\n", (long)size);
    while (!igraph_es_end(&g, &it)) {
        igraph_vector_push_back(&from, igraph_es_from(&g, &it));
        igraph_vector_push_back(&to, igraph_es_to(&g, &it));
        igraph_es_next(&g, &it);
    }
    igraph_vector_sort(&from);
    igraph_vector_sort(&to);
    igraph_vector_print(&from);
    igraph_vector_print(&to);

    igraph_es_destroy(&it);

    igraph_vector_destroy(&from);
    igraph_vector_destroy(&to);

    return 0;
}
