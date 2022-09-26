/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2010-2012  Gabor Csardi <csardi.gabor@gmail.com>
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
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include <igraph.h>
#include <stdio.h>

int main(void) {
    igraph_t graph;
    igraph_vector_int_list_t separators;
    igraph_integer_t i, n;

    igraph_famous(&graph, "zachary");
    igraph_vector_int_list_init(&separators, 0);
    igraph_all_minimal_st_separators(&graph, &separators);

    n = igraph_vector_int_list_size(&separators);
    for (i = 0; i < n; i++) {
        igraph_bool_t res;
        igraph_vector_int_t *sep = igraph_vector_int_list_get_ptr(&separators, i);

        igraph_is_separator(&graph, igraph_vss_vector(sep), &res);
        if (!res) {
            printf("Vertex set %" IGRAPH_PRId " is not a separator!\n", i);
            igraph_vector_int_print(sep);
            return 1;
        }
    }

    igraph_destroy(&graph);
    igraph_vector_int_list_destroy(&separators);

    return 0;
}
