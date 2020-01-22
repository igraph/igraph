/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2012  Gabor Csardi <csardi.gabor@gmail.com>
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
    igraph_t small, big;
    igraph_matrix_t small_coords, big_coords, merged_coords;
    igraph_vector_ptr_t graph_ptr, coords_ptr;
    igraph_arpack_options_t arpack_opts;

    /* To make things reproducible */
    igraph_rng_seed(igraph_rng_default(), 42);

    igraph_small(&big, 10, IGRAPH_UNDIRECTED,
                 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 0,
                 -1);

    igraph_small(&small, 3, IGRAPH_UNDIRECTED,
                 0, 1, 1, 2, 2, 0,
                 -1);

    igraph_arpack_options_init(&arpack_opts);

    igraph_matrix_init(&big_coords, 0, 0);
    igraph_layout_mds(&big, &big_coords, /*dist=*/ 0, /*dim=*/ 2,
                      &arpack_opts);

    igraph_matrix_init(&small_coords, 0, 0);
    igraph_layout_mds(&small, &small_coords, /*dist=*/ 0, /*dim=*/ 2,
                      &arpack_opts);

    igraph_vector_ptr_init(&graph_ptr, 2);
    igraph_vector_ptr_init(&coords_ptr, 2);
    igraph_matrix_init(&merged_coords, 0, 0);
    VECTOR(graph_ptr)[0] = &big;
    VECTOR(graph_ptr)[1] = &small;
    VECTOR(coords_ptr)[0] = &big_coords;
    VECTOR(coords_ptr)[1] = &small_coords;

    igraph_layout_merge_dla(&graph_ptr, &coords_ptr, &merged_coords);

    igraph_matrix_print(&merged_coords);

    igraph_matrix_destroy(&merged_coords);
    igraph_matrix_destroy(&small_coords);
    igraph_matrix_destroy(&big_coords);
    igraph_vector_ptr_destroy(&graph_ptr);
    igraph_vector_ptr_destroy(&coords_ptr);
    igraph_destroy(&small);
    igraph_destroy(&big);

    return 0;
}
