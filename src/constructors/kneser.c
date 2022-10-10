/*
   IGraph library.
   Copyright (C) 2022  The igraph development team

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

#include "igraph_vector_list.h"

igraph_vector_int_list_t vertices;

static igraph_integer_t igraph_i_number_of_subsets(igraph_integer_t n, igraph_integer_t k){
    if (n < k) {
        return 0;
    }
    if (k == 0) {
        return 1;
    }

    igraph_integer_t nCk = igraph_number_of_subsets(n - 1, k) + igraph_number_of_subsets(n - 1, k - 1);
    return nCk;
}

static igraph_i_generate_subsets(igraph_integer_t n, igraph_integer_t k, igraph_integer_t vertex_id, igraph_integer_t size, igraph_integer_t *subset_index, igraph_vector_int_t *subset){
    if (size == k){
        igraph_vector_list_insert_copy(&vertices, *subset_index, subset);
        (*subset_index)++;
        return;
    }
    if (vertex_id > n) return;
    VECTOR(*subset)[size] = vertex_id;
    genSubset(n, k, vertex_id + 1, size + 1, subset_index, subset);
    genSubset(n, k, vertex_id + 1, size, subset_index, subset);
}