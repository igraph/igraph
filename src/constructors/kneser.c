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

static igraph_integer_t igraph_i_number_of_subsets(igraph_integer_t n, igraph_integer_t k) {
    if (n < k) {
        return 0;
    }
    if (k == 0) {
        return 1;
    }

    return igraph_i_number_of_subsets(n - 1, k) + igraph_i_number_of_subsets(n - 1, k - 1);
}

static igraph_error_t igraph_i_generate_subsets(igraph_integer_t n, igraph_integer_t k, igraph_integer_t vertex_id, igraph_integer_t size,
                              igraph_integer_t *subset_index, igraph_vector_int_t *subset, igraph_vector_int_list_t *vertices) {
    if (size == k) {
        IGRAPH_CHECK(igraph_vector_list_insert_copy(vertices, *subset_index, subset));
        (*subset_index)++;
        return IGRAPH_SUCCESS;
    }
    if (vertex_id > n) {
        return IGRAPH_SUCCESS;
    }

    VECTOR(*subset)[size] = vertex_id;
    IGRAPH_CHECK(igraph_i_generate_subsets(n, k, vertex_id + 1, size + 1, subset_index, subset, vertices));
    IGRAPH_CHECK(igraph_i_generate_subsets(n, k, vertex_id + 1, size, subset_index, subset, vertices));
    return IGRAPH_SUCCESS;
}