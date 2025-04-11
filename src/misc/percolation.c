/*
    IGraph library.
    Copyright (C) 2022  The igraph development team <igraph@igraph.org>

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

#include "igraph_components.h"
#include "igraph_error.h"
#include "igraph_structural.h"

static igraph_error_t igraph_i_percolate_edge(igraph_vector_int_t *links, igraph_vector_int_t *sizes, igraph_integer_t* biggest, igraph_integer_t a, igraph_integer_t b) {
    // find head of each
    // TODO: Path compression
    while (VECTOR(*links)[a] != a) {
        a = VECTOR(*links)[a];
    }
    while (VECTOR(*links)[b] != b) {
        b = VECTOR(*links)[b];
    }
    // make a child of b
    VECTOR(*links)[a] = b;
    VECTOR(*sizes)[b] += VECTOR(*sizes)[a];
    // if made new biggest component, update size
    if (VECTOR(*sizes)[b] >= *biggest) *biggest = VECTOR(*sizes)[b];
    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_edge_list_percolation(const igraph_vector_int_t *edges, igraph_vector_int_t* output) {
    igraph_integer_t biggest = 0;
    igraph_integer_t vert_count = igraph_vector_int_max(edges);
    igraph_vector_int_t sizes;
    IGRAPH_CHECK(igraph_vector_int_init(&sizes, vert_count));
    igraph_vector_int_t links;
    IGRAPH_CHECK(igraph_vector_int_init(&links, vert_count));

    int edge_count = igraph_vector_int_size(edges) / 2;
    IGRAPH_CHECK(igraph_vector_int_resize(output, edge_count));

    for (igraph_integer_t i = 0; i < edge_count; i++) {
        igraph_i_percolate_edge(&links, &sizes, &biggest, VECTOR(*edges)[2*i], VECTOR(*edges)[2*i+1]);
        VECTOR(*output)[i] = biggest;
    }
    igraph_vector_int_destroy(&sizes);
    igraph_vector_int_destroy(&links);

    return IGRAPH_SUCCESS;
}

IGRAPH_EXPORT igraph_error_t igraph_bond_percolation(const igraph_t *graph, igraph_vector_t * output) {
    return IGRAPH_SUCCESS;
}

IGRAPH_EXPORT igraph_error_t igraph_site_percolation(const igraph_t *graph, igraph_vector_t * output) {
    return IGRAPH_SUCCESS;
}

