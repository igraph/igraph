/*
   igraph library.
   Copyright (C) 2026 The igraph development team <igraph@igraph.org>

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include "igraph_layout.h"

#include "igraph_interface.h"
#include "math/safe_intop.h"


igraph_error_t triangular_lattice_layout_rectangle_shape(const igraph_t *graph, igraph_matrix_t *res, const igraph_vector_int_t *dims) {
    igraph_int_t i, no_of_nodes;
    igraph_int_t x, y;

    if (igraph_vector_int_size(dims) != 2) {
        IGRAPH_ERROR("Rectangular lattice shape layout requires a dimension vector of length 2.",
                IGRAPH_EINVAL);
    }

    if (igraph_vector_int_any_smaller(dims, 0)) {
        IGRAPH_ERROR("Invalid dimension vector.", IGRAPH_EINVAL);
    }

   IGRAPH_CHECK(igraph_i_safe_vector_int_prod(dims, &no_of_nodes));

    if (graph != NULL) {
        if (no_of_nodes != igraph_vcount(graph)) {
            IGRAPH_ERROR("Number of graph vertices does not equal to number of nodes "
                         "given lattice dimensionality.", IGRAPH_EINVAL);
        }
    }

    IGRAPH_CHECK(igraph_matrix_resize(res, no_of_nodes, 2));

    const igraph_int_t n_elem_first_dim = VECTOR(*dims)[0];
    x = y = 0;
    for (i = 0; i < no_of_nodes; i++) {
        MATRIX(*res, i, 0) = x + 0.5 * y;
        MATRIX(*res, i, 1) = y * sqrt(3.0) / 2.0;
        x++;
        if (x == n_elem_first_dim) {
            x = 0;
            y++;
        }
    }

    return IGRAPH_SUCCESS;
}



igraph_error_t igraph_layout_triangular(igraph_t *graph, igraph_matrix_t *res, const igraph_vector_int_t *dims) {

    igraph_int_t num_dims = igraph_vector_int_size(dims);

    if (igraph_vector_int_any_smaller(dims, 0)) {
        IGRAPH_ERROR("Invalid dimension vector.", IGRAPH_EINVAL);
    }
    /* If a coordinate of dims is 0 the result is an empty graph.
    Should this be part of the function as well? If so, 
    should igraph_bool_t directed be added to function parameters? */
    if (igraph_vector_int_contains(dims, 0)) {
        return igraph_empty(graph, 0, false);
    }

    switch (num_dims) {
    case 1:
        // IGRAPH_CHECK(triangular_lattice_layout_triangle_shape(graph, res, dims));
        break;
    case 2:
        IGRAPH_CHECK(triangular_lattice_layout_rectangle_shape(graph, res, dims));
        break;
    case 3:
        // IGRAPH_CHECK(triangular_lattice_layout_hex_shape(graph, res, dims));
        break;
    default:
        IGRAPH_ERRORF(
            "The size of the dimension vector must be 1, 2 or 3, got %" IGRAPH_PRId ".",
            IGRAPH_EINVAL, num_dims);
    }

    return IGRAPH_SUCCESS;
}


