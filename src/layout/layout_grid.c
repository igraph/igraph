/*
   igraph library.
   Copyright (C) 2003-2020  The igraph development team

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

#include "igraph_layout.h"

#include "igraph_interface.h"

/**
 * \ingroup layout
 * \function igraph_layout_grid
 * \brief Places the vertices on a regular grid on the plane.
 *
 * \param graph Pointer to an initialized graph object.
 * \param res Pointer to an initialized matrix object. This will
 *        contain the result and will be resized as needed.
 * \param width The number of vertices in a single row of the grid.
 *        When zero or negative, the width of the grid will be the
 *        square root of the number of vertices, rounded up if needed.
 * \return Error code. The current implementation always returns with
 *         success.
 *
 * Time complexity: O(|V|), the number of vertices.
 */
igraph_error_t igraph_layout_grid(const igraph_t *graph, igraph_matrix_t *res, igraph_int_t width) {
    igraph_int_t i, no_of_nodes = igraph_vcount(graph);
    igraph_real_t x, y;

    IGRAPH_CHECK(igraph_matrix_resize(res, no_of_nodes, 2));

    if (width <= 0) {
        width = ceil(sqrt(no_of_nodes));
    }

    x = y = 0;
    for (i = 0; i < no_of_nodes; i++) {
        MATRIX(*res, i, 0) = x++;
        MATRIX(*res, i, 1) = y;
        if (x == width) {
            x = 0; y++;
        }
    }

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup layout
 * \function igraph_layout_grid_3d
 * \brief Places the vertices on a regular grid in the 3D space.
 *
 * \param graph Pointer to an initialized graph object.
 * \param res Pointer to an initialized matrix object. This will
 *        contain the result and will be resized as needed.
 * \param width  The number of vertices in a single row of the grid. When
 *               zero or negative, the width is determined automatically.
 * \param height The number of vertices in a single column of the grid. When
 *               zero or negative, the height is determined automatically.
 *
 * \return Error code. The current implementation always returns with
 *         success.
 *
 * Time complexity: O(|V|), the number of vertices.
 */
igraph_error_t igraph_layout_grid_3d(const igraph_t *graph, igraph_matrix_t *res,
                          igraph_int_t width, igraph_int_t height) {
    igraph_int_t i, no_of_nodes = igraph_vcount(graph);
    igraph_real_t x, y, z;

    IGRAPH_CHECK(igraph_matrix_resize(res, no_of_nodes, 3));

    if (width <= 0 && height <= 0) {
        width = height = ceil(pow(no_of_nodes, 1.0 / 3));
    } else if (width <= 0) {
        width = ceil(sqrt(no_of_nodes / (double)height));
    } else if (height <= 0) {
        height = ceil(sqrt(no_of_nodes / (double)width));
    }

    x = y = z = 0;
    for (i = 0; i < no_of_nodes; i++) {
        MATRIX(*res, i, 0) = x++;
        MATRIX(*res, i, 1) = y;
        MATRIX(*res, i, 2) = z;
        if (x == width) {
            x = 0; y++;
            if (y == height) {
                y = 0; z++;
            }
        }
    }

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup layout
 * \function igraph_layout_square
 * \brief Places the vertices on a regular grid in the d-dimensional space.
 *
 * \param graph Pointer to an initialized graph object.
 * \param res Pointer to an initialized matrix object. This will
 *        contain the result and will be resized as needed.
 *
 * \return Error code. The current implementation always returns with
 *         success.
 *
 * Time complexity: TODO: ADD TIMECOMPLEXITY
 */
igraph_error_t igraph_layout_square(const igraph_t *graph, igraph_matrix_t *res, igraph_vector_int_t *dimvector) {

    igraph_int_t i, j, no_of_nodes;
    igraph_int_t no_of_dims = igraph_vector_int_size(dimvector);
    igraph_vector_t ith_node_coords;

    no_of_nodes = 1;
    for (i = 0; i < no_of_dims; i++) {
        no_of_nodes *= VECTOR(*dimvector)[i];
    }

    if (graph != NULL) {
        if (no_of_nodes != igraph_vcount(graph)) {
            IGRAPH_ERROR("Number of graph vertices does not equal to number of nodes " 
                         "given lattice dimensionality", IGRAPH_EINVAL);
        }
    }

    IGRAPH_CHECK(igraph_matrix_resize(res, no_of_nodes, no_of_dims));

    igraph_vector_init(&ith_node_coords, no_of_dims);
    for (i = 0; i < no_of_nodes; i++) {
        igraph_matrix_set_row(res, &ith_node_coords, i);
        for (j = 0; j < no_of_dims; j++) {
            VECTOR(ith_node_coords)[j]++;
            if (VECTOR(ith_node_coords)[j] < VECTOR(*dimvector)[j]) {
                break;
            }
            VECTOR(ith_node_coords)[j] = 0;
        }
    }

    return IGRAPH_SUCCESS;
}