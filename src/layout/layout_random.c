/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
/*
   IGraph library.
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
#include "igraph_random.h"

/**
 * \ingroup layout
 * \function igraph_layout_random
 * \brief Places the vertices uniform randomly on a plane.
 *
 * \param graph Pointer to an initialized graph object.
 * \param res Pointer to an initialized matrix object. This will
 *        contain the result and will be resized as needed.
 * \return Error code. The current implementation always returns with
 * success.
 *
 * Time complexity: O(|V|), the
 * number of vertices.
 */
int igraph_layout_random(const igraph_t *graph, igraph_matrix_t *res) {

    long int no_of_nodes = igraph_vcount(graph);
    long int i;

    IGRAPH_CHECK(igraph_matrix_resize(res, no_of_nodes, 2));

    RNG_BEGIN();

    for (i = 0; i < no_of_nodes; i++) {
        MATRIX(*res, i, 0) = RNG_UNIF(-1, 1);
        MATRIX(*res, i, 1) = RNG_UNIF(-1, 1);
    }

    RNG_END();

    return 0;
}

/**
 * \function igraph_layout_random_3d
 * \brief Places the vertices uniform randomly in a cube.
 *
 * </para><para>
 * Vertex coordinates range from -1 to 1, and are placed in 3 columns
 * of a matrix, with a row for each vertex.
 *
 * \param graph The graph to place.
 * \param res Pointer to an initialized matrix object. It will be
 * resized to hold the result.
 * \return Error code. The current implementation always returns with
 * success.
 *
 * Added in version 0.2.</para><para>
 *
 * Time complexity: O(|V|), the number of vertices.
 */
int igraph_layout_random_3d(const igraph_t *graph, igraph_matrix_t *res) {

    long int no_of_nodes = igraph_vcount(graph);
    long int i;

    IGRAPH_CHECK(igraph_matrix_resize(res, no_of_nodes, 3));

    RNG_BEGIN();

    for (i = 0; i < no_of_nodes; i++) {
        MATRIX(*res, i, 0) = RNG_UNIF(-1, 1);
        MATRIX(*res, i, 1) = RNG_UNIF(-1, 1);
        MATRIX(*res, i, 2) = RNG_UNIF(-1, 1);
    }

    RNG_END();

    return IGRAPH_SUCCESS;
}
