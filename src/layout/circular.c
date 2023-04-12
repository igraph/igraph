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

#include "core/interruption.h"
#include "core/math.h"

/**
 * \ingroup layout
 * \function igraph_layout_circle
 * \brief Places the vertices uniformly on a circle in arbitrary order.
 *
 * \param graph Pointer to an initialized graph object.
 * \param res Pointer to an initialized matrix object. This will
 *        contain the result and will be resized as needed.
 * \param order The order of the vertices on the circle. The vertices
 *        not included here, will be placed at (0,0). Supply
 *        \ref igraph_vss_all() here to place vertices in the
 *        order of their vertex IDs.
 * \return Error code.
 *
 * Time complexity: O(|V|), the number of vertices.
 */
igraph_error_t igraph_layout_circle(const igraph_t *graph, igraph_matrix_t *res,
                         igraph_vs_t order) {

    const igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t vs_size;
    igraph_vit_t vit;

    IGRAPH_CHECK(igraph_vs_size(graph, &order, &vs_size));

    IGRAPH_CHECK(igraph_matrix_resize(res, no_of_nodes, 2));
    igraph_matrix_null(res);

    IGRAPH_CHECK(igraph_vit_create(graph, order, &vit));
    for (igraph_integer_t i = 0; !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit), i++) {
        igraph_real_t phi = 2 * M_PI / vs_size * i;
        igraph_integer_t idx = IGRAPH_VIT_GET(vit);
        MATRIX(*res, idx, 0) = cos(phi);
        MATRIX(*res, idx, 1) = sin(phi);
    }
    igraph_vit_destroy(&vit);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_layout_star
 * \brief Generates a star-like layout.
 *
 * \param graph The input graph. Its edges are ignored by this function.
 * \param res Pointer to an initialized matrix object. This will
 *        contain the result and will be resized as needed.
 * \param center The id of the vertex to put in the center. You can set it to
 *        any arbitrary value for the special case when the input graph has no
 *        vertices; otherwise it must be between 0 and the number of vertices
 *        minus 1.
 * \param order A numeric vector giving the order of the vertices
 *      (including the center vertex!). If a null pointer, then the
 *      vertices are placed in increasing vertex ID order.
 * \return Error code.
 *
 * Time complexity: O(|V|), linear in the number of vertices.
 *
 * \sa \ref igraph_layout_circle() and other layout generators.
 */
igraph_error_t igraph_layout_star(const igraph_t *graph, igraph_matrix_t *res,
                       igraph_integer_t center, const igraph_vector_int_t *order) {

    const igraph_integer_t no_of_nodes = igraph_vcount(graph);

    if (no_of_nodes > 0 && (center < 0 || center >= no_of_nodes)) {
        IGRAPH_ERROR("The given center is not a vertex of the graph.", IGRAPH_EINVAL);
    }
    if (order && igraph_vector_int_size(order) != no_of_nodes) {
        IGRAPH_ERROR("Invalid order vector length.", IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_matrix_resize(res, no_of_nodes, 2));

    if (no_of_nodes == 1) {
        MATRIX(*res, 0, 0) = MATRIX(*res, 0, 1) = 0.0;
    } else if (no_of_nodes > 1) {
        igraph_real_t phi = 0.0;
        for (igraph_integer_t i = 0; i < no_of_nodes; i++) {
            igraph_integer_t node = order ? VECTOR(*order)[i] : i;
            if (order && (node < 0 || node >= no_of_nodes)) {
                IGRAPH_ERROR("Elements in the order vector are not all vertices of the graph.", IGRAPH_EINVAL);
            }
            if (node != center) {
                MATRIX(*res, node, 0) = cos(phi);
                MATRIX(*res, node, 1) = sin(phi);
                phi += 2.0 * M_PI / (no_of_nodes - 1);
            } else {
                MATRIX(*res, node, 0) = MATRIX(*res, node, 1) = 0.0;
            }
        }
    }

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_layout_sphere
 * \brief Places vertices (more or less) uniformly on a sphere.
 *
 * The vertices are placed with approximately equal spacing on a spiral
 * wrapped around a sphere, in the order of their vertex IDs. Vertices
 * with consecutive vertex IDs are placed near each other.
 *
 * </para><para>
 * The algorithm was described in the following paper:
 *
 * </para><para>
 * Distributing many points on a sphere by E.B. Saff and
 * A.B.J. Kuijlaars, \emb Mathematical Intelligencer \eme 19.1 (1997)
 * 5--11. https://doi.org/10.1007/BF03024331
 *
 * \param graph Pointer to an initialized graph object.
 * \param res Pointer to an initialized matrix object. This will
 *        contain the result and will be resized as needed.
 * \return Error code. The current implementation always returns with
 * success.
 *
 * Added in version 0.2.</para><para>
 *
 * Time complexity: O(|V|), the number of vertices in the graph.
 */
igraph_error_t igraph_layout_sphere(const igraph_t *graph, igraph_matrix_t *res) {

    const igraph_integer_t no_of_nodes = igraph_vcount(graph);
    const igraph_real_t sqrt_no_of_nodes = sqrt(no_of_nodes);

    IGRAPH_CHECK(igraph_matrix_resize(res, no_of_nodes, 3));

    igraph_real_t phi = 0;
    for (igraph_integer_t i = 0; i < no_of_nodes; i++) {
        igraph_real_t r, z;

        /* The first and last point are handled separately to avoid
         * division by zero or 1-z*z becoming slightly negative due
         * to roundoff errors. */
        if (i == 0) {
            z = -1; r = 0;
        } else if (i == no_of_nodes-1) {
            z = 1; r = 0;
        } else {
            z = -1.0 + 2.0 * i / (no_of_nodes - 1);
            r = sqrt(1 - z*z);
            phi += 3.6 / (sqrt_no_of_nodes*r);
        }

        igraph_real_t x = r*cos(phi);
        igraph_real_t y = r*sin(phi);

        MATRIX(*res, i, 0) = x;
        MATRIX(*res, i, 1) = y;
        MATRIX(*res, i, 2) = z;

        IGRAPH_ALLOW_INTERRUPTION();
    }

    return IGRAPH_SUCCESS;
}
