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
 * \brief Places the vertices uniformly on a circle, in the order of vertex ids.
 *
 * \param graph Pointer to an initialized graph object.
 * \param res Pointer to an initialized matrix object. This will
 *        contain the result and will be resized as needed.
 * \param order The order of the vertices on the circle. The vertices
 *        not included here, will be placed at (0,0). Supply
 *        \ref igraph_vss_all() here for all vertices, in the order of
 *        their vertex ids.
 * \return Error code.
 *
 * Time complexity: O(|V|), the
 * number of vertices.
 */
int igraph_layout_circle(const igraph_t *graph, igraph_matrix_t *res,
                         igraph_vs_t order) {

    long int no_of_nodes = igraph_vcount(graph);
    igraph_integer_t vs_size;
    long int i;
    igraph_vit_t vit;

    IGRAPH_CHECK(igraph_vs_size(graph, &order, &vs_size));

    IGRAPH_CHECK(igraph_matrix_resize(res, no_of_nodes, 2));
    igraph_matrix_null(res);

    igraph_vit_create(graph, order, &vit);
    for (i = 0; !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit), i++) {
        igraph_real_t phi = 2 * M_PI / vs_size * i;
        int idx = IGRAPH_VIT_GET(vit);
        MATRIX(*res, idx, 0) = cos(phi);
        MATRIX(*res, idx, 1) = sin(phi);
    }
    igraph_vit_destroy(&vit);

    return 0;
}

/**
 * \function igraph_layout_star
 * \brief Generates a star-like layout.
 *
 * \param graph The input graph. Its edges are ignored by this function.
 * \param res Pointer to an initialized matrix object. This will
 *        contain the result and will be resized as needed.
 * \param center The id of the vertex to put in the center.
 * \param order A numeric vector giving the order of the vertices
 *      (including the center vertex!). If a null pointer, then the
 *      vertices are placed in increasing vertex id order.
 * \return Error code.
 *
 * Time complexity: O(|V|), linear in the number of vertices.
 *
 * \sa \ref igraph_layout_circle() and other layout generators.
 */
int igraph_layout_star(const igraph_t *graph, igraph_matrix_t *res,
                       igraph_integer_t center, const igraph_vector_t *order) {

    long int no_of_nodes = igraph_vcount(graph);
    long int c = center;
    long int i;
    igraph_real_t step;
    igraph_real_t phi;

    if (center < 0 || center >= no_of_nodes) {
        IGRAPH_ERROR("The given center is not a vertex of the graph.", IGRAPH_EINVAL);
    }
    if (order && igraph_vector_size(order) != no_of_nodes) {
        IGRAPH_ERROR("Invalid order vector length.", IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_matrix_resize(res, no_of_nodes, 2));

    if (no_of_nodes == 1) {
        MATRIX(*res, 0, 0) = MATRIX(*res, 0, 1) = 0.0;
    } else {
        for (i = 0, step = 2 * M_PI / (no_of_nodes - 1), phi = 0;
             i < no_of_nodes; i++) {
            long int node = order ? (long int) VECTOR(*order)[i] : i;
            if (order && (node < 0 || node >= no_of_nodes)) {
                IGRAPH_ERROR("Elements in the order vector are not all vertices of the graph.", IGRAPH_EINVAL);
            }
            if (node != c) {
                MATRIX(*res, node, 0) = cos(phi);
                MATRIX(*res, node, 1) = sin(phi);
                phi += step;
            } else {
                MATRIX(*res, node, 0) = MATRIX(*res, node, 1) = 0.0;
            }
        }
    }

    return 0;
}

/**
 * \function igraph_layout_sphere
 * \brief Places vertices (more or less) uniformly on a sphere.
 *
 * </para><para>
 * The algorithm was described in the following paper:
 * Distributing many points on a sphere by E.B. Saff and
 * A.B.J. Kuijlaars, \emb Mathematical Intelligencer \eme 19.1 (1997)
 * 5--11.
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
int igraph_layout_sphere(const igraph_t *graph, igraph_matrix_t *res) {

    long int no_of_nodes = igraph_vcount(graph);
    long int i;
    igraph_real_t h;

    IGRAPH_CHECK(igraph_matrix_resize(res, no_of_nodes, 3));

    if (no_of_nodes != 0) {
        MATRIX(*res, 0, 0) = M_PI;
        MATRIX(*res, 0, 1) = 0;
    }
    for (i = 1; i < no_of_nodes - 1; i++) {
        h = -1 + 2 * i / (double)(no_of_nodes - 1);
        MATRIX(*res, i, 0) = acos(h);
        MATRIX(*res, i, 1) = fmod((MATRIX(*res, i - 1, 1) +
                                   3.6 / sqrt(no_of_nodes * (1 - h * h))), 2 * M_PI);
        IGRAPH_ALLOW_INTERRUPTION();
    }
    if (no_of_nodes >= 2) {
        MATRIX(*res, no_of_nodes - 1, 0) = 0;
        MATRIX(*res, no_of_nodes - 1, 1) = 0;
    }

    for (i = 0; i < no_of_nodes; i++) {
        igraph_real_t x = cos(MATRIX(*res, i, 1)) * sin(MATRIX(*res, i, 0));
        igraph_real_t y = sin(MATRIX(*res, i, 1)) * sin(MATRIX(*res, i, 0));
        igraph_real_t z = cos(MATRIX(*res, i, 0));
        MATRIX(*res, i, 0) = x;
        MATRIX(*res, i, 1) = y;
        MATRIX(*res, i, 2) = z;
        IGRAPH_ALLOW_INTERRUPTION();
    }

    return 0;
}
