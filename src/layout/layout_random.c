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

#include "layout/layout_internal.h"

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
igraph_error_t igraph_layout_random(const igraph_t *graph, igraph_matrix_t *res) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t i;

    IGRAPH_CHECK(igraph_matrix_resize(res, no_of_nodes, 2));

    RNG_BEGIN();

    for (i = 0; i < no_of_nodes; i++) {
        MATRIX(*res, i, 0) = RNG_UNIF(-1, 1);
        MATRIX(*res, i, 1) = RNG_UNIF(-1, 1);
    }

    RNG_END();

    return IGRAPH_SUCCESS;
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
igraph_error_t igraph_layout_random_3d(const igraph_t *graph, igraph_matrix_t *res) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t i;

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


/* The following functions generate suitable initial random layouts for
 * the Fruchterman-Reingold and Kamada-Kawai algorithms. */

igraph_error_t igraph_i_layout_random_bounded(
        const igraph_t *graph,
        igraph_matrix_t *res,
        const igraph_vector_t *minx, const igraph_vector_t *maxx,
        const igraph_vector_t *miny, const igraph_vector_t *maxy) {

    const igraph_integer_t no_nodes = igraph_vcount(graph);
    const igraph_real_t width = sqrt(no_nodes), height = width;

    igraph_real_t dminx = -width/2,  dmaxx = width/2,
                  dminy = -height/2, dmaxy = height/2; /* default values */

    /* Caller should ensure that minx, etc. do not contain NaN. */

    if (minx && !igraph_vector_empty(minx)) {
        igraph_real_t m = igraph_vector_max(minx);
        if (m == IGRAPH_POSINFINITY) {
            IGRAPH_ERROR("Infinite lower coordinate bound for graph layout.", IGRAPH_EINVAL);
        }
        if (m > dmaxx) {
            dmaxx += m;
        }
    }
    if (maxx && !igraph_vector_empty(maxx)) {
        igraph_real_t m = igraph_vector_min(maxx);
        if (m == IGRAPH_NEGINFINITY) {
            IGRAPH_ERROR("Negative infinite upper coordinate bound for graph layout.", IGRAPH_EINVAL);
        }
        if (m < dminx) {
            dminx -= m;
        }
    }
    if (miny && !igraph_vector_empty(miny)) {
        igraph_real_t m = igraph_vector_max(miny);
        if (m == IGRAPH_POSINFINITY) {
            IGRAPH_ERROR("Infinite lower coordinate bound for graph layout.", IGRAPH_EINVAL);
        }
        if (m > dmaxy) {
            dmaxy += m;
        }
    }
    if (maxy && !igraph_vector_empty(maxy)) {
        igraph_real_t m = igraph_vector_min(maxy);
        if (m == IGRAPH_NEGINFINITY) {
            IGRAPH_ERROR("Negative infinite upper coordinate bound for graph layout.", IGRAPH_EINVAL);
        }
        if (m < dminy) {
            dminy -= m;
        }
    }

    RNG_BEGIN();
    IGRAPH_CHECK(igraph_matrix_resize(res, no_nodes, 2));
    for (igraph_integer_t i = 0; i < no_nodes; i++) {
        igraph_real_t x1 = minx ? VECTOR(*minx)[i] : dminx;
        igraph_real_t x2 = maxx ? VECTOR(*maxx)[i] : dmaxx;
        igraph_real_t y1 = miny ? VECTOR(*miny)[i] : dminy;
        igraph_real_t y2 = maxy ? VECTOR(*maxy)[i] : dmaxy;
        if (!isfinite(x1)) {
            x1 = -width / 2;
        }
        if (!isfinite(x2)) {
            x2 =  width / 2;
        }
        if (!isfinite(y1)) {
            y1 = -height / 2;
        }
        if (!isfinite(y2)) {
            y2 =  height / 2;
        }
        MATRIX(*res, i, 0) = RNG_UNIF(x1, x2);
        MATRIX(*res, i, 1) = RNG_UNIF(y1, y2);
    }
    RNG_END();

    return IGRAPH_SUCCESS;
}

igraph_error_t igraph_i_layout_random_bounded_3d(
        const igraph_t *graph, igraph_matrix_t *res,
        const igraph_vector_t *minx, const igraph_vector_t *maxx,
        const igraph_vector_t *miny, const igraph_vector_t *maxy,
        const igraph_vector_t *minz, const igraph_vector_t *maxz) {

    const igraph_integer_t no_nodes = igraph_vcount(graph);
    const igraph_real_t width = sqrt(no_nodes), height = width, depth = width;

    igraph_real_t dminx = -width/2,  dmaxx = width/2,
                  dminy = -height/2, dmaxy = height/2,
                  dminz = -depth/2,  dmaxz = depth/2; /* default values */

    /* Caller should ensure that minx, etc. do not contain NaN. */

    if (minx && !igraph_vector_empty(minx)) {
        igraph_real_t m = igraph_vector_max(minx);
        if (m == IGRAPH_POSINFINITY) {
            IGRAPH_ERROR("Infinite lower coordinate bound for graph layout.", IGRAPH_EINVAL);
        }
        if (m > dmaxx) {
            dmaxx += m;
        }
    }
    if (maxx && !igraph_vector_empty(maxx)) {
        igraph_real_t m = igraph_vector_min(maxx);
        if (m == IGRAPH_NEGINFINITY) {
            IGRAPH_ERROR("Negative infinite upper coordinate bound for graph layout.", IGRAPH_EINVAL);
        }
        if (m < dminx) {
            dminx -= m;
        }
    }
    if (miny && !igraph_vector_empty(miny)) {
        igraph_real_t m = igraph_vector_max(miny);
        if (m == IGRAPH_POSINFINITY) {
            IGRAPH_ERROR("Infinite lower coordinate bound for graph layout.", IGRAPH_EINVAL);
        }
        if (m > dmaxy) {
            dmaxy += m;
        }
    }
    if (maxy && !igraph_vector_empty(maxy)) {
        igraph_real_t m = igraph_vector_min(maxy);
        if (m == IGRAPH_NEGINFINITY) {
            IGRAPH_ERROR("Negative infinite upper coordinate bound for graph layout.", IGRAPH_EINVAL);
        }
        if (m < dminy) {
            dminy -= m;
        }
    }
    if (minz && !igraph_vector_empty(minz)) {
        igraph_real_t m = igraph_vector_max(minz);
        if (m == IGRAPH_POSINFINITY) {
            IGRAPH_ERROR("Infinite lower coordinate bound for graph layout.", IGRAPH_EINVAL);
        }
        if (m > dmaxz) {
            dmaxz += m;
        }
    }
    if (maxz && !igraph_vector_empty(maxz)) {
        igraph_real_t m = igraph_vector_min(maxz);
        if (m == IGRAPH_NEGINFINITY) {
            IGRAPH_ERROR("Negative infinite upper coordinate bound for graph layout.", IGRAPH_EINVAL);
        }
        if (m < dminz) {
            dminz -= m;
        }
    }

    RNG_BEGIN();
    IGRAPH_CHECK(igraph_matrix_resize(res, no_nodes, 3));
    for (igraph_integer_t i = 0; i < no_nodes; i++) {
        igraph_real_t x1 = minx ? VECTOR(*minx)[i] : dminx;
        igraph_real_t x2 = maxx ? VECTOR(*maxx)[i] : dmaxx;
        igraph_real_t y1 = miny ? VECTOR(*miny)[i] : dminy;
        igraph_real_t y2 = maxy ? VECTOR(*maxy)[i] : dmaxy;
        igraph_real_t z1 = minz ? VECTOR(*minz)[i] : dminz;
        igraph_real_t z2 = maxz ? VECTOR(*maxz)[i] : dmaxz;
        if (!isfinite(x1)) {
            x1 = -width / 2;
        }
        if (!isfinite(x2)) {
            x2 =  width / 2;
        }
        if (!isfinite(y1)) {
            y1 = -height / 2;
        }
        if (!isfinite(y2)) {
            y2 =  height / 2;
        }
        if (!isfinite(z1)) {
            z1 = -depth / 2;
        }
        if (!isfinite(z2)) {
            z2 =  depth / 2;
        }
        MATRIX(*res, i, 0) = RNG_UNIF(x1, x2);
        MATRIX(*res, i, 1) = RNG_UNIF(y1, y2);
        MATRIX(*res, i, 2) = RNG_UNIF(z1, z2);
    }
    RNG_END();

    return IGRAPH_SUCCESS;
}
