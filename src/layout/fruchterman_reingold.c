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

#include "igraph_random.h"
#include "igraph_interface.h"
#include "igraph_components.h"

#include "core/grid.h"

static int igraph_layout_i_fr(const igraph_t *graph,
                              igraph_matrix_t *res,
                              igraph_bool_t use_seed,
                              igraph_integer_t niter,
                              igraph_real_t start_temp,
                              const igraph_vector_t *weight,
                              const igraph_vector_t *minx,
                              const igraph_vector_t *maxx,
                              const igraph_vector_t *miny,
                              const igraph_vector_t *maxy) {

    igraph_integer_t no_nodes = igraph_vcount(graph);
    igraph_integer_t no_edges = igraph_ecount(graph);
    igraph_integer_t i;
    igraph_vector_float_t dispx, dispy;
    igraph_real_t temp = start_temp;
    igraph_real_t difftemp = start_temp / niter;
    float width = sqrtf(no_nodes), height = width;
    igraph_bool_t conn = 1;
    float C = 0;

    igraph_is_connected(graph, &conn, IGRAPH_WEAK);
    if (!conn) {
        C = no_nodes * sqrtf(no_nodes);
    }

    RNG_BEGIN();

    if (!use_seed) {
        IGRAPH_CHECK(igraph_matrix_resize(res, no_nodes, 2));
        for (i = 0; i < no_nodes; i++) {
            igraph_real_t x1 = minx ? VECTOR(*minx)[i] : -width / 2;
            igraph_real_t x2 = maxx ? VECTOR(*maxx)[i] :  width / 2;
            igraph_real_t y1 = miny ? VECTOR(*miny)[i] : -height / 2;
            igraph_real_t y2 = maxy ? VECTOR(*maxy)[i] :  height / 2;
            if (!igraph_finite(x1)) {
                x1 = -sqrt(no_nodes) / 2;
            }
            if (!igraph_finite(x2)) {
                x2 =  sqrt(no_nodes) / 2;
            }
            if (!igraph_finite(y1)) {
                y1 = -sqrt(no_nodes) / 2;
            }
            if (!igraph_finite(y2)) {
                y2 =  sqrt(no_nodes) / 2;
            }
            MATRIX(*res, i, 0) = RNG_UNIF(x1, x2);
            MATRIX(*res, i, 1) = RNG_UNIF(y1, y2);
        }
    }

    IGRAPH_CHECK(igraph_vector_float_init(&dispx, no_nodes));
    IGRAPH_FINALLY(igraph_vector_float_destroy, &dispx);
    IGRAPH_CHECK(igraph_vector_float_init(&dispy, no_nodes));
    IGRAPH_FINALLY(igraph_vector_float_destroy, &dispy);

    for (i = 0; i < niter; i++) {
        igraph_integer_t v, u, e;

        /* calculate repulsive forces, we have a special version
           for unconnected graphs */
        igraph_vector_float_null(&dispx);
        igraph_vector_float_null(&dispy);
        if (conn) {
            for (v = 0; v < no_nodes; v++) {
                for (u = v + 1; u < no_nodes; u++) {
                    float dx = MATRIX(*res, v, 0) - MATRIX(*res, u, 0);
                    float dy = MATRIX(*res, v, 1) - MATRIX(*res, u, 1);
                    float dlen = dx * dx + dy * dy;

                    if (dlen == 0) {
                        dx = RNG_UNIF01() * 1e-9;
                        dy = RNG_UNIF01() * 1e-9;
                        dlen = dx * dx + dy * dy;
                    }

                    VECTOR(dispx)[v] += dx / dlen;
                    VECTOR(dispy)[v] += dy / dlen;
                    VECTOR(dispx)[u] -= dx / dlen;
                    VECTOR(dispy)[u] -= dy / dlen;
                }
            }
        } else {
            for (v = 0; v < no_nodes; v++) {
                for (u = v + 1; u < no_nodes; u++) {
                    float dx = MATRIX(*res, v, 0) - MATRIX(*res, u, 0);
                    float dy = MATRIX(*res, v, 1) - MATRIX(*res, u, 1);
                    float dlen, rdlen;

                    dlen = dx * dx + dy * dy;
                    if (dlen == 0) {
                        dx = RNG_UNIF(0, 1e-6);
                        dy = RNG_UNIF(0, 1e-6);
                        dlen = dx * dx + dy * dy;
                    }

                    rdlen = sqrt(dlen);

                    VECTOR(dispx)[v] += dx * (C - dlen * rdlen) / (dlen * C);
                    VECTOR(dispy)[v] += dy * (C - dlen * rdlen) / (dlen * C);
                    VECTOR(dispx)[u] -= dx * (C - dlen * rdlen) / (dlen * C);
                    VECTOR(dispy)[u] -= dy * (C - dlen * rdlen) / (dlen * C);
                }
            }
        }

        /* calculate attractive forces */
        for (e = 0; e < no_edges; e++) {
            /* each edges is an ordered pair of vertices v and u */
            igraph_integer_t v = IGRAPH_FROM(graph, e);
            igraph_integer_t u = IGRAPH_TO(graph, e);
            igraph_real_t dx = MATRIX(*res, v, 0) - MATRIX(*res, u, 0);
            igraph_real_t dy = MATRIX(*res, v, 1) - MATRIX(*res, u, 1);
            igraph_real_t w = weight ? VECTOR(*weight)[e] : 1.0;
            igraph_real_t dlen = sqrt(dx * dx + dy * dy) * w;
            VECTOR(dispx)[v] -= (dx * dlen);
            VECTOR(dispy)[v] -= (dy * dlen);
            VECTOR(dispx)[u] += (dx * dlen);
            VECTOR(dispy)[u] += (dy * dlen);
        }

        /* limit max displacement to temperature t and prevent from
           displacement outside frame */
        for (v = 0; v < no_nodes; v++) {
            igraph_real_t dx = VECTOR(dispx)[v] + RNG_UNIF01() * 1e-9;
            igraph_real_t dy = VECTOR(dispy)[v] + RNG_UNIF01() * 1e-9;
            igraph_real_t displen = sqrt(dx * dx + dy * dy);
            igraph_real_t mx = fabs(dx) < temp ? dx : temp;
            igraph_real_t my = fabs(dy) < temp ? dy : temp;
            if (displen > 0) {
                MATRIX(*res, v, 0) += (dx / displen) * mx;
                MATRIX(*res, v, 1) += (dy / displen) * my;
            }
            if (minx && MATRIX(*res, v, 0) < VECTOR(*minx)[v]) {
                MATRIX(*res, v, 0) = VECTOR(*minx)[v];
            }
            if (maxx && MATRIX(*res, v, 0) > VECTOR(*maxx)[v]) {
                MATRIX(*res, v, 0) = VECTOR(*maxx)[v];
            }
            if (miny && MATRIX(*res, v, 1) < VECTOR(*miny)[v]) {
                MATRIX(*res, v, 1) = VECTOR(*miny)[v];
            }
            if (maxy && MATRIX(*res, v, 1) > VECTOR(*maxy)[v]) {
                MATRIX(*res, v, 1) = VECTOR(*maxy)[v];
            }
        }

        temp -= difftemp;
    }

    RNG_END();

    igraph_vector_float_destroy(&dispx);
    igraph_vector_float_destroy(&dispy);
    IGRAPH_FINALLY_CLEAN(2);

    return 0;
}

static int igraph_layout_i_grid_fr(
        const igraph_t *graph,
        igraph_matrix_t *res, igraph_bool_t use_seed,
        igraph_integer_t niter, igraph_real_t start_temp,
        const igraph_vector_t *weight, const igraph_vector_t *minx,
        const igraph_vector_t *maxx, const igraph_vector_t *miny,
        const igraph_vector_t *maxy) {

    igraph_integer_t no_nodes = igraph_vcount(graph);
    igraph_integer_t no_edges = igraph_ecount(graph);
    float width = sqrtf(no_nodes), height = width;
    igraph_2dgrid_t grid;
    igraph_vector_float_t dispx, dispy;
    igraph_real_t temp = start_temp;
    igraph_real_t difftemp = start_temp / niter;
    igraph_2dgrid_iterator_t vidit;
    igraph_integer_t i;
    const float cellsize = 2.0;

    RNG_BEGIN();

    if (!use_seed) {
        IGRAPH_CHECK(igraph_matrix_resize(res, no_nodes, 2));
        for (i = 0; i < no_nodes; i++) {
            igraph_real_t x1 = minx ? VECTOR(*minx)[i] : -width / 2;
            igraph_real_t x2 = maxx ? VECTOR(*maxx)[i] :  width / 2;
            igraph_real_t y1 = miny ? VECTOR(*miny)[i] : -height / 2;
            igraph_real_t y2 = maxy ? VECTOR(*maxy)[i] :  height / 2;
            if (!igraph_finite(x1)) {
                x1 = -sqrt(no_nodes) / 2;
            }
            if (!igraph_finite(x2)) {
                x2 =  sqrt(no_nodes) / 2;
            }
            if (!igraph_finite(y1)) {
                y1 = -sqrt(no_nodes) / 2;
            }
            if (!igraph_finite(y2)) {
                y2 =  sqrt(no_nodes) / 2;
            }
            MATRIX(*res, i, 0) = RNG_UNIF(x1, x2);
            MATRIX(*res, i, 1) = RNG_UNIF(y1, y2);
        }
    }

    /* make grid */
    IGRAPH_CHECK(igraph_2dgrid_init(&grid, res, -width / 2, width / 2, cellsize,
                                    -height / 2, height / 2, cellsize));
    IGRAPH_FINALLY(igraph_2dgrid_destroy, &grid);

    /* place vertices on grid */
    for (i = 0; i < no_nodes; i++) {
        igraph_2dgrid_add2(&grid, i);
    }

    IGRAPH_CHECK(igraph_vector_float_init(&dispx, no_nodes));
    IGRAPH_FINALLY(igraph_vector_float_destroy, &dispx);
    IGRAPH_CHECK(igraph_vector_float_init(&dispy, no_nodes));
    IGRAPH_FINALLY(igraph_vector_float_destroy, &dispy);

    for (i = 0; i < niter; i++) {
        igraph_integer_t v, u, e;

        igraph_vector_float_null(&dispx);
        igraph_vector_float_null(&dispy);

        /* repulsion */
        igraph_2dgrid_reset(&grid, &vidit);
        while ( (v = igraph_2dgrid_next(&grid, &vidit) - 1) != -1) {
            while ( (u = igraph_2dgrid_next_nei(&grid, &vidit) - 1) != -1) {
                float dx = MATRIX(*res, v, 0) - MATRIX(*res, u, 0);
                float dy = MATRIX(*res, v, 1) - MATRIX(*res, u, 1);
                float dlen = dx * dx + dy * dy;
                if (dlen < cellsize * cellsize) {
                    VECTOR(dispx)[v] += dx / dlen;
                    VECTOR(dispy)[v] += dy / dlen;
                    VECTOR(dispx)[u] -= dx / dlen;
                    VECTOR(dispy)[u] -= dy / dlen;
                }
            }
        }

        /* attraction */
        for (e = 0; e < no_edges; e++) {
            igraph_integer_t v = IGRAPH_FROM(graph, e);
            igraph_integer_t u = IGRAPH_TO(graph, e);
            igraph_real_t dx = MATRIX(*res, v, 0) - MATRIX(*res, u, 0);
            igraph_real_t dy = MATRIX(*res, v, 1) - MATRIX(*res, u, 1);
            igraph_real_t w = weight ? VECTOR(*weight)[e] : 1.0;
            igraph_real_t dlen = sqrt(dx * dx + dy * dy) * w;
            VECTOR(dispx)[v] -= (dx * dlen);
            VECTOR(dispy)[v] -= (dy * dlen);
            VECTOR(dispx)[u] += (dx * dlen);
            VECTOR(dispy)[u] += (dy * dlen);
        }

        /* update */
        for (v = 0; v < no_nodes; v++) {
            igraph_real_t dx = VECTOR(dispx)[v] + RNG_UNIF01() * 1e-9;
            igraph_real_t dy = VECTOR(dispy)[v] + RNG_UNIF01() * 1e-9;
            igraph_real_t displen = sqrt(dx * dx + dy * dy);
            igraph_real_t mx = fabs(dx) < temp ? dx : temp;
            igraph_real_t my = fabs(dy) < temp ? dy : temp;
            if (displen > 0) {
                MATRIX(*res, v, 0) += (dx / displen) * mx;
                MATRIX(*res, v, 1) += (dy / displen) * my;
            }
            if (minx && MATRIX(*res, v, 0) < VECTOR(*minx)[v]) {
                MATRIX(*res, v, 0) = VECTOR(*minx)[v];
            }
            if (maxx && MATRIX(*res, v, 0) > VECTOR(*maxx)[v]) {
                MATRIX(*res, v, 0) = VECTOR(*maxx)[v];
            }
            if (miny && MATRIX(*res, v, 1) < VECTOR(*miny)[v]) {
                MATRIX(*res, v, 1) = VECTOR(*miny)[v];
            }
            if (maxy && MATRIX(*res, v, 1) > VECTOR(*maxy)[v]) {
                MATRIX(*res, v, 1) = VECTOR(*maxy)[v];
            }
        }

        temp -= difftemp;
    }

    igraph_vector_float_destroy(&dispx);
    igraph_vector_float_destroy(&dispy);
    igraph_2dgrid_destroy(&grid);
    IGRAPH_FINALLY_CLEAN(3);
    return 0;
}

/**
 * \ingroup layout
 * \function igraph_layout_fruchterman_reingold
 * \brief Places the vertices on a plane according to the Fruchterman-Reingold algorithm.
 *
 * </para><para>
 * This is a force-directed layout, see Fruchterman, T.M.J. and
 * Reingold, E.M.: Graph Drawing by Force-directed Placement.
 * Software -- Practice and Experience, 21/11, 1129--1164,
 * 1991.
 * \param graph Pointer to an initialized graph object.
 * \param res Pointer to an initialized matrix object. This will
 *        contain the result and will be resized as needed.
 * \param use_seed Logical, if true the supplied values in the
 *        \p res argument are used as an initial layout, if
 *        false a random initial layout is used.
 * \param niter The number of iterations to do. A reasonable
 *        default value is 500.
 * \param start_temp Start temperature. This is the maximum amount
 *        of movement allowed along one axis, within one step, for a
 *        vertex. Currently it is decreased linearly to zero during
 *        the iteration.
 * \param grid Whether to use the (fast but less accurate) grid based
 *        version of the algorithm. Possible values: \c
 *        IGRAPH_LAYOUT_GRID, \c IGRAPH_LAYOUT_NOGRID, \c
 *        IGRAPH_LAYOUT_AUTOGRID. The last one uses the grid based
 *        version only for large graphs, currently the ones with
 *        more than 1000 vertices.
 * \param weight Pointer to a vector containing edge weights,
 *        the attraction along the edges will be multiplied by these.
 *        It will be ignored if it is a null-pointer.
 * \param minx Pointer to a vector, or a \c NULL pointer. If not a
 *        \c NULL pointer then the vector gives the minimum
 *        \quote x \endquote coordinate for every vertex.
 * \param maxx Same as \p minx, but the maximum \quote x \endquote
 *        coordinates.
 * \param miny Pointer to a vector, or a \c NULL pointer. If not a
 *        \c NULL pointer then the vector gives the minimum
 *        \quote y \endquote coordinate for every vertex.
 * \param maxy Same as \p miny, but the maximum \quote y \endquote
 *        coordinates.
 * \return Error code.
 *
 * Time complexity: O(|V|^2) in each
 * iteration, |V| is the number of
 * vertices in the graph.
 */

int igraph_layout_fruchterman_reingold(const igraph_t *graph,
                                       igraph_matrix_t *res,
                                       igraph_bool_t use_seed,
                                       igraph_integer_t niter,
                                       igraph_real_t start_temp,
                                       igraph_layout_grid_t grid,
                                       const igraph_vector_t *weight,
                                       const igraph_vector_t *minx,
                                       const igraph_vector_t *maxx,
                                       const igraph_vector_t *miny,
                                       const igraph_vector_t *maxy) {

    igraph_integer_t no_nodes = igraph_vcount(graph);

    if (niter < 0) {
        IGRAPH_ERROR("Number of iterations must be non-negative in "
                     "Fruchterman-Reingold layout.", IGRAPH_EINVAL);
    }

    if (use_seed && (igraph_matrix_nrow(res) != no_nodes ||
                     igraph_matrix_ncol(res) != 2)) {
        IGRAPH_ERROR("Invalid start position matrix size in "
                     "Fruchterman-Reingold layout.", IGRAPH_EINVAL);
    }

    if (weight && igraph_vector_size(weight) != igraph_ecount(graph)) {
        IGRAPH_ERROR("Invalid weight vector length.", IGRAPH_EINVAL);
    }

    if (minx && igraph_vector_size(minx) != no_nodes) {
        IGRAPH_ERROR("Invalid minx vector length.", IGRAPH_EINVAL);
    }
    if (maxx && igraph_vector_size(maxx) != no_nodes) {
        IGRAPH_ERROR("Invalid maxx vector length.", IGRAPH_EINVAL);
    }
    if (minx && maxx && !igraph_vector_all_le(minx, maxx)) {
        IGRAPH_ERROR("minx must not be greater than maxx.", IGRAPH_EINVAL);
    }
    if (miny && igraph_vector_size(miny) != no_nodes) {
        IGRAPH_ERROR("Invalid miny vector length.", IGRAPH_EINVAL);
    }
    if (maxy && igraph_vector_size(maxy) != no_nodes) {
        IGRAPH_ERROR("Invalid maxy vector length.", IGRAPH_EINVAL);
    }
    if (miny && maxy && !igraph_vector_all_le(miny, maxy)) {
        IGRAPH_ERROR("miny must not be greater than maxy.", IGRAPH_EINVAL);
    }

    if (grid == IGRAPH_LAYOUT_AUTOGRID) {
        if (no_nodes > 1000) {
            grid = IGRAPH_LAYOUT_GRID;
        } else {
            grid = IGRAPH_LAYOUT_NOGRID;
        }
    }

    if (grid == IGRAPH_LAYOUT_GRID) {
        return igraph_layout_i_grid_fr(graph, res, use_seed, niter, start_temp,
                                       weight, minx, maxx, miny, maxy);
    } else {
        return igraph_layout_i_fr(graph, res, use_seed, niter, start_temp,
                                  weight, minx, maxx, miny, maxy);
    }
}

/**
 * \function igraph_layout_fruchterman_reingold_3d
 * \brief 3D Fruchterman-Reingold algorithm.
 *
 * This is the 3D version of the force based
 * Fruchterman-Reingold layout (see \ref
 * igraph_layout_fruchterman_reingold for the 2D version
 *
 * \param graph Pointer to an initialized graph object.
 * \param res Pointer to an initialized matrix object. This will
 *        contain the result and will be resized as needed.
 * \param use_seed Logical, if true the supplied values in the
 *        \p res argument are used as an initial layout, if
 *        false a random initial layout is used.
 * \param niter The number of iterations to do. A reasonable
 *        default value is 500.
 * \param start_temp Start temperature. This is the maximum amount
 *        of movement alloved along one axis, within one step, for a
 *        vertex. Currently it is decreased linearly to zero during
 *        the iteration.
 * \param weight Pointer to a vector containing edge weights,
 *        the attraction along the edges will be multiplied by these.
 *        It will be ignored if it is a null-pointer.
 * \param minx Pointer to a vector, or a \c NULL pointer. If not a
 *        \c NULL pointer then the vector gives the minimum
 *        \quote x \endquote coordinate for every vertex.
 * \param maxx Same as \p minx, but the maximum \quote x \endquote
 *        coordinates.
 * \param miny Pointer to a vector, or a \c NULL pointer. If not a
 *        \c NULL pointer then the vector gives the minimum
 *        \quote y \endquote coordinate for every vertex.
 * \param maxy Same as \p miny, but the maximum \quote y \endquote
 *        coordinates.
 * \param minz Pointer to a vector, or a \c NULL pointer. If not a
 *        \c NULL pointer then the vector gives the minimum
 *        \quote z \endquote coordinate for every vertex.
 * \param maxz Same as \p minz, but the maximum \quote z \endquote
 *        coordinates.
 * \return Error code.
 *
 * Added in version 0.2.</para><para>
 *
 * Time complexity: O(|V|^2) in each
 * iteration, |V| is the number of
 * vertices in the graph.
 *
 */

int igraph_layout_fruchterman_reingold_3d(const igraph_t *graph,
        igraph_matrix_t *res,
        igraph_bool_t use_seed,
        igraph_integer_t niter,
        igraph_real_t start_temp,
        const igraph_vector_t *weight,
        const igraph_vector_t *minx,
        const igraph_vector_t *maxx,
        const igraph_vector_t *miny,
        const igraph_vector_t *maxy,
        const igraph_vector_t *minz,
        const igraph_vector_t *maxz) {

    igraph_integer_t no_nodes = igraph_vcount(graph);
    igraph_integer_t no_edges = igraph_ecount(graph);
    igraph_integer_t i;
    igraph_vector_float_t dispx, dispy, dispz;
    igraph_real_t temp = start_temp;
    igraph_real_t difftemp = start_temp / niter;
    float width = sqrtf(no_nodes), height = width, depth = width;
    igraph_bool_t conn = 1;
    float C = 0;

    if (niter < 0) {
        IGRAPH_ERROR("Number of iterations must be non-negative in "
                     "Fruchterman-Reingold layout", IGRAPH_EINVAL);
    }

    if (use_seed && (igraph_matrix_nrow(res) != no_nodes ||
                     igraph_matrix_ncol(res) != 3)) {
        IGRAPH_ERROR("Invalid start position matrix size in "
                     "Fruchterman-Reingold layout", IGRAPH_EINVAL);
    }

    if (weight && igraph_vector_size(weight) != igraph_ecount(graph)) {
        IGRAPH_ERROR("Invalid weight vector length", IGRAPH_EINVAL);
    }

    if (minx && igraph_vector_size(minx) != no_nodes) {
        IGRAPH_ERROR("Invalid minx vector length", IGRAPH_EINVAL);
    }
    if (maxx && igraph_vector_size(maxx) != no_nodes) {
        IGRAPH_ERROR("Invalid maxx vector length", IGRAPH_EINVAL);
    }
    if (minx && maxx && !igraph_vector_all_le(minx, maxx)) {
        IGRAPH_ERROR("minx must not be greater than maxx", IGRAPH_EINVAL);
    }
    if (miny && igraph_vector_size(miny) != no_nodes) {
        IGRAPH_ERROR("Invalid miny vector length", IGRAPH_EINVAL);
    }
    if (maxy && igraph_vector_size(maxy) != no_nodes) {
        IGRAPH_ERROR("Invalid maxy vector length", IGRAPH_EINVAL);
    }
    if (miny && maxy && !igraph_vector_all_le(miny, maxy)) {
        IGRAPH_ERROR("miny must not be greater than maxy", IGRAPH_EINVAL);
    }
    if (minz && igraph_vector_size(minz) != no_nodes) {
        IGRAPH_ERROR("Invalid minz vector length", IGRAPH_EINVAL);
    }
    if (maxz && igraph_vector_size(maxz) != no_nodes) {
        IGRAPH_ERROR("Invalid maxz vector length", IGRAPH_EINVAL);
    }
    if (minz && maxz && !igraph_vector_all_le(minz, maxz)) {
        IGRAPH_ERROR("minz must not be greater than maxz", IGRAPH_EINVAL);
    }

    igraph_is_connected(graph, &conn, IGRAPH_WEAK);
    if (!conn) {
        C = no_nodes * sqrtf(no_nodes);
    }

    RNG_BEGIN();

    if (!use_seed) {
        IGRAPH_CHECK(igraph_matrix_resize(res, no_nodes, 3));
        for (i = 0; i < no_nodes; i++) {
            igraph_real_t x1 = minx ? VECTOR(*minx)[i] : -width / 2;
            igraph_real_t x2 = maxx ? VECTOR(*maxx)[i] :  width / 2;
            igraph_real_t y1 = miny ? VECTOR(*miny)[i] : -height / 2;
            igraph_real_t y2 = maxy ? VECTOR(*maxy)[i] :  height / 2;
            igraph_real_t z1 = minz ? VECTOR(*minz)[i] : -depth / 2;
            igraph_real_t z2 = maxz ? VECTOR(*maxz)[i] :  depth / 2;
            MATRIX(*res, i, 0) = RNG_UNIF(x1, x2);
            MATRIX(*res, i, 1) = RNG_UNIF(y1, y2);
            MATRIX(*res, i, 2) = RNG_UNIF(z1, z2);
        }
    }

    IGRAPH_CHECK(igraph_vector_float_init(&dispx, no_nodes));
    IGRAPH_FINALLY(igraph_vector_float_destroy, &dispx);
    IGRAPH_CHECK(igraph_vector_float_init(&dispy, no_nodes));
    IGRAPH_FINALLY(igraph_vector_float_destroy, &dispy);
    IGRAPH_CHECK(igraph_vector_float_init(&dispz, no_nodes));
    IGRAPH_FINALLY(igraph_vector_float_destroy, &dispz);

    for (i = 0; i < niter; i++) {
        igraph_integer_t v, u, e;

        /* calculate repulsive forces, we have a special version
           for unconnected graphs */
        igraph_vector_float_null(&dispx);
        igraph_vector_float_null(&dispy);
        igraph_vector_float_null(&dispz);
        if (conn) {
            for (v = 0; v < no_nodes; v++) {
                for (u = v + 1; u < no_nodes; u++) {
                    float dx = MATRIX(*res, v, 0) - MATRIX(*res, u, 0);
                    float dy = MATRIX(*res, v, 1) - MATRIX(*res, u, 1);
                    float dz = MATRIX(*res, v, 2) - MATRIX(*res, u, 2);
                    float dlen = dx * dx + dy * dy + dz * dz;

                    if (dlen == 0) {
                        dx = RNG_UNIF01() * 1e-9;
                        dy = RNG_UNIF01() * 1e-9;
                        dz = RNG_UNIF01() * 1e-9;
                        dlen = dx * dx + dy * dy + dz * dz;
                    }

                    VECTOR(dispx)[v] += dx / dlen;
                    VECTOR(dispy)[v] += dy / dlen;
                    VECTOR(dispz)[v] += dz / dlen;
                    VECTOR(dispx)[u] -= dx / dlen;
                    VECTOR(dispy)[u] -= dy / dlen;
                    VECTOR(dispz)[u] -= dz / dlen;
                }
            }
        } else {
            for (v = 0; v < no_nodes; v++) {
                for (u = v + 1; u < no_nodes; u++) {
                    float dx = MATRIX(*res, v, 0) - MATRIX(*res, u, 0);
                    float dy = MATRIX(*res, v, 1) - MATRIX(*res, u, 1);
                    float dz = MATRIX(*res, v, 2) - MATRIX(*res, u, 2);
                    float dlen, rdlen;

                    dlen = dx * dx + dy * dy + dz * dz;
                    if (dlen == 0) {
                        dx = RNG_UNIF01() * 1e-9;
                        dy = RNG_UNIF01() * 1e-9;
                        dz = RNG_UNIF01() * 1e-9;
                        dlen = dx * dx + dy * dy + dz * dz;
                    }

                    rdlen = sqrt(dlen);

                    VECTOR(dispx)[v] += dx * (C - dlen * rdlen) / (dlen * C);
                    VECTOR(dispy)[v] += dy * (C - dlen * rdlen) / (dlen * C);
                    VECTOR(dispy)[v] += dz * (C - dlen * rdlen) / (dlen * C);
                    VECTOR(dispx)[u] -= dx * (C - dlen * rdlen) / (dlen * C);
                    VECTOR(dispy)[u] -= dy * (C - dlen * rdlen) / (dlen * C);
                    VECTOR(dispz)[u] -= dz * (C - dlen * rdlen) / (dlen * C);
                }
            }
        }

        /* calculate attractive forces */
        for (e = 0; e < no_edges; e++) {
            /* each edges is an ordered pair of vertices v and u */
            igraph_integer_t v = IGRAPH_FROM(graph, e);
            igraph_integer_t u = IGRAPH_TO(graph, e);
            igraph_real_t dx = MATRIX(*res, v, 0) - MATRIX(*res, u, 0);
            igraph_real_t dy = MATRIX(*res, v, 1) - MATRIX(*res, u, 1);
            igraph_real_t dz = MATRIX(*res, v, 2) - MATRIX(*res, u, 2);
            igraph_real_t w = weight ? VECTOR(*weight)[e] : 1.0;
            igraph_real_t dlen = sqrt(dx * dx + dy * dy + dz * dz) * w;
            VECTOR(dispx)[v] -= (dx * dlen);
            VECTOR(dispy)[v] -= (dy * dlen);
            VECTOR(dispz)[v] -= (dz * dlen);
            VECTOR(dispx)[u] += (dx * dlen);
            VECTOR(dispy)[u] += (dy * dlen);
            VECTOR(dispz)[u] += (dz * dlen);
        }

        /* limit max displacement to temperature t and prevent from
           displacement outside frame */
        for (v = 0; v < no_nodes; v++) {
            igraph_real_t dx = VECTOR(dispx)[v] + RNG_UNIF01() * 1e-9;
            igraph_real_t dy = VECTOR(dispy)[v] + RNG_UNIF01() * 1e-9;
            igraph_real_t dz = VECTOR(dispz)[v] + RNG_UNIF01() * 1e-9;
            igraph_real_t displen = sqrt(dx * dx + dy * dy + dz * dz);
            igraph_real_t mx = fabs(dx) < temp ? dx : temp;
            igraph_real_t my = fabs(dy) < temp ? dy : temp;
            igraph_real_t mz = fabs(dz) < temp ? dz : temp;
            if (displen > 0) {
                MATRIX(*res, v, 0) += (dx / displen) * mx;
                MATRIX(*res, v, 1) += (dy / displen) * my;
                MATRIX(*res, v, 2) += (dz / displen) * mz;
            }
            if (minx && MATRIX(*res, v, 0) < VECTOR(*minx)[v]) {
                MATRIX(*res, v, 0) = VECTOR(*minx)[v];
            }
            if (maxx && MATRIX(*res, v, 0) > VECTOR(*maxx)[v]) {
                MATRIX(*res, v, 0) = VECTOR(*maxx)[v];
            }
            if (miny && MATRIX(*res, v, 1) < VECTOR(*miny)[v]) {
                MATRIX(*res, v, 1) = VECTOR(*miny)[v];
            }
            if (maxy && MATRIX(*res, v, 1) > VECTOR(*maxy)[v]) {
                MATRIX(*res, v, 1) = VECTOR(*maxy)[v];
            }
            if (minz && MATRIX(*res, v, 2) < VECTOR(*minz)[v]) {
                MATRIX(*res, v, 2) = VECTOR(*minz)[v];
            }
            if (maxz && MATRIX(*res, v, 2) > VECTOR(*maxz)[v]) {
                MATRIX(*res, v, 2) = VECTOR(*maxz)[v];
            }
        }

        temp -= difftemp;
    }

    RNG_END();

    igraph_vector_float_destroy(&dispx);
    igraph_vector_float_destroy(&dispy);
    igraph_vector_float_destroy(&dispz);
    IGRAPH_FINALLY_CLEAN(3);

    return 0;
}
