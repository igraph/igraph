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
#include "core/interruption.h"
#include "layout/layout_internal.h"

static igraph_error_t igraph_layout_i_fr(const igraph_t *graph,
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
    igraph_vector_t dispx, dispy;
    igraph_real_t temp = start_temp;
    igraph_real_t difftemp = start_temp / niter;
    igraph_bool_t conn = true;
    igraph_real_t C = 0;

    IGRAPH_CHECK(igraph_is_connected(graph, &conn, IGRAPH_WEAK));
    if (!conn) {
        C = no_nodes * sqrt(no_nodes);
    }

    if (!use_seed) {
        igraph_i_layout_random_bounded(graph, res, minx, maxx, miny, maxy);
    }

    IGRAPH_VECTOR_INIT_FINALLY(&dispx, no_nodes);
    IGRAPH_VECTOR_INIT_FINALLY(&dispy, no_nodes);

    RNG_BEGIN();
    for (i = 0; i < niter; i++) {
        igraph_integer_t v, u, e;

        IGRAPH_ALLOW_INTERRUPTION();

        /* calculate repulsive forces, we have a special version
           for unconnected graphs */
        igraph_vector_null(&dispx);
        igraph_vector_null(&dispy);
        if (conn) {
            for (v = 0; v < no_nodes; v++) {
                for (u = v + 1; u < no_nodes; u++) {
                    igraph_real_t dx = MATRIX(*res, v, 0) - MATRIX(*res, u, 0);
                    igraph_real_t dy = MATRIX(*res, v, 1) - MATRIX(*res, u, 1);
                    igraph_real_t dlen = dx * dx + dy * dy;

                    while (dlen == 0) {
                        dx = RNG_UNIF(-1e-9, 1e-9);
                        dy = RNG_UNIF(-1e-9, 1e-9);
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
                    igraph_real_t dx = MATRIX(*res, v, 0) - MATRIX(*res, u, 0);
                    igraph_real_t dy = MATRIX(*res, v, 1) - MATRIX(*res, u, 1);
                    igraph_real_t dlen, rdlen;

                    dlen = dx * dx + dy * dy;
                    while (dlen == 0) {
                        dx = RNG_UNIF(-1e-9, 1e-9);
                        dy = RNG_UNIF(-1e-9, 1e-9);
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
            /* each edge is an ordered pair of vertices v and u */
            igraph_integer_t v = IGRAPH_FROM(graph, e);
            igraph_integer_t u = IGRAPH_TO(graph, e);
            igraph_real_t dx = MATRIX(*res, v, 0) - MATRIX(*res, u, 0);
            igraph_real_t dy = MATRIX(*res, v, 1) - MATRIX(*res, u, 1);
            igraph_real_t w = weight ? VECTOR(*weight)[e] : 1.0;
            igraph_real_t dlen = sqrt(dx*dx + dy*dy) * w;
            VECTOR(dispx)[v] -= (dx * dlen);
            VECTOR(dispy)[v] -= (dy * dlen);
            VECTOR(dispx)[u] += (dx * dlen);
            VECTOR(dispy)[u] += (dy * dlen);
        }

        /* limit max displacement to temperature t and prevent from
           displacement outside frame */
        for (v = 0; v < no_nodes; v++) {
            igraph_real_t dx = VECTOR(dispx)[v] + RNG_UNIF(-1e-9, 1e-9);
            igraph_real_t dy = VECTOR(dispy)[v] + RNG_UNIF(-1e-9, 1e-9);
            igraph_real_t displen = sqrt(dx * dx + dy * dy);

            if (displen > temp) {
                dx *= temp/displen;
                dy *= temp/displen;
            }

            if (displen > 0) {
                MATRIX(*res, v, 0) += dx;
                MATRIX(*res, v, 1) += dy;
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

    igraph_vector_destroy(&dispx);
    igraph_vector_destroy(&dispy);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_layout_i_grid_fr(
        const igraph_t *graph,
        igraph_matrix_t *res, igraph_bool_t use_seed,
        igraph_integer_t niter, igraph_real_t start_temp,
        const igraph_vector_t *weight, const igraph_vector_t *minx,
        const igraph_vector_t *maxx, const igraph_vector_t *miny,
        const igraph_vector_t *maxy) {

    igraph_integer_t no_nodes = igraph_vcount(graph);
    igraph_integer_t no_edges = igraph_ecount(graph);
    igraph_real_t width = sqrt(no_nodes), height = width;
    igraph_2dgrid_t grid;
    igraph_vector_t dispx, dispy;
    igraph_real_t temp = start_temp;
    igraph_real_t difftemp = start_temp / niter;
    igraph_2dgrid_iterator_t vidit;
    igraph_integer_t i;
    const igraph_real_t cellsize = 2.0;

    if (!use_seed) {
        igraph_i_layout_random_bounded(graph, res, minx, maxx, miny, maxy);
    }

    /* make grid */
    IGRAPH_CHECK(igraph_2dgrid_init(&grid, res, -width / 2, width / 2, cellsize,
                                    -height / 2, height / 2, cellsize));
    IGRAPH_FINALLY(igraph_2dgrid_destroy, &grid);

    /* place vertices on grid */
    for (i = 0; i < no_nodes; i++) {
        igraph_2dgrid_add2(&grid, i);
    }

    IGRAPH_VECTOR_INIT_FINALLY(&dispx, no_nodes);
    IGRAPH_VECTOR_INIT_FINALLY(&dispy, no_nodes);

    RNG_BEGIN();
    for (i = 0; i < niter; i++) {
        igraph_integer_t v, u, e;

        IGRAPH_ALLOW_INTERRUPTION();

        igraph_vector_null(&dispx);
        igraph_vector_null(&dispy);

        /* repulsion */
        igraph_2dgrid_reset(&grid, &vidit);
        while ( (v = igraph_2dgrid_next(&grid, &vidit) - 1) != -1) {
            while ( (u = igraph_2dgrid_next_nei(&grid, &vidit) - 1) != -1) {
                igraph_real_t dx = MATRIX(*res, v, 0) - MATRIX(*res, u, 0);
                igraph_real_t dy = MATRIX(*res, v, 1) - MATRIX(*res, u, 1);
                igraph_real_t dlen = dx * dx + dy * dy;
                while (dlen == 0) {
                    dx = RNG_UNIF(-1e-9, 1e-9);
                    dy = RNG_UNIF(-1e-9, 1e-9);
                    dlen = dx * dx + dy * dy;
                }
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
            igraph_real_t dlen = sqrt(dx*dx + dy*dy) * w;
            VECTOR(dispx)[v] -= (dx * dlen);
            VECTOR(dispy)[v] -= (dy * dlen);
            VECTOR(dispx)[u] += (dx * dlen);
            VECTOR(dispy)[u] += (dy * dlen);
        }

        /* update */
        for (v = 0; v < no_nodes; v++) {
            igraph_real_t dx = VECTOR(dispx)[v] + RNG_UNIF(-1e-9, 1e-9);
            igraph_real_t dy = VECTOR(dispy)[v] + RNG_UNIF(-1e-9, 1e-9);
            igraph_real_t displen = sqrt(dx * dx + dy * dy);

            if (displen > temp) {
                dx *= temp/displen;
                dy *= temp/displen;
            }

            if (displen > 0) {
                MATRIX(*res, v, 0) += dx;
                MATRIX(*res, v, 1) += dy;
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

    igraph_vector_destroy(&dispx);
    igraph_vector_destroy(&dispy);
    igraph_2dgrid_destroy(&grid);
    IGRAPH_FINALLY_CLEAN(3);
    return IGRAPH_SUCCESS;
}

/**
 * \ingroup layout
 * \function igraph_layout_fruchterman_reingold
 * \brief Places the vertices on a plane according to the Fruchterman-Reingold algorithm.
 *
 * </para><para>
 * This is a force-directed layout that simulates an attractive force \c f_a between
 * connected vertex pairs and a repulsive force \c f_r between all vertex pairs.
 * The forces are computed as a function of the distance \c d between the two vertices as
 *
 * </para><para>
 * <code>f_a(d) = -w * d^2</code> and <code>f_r(d) = 1/d</code>,
 *
 * </para><para>
 * where \c w represents the edge weight. The equilibrium distance of two connected
 * vertices is thus <code>1/w^3</code>, assuming no other forces acting on them.
 *
 * </para><para>
 * In disconnected graphs, igraph effectively inserts a weak connection of weight
 * <code>n^(-3/2)</code> between all pairs of vertices, where \c n is the vertex count.
 * This ensures that components are kept near each other.
 *
 * </para><para>
 * Reference:
 *
 * </para><para>
 * Fruchterman, T.M.J. and Reingold, E.M.:
 * Graph Drawing by Force-directed Placement.
 * Software -- Practice and Experience, 21/11, 1129--1164,
 * 1991. https://doi.org/10.1002/spe.4380211102
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
 *        of movement allowed along one axis, within one step, for a
 *        vertex. Currently it is decreased linearly to zero during
 *        the iteration.
 * \param grid Whether to use the (fast but less accurate) grid based
 *        version of the algorithm. Possible values: \c
 *        IGRAPH_LAYOUT_GRID, \c IGRAPH_LAYOUT_NOGRID, \c
 *        IGRAPH_LAYOUT_AUTOGRID. The last one uses the grid based
 *        version only for large graphs, currently the ones with
 *        more than 1000 vertices.
 * \param weights Pointer to a vector containing edge weights,
 *        the attraction along the edges will be multiplied by these.
 *        Weights must be positive.
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

igraph_error_t igraph_layout_fruchterman_reingold(const igraph_t *graph,
                                       igraph_matrix_t *res,
                                       igraph_bool_t use_seed,
                                       igraph_integer_t niter,
                                       igraph_real_t start_temp,
                                       igraph_layout_grid_t grid,
                                       const igraph_vector_t *weights,
                                       const igraph_vector_t *minx,
                                       const igraph_vector_t *maxx,
                                       const igraph_vector_t *miny,
                                       const igraph_vector_t *maxy) {

    igraph_integer_t no_nodes = igraph_vcount(graph);
    igraph_integer_t no_edges = igraph_ecount(graph);

    if (niter < 0) {
        IGRAPH_ERROR("Number of iterations must be non-negative in "
                     "Fruchterman-Reingold layout.", IGRAPH_EINVAL);
    }

    if (use_seed && (igraph_matrix_nrow(res) != no_nodes ||
                     igraph_matrix_ncol(res) != 2)) {
        IGRAPH_ERROR("Invalid start position matrix size in "
                     "Fruchterman-Reingold layout.", IGRAPH_EINVAL);
    }

    if (weights && igraph_vector_size(weights) != no_edges) {
        IGRAPH_ERROR("Invalid weight vector length.", IGRAPH_EINVAL);
    }
    if (weights && no_edges > 0 && igraph_vector_min(weights) <= 0) {
        IGRAPH_ERROR("Weights must be positive for Fruchterman-Reingold layout.", IGRAPH_EINVAL);
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
                                       weights, minx, maxx, miny, maxy);
    } else {
        return igraph_layout_i_fr(graph, res, use_seed, niter, start_temp,
                                  weights, minx, maxx, miny, maxy);
    }
}

/**
 * \function igraph_layout_fruchterman_reingold_3d
 * \brief 3D Fruchterman-Reingold algorithm.
 *
 * This is the 3D version of the force based Fruchterman-Reingold layout.
 * See \ref igraph_layout_fruchterman_reingold() for the 2D version.
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
 * \param weights Pointer to a vector containing edge weights,
 *        the attraction along the edges will be multiplied by these.
 *        Weights must be positive.
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

igraph_error_t igraph_layout_fruchterman_reingold_3d(const igraph_t *graph,
        igraph_matrix_t *res,
        igraph_bool_t use_seed,
        igraph_integer_t niter,
        igraph_real_t start_temp,
        const igraph_vector_t *weights,
        const igraph_vector_t *minx,
        const igraph_vector_t *maxx,
        const igraph_vector_t *miny,
        const igraph_vector_t *maxy,
        const igraph_vector_t *minz,
        const igraph_vector_t *maxz) {

    const igraph_integer_t no_nodes = igraph_vcount(graph);
    const igraph_integer_t no_edges = igraph_ecount(graph);
    igraph_integer_t i;
    igraph_vector_t dispx, dispy, dispz;
    igraph_real_t temp = start_temp;
    igraph_real_t difftemp = start_temp / niter;
    igraph_bool_t conn = true;
    igraph_real_t C = 0;

    if (niter < 0) {
        IGRAPH_ERROR("Number of iterations must be non-negative in "
                     "Fruchterman-Reingold layout", IGRAPH_EINVAL);
    }

    if (use_seed && (igraph_matrix_nrow(res) != no_nodes ||
                     igraph_matrix_ncol(res) != 3)) {
        IGRAPH_ERROR("Invalid start position matrix size in "
                     "Fruchterman-Reingold layout", IGRAPH_EINVAL);
    }

    if (weights && igraph_vector_size(weights) != igraph_ecount(graph)) {
        IGRAPH_ERROR("Invalid weight vector length", IGRAPH_EINVAL);
    }
    if (weights && no_edges > 0 && igraph_vector_min(weights) <= 0) {
        IGRAPH_ERROR("Weights must be positive for Fruchterman-Reingold layout.", IGRAPH_EINVAL);
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

    IGRAPH_CHECK(igraph_is_connected(graph, &conn, IGRAPH_WEAK));
    if (!conn) {
        C = no_nodes * sqrt(no_nodes);
    }

    if (!use_seed) {
        igraph_i_layout_random_bounded_3d(graph, res, minx, maxx, miny, maxy, minz, maxz);
    }

    IGRAPH_VECTOR_INIT_FINALLY(&dispx, no_nodes);
    IGRAPH_VECTOR_INIT_FINALLY(&dispy, no_nodes);
    IGRAPH_VECTOR_INIT_FINALLY(&dispz, no_nodes);

    RNG_BEGIN();
    for (i = 0; i < niter; i++) {
        igraph_integer_t v, u, e;

        IGRAPH_ALLOW_INTERRUPTION();

        /* calculate repulsive forces, we have a special version
           for unconnected graphs */
        igraph_vector_null(&dispx);
        igraph_vector_null(&dispy);
        igraph_vector_null(&dispz);
        if (conn) {
            for (v = 0; v < no_nodes; v++) {
                for (u = v + 1; u < no_nodes; u++) {
                    igraph_real_t dx = MATRIX(*res, v, 0) - MATRIX(*res, u, 0);
                    igraph_real_t dy = MATRIX(*res, v, 1) - MATRIX(*res, u, 1);
                    igraph_real_t dz = MATRIX(*res, v, 2) - MATRIX(*res, u, 2);
                    igraph_real_t dlen = dx * dx + dy * dy + dz * dz;

                    while (dlen == 0) {
                        dx = RNG_UNIF(-1e-9, 1e-9);
                        dy = RNG_UNIF(-1e-9, 1e-9);
                        dz = RNG_UNIF(-1e-9, 1e-9);
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
                    igraph_real_t dx = MATRIX(*res, v, 0) - MATRIX(*res, u, 0);
                    igraph_real_t dy = MATRIX(*res, v, 1) - MATRIX(*res, u, 1);
                    igraph_real_t dz = MATRIX(*res, v, 2) - MATRIX(*res, u, 2);
                    igraph_real_t dlen, rdlen;

                    dlen = dx * dx + dy * dy + dz * dz;
                    while (dlen == 0) {
                        dx = RNG_UNIF(-1e-9, 1e-9);
                        dy = RNG_UNIF(-1e-9, 1e-9);
                        dz = RNG_UNIF(-1e-9, 1e-9);
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
            igraph_real_t w = weights ? VECTOR(*weights)[e] : 1.0;
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
            igraph_real_t dx = VECTOR(dispx)[v] + RNG_UNIF(-1e-9, 1e-9);
            igraph_real_t dy = VECTOR(dispy)[v] + RNG_UNIF(-1e-9, 1e-9);
            igraph_real_t dz = VECTOR(dispz)[v] + RNG_UNIF(-1e-9, 1e-9);
            igraph_real_t displen = sqrt(dx * dx + dy * dy + dz * dz);

            if (displen > temp) {
                dx *= temp/displen;
                dy *= temp/displen;
                dz *= temp/displen;
            }

            if (displen > 0) {
                MATRIX(*res, v, 0) += dx;
                MATRIX(*res, v, 1) += dy;
                MATRIX(*res, v, 2) += dz;
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

    igraph_vector_destroy(&dispx);
    igraph_vector_destroy(&dispy);
    igraph_vector_destroy(&dispz);
    IGRAPH_FINALLY_CLEAN(3);

    return IGRAPH_SUCCESS;
}
