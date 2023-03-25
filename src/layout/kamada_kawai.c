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
#include "igraph_paths.h"

#include "core/interruption.h"
#include "layout/layout_internal.h"

/* Energy gradient values below this threshold are considered to be zero. */
#define KK_EPS 1e-13

/**
 * \ingroup layout
 * \function igraph_layout_kamada_kawai
 * \brief Places the vertices on a plane according to the Kamada-Kawai algorithm.
 *
 * This is a force-directed layout. A spring is inserted between all pairs
 * of vertices, both those which are directly connected and those that are not.
 * The unstretched length of springs is chosen based on the undirected graph distance
 * between the corresponding pair of vertices. Thus, in a weighted graph, increasing
 * the weight between two vertices pushes them apart. The Young modulus of springs
 * is inversely proportional to the graph distance, ensuring that springs between
 * far-apart veritces will have a smaller effect on the layout.
 *
 * </para><para>
 * Disconnected graphs are handled by assuming that the graph distance between
 * vertices in different components is the same as the graph diameter.
 *
 * </para><para>
 * This layout works particularly well for locally connected spatial networks
 * such as lattices.
 *
 * </para><para>
 * This layout algorithm is not suitable for large graphs. The memory
 * requirements are of the order O(|V|^2).
 *
 * </para><para>
 * Reference:
 *
 * </para><para>
 * Kamada, T. and Kawai, S.:
 * An Algorithm for Drawing General Undirected Graphs.
 * Information Processing Letters, 31/1, 7--15, 1989.
 * https://doi.org/10.1016/0020-0190(89)90102-6
 *
 * \param graph A graph object.
 * \param res Pointer to an initialized matrix object. This will
 *        contain the result (x-positions in column zero and
 *        y-positions in column one) and will be resized if needed.
 * \param use_seed Boolean, whether to use the values supplied in the
 *        \p res argument as the initial configuration. If zero and there
 *        are any limits on the X or Y coordinates, then a random initial
 *        configuration is used. Otherwise the vertices are placed on a
 *        circle of radius 1 as the initial configuration.
 * \param maxiter The maximum number of iterations to perform. A reasonable
 *        default value is at least ten (or more) times the number of
 *        vertices.
 * \param epsilon Stop the iteration, if the maximum delta value of the
 *        algorithm is smaller than still. It is safe to leave it at zero,
 *        and then \p maxiter iterations are performed.
 * \param kkconst The Kamada-Kawai vertex attraction constant.
 *        Typical value: number of vertices.
 * \param weights Edge weights, larger values will result longer edges.
 *        Weights must be positive. Pass \c NULL to assume unit weights
 *        for all edges.
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
 * Time complexity: O(|V|) for each iteration, after an O(|V|^2
 * log|V|) initialization step. |V| is the number of vertices in the
 * graph.
 */

igraph_error_t igraph_layout_kamada_kawai(const igraph_t *graph, igraph_matrix_t *res,
                               igraph_bool_t use_seed, igraph_integer_t maxiter,
                               igraph_real_t epsilon, igraph_real_t kkconst,
                               const igraph_vector_t *weights,
                               const igraph_vector_t *minx, const igraph_vector_t *maxx,
                               const igraph_vector_t *miny, const igraph_vector_t *maxy) {

    igraph_integer_t no_nodes = igraph_vcount(graph);
    igraph_integer_t no_edges = igraph_ecount(graph);
    igraph_real_t L, L0 = sqrt(no_nodes);
    igraph_matrix_t dij, lij, kij;
    igraph_real_t max_dij;
    igraph_vector_t D1, D2;
    igraph_integer_t i, j, m;

    if (maxiter < 0) {
        IGRAPH_ERROR("Number of iterations must be non-negative in "
                     "Kamada-Kawai layout.", IGRAPH_EINVAL);
    }
    if (kkconst <= 0) {
        IGRAPH_ERROR("`K' constant must be positive in Kamada-Kawai layout.",
                     IGRAPH_EINVAL);
    }

    if (use_seed && (igraph_matrix_nrow(res) != no_nodes ||
                     igraph_matrix_ncol(res) != 2)) {
        IGRAPH_ERROR("Invalid start position matrix size in "
                     "Kamada-Kawai layout.", IGRAPH_EINVAL);
    }
    if (weights && igraph_vector_size(weights) != no_edges) {
        IGRAPH_ERROR("Invalid weight vector length.", IGRAPH_EINVAL);
    }
    if (weights && no_edges > 0 && igraph_vector_min(weights) <= 0) {
        IGRAPH_ERROR("Weights must be positive for Kamada-Kawai layout.", IGRAPH_EINVAL);
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

    if (!use_seed) {
        if (minx || maxx || miny || maxy) {
            igraph_i_layout_random_bounded(graph, res, minx, maxx, miny, maxy);
        } else {
            igraph_layout_circle(graph, res, /* order= */ igraph_vss_all());
            /* The original paper recommends using a radius of 0.5*L0 here.
             * The coefficient of 0.36 was chosen empirically so that this initial
             * layout would be as close as possible to the equilibrium layout
             * when the graph is a cycle graph. */
            igraph_matrix_scale(res, 0.36 * L0);
        }
    }

    if (no_nodes <= 1) {
        return IGRAPH_SUCCESS;
    }

    IGRAPH_MATRIX_INIT_FINALLY(&dij, no_nodes, no_nodes);
    IGRAPH_MATRIX_INIT_FINALLY(&kij, no_nodes, no_nodes);
    IGRAPH_MATRIX_INIT_FINALLY(&lij, no_nodes, no_nodes);

    IGRAPH_CHECK(igraph_distances_dijkstra(graph, &dij, igraph_vss_all(),
                 igraph_vss_all(), weights, IGRAPH_ALL));

    /* Find largest finite distance */
    max_dij = 0.0;
    for (i = 0; i < no_nodes; i++) {
        for (j = i + 1; j < no_nodes; j++) {
            if (!isfinite(MATRIX(dij, i, j))) {
                continue;
            }
            if (MATRIX(dij, i, j) > max_dij) {
                max_dij = MATRIX(dij, i, j);
            }
        }
    }

    /* Replace infinite distances by the largest finite distance,
     * effectively making the graph connected. */
    for (i = 0; i < no_nodes; i++) {
        for (j = 0; j < no_nodes; j++) {
            if (MATRIX(dij, i, j) > max_dij) {
                MATRIX(dij, i, j) = max_dij;
            }
        }
    }

    L = L0 / max_dij;
    for (i = 0; i < no_nodes; i++) {
        for (j = 0; j < no_nodes; j++) {
            igraph_real_t tmp = MATRIX(dij, i, j) * MATRIX(dij, i, j);
            if (i == j) {
                continue;
            }
            MATRIX(kij, i, j) = kkconst / tmp;
            MATRIX(lij, i, j) = L * MATRIX(dij, i, j);
        }
    }

    /* Initialize delta */
    IGRAPH_VECTOR_INIT_FINALLY(&D1, no_nodes);
    IGRAPH_VECTOR_INIT_FINALLY(&D2, no_nodes);
    for (m = 0; m < no_nodes; m++) {
        igraph_real_t myD1 = 0.0, myD2 = 0.0;
        for (i = 0; i < no_nodes; i++) {
            igraph_real_t dx, dy, mi_dist;
            if (i == m) {
                continue;
            }
            dx = MATRIX(*res, m, 0) - MATRIX(*res, i, 0);
            dy = MATRIX(*res, m, 1) - MATRIX(*res, i, 1);
            mi_dist = sqrt(dx*dx + dy*dy);
            myD1 += MATRIX(kij, m, i) * (dx - MATRIX(lij, m, i) * dx / mi_dist);
            myD2 += MATRIX(kij, m, i) * (dy - MATRIX(lij, m, i) * dy / mi_dist);
        }
        VECTOR(D1)[m] = myD1;
        VECTOR(D2)[m] = myD2;
    }

    for (j = 0; j < maxiter; j++) {
        igraph_real_t myD1, myD2, A, B, C;
        igraph_real_t max_delta, delta_x, delta_y;
        igraph_real_t old_x, old_y, new_x, new_y;

        IGRAPH_ALLOW_INTERRUPTION();

        myD1 = 0.0, myD2 = 0.0, A = 0.0, B = 0.0, C = 0.0;

        /* Select maximal delta */
        m = 0; max_delta = -1;
        for (i = 0; i < no_nodes; i++) {
            igraph_real_t delta = (VECTOR(D1)[i] * VECTOR(D1)[i] +
                                   VECTOR(D2)[i] * VECTOR(D2)[i]);
            if (delta > max_delta) {
                m = i; max_delta = delta;
            }
        }
        if (max_delta < epsilon) {
            break;
        }
        old_x = MATRIX(*res, m, 0);
        old_y = MATRIX(*res, m, 1);

        /* Calculate D1 and D2, A, B, C */
        for (i = 0; i < no_nodes; i++) {
            igraph_real_t dx, dy, dist, den;
            if (i == m) {
                continue;
            }
            dx = old_x - MATRIX(*res, i, 0);
            dy = old_y - MATRIX(*res, i, 1);
            dist = sqrt(dx*dx + dy*dy);
            den = dist * (dx * dx + dy * dy);
            A += MATRIX(kij, m, i) * (1 - MATRIX(lij, m, i) * dy * dy / den);
            B += MATRIX(kij, m, i) * MATRIX(lij, m, i) * dx * dy / den;
            C += MATRIX(kij, m, i) * (1 - MATRIX(lij, m, i) * dx * dx / den);
        }
        myD1 = VECTOR(D1)[m];
        myD2 = VECTOR(D2)[m];

        /* We need to solve the following linear equations, corresponding to
         * eqs. (11) and (12) in the paper.
         *
         * A * delta_x + B * delta_y == myD1
         * B * delta_x + C * delta_y == myD2
         *
         * We special-case the equilibrium case, i.e. when the energy gradient
         * is zero and no displacement is necessary. This is important for the
         * case of path graphs, where the determinant of the LHS will be
         * zero in equilibrium, causing numerical problems.
         */
        if (myD1*myD1 + myD2*myD2 < KK_EPS*KK_EPS) {
            delta_x = 0;
            delta_y = 0;
        } else {
            igraph_real_t det = C * A - B * B;
            delta_y = (B * myD1 - A * myD2) / det;
            delta_x = (B * myD2 - C * myD1) / det;
        }

        new_x = old_x + delta_x;
        new_y = old_y + delta_y;

        /* Limits, if given */
        if (minx && new_x < VECTOR(*minx)[m]) {
            new_x = VECTOR(*minx)[m];
        }
        if (maxx && new_x > VECTOR(*maxx)[m]) {
            new_x = VECTOR(*maxx)[m];
        }
        if (miny && new_y < VECTOR(*miny)[m]) {
            new_y = VECTOR(*miny)[m];
        }
        if (maxy && new_y > VECTOR(*maxy)[m]) {
            new_y = VECTOR(*maxy)[m];
        }

        /* Update delta, only with/for the affected node */
        VECTOR(D1)[m] = VECTOR(D2)[m] = 0.0;
        for (i = 0; i < no_nodes; i++) {
            igraph_real_t old_dx, old_dy, new_dx, new_dy, new_mi_dist, old_mi_dist;
            if (i == m) {
                continue;
            }
            old_dx = old_x - MATRIX(*res, i, 0);
            old_dy = old_y - MATRIX(*res, i, 1);
            old_mi_dist = sqrt(old_dx*old_dx + old_dy*old_dy);
            new_dx = new_x - MATRIX(*res, i, 0);
            new_dy = new_y - MATRIX(*res, i, 1);
            new_mi_dist = sqrt(new_dx*new_dx + new_dy*new_dy);

            VECTOR(D1)[i] -= MATRIX(kij, m, i) *
                             (-old_dx + MATRIX(lij, m, i) * old_dx / old_mi_dist);
            VECTOR(D2)[i] -= MATRIX(kij, m, i) *
                             (-old_dy + MATRIX(lij, m, i) * old_dy / old_mi_dist);
            VECTOR(D1)[i] += MATRIX(kij, m, i) *
                             (-new_dx + MATRIX(lij, m, i) * new_dx / new_mi_dist);
            VECTOR(D2)[i] += MATRIX(kij, m, i) *
                             (-new_dy + MATRIX(lij, m, i) * new_dy / new_mi_dist);

            VECTOR(D1)[m] += MATRIX(kij, m, i) *
                             (new_dx - MATRIX(lij, m, i) * new_dx / new_mi_dist);
            VECTOR(D2)[m] += MATRIX(kij, m, i) *
                             (new_dy - MATRIX(lij, m, i) * new_dy / new_mi_dist);
        }

        /* Update coordinates*/
        MATRIX(*res, m, 0) = new_x;
        MATRIX(*res, m, 1) = new_y;
    }

    igraph_vector_destroy(&D2);
    igraph_vector_destroy(&D1);
    igraph_matrix_destroy(&lij);
    igraph_matrix_destroy(&kij);
    igraph_matrix_destroy(&dij);
    IGRAPH_FINALLY_CLEAN(5);

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup layout
 * \function igraph_layout_kamada_kawai_3d
 * \brief 3D version of the Kamada-Kawai layout generator.
 *
 * This is the 3D version of \ref igraph_layout_kamada_kawai().
 * See the documentation of that function for more information.
 *
 * </para><para>
 * This layout algorithm is not suitable for large graphs. The memory
 * requirements are of the order O(|V|^2).
 *
 * \param graph A graph object.
 * \param res Pointer to an initialized matrix object. This will
 *        contain the result (x-, y- and z-positions in columns one
 *        through three) and will be resized if needed.
 * \param use_seed Boolean, whether to use the values supplied in the
 *        \p res argument as the initial configuration. If zero and there
 *        are any limits on the z, y or z coordinates, then a random initial
 *        configuration is used. Otherwise the vertices are placed uniformly
 *        on a sphere of radius 1 as the initial configuration.
 * \param maxiter The maximum number of iterations to perform. A reasonable
 *        default value is at least ten (or more) times the number of
 *        vertices.
 * \param epsilon Stop the iteration, if the maximum delta value of the
 *        algorithm is smaller than still. It is safe to leave it at zero,
 *        and then \p maxiter iterations are performed.
 * \param kkconst The Kamada-Kawai vertex attraction constant.
 *        Typical value: number of vertices.
 * \param weights Edge weights, larger values will result longer edges.
 *        Weights must be positive. Pass \c NULL to assume unit weights
 *        for all edges.
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
 * Time complexity: O(|V|) for each iteration, after an O(|V|^2
 * log|V|) initialization step. |V| is the number of vertices in the
 * graph.
 */

igraph_error_t igraph_layout_kamada_kawai_3d(const igraph_t *graph, igraph_matrix_t *res,
                                  igraph_bool_t use_seed, igraph_integer_t maxiter,
                                  igraph_real_t epsilon, igraph_real_t kkconst,
                                  const igraph_vector_t *weights,
                                  const igraph_vector_t *minx, const igraph_vector_t *maxx,
                                  const igraph_vector_t *miny, const igraph_vector_t *maxy,
                                  const igraph_vector_t *minz, const igraph_vector_t *maxz) {

    const igraph_integer_t no_nodes = igraph_vcount(graph);
    const igraph_integer_t no_edges = igraph_ecount(graph);
    igraph_real_t L, L0 = sqrt(no_nodes);
    igraph_matrix_t dij, lij, kij;
    igraph_real_t max_dij;
    igraph_vector_t D1, D2, D3;
    igraph_integer_t i, j, m;

    if (maxiter < 0) {
        IGRAPH_ERROR("Number of iterations must be non-negatice in "
                     "Kamada-Kawai layout", IGRAPH_EINVAL);
    }
    if (kkconst <= 0) {
        IGRAPH_ERROR("`K' constant must be positive in Kamada-Kawai layout",
                     IGRAPH_EINVAL);
    }

    if (use_seed && (igraph_matrix_nrow(res) != no_nodes ||
                     igraph_matrix_ncol(res) != 3)) {
        IGRAPH_ERROR("Invalid start position matrix size in "
                     "3d Kamada-Kawai layout", IGRAPH_EINVAL);
    }
    if (weights && igraph_vector_size(weights) != no_edges) {
        IGRAPH_ERROR("Invalid weight vector length", IGRAPH_EINVAL);
    }
    if (weights && no_edges > 0 && igraph_vector_min(weights) <= 0) {
        IGRAPH_ERROR("Weights must be positive for Kamada-Kawai layout.", IGRAPH_EINVAL);
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

    if (!use_seed) {
        if (minx || maxx || miny || maxy || minz || maxz) {
            igraph_i_layout_random_bounded_3d(graph, res, minx, maxx, miny, maxy, minz, maxz);
        } else {
            igraph_layout_sphere(graph, res);
            /* The coefficient of 0.36 was chosen empirically so that this initial layout
             * would be as close as possible to the equilibrium layout when the graph is
             * a Goldberg polyhedron, i.e. having a naturally spherical layout. */
            igraph_matrix_scale(res, 0.36*L0);
        }
    }

    if (no_nodes <= 1) {
        return IGRAPH_SUCCESS;
    }

    IGRAPH_MATRIX_INIT_FINALLY(&dij, no_nodes, no_nodes);
    IGRAPH_MATRIX_INIT_FINALLY(&kij, no_nodes, no_nodes);
    IGRAPH_MATRIX_INIT_FINALLY(&lij, no_nodes, no_nodes);
    IGRAPH_CHECK(igraph_distances_dijkstra(graph, &dij, igraph_vss_all(),
                 igraph_vss_all(), weights, IGRAPH_ALL));

    max_dij = 0.0;
    for (i = 0; i < no_nodes; i++) {
        for (j = i + 1; j < no_nodes; j++) {
            if (!isfinite(MATRIX(dij, i, j))) {
                continue;
            }
            if (MATRIX(dij, i, j) > max_dij) {
                max_dij = MATRIX(dij, i, j);
            }
        }
    }
    for (i = 0; i < no_nodes; i++) {
        for (j = 0; j < no_nodes; j++) {
            if (MATRIX(dij, i, j) > max_dij) {
                MATRIX(dij, i, j) = max_dij;
            }
        }
    }

    L = L0 / max_dij;
    for (i = 0; i < no_nodes; i++) {
        for (j = 0; j < no_nodes; j++) {
            igraph_real_t tmp = MATRIX(dij, i, j) * MATRIX(dij, i, j);
            if (i == j) {
                continue;
            }
            MATRIX(kij, i, j) = kkconst / tmp;
            MATRIX(lij, i, j) = L * MATRIX(dij, i, j);
        }
    }

    /* Initialize delta */
    IGRAPH_VECTOR_INIT_FINALLY(&D1, no_nodes);
    IGRAPH_VECTOR_INIT_FINALLY(&D2, no_nodes);
    IGRAPH_VECTOR_INIT_FINALLY(&D3, no_nodes);
    for (m = 0; m < no_nodes; m++) {
        igraph_real_t dx, dy, dz, mi_dist;
        igraph_real_t myD1 = 0.0, myD2 = 0.0, myD3 = 0.0;
        for (i = 0; i < no_nodes; i++) {
            if (i == m) {
                continue;
            }
            dx = MATRIX(*res, m, 0) - MATRIX(*res, i, 0);
            dy = MATRIX(*res, m, 1) - MATRIX(*res, i, 1);
            dz = MATRIX(*res, m, 2) - MATRIX(*res, i, 2);
            mi_dist = sqrt(dx * dx + dy * dy + dz * dz);
            myD1 += MATRIX(kij, m, i) * (dx - MATRIX(lij, m, i) * dx / mi_dist);
            myD2 += MATRIX(kij, m, i) * (dy - MATRIX(lij, m, i) * dy / mi_dist);
            myD3 += MATRIX(kij, m, i) * (dz - MATRIX(lij, m, i) * dz / mi_dist);
        }
        VECTOR(D1)[m] = myD1;
        VECTOR(D2)[m] = myD2;
        VECTOR(D3)[m] = myD3;
    }

    for (j = 0; j < maxiter; j++) {

        igraph_real_t Ax = 0.0, Ay = 0.0, Az = 0.0;
        igraph_real_t Axx = 0.0, Axy = 0.0, Axz = 0.0, Ayy = 0.0, Ayz = 0.0, Azz = 0.0;
        igraph_real_t max_delta, delta_x, delta_y, delta_z;
        igraph_real_t old_x, old_y, old_z, new_x, new_y, new_z;

        IGRAPH_ALLOW_INTERRUPTION();

        /* Select maximal delta */
        m = 0; max_delta = -1;
        for (i = 0; i < no_nodes; i++) {
            igraph_real_t delta = (VECTOR(D1)[i] * VECTOR(D1)[i] +
                                   VECTOR(D2)[i] * VECTOR(D2)[i] +
                                   VECTOR(D3)[i] * VECTOR(D3)[i]);
            if (delta > max_delta) {
                m = i; max_delta = delta;
            }
        }
        if (max_delta < epsilon) {
            break;
        }
        old_x = MATRIX(*res, m, 0);
        old_y = MATRIX(*res, m, 1);
        old_z = MATRIX(*res, m, 2);

        /* Calculate D1, D2 and D3, and other coefficients */
        for (i = 0; i < no_nodes; i++) {
            igraph_real_t dx, dy, dz, dist, den, k_mi, l_mi;
            if (i == m) {
                continue;
            }
            dx = old_x - MATRIX(*res, i, 0);
            dy = old_y - MATRIX(*res, i, 1);
            dz = old_z - MATRIX(*res, i, 2);
            dist = sqrt(dx * dx + dy * dy + dz * dz);
            den = dist * (dx * dx + dy * dy + dz * dz);

            k_mi = MATRIX(kij, m, i);
            l_mi = MATRIX(lij, m, i);
            Axx += k_mi * (1 - l_mi * (dy * dy + dz * dz) / den);
            Ayy += k_mi * (1 - l_mi * (dx * dx + dz * dz) / den);
            Azz += k_mi * (1 - l_mi * (dx * dx + dy * dy) / den);
            Axy += k_mi * l_mi * dx * dy / den;
            Axz += k_mi * l_mi * dx * dz / den;
            Ayz += k_mi * l_mi * dy * dz / den;
        }
        Ax = -VECTOR(D1)[m];
        Ay = -VECTOR(D2)[m];
        Az = -VECTOR(D3)[m];

        /* Need to solve some linear equations, we just use Cramer's rule */
#define DET(a,b,c,d,e,f,g,h,i) ((a*e*i+b*f*g+c*d*h)-(c*e*g+b*d*i+a*f*h))

        /* See comments in 2D version for the reason for this check */
        if (Ax*Ax + Ay*Ay + Az*Az < KK_EPS*KK_EPS) {
            delta_x = delta_y = delta_z = 0;
        } else {
            igraph_real_t detnum;
            detnum  = DET(Axx, Axy, Axz, Axy, Ayy, Ayz, Axz, Ayz, Azz);
            delta_x = DET(Ax, Ay, Az, Axy, Ayy, Ayz, Axz, Ayz, Azz) / detnum;
            delta_y = DET(Axx, Axy, Axz, Ax, Ay, Az, Axz, Ayz, Azz) / detnum;
            delta_z = DET(Axx, Axy, Axz, Axy, Ayy, Ayz, Ax, Ay, Az ) / detnum;
        }

        new_x = old_x + delta_x;
        new_y = old_y + delta_y;
        new_z = old_z + delta_z;

        /* Limits, if given */
        if (minx && new_x < VECTOR(*minx)[m]) {
            new_x = VECTOR(*minx)[m];
        }
        if (maxx && new_x > VECTOR(*maxx)[m]) {
            new_x = VECTOR(*maxx)[m];
        }
        if (miny && new_y < VECTOR(*miny)[m]) {
            new_y = VECTOR(*miny)[m];
        }
        if (maxy && new_y > VECTOR(*maxy)[m]) {
            new_y = VECTOR(*maxy)[m];
        }
        if (minz && new_z < VECTOR(*minz)[m]) {
            new_z = VECTOR(*minz)[m];
        }
        if (maxz && new_z > VECTOR(*maxz)[m]) {
            new_z = VECTOR(*maxz)[m];
        }

        /* Update delta, only with/for the affected node */
        VECTOR(D1)[m] = VECTOR(D2)[m] = VECTOR(D3)[m] = 0.0;
        for (i = 0; i < no_nodes; i++) {
            igraph_real_t old_dx, old_dy, old_dz, old_mi_dist, new_dx, new_dy, new_dz, new_mi_dist;
            if (i == m) {
                continue;
            }
            old_dx = old_x - MATRIX(*res, i, 0);
            old_dy = old_y - MATRIX(*res, i, 1);
            old_dz = old_z - MATRIX(*res, i, 2);
            old_mi_dist = sqrt(old_dx * old_dx + old_dy * old_dy +
                               old_dz * old_dz);
            new_dx = new_x - MATRIX(*res, i, 0);
            new_dy = new_y - MATRIX(*res, i, 1);
            new_dz = new_z - MATRIX(*res, i, 2);
            new_mi_dist = sqrt(new_dx * new_dx + new_dy * new_dy +
                               new_dz * new_dz);

            VECTOR(D1)[i] -= MATRIX(kij, m, i) *
                             (-old_dx + MATRIX(lij, m, i) * old_dx / old_mi_dist);
            VECTOR(D2)[i] -= MATRIX(kij, m, i) *
                             (-old_dy + MATRIX(lij, m, i) * old_dy / old_mi_dist);
            VECTOR(D3)[i] -= MATRIX(kij, m, i) *
                             (-old_dz + MATRIX(lij, m, i) * old_dz / old_mi_dist);

            VECTOR(D1)[i] += MATRIX(kij, m, i) *
                             (-new_dx + MATRIX(lij, m, i) * new_dx / new_mi_dist);
            VECTOR(D2)[i] += MATRIX(kij, m, i) *
                             (-new_dy + MATRIX(lij, m, i) * new_dy / new_mi_dist);
            VECTOR(D3)[i] += MATRIX(kij, m, i) *
                             (-new_dz + MATRIX(lij, m, i) * new_dz / new_mi_dist);

            VECTOR(D1)[m] += MATRIX(kij, m, i) *
                             (new_dx - MATRIX(lij, m, i) * new_dx / new_mi_dist);
            VECTOR(D2)[m] += MATRIX(kij, m, i) *
                             (new_dy - MATRIX(lij, m, i) * new_dy / new_mi_dist);
            VECTOR(D3)[m] += MATRIX(kij, m, i) *
                             (new_dz - MATRIX(lij, m, i) * new_dz / new_mi_dist);
        }

        /* Update coordinates*/
        MATRIX(*res, m, 0) = new_x;
        MATRIX(*res, m, 1) = new_y;
        MATRIX(*res, m, 2) = new_z;
    }

    igraph_vector_destroy(&D3);
    igraph_vector_destroy(&D2);
    igraph_vector_destroy(&D1);
    igraph_matrix_destroy(&lij);
    igraph_matrix_destroy(&kij);
    igraph_matrix_destroy(&dij);
    IGRAPH_FINALLY_CLEAN(6);

    return IGRAPH_SUCCESS;
}
