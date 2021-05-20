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

#include "core/math.h"
#include "core/interruption.h"

/**
 * \ingroup layout
 * \function igraph_layout_gem
 *
 * The GEM layout algorithm, as described in Arne Frick, Andreas Ludwig,
 * Heiko Mehldau: A Fast Adaptive Layout Algorithm for Undirected Graphs,
 * Proc. Graph Drawing 1994, LNCS 894, pp. 388-403, 1995.
 * \param graph The input graph. Edge directions are ignored in
 *        directed graphs.
 * \param res The result is stored here. If the \p use_seed argument
 *        is true (non-zero), then this matrix is also used as the
 *        starting point of the algorithm.
 * \param use_seed Boolean, whether to use the supplied coordinates in
 *        \p res as the starting point. If false (zero), then a
 *        uniform random starting point is used.
 * \param maxiter The maximum number of iterations to
 *        perform. Updating a single vertex counts as an iteration.
 *        A reasonable default is 40 * n * n, where n is the number of
 *        vertices. The original paper suggests 4 * n * n, but this
 *        usually only works if the other parameters are set up carefully.
 * \param temp_max The maximum allowed local temperature. A reasonable
 *        default is the number of vertices.
 * \param temp_min The global temperature at which the algorithm
 *        terminates (even before reaching \p maxiter iterations). A
 *        reasonable default is 1/10.
 * \param temp_init Initial local temperature of all vertices. A
 *        reasonable default is the square root of the number of
 *        vertices.
 * \return Error code.
 *
 * Time complexity: O(t * n * (n+e)), where n is the number of vertices,
 * e is the number of edges and t is the number of time steps
 * performed.
 */

int igraph_layout_gem(const igraph_t *graph, igraph_matrix_t *res,
                      igraph_bool_t use_seed, igraph_integer_t maxiter,
                      igraph_real_t temp_max, igraph_real_t temp_min,
                      igraph_real_t temp_init) {

    igraph_integer_t no_nodes = igraph_vcount(graph);
    igraph_vector_int_t perm;
    igraph_vector_float_t impulse_x, impulse_y, temp, skew_gauge;
    igraph_integer_t i;
    float temp_global;
    igraph_integer_t perm_pointer = 0;
    float barycenter_x = 0, barycenter_y = 0;
    igraph_vector_t phi;
    igraph_vector_t neis;
    const float elen_des2 = 128 * 128;
    const float gamma = 1 / 16.0f;
    const float alpha_o = (float)M_PI;
    const float alpha_r = (float)M_PI / 3.0f;
    const float sigma_o = 1.0f / 3.0f;
    const float sigma_r = 1.0f / 2.0f / no_nodes;

    if (maxiter < 0) {
        IGRAPH_ERROR("Number of iterations must be non-negative in GEM layout",
                     IGRAPH_EINVAL);
    }
    if (use_seed && (igraph_matrix_nrow(res) != no_nodes ||
                     igraph_matrix_ncol(res) != 2)) {
        IGRAPH_ERROR("Invalid start position matrix size in GEM layout",
                     IGRAPH_EINVAL);
    }
    if (temp_max <= 0) {
        IGRAPH_ERROR("Maximum temperature should be positive in GEM layout",
                     IGRAPH_EINVAL);
    }
    if (temp_min <= 0) {
        IGRAPH_ERROR("Minimum temperature should be positive in GEM layout",
                     IGRAPH_EINVAL);
    }
    if (temp_init <= 0) {
        IGRAPH_ERROR("Initial temperature should be positive in GEM layout",
                     IGRAPH_EINVAL);
    }
    if (temp_max < temp_init || temp_init < temp_min) {
        IGRAPH_ERROR("Minimum <= Initial <= Maximum temperature is required "
                     "in GEM layout", IGRAPH_EINVAL);
    }

    if (no_nodes == 0) {
        return 0;
    }

    IGRAPH_CHECK(igraph_vector_float_init(&impulse_x, no_nodes));
    IGRAPH_FINALLY(igraph_vector_float_destroy, &impulse_x);
    IGRAPH_CHECK(igraph_vector_float_init(&impulse_y, no_nodes));
    IGRAPH_FINALLY(igraph_vector_float_destroy, &impulse_y);
    IGRAPH_CHECK(igraph_vector_float_init(&temp, no_nodes));
    IGRAPH_FINALLY(igraph_vector_float_destroy, &temp);
    IGRAPH_CHECK(igraph_vector_float_init(&skew_gauge, no_nodes));
    IGRAPH_FINALLY(igraph_vector_float_destroy, &skew_gauge);
    IGRAPH_CHECK(igraph_vector_int_init_seq(&perm, 0, no_nodes - 1));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &perm);
    IGRAPH_VECTOR_INIT_FINALLY(&phi, no_nodes);
    IGRAPH_VECTOR_INIT_FINALLY(&neis, 10);

    RNG_BEGIN();

    /* Initialization */
    igraph_degree(graph, &phi, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS);
    if (!use_seed) {
        const igraph_real_t width_half = no_nodes * 100, height_half = width_half;
        IGRAPH_CHECK(igraph_matrix_resize(res, no_nodes, 2));
        for (i = 0; i < no_nodes; i++) {
            MATRIX(*res, i, 0) = RNG_UNIF(-width_half, width_half);
            MATRIX(*res, i, 1) = RNG_UNIF(-height_half, height_half);
            barycenter_x += MATRIX(*res, i, 0);
            barycenter_y += MATRIX(*res, i, 1);
            VECTOR(phi)[i] *= (VECTOR(phi)[i] / 2.0 + 1.0);
        }
    } else {
        for (i = 0; i < no_nodes; i++) {
            barycenter_x += MATRIX(*res, i, 0);
            barycenter_y += MATRIX(*res, i, 1);
            VECTOR(phi)[i] *= (VECTOR(phi)[i] / 2.0 + 1.0);
        }
    }
    igraph_vector_float_fill(&temp, temp_init);
    temp_global = temp_init * no_nodes;

    while (temp_global > temp_min * no_nodes && maxiter > 0) {
        igraph_integer_t u, v, nlen, j;
        float px, py, pvx, pvy;

        IGRAPH_ALLOW_INTERRUPTION();

        /* choose a vertex v to update */
        if (!perm_pointer) {
            igraph_vector_int_shuffle(&perm);
            perm_pointer = no_nodes - 1;
        }
        v = VECTOR(perm)[perm_pointer--];

        /* compute v's impulse */
        px = (barycenter_x / no_nodes - MATRIX(*res, v, 0)) * gamma * VECTOR(phi)[v];
        py = (barycenter_y / no_nodes - MATRIX(*res, v, 1)) * gamma * VECTOR(phi)[v];
        px += RNG_UNIF(-32.0, 32.0);
        py += RNG_UNIF(-32.0, 32.0);

        for (u = 0; u < no_nodes; u++) {
            float dx, dy, dist2;
            if (u == v) {
                continue;
            }
            dx = MATRIX(*res, v, 0) - MATRIX(*res, u, 0);
            dy = MATRIX(*res, v, 1) - MATRIX(*res, u, 1);
            dist2 = dx * dx + dy * dy;
            if (dist2 != 0) {
                px += dx * elen_des2 / dist2;
                py += dy * elen_des2 / dist2;
            }
        }

        IGRAPH_CHECK(igraph_neighbors(graph, &neis, v, IGRAPH_ALL));
        nlen = igraph_vector_size(&neis);
        for (j = 0; j < nlen; j++) {
            igraph_integer_t u = VECTOR(neis)[j];
            float dx = MATRIX(*res, v, 0) - MATRIX(*res, u, 0);
            float dy = MATRIX(*res, v, 1) - MATRIX(*res, u, 1);
            float dist2 = dx * dx + dy * dy;
            px -= dx * dist2 / (elen_des2 * VECTOR(phi)[v]);
            py -= dy * dist2 / (elen_des2 * VECTOR(phi)[v]);
        }

        /* update v's position and temperature */
        if (px != 0 || py != 0) {
            float plen = sqrtf(px * px + py * py);
            px *= VECTOR(temp)[v] / plen;
            py *= VECTOR(temp)[v] / plen;
            MATRIX(*res, v, 0) += px;
            MATRIX(*res, v, 1) += py;
            barycenter_x += px;
            barycenter_y += py;
        }

        pvx = VECTOR(impulse_x)[v]; pvy = VECTOR(impulse_y)[v];
        if (pvx != 0 || pvy != 0) {
            float beta = atan2f(pvy - py, pvx - px);
            float sin_beta = sinf(beta);
            float sign_sin_beta = (sin_beta > 0) ? 1 : ((sin_beta < 0) ? -1 : 0);
            float cos_beta = cosf(beta);
            float abs_cos_beta = fabsf(cos_beta);
            float old_temp = VECTOR(temp)[v];
            if (sin(beta) >= sin(M_PI_2 + alpha_r / 2.0)) {
                VECTOR(skew_gauge)[v] += sigma_r * sign_sin_beta;
            }
            if (abs_cos_beta >= cosf(alpha_o / 2.0)) {
                VECTOR(temp)[v] *= sigma_o * cos_beta;
            }
            VECTOR(temp)[v] *= (1 - fabsf(VECTOR(skew_gauge)[v]));
            if (VECTOR(temp)[v] > temp_max) {
                VECTOR(temp)[v] = temp_max;
            }
            VECTOR(impulse_x)[v] = px;
            VECTOR(impulse_y)[v] = py;
            temp_global += VECTOR(temp)[v] - old_temp;
        }

        maxiter--;

    } /* while temp && iter */


    RNG_END();

    igraph_vector_destroy(&neis);
    igraph_vector_destroy(&phi);
    igraph_vector_int_destroy(&perm);
    igraph_vector_float_destroy(&skew_gauge);
    igraph_vector_float_destroy(&temp);
    igraph_vector_float_destroy(&impulse_y);
    igraph_vector_float_destroy(&impulse_x);
    IGRAPH_FINALLY_CLEAN(7);

    return 0;
}
