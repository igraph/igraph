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
#include "layout/layout_internal.h"

#include <math.h>

/* not 'static', used in tests */
igraph_bool_t igraph_i_layout_segments_intersect(float p0_x, float p0_y,
        float p1_x, float p1_y,
        float p2_x, float p2_y,
        float p3_x, float p3_y) {
    float s1_x = p1_x - p0_x;
    float s1_y = p1_y - p0_y;
    float s2_x = p3_x - p2_x;
    float s2_y = p3_y - p2_y;

    float s1, s2, t1, t2, s, t;
    s1 = (-s1_y * (p0_x - p2_x) + s1_x * (p0_y - p2_y));
    s2 = (-s2_x * s1_y + s1_x * s2_y);
    if (s2 == 0) {
        return 0;
    }
    t1 = ( s2_x * (p0_y - p2_y) - s2_y * (p0_x - p2_x));
    t2 = (-s2_x * s1_y + s1_x * s2_y);
    s = s1 / s2;
    t = t1 / t2;

    return s >= 0 && s <= 1 && t >= 0 && t <= 1 ? 1 : 0;
}

/* not 'static', used in tests */
float igraph_i_layout_point_segment_dist2(float v_x, float v_y,
                                   float u1_x, float u1_y,
                                   float u2_x, float u2_y) {

    float dx = u2_x - u1_x;
    float dy = u2_y - u1_y;
    float l2 = dx * dx + dy * dy;
    float t, p_x, p_y;
    if (l2 == 0) {
        return (v_x - u1_x) * (v_x - u1_x) + (v_y - u1_y) * (v_y - u1_y);
    }
    t = ((v_x - u1_x) * dx + (v_y - u1_y) * dy) / l2;
    if (t < 0.0) {
        return (v_x - u1_x) * (v_x - u1_x) + (v_y - u1_y) * (v_y - u1_y);
    } else if (t > 1.0) {
        return (v_x - u2_x) * (v_x - u2_x) + (v_y - u2_y) * (v_y - u2_y);
    }
    p_x = u1_x + t * dx;
    p_y = u1_y + t * dy;
    return (v_x - p_x) * (v_x - p_x) + (v_y - p_y) * (v_y - p_y);
}

/**
 * \function igraph_layout_davidson_harel
 * Davidson-Harel layout algorithm
 *
 * This function implements the algorithm by Davidson and Harel,
 * see Ron Davidson, David Harel: Drawing Graphs Nicely Using
 * Simulated Annealing. ACM Transactions on Graphics 15(4),
 * pp. 301-331, 1996.
 *
 * </para><para>
 * The algorithm uses simulated annealing and a sophisticated
 * energy function, which is unfortunately hard to parameterize
 * for different graphs. The original publication did not disclose any
 * parameter values, and the ones below were determined by
 * experimentation.
 *
 * </para><para>
 * The algorithm consists of two phases, an annealing phase, and a
 * fine-tuning phase. There is no simulated annealing in the second
 * phase.
 *
 * </para><para>
 * Our implementation tries to follow the original publication, as
 * much as possible. The only major difference is that coordinates are
 * explicitly kept within the bounds of the rectangle of the layout.
 *
 * \param graph The input graph, edge directions are ignored.
 * \param res A matrix, the result is stored here. It can be used to
 *     supply start coordinates, see \p use_seed.
 * \param use_seed Boolean, whether to use the supplied \p res as
 *     start coordinates.
 * \param maxiter The maximum number of annealing iterations. A
 *     reasonable value for smaller graphs is 10.
 * \param fineiter The number of fine tuning iterations. A reasonable
 *     value is max(10, log2(n)) where n is the number of vertices.
 * \param cool_fact Cooling factor. A reasonable value is 0.75.
 * \param weight_node_dist Weight for the node-node distances
 *     component of the energy function. Reasonable value: 1.0.
 * \param weight_border Weight for the distance from the border
 *     component of the energy function. It can be set to zero, if
 *     vertices are allowed to sit on the border.
 * \param weight_edge_lengths Weight for the edge length component
 *     of the energy function, a reasonable value is the density of
 *     the graph divided by 10.
 * \param weight_edge_crossings Weight for the edge crossing component
 *     of the energy function, a reasonable default is 1 minus the
 *     square root of the density of the graph.
 * \param weight_node_edge_dist Weight for the node-edge distance
 *     component of the energy function. A reasonable value is
 *     1 minus the density, divided by 5.
 * \return Error code.
 *
 * Time complexity: one first phase iteration has time complexity
 * O(n^2+m^2), one fine tuning iteration has time complexity O(mn).
 * Time complexity might be smaller if some of the weights of the
 * components of the energy function are set to zero.
 *
 */

int igraph_layout_davidson_harel(const igraph_t *graph, igraph_matrix_t *res,
                                 igraph_bool_t use_seed, igraph_integer_t maxiter,
                                 igraph_integer_t fineiter, igraph_real_t cool_fact,
                                 igraph_real_t weight_node_dist, igraph_real_t weight_border,
                                 igraph_real_t weight_edge_lengths,
                                 igraph_real_t weight_edge_crossings,
                                 igraph_real_t weight_node_edge_dist) {

    igraph_integer_t no_nodes = igraph_vcount(graph);
    igraph_integer_t no_edges = igraph_ecount(graph);
    float width = sqrt(no_nodes) * 10, height = width;
    igraph_vector_int_t perm;
    igraph_bool_t fine_tuning = 0;
    igraph_integer_t round, i;
    igraph_vector_float_t try_x, try_y;
    igraph_vector_int_t try_idx;
    float move_radius = width / 2;
    float fine_tuning_factor = 0.01f;
    igraph_vector_t neis;
    float min_x = width / 2, max_x = -width / 2, min_y = height / 2, max_y = -height / 2;

    igraph_integer_t no_tries = 30;
    float w_node_dist = weight_node_dist ;          /* 1.0 */
    float w_borderlines = weight_border;            /* 0.0 */
    float w_edge_lengths = weight_edge_lengths;     /* 0.0001; */
    float w_edge_crossings = weight_edge_crossings; /* 1.0 */
    float w_node_edge_dist = weight_node_edge_dist; /* 0.2 */

    if (use_seed && (igraph_matrix_nrow(res) != no_nodes ||
                     igraph_matrix_ncol(res) != 2)) {
        IGRAPH_ERROR("Invalid start position matrix size in "
                     "Davidson-Harel layout", IGRAPH_EINVAL);
    }
    if (maxiter < 0) {
        IGRAPH_ERROR("Number of iterations must be non-negative in "
                     "Davidson-Harel layout", IGRAPH_EINVAL);
    }
    if (fineiter < 0) {
        IGRAPH_ERROR("Number of fine tuning iterations must be non-negative in "
                     "Davidson-Harel layout", IGRAPH_EINVAL);
    }
    if (cool_fact <= 0 || cool_fact >= 1) {
        IGRAPH_ERROR("Cooling factor must be in (0,1) in "
                     "Davidson-Harel layout", IGRAPH_EINVAL);
    }

    if (no_nodes == 0) {
        return 0;
    }

    IGRAPH_CHECK(igraph_vector_int_init_seq(&perm, 0, no_nodes - 1));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &perm);
    IGRAPH_CHECK(igraph_vector_float_init(&try_x, no_tries));
    IGRAPH_FINALLY(igraph_vector_float_destroy, &try_x);
    IGRAPH_CHECK(igraph_vector_float_init(&try_y, no_tries));
    IGRAPH_FINALLY(igraph_vector_float_destroy, &try_y);
    IGRAPH_CHECK(igraph_vector_int_init_seq(&try_idx, 0, no_tries - 1));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &try_idx);
    IGRAPH_VECTOR_INIT_FINALLY(&neis, 100);

    RNG_BEGIN();

    if (!use_seed) {
        IGRAPH_CHECK(igraph_matrix_resize(res, no_nodes, 2));
        for (i = 0; i < no_nodes; i++) {
            float x, y;
            x = MATRIX(*res, i, 0) = RNG_UNIF(-width / 2, width / 2);
            y = MATRIX(*res, i, 1) = RNG_UNIF(-height / 2, height / 2);
            if (x < min_x) {
                min_x = x;
            } else if (x > max_x) {
                max_x = x;
            }
            if (y < min_y) {
                min_y = y;
            } else if (y > max_y) {
                max_y = y;
            }
        }
    } else {
        min_x = IGRAPH_INFINITY; max_x = IGRAPH_NEGINFINITY;
        min_y = IGRAPH_INFINITY; max_y = IGRAPH_NEGINFINITY;
        for (i = 0; i < no_nodes; i++) {
            float x = MATRIX(*res, i, 0);
            float y = MATRIX(*res, i, 1);
            if (x < min_x) {
                min_x = x;
            } else if (x > max_x) {
                max_x = x;
            }
            if (y < min_y) {
                min_y = y;
            } else if (y > max_y) {
                max_y = y;
            }
        }
    }

    for (i = 0; i < no_tries; i++) {
        float phi = 2 * M_PI / no_tries * i;
        VECTOR(try_x)[i] = cos(phi);
        VECTOR(try_y)[i] = sin(phi);
    }

    for (round = 0; round < maxiter + fineiter; round++) {
        igraph_integer_t p;
        igraph_vector_int_shuffle(&perm);

        IGRAPH_ALLOW_INTERRUPTION();

        fine_tuning = round >= maxiter;
        if (fine_tuning) {
            float fx = fine_tuning_factor * (max_x - min_x);
            float fy = fine_tuning_factor * (max_y - min_y);
            move_radius = fx < fy ? fx : fy;
        }

        for (p = 0; p < no_nodes; p++) {
            igraph_integer_t t;
            igraph_integer_t v = VECTOR(perm)[p];
            igraph_vector_int_shuffle(&try_idx);

            for (t = 0; t < no_tries; t++) {
                float diff_energy = 0.0;
                int ti = VECTOR(try_idx)[t];

                /* Try moving it */
                float old_x = MATRIX(*res, v, 0);
                float old_y = MATRIX(*res, v, 1);
                float new_x = old_x + move_radius * VECTOR(try_x)[ti];
                float new_y = old_y + move_radius * VECTOR(try_y)[ti];

                if (new_x < -width / 2) {
                    new_x = -width / 2 - 1e-6;
                }
                if (new_x >  width / 2) {
                    new_x =  width / 2 - 1e-6;
                }
                if (new_y < -height / 2) {
                    new_y = -height / 2 - 1e-6;
                }
                if (new_y >  height / 2) {
                    new_y =  height / 2 - 1e-6;
                }

                if (w_node_dist != 0) {
                    igraph_integer_t u;
                    for (u = 0; u < no_nodes; u++) {
                        float odx, ody, odist2, dx, dy, dist2;
                        if (u == v) {
                            continue;
                        }
                        odx = old_x - MATRIX(*res, u, 0);
                        ody = old_y - MATRIX(*res, u, 1);
                        dx = new_x - MATRIX(*res, u, 0);
                        dy = new_y - MATRIX(*res, u, 1);
                        odist2 = odx * odx + ody * ody;
                        dist2 = dx * dx + dy * dy;
                        diff_energy += w_node_dist / dist2 - w_node_dist / odist2;
                    }
                }

                if (w_borderlines != 0) {
                    float odx1 = width / 2 - old_x, odx2 = old_x + width / 2;
                    float ody1 = height / 2 - old_y, ody2 = old_y + height / 2;
                    float dx1 = width / 2 - new_x, dx2 = new_x + width / 2;
                    float dy1 = height / 2 - new_y, dy2 = new_y + height / 2;
                    if (odx1 < 0) {
                        odx1 = 2;
                    } if (odx2 < 0) {
                        odx2 = 2;
                    }
                    if (ody1 < 0) {
                        ody1 = 2;
                    } if (ody2 < 0) {
                        ody2 = 2;
                    }
                    if (dx1 < 0) {
                        dx1 = 2;
                    } if (dx2 < 0) {
                        dx2 = 2;
                    }
                    if (dy1 < 0) {
                        dy1 = 2;
                    } if (dy2 < 0) {
                        dy2 = 2;
                    }
                    diff_energy -= w_borderlines *
                                   (1.0 / (odx1 * odx1) + 1.0 / (odx2 * odx2) +
                                    1.0 / (ody1 * ody1) + 1.0 / (ody2 * ody2));
                    diff_energy += w_borderlines *
                                   (1.0 / (dx1 * dx1) + 1.0 / (dx2 * dx2) +
                                    1.0 / (dy1 * dy1) + 1.0 / (dy2 * dy2));
                }

                if (w_edge_lengths != 0) {
                    igraph_integer_t len, j;
                    igraph_neighbors(graph, &neis, v, IGRAPH_ALL);
                    len = igraph_vector_size(&neis);
                    for (j = 0; j < len; j++) {
                        igraph_integer_t u = VECTOR(neis)[j];
                        float odx = old_x - MATRIX(*res, u, 0);
                        float ody = old_y - MATRIX(*res, u, 1);
                        float odist2 = odx * odx + ody * ody;
                        float dx = new_x - MATRIX(*res, u, 0);
                        float dy = new_y - MATRIX(*res, u, 1);
                        float dist2 = dx * dx + dy * dy;
                        diff_energy += w_edge_lengths * (dist2 - odist2);
                    }
                }

                if (w_edge_crossings != 0) {
                    igraph_integer_t len, j, no = 0;
                    igraph_neighbors(graph, &neis, v, IGRAPH_ALL);
                    len = igraph_vector_size(&neis);
                    for (j = 0; j < len; j++) {
                        igraph_integer_t u = VECTOR(neis)[j];
                        float u_x = MATRIX(*res, u, 0);
                        float u_y = MATRIX(*res, u, 1);
                        igraph_integer_t e;
                        for (e = 0; e < no_edges; e++) {
                            igraph_integer_t u1 = IGRAPH_FROM(graph, e);
                            igraph_integer_t u2 = IGRAPH_TO(graph, e);
                            float u1_x, u1_y, u2_x, u2_y;
                            if (u1 == v || u2 == v || u1 == u || u2 == u) {
                                continue;
                            }
                            u1_x = MATRIX(*res, u1, 0);
                            u1_y = MATRIX(*res, u1, 1);
                            u2_x = MATRIX(*res, u2, 0);
                            u2_y = MATRIX(*res, u2, 1);
                            no -= igraph_i_layout_segments_intersect(old_x, old_y, u_x, u_y,
                                                              u1_x, u1_y, u2_x, u2_y);
                            no += igraph_i_layout_segments_intersect(new_x, new_y, u_x, u_y,
                                                              u1_x, u1_y, u2_x, u2_y);
                        }
                    }
                    diff_energy += w_edge_crossings * no;
                }

                if (w_node_edge_dist != 0 && fine_tuning) {
                    igraph_integer_t e, no;

                    /* All non-incident edges from the moved 'v' */
                    for (e = 0; e < no_edges; e++) {
                        igraph_integer_t u1 = IGRAPH_FROM(graph, e);
                        igraph_integer_t u2 = IGRAPH_TO(graph, e);
                        float u1_x, u1_y, u2_x, u2_y, d_ev;
                        if (u1 == v || u2 == v) {
                            continue;
                        }
                        u1_x = MATRIX(*res, u1, 0);
                        u1_y = MATRIX(*res, u1, 1);
                        u2_x = MATRIX(*res, u2, 0);
                        u2_y = MATRIX(*res, u2, 1);
                        d_ev = igraph_i_layout_point_segment_dist2(old_x, old_y, u1_x, u1_y,
                                                            u2_x, u2_y);
                        diff_energy -= w_node_edge_dist / d_ev;
                        d_ev = igraph_i_layout_point_segment_dist2(new_x, new_y, u1_x, u1_y,
                                                            u2_x, u2_y);
                        diff_energy += w_node_edge_dist / d_ev;
                    }

                    /* All other nodes from all of v's incident edges */
                    igraph_incident(graph, &neis, v, IGRAPH_ALL);
                    no = igraph_vector_size(&neis);
                    for (e = 0; e < no; e++) {
                        igraph_integer_t mye = VECTOR(neis)[e];
                        igraph_integer_t u = IGRAPH_OTHER(graph, mye, v);
                        float u_x = MATRIX(*res, u, 0);
                        float u_y = MATRIX(*res, u, 1);
                        igraph_integer_t w;
                        for (w = 0; w < no_nodes; w++) {
                            float w_x, w_y, d_ev;
                            if (w == v || w == u) {
                                continue;
                            }
                            w_x = MATRIX(*res, w, 0);
                            w_y = MATRIX(*res, w, 1);
                            d_ev = igraph_i_layout_point_segment_dist2(w_x, w_y, old_x,
                                                                old_y, u_x, u_y);
                            diff_energy -= w_node_edge_dist / d_ev;
                            d_ev = igraph_i_layout_point_segment_dist2(w_x, w_y, new_x, new_y,
                                                                u_x, u_y);
                            diff_energy += w_node_edge_dist / d_ev;
                        }
                    }
                } /* w_node_edge_dist != 0 && fine_tuning */

                if (diff_energy < 0 ||
                    (!fine_tuning && RNG_UNIF01() < exp(-diff_energy / move_radius))) {
                    MATRIX(*res, v, 0) = new_x;
                    MATRIX(*res, v, 1) = new_y;
                    if (new_x < min_x) {
                        min_x = new_x;
                    } else if (new_x > max_x) {
                        max_x = new_x;
                    }
                    if (new_y < min_y) {
                        min_y = new_y;
                    } else if (new_y > max_y) {
                        max_y = new_y;
                    }
                }

            } /* t < no_tries */

        } /* p < no_nodes  */

        move_radius *= cool_fact;

    } /* round < maxiter */

    RNG_END();

    igraph_vector_destroy(&neis);
    igraph_vector_int_destroy(&try_idx);
    igraph_vector_float_destroy(&try_x);
    igraph_vector_float_destroy(&try_y);
    igraph_vector_int_destroy(&perm);
    IGRAPH_FINALLY_CLEAN(5);

    return 0;
}
