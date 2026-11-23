/*
   igraph library.
   Copyright (C) 2026  The igraph development team <igraph@igraph.org>

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include <igraph.h>
#include "test_utilities.h"

/* Points that lie exactly on the boundary of a closed exclusion region must
 * remove the edge, also when the coordinates are far from the origin compared
 * to the distances between the points. The tolerance that ensures this must
 * not break the containment of skeletons in each other, or connectivity. */

/* Checks that the graph has exactly the given edges. */
static void check_edges(const igraph_t *graph, const igraph_int_t *edges, igraph_int_t count) {
    IGRAPH_ASSERT(igraph_ecount(graph) == count);
    for (igraph_int_t i = 0; i < count; i++) {
        igraph_bool_t adj;
        igraph_are_adjacent(graph, edges[2 * i], edges[2 * i + 1], &adj);
        IGRAPH_ASSERT(adj);
    }
}

/* Checks that two graphs on the same vertices have the same edges, in any order. */
static void check_same_edges(const igraph_t *graph1, const igraph_t *graph2) {
    IGRAPH_ASSERT(igraph_ecount(graph1) == igraph_ecount(graph2));
    for (igraph_int_t e = 0; e < igraph_ecount(graph1); e++) {
        igraph_bool_t adj;
        igraph_are_adjacent(graph2, IGRAPH_FROM(graph1, e), IGRAPH_TO(graph1, e), &adj);
        IGRAPH_ASSERT(adj);
    }
}

/* Checks that the Gabriel graph, the beta = 1 skeletons and the weighted Gabriel
 * graph with several cutoffs are identical, and returns the Gabriel graph. */
static void gabriel_identity(const igraph_matrix_t *points, igraph_t *gabriel) {
    const igraph_real_t cutoffs[] = { 1, 1.5, 2, IGRAPH_INFINITY };
    igraph_t graph;
    igraph_vector_t weights;
    igraph_bool_t same;

    igraph_vector_init(&weights, 0);
    igraph_gabriel_graph(gabriel, points);

    igraph_lune_beta_skeleton(&graph, points, 1);
    igraph_is_same_graph(gabriel, &graph, &same);
    IGRAPH_ASSERT(same);
    igraph_destroy(&graph);

    if (igraph_matrix_ncol(points) == 2) {
        igraph_circle_beta_skeleton(&graph, points, 1);
        igraph_is_same_graph(gabriel, &graph, &same);
        IGRAPH_ASSERT(same);
        igraph_destroy(&graph);
    }

    for (size_t i = 0; i < sizeof(cutoffs) / sizeof(cutoffs[0]); i++) {
        igraph_beta_weighted_gabriel_graph(&graph, &weights, points, cutoffs[i]);
        igraph_is_same_graph(gabriel, &graph, &same);
        IGRAPH_ASSERT(same);
        igraph_destroy(&graph);
    }

    igraph_vector_destroy(&weights);
}

/* Checks that every edge of 'small' is also an edge of 'big'. */
static void check_subgraph(const igraph_t *small, const igraph_t *big) {
    for (igraph_int_t e = 0; e < igraph_ecount(small); e++) {
        igraph_bool_t adj;
        igraph_are_adjacent(big, IGRAPH_FROM(small, e), IGRAPH_TO(small, e), &adj);
        IGRAPH_ASSERT(adj);
    }
}

static void check_connected(const igraph_t *graph) {
    igraph_bool_t connected;
    igraph_is_connected(graph, &connected, IGRAPH_WEAK);
    IGRAPH_ASSERT(connected);
}

typedef enum { LUNE, CIRCLE, RNG } skeleton_t;

static void skeleton(igraph_t *graph, const igraph_matrix_t *points, skeleton_t type, igraph_real_t beta) {
    switch (type) {
    case LUNE: igraph_lune_beta_skeleton(graph, points, beta); break;
    case CIRCLE: igraph_circle_beta_skeleton(graph, points, beta); break;
    case RNG: igraph_relative_neighborhood_graph(graph, points); break;
    }
}

/* Checks that each skeleton of the sequence contains the next one, and that
 * the skeletons are connected where required. */
static void check_sequence(const igraph_matrix_t *points, const skeleton_t *types,
                           const igraph_real_t *betas, igraph_int_t count, igraph_real_t connected_below) {
    igraph_t prev, next;
    skeleton(&prev, points, types[0], betas[0]);
    for (igraph_int_t i = 0; i < count; i++) {
        if (betas[i] < connected_below) {
            check_connected(&prev);
        }
        if (i + 1 < count) {
            skeleton(&next, points, types[i + 1], betas[i + 1]);
            check_subgraph(&next, &prev);
            igraph_destroy(&prev);
            prev = next;
        }
    }
    igraph_destroy(&prev);
}

/* Checks that the lune-based and circle-based skeletons are nested as beta
 * increases, with the relative neighborhood graph between beta < 2 and beta = 2,
 * that the skeletons are connected where they must be, and in 2D that circle-based
 * skeletons are contained in lune-based ones. */
static void check_nesting(const igraph_matrix_t *points) {
    const skeleton_t lune_types[] = { LUNE, LUNE, LUNE, LUNE, LUNE, RNG, LUNE, LUNE };
    const igraph_real_t lune_betas[] = { 0.3, 0.7, 1, 1.5, 1.99, 1.999, 2, 3 };
    const skeleton_t circle_types[] = { CIRCLE, CIRCLE, CIRCLE, CIRCLE, CIRCLE, CIRCLE, CIRCLE };
    const igraph_real_t circle_betas[] = { 0.3, 0.7, 1, 1.05, 2 / sqrt(3.0), 1.5, 3 };

    if (igraph_matrix_ncol(points) == 2) {
        igraph_t lune, circle;

        check_sequence(points, lune_types, lune_betas, 8, 2);
        check_sequence(points, circle_types, circle_betas, 7, 1.01);
        for (igraph_int_t i = 0; i < 7; i++) {
            skeleton(&lune, points, LUNE, circle_betas[i]);
            skeleton(&circle, points, CIRCLE, circle_betas[i]);
            check_subgraph(&circle, &lune);
            if (circle_betas[i] <= 1) {
                check_subgraph(&lune, &circle);
            }
            igraph_destroy(&circle);
            igraph_destroy(&lune);
        }
    } else {
        /* Lune-based skeletons with beta < 1 are only available in 2D. */
        check_sequence(points, lune_types + 2, lune_betas + 2, 6, 2);
    }
}

/* Checks that all edges have length 'spacing', and that there are 'count' of them. */
static void check_unit_edges(const igraph_t *graph, const igraph_matrix_t *points, igraph_real_t spacing, igraph_int_t count) {
    IGRAPH_ASSERT(igraph_ecount(graph) == count);
    for (igraph_int_t e = 0; e < count; e++) {
        igraph_real_t sqr_length = 0;
        for (igraph_int_t k = 0; k < igraph_matrix_ncol(points); k++) {
            const igraph_real_t d = MATRIX(*points, IGRAPH_FROM(graph, e), k) - MATRIX(*points, IGRAPH_TO(graph, e), k);
            sqr_length += d * d;
        }
        IGRAPH_ASSERT(fabs(sqrt(sqr_length) / spacing - 1) < 1e-6);
    }
}

int main(void) {
    igraph_t graph;
    igraph_matrix_t points;
    igraph_vector_t weights;

    igraph_vector_init(&weights, 0);

    {
        /* Right angle at point 2 for the stored coordinates, at two offsets. */
        const igraph_real_t right_far[] = { 1e6, 1e6,  1e6 + 0.3, 1e6 + 0.6,  1e6 + 0.3, 1e6 };
        const igraph_real_t right_near[] = { 19, 19,  19.1, 19.2,  19.1, 19 };
        const igraph_int_t expected[] = { 0, 2,  1, 2 };

        igraph_matrix_init_array(&points, right_far, 3, 2, IGRAPH_ROW_MAJOR);
        gabriel_identity(&points, &graph);
        check_edges(&graph, expected, 2);
        igraph_destroy(&graph);
        igraph_matrix_destroy(&points);

        igraph_matrix_init_array(&points, right_near, 3, 2, IGRAPH_ROW_MAJOR);
        gabriel_identity(&points, &graph);
        check_edges(&graph, expected, 2);
        igraph_destroy(&graph);
        igraph_matrix_destroy(&points);
    }

    {
        /* Four cells of a 30 arc-second grid in longitude and latitude. Point 0 is
         * the right-angle corner of points 2 and 3, so 2-3 is not a Gabriel edge. */
        const igraph_real_t grid[] = {
            19.058333333333334, 47.50833333333333,
            19.05, 47.516666666666666,
            19.058333333333334, 47.516666666666666,
            19.066666666666666, 47.50833333333333
        };
        const igraph_int_t expected[] = { 0, 2,  0, 3,  1, 2 };

        igraph_matrix_init_array(&points, grid, 4, 2, IGRAPH_ROW_MAJOR);
        gabriel_identity(&points, &graph);
        check_edges(&graph, expected, 3);
        igraph_destroy(&graph);
        igraph_matrix_destroy(&points);
    }

    {
        /* An exactly equilateral triangle in 3D: the relative neighborhood graph
         * keeps all edges, the closed beta = 2 lune removes all of them. */
        const igraph_real_t tri[] = {
            1e6, 1e6, 1e6,
            1e6 + 0.1, 1e6 + 0.1, 1e6,
            1e6 + 0.1, 1e6, 1e6 + 0.1
        };

        igraph_matrix_init_array(&points, tri, 3, 3, IGRAPH_ROW_MAJOR);
        igraph_relative_neighborhood_graph(&graph, &points);
        IGRAPH_ASSERT(igraph_ecount(&graph) == 3);
        igraph_destroy(&graph);
        igraph_lune_beta_skeleton(&graph, &points, 2);
        IGRAPH_ASSERT(igraph_ecount(&graph) == 0);
        igraph_destroy(&graph);
        igraph_matrix_destroy(&points);
    }

    {
        /* Point 2 is on both spheres of the beta = 5 lune of edge 0-1, so the
         * edge disappears at exactly beta = 5. The other edges disappear at 2. */
        const igraph_real_t apex[] = { 1e6 - 1, 1e6,  1e6 + 1, 1e6,  1e6, 1e6 + 3 };
        const igraph_int_t gabriel[] = { 0, 1,  0, 2,  1, 2 };
        igraph_int_t eid;

        igraph_matrix_init_array(&points, apex, 3, 2, IGRAPH_ROW_MAJOR);

        igraph_lune_beta_skeleton(&graph, &points, 4.9);
        check_edges(&graph, gabriel, 1);
        igraph_destroy(&graph);

        igraph_lune_beta_skeleton(&graph, &points, 5);
        IGRAPH_ASSERT(igraph_ecount(&graph) == 0);
        igraph_destroy(&graph);

        /* The weight of 0-1 is the threshold 5, also when the cutoff equals it. */
        igraph_beta_weighted_gabriel_graph(&graph, &weights, &points, 5);
        check_edges(&graph, gabriel, 3);
        igraph_get_eid(&graph, &eid, 0, 1, IGRAPH_UNDIRECTED, true);
        IGRAPH_ASSERT(VECTOR(weights)[eid] == 5);
        igraph_get_eid(&graph, &eid, 0, 2, IGRAPH_UNDIRECTED, true);
        IGRAPH_ASSERT(VECTOR(weights)[eid] == 2);
        igraph_destroy(&graph);

        igraph_beta_weighted_gabriel_graph(&graph, &weights, &points, 4.9);
        igraph_get_eid(&graph, &eid, 0, 1, IGRAPH_UNDIRECTED, true);
        IGRAPH_ASSERT(VECTOR(weights)[eid] == IGRAPH_INFINITY);
        igraph_destroy(&graph);

        igraph_matrix_destroy(&points);
    }

    {
        /* Unit square: the diagonal from the Delaunay triangulation has two points
         * on its Gabriel circle. All weights are infinite with cutoff 1. */
        const igraph_real_t square[] = { 0, 0,  1, 0,  1, 1,  0, 1 };

        igraph_matrix_init_array(&points, square, 4, 2, IGRAPH_ROW_MAJOR);
        gabriel_identity(&points, &graph);
        IGRAPH_ASSERT(igraph_ecount(&graph) == 4);
        igraph_destroy(&graph);

        igraph_beta_weighted_gabriel_graph(&graph, &weights, &points, 1);
        IGRAPH_ASSERT(igraph_vector_size(&weights) == 4);
        IGRAPH_ASSERT(igraph_vector_min(&weights) == IGRAPH_INFINITY);
        igraph_destroy(&graph);

        CHECK_ERROR(igraph_beta_weighted_gabriel_graph(&graph, &weights, &points, 0.5), IGRAPH_EINVAL);
        CHECK_ERROR(igraph_beta_weighted_gabriel_graph(&graph, &weights, &points, IGRAPH_NAN), IGRAPH_EINVAL);

        igraph_matrix_destroy(&points);
    }

    {
        /* Point 2 is near endpoint 0, and outside the Gabriel circle of 0-1 by
         * less than the tolerance. Whatever the decision, the functions agree. */
        /* Radius 1 + 32 * DBL_EPSILON, angle pi - 0.1. */
        const igraph_real_t near_circle[] = { -1, 0,  1, 0,  -0.9950041652780328, 0.09983341664682907 };

        igraph_matrix_init_array(&points, near_circle, 3, 2, IGRAPH_ROW_MAJOR);
        gabriel_identity(&points, &graph);
        igraph_destroy(&graph);
        igraph_matrix_destroy(&points);
    }

    {
        /* Two points: a single edge that persists for all beta. */
        const igraph_real_t two[] = { 0, 0,  1, 0 };

        igraph_matrix_init_array(&points, two, 2, 2, IGRAPH_ROW_MAJOR);
        igraph_beta_weighted_gabriel_graph(&graph, &weights, &points, IGRAPH_INFINITY);
        IGRAPH_ASSERT(igraph_ecount(&graph) == 1);
        IGRAPH_ASSERT(VECTOR(weights)[0] == IGRAPH_INFINITY);
        igraph_destroy(&graph);
        igraph_matrix_destroy(&points);
    }

    {
        /* Random cells of a 30 arc-second grid in longitude and latitude. */
        igraph_vector_int_t cells;

        igraph_rng_seed(igraph_rng_default(), 42);
        igraph_vector_int_init_range(&cells, 0, 100);
        igraph_vector_int_shuffle(&cells);
        igraph_matrix_init(&points, 30, 2);
        for (igraph_int_t i = 0; i < 30; i++) {
            MATRIX(points, i, 0) = 19.05 + (VECTOR(cells)[i] % 10) / 120.0;
            MATRIX(points, i, 1) = 47.5 + (VECTOR(cells)[i] / 10) / 120.0;
        }
        gabriel_identity(&points, &graph);
        igraph_destroy(&graph);
        igraph_matrix_destroy(&points);
        igraph_vector_int_destroy(&cells);
    }

    {
        /* Triangular lattices far from the origin. The rounded coordinates are no
         * longer exactly equilateral, but are treated as such: the relative
         * neighborhood graph keeps the triangles, the beta = 2 skeleton has no edges. */
        const igraph_real_t offsets[][3] = {
            /* x, y, spacing */
            { 0, 0, 1 }, { 47.5, 19.05, 1e-3 }, { 47.5, 19.05, 1e-4 }, { 5e6, 3e5, 10 }
        };
        const igraph_int_t k = 6;

        igraph_matrix_init(&points, k * k, 2);
        for (size_t c = 0; c < sizeof(offsets) / sizeof(offsets[0]); c++) {
            const igraph_real_t h = offsets[c][2];
            for (igraph_int_t i = 0; i < k; i++) {
                for (igraph_int_t j = 0; j < k; j++) {
                    MATRIX(points, i * k + j, 0) = offsets[c][0] + (i + 0.5 * (j % 2)) * h;
                    MATRIX(points, i * k + j, 1) = offsets[c][1] + j * h * sqrt(3.0) / 2;
                }
            }

            igraph_delaunay_graph(&graph, &points);
            IGRAPH_ASSERT(igraph_ecount(&graph) == 89);
            igraph_destroy(&graph);

            gabriel_identity(&points, &graph);
            IGRAPH_ASSERT(igraph_ecount(&graph) == 85);
            igraph_destroy(&graph);

            igraph_relative_neighborhood_graph(&graph, &points);
            IGRAPH_ASSERT(igraph_ecount(&graph) == 85);
            igraph_destroy(&graph);

            igraph_lune_beta_skeleton(&graph, &points, 2);
            IGRAPH_ASSERT(igraph_ecount(&graph) == 0);
            igraph_destroy(&graph);
        }
        igraph_matrix_destroy(&points);
    }

    {
        /* Random cells of a grid, in longitude and latitude and in integer
         * coordinates. The skeletons must not depend on the coordinates. */
        const igraph_real_t betas[] = { 1, 1.5, 2, 3 };
        igraph_vector_int_t cells;
        igraph_matrix_t int_points;
        igraph_t graph2;

        igraph_rng_seed(igraph_rng_default(), 137);
        igraph_vector_int_init_range(&cells, 0, 400);
        igraph_vector_int_shuffle(&cells);
        igraph_matrix_init(&points, 100, 2);
        igraph_matrix_init(&int_points, 100, 2);
        for (igraph_int_t i = 0; i < 100; i++) {
            MATRIX(int_points, i, 0) = VECTOR(cells)[i] % 20;
            MATRIX(int_points, i, 1) = VECTOR(cells)[i] / 20;
            MATRIX(points, i, 0) = 19.05 + MATRIX(int_points, i, 0) / 120;
            MATRIX(points, i, 1) = 47.5 + MATRIX(int_points, i, 1) / 120;
        }

        for (size_t i = 0; i < sizeof(betas) / sizeof(betas[0]); i++) {
            igraph_lune_beta_skeleton(&graph, &points, betas[i]);
            igraph_lune_beta_skeleton(&graph2, &int_points, betas[i]);
            check_same_edges(&graph, &graph2);
            igraph_destroy(&graph2);
            igraph_destroy(&graph);

            igraph_circle_beta_skeleton(&graph, &points, betas[i]);
            igraph_circle_beta_skeleton(&graph2, &int_points, betas[i]);
            check_same_edges(&graph, &graph2);
            igraph_destroy(&graph2);
            igraph_destroy(&graph);
        }

        igraph_relative_neighborhood_graph(&graph, &points);
        igraph_relative_neighborhood_graph(&graph2, &int_points);
        check_same_edges(&graph, &graph2);
        igraph_destroy(&graph2);
        igraph_destroy(&graph);

        igraph_matrix_destroy(&int_points);
        igraph_matrix_destroy(&points);
        igraph_vector_int_destroy(&cells);
    }

    {
        /* A short baseline 0-2 and a distant point 1, at a right angle at point 0.
         * Point 2 is within the tolerance of the Gabriel circle of 0-1, but it is
         * farther from 1 than 0 is, so it must not remove 0-1, an edge of the
         * minimum spanning tree. */
        const igraph_real_t mm_km[] = { 1e6, 1e6,  1e6 + 1000, 1e6,  1e6, 1e6 + 0.001 };
        const igraph_real_t cm_100km[] = { 5e6, 5e6,  5.1e6, 5e6,  5e6, 5e6 + 0.01 };
        const igraph_real_t *cases[] = { mm_km, cm_100km };
        const igraph_real_t cutoffs[] = { 1, 1.99, 2, IGRAPH_INFINITY };
        igraph_int_t eid;

        for (size_t c = 0; c < sizeof(cases) / sizeof(cases[0]); c++) {
            igraph_matrix_init_array(&points, cases[c], 3, 2, IGRAPH_ROW_MAJOR);
            gabriel_identity(&points, &graph);
            check_connected(&graph);
            igraph_destroy(&graph);
            check_nesting(&points);

            /* Point 2 gives no finite threshold for 0-1. Point 0 is on the Gabriel
             * circle of 1-2, but it is not closer to 1 than 2 is, up to the
             * tolerance, so its threshold for 1-2 is 2, not 1. */
            for (size_t i = 0; i < sizeof(cutoffs) / sizeof(cutoffs[0]); i++) {
                igraph_beta_weighted_gabriel_graph(&graph, &weights, &points, cutoffs[i]);
                igraph_get_eid(&graph, &eid, 0, 1, IGRAPH_UNDIRECTED, true);
                IGRAPH_ASSERT(VECTOR(weights)[eid] == IGRAPH_INFINITY);
                igraph_get_eid(&graph, &eid, 1, 2, IGRAPH_UNDIRECTED, true);
                IGRAPH_ASSERT(VECTOR(weights)[eid] == (cutoffs[i] < 2 ? IGRAPH_INFINITY : 2));
                igraph_destroy(&graph);
            }
            igraph_matrix_destroy(&points);
        }
    }

    {
        /* As above, with point 2 slightly farther from the Gabriel circle. The
         * weight of 0-1 is infinite, and the skeletons for beta < 2 keep 0-1. */
        const igraph_real_t far[] = { 1e6, 1e6,  1e6 + 1000, 1e6,  1e6, 1e6 + 0.0013 };
        igraph_bool_t adj;

        igraph_matrix_init_array(&points, far, 3, 2, IGRAPH_ROW_MAJOR);
        igraph_lune_beta_skeleton(&graph, &points, 1.5);
        igraph_are_adjacent(&graph, 0, 1, &adj);
        IGRAPH_ASSERT(adj);
        igraph_destroy(&graph);
        check_nesting(&points);
        igraph_matrix_destroy(&points);
    }

    {
        /* Edge 0-1 with a point close to endpoint 0, outside the regions of 0-1
         * for small beta. The tolerance of these regions must not grow as beta
         * decreases: otherwise a smaller beta removes 0-1 while a larger beta keeps it. */
        const igraph_real_t cases[][6] = {
            { 1e6, 1e6,  1e6 + 1, 1e6,  1e6 + 2.5e-9, 1e6 + 5e-9 },
            { 1e6, 1e6,  1e6 + 1, 1e6,  1e6 + 1.25e-5, 1e6 + 5e-7 },
            { 1e6, 1e6,  1e6 + 1e-5, 1e6,  1e6 + 0.5e-5, 1e6 + 0.6e-5 },
            { 0, 0,  1, 0,  1e-9, 1.2e-11 }
        };
        const skeleton_t types[] = { LUNE, LUNE, LUNE, LUNE, LUNE, LUNE, LUNE, LUNE, LUNE };
        const igraph_real_t betas[] = { 1e-4, 1e-3, 2e-3, 0.01, 0.02, 0.1, 0.2, 0.5, 1 };
        igraph_bool_t adj;

        for (size_t c = 0; c < sizeof(cases) / sizeof(cases[0]); c++) {
            igraph_matrix_init_array(&points, cases[c], 3, 2, IGRAPH_ROW_MAJOR);
            check_sequence(&points, types, betas, sizeof(betas) / sizeof(betas[0]), 2);
            igraph_lune_beta_skeleton(&graph, &points, 1e-4);
            igraph_are_adjacent(&graph, 0, 1, &adj);
            IGRAPH_ASSERT(adj);
            igraph_destroy(&graph);
            igraph_matrix_destroy(&points);
        }
    }

    {
        /* A thin rectangle in 4D, which uses all pairs of points as candidates. */
        const igraph_real_t rect[] = {
            0, 0, 0, 0,
            1, 0, 0, 0,
            1, 1.5e-7, 0, 0,
            0, 1.5e-7, 0, 0
        };

        igraph_matrix_init_array(&points, rect, 4, 4, IGRAPH_ROW_MAJOR);
        check_nesting(&points);
        igraph_matrix_destroy(&points);
    }

    {
        /* A grid in projected coordinates with a point 1 mm north of one site.
         * That point removes the edge to the next site to the north, but no
         * other edge between grid sites. */
        const igraph_int_t k = 5, centre = 2 * 5 + 2, extra = 5 * 5;
        igraph_int_t grid_edges = 0;
        igraph_bool_t adj;

        igraph_matrix_init(&points, k * k + 1, 2);
        for (igraph_int_t i = 0; i < k * k; i++) {
            MATRIX(points, i, 0) = 5e5 + 1000 * (i / k);
            MATRIX(points, i, 1) = 5e6 + 1000 * (i % k);
        }
        MATRIX(points, extra, 0) = MATRIX(points, centre, 0);
        MATRIX(points, extra, 1) = MATRIX(points, centre, 1) + 0.001;

        igraph_gabriel_graph(&graph, &points);
        for (igraph_int_t i = 0; i < k * k; i++) {
            if (i % k + 1 < k) {
                igraph_are_adjacent(&graph, i, i + 1, &adj);
                IGRAPH_ASSERT(adj == (i != centre));
            }
            if (i + k < k * k) {
                igraph_are_adjacent(&graph, i, i + k, &adj);
                IGRAPH_ASSERT(adj);
            }
        }
        for (igraph_int_t e = 0; e < igraph_ecount(&graph); e++) {
            if (IGRAPH_FROM(&graph, e) != extra && IGRAPH_TO(&graph, e) != extra) {
                grid_edges++;
            }
        }
        IGRAPH_ASSERT(grid_edges == 2 * k * (k - 1) - 1);
        igraph_are_adjacent(&graph, extra, centre, &adj);
        IGRAPH_ASSERT(adj);
        igraph_are_adjacent(&graph, extra, centre + 1, &adj);
        IGRAPH_ASSERT(adj);
        igraph_destroy(&graph);

        check_nesting(&points);
        igraph_matrix_destroy(&points);
    }

    {
        /* Random points far from the origin, in 2D and 3D. */
        const igraph_real_t offsets[][3] = { { 0, 0, 0 }, { 19.05, 47.5, 0 }, { 5e5, 5e6, 100 } };
        const igraph_real_t scales[] = { 1, 0.01, 1000 };

        igraph_rng_seed(igraph_rng_default(), 4242);
        for (igraph_int_t dim = 2; dim <= 3; dim++) {
            for (size_t c = 0; c < sizeof(scales) / sizeof(scales[0]); c++) {
                igraph_matrix_init(&points, 60, dim);
                for (igraph_int_t i = 0; i < 60; i++) {
                    for (igraph_int_t k = 0; k < dim; k++) {
                        MATRIX(points, i, k) = offsets[c][k] + scales[c] * RNG_UNIF01();
                    }
                }
                check_nesting(&points);
                igraph_matrix_destroy(&points);
            }
        }
    }

    {
        /* Rotated lattices far from the origin. The Gabriel graph of a square or
         * cubic lattice has no diagonals. The relative neighborhood graph of a
         * triangular lattice keeps all triangles, the beta = 2 skeleton has no edges. */
        const igraph_int_t k = 6;
        const igraph_real_t c = cos(0.3), s = sin(0.3), c2 = cos(0.7), s2 = sin(0.7);

        /* Square lattice in projected coordinates. */
        igraph_matrix_init(&points, k * k, 2);
        for (igraph_int_t i = 0; i < k * k; i++) {
            const igraph_real_t x = 30.0 * (i / k), y = 30.0 * (i % k);
            MATRIX(points, i, 0) = 5e5 + c * x - s * y;
            MATRIX(points, i, 1) = 5e6 + s * x + c * y;
        }
        gabriel_identity(&points, &graph);
        check_unit_edges(&graph, &points, 30, 2 * k * (k - 1));
        igraph_destroy(&graph);
        check_nesting(&points);

        /* Triangular lattice in geographic coordinates. */
        for (igraph_int_t i = 0; i < k * k; i++) {
            const igraph_real_t x = (i / k + 0.5 * (i % k)) / 120, y = (i % k) * sqrt(3.0) / 2 / 120;
            MATRIX(points, i, 0) = 19.05 + c * x - s * y;
            MATRIX(points, i, 1) = 47.5 + s * x + c * y;
        }
        igraph_relative_neighborhood_graph(&graph, &points);
        check_unit_edges(&graph, &points, 1.0 / 120, 3 * k * k - 4 * k + 1);
        igraph_destroy(&graph);
        igraph_lune_beta_skeleton(&graph, &points, 1.99);
        check_unit_edges(&graph, &points, 1.0 / 120, 3 * k * k - 4 * k + 1);
        igraph_destroy(&graph);
        igraph_lune_beta_skeleton(&graph, &points, 2);
        IGRAPH_ASSERT(igraph_ecount(&graph) == 0);
        igraph_destroy(&graph);
        check_nesting(&points);
        igraph_matrix_destroy(&points);

        /* Cubic lattice in projected coordinates. */
        igraph_matrix_init(&points, 4 * 4 * 4, 3);
        for (igraph_int_t i = 0; i < 4 * 4 * 4; i++) {
            const igraph_real_t x = i / 16, y = (i / 4) % 4, z = i % 4;
            const igraph_real_t x1 = c * x - s * y, y1 = s * x + c * y;
            MATRIX(points, i, 0) = 5e5 + x1;
            MATRIX(points, i, 1) = 5e6 + c2 * y1 - s2 * z;
            MATRIX(points, i, 2) = 100 + s2 * y1 + c2 * z;
        }
        gabriel_identity(&points, &graph);
        check_unit_edges(&graph, &points, 1, 3 * 4 * 4 * 3);
        igraph_destroy(&graph);
        check_nesting(&points);
        igraph_matrix_destroy(&points);
    }

    {
        /* Circle-based skeletons. An equilateral triangle loses all its edges
         * exactly at beta = 2 / sqrt(3), where its points are on the boundaries
         * of the regions. Above beta = 1, a short perpendicular leg 0-2 removes
         * the long edge 0-1, and point 0 removes 1-2, so the graph is not connected. */
        const igraph_real_t equilateral[] = { 0, 0,  1, 0,  0.5, 0.8660254037844386 };
        const igraph_real_t right[] = { 0, 0,  1, 0,  0, 0.1 };
        const igraph_int_t short_leg[] = { 0, 2 };
        const igraph_int_t legs[] = { 0, 1,  0, 2 };

        igraph_matrix_init_array(&points, equilateral, 3, 2, IGRAPH_ROW_MAJOR);
        igraph_circle_beta_skeleton(&graph, &points, 1.1);
        IGRAPH_ASSERT(igraph_ecount(&graph) == 3);
        igraph_destroy(&graph);
        igraph_circle_beta_skeleton(&graph, &points, 2 / sqrt(3.0));
        IGRAPH_ASSERT(igraph_ecount(&graph) == 0);
        igraph_destroy(&graph);
        igraph_circle_beta_skeleton(&graph, &points, 1.2);
        IGRAPH_ASSERT(igraph_ecount(&graph) == 0);
        igraph_destroy(&graph);
        igraph_matrix_destroy(&points);

        igraph_matrix_init_array(&points, right, 3, 2, IGRAPH_ROW_MAJOR);
        igraph_circle_beta_skeleton(&graph, &points, 1.01);
        check_edges(&graph, short_leg, 1);
        igraph_destroy(&graph);
        igraph_lune_beta_skeleton(&graph, &points, 1.01);
        check_edges(&graph, legs, 2);
        igraph_destroy(&graph);
        check_nesting(&points);
        igraph_matrix_destroy(&points);
    }

    igraph_vector_destroy(&weights);

    VERIFY_FINALLY_STACK();
    return 0;
}
