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
 * to the distances between the points. */

/* Checks that the graph has exactly the given edges. */
static void check_edges(const igraph_t *graph, const igraph_int_t *edges, igraph_int_t count) {
    IGRAPH_ASSERT(igraph_ecount(graph) == count);
    for (igraph_int_t i = 0; i < count; i++) {
        igraph_bool_t adj;
        igraph_are_adjacent(graph, edges[2 * i], edges[2 * i + 1], &adj);
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

    igraph_vector_destroy(&weights);

    VERIFY_FINALLY_STACK();
    return 0;
}
