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
#include <math.h>
#include <stdio.h>

#include "test_utilities.h"
#include "core/barnes_hut.h"

/* --- Force Callbacks for Testing --- */

/* Simple repulsive force: Inverse square law F = m1*m2 / d^2 */
static void test_repulsive_force(
    const igraph_bh_point_t *p1,
    const igraph_bh_point_t *p2,
    igraph_real_t *force,
    void *user_data
) {
    IGRAPH_UNUSED(user_data);
    igraph_real_t dx = p1->coord[0] - p2->coord[0];
    igraph_real_t dy = p1->coord[1] - p2->coord[1];
    igraph_real_t dz = p1->coord[2] - p2->coord[2]; // Will be 0.0 in 2D

    igraph_real_t dist_sq = dx*dx + dy*dy + dz*dz;
    if (dist_sq < 1e-6) return; // Avoid division by zero

    igraph_real_t dist = sqrt(dist_sq);
    igraph_real_t f = (p1->mass * p2->mass) / dist_sq;

    force[0] += f * (dx / dist);
    force[1] += f * (dy / dist);
    force[2] += f * (dz / dist);
}

/* Simple attractive force: Hooke's Law (springs) F = k * d */
static void test_attractive_force(
    const igraph_bh_point_t *p1,
    const igraph_bh_point_t *p2,
    igraph_real_t *force,
    void *user_data
) {
    IGRAPH_UNUSED(user_data);
    // Vector pointing from p1 to p2
    igraph_real_t dx = p2->coord[0] - p1->coord[0];
    igraph_real_t dy = p2->coord[1] - p1->coord[1];
    igraph_real_t dz = p2->coord[2] - p1->coord[2];

    // F = 1.0 * distance vector
    force[0] += dx;
    force[1] += dy;
    force[2] += dz;
}


/* --- Unit Tests --- */

int test_init_and_destroy(void) {
    igraph_bh_tree_t tree;

    /* Initialize and verify defaults/parameters */
    IGRAPH_ASSERT(igraph_bh_tree_init(&tree, 2, 0.5, 15, 2) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(tree.dim == 2);
    IGRAPH_ASSERT(tree.bh_theta == 0.5);
    IGRAPH_ASSERT(tree.max_level == 15);
    IGRAPH_ASSERT(tree.leaf_capacity == 2);
    IGRAPH_ASSERT(tree.nodes == NULL);
    IGRAPH_ASSERT(tree.points == NULL);
    IGRAPH_ASSERT(tree.node_count == 0);

    igraph_bh_tree_destroy(&tree);

    return 0;
}

int test_empty_tree(void) {
    igraph_bh_tree_t tree;
    igraph_matrix_t coords;

    igraph_matrix_init(&coords, 0, 2); /* Empty matrix */

    IGRAPH_ASSERT(igraph_bh_tree_init(&tree, 2, 0.6, 10, 1) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_bh_tree_build(&tree, &coords, NULL) == IGRAPH_SUCCESS);

    IGRAPH_ASSERT(tree.point_count == 0);
    IGRAPH_ASSERT(tree.node_count == 0);

    igraph_bh_tree_destroy(&tree);
    igraph_matrix_destroy(&coords);
    return 0;
}

int test_2d_tree_build(void) {
    igraph_bh_tree_t tree;
    igraph_matrix_t coords;

    /* Place 4 points nicely in 4 separate quadrants to force a clean split */
    igraph_matrix_init(&coords, 4, 2);
    MATRIX(coords, 0, 0) = 0.0; MATRIX(coords, 0, 1) = 0.0;
    MATRIX(coords, 1, 0) = 10.0; MATRIX(coords, 1, 1) = 0.0;
    MATRIX(coords, 2, 0) = 0.0; MATRIX(coords, 2, 1) = 10.0;
    MATRIX(coords, 3, 0) = 10.0; MATRIX(coords, 3, 1) = 10.0;

    IGRAPH_ASSERT(igraph_bh_tree_init(&tree, 2, 0.5, 10, 1) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_bh_tree_build(&tree, &coords, NULL) == IGRAPH_SUCCESS);

    IGRAPH_ASSERT(tree.point_count == 4);
    IGRAPH_ASSERT(tree.node_count > 1); /* Should have root + 4 children */

    /* Verify Root Node stats */
    igraph_bh_node_t *root = &tree.nodes[0];
    IGRAPH_ASSERT(root->mass == 4.0); /* Default mass is 1.0 per point */
    IGRAPH_ASSERT(root->center[0] == 5.0); /* CoM should be strictly center */
    IGRAPH_ASSERT(root->center[1] == 5.0);
    IGRAPH_ASSERT(root->is_leaf == 0);
    IGRAPH_ASSERT(root->point_count == 4);

    igraph_bh_tree_destroy(&tree);
    igraph_matrix_destroy(&coords);
    return 0;
}

int test_3d_coincident_points(void) {
    igraph_bh_tree_t tree;
    igraph_matrix_t coords;

    /* Test identical coordinates. This tests if `max_level` gracefully stops infinite subdivision */
    igraph_matrix_init(&coords, 3, 3);
    for (igraph_integer_t i = 0; i < 3; i++) {
        MATRIX(coords, i, 0) = 5.0;
        MATRIX(coords, i, 1) = -2.0;
        MATRIX(coords, i, 2) = 1.0;
    }

    IGRAPH_ASSERT(igraph_bh_tree_init(&tree, 3, 0.5, 5, 1) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_bh_tree_build(&tree, &coords, NULL) == IGRAPH_SUCCESS);

    /* Should reach max_level 5, so depth shouldn't explode */
    IGRAPH_ASSERT(tree.node_count > 1);
    IGRAPH_ASSERT(tree.nodes[0].mass == 3.0);

    igraph_bh_tree_destroy(&tree);
    igraph_matrix_destroy(&coords);
    return 0;
}

int test_repulsive_force_calculation(void) {
    igraph_bh_tree_t tree;
    igraph_matrix_t coords;
    igraph_matrix_t forces;

    igraph_matrix_init(&coords, 2, 2);
    /* Two points separated by 2.0 units on X-axis */
    MATRIX(coords, 0, 0) = -1.0; MATRIX(coords, 0, 1) = 0.0;
    MATRIX(coords, 1, 0) = 1.0;  MATRIX(coords, 1, 1) = 0.0;

    igraph_matrix_init(&forces, 2, 2);
    igraph_matrix_fill(&forces, 0.0);

    IGRAPH_ASSERT(igraph_bh_tree_init(&tree, 2, 0.6, 10, 1) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_bh_tree_build(&tree, &coords, NULL) == IGRAPH_SUCCESS);

    IGRAPH_ASSERT(igraph_bh_calculate_repulsive_forces(&tree, &forces, test_repulsive_force, NULL) == IGRAPH_SUCCESS);

    /* * Distance is 2.0. Inverse square law gives f = 1 / (2^2) = 0.25
     * Point 0 is at -1.0, so it gets pushed away in the negative X direction.
     * Point 1 is at 1.0, so it gets pushed in the positive X direction.
     */
    IGRAPH_ASSERT(MATRIX(forces, 0, 0) < -0.24 && MATRIX(forces, 0, 0) > -0.26);
    IGRAPH_ASSERT(MATRIX(forces, 1, 0) > 0.24 && MATRIX(forces, 1, 0) < 0.26);

    /* Y forces should be zero */
    IGRAPH_ASSERT(MATRIX(forces, 0, 1) == 0.0);
    IGRAPH_ASSERT(MATRIX(forces, 1, 1) == 0.0);

    igraph_bh_tree_destroy(&tree);
    igraph_matrix_destroy(&coords);
    igraph_matrix_destroy(&forces);
    return 0;
}

int test_attractive_force_calculation(void) {
    igraph_bh_tree_t tree;
    igraph_matrix_t coords;
    igraph_matrix_t forces;
    igraph_vector_int_t from, to;
    igraph_vector_t weights;

    igraph_matrix_init(&coords, 2, 2);
    /* Two points separated by 5.0 units on Y-axis */
    MATRIX(coords, 0, 0) = 0.0; MATRIX(coords, 0, 1) = 0.0;
    MATRIX(coords, 1, 0) = 0.0; MATRIX(coords, 1, 1) = 5.0;

    igraph_matrix_init(&forces, 2, 2);
    igraph_matrix_fill(&forces, 0.0);

    igraph_vector_int_init(&from, 1);
    igraph_vector_int_init(&to, 1);
    VECTOR(from)[0] = 0;
    VECTOR(to)[0] = 1;

    igraph_vector_init(&weights, 1);
    VECTOR(weights)[0] = 2.0; /* Edge weight multipliers */

    IGRAPH_ASSERT(igraph_bh_tree_init(&tree, 2, 0.6, 10, 1) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_bh_tree_build(&tree, &coords, NULL) == IGRAPH_SUCCESS);

    IGRAPH_ASSERT(igraph_bh_calculate_attractive_forces(
        &tree, &from, &to, &weights, &forces, test_attractive_force, NULL
    ) == IGRAPH_SUCCESS);

    /*
     * Point 0 (at 0) is pulled to Point 1 (at 5) with dist 5. F = 5. Weight = 2. Total = +10 on Y.
     * Point 1 is pulled to Point 0. Total = -10 on Y.
     */
    IGRAPH_ASSERT(MATRIX(forces, 0, 1) == 10.0);
    IGRAPH_ASSERT(MATRIX(forces, 1, 1) == -10.0);

    /* Clean up */
    igraph_bh_tree_destroy(&tree);
    igraph_matrix_destroy(&coords);
    igraph_matrix_destroy(&forces);
    igraph_vector_int_destroy(&from);
    igraph_vector_int_destroy(&to);
    igraph_vector_destroy(&weights);
    return 0;
}


int main(void) {
    igraph_set_error_handler(igraph_error_handler_ignore);

    RUN_TEST(test_init_and_destroy);
    RUN_TEST(test_empty_tree);
    RUN_TEST(test_2d_tree_build);
    RUN_TEST(test_3d_coincident_points);
    RUN_TEST(test_repulsive_force_calculation);
    RUN_TEST(test_attractive_force_calculation);

    return 0;
}
