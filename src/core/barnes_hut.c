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

#include "barnes_hut.h"
#include "igraph_memory.h"
#include "igraph_error.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#define IGRAPH_BH_DEFAULT_THETA 0.6
#define IGRAPH_BH_DEFAULT_MAX_LEVEL 10
#define IGRAPH_BH_DEFAULT_LEAF_CAPACITY 1
#define IGRAPH_BH_CHILDREN_COUNT(dim) (1 << (dim))

/* Helper to dynamically grow the contiguous node array */
static igraph_error_t ensure_node_capacity(igraph_bh_tree_t *tree, igraph_integer_t needed) {
    if (needed <= tree->capacity) {
        return IGRAPH_SUCCESS;
    }

    igraph_integer_t new_cap = tree->capacity == 0 ? 1024 : tree->capacity * 2;
    while (new_cap < needed) {
        new_cap *= 2;
    }

    igraph_bh_node_t *new_nodes = IGRAPH_REALLOC(tree->nodes, new_cap, igraph_bh_node_t);
    if (!new_nodes) {
        IGRAPH_ERROR("Failed to allocate memory for Barnes-Hut tree nodes.", IGRAPH_ENOMEM);
    }

    tree->nodes = new_nodes;
    tree->capacity = new_cap;
    return IGRAPH_SUCCESS;
}

static igraph_error_t build_bh_node_recursive(
    igraph_bh_tree_t *tree,
    igraph_integer_t node_idx,
    igraph_integer_t start_idx,
    igraph_integer_t count,
    igraph_real_t center[3],
    igraph_real_t size,
    igraph_integer_t level
) {
    tree->nodes[node_idx].size = size;
    tree->nodes[node_idx].point_count = count;

    if (count == 0) {
        tree->nodes[node_idx].is_leaf = 1;
        tree->nodes[node_idx].data.first_point_idx = start_idx;
        tree->nodes[node_idx].mass = 0;
        tree->nodes[node_idx].center[0] = center[0];
        tree->nodes[node_idx].center[1] = center[1];
        tree->nodes[node_idx].center[2] = center[2];
        return IGRAPH_SUCCESS;
    }

    igraph_real_t total_mass = 0.0;
    igraph_real_t cx = 0.0, cy = 0.0, cz = 0.0;

    for (igraph_integer_t i = 0; i < count; i++) {
        igraph_integer_t pid = tree->point_indices[start_idx + i];
        igraph_bh_point_t *p = &tree->points[pid];
        igraph_real_t m = p->mass > 0 ? p->mass : 1.0;

        total_mass += m;
        cx += p->coord[0] * m;
        cy += p->coord[1] * m;
        if (tree->dim == 3) {
            cz += p->coord[2] * m;
        }
    }

    tree->nodes[node_idx].mass = total_mass;
    if (total_mass > 0) {
        tree->nodes[node_idx].center[0] = cx / total_mass;
        tree->nodes[node_idx].center[1] = cy / total_mass;
        tree->nodes[node_idx].center[2] = (tree->dim == 3) ? (cz / total_mass) : 0.0;
    } else {
        tree->nodes[node_idx].center[0] = center[0];
        tree->nodes[node_idx].center[1] = center[1];
        tree->nodes[node_idx].center[2] = center[2];
    }

    /* Base case: node is a leaf if under capacity or reached max tree depth */
    if (count <= tree->leaf_capacity || level >= tree->max_level) {
        tree->nodes[node_idx].is_leaf = 1;
        tree->nodes[node_idx].data.first_point_idx = start_idx;
        return IGRAPH_SUCCESS;
    }

    tree->nodes[node_idx].is_leaf = 0;

    int n_children = IGRAPH_BH_CHILDREN_COUNT(tree->dim);
    igraph_integer_t *counts = IGRAPH_CALLOC(n_children, igraph_integer_t);
    if (!counts) {
        IGRAPH_ERROR("Memory allocation failed for Barnes-Hut counts.", IGRAPH_ENOMEM);
    }

    /* Assign points to quadrants based on strict geometric center, NOT center of mass */
    for (igraph_integer_t i = 0; i < count; i++) {
        igraph_integer_t pid = tree->point_indices[start_idx + i];
        igraph_bh_point_t *p = &tree->points[pid];
        int child_idx = 0;
        if (p->coord[0] >= center[0]) child_idx |= 1;
        if (p->coord[1] >= center[1]) child_idx |= 2;
        if (tree->dim == 3 && p->coord[2] >= center[2]) child_idx |= 4;
        counts[child_idx]++;
    }

    igraph_integer_t *offsets = IGRAPH_CALLOC(n_children + 1, igraph_integer_t);
    if (!offsets) {
        IGRAPH_FREE(counts);
        IGRAPH_ERROR("Memory allocation failed for Barnes-Hut offsets.", IGRAPH_ENOMEM);
    }
    for (int i = 1; i <= n_children; i++) {
        offsets[i] = offsets[i-1] + counts[i-1];
    }

    igraph_integer_t *temp = IGRAPH_MALLOC(count * sizeof(igraph_integer_t));
    if (!temp) {
        IGRAPH_FREE(offsets);
        IGRAPH_FREE(counts);
        IGRAPH_ERROR("Memory allocation failed for Barnes-Hut temp array.", IGRAPH_ENOMEM);
    }

    igraph_integer_t *cursors = IGRAPH_CALLOC(n_children, igraph_integer_t);
    if (!cursors) {
        IGRAPH_FREE(temp);
        IGRAPH_FREE(offsets);
        IGRAPH_FREE(counts);
        IGRAPH_ERROR("Memory allocation failed for Barnes-Hut cursors.", IGRAPH_ENOMEM);
    }
    for (int i = 0; i < n_children; i++) {
        cursors[i] = offsets[i];
    }

    /* Permute point indices so children's points sit contiguously in memory */
    for (igraph_integer_t i = 0; i < count; i++) {
        igraph_integer_t pid = tree->point_indices[start_idx + i];
        igraph_bh_point_t *p = &tree->points[pid];
        int child_idx = 0;
        if (p->coord[0] >= center[0]) child_idx |= 1;
        if (p->coord[1] >= center[1]) child_idx |= 2;
        if (tree->dim == 3 && p->coord[2] >= center[2]) child_idx |= 4;
        temp[cursors[child_idx]++] = pid;
    }

    for (igraph_integer_t i = 0; i < count; i++) {
        tree->point_indices[start_idx + i] = temp[i];
    }

    IGRAPH_FREE(cursors);
    IGRAPH_FREE(temp);

    /* Allocate all children contiguously at the end of the current node array */
    IGRAPH_CHECK(ensure_node_capacity(tree, tree->node_count + n_children));
    igraph_integer_t first_child_idx = tree->node_count;
    tree->node_count += n_children;

    /* Re-fetch the node pointer in case ensure_node_capacity invoked a realloc */
    tree->nodes[node_idx].data.first_child_idx = first_child_idx;

    igraph_real_t half_size = size / 2.0;
    for (int i = 0; i < n_children; i++) {
        igraph_real_t child_center[3];
        child_center[0] = center[0] + ((i & 1) ? half_size : -half_size);
        child_center[1] = center[1] + ((i & 2) ? half_size : -half_size);
        child_center[2] = (tree->dim == 3) ? (center[2] + ((i & 4) ? half_size : -half_size)) : 0.0;

        igraph_error_t err = build_bh_node_recursive(
            tree, first_child_idx + i,
            start_idx + offsets[i], counts[i],
            child_center, half_size, level + 1
        );

        if (err != IGRAPH_SUCCESS) {
            IGRAPH_FREE(offsets);
            IGRAPH_FREE(counts);
            return err;
        }
    }

    IGRAPH_FREE(offsets);
    IGRAPH_FREE(counts);

    return IGRAPH_SUCCESS;
}


igraph_error_t igraph_bh_tree_init(
    igraph_bh_tree_t *tree,
    igraph_integer_t dim,
    igraph_real_t theta,
    igraph_integer_t max_level,
    igraph_integer_t leaf_capacity
) {
    tree->nodes = NULL;
    tree->points = NULL;
    tree->point_indices = NULL;
    tree->node_count = 0;
    tree->capacity = 0;
    tree->point_count = 0;

    tree->dim = dim;
    tree->bh_theta = theta > 0 ? theta : IGRAPH_BH_DEFAULT_THETA;
    tree->max_level = max_level > 0 ? max_level : IGRAPH_BH_DEFAULT_MAX_LEVEL;
    tree->leaf_capacity = leaf_capacity > 0 ? leaf_capacity : IGRAPH_BH_DEFAULT_LEAF_CAPACITY;

    return IGRAPH_SUCCESS;
}

void igraph_bh_tree_destroy(igraph_bh_tree_t *tree) {
    if (tree->points) {
        IGRAPH_FREE(tree->points);
    }
    if (tree->nodes) {
        IGRAPH_FREE(tree->nodes);
    }
    if (tree->point_indices) {
        IGRAPH_FREE(tree->point_indices);
    }

    tree->points = NULL;
    tree->nodes = NULL;
    tree->point_indices = NULL;
    tree->node_count = 0;
    tree->point_count = 0;
    tree->capacity = 0;
}

igraph_error_t igraph_bh_tree_build(
    igraph_bh_tree_t *tree,
    const igraph_matrix_t *coords,
    const igraph_vector_t *masses
) {
    /* Cleanup existing data if reusing the structure */
    if (tree->points) IGRAPH_FREE(tree->points);
    if (tree->point_indices) IGRAPH_FREE(tree->point_indices);
    tree->node_count = 0;

    igraph_integer_t n = igraph_matrix_nrow(coords);

    if (n == 0) {
        tree->point_count = 0;
        return IGRAPH_SUCCESS;
    }

    tree->point_count = n;

    tree->points = IGRAPH_CALLOC(n, igraph_bh_point_t);
    if (!tree->points) return IGRAPH_ENOMEM;

    tree->point_indices = IGRAPH_MALLOC(n * sizeof(igraph_integer_t));
    if (!tree->point_indices) {
        IGRAPH_FREE(tree->points);
        tree->points = NULL;
        return IGRAPH_ENOMEM;
    }

    igraph_real_t min_x = MATRIX(*coords, 0, 0);
    igraph_real_t max_x = min_x;
    igraph_real_t min_y = MATRIX(*coords, 0, 1);
    igraph_real_t max_y = min_y;
    igraph_real_t min_z = (tree->dim == 3) ? MATRIX(*coords, 0, 2) : 0.0;
    igraph_real_t max_z = min_z;

    for (igraph_integer_t i = 0; i < n; i++) {
        tree->point_indices[i] = i;

        tree->points[i].id = i;
        tree->points[i].coord[0] = MATRIX(*coords, i, 0);
        tree->points[i].coord[1] = MATRIX(*coords, i, 1);
        tree->points[i].coord[2] = (tree->dim == 3) ? MATRIX(*coords, i, 2) : 0.0;
        tree->points[i].mass = masses ? VECTOR(*masses)[i] : 1.0;
        tree->points[i].data = NULL;

        if (tree->points[i].coord[0] < min_x) min_x = tree->points[i].coord[0];
        if (tree->points[i].coord[0] > max_x) max_x = tree->points[i].coord[0];
        if (tree->points[i].coord[1] < min_y) min_y = tree->points[i].coord[1];
        if (tree->points[i].coord[1] > max_y) max_y = tree->points[i].coord[1];
        if (tree->dim == 3) {
            if (tree->points[i].coord[2] < min_z) min_z = tree->points[i].coord[2];
            if (tree->points[i].coord[2] > max_z) max_z = tree->points[i].coord[2];
        }
    }

    igraph_real_t cx = (min_x + max_x) / 2.0;
    igraph_real_t cy = (min_y + max_y) / 2.0;
    igraph_real_t cz = (tree->dim == 3) ? ((min_z + max_z) / 2.0) : 0.0;
    igraph_real_t size = max_x - min_x;
    if (max_y - min_y > size) size = max_y - min_y;
    if (tree->dim == 3 && max_z - min_z > size) size = max_z - min_z;
    size *= 1.01;

    IGRAPH_CHECK(ensure_node_capacity(tree, 1));
    igraph_integer_t root_idx = tree->node_count++;

    igraph_real_t root_center[3] = {cx, cy, cz};

    IGRAPH_CHECK(build_bh_node_recursive(
        tree, root_idx, 0, n,
        root_center, size, 0
    ));

    return IGRAPH_SUCCESS;
}

static void calculate_force_node_to_tree(
    const igraph_bh_tree_t *tree,
    igraph_integer_t node_idx,
    const igraph_bh_point_t *point,
    igraph_real_t *force,
    igraph_bh_force_func_t force_func,
    void *user_data
) {
    const igraph_bh_node_t *node = &tree->nodes[node_idx];

    /* Skip empty nodes */
    if (node->point_count == 0) return;

    igraph_real_t dx = point->coord[0] - node->center[0];
    igraph_real_t dy = point->coord[1] - node->center[1];
    igraph_real_t dz = (tree->dim == 3) ? (point->coord[2] - node->center[2]) : 0.0;
    igraph_real_t dist_sq = dx*dx + dy*dy + dz*dz;
    igraph_real_t dist = sqrt(dist_sq);

    if (node->is_leaf) {
        for (igraph_integer_t i = 0; i < node->point_count; i++) {
            igraph_integer_t pid = tree->point_indices[node->data.first_point_idx + i];
            const igraph_bh_point_t *other = &tree->points[pid];
            if (other->id != point->id) {
                igraph_real_t f[3] = {0, 0, 0};
                force_func(point, other, f, user_data);
                force[0] += f[0];
                force[1] += f[1];
                if (tree->dim == 3) {
                    force[2] += f[2];
                }
            }
        }
    } else {
        /* Multipole Acceptance Criterion (MAC) */
        if (dist > 0 && node->size < tree->bh_theta * dist) {
            igraph_bh_point_t pseudo_point;
            pseudo_point.id = -1;
            pseudo_point.mass = node->mass;
            pseudo_point.coord[0] = node->center[0];
            pseudo_point.coord[1] = node->center[1];
            pseudo_point.coord[2] = node->center[2];
            pseudo_point.data = NULL;

            igraph_real_t f[3] = {0, 0, 0};
            force_func(point, &pseudo_point, f, user_data);
            force[0] += f[0];
            force[1] += f[1];
            if (tree->dim == 3) {
                force[2] += f[2];
            }
        } else {
            int n_children = IGRAPH_BH_CHILDREN_COUNT(tree->dim);
            for (int i = 0; i < n_children; i++) {
                calculate_force_node_to_tree(
                    tree, node->data.first_child_idx + i, point, force, force_func, user_data
                );
            }
        }
    }
}

igraph_error_t igraph_bh_calculate_repulsive_forces(
    const igraph_bh_tree_t *tree,
    igraph_matrix_t *forces,
    igraph_bh_force_func_t force_func,
    void *user_data
) {
    if (!tree || tree->point_count == 0 || tree->node_count == 0) {
        return IGRAPH_SUCCESS;
    }

    igraph_integer_t n = tree->point_count;

    for (igraph_integer_t i = 0; i < n; i++) {
        igraph_real_t force[3] = {0, 0, 0};

        /* The root is always at index 0 in the flat array */
        calculate_force_node_to_tree(tree, 0, &tree->points[i], force, force_func, user_data);

        MATRIX(*forces, i, 0) = force[0];
        MATRIX(*forces, i, 1) = force[1];
        if (tree->dim == 3) {
            MATRIX(*forces, i, 2) = force[2];
        }
    }

    return IGRAPH_SUCCESS;
}

igraph_error_t igraph_bh_calculate_attractive_forces(
    const igraph_bh_tree_t *tree,
    const igraph_vector_int_t *from,
    const igraph_vector_int_t *to,
    const igraph_vector_t *weights,
    igraph_matrix_t *forces,
    igraph_bh_force_func_t force_func,
    void *user_data
) {
    if (!from || !to || igraph_vector_int_size(from) == 0) {
        return IGRAPH_SUCCESS;
    }

    igraph_integer_t n_edges = igraph_vector_int_size(from);

    for (igraph_integer_t e = 0; e < n_edges; e++) {
        igraph_integer_t from_idx = VECTOR(*from)[e];
        igraph_integer_t to_idx = VECTOR(*to)[e];
        igraph_real_t w = weights ? VECTOR(*weights)[e] : 1.0;

        igraph_bh_point_t *p1 = &tree->points[from_idx];
        igraph_bh_point_t *p2 = &tree->points[to_idx];

        igraph_real_t f[3] = {0, 0, 0};
        force_func(p1, p2, f, user_data);

        f[0] *= w;
        f[1] *= w;
        if (tree->dim == 3) {
            f[2] *= w;
        }

        MATRIX(*forces, from_idx, 0) += f[0];
        MATRIX(*forces, from_idx, 1) += f[1];
        if (tree->dim == 3) {
            MATRIX(*forces, from_idx, 2) += f[2];
        }

        MATRIX(*forces, to_idx, 0) -= f[0];
        MATRIX(*forces, to_idx, 1) -= f[1];
        if (tree->dim == 3) {
            MATRIX(*forces, to_idx, 2) -= f[2];
        }
    }

    return IGRAPH_SUCCESS;
}

void igraph_bh_tree_print_stats(const igraph_bh_tree_t *tree) {
    printf("Barnes-Hut Tree Stats:\n");
    printf("  Points: %d\n", (int)tree->point_count);
    printf("  Nodes: %d (Allocated Capacity: %d)\n", (int)tree->node_count, (int)tree->capacity);
    printf("  Dimension: %d\n", (int)tree->dim);
    printf("  Theta: %f\n", (double)tree->bh_theta);
    printf("  Max level: %d\n", (int)tree->max_level);
    printf("  Leaf Capacity: %d\n", (int)tree->leaf_capacity);
}
