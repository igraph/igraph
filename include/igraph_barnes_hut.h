/*
   igraph library.
   Copyright (C) 2026  The igraph development team <igraph@igraph.org>

   This program is free software; you can redistribute it and/or modify ...
   [License block truncated for brevity]
 */

#ifndef IGRAPH_BARNES_HUT_H
#define IGRAPH_BARNES_HUT_H

#include "igraph_decls.h"
#include "igraph_types.h"
#include "igraph_matrix.h"
#include "igraph_vector.h"

IGRAPH_BEGIN_C_DECLS

/* Represents a physical entity in the space */
typedef struct igraph_bh_point_t {
    igraph_real_t coord[3];
    igraph_real_t mass;
    igraph_integer_t id;     /* Typically the graph vertex ID */
    void *data;              /* Optional user payload */
} igraph_bh_point_t;

/* * Flat node structure.
 * Pointers have been removed in favor of integer offsets into global arrays
 * to guarantee optimal CPU cache locality and zero per-node allocations.
 */
typedef struct igraph_bh_node_t {
    igraph_real_t mass;
    igraph_real_t center[3];
    igraph_real_t size;

    igraph_integer_t point_count;
    igraph_bool_t is_leaf;

    union {
        /* If is_leaf == true: offset into tree->point_indices */
        igraph_integer_t first_point_idx;

        /* If is_leaf == false: offset into tree->nodes for the first of 2^dim children */
        igraph_integer_t first_child_idx;
    } data;
} igraph_bh_node_t;

/* The main tree context */
typedef struct {
    igraph_bh_node_t *nodes;             /* Contiguous array of all nodes */
    igraph_integer_t node_count;
    igraph_integer_t capacity;           /* Allocated size of 'nodes' array */

    igraph_bh_point_t *points;           /* Contiguous array of point data */
    igraph_integer_t *point_indices;     /* Permuted indices mapped to leaf nodes */
    igraph_integer_t point_count;

    igraph_integer_t dim;                /* 2 or 3 */
    igraph_integer_t max_level;          /* Maximum recursion depth */
    igraph_integer_t leaf_capacity;      /* Max points allowed in a leaf before subdividing */
    igraph_real_t bh_theta;              /* Barnes-Hut MAC threshold */
} igraph_bh_tree_t;

/* =========================================================================
 * PHYSICS PLUGINS (CALLBACK KERNELS)
 * ========================================================================= */

/**
 * @brief Repulsion Kernel (O(N log N) Spatial Domain)
 * Computes the force exerted ON p1 BY p2 (which may be a macroscopic pseudo-node).
 */
typedef void (*igraph_bh_repulsion_kernel_t)(
    const igraph_bh_point_t *p1,
    const igraph_bh_point_t *p2,
    igraph_real_t dx,       /* Pre-computed: p1->coord - p2->coord */
    igraph_real_t dy,
    igraph_real_t dz,
    igraph_real_t dist_sq,  /* Pre-computed, guaranteed > 0 by BH Engine */
    igraph_real_t force[3], /* OUTPUT: The force vector to apply to p1 */
    void *user_data         /* Context for layout hyperparameters */
);

/**
 * @brief Attraction Kernel (O(E) Topological Domain)
 * Computes the forces exerted across a specific directed/undirected edge.
 */
typedef void (*igraph_bh_attraction_kernel_t)(
    const igraph_bh_point_t *p1,    /* Source */
    const igraph_bh_point_t *p2,    /* Target */
    igraph_real_t dx,               /* Pre-computed: p1->coord - p2->coord */
    igraph_real_t dy,
    igraph_real_t dz,
    igraph_real_t dist_sq,          /* Pre-computed, guaranteed > 0 by BH Engine */
    igraph_real_t weight,           /* The edge weight */
    igraph_real_t force_p1[3],      /* OUTPUT: Force applied to p1 */
    igraph_real_t force_p2[3],      /* OUTPUT: Force applied to p2 */
    void *user_data                 /* Context for layout hyperparameters */
);

/* =========================================================================
 * LIFECYCLE & ITERATORS
 * ========================================================================= */

igraph_error_t igraph_bh_tree_init(
    igraph_bh_tree_t *tree,
    igraph_integer_t dim,
    igraph_real_t theta,
    igraph_integer_t max_level,
    igraph_integer_t leaf_capacity
);

void igraph_bh_tree_get_scaling_params(
    igraph_integer_t n_points,
    igraph_integer_t dim,
    igraph_integer_t *max_level,
    igraph_integer_t *leaf_capacity
);

void igraph_bh_tree_destroy(igraph_bh_tree_t *tree);

/**
 * Build the tree layout using the provided graph coordinates and masses.
 */
igraph_error_t igraph_bh_tree_build(
    igraph_bh_tree_t *tree,
    const igraph_matrix_t *coords,
    const igraph_vector_t *masses
);

/* Applies the repulsion kernel using the O(N log N) spatial tree */
igraph_error_t igraph_bh_apply_repulsion_from_tree(
    const igraph_bh_tree_t *tree,
    igraph_matrix_t *forces,
    igraph_bh_repulsion_kernel_t kernel,
    void *user_data
);

/* Applies the attraction kernel using the O(E) topological edge list */
igraph_error_t igraph_bh_apply_attraction_from_edges(
    const igraph_bh_tree_t *tree, /* Passed so callbacks can access point attributes (mass, id) */
    const igraph_vector_int_t *from,
    const igraph_vector_int_t *to,
    const igraph_vector_t *weights,
    igraph_matrix_t *forces,
    igraph_bh_attraction_kernel_t kernel,
    void *user_data
);

void igraph_bh_tree_print_stats(const igraph_bh_tree_t *tree);

IGRAPH_END_C_DECLS

#endif
