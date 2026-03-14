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

/* Callback for calculating pairwise forces */
typedef void (*igraph_bh_force_func_t)(
    const igraph_bh_point_t *p1,
    const igraph_bh_point_t *p2,
    igraph_real_t *force,
    void *user_data
);

/**
 * Initialize the tree parameters. Does NOT allocate point/node memory yet.
 */
igraph_error_t igraph_bh_tree_init(
    igraph_bh_tree_t *tree,
    igraph_integer_t dim,
    igraph_real_t theta,
    igraph_integer_t max_level,
    igraph_integer_t leaf_capacity
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

igraph_error_t igraph_bh_calculate_repulsive_forces(
    const igraph_bh_tree_t *tree,
    igraph_matrix_t *forces,
    igraph_bh_force_func_t force_func,
    void *user_data
);

igraph_error_t igraph_bh_calculate_attractive_forces(
    const igraph_bh_tree_t *tree,
    const igraph_vector_int_t *from,
    const igraph_vector_int_t *to,
    const igraph_vector_t *weights,
    igraph_matrix_t *forces,
    igraph_bh_force_func_t force_func,
    void *user_data
);

void igraph_bh_tree_print_stats(const igraph_bh_tree_t *tree);

IGRAPH_END_C_DECLS

#endif
