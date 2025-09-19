/*
   igraph library.
   Copyright (C) 2009-2025  The igraph development team <igraph@igraph.org>

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

#ifndef IGRAPH_CONSTANTS_H
#define IGRAPH_CONSTANTS_H

#include "igraph_decls.h"

IGRAPH_BEGIN_C_DECLS

/* -------------------------------------------------- */
/* Constants                                          */
/* -------------------------------------------------- */

/**
 * \define IGRAPH_UNLIMITED
 * \brief Constant for "do not limit results".
 *
 * A constant signifying that no limitation should be used with various cutoff,
 * size limit or result set size parameters, such as minimum or maximum clique
 * size, number of results returned, cutoff for path lengths, etc. Currently
 * defined to <code>-1</code>.
 */
#define IGRAPH_UNLIMITED (-1)
/* Note to maintainers: IGRAPH_UNLIMITED is intended to support readability
 * when *negative* parameter values indicate "no limit". Do not test
 * directly against IGRAPH_UNLIMITED in implementations, test for negative
 * values instead. */

/* These constants are meant to be used for sake of readability */
enum { IGRAPH_UNDIRECTED = 0, IGRAPH_DIRECTED = 1 };
enum { IGRAPH_NO_MULTIPLE = 0, IGRAPH_MULTIPLE = 1 };
enum { IGRAPH_EDGE_UNLABELED = 0, IGRAPH_EDGE_LABELED = 1 };

/**
 * \typedef igraph_loops_t
 * \brief How to interpret self-loops in undirected graphs?
 *
 * Controls the interpretation of self-loops in undirected graphs, typically
 * in the context of adjacency matrices or degrees.
 *
 * </para><para>These constants are also used to improve readability in
 * boolean contexts, with \c IGRAPH_NO_LOOPS, equivalent to \c false,
 * signifying that loops should be ignored and \c IGRAPH_LOOPS, equivalent
 * to \c true, that loops should be considered.
 *
 * \enumval IGRAPH_NO_LOOPS Self-loops are ignored.
 * \enumval IGRAPH_LOOPS_TWICE Self-loops are considered, and counted twice
 *    in undirected graphs. For example, a self-loop contributes two to the
 *    degree of a vertex and to diagonal entries of adjacency matrices. This
 *    is the standard interpretation in graph theory, thus \c IGRAPH_LOOPS
 *    serves as an alias for this option.
 * \enumval IGRAPH_LOOPS_ONCE Self-loops are considered, and counted only
 *    once in undirected graphs.
 */
typedef enum {
    IGRAPH_NO_LOOPS = 0,
    IGRAPH_LOOPS_TWICE = 1,
    IGRAPH_LOOPS_ONCE = 2,
    IGRAPH_LOOPS = IGRAPH_LOOPS_TWICE
} igraph_loops_t;
/* Note for the enum above: yes, IGRAPH_LOOPS_TWICE is 1, and IGRAPH_LOOPS_ONCE
 * is 2. This is intentional, for the sake of backwards compatibility with
 * earlier versions where we only had IGRAPH_LOOPS and it meant
 * IGRAPH_LOOPS_TWICE */

typedef enum { IGRAPH_ASCENDING = 0, IGRAPH_DESCENDING = 1 } igraph_order_t;

/**
 * \typedef igraph_neimode_t
 * \brief How to interpret edge directions in directed graphs?
 *
 * These "neighbor mode" constants are typically used to specify the treatment
 * of edge directions in directed graphs, or which vertices to consider as
 * adjacent to (i.e. neighbor of) a vertex. It is typically ignored in undirected
 * graphs.
 *
 * \enumval IGRAPH_OUT Follow edge directions in directed graphs, or consider
 *    out-neighbors of vertices.
 * \enumval IGRAPH_IN Follow edges in the reverse direction in directed graphs,
 *    or consider in-neighbors of vertices.
 * \enumval IGRAPH_ALL Ignore edge directions in directed graphs, or consider
 *    all neighbours (both out and in-neighbors) of vertices.
 */
typedef enum {
    IGRAPH_OUT = 1,
    IGRAPH_IN = 2,
    IGRAPH_ALL = 3
} igraph_neimode_t;
/* Do not renumber the vaues above! Some internal code treats them as bitmasks
 * and assumes that IGRAPH_ALL == IGRAPH_IN | IGRAPH_OUT and IGRAPH_IN & IGRAPH_OUT == 0. */

/* Reverse IGRAPH_OUT to IGRAPH_IN and vice versa. Leave other values alone. */
#define IGRAPH_REVERSE_MODE(mode) \
    ((mode) == IGRAPH_IN ? IGRAPH_OUT : ((mode) == IGRAPH_OUT ? IGRAPH_IN : (mode)))

typedef enum { IGRAPH_WEAK = 1, IGRAPH_STRONG = 2 } igraph_connectedness_t;

typedef enum { IGRAPH_RECIPROCITY_DEFAULT = 0,
               IGRAPH_RECIPROCITY_RATIO = 1
             } igraph_reciprocity_t;

typedef enum { IGRAPH_ADJ_DIRECTED = 0,
               IGRAPH_ADJ_UNDIRECTED,
               IGRAPH_ADJ_UPPER, IGRAPH_ADJ_LOWER, IGRAPH_ADJ_MIN,
               IGRAPH_ADJ_PLUS,
               IGRAPH_ADJ_MAX,
             } igraph_adjacency_t;

typedef enum { IGRAPH_STAR_OUT = 0, IGRAPH_STAR_IN,
               IGRAPH_STAR_UNDIRECTED,
               IGRAPH_STAR_MUTUAL
             } igraph_star_mode_t;

typedef enum { IGRAPH_WHEEL_OUT = 0, IGRAPH_WHEEL_IN,
               IGRAPH_WHEEL_UNDIRECTED,
               IGRAPH_WHEEL_MUTUAL
             } igraph_wheel_mode_t;

typedef enum { IGRAPH_TREE_OUT = 0, IGRAPH_TREE_IN,
               IGRAPH_TREE_UNDIRECTED
             } igraph_tree_mode_t;

typedef enum { IGRAPH_GET_ADJACENCY_UPPER = 0,
               IGRAPH_GET_ADJACENCY_LOWER,
               IGRAPH_GET_ADJACENCY_BOTH
             } igraph_get_adjacency_t;

typedef enum { IGRAPH_DEGSEQ_CONFIGURATION = 0,     /* Configuration model, allowing non-simple graphs */
               IGRAPH_DEGSEQ_VL,                    /* Viger-Latapy, generates simple connected graphs */
               IGRAPH_DEGSEQ_FAST_HEUR_SIMPLE,      /* Fast heuristic, generates simple graphs */
               IGRAPH_DEGSEQ_CONFIGURATION_SIMPLE,  /* Configuration model, generates simple graphs */
               IGRAPH_DEGSEQ_EDGE_SWITCHING_SIMPLE, /* Edge-switching MCMC, generates simple graphs */
             } igraph_degseq_t;

typedef enum { IGRAPH_REALIZE_DEGSEQ_SMALLEST = 0,
               IGRAPH_REALIZE_DEGSEQ_LARGEST,
               IGRAPH_REALIZE_DEGSEQ_INDEX
             } igraph_realize_degseq_t;

typedef enum { IGRAPH_RANDOM_TREE_PRUFER = 0,
               IGRAPH_RANDOM_TREE_LERW
             } igraph_random_tree_t;

typedef enum { IGRAPH_EDGEORDER_ID = 0,
               IGRAPH_EDGEORDER_FROM,
               IGRAPH_EDGEORDER_TO
             } igraph_edgeorder_type_t;

typedef enum { IGRAPH_TO_DIRECTED_ARBITRARY = 0,
               IGRAPH_TO_DIRECTED_MUTUAL,
               IGRAPH_TO_DIRECTED_RANDOM,
               IGRAPH_TO_DIRECTED_ACYCLIC
             } igraph_to_directed_t;

typedef enum { IGRAPH_TO_UNDIRECTED_EACH = 0,
               IGRAPH_TO_UNDIRECTED_COLLAPSE,
               IGRAPH_TO_UNDIRECTED_MUTUAL
             } igraph_to_undirected_t;

typedef enum { IGRAPH_VCONN_NEI_ERROR = 0,
               IGRAPH_VCONN_NEI_NUMBER_OF_NODES,
               IGRAPH_VCONN_NEI_IGNORE,
               IGRAPH_VCONN_NEI_NEGATIVE
             } igraph_vconn_nei_t;

typedef enum { IGRAPH_SPINCOMM_UPDATE_SIMPLE = 0,
               IGRAPH_SPINCOMM_UPDATE_CONFIG
             } igraph_spincomm_update_t;

typedef enum { IGRAPH_TRANSITIVITY_NAN = 0,
               IGRAPH_TRANSITIVITY_ZERO
             } igraph_transitivity_mode_t;

typedef enum { IGRAPH_SPINCOMM_IMP_ORIG = 0,
               IGRAPH_SPINCOMM_IMP_NEG
             } igraph_spinglass_implementation_t;

typedef enum { IGRAPH_COMMCMP_VI = 0,
               IGRAPH_COMMCMP_NMI,
               IGRAPH_COMMCMP_SPLIT_JOIN,
               IGRAPH_COMMCMP_RAND,
               IGRAPH_COMMCMP_ADJUSTED_RAND
             } igraph_community_comparison_t;

typedef enum { IGRAPH_ADD_WEIGHTS_NO = 0,
               IGRAPH_ADD_WEIGHTS_YES,
               IGRAPH_ADD_WEIGHTS_IF_PRESENT
             } igraph_add_weights_t;

typedef enum { IGRAPH_BARABASI_BAG = 0,
               IGRAPH_BARABASI_PSUMTREE,
               IGRAPH_BARABASI_PSUMTREE_MULTIPLE
             } igraph_barabasi_algorithm_t;

typedef enum { IGRAPH_FAS_EXACT_IP = 0,
               IGRAPH_FAS_APPROX_EADES,
               IGRAPH_FAS_EXACT_IP_CG,
               IGRAPH_FAS_EXACT_IP_TI
             } igraph_fas_algorithm_t;

typedef enum { IGRAPH_FVS_EXACT_IP = 0
             } igraph_fvs_algorithm_t;

typedef enum { IGRAPH_SUBGRAPH_AUTO = 0,
               IGRAPH_SUBGRAPH_COPY_AND_DELETE,
               IGRAPH_SUBGRAPH_CREATE_FROM_SCRATCH
             } igraph_subgraph_implementation_t;

typedef enum { IGRAPH_LAYOUT_GRID = 0,
               IGRAPH_LAYOUT_NOGRID,
               IGRAPH_LAYOUT_AUTOGRID
             } igraph_layout_grid_t;

typedef enum { IGRAPH_RANDOM_WALK_STUCK_ERROR = 0,
               IGRAPH_RANDOM_WALK_STUCK_RETURN
             } igraph_random_walk_stuck_t;

typedef enum { IGRAPH_VORONOI_FIRST = 0,
               IGRAPH_VORONOI_LAST,
               IGRAPH_VORONOI_RANDOM
             } igraph_voronoi_tiebreaker_t;

typedef enum { IGRAPH_CHUNG_LU_ORIGINAL = 0,
               IGRAPH_CHUNG_LU_MAXENT,
               IGRAPH_CHUNG_LU_NR
             } igraph_chung_lu_t;

typedef enum { IGRAPH_ROW_MAJOR = 0,
               IGRAPH_COLUMN_MAJOR = 1
             } igraph_matrix_storage_t;

typedef enum { IGRAPH_MST_AUTOMATIC = 0,
               IGRAPH_MST_UNWEIGHTED,
               IGRAPH_MST_PRIM,
               IGRAPH_MST_KRUSKAL
             } igraph_mst_algorithm_t;

typedef enum { IGRAPH_PRODUCT_CARTESIAN = 0,
               IGRAPH_PRODUCT_LEXICOGRAPHIC,
               IGRAPH_PRODUCT_STRONG,
               IGRAPH_PRODUCT_TENSOR,
               IGRAPH_PRODUCT_MODULAR
             } igraph_product_t;

/**
 * \typedef igraph_lpa_variant_t
 * \brief Label propagation algorithm variants of implementation
 *
 * Variants to run the label propagation algorithm.
 * \enumval IGRAPH_LPA_DOMINANCE Check for dominance of all nodes after each iteration
 * \enumval IGRAPH_LPA_RETENTION Keep current label if among dominant labels, only check if labels changed
 * \enumval IGRAPH_LPA_FAST Sample from dominant labels, only check neighbors
 */
typedef enum {
    IGRAPH_LPA_DOMINANCE = 0, /* Sample from dominant labels, check for dominance after each iteration. */
    IGRAPH_LPA_RETENTION,     /* Keep current label if among dominant labels, only check if labels changed. */
    IGRAPH_LPA_FAST           /* Sample from dominant labels, only check neighbors. */
} igraph_lpa_variant_t;

IGRAPH_END_C_DECLS

#endif
