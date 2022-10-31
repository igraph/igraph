/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2009-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge, MA 02139 USA

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
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#ifndef IGRAPH_CONSTANTS_H
#define IGRAPH_CONSTANTS_H

#include "igraph_config.h"
#include "igraph_decls.h"

__BEGIN_DECLS

/* -------------------------------------------------- */
/* Constants                                          */
/* -------------------------------------------------- */

typedef enum { IGRAPH_UNDIRECTED = 0, IGRAPH_DIRECTED = 1 } igraph_i_directed_t;

/* Note for the enum below: yes, IGRAPH_LOOPS_TWICE is 1, and IGRAPH_LOOPS_ONCE
 * is 2. This is intentional, for the sake of backwards compatibility with
 * earlier versions where we only had IGRAPH_LOOPS and it meant
 * IGRAPH_LOOPS_TWICE */
typedef enum { IGRAPH_NO_LOOPS = 0, IGRAPH_LOOPS = 1, IGRAPH_LOOPS_TWICE = 1, IGRAPH_LOOPS_ONCE = 2 } igraph_loops_t;

typedef enum { IGRAPH_NO_MULTIPLE = 0, IGRAPH_MULTIPLE = 1 } igraph_multiple_t;

typedef enum { IGRAPH_ASCENDING = 0, IGRAPH_DESCENDING = 1 } igraph_order_t;

typedef enum { IGRAPH_MINIMUM = 0, IGRAPH_MAXIMUM = 1 } igraph_optimal_t;

/* Do not renumber the following values! Some internal code treats them as bitmasks
 * and assumes that IGRAPH_ALL == IGRAPH_IN | IGRAPH_OUT and IGRAPH_IN & IGRAPH_OUT == 0. */
typedef enum { IGRAPH_OUT = 1, IGRAPH_IN = 2, IGRAPH_ALL = 3 } igraph_neimode_t;

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

typedef enum { IGRAPH_ERDOS_RENYI_GNP = 0,
               IGRAPH_ERDOS_RENYI_GNM
             } igraph_erdos_renyi_t;

typedef enum { IGRAPH_GET_ADJACENCY_UPPER = 0,
               IGRAPH_GET_ADJACENCY_LOWER,
               IGRAPH_GET_ADJACENCY_BOTH
             } igraph_get_adjacency_t;

typedef enum { IGRAPH_DEGSEQ_CONFIGURATION = 0,     /* Configuration model, allowing non-simple graphs */
               IGRAPH_DEGSEQ_VL,                    /* Viger-Latapy, generates simple connected graphs */
               IGRAPH_DEGSEQ_FAST_HEUR_SIMPLE,      /* Fast heuristic, generates simple graphs */
               IGRAPH_DEGSEQ_CONFIGURATION_SIMPLE,  /* Configuration model, generates simple graphs */
               IGRAPH_DEGSEQ_EDGE_SWITCHING_SIMPLE, /* Edge-switching MCMC, generates simple graphs */

               /* Deprecated, kept for backwards compatibility: */
               IGRAPH_DEGSEQ_SIMPLE IGRAPH_DEPRECATED_ENUMVAL = IGRAPH_DEGSEQ_CONFIGURATION,
               IGRAPH_DEGSEQ_SIMPLE_NO_MULTIPLE IGRAPH_DEPRECATED_ENUMVAL = IGRAPH_DEGSEQ_FAST_HEUR_SIMPLE,
               IGRAPH_DEGSEQ_SIMPLE_NO_MULTIPLE_UNIFORM IGRAPH_DEPRECATED_ENUMVAL = IGRAPH_DEGSEQ_CONFIGURATION_SIMPLE
             } igraph_degseq_t;

typedef enum { IGRAPH_REALIZE_DEGSEQ_SMALLEST = 0,
               IGRAPH_REALIZE_DEGSEQ_LARGEST,
               IGRAPH_REALIZE_DEGSEQ_INDEX
             } igraph_realize_degseq_t;

typedef enum { IGRAPH_RANDOM_TREE_PRUFER = 0,
               IGRAPH_RANDOM_TREE_LERW
             } igraph_random_tree_t;

typedef enum { IGRAPH_FILEFORMAT_EDGELIST = 0,
               IGRAPH_FILEFORMAT_NCOL,
               IGRAPH_FILEFORMAT_PAJEK,
               IGRAPH_FILEFORMAT_LGL,
               IGRAPH_FILEFORMAT_GRAPHML
             } igraph_fileformat_type_t;

typedef enum { IGRAPH_REWIRING_SIMPLE = 0,
               IGRAPH_REWIRING_SIMPLE_LOOPS
             } igraph_rewiring_t;

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

typedef enum { IGRAPH_DONT_SIMPLIFY = 0,
               IGRAPH_SIMPLIFY
             } igraph_lazy_adlist_simplify_t;

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
               IGRAPH_FAS_APPROX_EADES
             } igraph_fas_algorithm_t;

typedef enum { IGRAPH_SUBGRAPH_AUTO = 0,
               IGRAPH_SUBGRAPH_COPY_AND_DELETE,
               IGRAPH_SUBGRAPH_CREATE_FROM_SCRATCH
             } igraph_subgraph_implementation_t;

typedef enum { IGRAPH_IMITATE_AUGMENTED = 0,
               IGRAPH_IMITATE_BLIND,
               IGRAPH_IMITATE_CONTRACTED
             } igraph_imitate_algorithm_t;

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

typedef enum { IGRAPH_ROW_MAJOR = 0,
               IGRAPH_COLUMN_MAJOR = 1
             } igraph_matrix_storage_t;

__END_DECLS

#endif
