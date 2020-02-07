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

#include "igraph_decls.h"
#include "igraph_types.h"
#include "igraph_datatype.h"

__BEGIN_DECLS

/* -------------------------------------------------- */
/* Constants                                          */
/* -------------------------------------------------- */

typedef enum { IGRAPH_UNDIRECTED = 0, IGRAPH_DIRECTED = 1 } igraph_i_directed_t;

typedef enum { IGRAPH_NO_LOOPS = 0, IGRAPH_LOOPS = 1 } igraph_i_loops_t;

typedef enum { IGRAPH_NO_MULTIPLE = 0, IGRAPH_MULTIPLE = 1 } igraph_i_multiple_t;

typedef enum { IGRAPH_ASCENDING = 0, IGRAPH_DESCENDING = 1 } igraph_order_t;

typedef enum { IGRAPH_MINIMUM = 0, IGRAPH_MAXIMUM = 1 } igraph_optimal_t;

typedef enum { IGRAPH_OUT = 1, IGRAPH_IN = 2, IGRAPH_ALL = 3,
               IGRAPH_TOTAL = 3
             } igraph_neimode_t;

typedef enum { IGRAPH_WEAK = 1, IGRAPH_STRONG = 2 } igraph_connectedness_t;

typedef enum { IGRAPH_RECIPROCITY_DEFAULT = 0,
               IGRAPH_RECIPROCITY_RATIO = 1
             } igraph_reciprocity_t;

typedef enum { IGRAPH_ADJ_DIRECTED = 0,
               IGRAPH_ADJ_UNDIRECTED = 1, IGRAPH_ADJ_MAX = 1,
               IGRAPH_ADJ_UPPER, IGRAPH_ADJ_LOWER, IGRAPH_ADJ_MIN,
               IGRAPH_ADJ_PLUS
             } igraph_adjacency_t;

typedef enum { IGRAPH_STAR_OUT = 0, IGRAPH_STAR_IN,
               IGRAPH_STAR_UNDIRECTED,
               IGRAPH_STAR_MUTUAL
             } igraph_star_mode_t;

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

typedef enum { IGRAPH_DEGSEQ_SIMPLE = 0,
               IGRAPH_DEGSEQ_VL,
               IGRAPH_DEGSEQ_SIMPLE_NO_MULTIPLE,
               IGRAPH_DEGSEQ_SIMPLE_NO_MULTIPLE_UNIFORM
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
               IGRAPH_TO_DIRECTED_MUTUAL
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

typedef igraph_real_t  igraph_scalar_function_t(const igraph_vector_t *var,
        const igraph_vector_t *par,
        void* extra);
typedef void igraph_vector_function_t(const igraph_vector_t *var,
                                      const igraph_vector_t *par,
                                      igraph_vector_t* res, void* extra);

typedef enum { IGRAPH_LAYOUT_GRID = 0,
               IGRAPH_LAYOUT_NOGRID,
               IGRAPH_LAYOUT_AUTOGRID
             } igraph_layout_grid_t;

typedef enum { IGRAPH_RANDOM_WALK_STUCK_ERROR = 0,
               IGRAPH_RANDOM_WALK_STUCK_RETURN
             } igraph_random_walk_stuck_t;


__END_DECLS

#endif
