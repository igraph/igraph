/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2003, 2004, 2005  Gabor Csardi <csardi@rmki.kfki.hu>
   MTA RMKI, Konkoly-Thege Miklos st. 29-33, Budapest 1121, Hungary
   
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

#ifndef RESTGAME_H
#define RESTGAME_H

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

#ifndef _GNU_SOURCE
# define _GNU_SOURCE
#endif

#include "igraph_types.h"
#include "igraph_error.h"
#include "igraph_interrupt.h"
#include "igraph_arpack.h"

#include <stdio.h> 		/* FILE */

/* -------------------------------------------------- */
/* Constants                                          */
/* -------------------------------------------------- */

typedef enum { IGRAPH_UNDIRECTED=0, IGRAPH_DIRECTED=1 } igraph_i_directed_t;

typedef enum { IGRAPH_NO_LOOPS=0, IGRAPH_LOOPS=1 } igraph_i_loops_t;

typedef enum { IGRAPH_OUT=1, IGRAPH_IN=2, IGRAPH_ALL=3,
	       IGRAPH_TOTAL=3 } igraph_neimode_t;

typedef enum { IGRAPH_WEAK=1, IGRAPH_STRONG=2 } igraph_connectedness_t;

typedef enum { IGRAPH_ADJ_DIRECTED=0, 
	       IGRAPH_ADJ_UNDIRECTED=1, IGRAPH_ADJ_MAX=1,
               IGRAPH_ADJ_UPPER, IGRAPH_ADJ_LOWER, IGRAPH_ADJ_MIN,
	       IGRAPH_ADJ_PLUS } igraph_adjacency_t;

typedef enum { IGRAPH_STAR_OUT=0, IGRAPH_STAR_IN,
	       IGRAPH_STAR_UNDIRECTED } igraph_star_mode_t;

typedef enum { IGRAPH_TREE_OUT=0, IGRAPH_TREE_IN,
	       IGRAPH_TREE_UNDIRECTED } igraph_tree_mode_t;

typedef enum { IGRAPH_ERDOS_RENYI_GNP=0, 
	       IGRAPH_ERDOS_RENYI_GNM } igraph_erdos_renyi_t;

typedef enum { IGRAPH_GET_ADJACENCY_UPPER=0,
	       IGRAPH_GET_ADJACENCY_LOWER,
	       IGRAPH_GET_ADJACENCY_BOTH } igraph_get_adjacency_t;

typedef enum { IGRAPH_DEGSEQ_SIMPLE=0,
	       IGRAPH_DEGSEQ_VL } igraph_degseq_t;

typedef enum { IGRAPH_FILEFORMAT_EDGELIST=0,
	       IGRAPH_FILEFORMAT_NCOL,
	       IGRAPH_FILEFORMAT_PAJEK,
               IGRAPH_FILEFORMAT_LGL,
               IGRAPH_FILEFORMAT_GRAPHML } igraph_fileformat_type_t;

typedef enum { IGRAPH_REWIRING_SIMPLE=0 } igraph_rewiring_t;

typedef enum { IGRAPH_EDGEORDER_ID=0,
	       IGRAPH_EDGEORDER_FROM,
	       IGRAPH_EDGEORDER_TO } igraph_edgeorder_type_t;

typedef enum { IGRAPH_TO_DIRECTED_ARBITRARY=0,
	       IGRAPH_TO_DIRECTED_MUTUAL } igraph_to_directed_t;

typedef enum { IGRAPH_TO_UNDIRECTED_EACH=0,
	       IGRAPH_TO_UNDIRECTED_COLLAPSE } igraph_to_undirected_t;

typedef enum { IGRAPH_VCONN_NEI_ERROR=0,
	       IGRAPH_VCONN_NEI_INFINITY,
	       IGRAPH_VCONN_NEI_IGNORE }igraph_vconn_nei_t;

typedef enum { IGRAPH_SPINCOMM_UPDATE_SIMPLE=0,
	       IGRAPH_SPINCOMM_UPDATE_CONFIG } igraph_spincomm_update_t; 

typedef enum { IGRAPH_DONT_SIMPLIFY=0,
	       IGRAPH_SIMPLIFY } igraph_lazy_adlist_simplify_t;

typedef enum { IGRAPH_TRANSITIVITY_NAN=0,
               IGRAPH_TRANSITIVITY_ZERO } igraph_transitivity_mode_t;

typedef enum { IGRAPH_SPINCOMM_IMP_ORIG=0, 
	       IGRAPH_SPINCOMM_IMP_NEG } igraph_spinglass_implementation_t;

typedef igraph_real_t  igraph_scalar_function_t(const igraph_vector_t *var, 
						const igraph_vector_t *par,
						void* extra);
typedef void igraph_vector_function_t(const igraph_vector_t *var, 
				      const igraph_vector_t *par,
				      igraph_vector_t* res, void* extra);


#include "igraph_datatype.h"
#include "igraph_iterators.h"
#include "igraph_interface.h"
#include "igraph_constructors.h"
#include "igraph_games.h"
#include "igraph_centrality.h"
#include "igraph_paths.h"
#include "igraph_components.h"
#include "igraph_structural.h"
#include "igraph_transitivity.h"
#include "igraph_neighborhood.h"
#include "igraph_topology.h"
#include "igraph_bipartite.h"
#include "igraph_cliques.h"
#include "igraph_layout.h"
#include "igraph_visitor.h"
#include "igraph_community.h"
#include "igraph_conversion.h"
#include "igraph_foreign.h"
#include "igraph_motifs.h"
#include "igraph_operators.h"
#include "igraph_flow.h"
#include "igraph_revolver.h"
#include "igraph_nongraph.h"
#include "igraph_cocitation.h"
#include "igraph_adjlist.h"
#include "igraph_progress.h"

/* -------------------------------------------------- */
/* Eigenvectors and eigenvalues                       */
/* -------------------------------------------------- */

int igraph_eigen_tred2(const igraph_matrix_t *A,
		       igraph_vector_t *D,
		       igraph_vector_t *E,
		       igraph_matrix_t *Z);

int igraph_eigen_tql2(igraph_vector_t *D,
		      igraph_vector_t *E,
		      igraph_matrix_t *Z);

int igraph_eigen_tred1(const igraph_matrix_t *A,
		       igraph_vector_t *D,
		       igraph_vector_t *E2);

int igraph_eigen_tqlrat(igraph_vector_t *D,
			igraph_vector_t *E2);

int igraph_eigen_rs(const igraph_matrix_t *A,
		    igraph_vector_t *values,
		    igraph_matrix_t *vectors);


extern unsigned int igraph_i_isoclass_3[];
extern unsigned int igraph_i_isoclass_4[];
extern unsigned int igraph_i_isoclass_3u[];
extern unsigned int igraph_i_isoclass_4u[];
extern unsigned int igraph_i_isoclass2_3[];
extern unsigned int igraph_i_isoclass2_4[];
extern unsigned int igraph_i_isoclass2_3u[];
extern unsigned int igraph_i_isoclass2_4u[];
extern unsigned int igraph_i_isoclass_3_idx[];
extern unsigned int igraph_i_isoclass_4_idx[];
extern unsigned int igraph_i_isoclass_3u_idx[];
extern unsigned int igraph_i_isoclass_4u_idx[];

#include "igraph_attributes.h"

__END_DECLS
  
#endif
