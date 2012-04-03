/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2009  Gabor Csardi <csardi@rmki.kfki.hu>
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

#ifndef IGRAPH_TOPOLOGY_H
#define IGRAPH_TOPOLOGY_H

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

#include "igraph_constants.h"
#include "igraph_datatype.h"
#include "igraph_types.h"
#include "igraph_vector_ptr.h"

__BEGIN_DECLS

int igraph_topological_sorting(const igraph_t *graph, igraph_vector_t *res,
			       igraph_neimode_t mode);
int igraph_is_dag(const igraph_t *graph, igraph_bool_t *res);

/* -------------------------------------------------- */
/* Graph isomorphisms                                 */
/* -------------------------------------------------- */

/* Common functions */
int igraph_permute_vertices(const igraph_t *graph, igraph_t *res,
			    const igraph_vector_t *permutation);

/* Generic interface */
int igraph_isomorphic(const igraph_t *graph1, const igraph_t *graph2,
		      igraph_bool_t *iso);
int igraph_subisomorphic(const igraph_t *graph1, const igraph_t *graph2,
			 igraph_bool_t *iso);

/* VF2 family*/
/** 
 * \typedef igraph_isohandler_t
 * Callback type, called when an isomorphism was found
 * 
 * See the details at the documentation of \ref
 * igraph_isomorphic_function_vf2().
 * \param map12 The mapping from the first graph to the second.
 * \param map21 The mapping from the second graph to the first, the
 *   inverse of \p map12 basically.
 * \param arg This extra argument was passed to \ref
 *   igraph_isomorphic_function_vf2() when it was called.
 * \return Boolean, whether to continue with the isomorphism search.
 */


typedef igraph_bool_t igraph_isohandler_t(const igraph_vector_t *map12, 
					  const igraph_vector_t *map21, void *arg);

typedef igraph_bool_t igraph_isocompat_t(const igraph_t *graph1,
					 const igraph_t *graph2,
					 const igraph_integer_t g1_num,
					 const igraph_integer_t g2_num,
					 void *arg);

/**
 * \struct igraph_isomorphic_function_vf2_args_t 
 * Information about a BLISS run
 * Parameter bundle to be used for different subfunctions of vf2 algorithms.
 * 
 * \member isohandler_fn_arg Argument passed through to isohandler_fn
 * \member node_compat_fn_arg Argument passed through to node_compat_fn
 * \member edge_compat_fn_arg Argument passed through to edge_compat_fn
 */

typedef struct{
	void *isohandler_fn_arg;
	void *node_compat_fn_arg;
	void *edge_compat_fn_arg;
} igraph_isomorphic_function_vf2_args_t;

int igraph_isomorphic_vf2(const igraph_t *graph1, const igraph_t *graph2, 
			  const igraph_vector_int_t *vertex_color1,
			  const igraph_vector_int_t *vertex_color2,
			  const igraph_vector_int_t *edge_color1,
			  const igraph_vector_int_t *edge_color2,
			  igraph_bool_t *iso,
			  igraph_vector_t *map12, 
			  igraph_vector_t *map21,
			  igraph_isocompat_t *node_compat_fn,
			  igraph_isocompat_t *edge_compat_fn,
			  igraph_isomorphic_function_vf2_args_t *args);
int igraph_isomorphic_function_vf2(const igraph_t *graph1, const igraph_t *graph2,
				   const igraph_vector_int_t *vertex_color1,
				   const igraph_vector_int_t *vertex_color2,
				   const igraph_vector_int_t *edge_color1,
				   const igraph_vector_int_t *edge_color2,
				   igraph_vector_t *map12, igraph_vector_t *map21,
				   igraph_isohandler_t *isohandler_fn,
				   igraph_isocompat_t *node_compat_fn,
				   igraph_isocompat_t *edge_compat_fn,
				   igraph_isomorphic_function_vf2_args_t *args);
int igraph_count_isomorphisms_vf2(const igraph_t *graph1, const igraph_t *graph2, 
				  const igraph_vector_int_t *vertex_color1,
				  const igraph_vector_int_t *vertex_color2,
				  const igraph_vector_int_t *edge_color1,
				  const igraph_vector_int_t *edge_color2,
				  igraph_integer_t *count,
				  igraph_isocompat_t *node_compat_fn,
				  igraph_isocompat_t *edge_compat_fn,
				  igraph_isomorphic_function_vf2_args_t *args);
int igraph_get_isomorphisms_vf2(const igraph_t *graph1,
				const igraph_t *graph2,
				const igraph_vector_int_t *vertex_color1,
				const igraph_vector_int_t *vertex_color2,
				const igraph_vector_int_t *edge_color1,
				const igraph_vector_int_t *edge_color2,
				igraph_vector_ptr_t *maps,
				igraph_isocompat_t *node_compat_fn,
				igraph_isocompat_t *edge_compat_fn,
				igraph_isomorphic_function_vf2_args_t *args);

int igraph_subisomorphic_vf2(const igraph_t *graph1, const igraph_t *graph2, 
			     const igraph_vector_int_t *vertex_color1,
			     const igraph_vector_int_t *vertex_color2,
			     const igraph_vector_int_t *edge_color1,
			     const igraph_vector_int_t *edge_color2,
			     igraph_bool_t *iso,
			     igraph_vector_t *map12, 
			     igraph_vector_t *map21,
			     igraph_isocompat_t *node_compat_fn,
			     igraph_isocompat_t *edge_compat_fn,
				 igraph_isomorphic_function_vf2_args_t *args);
int igraph_subisomorphic_function_vf2(const igraph_t *graph1, 
				      const igraph_t *graph2,
				      const igraph_vector_int_t *vertex_color1,
				      const igraph_vector_int_t *vertex_color2,
				      const igraph_vector_int_t *edge_color1,
				      const igraph_vector_int_t *edge_color2,
				      igraph_vector_t *map12,
				      igraph_vector_t *map21,
				      igraph_isohandler_t *isohandler_fn,
				      igraph_isocompat_t *node_compat_fn,
				      igraph_isocompat_t *edge_compat_fn,
					  igraph_isomorphic_function_vf2_args_t *args);
int igraph_count_subisomorphisms_vf2(const igraph_t *graph1, const igraph_t *graph2, 
				     const igraph_vector_int_t *vertex_color1,
				     const igraph_vector_int_t *vertex_color2,
				     const igraph_vector_int_t *edge_color1,
				     const igraph_vector_int_t *edge_color2,
				     igraph_integer_t *count,
				     igraph_isocompat_t *node_compat_fn,
				     igraph_isocompat_t *edge_compat_fn,
					 igraph_isomorphic_function_vf2_args_t *args);
int igraph_get_subisomorphisms_vf2(const igraph_t *graph1,
				   const igraph_t *graph2,
				   const igraph_vector_int_t *vertex_color1,
				   const igraph_vector_int_t *vertex_color2,
				   const igraph_vector_int_t *edge_color1,
				   const igraph_vector_int_t *edge_color2,
				   igraph_vector_ptr_t *maps,
				   igraph_isocompat_t *node_compat_fn,
				   igraph_isocompat_t *edge_compat_fn,
				   igraph_isomorphic_function_vf2_args_t *args);

/* BLISS family */
/**
 * \struct igraph_bliss_info_t 
 * Information about a BLISS run
 * 
 * Some secondary information found by the BLISS algorithm is stored
 * here. It is useful if you wany to study the internal working of the
 * algorithm.
 * \member nof_nodes The number of nodes in the search tree.
 * \member nof_leaf_nodes The number of leaf nodes in the search tree.
 * \member nof_bad_nodes Number of bad nodes.
 * \member nof_canupdates Number of canrep updates.
 * \member max_level Maximum level.
 * \member group_size The size of the automorphism group of the graph,
 *    given as a string. It should be deallocated via
 *    <function>free()</function> if not needed any more.
 * 
 * See http://www.tcs.hut.fi/Software/bliss/index.html
 * for details about the algorithm and these parameters.
 */
typedef struct igraph_bliss_info_t {
  unsigned long nof_nodes;
  unsigned long nof_leaf_nodes;
  unsigned long nof_bad_nodes;
  unsigned long nof_canupdates;
  unsigned long max_level;
  char *group_size;
} igraph_bliss_info_t;

/**
 * \typedef igraph_bliss_sh_t
 * Splitting heuristics for BLISS
 * 
 * \enumval IGRAPH_BLISS_F First non-singleton cell.
 * \enumval IGRAPH_BLISS_FL First largest non-singleton cell.
 * \enumval IGRAPH_BLISS_FS First smallest non-singleton cell.
 * \enumval IGRAPH_BLISS_FM First maximally non-trivially connected
 *      non-singleton cell.
 * \enumval IGRAPH_BLISS_FLM Largest maximally non-trivially connected
 *      non-singleton cell.
 * \enumval IGRAPH_BLISS_FSM Smallest maximally non-trivially
 *      connected non-singletion cell.
 */

typedef enum { IGRAPH_BLISS_F=0, IGRAPH_BLISS_FL, 
	       IGRAPH_BLISS_FS, IGRAPH_BLISS_FM, 
	       IGRAPH_BLISS_FLM, IGRAPH_BLISS_FSM } igraph_bliss_sh_t;

int igraph_canonical_permutation(const igraph_t *graph, igraph_vector_t *labeling, 
				 igraph_bliss_sh_t sh, igraph_bliss_info_t *info);
int igraph_isomorphic_bliss(const igraph_t *graph1, const igraph_t *graph2,
			    igraph_bool_t *iso, igraph_vector_t *map12, 
			    igraph_vector_t *map21,
			    igraph_bliss_sh_t sh1, igraph_bliss_sh_t sh2, 
			    igraph_bliss_info_t *info1, igraph_bliss_info_t *info2);

int igraph_automorphisms(const igraph_t *graph,
			 igraph_bliss_sh_t sh, igraph_bliss_info_t *info);

/* Functions for 3-4 graphs */
int igraph_isomorphic_34(const igraph_t *graph1, const igraph_t *graph2, 
			 igraph_bool_t *iso);
int igraph_isoclass(const igraph_t *graph, igraph_integer_t *isoclass);
int igraph_isoclass_subgraph(const igraph_t *graph, igraph_vector_t *vids,
			     igraph_integer_t *isoclass);
int igraph_isoclass_create(igraph_t *graph, igraph_integer_t size,
			   igraph_integer_t number, igraph_bool_t directed);




__END_DECLS

#endif
