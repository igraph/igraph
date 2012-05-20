/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2003-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#ifndef IGRAPH_FLOW_H
#define IGRAPH_FLOW_H

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
#include "igraph_types.h"
#include "igraph_datatype.h"
#include "igraph_vector_ptr.h"

__BEGIN_DECLS

/* -------------------------------------------------- */
/* MAximum flows, minimum cuts & such                 */
/* -------------------------------------------------- */

int igraph_maxflow(const igraph_t *graph, igraph_real_t *value,
		   igraph_vector_t *flow, igraph_vector_t *cut,
		   igraph_vector_t *partition, igraph_vector_t *partition2,
		   igraph_integer_t source, igraph_integer_t target,
		   const igraph_vector_t *capacity);
int igraph_maxflow_value(const igraph_t *graph, igraph_real_t *value,
			 igraph_integer_t source, igraph_integer_t target,
			 const igraph_vector_t *capacity);

int igraph_st_mincut(const igraph_t *graph, igraph_real_t *value,
		     igraph_vector_t *cut, igraph_vector_t *partition,
		     igraph_vector_t *partition2,
		     igraph_integer_t source, igraph_integer_t target,
		     const igraph_vector_t *capacity);
int igraph_st_mincut_value(const igraph_t *graph, igraph_real_t *res,
                           igraph_integer_t source, igraph_integer_t target,
			   const igraph_vector_t *capacity);

int igraph_mincut_value(const igraph_t *graph, igraph_real_t *res, 
			const igraph_vector_t *capacity);
int igraph_mincut(const igraph_t *graph,
		  igraph_real_t *value,
		  igraph_vector_t *partition,
		  igraph_vector_t *partition2,
		  igraph_vector_t *cut,
		  const igraph_vector_t *capacity);

int igraph_st_vertex_connectivity(const igraph_t *graph, 
				  igraph_integer_t *res,
				  igraph_integer_t source,
				  igraph_integer_t target,
				  igraph_vconn_nei_t neighbors);
int igraph_vertex_connectivity(const igraph_t *graph, igraph_integer_t *res,
			       igraph_bool_t checks);

int igraph_st_edge_connectivity(const igraph_t *graph, igraph_integer_t *res,
				igraph_integer_t source, 
				igraph_integer_t target);
int igraph_edge_connectivity(const igraph_t *graph, igraph_integer_t *res,
			     igraph_bool_t checks);

int igraph_edge_disjoint_paths(const igraph_t *graph, igraph_integer_t *res,
			       igraph_integer_t source, 
			       igraph_integer_t target);
int igraph_vertex_disjoint_paths(const igraph_t *graph, igraph_integer_t *res,
				 igraph_integer_t source,
				 igraph_integer_t target);

int igraph_adhesion(const igraph_t *graph, igraph_integer_t *res,
		    igraph_bool_t checks);
int igraph_cohesion(const igraph_t *graph, igraph_integer_t *res,
		    igraph_bool_t checks);

/* s-t cut listing related stuff */

int igraph_even_tarjan_reduction(const igraph_t *graph, igraph_t *graphbar,
				 igraph_vector_t *capacity);

int igraph_residual_graph(const igraph_t *graph,
			  const igraph_vector_t *capacity,
			  igraph_t *residual,
			  igraph_vector_t *residual_capacity,
			  const igraph_vector_t *flow);
int igraph_i_residual_graph(const igraph_t *graph,
			    const igraph_vector_t *capacity,
			    igraph_t *residual,
			    igraph_vector_t *residual_capacity,
			    const igraph_vector_t *flow, 
			    igraph_vector_t *tmp);

int igraph_i_reverse_residual_graph(const igraph_t *graph,
				    const igraph_vector_t *capacity,
				    igraph_t *residual,
				    const igraph_vector_t *flow,
				    igraph_vector_t *tmp);
int igraph_reverse_residual_graph(const igraph_t *graph,
				  const igraph_vector_t *capacity,
				  igraph_t *residual,
				  const igraph_vector_t *flow);

int igraph_dominator_tree(const igraph_t *graph,
			  igraph_integer_t root,
			  igraph_vector_t *dom,
			  igraph_t *domtree,
			  igraph_vector_t *leftout,
			  igraph_neimode_t mode);

int igraph_all_st_cuts(const igraph_t *graph,
		       igraph_vector_ptr_t *cuts,
		       igraph_vector_ptr_t *partition1s,
		       igraph_integer_t source,
		       igraph_integer_t target);

int igraph_all_st_mincuts(const igraph_t *graph, igraph_real_t *value,
			  igraph_vector_ptr_t *cuts,
			  igraph_vector_ptr_t *partition1s,
			  igraph_integer_t source,
			  igraph_integer_t target,
			  const igraph_vector_t *capacity);

__END_DECLS

#endif
