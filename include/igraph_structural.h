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

#ifndef IGRAPH_STRUCTURAL_H
#define IGRAPH_STRUCTURAL_H

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
#include "igraph_iterators.h"

__BEGIN_DECLS

/* -------------------------------------------------- */
/* Basic query functions                              */
/* -------------------------------------------------- */

int igraph_are_connected(const igraph_t *graph, igraph_integer_t v1, igraph_integer_t v2, igraph_bool_t *res);

/* -------------------------------------------------- */
/* Structural properties                              */
/* -------------------------------------------------- */

int igraph_minimum_spanning_tree_unweighted(const igraph_t *graph, 
					    igraph_t *mst);
int igraph_minimum_spanning_tree_prim(const igraph_t *graph, igraph_t *mst,
				      const igraph_vector_t *weights);
int igraph_subcomponent(const igraph_t *graph, igraph_vector_t *res, igraph_real_t vid, 
			igraph_neimode_t mode);	
int igraph_rewire(igraph_t *graph, igraph_integer_t n, igraph_rewiring_t mode);
int igraph_subgraph(const igraph_t *graph, igraph_t *res, 
		    const igraph_vs_t vids);
int igraph_simplify(igraph_t *graph, igraph_bool_t multiple, igraph_bool_t loops);
int igraph_reciprocity(const igraph_t *graph, igraph_real_t *res,
		       igraph_bool_t ignore_loops);

int igraph_maxdegree(const igraph_t *graph, igraph_integer_t *res,
		     igraph_vs_t vids, igraph_neimode_t mode, 
		     igraph_bool_t loops);
int igraph_density(const igraph_t *graph, igraph_real_t *res, 
		   igraph_bool_t loops);

int igraph_is_loop(const igraph_t *graph, igraph_vector_bool_t *res, 
		   igraph_es_t es);
int igraph_is_simple(const igraph_t *graph, igraph_bool_t *res);
int igraph_is_multiple(const igraph_t *graph, igraph_vector_bool_t *res, 
		       igraph_es_t es);
int igraph_count_multiple(const igraph_t *graph, igraph_vector_t *res, igraph_es_t es);
int igraph_girth(const igraph_t *graph, igraph_integer_t *girth, 
		 igraph_vector_t *circle);
int igraph_add_edge(igraph_t *graph, igraph_integer_t from, igraph_integer_t to);

int igraph_unfold_tree(const igraph_t *graph, igraph_t *tree,
		       igraph_neimode_t mode, const igraph_vector_t *roots,
		       igraph_vector_t *vertex_index);

int igraph_is_mutual(igraph_t *graph, igraph_vector_bool_t *res, igraph_es_t es);

int igraph_maximum_cardinality_search(const igraph_t *graph,
				      igraph_vector_t *alpha,
				      igraph_vector_t *alpham1);
int igraph_is_chordal(const igraph_t *graph,
		      const igraph_vector_t *alpha,
		      const igraph_vector_t *alpham1,
		      igraph_bool_t *chordal,
		      igraph_vector_t *fill_in,
		      igraph_t *newgraph);
int igraph_avg_nearest_neighbor_degree(const igraph_t *graph,
				       igraph_vs_t vids,
				       igraph_vector_t *knn,
				       igraph_vector_t *knnk, 
				       const igraph_vector_t *weights);

/* -------------------------------------------------- */
/* Spectral Properties                                */
/* -------------------------------------------------- */

int igraph_laplacian(const igraph_t *graph, igraph_matrix_t *res,
		     igraph_bool_t normalized, 
		     const igraph_vector_t *weights);


__END_DECLS

#endif
