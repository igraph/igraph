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

#ifndef IGRAPH_CONVERSION_H
#define IGRAPH_CONVERSION_H

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
#include "igraph_spmatrix.h"
#include "igraph_matrix.h"
#include "igraph_sparsemat.h"
#include "igraph_attributes.h"

__BEGIN_DECLS

/* -------------------------------------------------- */
/* Conversion                                         */
/* -------------------------------------------------- */

int igraph_get_adjacency(const igraph_t *graph, igraph_matrix_t *res,
			 igraph_get_adjacency_t type, igraph_bool_t eids);
int igraph_get_adjacency_sparse(const igraph_t *graph, igraph_spmatrix_t *res,
			        igraph_get_adjacency_t type);

int igraph_get_stochastic(const igraph_t *graph, 
			  igraph_matrix_t *matrix,
			  igraph_bool_t column_wise);

int igraph_get_stochastic_sparsemat(const igraph_t *graph, 
				    igraph_sparsemat_t *sparsemat,
				    igraph_bool_t column_wise);

int igraph_get_edgelist(const igraph_t *graph, igraph_vector_t *res, igraph_bool_t bycol);

int igraph_to_directed(igraph_t *graph, 
		       igraph_to_directed_t flags);
int igraph_to_undirected(igraph_t *graph,
			 igraph_to_undirected_t flags,
			 const igraph_attribute_combination_t *edge_comb);

__END_DECLS

#endif
