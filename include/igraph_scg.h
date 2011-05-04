/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2011  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard st, Cambridge, MA, 02138 USA
   
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
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
   02110-1301 USA

*/

#ifndef IGRAPH_SCG_H
#define IGRAPH_SCG_H

#include "igraph_types.h"
#include "igraph_vector.h"
#include "igraph_matrix.h"
#include "igraph_sparsemat.h"

typedef enum { IGRAPH_SCG_SYMMETRIC=1, IGRAPH_SCG_LAPLACIAN=2,
	       IGRAPH_SCG_STOCHASTIC=3 } igraph_scg_matrix_t;

typedef enum { IGRAPH_SCG_OPTIMUM=1, IGRAPH_SCG_INTERV_KM=2,
	       IGRAPH_SCG_INTERV=3, IGRAPH_SCG_EXACT=4 } 
  igraph_scg_algorithm_t;

typedef enum { IGRAPH_SCG_NORM_ROW=1, IGRAPH_SCG_NORM_COL=2 }
  igraph_scg_norm_t;

typedef enum { IGRAPH_SCG_DIRECTION_DEFAULT=1, 
	       IGRAPH_SCG_DIRECTION_LEFT=2,
	       IGRAPH_SCG_DIRECTION_RIGHT=3 } igraph_scg_direction_t;

int igraph_scg_grouping(const igraph_matrix_t *V, 
			igraph_vector_t *groups,
			igraph_integer_t intervals,
			const igraph_vector_t *intervals_vector,
			igraph_scg_matrix_t matrix_type,
			igraph_scg_algorithm_t algorithm,
			const igraph_vector_t *p,
			igraph_integer_t maxiter);

int igraph_scg_semiprojectors(const igraph_vector_t *groups,
			      igraph_scg_matrix_t matrix_type,
			      igraph_matrix_t *L,
			      igraph_matrix_t *R,
			      igraph_sparsemat_t *Lsparse,
			      igraph_sparsemat_t *Rsparse, 
			      const igraph_vector_t *p,
			      igraph_scg_norm_t norm);

int igraph_scg_norm_eps(const igraph_matrix_t *V,
			const igraph_vector_t *groups,
			igraph_vector_t *eps,
			igraph_scg_matrix_t matrix_type,
			const igraph_vector_t *p,
			igraph_scg_norm_t norm);

int igraph_scg(const igraph_t *graph,
	       const igraph_vector_t *ev,
	       igraph_integer_t intervals, 
	       const igraph_vector_t *intervals_vector,
	       igraph_scg_matrix_t matrix_type,
	       igraph_scg_algorithm_t algorithm,
	       igraph_scg_norm_t norm, 
	       igraph_scg_direction_t direction,
	       const igraph_matrix_t *evec,
	       const igraph_matrix_complex_t *evec_cplx,
	       const igraph_vector_t *given_p,
	       const igraph_vector_t *given_groups,
	       igraph_bool_t use_arpack,
	       igraph_integer_t maxiter,
	       igraph_t *scg_graph,
	       igraph_matrix_t *scg_matrix,
	       igraph_sparsemat_t *scg_sparsemat,
	       igraph_matrix_t *L,
	       igraph_matrix_t *R,
	       igraph_sparsemat_t *Lsparse,
	       igraph_sparsemat_t *Rsparse,
	       igraph_vector_t *eigenvalues,
	       igraph_vector_complex_t *eigenvalues_cplx,
	       igraph_matrix_t *eigenvectors,
	       igraph_matrix_complex_t *eigenvectors_cplx,
	       igraph_vector_t *groups,
	       igraph_vector_t *p);

int igraph_scg_matrix(const igraph_matrix_t *matrix,
		      const igraph_sparsemat_t *sparsemat,
		      const igraph_vector_t *ev,
		      igraph_integer_t intervals, 
		      const igraph_vector_t *intervals_vector,
		      igraph_scg_matrix_t matrix_type,
		      igraph_scg_algorithm_t algorithm,
		      igraph_scg_norm_t norm, 
		      igraph_scg_direction_t direction,
		      const igraph_matrix_t *evec,
		      const igraph_matrix_complex_t *evec_cplx,
		      const igraph_vector_t *given_p,
		      const igraph_vector_t *given_groups,
		      igraph_bool_t use_arpack,
		      igraph_integer_t maxiter,
		      igraph_t *scg_graph,
		      igraph_matrix_t *scg_matrix,
		      igraph_sparsemat_t *scg_sparsemat,
		      igraph_matrix_t *L,
		      igraph_matrix_t *R,
		      igraph_sparsemat_t *Lsparse,
		      igraph_sparsemat_t *Rsparse,
		      igraph_vector_t *eigenvalues,
		      igraph_vector_complex_t *eigenvalues_cplx,
		      igraph_matrix_t *eigenvectors,
		      igraph_matrix_complex_t *eigenvectors_cplx,
		      igraph_vector_t *groups,
		      igraph_vector_t *p);

#endif
