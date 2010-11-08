/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2010  Gabor Csardi <csardi.gabor@gmail.com>
   Rue de l'Industrie 5, Lausanne 1005, Switzerland
   
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

#include "igraph_eigen.h"

/** 
 * \function igraph_eigen_matrix_symmetric
 * 
 */ 

int igraph_eigen_matrix_symmetric(const igraph_matrix_t *A,
				  const igraph_sparsemat_t *sA,
				  const igraph_arpack_function_t *fun, 
				  void *extra,
				  igraph_eigen_algorithm_t algorithm,
				  const igraph_eigen_which_t *which,
				  igraph_arpack_options_t *options,
				  igraph_arpack_storage_t *storage,
				  igraph_vector_t *values, 
				  igraph_matrix_t *vectors) {
  /* TODO */
  return 0;
}

/** 
 * \function igraph_eigen_matrix
 * 
 */ 

int igraph_eigen_matrix(const igraph_matrix_t *A,
			const igraph_sparsemat_t *sA,
			const igraph_arpack_function_t *fun,
			void *extra,
			igraph_eigen_algorithm_t algorithm,
			const igraph_eigen_which_t *which,
			igraph_arpack_options_t *options,
			igraph_arpack_storage_t *storage,
			igraph_vector_complex_t *values,
			igraph_matrix_complex_t *vectors) {
  /* TODO */
  return 0;
}

/** 
 * \function igraph_eigen_adjacency
 *
 */ 

int igraph_eigen_adjacency(const igraph_t *graph,
			   igraph_eigen_algorithm_t algorithm,
			   const igraph_eigen_which_t *which,
			   igraph_arpack_options_t *options,
			   igraph_arpack_storage_t *storage,
			   igraph_vector_t *values,
			   igraph_matrix_t *vectors,
			   igraph_vector_complex_t *cmplxvalues,
			   igraph_matrix_complex_t *cmplxvectors) {
  /* TODO */
  return 0;
}

/** 
 * \function igraph_eigen_laplacian
 *
 */ 

int igraph_eigen_laplacian(const igraph_t *graph,
			   igraph_eigen_algorithm_t algorithm,
			   const igraph_eigen_which_t *which,
			   igraph_arpack_options_t *options,
			   igraph_arpack_storage_t *storage,
			   igraph_vector_t *values,
			   igraph_matrix_t *vectors,
			   igraph_vector_complex_t *cmplxvalues,
			   igraph_matrix_complex_t *cmplxvectors) {
  /* TODO */
  return 0;
}
