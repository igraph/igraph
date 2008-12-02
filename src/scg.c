/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2008  Gabor Csardi <Gabor.Csardi@unil.ch>
   University of Lausanne, Rue de Bugnon 27, CH-1005 Lausanne, Switzerland
   
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

#include "igraph.h"

int igraph_scg(const igraph_t *graph, igraph_integer_t ev,
	       igraph_scg_matrix_t matrix_type, 
	       igraph_scg_algorithm_t algo,
	       igraph_scg_norm_t norm,
	       igraph_scg_direction_t direction,
	       const igraph_matrix_t *evec,
	       const igraph_vector_t *markovp,
	       const igraph_vector_t *group,
	       /* igraph_bool_t use_arpack, */
	       igraph_integer_t maxiter
	       /* igraph_?_t *semproj, */
	       /* igraph_?_t *epairs, */
	       /* igraph_?_t *c_markovp, */) {
  /* TODO */
  return 0;
}

int igraph_scg_matrix(const igraph_matrix_t *matrix, igraph_integer_t ev,
		      igraph_scg_algorithm_t algo,
		      igraph_scg_norm_t norm,
		      igraph_scg_direction_t direction,
		      const igraph_matrix_t *evec,
		      const igraph_vector_t *markovp,
		      const igraph_vector_t *group,
		      /* igraph_bool_t use_arpack, */
		      igraph_integer_t maxiter
		      /* igraph_?_t *semproj, */
		      /* igraph_?_t *epairs, */
		      /* igraph_?_t *c_markovp, */) {
  /* TODO */
  return 0;
}

