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

#ifndef LAPACK_H
#define LAPACK_H

#include "igraph_vector.h"
#include "igraph_matrix.h"

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

int igraph_lapack_dgetrf(igraph_matrix_t *a, igraph_vector_int_t *ipiv,
			 int *info);
int igraph_lapack_dgetrs(igraph_bool_t transpose, const igraph_matrix_t *a,
			 igraph_vector_int_t *ipiv, igraph_matrix_t *b,
			 int *info);
int igraph_lapack_dgesv(igraph_matrix_t *a, igraph_vector_int_t *ipiv,
			igraph_matrix_t *b, int *info);

typedef enum { IGRAPH_LAPACK_DSYEV_ALL,
	       IGRAPH_LAPACK_DSYEV_INTERVAL,
	       IGRAPH_LAPACK_DSYEV_SELECT } igraph_lapack_dsyev_which_t;

int igraph_lapack_dsyevr(const igraph_matrix_t *A, 
			 igraph_lapack_dsyev_which_t which,
			 igraph_real_t vl, igraph_real_t vu, int vestimate, 
			 int il, int iu, igraph_real_t abstol,
			 igraph_vector_t *values, igraph_matrix_t *vectors,
			 igraph_vector_int_t *support);

int igraph_lapack_dgeev(igraph_bool_t leftvec, igraph_bool_t rightvec,
			const igraph_matrix_t *A, 
			igraph_vector_t *valuesreal,
			igraph_vector_t *valuesimag, 
			igraph_matrix_t *vectorsleft,
			igraph_matrix_t *vectorsright, int *info);
			
__END_DECLS

#endif
