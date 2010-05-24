/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2007  Gabor Csardi <csardi@rmki.kfki.hu>
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
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
   02110-1301 USA

*/

#ifndef BLAS_H
#define BLAS_H

#include "igraph_types.h"
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

/**
 * \section about_blas About the BLAS interface in igraph
 * 
 * <para>
 * TODO
 * </para>
 * 
 * <para>
 * igraph does not contain all BLAS routines, only the ones
 * dealing with symmetric and non-symmetric matrices and vectors
 * using double precision real numbers.
 * </para>
 */

void igraph_blas_dgemv(igraph_bool_t transpose, igraph_real_t alpha,
        const igraph_matrix_t* a, const igraph_vector_t* x,
        igraph_real_t beta, igraph_vector_t* y);
void igraph_blas_dgemv_array(igraph_bool_t transpose, igraph_real_t alpha,
        const igraph_matrix_t* a, const igraph_real_t* x,
        igraph_real_t beta, igraph_real_t* y);

__END_DECLS

#endif
