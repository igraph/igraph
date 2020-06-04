/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2007-2012  Gabor Csardi <csardi.gabor@gmail.com>
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
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#ifndef IGRAPH_BLAS_H
#define IGRAPH_BLAS_H

#include "igraph_types.h"
#include "igraph_vector.h"
#include "igraph_matrix.h"
#include "igraph_decls.h"

__BEGIN_DECLS

/**
 * \section about_blas BLAS interface in igraph
 *
 * <para>
 * BLAS is a highly optimized library for basic linear algebra operations
 * such as vector-vector, matrix-vector and matrix-matrix product.
 * Please see http://www.netlib.org/blas/ for details and a reference
 * implementation in Fortran. igraph contains some wrapper functions
 * that can be used to call BLAS routines in a somewhat more
 * user-friendly way. Not all BLAS routines are included in igraph,
 * and even those which are included might not have wrappers;
 * the extension of the set of wrapped functions will probably be driven
 * by igraph's internal requirements. The wrapper functions usually
 * substitute double-precision floating point arrays used by BLAS with
 * \type igraph_vector_t and \type igraph_matrix_t instances and also
 * remove those parameters (such as the number of rows/columns) that
 * can be inferred from the passed arguments directly.
 * </para>
 */

DECLDIR void igraph_blas_dgemv(igraph_bool_t transpose, igraph_real_t alpha,
                               const igraph_matrix_t* a, const igraph_vector_t* x,
                               igraph_real_t beta, igraph_vector_t* y);
DECLDIR void igraph_blas_dgemv_array(igraph_bool_t transpose, igraph_real_t alpha,
                                     const igraph_matrix_t* a, const igraph_real_t* x,
                                     igraph_real_t beta, igraph_real_t* y);

DECLDIR igraph_real_t igraph_blas_dnrm2(const igraph_vector_t *v);

__END_DECLS

#endif
