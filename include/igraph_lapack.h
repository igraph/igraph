/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2010-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#ifndef IGRAPH_LAPACK_H
#define IGRAPH_LAPACK_H

#include "igraph_decls.h"
#include "igraph_vector.h"
#include "igraph_matrix.h"

__BEGIN_DECLS

/**
 * \section about_lapack LAPACK interface in igraph
 *
 * <para>
 * LAPACK is written in Fortran90 and provides routines for solving
 * systems of simultaneous linear equations, least-squares solutions
 * of linear systems of equations, eigenvalue problems, and singular
 * value problems. The associated matrix factorizations (LU, Cholesky,
 * QR, SVD, Schur, generalized Schur) are also provided, as are
 * related computations such as reordering of the Schur factorizations
 * and estimating condition numbers. Dense and banded matrices are
 * handled, but not general sparse matrices. In all areas, similar
 * functionality is provided for real and complex matrices, in both
 * single and double precision.
 * </para>
 *
 * <para>
 * igraph provides an interface to a very limited set of LAPACK
 * functions, using the regular igraph data structures.
 * </para>
 *
 * <para>
 * See more about LAPACK at http://www.netlib.org/lapack/
 * </para>
 */

IGRAPH_EXPORT int igraph_lapack_dgetrf(igraph_matrix_t *a, igraph_vector_int_t *ipiv,
                                       int *info);
IGRAPH_EXPORT int igraph_lapack_dgetrs(igraph_bool_t transpose, const igraph_matrix_t *a,
                                       const igraph_vector_int_t *ipiv, igraph_matrix_t *b);
IGRAPH_EXPORT int igraph_lapack_dgesv(igraph_matrix_t *a, igraph_vector_int_t *ipiv,
                                      igraph_matrix_t *b, int *info);

typedef enum { IGRAPH_LAPACK_DSYEV_ALL,
               IGRAPH_LAPACK_DSYEV_INTERVAL,
               IGRAPH_LAPACK_DSYEV_SELECT
             } igraph_lapack_dsyev_which_t;

IGRAPH_EXPORT int igraph_lapack_dsyevr(const igraph_matrix_t *A,
                                       igraph_lapack_dsyev_which_t which,
                                       igraph_real_t vl, igraph_real_t vu, int vestimate,
                                       int il, int iu, igraph_real_t abstol,
                                       igraph_vector_t *values, igraph_matrix_t *vectors,
                                       igraph_vector_int_t *support);

/* TODO: should we use complex vectors/matrices? */

IGRAPH_EXPORT int igraph_lapack_dgeev(const igraph_matrix_t *A,
                                      igraph_vector_t *valuesreal,
                                      igraph_vector_t *valuesimag,
                                      igraph_matrix_t *vectorsleft,
                                      igraph_matrix_t *vectorsright, int *info);

typedef enum { IGRAPH_LAPACK_DGEEVX_BALANCE_NONE = 0,
               IGRAPH_LAPACK_DGEEVX_BALANCE_PERM,
               IGRAPH_LAPACK_DGEEVX_BALANCE_SCALE,
               IGRAPH_LAPACK_DGEEVX_BALANCE_BOTH
             }
igraph_lapack_dgeevx_balance_t;

IGRAPH_EXPORT int igraph_lapack_dgeevx(igraph_lapack_dgeevx_balance_t balance,
                                       const igraph_matrix_t *A,
                                       igraph_vector_t *valuesreal,
                                       igraph_vector_t *valuesimag,
                                       igraph_matrix_t *vectorsleft,
                                       igraph_matrix_t *vectorsright,
                                       int *ilo, int *ihi, igraph_vector_t *scale,
                                       igraph_real_t *abnrm,
                                       igraph_vector_t *rconde,
                                       igraph_vector_t *rcondv,
                                       int *info);

IGRAPH_EXPORT int igraph_lapack_dgehrd(const igraph_matrix_t *A,
                                       int ilo, int ihi,
                                       igraph_matrix_t *result);

__END_DECLS

#endif
