/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
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

#ifndef BLAS_INTERNAL_H
#define BLAS_INTERNAL_H

/* Note: only files calling the BLAS routines directly need to
   include this header.
*/

#include "igraph_types.h"
#include "config.h"

#ifndef INTERNAL_BLAS
    #define igraphdaxpy_    daxpy_
    #define igraphdger_     dger_
    #define igraphdcopy_    dcopy_
    #define igraphdscal_    dscal_
    #define igraphdswap_    dswap_
    #define igraphdgemm_    dgemm_
    #define igraphdgemv_    dgemv_
    #define igraphddot_     ddot_
    #define igraphdnrm2_    dnrm2_
    #define igraphlsame_    lsame_
    #define igraphdrot_     drot_
    #define igraphidamax_   idamax_
    #define igraphdtrmm_    dtrmm_
    #define igraphdasum_    dasum_
    #define igraphdtrsm_    dtrsm_
    #define igraphdtrsv_    dtrsv_
    #define igraphdnrm2_    dnrm2_
    #define igraphdsymv_    dsymv_
    #define igraphdsyr2_    dsyr2_
    #define igraphdsyr2k_   dsyr2k_
    #define igraphdtrmv_    dtrmv_
    #define igraphdsyrk_    dsyrk_
#endif

int igraphdgemv_(char *trans, int *m, int *n, igraph_real_t *alpha,
                 igraph_real_t *a, int *lda, igraph_real_t *x, int *incx,
                 igraph_real_t *beta, igraph_real_t *y, int *incy);

int igraphdgemm_(char *transa, char *transb, int *m, int *n, int *k,
                 double *alpha, double *a, int *lda, double *b, int *ldb,
                 double *beta, double *c__, int *ldc);

double igraphdnrm2_(int *n, double *x, int *incx);

double igraphddot_(int *n, double *dx, int *incx, double *dy, int *incy);

#endif
