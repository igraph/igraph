/* spxat.h */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2015 Free Software Foundation, Inc.
*  Written by Andrew Makhorin <mao@gnu.org>.
*
*  GLPK is free software: you can redistribute it and/or modify it
*  under the terms of the GNU General Public License as published by
*  the Free Software Foundation, either version 3 of the License, or
*  (at your option) any later version.
*
*  GLPK is distributed in the hope that it will be useful, but WITHOUT
*  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
*  or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public
*  License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with GLPK. If not, see <http://www.gnu.org/licenses/>.
***********************************************************************/

#ifndef SPXAT_H
#define SPXAT_H

#include "spxlp.h"

typedef struct SPXAT SPXAT;

struct SPXAT
{     /* mxn-matrix A of constraint coefficients in sparse row-wise
       * format */
      int *ptr; /* int ptr[1+m+1]; */
      /* ptr[0] is not used;
       * ptr[i], 1 <= i <= m, is starting position of i-th row in
       * arrays ind and val; note that ptr[1] is always 1;
       * ptr[m+1] indicates the position after the last element in
       * arrays ind and val, i.e. ptr[m+1] = nnz+1, where nnz is the
       * number of non-zero elements in matrix A;
       * the length of i-th row (the number of non-zero elements in
       * that row) can be calculated as ptr[i+1] - ptr[i] */
      int *ind; /* int ind[1+nnz]; */
      /* column indices */
      double *val; /* double val[1+nnz]; */
      /* non-zero element values */
      double *work; /* double work[1+n]; */
      /* working array */
};

#define spx_alloc_at _glp_spx_alloc_at
void spx_alloc_at(SPXLP *lp, SPXAT *at);
/* allocate constraint matrix in sparse row-wise format */

#define spx_build_at _glp_spx_build_at
void spx_build_at(SPXLP *lp, SPXAT *at);
/* build constraint matrix in sparse row-wise format */

#define spx_at_prod _glp_spx_at_prod
void spx_at_prod(SPXLP *lp, SPXAT *at, double y[/*1+n*/], double s,
      const double x[/*1+m*/]);
/* compute product y := y + s * A'* x */

#define spx_nt_prod1 _glp_spx_nt_prod1
void spx_nt_prod1(SPXLP *lp, SPXAT *at, double y[/*1+n-m*/], int ign,
      double s, const double x[/*1+m*/]);
/* compute product y := y + s * N'* x */

#define spx_eval_trow1 _glp_spx_eval_trow1
void spx_eval_trow1(SPXLP *lp, SPXAT *at, const double rho[/*1+m*/],
      double trow[/*1+n-m*/]);
/* compute i-th row of simplex table */

#define spx_free_at _glp_spx_free_at
void spx_free_at(SPXLP *lp, SPXAT *at);
/* deallocate constraint matrix in sparse row-wise format */

#endif

/* eof */
