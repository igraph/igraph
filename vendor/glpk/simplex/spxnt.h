/* spxnt.h */

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

#ifndef SPXNT_H
#define SPXNT_H

#include "spxlp.h"

typedef struct SPXNT SPXNT;

struct SPXNT
{     /* mx(n-m)-matrix N composed of non-basic columns of constraint
       * matrix A, in sparse row-wise format */
      int *ptr; /* int ptr[1+m]; */
      /* ptr[0] is not used;
       * ptr[i], 1 <= i <= m, is starting position of i-th row in
       * arrays ind and val; note that ptr[1] is always 1;
       * these starting positions are set up *once* as if they would
       * correspond to rows of matrix A stored without gaps, i.e.
       * ptr[i+1] - ptr[i] is the number of non-zeros in i-th (i < m)
       * row of matrix A, and (nnz+1) - ptr[m] is the number of
       * non-zero in m-th (last) row of matrix A, where nnz is the
       * total number of non-zeros in matrix A */
      int *len; /* int len[1+m]; */
      /* len[0] is not used;
       * len[i], 1 <= i <= m, is the number of non-zeros in i-th row
       * of current matrix N */
      int *ind; /* int ind[1+nnz]; */
      /* column indices */
      double *val; /* double val[1+nnz]; */
      /* non-zero element values */
};

#define spx_alloc_nt _glp_spx_alloc_nt
void spx_alloc_nt(SPXLP *lp, SPXNT *nt);
/* allocate matrix N in sparse row-wise format */

#define spx_init_nt _glp_spx_init_nt
void spx_init_nt(SPXLP *lp, SPXNT *nt);
/* initialize row pointers for matrix N */

#define spx_nt_add_col _glp_spx_nt_add_col
void spx_nt_add_col(SPXLP *lp, SPXNT *nt, int j, int k);
/* add column N[j] = A[k] */

#define spx_build_nt _glp_spx_build_nt
void spx_build_nt(SPXLP *lp, SPXNT *nt);
/* build matrix N for current basis */

#define spx_nt_del_col _glp_spx_nt_del_col
void spx_nt_del_col(SPXLP *lp, SPXNT *nt, int j, int k);
/* remove column N[j] = A[k] from matrix N */

#define spx_update_nt _glp_spx_update_nt
void spx_update_nt(SPXLP *lp, SPXNT *nt, int p, int q);
/* update matrix N for adjacent basis */

#define spx_nt_prod _glp_spx_nt_prod
void spx_nt_prod(SPXLP *lp, SPXNT *nt, double y[/*1+n-m*/], int ign,
      double s, const double x[/*1+m*/]);
/* compute product y := y + s * N'* x */

#if 1 /* 31/III-2016 */
#define spx_nt_prod_s _glp_spx_nt_prod_s
void spx_nt_prod_s(SPXLP *lp, SPXNT *nt, FVS *y, int ign, double s,
      const FVS *x, double eps);
/* sparse version of spx_nt_prod */
#endif

#define spx_free_nt _glp_spx_free_nt
void spx_free_nt(SPXLP *lp, SPXNT *nt);
/* deallocate matrix N in sparse row-wise format */

#endif

/* eof */
