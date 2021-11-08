/* fhv.h (sparse updatable FHV-factorization) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2012-2013 Free Software Foundation, Inc.
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

#ifndef FHV_H
#define FHV_H

#include "luf.h"

/***********************************************************************
*  The structure FHV describes sparse updatable FHV-factorization.
*
*  The FHV-factorization has the following format:
*
*     A = F * H * V,                                                 (1)
*
*     F = P0 * L * P0',                                              (2)
*
*     H = H[1] * H[2] * ... * H[nfs],                                (3)
*
*     V = P * U * Q,                                                 (4)
*
*  where: A is a given (unsymmetric) square matrix; F, H, V are matrix
*  factors actually computed; L is a lower triangular matrix with unity
*  diagonal; U is an upper tringular matrix; H[k], k = 1, 2, ..., nfs,
*  is a row-like factor, which differs from unity matrix only in one
*  row called a non-trivial row; P0, P, Q are permutation matrices; and
*  P0' is a matrix transposed to P0.
*
*  Matrices F, V, P, Q are stored in the underlying LUF object.
*
*  Non-trivial rows of factors H[k] are stored as sparse vectors in the
*  right (static) part of the sparse vector area (SVA). Note that unity
*  diagonal elements of non-trivial rows are not stored.
*
*  Matrix P0 is stored in the same way as matrix P.
*
*  Matrices L and U are completely defined by matrices F, V, P, and Q,
*  and therefore not stored explicitly. */

typedef struct FHV FHV;

struct FHV
{     /* FHV-factorization */
      LUF *luf;
      /* LU-factorization (contains matrices F, V, P, Q) */
      /*--------------------------------------------------------------*/
      /* matrix H in the form of eta file */
      int nfs_max;
      /* maximal number of row-like factors (this limits the number of
       * updates of the factorization) */
      int nfs;
      /* current number of row-like factors, 0 <= nfs <= nfs_max */
      int *hh_ind; /* int hh_ind[1+nfs_max]; */
      /* hh_ind[0] is not used;
       * hh_ind[k], 1 <= k <= nfs, is number of non-trivial row of
       * factor H[k] */
      int hh_ref;
      /* reference number of sparse vector in SVA, which is non-trivial
       * row of factor H[1] */
#if 0 + 0
      int *hh_ptr = &sva->ptr[hh_ref-1];
      /* hh_ptr[0] is not used;
       * hh_ptr[k], 1 <= k <= nfs, is pointer to non-trivial row of
       * factor H[k] */
      int *hh_len = &sva->len[hh_ref-1];
      /* hh_len[0] is not used;
       * hh_len[k], 1 <= k <= nfs, is number of non-zero elements in
       * non-trivial row of factor H[k] */
#endif
      /*--------------------------------------------------------------*/
      /* matrix P0 */
      int *p0_ind; /* int p0_ind[1+n]; */
      /* p0_ind[i] = j means that P0[i,j] = 1 */
      int *p0_inv; /* int p0_inv[1+n]; */
      /* p0_inv[j] = i means that P0[i,j] = 1 */
};

#define fhv_ft_update _glp_fhv_ft_update
int fhv_ft_update(FHV *fhv, int q, int aq_len, const int aq_ind[],
      const double aq_val[], int ind[/*1+n*/], double val[/*1+n*/],
      double work[/*1+n*/]);
/* update FHV-factorization (Forrest-Tomlin) */

#define fhv_h_solve _glp_fhv_h_solve
void fhv_h_solve(FHV *fhv, double x[/*1+n*/]);
/* solve system H * x = b */

#define fhv_ht_solve _glp_fhv_ht_solve
void fhv_ht_solve(FHV *fhv, double x[/*1+n*/]);
/* solve system H' * x = b */

#endif

/* eof */
