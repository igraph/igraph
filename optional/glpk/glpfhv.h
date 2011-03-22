/* glpfhv.h (LP basis factorization, FHV eta file version) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*
*  Copyright (C) 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008,
*  2009, 2010 Andrew Makhorin, Department for Applied Informatics,
*  Moscow Aviation Institute, Moscow, Russia. All rights reserved.
*  E-mail: <mao@gnu.org>.
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

#ifndef GLPFHV_H
#define GLPFHV_H

#include "glpluf.h"

/***********************************************************************
*  The structure FHV defines the factorization of the basis mxm-matrix
*  B, where m is the number of rows in corresponding problem instance.
*
*  This factorization is the following sextet:
*
*     [B] = (F, H, V, P0, P, Q),                                     (1)
*
*  where F, H, and V are such matrices that
*
*     B = F * H * V,                                                 (2)
*
*  and P0, P, and Q are such permutation matrices that the matrix
*
*     L = P0 * F * inv(P0)                                           (3)
*
*  is lower triangular with unity diagonal, and the matrix
*
*     U = P * V * Q                                                  (4)
*
*  is upper triangular. All the matrices have the same order m, which
*  is the order of the basis matrix B.
*
*  The matrices F, V, P, and Q are stored in the structure LUF (see the
*  module GLPLUF), which is a member of the structure FHV.
*
*  The matrix H is stored in the form of eta file using row-like format
*  as follows:
*
*     H = H[1] * H[2] * ... * H[nfs],                                (5)
*
*  where H[k], k = 1, 2, ..., nfs, is a row-like factor, which differs
*  from the unity matrix only by one row, nfs is current number of row-
*  like factors. After the factorization has been built for some given
*  basis matrix B the matrix H has no factors and thus it is the unity
*  matrix. Then each time when the factorization is recomputed for an
*  adjacent basis matrix, the next factor H[k], k = 1, 2, ... is built
*  and added to the end of the eta file H.
*
*  Being sparse vectors non-trivial rows of the factors H[k] are stored
*  in the right part of the sparse vector area (SVA) in the same manner
*  as rows and columns of the matrix F.
*
*  For more details see the program documentation. */

typedef struct FHV FHV;

struct FHV
{     /* LP basis factorization */
      int m_max;
      /* maximal value of m (increased automatically, if necessary) */
      int m;
      /* the order of matrices B, F, H, V, P0, P, Q */
      int valid;
      /* the factorization is valid only if this flag is set */
      LUF *luf;
      /* LU-factorization (contains the matrices F, V, P, Q) */
      /*--------------------------------------------------------------*/
      /* matrix H in the form of eta file */
      int hh_max;
      /* maximal number of row-like factors (which limits the number of
         updates of the factorization) */
      int hh_nfs;
      /* current number of row-like factors (0 <= hh_nfs <= hh_max) */
      int *hh_ind; /* int hh_ind[1+hh_max]; */
      /* hh_ind[k], k = 1, ..., nfs, is the number of a non-trivial row
         of factor H[k] */
      int *hh_ptr; /* int hh_ptr[1+hh_max]; */
      /* hh_ptr[k], k = 1, ..., nfs, is a pointer to the first element
         of the non-trivial row of factor H[k] in the SVA */
      int *hh_len; /* int hh_len[1+hh_max]; */
      /* hh_len[k], k = 1, ..., nfs, is the number of non-zero elements
         in the non-trivial row of factor H[k] */
      /*--------------------------------------------------------------*/
      /* matrix P0 */
      int *p0_row; /* int p0_row[1+m_max]; */
      /* p0_row[i] = j means that p0[i,j] = 1 */
      int *p0_col; /* int p0_col[1+m_max]; */
      /* p0_col[j] = i means that p0[i,j] = 1 */
      /* if i-th row or column of the matrix F corresponds to i'-th row
         or column of the matrix L = P0*F*inv(P0), then p0_row[i'] = i
         and p0_col[i] = i' */
      /*--------------------------------------------------------------*/
      /* working arrays */
      int *cc_ind; /* int cc_ind[1+m_max]; */
      /* integer working array */
      double *cc_val; /* double cc_val[1+m_max]; */
      /* floating-point working array */
      /*--------------------------------------------------------------*/
      /* control parameters */
      double upd_tol;
      /* update tolerance; if after updating the factorization absolute
         value of some diagonal element u[k,k] of matrix U = P*V*Q is
         less than upd_tol * max(|u[k,*]|, |u[*,k]|), the factorization
         is considered as inaccurate */
      /*--------------------------------------------------------------*/
      /* some statistics */
      int nnz_h;
      /* current number of non-zeros in all factors of matrix H */
};

/* return codes: */
#define FHV_ESING    1  /* singular matrix */
#define FHV_ECOND    2  /* ill-conditioned matrix */
#define FHV_ECHECK   3  /* insufficient accuracy */
#define FHV_ELIMIT   4  /* update limit reached */
#define FHV_EROOM    5  /* SVA overflow */

#define fhv_create_it _glp_fhv_create_it
FHV *fhv_create_it(void);
/* create LP basis factorization */

#define fhv_factorize _glp_fhv_factorize
int fhv_factorize(FHV *fhv, int m, int (*col)(void *info, int j,
      int ind[], double val[]), void *info);
/* compute LP basis factorization */

#define fhv_h_solve _glp_fhv_h_solve
void fhv_h_solve(FHV *fhv, int tr, double x[]);
/* solve system H*x = b or H'*x = b */

#define fhv_ftran _glp_fhv_ftran
void fhv_ftran(FHV *fhv, double x[]);
/* perform forward transformation (solve system B*x = b) */

#define fhv_btran _glp_fhv_btran
void fhv_btran(FHV *fhv, double x[]);
/* perform backward transformation (solve system B'*x = b) */

#define fhv_update_it _glp_fhv_update_it
int fhv_update_it(FHV *fhv, int j, int len, const int ind[],
      const double val[]);
/* update LP basis factorization */

#define fhv_delete_it _glp_fhv_delete_it
void fhv_delete_it(FHV *fhv);
/* delete LP basis factorization */

#endif

/* eof */
