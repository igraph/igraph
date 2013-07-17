/* glpfhv.c (LP basis factorization, FHV eta file version) */

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

#ifdef __clang__
#pragma clang diagnostic ignored "-Wsign-conversion"
#endif

#include "glpfhv.h"
#include "glpenv.h"
#define xfault xerror

/* CAUTION: DO NOT CHANGE THE LIMIT BELOW */

#define M_MAX 100000000 /* = 100*10^6 */
/* maximal order of the basis matrix */

/***********************************************************************
*  NAME
*
*  fhv_create_it - create LP basis factorization
*
*  SYNOPSIS
*
*  #include "glpfhv.h"
*  FHV *fhv_create_it(void);
*
*  DESCRIPTION
*
*  The routine fhv_create_it creates a program object, which represents
*  a factorization of LP basis.
*
*  RETURNS
*
*  The routine fhv_create_it returns a pointer to the object created. */

FHV *fhv_create_it(void)
{     FHV *fhv;
      fhv = xmalloc(sizeof(FHV));
      fhv->m_max = fhv->m = 0;
      fhv->valid = 0;
      fhv->luf = luf_create_it();
      fhv->hh_max = 50;
      fhv->hh_nfs = 0;
      fhv->hh_ind = fhv->hh_ptr = fhv->hh_len = NULL;
      fhv->p0_row = fhv->p0_col = NULL;
      fhv->cc_ind = NULL;
      fhv->cc_val = NULL;
      fhv->upd_tol = 1e-6;
      fhv->nnz_h = 0;
      return fhv;
}

/***********************************************************************
*  NAME
*
*  fhv_factorize - compute LP basis factorization
*
*  SYNOPSIS
*
*  #include "glpfhv.h"
*  int fhv_factorize(FHV *fhv, int m, int (*col)(void *info, int j,
*     int ind[], double val[]), void *info);
*
*  DESCRIPTION
*
*  The routine fhv_factorize computes the factorization of the basis
*  matrix B specified by the routine col.
*
*  The parameter fhv specified the basis factorization data structure
*  created by the routine fhv_create_it.
*
*  The parameter m specifies the order of B, m > 0.
*
*  The formal routine col specifies the matrix B to be factorized. To
*  obtain j-th column of A the routine fhv_factorize calls the routine
*  col with the parameter j (1 <= j <= n). In response the routine col
*  should store row indices and numerical values of non-zero elements
*  of j-th column of B to locations ind[1,...,len] and val[1,...,len],
*  respectively, where len is the number of non-zeros in j-th column
*  returned on exit. Neither zero nor duplicate elements are allowed.
*
*  The parameter info is a transit pointer passed to the routine col.
*
*  RETURNS
*
*  0  The factorization has been successfully computed.
*
*  FHV_ESING
*     The specified matrix is singular within the working precision.
*
*  FHV_ECOND
*     The specified matrix is ill-conditioned.
*
*  For more details see comments to the routine luf_factorize.
*
*  ALGORITHM
*
*  The routine fhv_factorize calls the routine luf_factorize (see the
*  module GLPLUF), which actually computes LU-factorization of the basis
*  matrix B in the form
*
*     [B] = (F, V, P, Q),
*
*  where F and V are such matrices that
*
*     B = F * V,
*
*  and P and Q are such permutation matrices that the matrix
*
*     L = P * F * inv(P)
*
*  is lower triangular with unity diagonal, and the matrix
*
*     U = P * V * Q
*
*  is upper triangular.
*
*  In order to build the complete representation of the factorization
*  (see formula (1) in the file glpfhv.h) the routine fhv_factorize just
*  additionally sets H = I and P0 = P. */

int fhv_factorize(FHV *fhv, int m, int (*col)(void *info, int j,
      int ind[], double val[]), void *info)
{     int ret;
      if (m < 1)
         xfault("fhv_factorize: m = %d; invalid parameter\n", m);
      if (m > M_MAX)
         xfault("fhv_factorize: m = %d; matrix too big\n", m);
      fhv->m = m;
      /* invalidate the factorization */
      fhv->valid = 0;
      /* allocate/reallocate arrays, if necessary */
      if (fhv->hh_ind == NULL)
         fhv->hh_ind = xcalloc(1+fhv->hh_max, sizeof(int));
      if (fhv->hh_ptr == NULL)
         fhv->hh_ptr = xcalloc(1+fhv->hh_max, sizeof(int));
      if (fhv->hh_len == NULL)
         fhv->hh_len = xcalloc(1+fhv->hh_max, sizeof(int));
      if (fhv->m_max < m)
      {  if (fhv->p0_row != NULL) xfree(fhv->p0_row);
         if (fhv->p0_col != NULL) xfree(fhv->p0_col);
         if (fhv->cc_ind != NULL) xfree(fhv->cc_ind);
         if (fhv->cc_val != NULL) xfree(fhv->cc_val);
         fhv->m_max = m + 100;
         fhv->p0_row = xcalloc(1+fhv->m_max, sizeof(int));
         fhv->p0_col = xcalloc(1+fhv->m_max, sizeof(int));
         fhv->cc_ind = xcalloc(1+fhv->m_max, sizeof(int));
         fhv->cc_val = xcalloc(1+fhv->m_max, sizeof(double));
      }
      /* try to factorize the basis matrix */
      switch (luf_factorize(fhv->luf, m, col, info))
      {  case 0:
            break;
         case LUF_ESING:
            ret = FHV_ESING;
            goto done;
         case LUF_ECOND:
            ret = FHV_ECOND;
            goto done;
         default:
            xassert(fhv != fhv);
      }
      /* the basis matrix has been successfully factorized */
      fhv->valid = 1;
      /* H := I */
      fhv->hh_nfs = 0;
      /* P0 := P */
      memcpy(&fhv->p0_row[1], &fhv->luf->pp_row[1], sizeof(int) * m);
      memcpy(&fhv->p0_col[1], &fhv->luf->pp_col[1], sizeof(int) * m);
      /* currently H has no factors */
      fhv->nnz_h = 0;
      ret = 0;
done: /* return to the calling program */
      return ret;
}

/***********************************************************************
*  NAME
*
*  fhv_h_solve - solve system H*x = b or H'*x = b
*
*  SYNOPSIS
*
*  #include "glpfhv.h"
*  void fhv_h_solve(FHV *fhv, int tr, double x[]);
*
*  DESCRIPTION
*
*  The routine fhv_h_solve solves either the system H*x = b (if the
*  flag tr is zero) or the system H'*x = b (if the flag tr is non-zero),
*  where the matrix H is a component of the factorization specified by
*  the parameter fhv, H' is a matrix transposed to H.
*
*  On entry the array x should contain elements of the right-hand side
*  vector b in locations x[1], ..., x[m], where m is the order of the
*  matrix H. On exit this array will contain elements of the solution
*  vector x in the same locations. */

void fhv_h_solve(FHV *fhv, int tr, double x[])
{     int nfs = fhv->hh_nfs;
      int *hh_ind = fhv->hh_ind;
      int *hh_ptr = fhv->hh_ptr;
      int *hh_len = fhv->hh_len;
      int *sv_ind = fhv->luf->sv_ind;
      double *sv_val = fhv->luf->sv_val;
      int i, k, beg, end, ptr;
      double temp;
      if (!fhv->valid)
         xfault("fhv_h_solve: the factorization is not valid\n");
      if (!tr)
      {  /* solve the system H*x = b */
         for (k = 1; k <= nfs; k++)
         {  i = hh_ind[k];
            temp = x[i];
            beg = hh_ptr[k];
            end = beg + hh_len[k] - 1;
            for (ptr = beg; ptr <= end; ptr++)
               temp -= sv_val[ptr] * x[sv_ind[ptr]];
            x[i] = temp;
         }
      }
      else
      {  /* solve the system H'*x = b */
         for (k = nfs; k >= 1; k--)
         {  i = hh_ind[k];
            temp = x[i];
            if (temp == 0.0) continue;
            beg = hh_ptr[k];
            end = beg + hh_len[k] - 1;
            for (ptr = beg; ptr <= end; ptr++)
               x[sv_ind[ptr]] -= sv_val[ptr] * temp;
         }
      }
      return;
}

/***********************************************************************
*  NAME
*
*  fhv_ftran - perform forward transformation (solve system B*x = b)
*
*  SYNOPSIS
*
*  #include "glpfhv.h"
*  void fhv_ftran(FHV *fhv, double x[]);
*
*  DESCRIPTION
*
*  The routine fhv_ftran performs forward transformation, i.e. solves
*  the system B*x = b, where B is the basis matrix, x is the vector of
*  unknowns to be computed, b is the vector of right-hand sides.
*
*  On entry elements of the vector b should be stored in dense format
*  in locations x[1], ..., x[m], where m is the number of rows. On exit
*  the routine stores elements of the vector x in the same locations. */

void fhv_ftran(FHV *fhv, double x[])
{     int *pp_row = fhv->luf->pp_row;
      int *pp_col = fhv->luf->pp_col;
      int *p0_row = fhv->p0_row;
      int *p0_col = fhv->p0_col;
      if (!fhv->valid)
         xfault("fhv_ftran: the factorization is not valid\n");
      /* B = F*H*V, therefore inv(B) = inv(V)*inv(H)*inv(F) */
      fhv->luf->pp_row = p0_row;
      fhv->luf->pp_col = p0_col;
      luf_f_solve(fhv->luf, 0, x);
      fhv->luf->pp_row = pp_row;
      fhv->luf->pp_col = pp_col;
      fhv_h_solve(fhv, 0, x);
      luf_v_solve(fhv->luf, 0, x);
      return;
}

/***********************************************************************
*  NAME
*
*  fhv_btran - perform backward transformation (solve system B'*x = b)
*
*  SYNOPSIS
*
*  #include "glpfhv.h"
*  void fhv_btran(FHV *fhv, double x[]);
*
*  DESCRIPTION
*
*  The routine fhv_btran performs backward transformation, i.e. solves
*  the system B'*x = b, where B' is a matrix transposed to the basis
*  matrix B, x is the vector of unknowns to be computed, b is the vector
*  of right-hand sides.
*
*  On entry elements of the vector b should be stored in dense format
*  in locations x[1], ..., x[m], where m is the number of rows. On exit
*  the routine stores elements of the vector x in the same locations. */

void fhv_btran(FHV *fhv, double x[])
{     int *pp_row = fhv->luf->pp_row;
      int *pp_col = fhv->luf->pp_col;
      int *p0_row = fhv->p0_row;
      int *p0_col = fhv->p0_col;
      if (!fhv->valid)
         xfault("fhv_btran: the factorization is not valid\n");
      /* B = F*H*V, therefore inv(B') = inv(F')*inv(H')*inv(V') */
      luf_v_solve(fhv->luf, 1, x);
      fhv_h_solve(fhv, 1, x);
      fhv->luf->pp_row = p0_row;
      fhv->luf->pp_col = p0_col;
      luf_f_solve(fhv->luf, 1, x);
      fhv->luf->pp_row = pp_row;
      fhv->luf->pp_col = pp_col;
      return;
}

/***********************************************************************
*  NAME
*
*  fhv_update_it - update LP basis factorization
*
*  SYNOPSIS
*
*  #include "glpfhv.h"
*  int fhv_update_it(FHV *fhv, int j, int len, const int ind[],
*     const double val[]);
*
*  DESCRIPTION
*
*  The routine fhv_update_it updates the factorization of the basis
*  matrix B after replacing its j-th column by a new vector.
*
*  The parameter j specifies the number of column of B, which has been
*  replaced, 1 <= j <= m, where m is the order of B.
*
*  Row indices and numerical values of non-zero elements of the new
*  column of B should be placed in locations ind[1], ..., ind[len] and
*  val[1], ..., val[len], resp., where len is the number of non-zeros
*  in the column. Neither zero nor duplicate elements are allowed.
*
*  RETURNS
*
*  0  The factorization has been successfully updated.
*
*  FHV_ESING
*     The adjacent basis matrix is structurally singular, since after
*     changing j-th column of matrix V by the new column (see algorithm
*     below) the case k1 > k2 occured.
*
*  FHV_ECHECK
*     The factorization is inaccurate, since after transforming k2-th
*     row of matrix U = P*V*Q, its diagonal element u[k2,k2] is zero or
*     close to zero,
*
*  FHV_ELIMIT
*     Maximal number of H factors has been reached.
*
*  FHV_EROOM
*     Overflow of the sparse vector area.
*
*  In case of non-zero return code the factorization becomes invalid.
*  It should not be used until it has been recomputed with the routine
*  fhv_factorize.
*
*  ALGORITHM
*
*  The routine fhv_update_it is based on the transformation proposed by
*  Forrest and Tomlin.
*
*  Let j-th column of the basis matrix B have been replaced by new
*  column B[j]. In order to keep the equality B = F*H*V j-th column of
*  matrix V should be replaced by the column inv(F*H)*B[j].
*
*  From the standpoint of matrix U = P*V*Q, replacement of j-th column
*  of matrix V is equivalent to replacement of k1-th column of matrix U,
*  where k1 is determined by permutation matrix Q. Thus, matrix U loses
*  its upper triangular form and becomes the following:
*
*         1   k1       k2   m
*     1   x x * x x x x x x x
*         . x * x x x x x x x
*     k1  . . * x x x x x x x
*         . . * x x x x x x x
*         . . * . x x x x x x
*         . . * . . x x x x x
*         . . * . . . x x x x
*     k2  . . * . . . . x x x
*         . . . . . . . . x x
*     m   . . . . . . . . . x
*
*  where row index k2 corresponds to the lowest non-zero element of
*  k1-th column.
*
*  The routine moves rows and columns k1+1, k1+2, ..., k2 of matrix U
*  by one position to the left and upwards and moves k1-th row and k1-th
*  column to position k2. As the result of such symmetric permutations
*  matrix U becomes the following:
*
*         1   k1       k2   m
*     1   x x x x x x x * x x
*         . x x x x x x * x x
*     k1  . . x x x x x * x x
*         . . . x x x x * x x
*         . . . . x x x * x x
*         . . . . . x x * x x
*         . . . . . . x * x x
*     k2  . . x x x x x * x x
*         . . . . . . . . x x
*     m   . . . . . . . . . x
*
*  Then the routine performs gaussian elimination to eliminate elements
*  u[k2,k1], u[k2,k1+1], ..., u[k2,k2-1] using diagonal elements
*  u[k1,k1], u[k1+1,k1+1], ..., u[k2-1,k2-1] as pivots in the same way
*  as described in comments to the routine luf_factorize (see the module
*  GLPLUF). Note that actually all operations are performed on matrix V,
*  not on matrix U. During the elimination process the routine permutes
*  neither rows nor columns, so only k2-th row of matrix U is changed.
*
*  To keep the main equality B = F*H*V, each time when the routine
*  applies elementary gaussian transformation to the transformed row of
*  matrix V (which corresponds to k2-th row of matrix U), it also adds
*  a new element (gaussian multiplier) to the current row-like factor
*  of matrix H, which corresponds to the transformed row of matrix V. */

int fhv_update_it(FHV *fhv, int j, int len, const int ind[],
      const double val[])
{     int m = fhv->m;
      LUF *luf = fhv->luf;
      int *vr_ptr = luf->vr_ptr;
      int *vr_len = luf->vr_len;
      int *vr_cap = luf->vr_cap;
      double *vr_piv = luf->vr_piv;
      int *vc_ptr = luf->vc_ptr;
      int *vc_len = luf->vc_len;
      int *vc_cap = luf->vc_cap;
      int *pp_row = luf->pp_row;
      int *pp_col = luf->pp_col;
      int *qq_row = luf->qq_row;
      int *qq_col = luf->qq_col;
      int *sv_ind = luf->sv_ind;
      double *sv_val = luf->sv_val;
      double *work = luf->work;
      double eps_tol = luf->eps_tol;
      int *hh_ind = fhv->hh_ind;
      int *hh_ptr = fhv->hh_ptr;
      int *hh_len = fhv->hh_len;
      int *p0_row = fhv->p0_row;
      int *p0_col = fhv->p0_col;
      int *cc_ind = fhv->cc_ind;
      double *cc_val = fhv->cc_val;
      double upd_tol = fhv->upd_tol;
      int i, i_beg, i_end, i_ptr, j_beg, j_end, j_ptr, k, k1, k2, p, q,
         p_beg, p_end, p_ptr, ptr, ret;
      double f, temp;
      if (!fhv->valid)
         xfault("fhv_update_it: the factorization is not valid\n");
      if (!(1 <= j && j <= m))
         xfault("fhv_update_it: j = %d; column number out of range\n",
            j);
      /* check if the new factor of matrix H can be created */
      if (fhv->hh_nfs == fhv->hh_max)
      {  /* maximal number of updates has been reached */
         fhv->valid = 0;
         ret = FHV_ELIMIT;
         goto done;
      }
      /* convert new j-th column of B to dense format */
      for (i = 1; i <= m; i++)
         cc_val[i] = 0.0;
      for (k = 1; k <= len; k++)
      {  i = ind[k];
         if (!(1 <= i && i <= m))
            xfault("fhv_update_it: ind[%d] = %d; row number out of rang"
               "e\n", k, i);
         if (cc_val[i] != 0.0)
            xfault("fhv_update_it: ind[%d] = %d; duplicate row index no"
               "t allowed\n", k, i);
         if (val[k] == 0.0)
            xfault("fhv_update_it: val[%d] = %g; zero element not allow"
               "ed\n", k, val[k]);
         cc_val[i] = val[k];
      }
      /* new j-th column of V := inv(F * H) * (new B[j]) */
      fhv->luf->pp_row = p0_row;
      fhv->luf->pp_col = p0_col;
      luf_f_solve(fhv->luf, 0, cc_val);
      fhv->luf->pp_row = pp_row;
      fhv->luf->pp_col = pp_col;
      fhv_h_solve(fhv, 0, cc_val);
      /* convert new j-th column of V to sparse format */
      len = 0;
      for (i = 1; i <= m; i++)
      {  temp = cc_val[i];
         if (temp == 0.0 || fabs(temp) < eps_tol) continue;
         len++, cc_ind[len] = i, cc_val[len] = temp;
      }
      /* clear old content of j-th column of matrix V */
      j_beg = vc_ptr[j];
      j_end = j_beg + vc_len[j] - 1;
      for (j_ptr = j_beg; j_ptr <= j_end; j_ptr++)
      {  /* get row index of v[i,j] */
         i = sv_ind[j_ptr];
         /* find v[i,j] in the i-th row */
         i_beg = vr_ptr[i];
         i_end = i_beg + vr_len[i] - 1;
         for (i_ptr = i_beg; sv_ind[i_ptr] != j; i_ptr++) /* nop */;
         xassert(i_ptr <= i_end);
         /* remove v[i,j] from the i-th row */
         sv_ind[i_ptr] = sv_ind[i_end];
         sv_val[i_ptr] = sv_val[i_end];
         vr_len[i]--;
      }
      /* now j-th column of matrix V is empty */
      luf->nnz_v -= vc_len[j];
      vc_len[j] = 0;
      /* add new elements of j-th column of matrix V to corresponding
         row lists; determine indices k1 and k2 */
      k1 = qq_row[j], k2 = 0;
      for (ptr = 1; ptr <= len; ptr++)
      {  /* get row index of v[i,j] */
         i = cc_ind[ptr];
         /* at least one unused location is needed in i-th row */
         if (vr_len[i] + 1 > vr_cap[i])
         {  if (luf_enlarge_row(luf, i, vr_len[i] + 10))
            {  /* overflow of the sparse vector area */
               fhv->valid = 0;
               luf->new_sva = luf->sv_size + luf->sv_size;
               xassert(luf->new_sva > luf->sv_size);
               ret = FHV_EROOM;
               goto done;
            }
         }
         /* add v[i,j] to i-th row */
         i_ptr = vr_ptr[i] + vr_len[i];
         sv_ind[i_ptr] = j;
         sv_val[i_ptr] = cc_val[ptr];
         vr_len[i]++;
         /* adjust index k2 */
         if (k2 < pp_col[i]) k2 = pp_col[i];
      }
      /* capacity of j-th column (which is currently empty) should be
         not less than len locations */
      if (vc_cap[j] < len)
      {  if (luf_enlarge_col(luf, j, len))
         {  /* overflow of the sparse vector area */
            fhv->valid = 0;
            luf->new_sva = luf->sv_size + luf->sv_size;
            xassert(luf->new_sva > luf->sv_size);
            ret = FHV_EROOM;
            goto done;
         }
      }
      /* add new elements of matrix V to j-th column list */
      j_ptr = vc_ptr[j];
      memmove(&sv_ind[j_ptr], &cc_ind[1], len * sizeof(int));
      memmove(&sv_val[j_ptr], &cc_val[1], len * sizeof(double));
      vc_len[j] = len;
      luf->nnz_v += len;
      /* if k1 > k2, diagonal element u[k2,k2] of matrix U is zero and
         therefore the adjacent basis matrix is structurally singular */
      if (k1 > k2)
      {  fhv->valid = 0;
         ret = FHV_ESING;
         goto done;
      }
      /* perform implicit symmetric permutations of rows and columns of
         matrix U */
      i = pp_row[k1], j = qq_col[k1];
      for (k = k1; k < k2; k++)
      {  pp_row[k] = pp_row[k+1], pp_col[pp_row[k]] = k;
         qq_col[k] = qq_col[k+1], qq_row[qq_col[k]] = k;
      }
      pp_row[k2] = i, pp_col[i] = k2;
      qq_col[k2] = j, qq_row[j] = k2;
      /* now i-th row of the matrix V is k2-th row of matrix U; since
         no pivoting is used, only this row will be transformed */
      /* copy elements of i-th row of matrix V to the working array and
         remove these elements from matrix V */
      for (j = 1; j <= m; j++) work[j] = 0.0;
      i_beg = vr_ptr[i];
      i_end = i_beg + vr_len[i] - 1;
      for (i_ptr = i_beg; i_ptr <= i_end; i_ptr++)
      {  /* get column index of v[i,j] */
         j = sv_ind[i_ptr];
         /* store v[i,j] to the working array */
         work[j] = sv_val[i_ptr];
         /* find v[i,j] in the j-th column */
         j_beg = vc_ptr[j];
         j_end = j_beg + vc_len[j] - 1;
         for (j_ptr = j_beg; sv_ind[j_ptr] != i; j_ptr++) /* nop */;
         xassert(j_ptr <= j_end);
         /* remove v[i,j] from the j-th column */
         sv_ind[j_ptr] = sv_ind[j_end];
         sv_val[j_ptr] = sv_val[j_end];
         vc_len[j]--;
      }
      /* now i-th row of matrix V is empty */
      luf->nnz_v -= vr_len[i];
      vr_len[i] = 0;
      /* create the next row-like factor of the matrix H; this factor
         corresponds to i-th (transformed) row */
      fhv->hh_nfs++;
      hh_ind[fhv->hh_nfs] = i;
      /* hh_ptr[] will be set later */
      hh_len[fhv->hh_nfs] = 0;
      /* up to (k2 - k1) free locations are needed to add new elements
         to the non-trivial row of the row-like factor */
      if (luf->sv_end - luf->sv_beg < k2 - k1)
      {  luf_defrag_sva(luf);
         if (luf->sv_end - luf->sv_beg < k2 - k1)
         {  /* overflow of the sparse vector area */
            fhv->valid = luf->valid = 0;
            luf->new_sva = luf->sv_size + luf->sv_size;
            xassert(luf->new_sva > luf->sv_size);
            ret = FHV_EROOM;
            goto done;
         }
      }
      /* eliminate subdiagonal elements of matrix U */
      for (k = k1; k < k2; k++)
      {  /* v[p,q] = u[k,k] */
         p = pp_row[k], q = qq_col[k];
         /* this is the crucial point, where even tiny non-zeros should
            not be dropped */
         if (work[q] == 0.0) continue;
         /* compute gaussian multiplier f = v[i,q] / v[p,q] */
         f = work[q] / vr_piv[p];
         /* perform gaussian transformation:
            (i-th row) := (i-th row) - f * (p-th row)
            in order to eliminate v[i,q] = u[k2,k] */
         p_beg = vr_ptr[p];
         p_end = p_beg + vr_len[p] - 1;
         for (p_ptr = p_beg; p_ptr <= p_end; p_ptr++)
            work[sv_ind[p_ptr]] -= f * sv_val[p_ptr];
         /* store new element (gaussian multiplier that corresponds to
            p-th row) in the current row-like factor */
         luf->sv_end--;
         sv_ind[luf->sv_end] = p;
         sv_val[luf->sv_end] = f;
         hh_len[fhv->hh_nfs]++;
      }
      /* set pointer to the current row-like factor of the matrix H
         (if no elements were added to this factor, it is unity matrix
         and therefore can be discarded) */
      if (hh_len[fhv->hh_nfs] == 0)
         fhv->hh_nfs--;
      else
      {  hh_ptr[fhv->hh_nfs] = luf->sv_end;
         fhv->nnz_h += hh_len[fhv->hh_nfs];
      }
      /* store new pivot which corresponds to u[k2,k2] */
      vr_piv[i] = work[qq_col[k2]];
      /* new elements of i-th row of matrix V (which are non-diagonal
         elements u[k2,k2+1], ..., u[k2,m] of matrix U = P*V*Q) now are
         contained in the working array; add them to matrix V */
      len = 0;
      for (k = k2+1; k <= m; k++)
      {  /* get column index and value of v[i,j] = u[k2,k] */
         j = qq_col[k];
         temp = work[j];
         /* if v[i,j] is close to zero, skip it */
         if (fabs(temp) < eps_tol) continue;
         /* at least one unused location is needed in j-th column */
         if (vc_len[j] + 1 > vc_cap[j])
         {  if (luf_enlarge_col(luf, j, vc_len[j] + 10))
            {  /* overflow of the sparse vector area */
               fhv->valid = 0;
               luf->new_sva = luf->sv_size + luf->sv_size;
               xassert(luf->new_sva > luf->sv_size);
               ret = FHV_EROOM;
               goto done;
            }
         }
         /* add v[i,j] to j-th column */
         j_ptr = vc_ptr[j] + vc_len[j];
         sv_ind[j_ptr] = i;
         sv_val[j_ptr] = temp;
         vc_len[j]++;
         /* also store v[i,j] to the auxiliary array */
         len++, cc_ind[len] = j, cc_val[len] = temp;
      }
      /* capacity of i-th row (which is currently empty) should be not
         less than len locations */
      if (vr_cap[i] < len)
      {  if (luf_enlarge_row(luf, i, len))
         {  /* overflow of the sparse vector area */
            fhv->valid = 0;
            luf->new_sva = luf->sv_size + luf->sv_size;
            xassert(luf->new_sva > luf->sv_size);
            ret = FHV_EROOM;
            goto done;
         }
      }
      /* add new elements to i-th row list */
      i_ptr = vr_ptr[i];
      memmove(&sv_ind[i_ptr], &cc_ind[1], len * sizeof(int));
      memmove(&sv_val[i_ptr], &cc_val[1], len * sizeof(double));
      vr_len[i] = len;
      luf->nnz_v += len;
      /* updating is finished; check that diagonal element u[k2,k2] is
         not very small in absolute value among other elements in k2-th
         row and k2-th column of matrix U = P*V*Q */
      /* temp = max(|u[k2,*]|, |u[*,k2]|) */
      temp = 0.0;
      /* walk through k2-th row of U which is i-th row of V */
      i = pp_row[k2];
      i_beg = vr_ptr[i];
      i_end = i_beg + vr_len[i] - 1;
      for (i_ptr = i_beg; i_ptr <= i_end; i_ptr++)
         if (temp < fabs(sv_val[i_ptr])) temp = fabs(sv_val[i_ptr]);
      /* walk through k2-th column of U which is j-th column of V */
      j = qq_col[k2];
      j_beg = vc_ptr[j];
      j_end = j_beg + vc_len[j] - 1;
      for (j_ptr = j_beg; j_ptr <= j_end; j_ptr++)
         if (temp < fabs(sv_val[j_ptr])) temp = fabs(sv_val[j_ptr]);
      /* check that u[k2,k2] is not very small */
      if (fabs(vr_piv[i]) < upd_tol * temp)
      {  /* the factorization seems to be inaccurate and therefore must
            be recomputed */
         fhv->valid = 0;
         ret = FHV_ECHECK;
         goto done;
      }
      /* the factorization has been successfully updated */
      ret = 0;
done: /* return to the calling program */
      return ret;
}

/***********************************************************************
*  NAME
*
*  fhv_delete_it - delete LP basis factorization
*
*  SYNOPSIS
*
*  #include "glpfhv.h"
*  void fhv_delete_it(FHV *fhv);
*
*  DESCRIPTION
*
*  The routine fhv_delete_it deletes LP basis factorization specified
*  by the parameter fhv and frees all memory allocated to this program
*  object. */

void fhv_delete_it(FHV *fhv)
{     luf_delete_it(fhv->luf);
      if (fhv->hh_ind != NULL) xfree(fhv->hh_ind);
      if (fhv->hh_ptr != NULL) xfree(fhv->hh_ptr);
      if (fhv->hh_len != NULL) xfree(fhv->hh_len);
      if (fhv->p0_row != NULL) xfree(fhv->p0_row);
      if (fhv->p0_col != NULL) xfree(fhv->p0_col);
      if (fhv->cc_ind != NULL) xfree(fhv->cc_ind);
      if (fhv->cc_val != NULL) xfree(fhv->cc_val);
      xfree(fhv);
      return;
}

/* eof */
