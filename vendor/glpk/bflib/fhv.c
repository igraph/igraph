/* fhv.c (sparse updatable FHV-factorization) */

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

#include "env.h"
#include "fhv.h"

/***********************************************************************
*  fhv_ft_update - update FHV-factorization (Forrest-Tomlin)
*
*  This routine updates FHV-factorization of the original matrix A
*  after replacing its j-th column by a new one. The routine is based
*  on the method proposed by Forrest and Tomlin [1].
*
*  The parameter q specifies the number of column of A, which has been
*  replaced, 1 <= q <= n, where n is the order of A.
*
*  Row indices and numerical values of non-zero elements of the new
*  j-th column of A should be placed in locations aq_ind[1], ...,
*  aq_ind[aq_len] and aq_val[1], ..., aq_val[aq_len], respectively,
*  where aq_len is the number of non-zeros. Neither zero nor duplicate
*  elements are allowed.
*
*  The working arrays ind, val, and work should have at least 1+n
*  elements (0-th elements are not used).
*
*  RETURNS
*
*  0  The factorization has been successfully updated.
*
*  1  New matrix U = P'* V * Q' is upper triangular with zero diagonal
*     element u[s,s]. (Elimination was not performed.)
*
*  2  New matrix U = P'* V * Q' is upper triangular, and its diagonal
*     element u[s,s] or u[t,t] is too small in magnitude. (Elimination
*     was not performed.)
*
*  3  The same as 2, but after performing elimination.
*
*  4  The factorization has not been updated, because maximal number of
*     updates has been reached.
*
*  5  Accuracy test failed for the updated factorization.
*
*  BACKGROUND
*
*  The routine is based on the updating method proposed by Forrest and
*  Tomlin [1].
*
*  Let q-th column of the original matrix A have been replaced by new
*  column A[q]. Then, to keep the equality A = F * H * V, q-th column
*  of matrix V should be replaced by column V[q] = inv(F * H) * A[q].
*  From the standpoint of matrix U = P'* V * Q' such replacement is
*  equivalent to replacement of s-th column of matrix U, where s is
*  determined from q by permutation matrix Q. Thus, matrix U loses its
*  upper triangular form and becomes the following:
*
*        1   s       t   n
*     1  x x * x x x x x x
*        . x * x x x x x x
*     s  . . * x x x x x x
*        . . * x x x x x x
*        . . * . x x x x x
*        . . * . . x x x x
*     t  . . * . . . x x x
*        . . . . . . . x x
*     n  . . . . . . . . x
*
*  where t is largest row index of a non-zero element in s-th column.
*
*  The routine makes matrix U upper triangular as follows. First, it
*  moves rows and columns s+1, ..., t by one position to the left and
*  upwards, resp., and moves s-th row and s-th column to position t.
*  Due to such symmetric permutations matrix U becomes the following
*  (note that all diagonal elements remain on the diagonal, and element
*  u[s,s] becomes u[t,t]):
*
*        1   s       t   n
*     1  x x x x x x * x x
*        . x x x x x * x x
*     s  . . x x x x * x x
*        . . . x x x * x x
*        . . . . x x * x x
*        . . . . . x * x x
*     t  . . x x x x * x x
*        . . . . . . . x x
*     n  . . . . . . . . x
*
*  Then the routine performs gaussian elimination to eliminate
*  subdiagonal elements u[t,s], ..., u[t,t-1] using diagonal elements
*  u[s,s], ..., u[t-1,t-1] as pivots. During the elimination process
*  the routine permutes neither rows nor columns, so only t-th row is
*  changed. Should note that actually all operations are performed on
*  matrix V = P * U * Q, since matrix U is not stored.
*
*  To keep the equality A = F * H * V, the routine appends new row-like
*  factor H[k] to matrix H, and every time it applies elementary
*  gaussian transformation to eliminate u[t,j'] = v[p,j] using pivot
*  u[j',j'] = v[i,j], it also adds new element f[p,j] = v[p,j] / v[i,j]
*  (gaussian multiplier) to factor H[k], which initially is a unity
*  matrix. At the end of elimination process the row-like factor H[k]
*  may look as follows:
*
*        1               n          1   s       t   n
*     1  1 . . . . . . . .       1  1 . . . . . . . .
*        . 1 . . . . . . .          . 1 . . . . . . .
*        . . 1 . . . . . .       s  . . 1 . . . . . .
*     p  . x x 1 . x . x .          . . . 1 . . . . .
*        . . . . 1 . . . .          . . . . 1 . . . .
*        . . . . . 1 . . .          . . . . . 1 . . .
*        . . . . . . 1 . .       t  . . x x x x 1 . .
*        . . . . . . . 1 .          . . . . . . . 1 .
*     n  . . . . . . . . 1       n  . . . . . . . . 1
*
*              H[k]                 inv(P) * H[k] * P
*
*  If, however, s = t, no elimination is needed, in which case no new
*  row-like factor is created.
*
*  REFERENCES
*
*  1. J.J.H.Forrest and J.A.Tomlin, "Updated triangular factors of the
*     basis to maintain sparsity in the product form simplex method,"
*     Math. Prog. 2 (1972), pp. 263-78. */

int fhv_ft_update(FHV *fhv, int q, int aq_len, const int aq_ind[],
      const double aq_val[], int ind[/*1+n*/], double val[/*1+n*/],
      double work[/*1+n*/])
{     LUF *luf = fhv->luf;
      int n = luf->n;
      SVA *sva = luf->sva;
      int *sv_ind = sva->ind;
      double *sv_val = sva->val;
      int vr_ref = luf->vr_ref;
      int *vr_ptr = &sva->ptr[vr_ref-1];
      int *vr_len = &sva->len[vr_ref-1];
      int *vr_cap = &sva->cap[vr_ref-1];
      double *vr_piv = luf->vr_piv;
      int vc_ref = luf->vc_ref;
      int *vc_ptr = &sva->ptr[vc_ref-1];
      int *vc_len = &sva->len[vc_ref-1];
      int *vc_cap = &sva->cap[vc_ref-1];
      int *pp_ind = luf->pp_ind;
      int *pp_inv = luf->pp_inv;
      int *qq_ind = luf->qq_ind;
      int *qq_inv = luf->qq_inv;
      int *hh_ind = fhv->hh_ind;
      int hh_ref = fhv->hh_ref;
      int *hh_ptr = &sva->ptr[hh_ref-1];
      int *hh_len = &sva->len[hh_ref-1];
#if 1 /* FIXME */
      const double eps_tol = DBL_EPSILON;
      const double vpq_tol = 1e-5;
      const double err_tol = 1e-10;
#endif
      int end, i, i_end, i_ptr, j, j_end, j_ptr, k, len, nnz, p, p_end,
         p_ptr, ptr, q_end, q_ptr, s, t;
      double f, vpq, temp;
      /*--------------------------------------------------------------*/
      /* replace current q-th column of matrix V by new one           */
      /*--------------------------------------------------------------*/
      xassert(1 <= q && q <= n);
      /* convert new q-th column of matrix A to dense format */
      for (i = 1; i <= n; i++)
         val[i] = 0.0;
      xassert(0 <= aq_len && aq_len <= n);
      for (k = 1; k <= aq_len; k++)
      {  i = aq_ind[k];
         xassert(1 <= i && i <= n);
         xassert(val[i] == 0.0);
         xassert(aq_val[k] != 0.0);
         val[i] = aq_val[k];
      }
      /* compute new q-th column of matrix V:
       * new V[q] = inv(F * H) * (new A[q]) */
      luf->pp_ind = fhv->p0_ind;
      luf->pp_inv = fhv->p0_inv;
      luf_f_solve(luf, val);
      luf->pp_ind = pp_ind;
      luf->pp_inv = pp_inv;
      fhv_h_solve(fhv, val);
      /* q-th column of V = s-th column of U */
      s = qq_inv[q];
      /* determine row number of element v[p,q] that corresponds to
       * diagonal element u[s,s] */
      p = pp_inv[s];
      /* convert new q-th column of V to sparse format;
       * element v[p,q] = u[s,s] is not included in the element list
       * and stored separately */
      vpq = 0.0;
      len = 0;
      for (i = 1; i <= n; i++)
      {  temp = val[i];
#if 1 /* FIXME */
         if (-eps_tol < temp && temp < +eps_tol)
#endif
            /* nop */;
         else if (i == p)
            vpq = temp;
         else
         {  ind[++len] = i;
            val[len] = temp;
         }
      }
      /* clear q-th column of matrix V */
      for (q_end = (q_ptr = vc_ptr[q]) + vc_len[q];
         q_ptr < q_end; q_ptr++)
      {  /* get row index of v[i,q] */
         i = sv_ind[q_ptr];
         /* find and remove v[i,q] from i-th row */
         for (i_end = (i_ptr = vr_ptr[i]) + vr_len[i];
            sv_ind[i_ptr] != q; i_ptr++)
            /* nop */;
         xassert(i_ptr < i_end);
         sv_ind[i_ptr] = sv_ind[i_end-1];
         sv_val[i_ptr] = sv_val[i_end-1];
         vr_len[i]--;
      }
      /* now q-th column of matrix V is empty */
      vc_len[q] = 0;
      /* put new q-th column of V (except element v[p,q] = u[s,s]) in
       * column-wise format */
      if (len > 0)
      {  if (vc_cap[q] < len)
         {  if (sva->r_ptr - sva->m_ptr < len)
            {  sva_more_space(sva, len);
               sv_ind = sva->ind;
               sv_val = sva->val;
            }
            sva_enlarge_cap(sva, vc_ref-1+q, len, 0);
         }
         ptr = vc_ptr[q];
         memcpy(&sv_ind[ptr], &ind[1], len * sizeof(int));
         memcpy(&sv_val[ptr], &val[1], len * sizeof(double));
         vc_len[q] = len;
      }
      /* put new q-th column of V (except element v[p,q] = u[s,s]) in
       * row-wise format, and determine largest row number t such that
       * u[s,t] != 0 */
      t = (vpq == 0.0 ? 0 : s);
      for (k = 1; k <= len; k++)
      {  /* get row index of v[i,q] */
         i = ind[k];
         /* put v[i,q] to i-th row */
         if (vr_cap[i] == vr_len[i])
         {  /* reserve extra locations in i-th row to reduce further
             * relocations of that row */
#if 1 /* FIXME */
            int need = vr_len[i] + 5;
#endif
            if (sva->r_ptr - sva->m_ptr < need)
            {  sva_more_space(sva, need);
               sv_ind = sva->ind;
               sv_val = sva->val;
            }
            sva_enlarge_cap(sva, vr_ref-1+i, need, 0);
         }
         sv_ind[ptr = vr_ptr[i] + (vr_len[i]++)] = q;
         sv_val[ptr] = val[k];
         /* v[i,q] is non-zero; increase t */
         if (t < pp_ind[i])
            t = pp_ind[i];
      }
      /*--------------------------------------------------------------*/
      /* check if matrix U is already upper triangular                */
      /*--------------------------------------------------------------*/
      /* check if there is a spike in s-th column of matrix U, which
       * is q-th column of matrix V */
      if (s >= t)
      {  /* no spike; matrix U is already upper triangular */
         /* store its diagonal element u[s,s] = v[p,q] */
         vr_piv[p] = vpq;
         if (s > t)
         {  /* matrix U is structurally singular, because its diagonal
             * element u[s,s] = v[p,q] is exact zero */
            xassert(vpq == 0.0);
            return 1;
         }
#if 1 /* FIXME */
         else if (-vpq_tol < vpq && vpq < +vpq_tol)
#endif
         {  /* matrix U is not well conditioned, because its diagonal
             * element u[s,s] = v[p,q] is too small in magnitude */
            return 2;
         }
         else
         {  /* normal case */
            return 0;
         }
      }
      /*--------------------------------------------------------------*/
      /* perform implicit symmetric permutations of rows and columns  */
      /* of matrix U                                                  */
      /*--------------------------------------------------------------*/
      /* currently v[p,q] = u[s,s] */
      xassert(p == pp_inv[s] && q == qq_ind[s]);
      for (k = s; k < t; k++)
      {  pp_ind[pp_inv[k] = pp_inv[k+1]] = k;
         qq_inv[qq_ind[k] = qq_ind[k+1]] = k;
      }
      /* now v[p,q] = u[t,t] */
      pp_ind[pp_inv[t] = p] = qq_inv[qq_ind[t] = q] = t;
      /*--------------------------------------------------------------*/
      /* check if matrix U is already upper triangular                */
      /*--------------------------------------------------------------*/
      /* check if there is a spike in t-th row of matrix U, which is
       * p-th row of matrix V */
      for (p_end = (p_ptr = vr_ptr[p]) + vr_len[p];
         p_ptr < p_end; p_ptr++)
      {  if (qq_inv[sv_ind[p_ptr]] < t)
            break; /* spike detected */
      }
      if (p_ptr == p_end)
      {  /* no spike; matrix U is already upper triangular */
         /* store its diagonal element u[t,t] = v[p,q] */
         vr_piv[p] = vpq;
#if 1 /* FIXME */
         if (-vpq_tol < vpq && vpq < +vpq_tol)
#endif
         {  /* matrix U is not well conditioned, because its diagonal
             * element u[t,t] = v[p,q] is too small in magnitude */
            return 2;
         }
         else
         {  /* normal case */
            return 0;
         }
      }
      /*--------------------------------------------------------------*/
      /* copy p-th row of matrix V, which is t-th row of matrix U, to */
      /* working array                                                */
      /*--------------------------------------------------------------*/
      /* copy p-th row of matrix V, including element v[p,q] = u[t,t],
       * to the working array in dense format and remove these elements
       * from matrix V; since no pivoting is used, only this row will
       * change during elimination */
      for (j = 1; j <= n; j++)
         work[j] = 0.0;
      work[q] = vpq;
      for (p_end = (p_ptr = vr_ptr[p]) + vr_len[p];
         p_ptr < p_end; p_ptr++)
      {  /* get column index of v[p,j] and store this element to the
          * working array */
         work[j = sv_ind[p_ptr]] = sv_val[p_ptr];
         /* find and remove v[p,j] from j-th column */
         for (j_end = (j_ptr = vc_ptr[j]) + vc_len[j];
            sv_ind[j_ptr] != p; j_ptr++)
            /* nop */;
         xassert(j_ptr < j_end);
         sv_ind[j_ptr] = sv_ind[j_end-1];
         sv_val[j_ptr] = sv_val[j_end-1];
         vc_len[j]--;
      }
      /* now p-th row of matrix V is temporarily empty */
      vr_len[p] = 0;
      /*--------------------------------------------------------------*/
      /* perform gaussian elimination                                 */
      /*--------------------------------------------------------------*/
      /* transform p-th row of matrix V stored in working array, which
       * is t-th row of matrix U, to eliminate subdiagonal elements
       * u[t,s], ..., u[t,t-1]; corresponding gaussian multipliers will
       * form non-trivial row of new row-like factor */
      nnz = 0; /* number of non-zero gaussian multipliers */
      for (k = s; k < t; k++)
      {  /* diagonal element u[k,k] = v[i,j] is used as pivot */
         i = pp_inv[k], j = qq_ind[k];
         /* take subdiagonal element u[t,k] = v[p,j] */
         temp = work[j];
#if 1 /* FIXME */
         if (-eps_tol < temp && temp < +eps_tol)
            continue;
#endif
         /* compute and save gaussian multiplier:
          * f := u[t,k] / u[k,k] = v[p,j] / v[i,j] */
         ind[++nnz] = i;
         val[nnz] = f = work[j] / vr_piv[i];
         /* gaussian transformation to eliminate u[t,k] = v[p,j]:
          * (p-th row of V) := (p-th row of V) - f * (i-th row of V) */
         for (i_end = (i_ptr = vr_ptr[i]) + vr_len[i];
            i_ptr < i_end; i_ptr++)
            work[sv_ind[i_ptr]] -= f * sv_val[i_ptr];
      }
      /* now matrix U is again upper triangular */
#if 1 /* FIXME */
      if (-vpq_tol < work[q] && work[q] < +vpq_tol)
#endif
      {  /* however, its new diagonal element u[t,t] = v[p,q] is too
          * small in magnitude */
         return 3;
      }
      /*--------------------------------------------------------------*/
      /* create new row-like factor H[k] and add to eta file H        */
      /*--------------------------------------------------------------*/
      /* (nnz = 0 means that all subdiagonal elements were too small
       * in magnitude) */
      if (nnz > 0)
      {  if (fhv->nfs == fhv->nfs_max)
         {  /* maximal number of row-like factors has been reached */
            return 4;
         }
         k = ++(fhv->nfs);
         hh_ind[k] = p;
         /* store non-trivial row of H[k] in right (dynamic) part of
          * SVA (diagonal unity element is not stored) */
         if (sva->r_ptr - sva->m_ptr < nnz)
         {  sva_more_space(sva, nnz);
            sv_ind = sva->ind;
            sv_val = sva->val;
         }
         sva_reserve_cap(sva, fhv->hh_ref-1+k, nnz);
         ptr = hh_ptr[k];
         memcpy(&sv_ind[ptr], &ind[1], nnz * sizeof(int));
         memcpy(&sv_val[ptr], &val[1], nnz * sizeof(double));
         hh_len[k] = nnz;
      }
      /*--------------------------------------------------------------*/
      /* copy transformed p-th row of matrix V, which is t-th row of  */
      /* matrix U, from working array back to matrix V                */
      /*--------------------------------------------------------------*/
      /* copy elements of transformed p-th row of matrix V, which are
       * non-diagonal elements u[t,t+1], ..., u[t,n] of matrix U, from
       * working array to corresponding columns of matrix V (note that
       * diagonal element u[t,t] = v[p,q] not copied); also transform
       * p-th row of matrix V to sparse format */
      len = 0;
      for (k = t+1; k <= n; k++)
      {  /* j-th column of V = k-th column of U */
         j = qq_ind[k];
         /* take non-diagonal element v[p,j] = u[t,k] */
         temp = work[j];
#if 1 /* FIXME */
         if (-eps_tol < temp && temp < +eps_tol)
            continue;
#endif
         /* add v[p,j] to j-th column of matrix V */
         if (vc_cap[j] == vc_len[j])
         {  /* reserve extra locations in j-th column to reduce further
             * relocations of that column */
#if 1 /* FIXME */
            int need = vc_len[j] + 5;
#endif
            if (sva->r_ptr - sva->m_ptr < need)
            {  sva_more_space(sva, need);
               sv_ind = sva->ind;
               sv_val = sva->val;
            }
            sva_enlarge_cap(sva, vc_ref-1+j, need, 0);
         }
         sv_ind[ptr = vc_ptr[j] + (vc_len[j]++)] = p;
         sv_val[ptr] = temp;
         /* store element v[p,j] = u[t,k] to working sparse vector */
         ind[++len] = j;
         val[len] = temp;
      }
      /* copy elements from working sparse vector to p-th row of matrix
       * V (this row is currently empty) */
      if (vr_cap[p] < len)
      {  if (sva->r_ptr - sva->m_ptr < len)
         {  sva_more_space(sva, len);
            sv_ind = sva->ind;
            sv_val = sva->val;
         }
         sva_enlarge_cap(sva, vr_ref-1+p, len, 0);
      }
      ptr = vr_ptr[p];
      memcpy(&sv_ind[ptr], &ind[1], len * sizeof(int));
      memcpy(&sv_val[ptr], &val[1], len * sizeof(double));
      vr_len[p] = len;
      /* store new diagonal element u[t,t] = v[p,q] */
      vr_piv[p] = work[q];
      /*--------------------------------------------------------------*/
      /* perform accuracy test (only if new H[k] was added)           */
      /*--------------------------------------------------------------*/
      if (nnz > 0)
      {  /* copy p-th (non-trivial) row of row-like factor H[k] (except
          * unity diagonal element) to working array in dense format */
         for (j = 1; j <= n; j++)
            work[j] = 0.0;
         k = fhv->nfs;
         for (end = (ptr = hh_ptr[k]) + hh_len[k]; ptr < end; ptr++)
            work[sv_ind[ptr]] = sv_val[ptr];
         /* compute inner product of p-th (non-trivial) row of matrix
          * H[k] and q-th column of matrix V */
         temp = vr_piv[p]; /* 1 * v[p,q] */
         ptr = vc_ptr[q];
         end = ptr + vc_len[q];
         for (; ptr < end; ptr++)
            temp += work[sv_ind[ptr]] * sv_val[ptr];
         /* inner product should be equal to element v[p,q] *before*
          * matrix V was transformed */
         /* compute relative error */
         temp = fabs(vpq - temp) / (1.0 + fabs(vpq));
#if 1 /* FIXME */
         if (temp > err_tol)
#endif
         {  /* relative error is too large */
            return 5;
         }
      }
      /* factorization has been successfully updated */
      return 0;
}

/***********************************************************************
*  fhv_h_solve - solve system H * x = b
*
*  This routine solves the system H * x = b, where the matrix H is the
*  middle factor of the sparse updatable FHV-factorization.
*
*  On entry the array x should contain elements of the right-hand side
*  vector b in locations x[1], ..., x[n], where n is the order of the
*  matrix H. On exit this array will contain elements of the solution
*  vector x in the same locations. */

void fhv_h_solve(FHV *fhv, double x[/*1+n*/])
{     SVA *sva = fhv->luf->sva;
      int *sv_ind = sva->ind;
      double *sv_val = sva->val;
      int nfs = fhv->nfs;
      int *hh_ind = fhv->hh_ind;
      int hh_ref = fhv->hh_ref;
      int *hh_ptr = &sva->ptr[hh_ref-1];
      int *hh_len = &sva->len[hh_ref-1];
      int i, k, end, ptr;
      double x_i;
      for (k = 1; k <= nfs; k++)
      {  x_i = x[i = hh_ind[k]];
         for (end = (ptr = hh_ptr[k]) + hh_len[k]; ptr < end; ptr++)
            x_i -= sv_val[ptr] * x[sv_ind[ptr]];
         x[i] = x_i;
      }
      return;
}

/***********************************************************************
*  fhv_ht_solve - solve system H' * x = b
*
*  This routine solves the system H' * x = b, where H' is a matrix
*  transposed to the matrix H, which is the middle factor of the sparse
*  updatable FHV-factorization.
*
*  On entry the array x should contain elements of the right-hand side
*  vector b in locations x[1], ..., x[n], where n is the order of the
*  matrix H. On exit this array will contain elements of the solution
*  vector x in the same locations. */

void fhv_ht_solve(FHV *fhv, double x[/*1+n*/])
{     SVA *sva = fhv->luf->sva;
      int *sv_ind = sva->ind;
      double *sv_val = sva->val;
      int nfs = fhv->nfs;
      int *hh_ind = fhv->hh_ind;
      int hh_ref = fhv->hh_ref;
      int *hh_ptr = &sva->ptr[hh_ref-1];
      int *hh_len = &sva->len[hh_ref-1];
      int k, end, ptr;
      double x_j;
      for (k = nfs; k >= 1; k--)
      {  if ((x_j = x[hh_ind[k]]) == 0.0)
            continue;
         for (end = (ptr = hh_ptr[k]) + hh_len[k]; ptr < end; ptr++)
            x[sv_ind[ptr]] -= sv_val[ptr] * x_j;
      }
      return;
}

/* eof */
