/* sgf.c (sparse Gaussian factorizer) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2012-2015 Free Software Foundation, Inc.
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
#include "sgf.h"

/***********************************************************************
*  sgf_reduce_nuc - initial reordering to minimize nucleus size
*
*  On entry to this routine it is assumed that V = A and F = P = Q = I,
*  where A is the original matrix to be factorized. It is also assumed
*  that matrix V = A is stored in both row- and column-wise formats.
*
*  This routine performs (implicit) non-symmetric permutations of rows
*  and columns of matrix U = P'* V * Q' to reduce it to the form:
*
*        1     k1    k2    n
*     1  x x x x x x x x x x
*        . x x x x x x x x x
*        . . x x x x x x x x
*     k1 . . . * * * * x x x
*        . . . * * * * x x x
*        . . . * * * * x x x
*     k2 . . . * * * * x x x
*        . . . . . . . x x x
*        . . . . . . . . x x
*     n  . . . . . . . . . x
*
*  where non-zeros in rows and columns k1, k1+1, ..., k2 constitute so
*  called nucleus ('*'), whose size is minimized by the routine.
*
*  The numbers k1 and k2 are returned by the routine on exit. Usually,
*  if the nucleus exists, 1 <= k1 < k2 <= n. However, if the resultant
*  matrix U is upper triangular (has no nucleus), k1 = n+1 and k2 = n.
*
*  Note that the routines sgf_choose_pivot and sgf_eliminate perform
*  exactly the same transformations (by processing row and columns
*  singletons), so preliminary minimization of the nucleus may not be
*  used. However, processing row and column singletons by the routines
*  sgf_minimize_nuc and sgf_singl_phase is more efficient. */

#if 1 /* 21/II-2016 */
/* Normally this routine returns zero. If the matrix is structurally
*  singular, the routine returns non-zero. */
#endif

int sgf_reduce_nuc(LUF *luf, int *k1_, int *k2_, int cnt[/*1+n*/],
      int list[/*1+n*/])
{     int n = luf->n;
      SVA *sva = luf->sva;
      int *sv_ind = sva->ind;
      int vr_ref = luf->vr_ref;
      int *vr_ptr = &sva->ptr[vr_ref-1];
      int *vr_len = &sva->len[vr_ref-1];
      int vc_ref = luf->vc_ref;
      int *vc_ptr = &sva->ptr[vc_ref-1];
      int *vc_len = &sva->len[vc_ref-1];
      int *pp_ind = luf->pp_ind;
      int *pp_inv = luf->pp_inv;
      int *qq_ind = luf->qq_ind;
      int *qq_inv = luf->qq_inv;
      int i, ii, j, jj, k1, k2, ns, ptr, end;
      /* initial nucleus is U = V = A */
      k1 = 1, k2 = n;
      /*--------------------------------------------------------------*/
      /* process column singletons                                    */
      /*--------------------------------------------------------------*/
      /* determine initial counts of columns of V and initialize list
       * of active column singletons */
      ns = 0; /* number of active column singletons */
      for (j = 1; j <= n; j++)
      {  if ((cnt[j] = vc_len[j]) == 1)
            list[++ns] = j;
      }
      /* process active column singletons */
      while (ns > 0)
      {  /* column singleton is in j-th column of V */
         j = list[ns--];
#if 1 /* 21/II-2016 */
         if (cnt[j] == 0)
         {  /* j-th column in the current nucleus is actually empty */
            /* this happened because on a previous step in the nucleus
             * there were two or more identical column singletons (that
             * means structural singularity), so removing one of them
             * from the nucleus made other columns empty */
            return 1;
         }
#endif
         /* find i-th row of V containing column singleton */
         ptr = vc_ptr[j];
         end = ptr + vc_len[j];
         for (; pp_ind[i = sv_ind[ptr]] < k1; ptr++)
            /* nop */;
         xassert(ptr < end);
         /* permute rows and columns of U to move column singleton to
          * position u[k1,k1] */
         ii = pp_ind[i];
         luf_swap_u_rows(k1, ii);
         jj = qq_inv[j];
         luf_swap_u_cols(k1, jj);
         /* nucleus size decreased */
         k1++;
         /* walk thru i-th row of V and decrease column counts; this
          * may cause new column singletons to appear */
         ptr = vr_ptr[i];
         end = ptr + vr_len[i];
         for (; ptr < end; ptr++)
         {  if (--(cnt[j = sv_ind[ptr]]) == 1)
               list[++ns] = j;
         }
      }
      /* nucleus begins at k1-th row/column of U */
      if (k1 > n)
      {  /* U is upper triangular; no nucleus exist */
         goto done;
      }
      /*--------------------------------------------------------------*/
      /* process row singletons                                       */
      /*--------------------------------------------------------------*/
      /* determine initial counts of rows of V and initialize list of
       * active row singletons */
      ns = 0; /* number of active row singletons */
      for (i = 1; i <= n; i++)
      {  if (pp_ind[i] < k1)
         {  /* corresponding row of U is above its k1-th row; set its
             * count to zero to prevent including it in active list */
            cnt[i] = 0;
         }
         else if ((cnt[i] = vr_len[i]) == 1)
            list[++ns] = i;
      }
      /* process active row singletons */
      while (ns > 0)
      {  /* row singleton is in i-th row of V */
         i = list[ns--];
#if 1 /* 21/II-2016 */
         if (cnt[i] == 0)
         {  /* i-th row in the current nucleus is actually empty */
            /* (see comments above for similar case of empty column) */
            return 2;
         }
#endif
         /* find j-th column of V containing row singleton */
         ptr = vr_ptr[i];
         end = ptr + vr_len[i];
         for (; qq_inv[j = sv_ind[ptr]] > k2; ptr++)
            /* nop */;
         xassert(ptr < end);
         /* permute rows and columns of U to move row singleton to
          * position u[k2,k2] */
         ii = pp_ind[i];
         luf_swap_u_rows(k2, ii);
         jj = qq_inv[j];
         luf_swap_u_cols(k2, jj);
         /* nucleus size decreased */
         k2--;
         /* walk thru j-th column of V and decrease row counts; this
          * may cause new row singletons to appear */
         ptr = vc_ptr[j];
         end = ptr + vc_len[j];
         for (; ptr < end; ptr++)
         {  if (--(cnt[i = sv_ind[ptr]]) == 1)
               list[++ns] = i;
         }
      }
      /* nucleus ends at k2-th row/column of U */
      xassert(k1 < k2);
done: *k1_ = k1, *k2_ = k2;
      return 0;
}

/***********************************************************************
*  sgf_singl_phase - compute LU-factorization (singleton phase)
*
*  It is assumed that on entry to the routine L = P'* F * P = F = I
*  and matrix U = P'* V * Q' has the following structure (provided by
*  the routine sgf_reduce_nuc):
*
*        1     k1    k2    n
*     1  a a a b b b b c c c
*        . a a b b b b c c c
*        . . a b b b b c c c
*     k1 . . . * * * * d d d
*        . . . * * * * d d d
*        . . . * * * * d d d
*     k2 . . . * * * * d d d
*        . . . . . . . e e e
*        . . . . . . . . e e
*     n  . . . . . . . . . e
*
*  First, the routine performs (implicit) symmetric permutations of
*  rows and columns of matrix U to place them in the following order:
*
*     1, 2, ..., k1-1; n, n-1, ..., k2+1; k1, k1+1, ..., k2
*
*  This changes the structure of matrix U as follows:
*
*        1     k1    k2'   n
*     1  a a a c c c b b b b
*        . a a c c c b b b b
*        . . a c c c b b b b
*     k1 . . . e . . . . . .
*        . . . e e . . . . .
*        . . . e e e . . . .
*     k2'. . . d d d * * * *
*        . . . d d d * * * *
*        . . . d d d * * * *
*     n  . . . d d d * * * *
*
*  where k2' = n - k2 + k1.
*
*  Then the routine performs elementary gaussian transformations to
*  eliminate subdiagonal elements in columns k1, ..., k2'-1 of U. The
*  effect is the same as if the routine sgf_eliminate would be called
*  for k = 1, ..., k2'-1 using diagonal elements u[k,k] as pivots.
*
*  After elimination matrices L and U becomes the following:
*
*        1     k1   k2'    n        1     k1   k2'    n
*     1  1 . . . . . . . . .     1  a a a c c c b b b b
*        . 1 . . . . . . . .        . a a c c c b b b b
*        . . 1 . . . . . . .        . . a c c c b b b b
*     k1 . . . 1 . . . . . .     k1 . . . e . . . . . .
*        . . . e'1 . . . . .        . . . . e . . . . .
*        . . . e'e'1 . . . .        . . . . . e . . . .
*     k2'. . . d'd'd'1 . . .     k2'. . . . . . * * * *
*        . . . d'd'd'. 1 . .        . . . . . . * * * *
*        . . . d'd'd'. . 1 .        . . . . . . * * * *
*     n  . . . d'd'd'. . . 1     n  . . . . . . * * * *
*
*             matrix L                   matrix U
*
*  where columns k1, ..., k2'-1 of L consist of subdiagonal elements
*  of initial matrix U divided by pivots u[k,k].
*
*  On exit the routine returns k2', the elimination step number, from
*  which computing of the factorization should be continued. Note that
*  k2' = n+1 means that matrix U is already upper triangular. */

int sgf_singl_phase(LUF *luf, int k1, int k2, int updat,
      int ind[/*1+n*/], double val[/*1+n*/])
{     int n = luf->n;
      SVA *sva = luf->sva;
      int *sv_ind = sva->ind;
      double *sv_val = sva->val;
      int fc_ref = luf->fc_ref;
      int *fc_ptr = &sva->ptr[fc_ref-1];
      int *fc_len = &sva->len[fc_ref-1];
      int vr_ref = luf->vr_ref;
      int *vr_ptr = &sva->ptr[vr_ref-1];
      int *vr_len = &sva->len[vr_ref-1];
      double *vr_piv = luf->vr_piv;
      int vc_ref = luf->vc_ref;
      int *vc_ptr = &sva->ptr[vc_ref-1];
      int *vc_len = &sva->len[vc_ref-1];
      int *pp_ind = luf->pp_ind;
      int *pp_inv = luf->pp_inv;
      int *qq_ind = luf->qq_ind;
      int *qq_inv = luf->qq_inv;
      int i, j, k, ptr, ptr1, end, len;
      double piv;
      /* (see routine sgf_reduce_nuc) */
      xassert((1 <= k1 && k1 < k2 && k2 <= n)
         || (k1 == n+1 && k2 == n));
      /* perform symmetric permutations of rows/columns of U */
      for (k = k1; k <= k2; k++)
         pp_ind[pp_inv[k]] = qq_inv[qq_ind[k]] = k - k2 + n;
      for (k = k2+1; k <= n; k++)
         pp_ind[pp_inv[k]] = qq_inv[qq_ind[k]] = n - k + k1;
      for (k = 1; k <= n; k++)
         pp_inv[pp_ind[k]] = qq_ind[qq_inv[k]] = k;
      /* determine k2' */
      k2 = n - k2 + k1;
      /* process rows and columns of V corresponding to rows and
       * columns 1, ..., k1-1 of U */
      for (k = 1; k < k1; k++)
      {  /* k-th row of U = i-th row of V */
         i = pp_inv[k];
         /* find pivot u[k,k] = v[i,j] in i-th row of V */
         ptr = vr_ptr[i];
         end = ptr + vr_len[i];
         for (; qq_inv[sv_ind[ptr]] != k; ptr++)
            /* nop */;
         xassert(ptr < end);
         /* store pivot */
         vr_piv[i] = sv_val[ptr];
         /* and remove it from i-th row of V */
         sv_ind[ptr] = sv_ind[end-1];
         sv_val[ptr] = sv_val[end-1];
         vr_len[i]--;
         /* clear column of V corresponding to k-th column of U */
         vc_len[qq_ind[k]] = 0;
      }
      /* clear rows of V corresponding to rows k1, ..., k2'-1 of U */
      for (k = k1; k < k2; k++)
         vr_len[pp_inv[k]] = 0;
      /* process rows and columns of V corresponding to rows and
       * columns k2', ..., n of U */
      for (k = k2; k <= n; k++)
      {  /* k-th row of U = i-th row of V */
         i = pp_inv[k];
         /* remove elements from i-th row of V that correspond to
          * elements u[k,k1], ..., u[k,k2'-1] */
         ptr = ptr1 = vr_ptr[i];
         end = ptr + vr_len[i];
         for (; ptr < end; ptr++)
         {  if (qq_inv[sv_ind[ptr]] >= k2)
            {  sv_ind[ptr1] = sv_ind[ptr];
               sv_val[ptr1] = sv_val[ptr];
               ptr1++;
            }
         }
         vr_len[i] = ptr1 - vr_ptr[i];
         /* k-th column of U = j-th column of V */
         j = qq_ind[k];
         /* remove elements from j-th column of V that correspond to
          * elements u[1,k], ..., u[k1-1,k] */
         ptr = ptr1 = vc_ptr[j];
         end = ptr + vc_len[j];
         for (; ptr < end; ptr++)
         {  if (pp_ind[sv_ind[ptr]] >= k2)
               /* element value is not needed in this case */
               sv_ind[ptr1++] = sv_ind[ptr];
         }
         vc_len[j] = ptr1 - vc_ptr[j];
      }
      /* process columns of V corresponding to columns k1, ..., k2'-1
       * of U, build columns of F */
      for (k = k1; k < k2; k++)
      {  /* k-th column of U = j-th column of V */
         j = qq_ind[k];
         /* remove elements from j-th column of V that correspond to
          * pivot (diagonal) element u[k,k] and subdiagonal elements
          * u[k+1,k], ..., u[n,k]; subdiagonal elements are stored for
          * further addition to matrix F */
         len = 0;
         piv = 0.0;
         ptr = vc_ptr[j];
         end = ptr + vc_len[j];
         for (; ptr < end; ptr++)
         {  i = sv_ind[ptr]; /* v[i,j] */
            if (pp_ind[i] == k)
            {  /* store pivot v[i,j] = u[k,k] */
               piv = vr_piv[i] = sv_val[ptr];
            }
            else if (pp_ind[i] > k)
            {  /* store subdiagonal element v[i,j] = u[i',k] */
               len++;
               ind[len] = i;
               val[len] = sv_val[ptr];
            }
         }
         /* clear j-th column of V = k-th column of U */
         vc_len[j] = 0;
         /* build k-th column of L = j-th column of F */
         j = pp_inv[k];
         xassert(piv != 0.0);
         if (len > 0)
         {  if (sva->r_ptr - sva->m_ptr < len)
            {  sva_more_space(sva, len);
               sv_ind = sva->ind;
               sv_val = sva->val;
            }
            sva_reserve_cap(sva, fc_ref-1+j, len);
            for (ptr = fc_ptr[j], ptr1 = 1; ptr1 <= len; ptr++, ptr1++)
            {  sv_ind[ptr] = ind[ptr1];
               sv_val[ptr] = val[ptr1] / piv;
            }
            fc_len[j] = len;
         }
      }
      /* if it is not planned to update matrix V, relocate all its
       * non-active rows corresponding to rows 1, ..., k2'-1 of U to
       * the right (static) part of SVA */
      if (!updat)
      {  for (k = 1; k < k2; k++)
         {  i = pp_inv[k];
            len = vr_len[i];
            if (sva->r_ptr - sva->m_ptr < len)
            {  sva_more_space(sva, len);
               sv_ind = sva->ind;
               sv_val = sva->val;
            }
            sva_make_static(sva, vr_ref-1+i);
         }
      }
      /* elimination steps 1, ..., k2'-1 have been performed */
      return k2;
}

/***********************************************************************
*  sgf_choose_pivot - choose pivot element v[p,q]
*
*  This routine chooses pivot element v[p,q], k <= p, q <= n, in the
*  active submatrix of matrix V = P * U * Q, where k is the number of
*  current elimination step, 1 <= k <= n.
*
*  It is assumed that on entry to the routine matrix U = P'* V * Q' has
*  the following partially triangularized form:
*
*        1       k         n
*     1  x x x x x x x x x x
*        . x x x x x x x x x
*        . . x x x x x x x x
*        . . . x x x x x x x
*     k  . . . . * * * * * *
*        . . . . * * * * * *
*        . . . . * * * * * *
*        . . . . * * * * * *
*        . . . . * * * * * *
*     n  . . . . * * * * * *
*
*  where rows and columns k, k+1, ..., n belong to the active submatrix
*  (its elements are marked by '*').
*
*  Since the matrix U is not stored, the routine works with the matrix
*  V = P * U * Q. It is assumed that the row-wise representation
*  corresponds to the matrix V, but the column-wise representation
*  corresponds to the active submatrix of the matrix V, i.e. elements,
*  which are not in the active submatrix, are not included in column
*  vectors. It is also assumed that each active row of the matrix V is
*  in the set R[len], where len is the number of non-zeros in the row,
*  and each active column of the matrix V is in the set C[len], where
*  len is the number of non-zeros in the column (in the latter case
*  only elements of the active submatrix are counted; such elements are
*  marked by '*' on the figure above).
*
*  For the reason of numerical stability the routine applies so called
*  threshold pivoting proposed by J.Reid. It is assumed that an element
*  v[i,j] can be selected as a pivot candidate if it is not very small
*  (in magnitude) among other elements in the same row, i.e. if it
*  satisfies to the stability condition |v[i,j]| >= tol * max|v[i,*]|,
*  where 0 < tol < 1 is a given tolerance.
*
*  In order to keep sparsity of the matrix V the routine uses Markowitz
*  strategy, trying to choose such element v[p,q], which satisfies to
*  the stability condition (see above) and has smallest Markowitz cost
*  (nr[p]-1) * (nc[q]-1), where nr[p] and nc[q] are, resp., numbers of
*  non-zeros in p-th row and q-th column of the active submatrix.
*
*  In order to reduce the search, i.e. not to walk through all elements
*  of the active submatrix, the routine uses a technique proposed by
*  I.Duff. This technique is based on using the sets R[len] and C[len]
*  of active rows and columns.
*
*  If the pivot element v[p,q] has been chosen, the routine stores its
*  indices to locations *p and *q and returns zero. Otherwise, non-zero
*  is returned. */

int sgf_choose_pivot(SGF *sgf, int *p_, int *q_)
{     LUF *luf = sgf->luf;
      int n = luf->n;
      SVA *sva = luf->sva;
      int *sv_ind = sva->ind;
      double *sv_val = sva->val;
      int vr_ref = luf->vr_ref;
      int *vr_ptr = &sva->ptr[vr_ref-1];
      int *vr_len = &sva->len[vr_ref-1];
      int vc_ref = luf->vc_ref;
      int *vc_ptr = &sva->ptr[vc_ref-1];
      int *vc_len = &sva->len[vc_ref-1];
      int *rs_head = sgf->rs_head;
      int *rs_next = sgf->rs_next;
      int *cs_head = sgf->cs_head;
      int *cs_prev = sgf->cs_prev;
      int *cs_next = sgf->cs_next;
      double *vr_max = sgf->vr_max;
      double piv_tol = sgf->piv_tol;
      int piv_lim = sgf->piv_lim;
      int suhl = sgf->suhl;
      int i, i_ptr, i_end, j, j_ptr, j_end, len, min_i, min_j, min_len,
         ncand, next_j, p, q;
      double best, big, cost, temp;
      /* no pivot candidate has been chosen so far */
      p = q = 0, best = DBL_MAX, ncand = 0;
      /* if the active submatrix contains a column having the only
       * non-zero element (column singleton), choose it as the pivot */
      j = cs_head[1];
      if (j != 0)
      {  xassert(vc_len[j] == 1);
         p = sv_ind[vc_ptr[j]], q = j;
         goto done;
      }
      /* if the active submatrix contains a row having the only
       * non-zero element (row singleton), choose it as the pivot */
      i = rs_head[1];
      if (i != 0)
      {  xassert(vr_len[i] == 1);
         p = i, q = sv_ind[vr_ptr[i]];
         goto done;
      }
      /* the active submatrix contains no singletons; walk thru its
       * other non-empty rows and columns */
      for (len = 2; len <= n; len++)
      {  /* consider active columns containing len non-zeros */
         for (j = cs_head[len]; j != 0; j = next_j)
         {  /* save the number of next column of the same length */
            next_j = cs_next[j];
            /* find an element in j-th column, which is placed in the
             * row with minimal number of non-zeros and satisfies to
             * the stability condition (such element may not exist) */
            min_i = min_j = 0, min_len = INT_MAX;
            for (j_end = (j_ptr = vc_ptr[j]) + vc_len[j];
               j_ptr < j_end; j_ptr++)
            {  /* get row index of v[i,j] */
               i = sv_ind[j_ptr];
               /* if i-th row is not shorter, skip v[i,j] */
               if (vr_len[i] >= min_len)
                  continue;
               /* big := max|v[i,*]| */
               if ((big = vr_max[i]) < 0.0)
               {  /* largest magnitude is unknown; compute it */
                  for (i_end = (i_ptr = vr_ptr[i]) + vr_len[i];
                     i_ptr < i_end; i_ptr++)
                  {  if ((temp = sv_val[i_ptr]) < 0.0)
                        temp = -temp;
                     if (big < temp)
                        big = temp;
                  }
                  xassert(big > 0.0);
                  vr_max[i] = big;
               }
               /* find v[i,j] in i-th row */
               for (i_end = (i_ptr = vr_ptr[i]) + vr_len[i];
                  sv_ind[i_ptr] != j; i_ptr++)
                  /* nop */;
               xassert(i_ptr < i_end);
               /* if |v[i,j]| < piv_tol * max|v[i,*]|, skip v[i,j] */
               if ((temp = sv_val[i_ptr]) < 0.0)
                  temp = -temp;
               if (temp < piv_tol * big)
                  continue;
               /* v[i,j] is a better candidate */
               min_i = i, min_j = j, min_len = vr_len[i];
               /* if Markowitz cost of v[i,j] is not greater than
                * (len-1)**2, v[i,j] can be chosen as the pivot right
                * now; this heuristic reduces the search and works well
                * in many cases */
               if (min_len <= len)
               {  p = min_i, q = min_j;
                  goto done;
               }
            }
            /* j-th column has been scanned */
            if (min_i != 0)
            {  /* element v[min_i,min_j] is a next pivot candidate */
               ncand++;
               /* compute its Markowitz cost */
               cost = (double)(min_len - 1) * (double)(len - 1);
               /* if this element is better, choose it as the pivot */
               if (cost < best)
                  p = min_i, q = min_j, best = cost;
               /* if piv_lim candidates were considered, terminate
                * the search, because it is doubtful that a much better
                * candidate will be found */
               if (ncand == piv_lim)
                  goto done;
            }
            else if (suhl)
            {  /* j-th column has no eligible elements that satisfy to
                * the stability criterion; Uwe Suhl suggests to exclude
                * such column from further considerations until it
                * becomes a column singleton; in hard cases this may
                * significantly reduce the time needed to choose the
                * pivot element */
               sgf_deactivate_col(j);
               cs_prev[j] = cs_next[j] = j;
            }
         }
         /* consider active rows containing len non-zeros */
         for (i = rs_head[len]; i != 0; i = rs_next[i])
         {  /* big := max|v[i,*]| */
            if ((big = vr_max[i]) < 0.0)
            {  /* largest magnitude is unknown; compute it */
               for (i_end = (i_ptr = vr_ptr[i]) + vr_len[i];
                  i_ptr < i_end; i_ptr++)
               {  if ((temp = sv_val[i_ptr]) < 0.0)
                     temp = -temp;
                  if (big < temp)
                     big = temp;
               }
               xassert(big > 0.0);
               vr_max[i] = big;
            }
            /* find an element in i-th row, which is placed in the
             * column with minimal number of non-zeros and satisfies to
             * the stability condition (such element always exists) */
            min_i = min_j = 0, min_len = INT_MAX;
            for (i_end = (i_ptr = vr_ptr[i]) + vr_len[i];
               i_ptr < i_end; i_ptr++)
            {  /* get column index of v[i,j] */
               j = sv_ind[i_ptr];
               /* if j-th column is not shorter, skip v[i,j] */
               if (vc_len[j] >= min_len)
                  continue;
               /* if |v[i,j]| < piv_tol * max|v[i,*]|, skip v[i,j] */
               if ((temp = sv_val[i_ptr]) < 0.0)
                  temp = -temp;
               if (temp < piv_tol * big)
                  continue;
               /* v[i,j] is a better candidate */
               min_i = i, min_j = j, min_len = vc_len[j];
               /* if Markowitz cost of v[i,j] is not greater than
                * (len-1)**2, v[i,j] can be chosen as the pivot right
                * now; this heuristic reduces the search and works well
                * in many cases */
               if (min_len <= len)
               {  p = min_i, q = min_j;
                  goto done;
               }
            }
            /* i-th row has been scanned */
            if (min_i != 0)
            {  /* element v[min_i,min_j] is a next pivot candidate */
               ncand++;
               /* compute its Markowitz cost */
               cost = (double)(len - 1) * (double)(min_len - 1);
               /* if this element is better, choose it as the pivot */
               if (cost < best)
                  p = min_i, q = min_j, best = cost;
               /* if piv_lim candidates were considered, terminate
                * the search, because it is doubtful that a much better
                * candidate will be found */
               if (ncand == piv_lim)
                  goto done;
            }
            else
            {  /* this can never be */
               xassert(min_i != min_i);
            }
         }
      }
done: /* report the pivot to the factorization routine */
      *p_ = p, *q_ = q;
      return (p == 0);
}

/***********************************************************************
*  sgf_eliminate - perform gaussian elimination
*
*  This routine performs elementary gaussian transformations in order
*  to eliminate subdiagonal elements in k-th column of matrix
*  U = P'* V * Q' using pivot element u[k,k], where k is the number of
*  current elimination step, 1 <= k <= n.
*
*  The parameters p and q specify, resp., row and column indices of the
*  pivot element v[p,q] = u[k,k].
*
*  On entry the routine assumes that partially triangularized matrices
*  L = P'* F * P and U = P'* V * Q' have the following structure:
*
*        1       k         n       1        k         n
*     1  1 . . . . . . . . .     1  x x x x x x x x x x
*        x 1 . . . . . . . .        . x x x x x x x x x
*        x x 1 . . . . . . .        . . x x x x x x x x
*        x x x 1 . . . . . .        . . . x x x x x x x
*     k  x x x x 1 . . . . .     k  . . . . * * * * * *
*        x x x x _ 1 . . . .        . . . . # * * * * *
*        x x x x _ . 1 . . .        . . . . # * * * * *
*        x x x x _ . . 1 . .        . . . . # * * * * *
*        x x x x _ . . . 1 .        . . . . # * * * * *
*     n  x x x x _ . . . . 1     n  . . . . # * * * * *
*
*             matrix L                   matrix U
*
*  where rows and columns k, k+1, ..., n of matrix U constitute the
*  active submatrix. Elements to be eliminated are marked by '#', and
*  other elements of the active submatrix are marked by '*'. May note
*  that each eliminated non-zero element u[i,k] of matrix U gives
*  corresponding non-zero element l[i,k] of matrix L (marked by '_').
*
*  Actually all operations are performed on matrix V. It is assumed
*  that the row-wise representation corresponds to matrix V, but the
*  column-wise representation corresponds to the active submatrix of
*  matrix V (or, more precisely, to its pattern, because only row
*  indices for columns of the active submatrix are used on this stage).
*
*  Let u[k,k] = v[p,q] be the pivot. In order to eliminate subdiagonal
*  elements u[i',k] = v[i,q], i'= k+1, k+2, ..., n, the routine applies
*  the following elementary gaussian transformations:
*
*     (i-th row of V) := (i-th row of V) - f[i,p] * (p-th row of V),
*
*  where f[i,p] = v[i,q] / v[p,q] is a gaussian multiplier stored to
*  p-th column of matrix F to keep the main equality A = F * V
*  (corresponding elements l[i',k] of matrix L are marked by '_' on the
*  figure above).
*
*  NOTE: On entry to the routine the working arrays flag and work
*        should contain zeros. This status is retained by the routine
*        on exit. */

int sgf_eliminate(SGF *sgf, int p, int q)
{     LUF *luf = sgf->luf;
      int n = luf->n;
      SVA *sva = luf->sva;
      int *sv_ind = sva->ind;
      double *sv_val = sva->val;
      int fc_ref = luf->fc_ref;
      int *fc_ptr = &sva->ptr[fc_ref-1];
      int *fc_len = &sva->len[fc_ref-1];
      int vr_ref = luf->vr_ref;
      int *vr_ptr = &sva->ptr[vr_ref-1];
      int *vr_len = &sva->len[vr_ref-1];
      int *vr_cap = &sva->cap[vr_ref-1];
      double *vr_piv = luf->vr_piv;
      int vc_ref = luf->vc_ref;
      int *vc_ptr = &sva->ptr[vc_ref-1];
      int *vc_len = &sva->len[vc_ref-1];
      int *vc_cap = &sva->cap[vc_ref-1];
      int *rs_head = sgf->rs_head;
      int *rs_prev = sgf->rs_prev;
      int *rs_next = sgf->rs_next;
      int *cs_head = sgf->cs_head;
      int *cs_prev = sgf->cs_prev;
      int *cs_next = sgf->cs_next;
      double *vr_max = sgf->vr_max;
      char *flag = sgf->flag;
      double *work = sgf->work;
      double eps_tol = sgf->eps_tol;
      int nnz_diff = 0;
      int fill, i, i_ptr, i_end, j, j_ptr, j_end, ptr, len, loc, loc1;
      double vpq, fip, vij;
      xassert(1 <= p && p <= n);
      xassert(1 <= q && q <= n);
      /* remove p-th row from the active set; this row will never
       * return there */
      sgf_deactivate_row(p);
      /* process p-th (pivot) row */
      ptr = 0;
      for (i_end = (i_ptr = vr_ptr[p]) + vr_len[p];
         i_ptr < i_end; i_ptr++)
      {  /* get column index of v[p,j] */
         j = sv_ind[i_ptr];
         if (j == q)
         {  /* save pointer to pivot v[p,q] */
            ptr = i_ptr;
         }
         else
         {  /* store v[p,j], j != q, to working array */
            flag[j] = 1;
            work[j] = sv_val[i_ptr];
         }
         /* remove j-th column from the active set; q-th column will
          * never return there while other columns will return to the
          * active set with new length */
         if (cs_next[j] == j)
         {  /* j-th column was marked by the pivoting routine according
             * to Uwe Suhl's suggestion and is already inactive */
            xassert(cs_prev[j] == j);
         }
         else
            sgf_deactivate_col(j);
         nnz_diff -= vc_len[j];
         /* find and remove v[p,j] from j-th column */
         for (j_end = (j_ptr = vc_ptr[j]) + vc_len[j];
            sv_ind[j_ptr] != p; j_ptr++)
            /* nop */;
         xassert(j_ptr < j_end);
         sv_ind[j_ptr] = sv_ind[j_end-1];
         vc_len[j]--;
      }
      /* save pivot v[p,q] and remove it from p-th row */
      xassert(ptr > 0);
      vpq = vr_piv[p] = sv_val[ptr];
      sv_ind[ptr] = sv_ind[i_end-1];
      sv_val[ptr] = sv_val[i_end-1];
      vr_len[p]--;
      /* if it is not planned to update matrix V, relocate p-th row to
       * the right (static) part of SVA */
      if (!sgf->updat)
      {  len = vr_len[p];
         if (sva->r_ptr - sva->m_ptr < len)
         {  sva_more_space(sva, len);
            sv_ind = sva->ind;
            sv_val = sva->val;
         }
         sva_make_static(sva, vr_ref-1+p);
      }
      /* copy the pattern (row indices) of q-th column of the active
       * submatrix (from which v[p,q] has been just removed) to p-th
       * column of matrix F (without unity diagonal element) */
      len = vc_len[q];
      if (len > 0)
      {  if (sva->r_ptr - sva->m_ptr < len)
         {  sva_more_space(sva, len);
            sv_ind = sva->ind;
            sv_val = sva->val;
         }
         sva_reserve_cap(sva, fc_ref-1+p, len);
         memcpy(&sv_ind[fc_ptr[p]], &sv_ind[vc_ptr[q]],
            len * sizeof(int));
         fc_len[p] = len;
      }
      /* make q-th column of the active submatrix empty */
      vc_len[q] = 0;
      /* transform non-pivot rows of the active submatrix */
      for (loc = fc_len[p]-1; loc >= 0; loc--)
      {  /* get row index of v[i,q] = row index of f[i,p] */
         i = sv_ind[fc_ptr[p] + loc];
         xassert(i != p); /* v[p,q] was removed */
         /* remove i-th row from the active set; this row will return
          * there with new length */
         sgf_deactivate_row(i);
         /* find v[i,q] in i-th row */
         for (i_end = (i_ptr = vr_ptr[i]) + vr_len[i];
            sv_ind[i_ptr] != q; i_ptr++)
            /* nop */;
         xassert(i_ptr < i_end);
         /* compute gaussian multiplier f[i,p] = v[i,q] / v[p,q] */
         fip = sv_val[fc_ptr[p] + loc] = sv_val[i_ptr] / vpq;
         /* remove v[i,q] from i-th row */
         sv_ind[i_ptr] = sv_ind[i_end-1];
         sv_val[i_ptr] = sv_val[i_end-1];
         vr_len[i]--;
         /* perform elementary gaussian transformation:
          * (i-th row) := (i-th row) - f[i,p] * (p-th row)
          * note that p-th row of V, which is in the working array,
          * doesn't contain pivot v[p,q], and i-th row of V doesn't
          * contain v[i,q] to be eliminated */
         /* walk thru i-th row and transform existing elements */
         fill = vr_len[p];
         for (i_end = (i_ptr = ptr = vr_ptr[i]) + vr_len[i];
            i_ptr < i_end; i_ptr++)
         {  /* get column index and value of v[i,j] */
            j = sv_ind[i_ptr];
            vij = sv_val[i_ptr];
            if (flag[j])
            {  /* v[p,j] != 0 */
               flag[j] = 0, fill--;
               /* v[i,j] := v[i,j] - f[i,p] * v[p,j] */
               vij -= fip * work[j];
               if (-eps_tol < vij && vij < +eps_tol)
               {  /* new v[i,j] is close to zero; remove it from the
                   * active submatrix, i.e. replace it by exact zero */
                  /* find and remove v[i,j] from j-th column */
                  for (j_end = (j_ptr = vc_ptr[j]) + vc_len[j];
                     sv_ind[j_ptr] != i; j_ptr++)
                     /* nop */;
                  xassert(j_ptr < j_end);
                  sv_ind[j_ptr] = sv_ind[j_end-1];
                  vc_len[j]--;
                  continue;
               }
            }
            /* keep new v[i,j] in i-th row */
            sv_ind[ptr] = j;
            sv_val[ptr] = vij;
            ptr++;
         }
         /* (new length of i-th row may decrease because of numerical
          * cancellation) */
         vr_len[i] = len = ptr - vr_ptr[i];
         /* now flag[*] is the pattern of the set v[p,*] \ v[i,*], and
          * fill is the number of non-zeros in this set */
         if (fill == 0)
         {  /* no fill-in occurs */
            /* walk thru p-th row and restore the column flags */
            for (i_end = (i_ptr = vr_ptr[p]) + vr_len[p];
               i_ptr < i_end; i_ptr++)
               flag[sv_ind[i_ptr]] = 1; /* v[p,j] != 0 */
            goto skip;
         }
         /* up to fill new non-zero elements may appear in i-th row due
          * to fill-in; reserve locations for these elements (note that
          * actual length of i-th row is currently stored in len) */
         if (vr_cap[i] < len + fill)
         {  if (sva->r_ptr - sva->m_ptr < len + fill)
            {  sva_more_space(sva, len + fill);
               sv_ind = sva->ind;
               sv_val = sva->val;
            }
            sva_enlarge_cap(sva, vr_ref-1+i, len + fill, 0);
         }
         vr_len[i] += fill;
         /* walk thru p-th row and add new elements to i-th row */
         for (loc1 = vr_len[p]-1; loc1 >= 0; loc1--)
         {  /* get column index of v[p,j] */
            j = sv_ind[vr_ptr[p] + loc1];
            if (!flag[j])
            {  /* restore j-th column flag */
               flag[j] = 1;
               /* v[i,j] was computed earlier on transforming existing
                * elements of i-th row */
               continue;
            }
            /* v[i,j] := 0 - f[i,p] * v[p,j] */
            vij = - fip * work[j];
            if (-eps_tol < vij && vij < +eps_tol)
            {  /* new v[i,j] is close to zero; do not add it to the
                * active submatrix, i.e. replace it by exact zero */
               continue;
            }
            /* add new v[i,j] to i-th row */
            sv_ind[ptr = vr_ptr[i] + (len++)] = j;
            sv_val[ptr] = vij;
            /* add new v[i,j] to j-th column */
            if (vc_cap[j] == vc_len[j])
            {  /* we reserve extra locations in j-th column to reduce
                * further relocations of that column */
#if 1 /* FIXME */
               /* use control parameter to specify the number of extra
                * locations reserved */
               int need = vc_len[j] + 10;
#endif
               if (sva->r_ptr - sva->m_ptr < need)
               {  sva_more_space(sva, need);
                  sv_ind = sva->ind;
                  sv_val = sva->val;
               }
               sva_enlarge_cap(sva, vc_ref-1+j, need, 1);
            }
            sv_ind[vc_ptr[j] + (vc_len[j]++)] = i;
         }
         /* set final length of i-th row just transformed */
         xassert(len <= vr_len[i]);
         vr_len[i] = len;
skip:    /* return i-th row to the active set with new length */
         sgf_activate_row(i);
         /* since i-th row has been changed, largest magnitude of its
          * elements becomes unknown */
         vr_max[i] = -1.0;
      }
      /* walk thru p-th (pivot) row */
      for (i_end = (i_ptr = vr_ptr[p]) + vr_len[p];
         i_ptr < i_end; i_ptr++)
      {  /* get column index of v[p,j] */
         j = sv_ind[i_ptr];
         xassert(j != q); /* v[p,q] was removed */
         /* return j-th column to the active set with new length */
         if (cs_next[j] == j && vc_len[j] != 1)
         {  /* j-th column was marked by the pivoting routine and it is
             * still not a column singleton, so leave it incative */
            xassert(cs_prev[j] == j);
         }
         else
            sgf_activate_col(j);
         nnz_diff += vc_len[j];
         /* restore zero content of the working arrays */
         flag[j] = 0;
         work[j] = 0.0;
      }
      /* return the difference between the numbers of non-zeros in the
       * active submatrix on entry and on exit, resp. */
      return nnz_diff;
}

/***********************************************************************
*  sgf_dense_lu - compute dense LU-factorization with full pivoting
*
*  This routine performs Gaussian elimination with full pivoting to
*  compute dense LU-factorization of the specified matrix A of order n
*  in the form:
*
*     A = P * L * U * Q,                                             (1)
*
*  where L is lower triangular matrix with unit diagonal, U is upper
*  triangular matrix, P and Q are permutation matrices.
*
*  On entry to the routine elements of matrix A = (a[i,j]) should be
*  placed in the array elements a[0], ..., a[n^2-1] in dense row-wise
*  format. On exit from the routine matrix A is replaced by factors L
*  and U as follows:
*
*       u[1,1]   u[1,2]  ...  u[1,n-1]   u[1,n]
*       l[2,1]   u[2,2]  ...  u[2,n-1]   u[2,n]
*        . . . . . . . . . . . . . .
*     l[n-1,1] l[n-1,2]     u[n-1,n-1] u[n-1,n]
*       l[n,1]   l[n,2]  ...  l[n,n-1]   u[n,n]
*
*  The unit diagonal elements of L are not stored.
*
*  Information on permutations of rows and columns of active submatrix
*  during factorization is accumulated by the routine as follows. Every
*  time the routine permutes rows i and i' or columns j and j', it also
*  permutes elements r[i-1] and r[i'-1] or c[j-1] and c[j'-1], resp.
*  Thus, on entry to the routine elements r[0], r[1], ..., r[n-1] and
*  c[0], c[1], ..., c[n-1] should be initialized by some integers that
*  identify rows and columns of the original matrix A.
*
*  If the factorization has been successfully computed, the routine
*  returns zero. Otherwise, if on k-th elimination step, 1 <= k <= n,
*  all elements of the active submatrix are close to zero, the routine
*  returns k, in which case a partial factorization is stored in the
*  array a. */

int sgf_dense_lu(int n, double a_[], int r[], int c[], double eps)
{     /* non-optimized version */
      int i, j, k, p, q, ref;
      double akk, big, temp;
#     define a(i,j) a_[(i)*n+(j)]
      /* initially U = A, L = P = Q = I */
      /* main elimination loop */
      for (k = 0; k < n; k++)
      {  /* choose pivot u[p,q], k <= p, q <= n */
         p = q = -1, big = eps;
         for (i = k; i < n; i++)
         {  for (j = k; j < n; j++)
            {  /* temp = |u[i,j]| */
               if ((temp = a(i,j)) < 0.0)
                  temp = -temp;
               if (big < temp)
                  p = i, q = j, big = temp;
            }
         }
         if (p < 0)
         {  /* k-th elimination step failed */
            return k+1;
         }
         /* permute rows k and p */
         if (k != p)
         {  for (j = 0; j < n; j++)
               temp = a(k,j), a(k,j) = a(p,j), a(p,j) = temp;
            ref = r[k], r[k] = r[p], r[p] = ref;
         }
         /* permute columns k and q */
         if (k != q)
         {  for (i = 0; i < n; i++)
               temp = a(i,k), a(i,k) = a(i,q), a(i,q) = temp;
            ref = c[k], c[k] = c[q], c[q] = ref;
         }
         /* now pivot is in position u[k,k] */
         akk = a(k,k);
         /* eliminate subdiagonal elements u[k+1,k], ..., u[n,k] */
         for (i = k+1; i < n; i++)
         {  if (a(i,k) != 0.0)
            {  /* gaussian multiplier l[i,k] := u[i,k] / u[k,k] */
               temp = (a(i,k) /= akk);
               /* (i-th row) := (i-th row) - l[i,k] * (k-th row) */
               for (j = k+1; j < n; j++)
                  a(i,j) -= temp * a(k,j);
            }
         }
      }
#     undef a
      return 0;
}

/***********************************************************************
*  sgf_dense_phase - compute LU-factorization (dense phase)
*
*  This routine performs dense phase of computing LU-factorization.
*
*  The aim is two-fold. First, the main factorization routine switches
*  to dense phase when the active submatrix is relatively dense, so
*  using dense format allows significantly reduces overheads needed to
*  maintain sparse data structures. And second, that is more important,
*  on dense phase full pivoting is used (rather than partial pivoting)
*  that allows improving numerical stability, since round-off errors
*  tend to increase on last steps of the elimination process.
*
*  On entry the routine assumes that elimination steps 1, 2, ..., k-1
*  have been performed, so partially transformed matrices L = P'* F * P
*  and U = P'* V * Q' have the following structure:
*
*        1       k         n       1        k         n
*     1  1 . . . . . . . . .     1  x x x x x x x x x x
*        x 1 . . . . . . . .        . x x x x x x x x x
*        x x 1 . . . . . . .        . . x x x x x x x x
*        x x x 1 . . . . . .        . . . x x x x x x x
*     k  x x x x 1 . . . . .     k  . . . . * * * * * *
*        x x x x . 1 . . . .        . . . . * * * * * *
*        x x x x . . 1 . . .        . . . . * * * * * *
*        x x x x . . . 1 . .        . . . . * * * * * *
*        x x x x . . . . 1 .        . . . . * * * * * *
*     n  x x x x . . . . . 1     n  . . . . * * * * * *
*
*             matrix L                   matrix U
*
*  where rows and columns k, k+1, ..., n of matrix U constitute the
*  active submatrix A~, whose elements are marked by '*'.
*
*  The routine copies the active submatrix A~ to a working array in
*  dense format, compute dense factorization A~ = P~* L~* U~* Q~ using
*  full pivoting, and then copies non-zero elements of factors L~ and
*  U~ back to factors L and U (more precisely, to factors F and V).
*
*  If the factorization has been successfully computed, the routine
*  returns zero. Otherwise, if on k-th elimination step, 1 <= k <= n,
*  all elements of the active submatrix are close to zero, the routine
*  returns k (information on linearly dependent rows/columns in this
*  case is provided by matrices P and Q). */

int sgf_dense_phase(LUF *luf, int k, int updat)
{     int n = luf->n;
      SVA *sva = luf->sva;
      int *sv_ind = sva->ind;
      double *sv_val = sva->val;
      int fc_ref = luf->fc_ref;
      int *fc_ptr = &sva->ptr[fc_ref-1];
      int *fc_len = &sva->len[fc_ref-1];
      int *fc_cap = &sva->cap[fc_ref-1];
      int vr_ref = luf->vr_ref;
      int *vr_ptr = &sva->ptr[vr_ref-1];
      int *vr_len = &sva->len[vr_ref-1];
      int *vr_cap = &sva->cap[vr_ref-1];
      double *vr_piv = luf->vr_piv;
      int vc_ref = luf->vc_ref;
      int *vc_len = &sva->len[vc_ref-1];
      int *pp_inv = luf->pp_inv;
      int *pp_ind = luf->pp_ind;
      int *qq_ind = luf->qq_ind;
      int *qq_inv = luf->qq_inv;
      int a_end, a_ptr, end, i, ia, ii, j, ja, jj, ka, len, na, ne,
         need, ptr;
      double *a_;
      xassert(1 <= k && k <= n);
      /* active columns of V are not longer needed; make them empty */
      for (jj = k; jj <= n; jj++)
      {  /* jj is number of active column of U = P'* V * Q' */
         vc_len[qq_ind[jj]] = 0;
      }
      /* determine order of active submatrix A~ of matrix U */
      na = n - k + 1;
      xassert(1 <= na && na <= n);
      /* determine number of elements in dense triangular factor (L~ or
       * U~), except diagonal elements */
      ne = na * (na - 1) / 2;
      /* we allocate active submatrix A~ in free (middle) part of SVA;
       * to avoid defragmentation that could destroy A~ we also should
       * reserve ne locations to build rows of V from rows of U~ and ne
       * locations to build columns of F from columns of L~ */
      need = na * na + ne + ne;
      if (sva->r_ptr - sva->m_ptr < need)
      {  sva_more_space(sva, need);
         sv_ind = sva->ind;
         sv_val = sva->val;
      }
      /* free (middle) part of SVA is structured as follows:
       * end of left (dynamic) part
       * ne free locations for new rows of V
       * na free locations for active submatrix A~
       * unused locations, if any
       * ne free locations for new columns of F
       * beginning of right (static) part */
      a_ptr = sva->m_ptr + ne;
      a_end = a_ptr + na * na;
      /* copy active submatrix A~ from matrix V to working array in
       * dense row-wise format */
      a_ = &sva->val[a_ptr];
#     define a(ia, ja) a_[((ia) - 1) * na + ((ja) - 1)]
      for (ia = 1; ia <= na; ia++)
      {  /* clear ia-th row of A~ */
         for (ja = 1; ja <= na; ja++)
            a(ia, ja) = 0.0;
         /* ia-th row of A~ = (k-1+ia)-th row of U = i-th row of V */
         i = pp_inv[k-1+ia];
         ptr = vr_ptr[i];
         end = ptr + vr_len[i];
         for (; ptr < end; ptr++)
            a(ia, qq_inv[sv_ind[ptr]]-k+1) = sv_val[ptr];
         /* i-th row of V is no longer needed; make it empty */
         vr_len[i] = 0;
      }
      /* compute dense factorization A~ = P~* L~* U~* Q~ */
#if 1 /* FIXME: epsilon tolerance */
      ka = sgf_dense_lu(na, &a(1, 1), &pp_inv[k], &qq_ind[k], 1e-20);
#endif
      /* rows of U with numbers pp_inv[k, k+1, ..., n] were permuted
       * due to row permutations of A~; update matrix P using P~ */
      for (ii = k; ii <= n; ii++)
         pp_ind[pp_inv[ii]] = ii;
      /* columns of U with numbers qq_ind[k, k+1, ..., n] were permuted
       * due to column permutations of A~; update matrix Q using Q~ */
      for (jj = k; jj <= n; jj++)
         qq_inv[qq_ind[jj]] = jj;
      /* check if dense factorization is complete */
      if (ka != 0)
      {  /* A~ is singular to working precision */
         /* information on linearly dependent rows/columns is provided
          * by matrices P and Q */
         xassert(1 <= ka && ka <= na);
         return k - 1 + ka;
      }
      /* build new rows of V from rows of U~ */
      for (ia = 1; ia <= na; ia++)
      {  /* ia-th row of U~ = (k-1+ia)-th row of U = i-th row of V */
         i = pp_inv[k-1+ia];
         xassert(vr_len[i] == 0);
         /* store diagonal element u~[ia,ia] */
         vr_piv[i] = a(ia, ia);
         /* determine number of non-zero non-diagonal elements in ia-th
          * row of U~ */
         len = 0;
         for (ja = ia+1; ja <= na; ja++)
         {  if (a(ia, ja) != 0.0)
               len++;
         }
         /* reserve len locations for i-th row of matrix V in left
          * (dynamic) part of SVA */
         if (vr_cap[i] < len)
         {  /* there should be enough room in free part of SVA */
            xassert(sva->r_ptr - sva->m_ptr >= len);
            sva_enlarge_cap(sva, vr_ref-1+i, len, 0);
            /* left part of SVA should not overlap matrix A~ */
            xassert(sva->m_ptr <= a_ptr);
         }
         /* copy non-zero non-diaginal elements of ia-th row of U~ to
          * i-th row of V */
         ptr = vr_ptr[i];
         for (ja = ia+1; ja <= na; ja++)
         {  if (a(ia, ja) != 0.0)
            {  sv_ind[ptr] = qq_ind[k-1+ja];
               sv_val[ptr] = a(ia, ja);
               ptr++;
            }
         }
         xassert(ptr - vr_ptr[i] == len);
         vr_len[i] = len;
      }
      /* build new columns of F from columns of L~ */
      for (ja = 1; ja <= na; ja++)
      {  /* ja-th column of L~ = (k-1+ja)-th column of L = j-th column
          * of F */
         j = pp_inv[k-1+ja];
         xassert(fc_len[j] == 0);
         xassert(fc_cap[j] == 0);
         /* determine number of non-zero non-diagonal elements in ja-th
          * column of L~ */
         len = 0;
         for (ia = ja+1; ia <= na; ia++)
         {  if (a(ia, ja) != 0.0)
               len++;
         }
         /* reserve len locations for j-th column of matrix F in right
          * (static) part of SVA */
         /* there should be enough room in free part of SVA */
         xassert(sva->r_ptr - sva->m_ptr >= len);
         if (len > 0)
            sva_reserve_cap(sva, fc_ref-1+j, len);
         /* right part of SVA should not overlap matrix A~ */
         xassert(a_end <= sva->r_ptr);
         /* copy non-zero non-diagonal elements of ja-th column of L~
          * to j-th column of F */
         ptr = fc_ptr[j];
         for (ia = ja+1; ia <= na; ia++)
         {  if (a(ia, ja) != 0.0)
            {  sv_ind[ptr] = pp_inv[k-1+ia];
               sv_val[ptr] = a(ia, ja);
               ptr++;
            }
         }
         xassert(ptr - fc_ptr[j] == len);
         fc_len[j] = len;
      }
      /* factors L~ and U~ are no longer needed */
#     undef a
      /* if it is not planned to update matrix V, relocate all its new
       * rows to the right (static) part of SVA */
      if (!updat)
      {  for (ia = 1; ia <= na; ia++)
         {  i = pp_inv[k-1+ia];
            len = vr_len[i];
            if (sva->r_ptr - sva->m_ptr < len)
            {  sva_more_space(sva, len);
               sv_ind = sva->ind;
               sv_val = sva->val;
            }
            sva_make_static(sva, vr_ref-1+i);
         }
      }
      return 0;
}

/***********************************************************************
*  sgf_factorize - compute LU-factorization (main routine)
*
*  This routine computes sparse LU-factorization of specified matrix A
*  using Gaussian elimination.
*
*  On entry to the routine matrix V = A should be stored in column-wise
*  format.
*
*  If the factorization has been successfully computed, the routine
*  returns zero. Otherwise, if on k-th elimination step, 1 <= k <= n,
*  all elements of the active submatrix are close to zero, the routine
*  returns k (information on linearly dependent rows/columns in this
*  case is provided by matrices P and Q). */

#if 1 /* 21/II-2016 */
/* If the matrix A is structurally singular, the routine returns -1.
*  NOTE: This case can be detected only if the singl flag is set. */
#endif

int sgf_factorize(SGF *sgf, int singl)
{     LUF *luf = sgf->luf;
      int n = luf->n;
      SVA *sva = luf->sva;
      int vr_ref = luf->vr_ref;
      int *vr_len = &sva->len[vr_ref-1];
      double *vr_piv = luf->vr_piv;
      int vc_ref = luf->vc_ref;
      int *vc_len = &sva->len[vc_ref-1];
      int *pp_ind = luf->pp_ind;
      int *pp_inv = luf->pp_inv;
      int *qq_ind = luf->qq_ind;
      int *qq_inv = luf->qq_inv;
      int *rs_head = sgf->rs_head;
      int *rs_prev = sgf->rs_prev;
      int *rs_next = sgf->rs_next;
      int *cs_head = sgf->cs_head;
      int *cs_prev = sgf->cs_prev;
      int *cs_next = sgf->cs_next;
      double *vr_max = sgf->vr_max;
      char *flag = sgf->flag;
      double *work = sgf->work;
      int i, j, k, k1, k2, p, q, nnz;
      /* build matrix V = A in row-wise format */
      luf_build_v_rows(luf, rs_prev);
      /* P := Q := I, so V = U = A, F = L = I */
      for (k = 1; k <= n; k++)
      {  vr_piv[k] = 0.0;
         pp_ind[k] = pp_inv[k] = qq_ind[k] = qq_inv[k] = k;
      }
#ifdef GLP_DEBUG
      sva_check_area(sva);
      luf_check_all(luf, 1);
#endif
      /* perform singleton phase, if required */
      if (!singl)
      {  /* assume that nucleus is entire matrix U */
         k2 = 1;
      }
      else
      {  /* minimize nucleus size */
#if 0 /* 21/II-2016 */
         sgf_reduce_nuc(luf, &k1, &k2, rs_prev, rs_next);
#else
         if (sgf_reduce_nuc(luf, &k1, &k2, rs_prev, rs_next))
            return -1;
#endif
#ifdef GLP_DEBUG
         xprintf("n = %d; k1 = %d; k2 = %d\n", n, k1, k2);
#endif
         /* perform singleton phase */
         k2 = sgf_singl_phase(luf, k1, k2, sgf->updat, rs_prev, work);
      }
#ifdef GLP_DEBUG
      sva_check_area(sva);
      luf_check_all(luf, k2);
#endif
      /* initialize working arrays */
      rs_head[0] = cs_head[0] = 0;
      for (k = 1; k <= n; k++)
      {  rs_head[k] = cs_head[k] = 0;
         vr_max[k] = -1.0;
         flag[k] = 0;
         work[k] = 0.0;
      }
      /* build lists of active rows and columns of matrix V; determine
       * number of non-zeros in initial active submatrix */
      nnz = 0;
      for (k = k2; k <= n; k++)
      {  i = pp_inv[k];
         sgf_activate_row(i);
         nnz += vr_len[i];
         j = qq_ind[k];
         sgf_activate_col(j);
      }
      /* main factorization loop */
      for (k = k2; k <= n; k++)
      {  int na;
         double den;
         /* calculate density of active submatrix */
         na = n - k + 1; /* order of active submatrix */
#if 0 /* 21/VIII-2014 */
         den = (double)nnz / (double)(na * na);
#else
         den = (double)nnz / ((double)(na) * (double)(na));
#endif
         /* if active submatrix is relatively dense, switch to dense
          * phase */
#if 1 /* FIXME */
         if (na >= 5 && den >= 0.71)
         {
#ifdef GLP_DEBUG
            xprintf("na = %d; nnz = %d; den = %g\n", na, nnz, den);
#endif
            break;
         }
#endif
         /* choose pivot v[p,q] */
         if (sgf_choose_pivot(sgf, &p, &q) != 0)
            return k; /* failure */
         /* u[i,j] = v[p,q], k <= i, j <= n */
         i = pp_ind[p];
         xassert(k <= i && i <= n);
         j = qq_inv[q];
         xassert(k <= j && j <= n);
         /* move u[i,j] to position u[k,k] by implicit permutations of
          * rows and columns of matrix U */
         luf_swap_u_rows(k, i);
         luf_swap_u_cols(k, j);
         /* perform gaussian elimination */
         nnz += sgf_eliminate(sgf, p, q);
      }
#if 1 /* FIXME */
      if (k <= n)
      {  /* continue computing factorization in dense mode */
#ifdef GLP_DEBUG
         sva_check_area(sva);
         luf_check_all(luf, k);
#endif
         k = sgf_dense_phase(luf, k, sgf->updat);
         if (k != 0)
            return k; /* failure */
      }
#endif
#ifdef GLP_DEBUG
      sva_check_area(sva);
      luf_check_all(luf, n+1);
#endif
      /* defragment SVA; currently all columns of V are empty, so they
       * will have zero capacity as required by luf_build_v_cols */
      sva_defrag_area(sva);
      /* build matrix F in row-wise format */
      luf_build_f_rows(luf, rs_head);
      /* build matrix V in column-wise format */
      luf_build_v_cols(luf, sgf->updat, rs_head);
      return 0;
}

/* eof */
