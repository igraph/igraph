/* btf.c (sparse block triangular LU-factorization) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2013-2014 Free Software Foundation, Inc.
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

#include "btf.h"
#include "env.h"
#include "luf.h"
#include "mc13d.h"
#include "mc21a.h"

/***********************************************************************
*  btf_store_a_cols - store pattern of matrix A in column-wise format
*
*  This routine stores the pattern (that is, only indices of non-zero
*  elements) of the original matrix A in column-wise format.
*
*  On exit the routine returns the number of non-zeros in matrix A. */

int btf_store_a_cols(BTF *btf, int (*col)(void *info, int j, int ind[],
      double val[]), void *info, int ind[], double val[])
{     int n = btf->n;
      SVA *sva = btf->sva;
      int *sv_ind = sva->ind;
      int ac_ref = btf->ac_ref;
      int *ac_ptr = &sva->ptr[ac_ref-1];
      int *ac_len = &sva->len[ac_ref-1];
      int j, len, ptr, nnz;
      nnz = 0;
      for (j = 1; j <= n; j++)
      {  /* get j-th column */
         len = col(info, j, ind, val);
         xassert(0 <= len && len <= n);
         /* reserve locations for j-th column */
         if (len > 0)
         {  if (sva->r_ptr - sva->m_ptr < len)
            {  sva_more_space(sva, len);
               sv_ind = sva->ind;
            }
            sva_reserve_cap(sva, ac_ref+(j-1), len);
         }
         /* store pattern of j-th column */
         ptr = ac_ptr[j];
         memcpy(&sv_ind[ptr], &ind[1], len * sizeof(int));
         ac_len[j] = len;
         nnz += len;
      }
      return nnz;
}

/***********************************************************************
*  btf_make_blocks - permutations to block triangular form
*
*  This routine analyzes the pattern of the original matrix A and
*  determines permutation matrices P and Q such that A = P * A~* Q,
*  where A~ is an upper block triangular matrix.
*
*  On exit the routine returns symbolic rank of matrix A. */

int btf_make_blocks(BTF *btf)
{     int n = btf->n;
      SVA *sva = btf->sva;
      int *sv_ind = sva->ind;
      int *pp_ind = btf->pp_ind;
      int *pp_inv = btf->pp_inv;
      int *qq_ind = btf->qq_ind;
      int *qq_inv = btf->qq_inv;
      int *beg = btf->beg;
      int ac_ref = btf->ac_ref;
      int *ac_ptr = &sva->ptr[ac_ref-1];
      int *ac_len = &sva->len[ac_ref-1];
      int i, j, rank, *iperm, *pr, *arp, *cv, *out, *ip, *lenr, *lowl,
         *numb, *prev;
      /* determine column permutation matrix M such that matrix A * M
       * has zero-free diagonal */
      iperm = qq_inv; /* matrix M */
      pr  = btf->p1_ind; /* working array */
      arp = btf->p1_inv; /* working array */
      cv  = btf->q1_ind; /* working array */
      out = btf->q1_inv; /* working array */
      rank = mc21a(n, sv_ind, ac_ptr, ac_len, iperm, pr, arp, cv, out);
      xassert(0 <= rank && rank <= n);
      if (rank < n)
      {  /* A is structurally singular (rank is its symbolic rank) */
         goto done;
      }
      /* build pattern of matrix A * M */
      ip   = pp_ind; /* working array */
      lenr = qq_ind; /* working array */
      for (j = 1; j <= n; j++)
      {  ip[j] = ac_ptr[iperm[j]];
         lenr[j] = ac_len[iperm[j]];
      }
      /* determine symmetric permutation matrix S such that matrix
       * S * (A * M) * S' = A~ is upper block triangular */
      lowl = btf->p1_ind; /* working array */
      numb = btf->p1_inv; /* working array */
      prev = btf->q1_ind; /* working array */
      btf->num =
         mc13d(n, sv_ind, ip, lenr, pp_inv, beg, lowl, numb, prev);
      xassert(beg[1] == 1);
      beg[btf->num+1] = n+1;
      /* A * M = S' * A~ * S ==> A = S' * A~ * (S * M') */
      /* determine permutation matrix P = S' */
      for (j = 1; j <= n; j++)
         pp_ind[pp_inv[j]] = j;
      /* determine permutation matrix Q = S * M' = P' * M' */
      for (i = 1; i <= n; i++)
         qq_ind[i] = iperm[pp_inv[i]];
      for (i = 1; i <= n; i++)
         qq_inv[qq_ind[i]] = i;
done: return rank;
}

/***********************************************************************
*  btf_check_blocks - check structure of matrix A~
*
*  This routine checks that structure of upper block triangular matrix
*  A~ is correct.
*
*  NOTE: For testing/debugging only. */

void btf_check_blocks(BTF *btf)
{     int n = btf->n;
      SVA *sva = btf->sva;
      int *sv_ind = sva->ind;
      int *pp_ind = btf->pp_ind;
      int *pp_inv = btf->pp_inv;
      int *qq_ind = btf->qq_ind;
      int *qq_inv = btf->qq_inv;
      int num = btf->num;
      int *beg = btf->beg;
      int ac_ref = btf->ac_ref;
      int *ac_ptr = &sva->ptr[ac_ref-1];
      int *ac_len = &sva->len[ac_ref-1];
      int i, ii, j, jj, k, size, ptr, end, diag;
      xassert(n > 0);
      /* check permutation matrices P and Q */
      for (k = 1; k <= n; k++)
      {  xassert(1 <= pp_ind[k] && pp_ind[k] <= n);
         xassert(pp_inv[pp_ind[k]] == k);
         xassert(1 <= qq_ind[k] && qq_ind[k] <= n);
         xassert(qq_inv[qq_ind[k]] == k);
      }
      /* check that matrix A~ is upper block triangular with non-zero
       * diagonal */
      xassert(1 <= num && num <= n);
      xassert(beg[1] == 1);
      xassert(beg[num+1] == n+1);
      /* walk thru blocks of A~ */
      for (k = 1; k <= num; k++)
      {  /* determine size of k-th block */
         size = beg[k+1] - beg[k];
         xassert(size >= 1);
         /* walk thru columns of k-th block */
         for (jj = beg[k]; jj < beg[k+1]; jj++)
         {  diag = 0;
            /* jj-th column of A~ = j-th column of A */
            j = qq_ind[jj];
            /* walk thru elements of j-th column of A */
            ptr = ac_ptr[j];
            end = ptr + ac_len[j];
            for (; ptr < end; ptr++)
            {  /* determine row index of a[i,j] */
               i = sv_ind[ptr];
               /* i-th row of A = ii-th row of A~ */
               ii = pp_ind[i];
               /* a~[ii,jj] should not be below k-th block */
               xassert(ii < beg[k+1]);
               if (ii == jj)
               {  /* non-zero diagonal element of A~ encountered */
                  diag = 1;
               }
            }
            xassert(diag);
         }
      }
      return;
}

/***********************************************************************
*  btf_build_a_rows - build matrix A in row-wise format
*
*  This routine builds the row-wise representation of matrix A in the
*  right part of SVA using its column-wise representation.
*
*  The working array len should have at least 1+n elements (len[0] is
*  not used). */

void btf_build_a_rows(BTF *btf, int len[/*1+n*/])
{     int n = btf->n;
      SVA *sva = btf->sva;
      int *sv_ind = sva->ind;
      double *sv_val = sva->val;
      int ar_ref = btf->ar_ref;
      int *ar_ptr = &sva->ptr[ar_ref-1];
      int *ar_len = &sva->len[ar_ref-1];
      int ac_ref = btf->ac_ref;
      int *ac_ptr = &sva->ptr[ac_ref-1];
      int *ac_len = &sva->len[ac_ref-1];
      int i, j, end, nnz, ptr, ptr1;
      /* calculate the number of non-zeros in each row of matrix A and
       * the total number of non-zeros */
      nnz = 0;
      for (i = 1; i <= n; i++)
         len[i] = 0;
      for (j = 1; j <= n; j++)
      {  nnz += ac_len[j];
         for (end = (ptr = ac_ptr[j]) + ac_len[j]; ptr < end; ptr++)
            len[sv_ind[ptr]]++;
      }
      /* we need at least nnz free locations in SVA */
      if (sva->r_ptr - sva->m_ptr < nnz)
      {  sva_more_space(sva, nnz);
         sv_ind = sva->ind;
         sv_val = sva->val;
      }
      /* reserve locations for rows of matrix A */
      for (i = 1; i <= n; i++)
      {  if (len[i] > 0)
            sva_reserve_cap(sva, ar_ref-1+i, len[i]);
         ar_len[i] = len[i];
      }
      /* walk thru columns of matrix A and build its rows */
      for (j = 1; j <= n; j++)
      {  for (end = (ptr = ac_ptr[j]) + ac_len[j]; ptr < end; ptr++)
         {  i = sv_ind[ptr];
            sv_ind[ptr1 = ar_ptr[i] + (--len[i])] = j;
            sv_val[ptr1] = sv_val[ptr];
         }
      }
      return;
}

/***********************************************************************
*  btf_a_solve - solve system A * x = b
*
*  This routine solves the system A * x = b, where A is the original
*  matrix.
*
*  On entry the array b should contain elements of the right-hand size
*  vector b in locations b[1], ..., b[n], where n is the order of the
*  matrix A. On exit the array x will contain elements of the solution
*  vector in locations x[1], ..., x[n]. Note that the array b will be
*  clobbered on exit.
*
*  The routine also uses locations [1], ..., [max_size] of two working
*  arrays w1 and w2, where max_size is the maximal size of diagonal
*  blocks in BT-factorization (max_size <= n). */

void btf_a_solve(BTF *btf, double b[/*1+n*/], double x[/*1+n*/],
      double w1[/*1+n*/], double w2[/*1+n*/])
{     SVA *sva = btf->sva;
      int *sv_ind = sva->ind;
      double *sv_val = sva->val;
      int *pp_inv = btf->pp_inv;
      int *qq_ind = btf->qq_ind;
      int num = btf->num;
      int *beg = btf->beg;
      int ac_ref = btf->ac_ref;
      int *ac_ptr = &sva->ptr[ac_ref-1];
      int *ac_len = &sva->len[ac_ref-1];
      double *bb = w1;
      double *xx = w2;
      LUF luf;
      int i, j, jj, k, beg_k, flag;
      double t;
      for (k = num; k >= 1; k--)
      {  /* determine order of diagonal block A~[k,k] */
         luf.n = beg[k+1] - (beg_k = beg[k]);
         if (luf.n == 1)
         {  /* trivial case */
            /* solve system A~[k,k] * X[k] = B[k] */
            t = x[qq_ind[beg_k]] =
               b[pp_inv[beg_k]] / btf->vr_piv[beg_k];
            /* substitute X[k] into other equations */
            if (t != 0.0)
            {  int ptr = ac_ptr[qq_ind[beg_k]];
               int end = ptr + ac_len[qq_ind[beg_k]];
               for (; ptr < end; ptr++)
                  b[sv_ind[ptr]] -= sv_val[ptr] * t;
            }
         }
         else
         {  /* general case */
            /* construct B[k] */
            flag = 0;
            for (i = 1; i <= luf.n; i++)
            {  if ((bb[i] = b[pp_inv[i + (beg_k-1)]]) != 0.0)
                  flag = 1;
            }
            /* solve system A~[k,k] * X[k] = B[k] */
            if (!flag)
            {  /* B[k] = 0, so X[k] = 0 */
               for (j = 1; j <= luf.n; j++)
                  x[qq_ind[j + (beg_k-1)]] = 0.0;
               continue;
            }
            luf.sva = sva;
            luf.fr_ref = btf->fr_ref + (beg_k-1);
            luf.fc_ref = btf->fc_ref + (beg_k-1);
            luf.vr_ref = btf->vr_ref + (beg_k-1);
            luf.vr_piv = btf->vr_piv + (beg_k-1);
            luf.vc_ref = btf->vc_ref + (beg_k-1);
            luf.pp_ind = btf->p1_ind + (beg_k-1);
            luf.pp_inv = btf->p1_inv + (beg_k-1);
            luf.qq_ind = btf->q1_ind + (beg_k-1);
            luf.qq_inv = btf->q1_inv + (beg_k-1);
            luf_f_solve(&luf, bb);
            luf_v_solve(&luf, bb, xx);
            /* store X[k] and substitute it into other equations */
            for (j = 1; j <= luf.n; j++)
            {  jj = j + (beg_k-1);
               t = x[qq_ind[jj]] = xx[j];
               if (t != 0.0)
               {  int ptr = ac_ptr[qq_ind[jj]];
                  int end = ptr + ac_len[qq_ind[jj]];
                  for (; ptr < end; ptr++)
                     b[sv_ind[ptr]] -= sv_val[ptr] * t;
               }
            }
         }
      }
      return;
}

/***********************************************************************
*  btf_at_solve - solve system A'* x = b
*
*  This routine solves the system A'* x = b, where A' is a matrix
*  transposed to the original matrix A.
*
*  On entry the array b should contain elements of the right-hand size
*  vector b in locations b[1], ..., b[n], where n is the order of the
*  matrix A. On exit the array x will contain elements of the solution
*  vector in locations x[1], ..., x[n]. Note that the array b will be
*  clobbered on exit.
*
*  The routine also uses locations [1], ..., [max_size] of two working
*  arrays w1 and w2, where max_size is the maximal size of diagonal
*  blocks in BT-factorization (max_size <= n). */

void btf_at_solve(BTF *btf, double b[/*1+n*/], double x[/*1+n*/],
      double w1[/*1+n*/], double w2[/*1+n*/])
{     SVA *sva = btf->sva;
      int *sv_ind = sva->ind;
      double *sv_val = sva->val;
      int *pp_inv = btf->pp_inv;
      int *qq_ind = btf->qq_ind;
      int num = btf->num;
      int *beg = btf->beg;
      int ar_ref = btf->ar_ref;
      int *ar_ptr = &sva->ptr[ar_ref-1];
      int *ar_len = &sva->len[ar_ref-1];
      double *bb = w1;
      double *xx = w2;
      LUF luf;
      int i, j, jj, k, beg_k, flag;
      double t;
      for (k = 1; k <= num; k++)
      {  /* determine order of diagonal block A~[k,k] */
         luf.n = beg[k+1] - (beg_k = beg[k]);
         if (luf.n == 1)
         {  /* trivial case */
            /* solve system A~'[k,k] * X[k] = B[k] */
            t = x[pp_inv[beg_k]] =
               b[qq_ind[beg_k]] / btf->vr_piv[beg_k];
            /* substitute X[k] into other equations */
            if (t != 0.0)
            {  int ptr = ar_ptr[pp_inv[beg_k]];
               int end = ptr + ar_len[pp_inv[beg_k]];
               for (; ptr < end; ptr++)
                  b[sv_ind[ptr]] -= sv_val[ptr] * t;
            }
         }
         else
         {  /* general case */
            /* construct B[k] */
            flag = 0;
            for (i = 1; i <= luf.n; i++)
            {  if ((bb[i] = b[qq_ind[i + (beg_k-1)]]) != 0.0)
                  flag = 1;
            }
            /* solve system A~'[k,k] * X[k] = B[k] */
            if (!flag)
            {  /* B[k] = 0, so X[k] = 0 */
               for (j = 1; j <= luf.n; j++)
                  x[pp_inv[j + (beg_k-1)]] = 0.0;
               continue;
            }
            luf.sva = sva;
            luf.fr_ref = btf->fr_ref + (beg_k-1);
            luf.fc_ref = btf->fc_ref + (beg_k-1);
            luf.vr_ref = btf->vr_ref + (beg_k-1);
            luf.vr_piv = btf->vr_piv + (beg_k-1);
            luf.vc_ref = btf->vc_ref + (beg_k-1);
            luf.pp_ind = btf->p1_ind + (beg_k-1);
            luf.pp_inv = btf->p1_inv + (beg_k-1);
            luf.qq_ind = btf->q1_ind + (beg_k-1);
            luf.qq_inv = btf->q1_inv + (beg_k-1);
            luf_vt_solve(&luf, bb, xx);
            luf_ft_solve(&luf, xx);
            /* store X[k] and substitute it into other equations */
            for (j = 1; j <= luf.n; j++)
            {  jj = j + (beg_k-1);
               t = x[pp_inv[jj]] = xx[j];
               if (t != 0.0)
               {  int ptr = ar_ptr[pp_inv[jj]];
                  int end = ptr + ar_len[pp_inv[jj]];
                  for (; ptr < end; ptr++)
                     b[sv_ind[ptr]] -= sv_val[ptr] * t;
               }
            }
         }
      }
      return;
}

/***********************************************************************
*  btf_at_solve1 - solve system A'* y = e' to cause growth in y
*
*  This routine is a special version of btf_at_solve. It solves the
*  system A'* y = e' = e + delta e, where A' is a matrix transposed to
*  the original matrix A, e is the specified right-hand side vector,
*  and delta e is a vector of +1 and -1 chosen to cause growth in the
*  solution vector y.
*
*  On entry the array e should contain elements of the right-hand size
*  vector e in locations e[1], ..., e[n], where n is the order of the
*  matrix A. On exit the array y will contain elements of the solution
*  vector in locations y[1], ..., y[n]. Note that the array e will be
*  clobbered on exit.
*
*  The routine also uses locations [1], ..., [max_size] of two working
*  arrays w1 and w2, where max_size is the maximal size of diagonal
*  blocks in BT-factorization (max_size <= n). */

void btf_at_solve1(BTF *btf, double e[/*1+n*/], double y[/*1+n*/],
      double w1[/*1+n*/], double w2[/*1+n*/])
{     SVA *sva = btf->sva;
      int *sv_ind = sva->ind;
      double *sv_val = sva->val;
      int *pp_inv = btf->pp_inv;
      int *qq_ind = btf->qq_ind;
      int num = btf->num;
      int *beg = btf->beg;
      int ar_ref = btf->ar_ref;
      int *ar_ptr = &sva->ptr[ar_ref-1];
      int *ar_len = &sva->len[ar_ref-1];
      double *ee = w1;
      double *yy = w2;
      LUF luf;
      int i, j, jj, k, beg_k, ptr, end;
      double e_k, y_k;
      for (k = 1; k <= num; k++)
      {  /* determine order of diagonal block A~[k,k] */
         luf.n = beg[k+1] - (beg_k = beg[k]);
         if (luf.n == 1)
         {  /* trivial case */
            /* determine E'[k] = E[k] + delta E[k] */
            e_k = e[qq_ind[beg_k]];
            e_k = (e_k >= 0.0 ? e_k + 1.0 : e_k - 1.0);
            /* solve system A~'[k,k] * Y[k] = E[k] */
            y_k = y[pp_inv[beg_k]] = e_k / btf->vr_piv[beg_k];
            /* substitute Y[k] into other equations */
            ptr = ar_ptr[pp_inv[beg_k]];
            end = ptr + ar_len[pp_inv[beg_k]];
            for (; ptr < end; ptr++)
               e[sv_ind[ptr]] -= sv_val[ptr] * y_k;
         }
         else
         {  /* general case */
            /* construct E[k] */
            for (i = 1; i <= luf.n; i++)
               ee[i] = e[qq_ind[i + (beg_k-1)]];
            /* solve system A~'[k,k] * Y[k] = E[k] + delta E[k] */
            luf.sva = sva;
            luf.fr_ref = btf->fr_ref + (beg_k-1);
            luf.fc_ref = btf->fc_ref + (beg_k-1);
            luf.vr_ref = btf->vr_ref + (beg_k-1);
            luf.vr_piv = btf->vr_piv + (beg_k-1);
            luf.vc_ref = btf->vc_ref + (beg_k-1);
            luf.pp_ind = btf->p1_ind + (beg_k-1);
            luf.pp_inv = btf->p1_inv + (beg_k-1);
            luf.qq_ind = btf->q1_ind + (beg_k-1);
            luf.qq_inv = btf->q1_inv + (beg_k-1);
            luf_vt_solve1(&luf, ee, yy);
            luf_ft_solve(&luf, yy);
            /* store Y[k] and substitute it into other equations */
            for (j = 1; j <= luf.n; j++)
            {  jj = j + (beg_k-1);
               y_k = y[pp_inv[jj]] = yy[j];
               ptr = ar_ptr[pp_inv[jj]];
               end = ptr + ar_len[pp_inv[jj]];
               for (; ptr < end; ptr++)
                  e[sv_ind[ptr]] -= sv_val[ptr] * y_k;
            }
         }
      }
      return;
}

/***********************************************************************
*  btf_estimate_norm - estimate 1-norm of inv(A)
*
*  This routine estimates 1-norm of inv(A) by one step of inverse
*  iteration for the small singular vector as described in [1]. This
*  involves solving two systems of equations:
*
*     A'* y = e,
*
*     A * z = y,
*
*  where A' is a matrix transposed to A, and e is a vector of +1 and -1
*  chosen to cause growth in y. Then
*
*     estimate 1-norm of inv(A) = (1-norm of z) / (1-norm of y)
*
*  REFERENCES
*
*  1. G.E.Forsythe, M.A.Malcolm, C.B.Moler. Computer Methods for
*     Mathematical Computations. Prentice-Hall, Englewood Cliffs, N.J.,
*     pp. 30-62 (subroutines DECOMP and SOLVE). */

double btf_estimate_norm(BTF *btf, double w1[/*1+n*/], double
      w2[/*1+n*/], double w3[/*1+n*/], double w4[/*1+n*/])
{     int n = btf->n;
      double *e = w1;
      double *y = w2;
      double *z = w1;
      int i;
      double y_norm, z_norm;
      /* compute y = inv(A') * e to cause growth in y */
      for (i = 1; i <= n; i++)
         e[i] = 0.0;
      btf_at_solve1(btf, e, y, w3, w4);
      /* compute 1-norm of y = sum |y[i]| */
      y_norm = 0.0;
      for (i = 1; i <= n; i++)
         y_norm += (y[i] >= 0.0 ? +y[i] : -y[i]);
      /* compute z = inv(A) * y */
      btf_a_solve(btf, y, z, w3, w4);
      /* compute 1-norm of z = sum |z[i]| */
      z_norm = 0.0;
      for (i = 1; i <= n; i++)
         z_norm += (z[i] >= 0.0 ? +z[i] : -z[i]);
      /* estimate 1-norm of inv(A) = (1-norm of z) / (1-norm of y) */
      return z_norm / y_norm;
}

/* eof */
