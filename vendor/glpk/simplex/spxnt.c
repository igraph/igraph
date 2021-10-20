/* spxnt.c */

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

#include "env.h"
#include "spxnt.h"

/***********************************************************************
*  spx_alloc_nt - allocate matrix N in sparse row-wise format
*
*  This routine allocates the memory for arrays needed to represent the
*  matrix N composed of non-basic columns of the constraint matrix A. */

void spx_alloc_nt(SPXLP *lp, SPXNT *nt)
{     int m = lp->m;
      int nnz = lp->nnz;
      nt->ptr = talloc(1+m, int);
      nt->len = talloc(1+m, int);
      nt->ind = talloc(1+nnz, int);
      nt->val = talloc(1+nnz, double);
      return;
}

/***********************************************************************
*  spx_init_nt - initialize row pointers for matrix N
*
*  This routine initializes (sets up) row pointers for the matrix N
*  using column-wise representation of the constraint matrix A.
*
*  This routine needs to be called only once. */

void spx_init_nt(SPXLP *lp, SPXNT *nt)
{     int m = lp->m;
      int n = lp->n;
      int nnz = lp->nnz;
      int *A_ptr = lp->A_ptr;
      int *A_ind = lp->A_ind;
      int *NT_ptr = nt->ptr;
      int *NT_len = nt->len;
      int i, k, ptr, end;
      /* calculate NT_len[i] = maximal number of non-zeros in i-th row
       * of N = number of non-zeros in i-th row of A */
      memset(&NT_len[1], 0, m * sizeof(int));
      for (k = 1; k <= n; k++)
      {  ptr = A_ptr[k];
         end = A_ptr[k+1];
         for (; ptr < end; ptr++)
            NT_len[A_ind[ptr]]++;
      }
      /* initialize row pointers NT_ptr[i], i = 1,...,n-m */
      NT_ptr[1] = 1;
      for (i = 2; i <= m; i++)
         NT_ptr[i] = NT_ptr[i-1] + NT_len[i-1];
      xassert(NT_ptr[m] + NT_len[m] == nnz+1);
      return;
}

/***********************************************************************
*  spx_nt_add_col - add column N[j] = A[k] to matrix N
*
*  This routine adds elements of column N[j] = A[k], 1 <= j <= n-m,
*  1 <= k <= n, to the row-wise represntation of the matrix N. It is
*  assumed (with no check) that elements of the specified column are
*  missing in the row-wise represntation of N. */

void spx_nt_add_col(SPXLP *lp, SPXNT *nt, int j, int k)
{     int m = lp->m;
      int n = lp->n;
      int nnz = lp->nnz;
      int *A_ptr = lp->A_ptr;
      int *A_ind = lp->A_ind;
      double *A_val = lp->A_val;
      int *NT_ptr = nt->ptr;
      int *NT_len = nt->len;
      int *NT_ind = nt->ind;
      double *NT_val = nt->val;
      int i, ptr, end, pos;
      xassert(1 <= j && j <= n-m);
      xassert(1 <= k && k <= n);
      ptr = A_ptr[k];
      end = A_ptr[k+1];
      for (; ptr < end; ptr++)
      {  i = A_ind[ptr];
         /* add element N[i,j] = A[i,k] to i-th row of matrix N */
         pos = NT_ptr[i] + (NT_len[i]++);
         if (i < m)
            xassert(pos < NT_ptr[i+1]);
         else
            xassert(pos <= nnz);
         NT_ind[pos] = j;
         NT_val[pos] = A_val[ptr];
      }
      return;
}

/***********************************************************************
*  spx_build_nt - build matrix N for current basis
*
*  This routine builds the row-wise represntation of the matrix N
*  for the current basis by adding columns of the constraint matrix A
*  corresponding to non-basic variables. */

void spx_build_nt(SPXLP *lp, SPXNT *nt)
{     int m = lp->m;
      int n = lp->n;
      int *head = lp->head;
      int *NT_len = nt->len;
      int j, k;
      /* N := 0 */
      memset(&NT_len[1], 0, m * sizeof(int));
      /* add non-basic columns N[j] = A[k] */
      for (j = 1; j <= n-m; j++)
      {  k = head[m+j]; /* x[k] = xN[j] */
         spx_nt_add_col(lp, nt, j, k);
      }
      return;
}

/***********************************************************************
*  spx_nt_del_col - remove column N[j] = A[k] from matrix N
*
*  This routine removes elements of column N[j] = A[k], 1 <= j <= n-m,
*  1 <= k <= n, from the row-wise representation of the matrix N. It is
*  assumed (with no check) that elements of the specified column are
*  present in the row-wise representation of N. */

void spx_nt_del_col(SPXLP *lp, SPXNT *nt, int j, int k)
{     int m = lp->m;
      int n = lp->n;
      int *A_ptr = lp->A_ptr;
      int *A_ind = lp->A_ind;
      int *NT_ptr = nt->ptr;
      int *NT_len = nt->len;
      int *NT_ind = nt->ind;
      double *NT_val = nt->val;
      int i, ptr, end, ptr1, end1;
      xassert(1 <= j && j <= n-m);
      xassert(1 <= k && k <= n);
      ptr = A_ptr[k];
      end = A_ptr[k+1];
      for (; ptr < end; ptr++)
      {  i = A_ind[ptr];
         /* find element N[i,j] = A[i,k] in i-th row of matrix N */
         ptr1 = NT_ptr[i];
         end1 = ptr1 + NT_len[i];
         for (; NT_ind[ptr1] != j; ptr1++)
            /* nop */;
         xassert(ptr1 < end1);
         /* and remove it from i-th row element list */
         NT_len[i]--;
         NT_ind[ptr1] = NT_ind[end1-1];
         NT_val[ptr1] = NT_val[end1-1];
      }
      return;
}

/***********************************************************************
*  spx_update_nt - update matrix N for adjacent basis
*
*  This routine updates the row-wise represntation of matrix N for
*  the adjacent basis, where column N[q], 1 <= q <= n-m, is replaced by
*  column B[p], 1 <= p <= m, of the current basis matrix B. */

void spx_update_nt(SPXLP *lp, SPXNT *nt, int p, int q)
{     int m = lp->m;
      int n = lp->n;
      int *head = lp->head;
      xassert(1 <= p && p <= m);
      xassert(1 <= q && q <= n-m);
      /* remove old column N[q] corresponding to variable xN[q] */
      spx_nt_del_col(lp, nt, q, head[m+q]);
      /* add new column N[q] corresponding to variable xB[p] */
      spx_nt_add_col(lp, nt, q, head[p]);
      return;
}

/***********************************************************************
*  spx_nt_prod - compute product y := y + s * N'* x
*
*  This routine computes the product:
*
*     y := y + s * N'* x,
*
*  where N' is a matrix transposed to the mx(n-m)-matrix N composed
*  from non-basic columns of the constraint matrix A, x is a m-vector,
*  s is a scalar, y is (n-m)-vector.
*
*  If the flag ign is non-zero, the routine ignores the input content
*  of the array y assuming that y = 0.
*
*  The routine uses the row-wise representation of the matrix N and
*  computes the product as a linear combination:
*
*     y := y + s * (N'[1] * x[1] + ... + N'[m] * x[m]),
*
*  where N'[i] is i-th row of N, 1 <= i <= m. */

void spx_nt_prod(SPXLP *lp, SPXNT *nt, double y[/*1+n-m*/], int ign,
      double s, const double x[/*1+m*/])
{     int m = lp->m;
      int n = lp->n;
      int *NT_ptr = nt->ptr;
      int *NT_len = nt->len;
      int *NT_ind = nt->ind;
      double *NT_val = nt->val;
      int i, j, ptr, end;
      double t;
      if (ign)
      {  /* y := 0 */
         for (j = 1; j <= n-m; j++)
            y[j] = 0.0;
      }
      for (i = 1; i <= m; i++)
      {  if (x[i] != 0.0)
         {  /* y := y + s * (i-th row of N) * x[i] */
            t = s * x[i];
            ptr = NT_ptr[i];
            end = ptr + NT_len[i];
            for (; ptr < end; ptr++)
               y[NT_ind[ptr]] += NT_val[ptr] * t;
         }
      }
      return;
}

#if 1 /* 31/III-2016 */
void spx_nt_prod_s(SPXLP *lp, SPXNT *nt, FVS *y, int ign, double s,
      const FVS *x, double eps)
{     /* sparse version of spx_nt_prod */
      int *NT_ptr = nt->ptr;
      int *NT_len = nt->len;
      int *NT_ind = nt->ind;
      double *NT_val = nt->val;
      int *x_ind = x->ind;
      double *x_vec = x->vec;
      int *y_ind = y->ind;
      double *y_vec = y->vec;
      int i, j, k, nnz, ptr, end;
      double t;
      xassert(x->n == lp->m);
      xassert(y->n == lp->n-lp->m);
      if (ign)
      {  /* y := 0 */
         fvs_clear_vec(y);
      }
      nnz = y->nnz;
      for (k = x->nnz; k >= 1; k--)
      {  i = x_ind[k];
         /* y := y + s * (i-th row of N) * x[i] */
         t = s * x_vec[i];
         ptr = NT_ptr[i];
         end = ptr + NT_len[i];
         for (; ptr < end; ptr++)
         {  j = NT_ind[ptr];
            if (y_vec[j] == 0.0)
               y_ind[++nnz] = j;
            y_vec[j] += NT_val[ptr] * t;
            /* don't forget about numeric cancellation */
            if (y_vec[j] == 0.0)
               y_vec[j] = DBL_MIN;
         }
      }
      y->nnz = nnz;
      fvs_adjust_vec(y, eps);
      return;
}
#endif

/***********************************************************************
*  spx_free_nt - deallocate matrix N in sparse row-wise format
*
*  This routine deallocates the memory used for arrays of the program
*  object nt. */

void spx_free_nt(SPXLP *lp, SPXNT *nt)
{     xassert(lp == lp);
      tfree(nt->ptr);
      tfree(nt->len);
      tfree(nt->ind);
      tfree(nt->val);
      return;
}

/* eof */
