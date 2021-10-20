/* spxat.c */

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
#include "spxat.h"

/***********************************************************************
*  spx_alloc_at - allocate constraint matrix in sparse row-wise format
*
*  This routine allocates the memory for arrays needed to represent the
*  constraint matrix in sparse row-wise format. */

void spx_alloc_at(SPXLP *lp, SPXAT *at)
{     int m = lp->m;
      int n = lp->n;
      int nnz = lp->nnz;
      at->ptr = talloc(1+m+1, int);
      at->ind = talloc(1+nnz, int);
      at->val = talloc(1+nnz, double);
      at->work = talloc(1+n, double);
      return;
}

/***********************************************************************
*  spx_build_at - build constraint matrix in sparse row-wise format
*
*  This routine builds sparse row-wise representation of the constraint
*  matrix A using its sparse column-wise representation stored in the
*  lp object, and stores the result in the at object. */

void spx_build_at(SPXLP *lp, SPXAT *at)
{     int m = lp->m;
      int n = lp->n;
      int nnz = lp->nnz;
      int *A_ptr = lp->A_ptr;
      int *A_ind = lp->A_ind;
      double *A_val = lp->A_val;
      int *AT_ptr = at->ptr;
      int *AT_ind = at->ind;
      double *AT_val = at->val;
      int i, k, ptr, end, pos;
      /* calculate AT_ptr[i] = number of non-zeros in i-th row */
      memset(&AT_ptr[1], 0, m * sizeof(int));
      for (k = 1; k <= n; k++)
      {  ptr = A_ptr[k];
         end = A_ptr[k+1];
         for (; ptr < end; ptr++)
            AT_ptr[A_ind[ptr]]++;
      }
      /* set AT_ptr[i] to position after last element in i-th row */
      AT_ptr[1]++;
      for (i = 2; i <= m; i++)
         AT_ptr[i] += AT_ptr[i-1];
      xassert(AT_ptr[m] == nnz+1);
      AT_ptr[m+1] = nnz+1;
      /* build row-wise representation and re-arrange AT_ptr[i] */
      for (k = n; k >= 1; k--)
      {  /* copy elements from k-th column to corresponding rows */
         ptr = A_ptr[k];
         end = A_ptr[k+1];
         for (; ptr < end; ptr++)
         {  pos = --AT_ptr[A_ind[ptr]];
            AT_ind[pos] = k;
            AT_val[pos] = A_val[ptr];
         }
      }
      xassert(AT_ptr[1] == 1);
      return;
}

/***********************************************************************
*  spx_at_prod - compute product y := y + s * A'* x
*
*  This routine computes the product:
*
*     y := y + s * A'* x,
*
*  where A' is a matrix transposed to the mxn-matrix A of constraint
*  coefficients, x is a m-vector, s is a scalar, y is a n-vector.
*
*  The routine uses the row-wise representation of the matrix A and
*  computes the product as a linear combination:
*
*     y := y + s * (A'[1] * x[1] + ... + A'[m] * x[m]),
*
*  where A'[i] is i-th row of A, 1 <= i <= m. */

void spx_at_prod(SPXLP *lp, SPXAT *at, double y[/*1+n*/], double s,
      const double x[/*1+m*/])
{     int m = lp->m;
      int *AT_ptr = at->ptr;
      int *AT_ind = at->ind;
      double *AT_val = at->val;
      int i, ptr, end;
      double t;
      for (i = 1; i <= m; i++)
      {  if (x[i] != 0.0)
         {  /* y := y + s * (i-th row of A) * x[i] */
            t = s * x[i];
            ptr = AT_ptr[i];
            end = AT_ptr[i+1];
            for (; ptr < end; ptr++)
               y[AT_ind[ptr]] += AT_val[ptr] * t;
         }
      }
      return;
}

/***********************************************************************
*  spx_nt_prod1 - compute product y := y + s * N'* x
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
*  of the array y assuming that y = 0. */

void spx_nt_prod1(SPXLP *lp, SPXAT *at, double y[/*1+n-m*/], int ign,
      double s, const double x[/*1+m*/])
{     int m = lp->m;
      int n = lp->n;
      int *head = lp->head;
      double *work = at->work;
      int j, k;
      for (k = 1; k <= n; k++)
         work[k] = 0.0;
      if (!ign)
      {  for (j = 1; j <= n-m; j++)
            work[head[m+j]] = y[j];
      }
      spx_at_prod(lp, at, work, s, x);
      for (j = 1; j <= n-m; j++)
         y[j] = work[head[m+j]];
      return;
}

/***********************************************************************
*  spx_eval_trow1 - compute i-th row of simplex table
*
*  This routine computes i-th row of the current simplex table
*  T = (T[i,j]) = - inv(B) * N, 1 <= i <= m, using representation of
*  the constraint matrix A in row-wise format.
*
*  The vector rho = (rho[j]), which is i-th row of the basis inverse
*  inv(B), should be previously computed with the routine spx_eval_rho.
*  It is assumed that elements of this vector are stored in the array
*  locations rho[1], ..., rho[m].
*
*  There exist two ways to compute the simplex table row.
*
*  1. T[i,j], j = 1,...,n-m, is computed as inner product:
*
*                    m
*        T[i,j] = - sum a[i,k] * rho[i],
*                   i=1
*
*  where N[j] = A[k] is a column of the constraint matrix corresponding
*  to non-basic variable xN[j]. The estimated number of operations in
*  this case is:
*
*        n1 = (n - m) * (nnz(A) / n),
*
*  (n - m) is the number of columns of N, nnz(A) / n is the average
*  number of non-zeros in one column of A and, therefore, of N.
*
*  2. The simplex table row is computed as part of a linear combination
*     of rows of A with coefficients rho[i] != 0. The estimated number
*     of operations in this case is:
*
*        n2 = nnz(rho) * (nnz(A) / m),
*
*     where nnz(rho) is the number of non-zeros in the vector rho,
*     nnz(A) / m is the average number of non-zeros in one row of A.
*
*  If n1 < n2, the routine computes the simples table row using the
*  first way (like the routine spx_eval_trow). Otherwise, the routine
*  uses the second way calling the routine spx_nt_prod1.
*
*  On exit components of the simplex table row are stored in the array
*  locations trow[1], ... trow[n-m]. */

void spx_eval_trow1(SPXLP *lp, SPXAT *at, const double rho[/*1+m*/],
      double trow[/*1+n-m*/])
{     int m = lp->m;
      int n = lp->n;
      int nnz = lp->nnz;
      int i, j, nnz_rho;
      double cnt1, cnt2;
      /* determine nnz(rho) */
      nnz_rho = 0;
      for (i = 1; i <= m; i++)
      {  if (rho[i] != 0.0)
            nnz_rho++;
      }
      /* estimate the number of operations for both ways */
      cnt1 = (double)(n - m) * ((double)nnz / (double)n);
      cnt2 = (double)nnz_rho * ((double)nnz / (double)m);
      /* compute i-th row of simplex table */
      if (cnt1 < cnt2)
      {  /* as inner products */
         int *A_ptr = lp->A_ptr;
         int *A_ind = lp->A_ind;
         double *A_val = lp->A_val;
         int *head = lp->head;
         int k, ptr, end;
         double tij;
         for (j = 1; j <= n-m; j++)
         {  k = head[m+j]; /* x[k] = xN[j] */
            /* compute t[i,j] = - N'[j] * pi */
            tij = 0.0;
            ptr = A_ptr[k];
            end = A_ptr[k+1];
            for (; ptr < end; ptr++)
               tij -= A_val[ptr] * rho[A_ind[ptr]];
            trow[j] = tij;
         }
      }
      else
      {  /* as linear combination */
         spx_nt_prod1(lp, at, trow, 1, -1.0, rho);
      }
      return;
}

/***********************************************************************
*  spx_free_at - deallocate constraint matrix in sparse row-wise format
*
*  This routine deallocates the memory used for arrays of the program
*  object at. */

void spx_free_at(SPXLP *lp, SPXAT *at)
{     xassert(lp == lp);
      tfree(at->ptr);
      tfree(at->ind);
      tfree(at->val);
      tfree(at->work);
      return;
}

/* eof */
