/* triang.c (find maximal triangular part of rectangular matrix) */

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
#include "triang.h"

/***********************************************************************
*  triang - find maximal triangular part of rectangular matrix
*
*  Given a mxn sparse matrix A this routine finds permutation matrices
*  P and Q such that matrix A' = P * A * Q has the following structure:
*
*        1         s         n
*     1  * . . . . . x x x x x
*        * * . . . . x x x x x
*        * * * . . . x x x x x
*        * * * * . . x x x x x
*        * * * * * . x x x x x
*     s  * * * * * * x x x x x
*        x x x x x x x x x x x
*        x x x x x x x x x x x
*     m  x x x x x x x x x x x
*
*  where '*' are elements of the triangular part, '.' are structural
*  zeros, 'x' are other elements.
*
*  The formal routine mat specifies the original matrix A in both row-
*  and column-wise format. If the routine mat is called with k = +i,
*  1 <= i <= m, it should store column indices and values of non-zero
*  elements of i-th row of A in locations ind[1], ..., ind[len] and
*  val[1], ..., val[len], resp., where len is the returned number of
*  non-zeros in the row, 0 <= len <= n. Similarly, if the routine mat
*  is called with k = -j, 1 <= j <= n, it should store row indices and
*  values of non-zero elements of j-th column of A and return len, the
*  number of non-zeros in the column, 0 <= len <= m. Should note that
*  duplicate indices are not allowed.
*
*  The parameter info is a transit pointer passed to the routine mat.
*
*  The parameter tol is a tolerance. The routine triang guarantees that
*  each diagonal element in the triangular part of matrix A' is not
*  less in magnitude than tol * max, where max is the maximal magnitude
*  of elements in corresponding column.
*
*  On exit the routine triang stores information on the triangular part
*  found in the arrays rn and cn. Elements rn[1], ..., rn[s] specify
*  row numbers and elements cn[1], ..., cn[s] specify column numbers
*  of the original matrix A, which correspond to rows/columns 1, ..., s
*  of matrix A', where s is the size of the triangular part returned by
*  the routine, 0 <= s <= min(m, n). The order of rows and columns that
*  are not included in the triangular part remains unspecified.
*
*  ALGORITHM
*
*  The routine triang uses a simple greedy heuristic.
*
*  At some step the matrix A' = P * A * Q has the following structure:
*
*        1                   n
*     1  * . . . . . . . x x x
*        * * . . . . . . x x x
*        * * * . . . . . x x x
*        * * * * . . . . x x x
*        x x x x # # # # x x x
*        x x x x # # # # x x x
*        x x x x # # # # x x x
*        x x x x # # # # x x x
*     m  x x x x # # # # x x x
*
*  where '#' are elements of active submatrix. Initially P = Q = I, so
*  the active submatrix is the original matrix A = A'.
*
*  If some row has exactly one non-zero in the active submatrix (row
*  singleton), the routine includes this row and corresponding column
*  in the triangular part, and removes the column from the active
*  submatrix. Otherwise, the routine simply removes a column having
*  maximal number of non-zeros from the active submatrix in the hope
*  that new row singleton(s) will appear.
*
*  COMPLEXITY
*
*  The time complexity of the routine triang is O(nnz), where nnz is
*  number of non-zeros in the original matrix A. */

int triang(int m, int n, int (*mat)(void *info, int k, int ind[],
      double val[]), void *info, double tol, int rn[], int cn[])
{     int head, i, j, jj, k, kk, ks, len, len2, next_j, ns, size;
      int *cind, *rind, *cnt, *ptr, *list, *prev, *next;
      double *cval, *rval, *big;
      char *flag;
      /* allocate working arrays */
      cind = talloc(1+m, int);
      cval = talloc(1+m, double);
      rind = talloc(1+n, int);
      rval = talloc(1+n, double);
      cnt = ptr = talloc(1+m, int);
      list = talloc(1+n, int);
      prev = talloc(1+n, int);
      next = talloc(1+n, int);
      big = talloc(1+n, double);
      flag = talloc(1+n, char);
      /*--------------------------------------------------------------*/
      /* build linked lists of columns having equal lengths           */
      /*--------------------------------------------------------------*/
      /* ptr[len], 0 <= len <= m, is number of first column of length
       * len;
       * next[j], 1 <= j <= n, is number of next column having the same
       * length as column j;
       * big[j], 1 <= j <= n, is maximal magnitude of elements in j-th
       * column */
      for (len = 0; len <= m; len++)
         ptr[len] = 0;
      for (j = 1; j <= n; j++)
      {  /* get j-th column */
         len = mat(info, -j, cind, cval);
         xassert(0 <= len && len <= m);
         /* add this column to beginning of list ptr[len] */
         next[j] = ptr[len];
         ptr[len] = j;
         /* determine maximal magnitude of elements in this column */
         big[j] = 0.0;
         for (k = 1; k <= len; k++)
         {  if (big[j] < fabs(cval[k]))
               big[j] = fabs(cval[k]);
         }
      }
      /*--------------------------------------------------------------*/
      /* build doubly linked list of columns ordered by decreasing    */
      /* column lengths                                               */
      /*--------------------------------------------------------------*/
      /* head is number of first column in the list;
       * prev[j], 1 <= j <= n, is number of column that precedes j-th
       * column in the list;
       * next[j], 1 <= j <= n, is number of column that follows j-th
       * column in the list */
      head = 0;
      for (len = 0; len <= m; len++)
      {  /* walk thru list of columns of length len */
         for (j = ptr[len]; j != 0; j = next_j)
         {  next_j = next[j];
            /* add j-th column to beginning of the column list */
            prev[j] = 0;
            next[j] = head;
            if (head != 0)
               prev[head] = j;
            head = j;
         }
      }
      /*--------------------------------------------------------------*/
      /* build initial singleton list                                 */
      /*--------------------------------------------------------------*/
      /* there are used two list of columns:
       * 1) doubly linked list of active columns, in which all columns
       *    are ordered by decreasing column lengths;
       * 2) singleton list; an active column is included in this list
       *    if it has at least one row singleton in active submatrix */
      /* flag[j], 1 <= j <= n, is a flag of j-th column:
       * 0  j-th column is inactive;
       * 1  j-th column is active;
       * 2  j-th column is active and has row singleton(s) */
      /* initially all columns are active */
      for (j = 1; j <= n; j++)
         flag[j] = 1;
      /* initialize row counts and build initial singleton list */
      /* cnt[i], 1 <= i <= m, is number of non-zeros, which i-th row
       * has in active submatrix;
       * ns is size of singleton list;
       * list[1], ..., list[ns] are numbers of active columns included
       * in the singleton list */
      ns = 0;
      for (i = 1; i <= m; i++)
      {  /* get i-th row */
         len = cnt[i] = mat(info, +i, rind, rval);
         xassert(0 <= len && len <= n);
         if (len == 1)
         {  /* a[i,j] is row singleton */
            j = rind[1];
            xassert(1 <= j && j <= n);
            if (flag[j] != 2)
            {  /* include j-th column in singleton list */
               flag[j] = 2;
               list[++ns] = j;
            }
         }
      }
      /*--------------------------------------------------------------*/
      /* main loop                                                    */
      /*--------------------------------------------------------------*/
      size = 0; /* size of triangular part */
      /* loop until active column list is non-empty, i.e. until the
       * active submatrix has at least one column */
      while (head != 0)
      {  if (ns == 0)
         {  /* singleton list is empty */
            /* remove from the active submatrix a column of maximal
             * length in the hope that some row singletons appear */
            j = head;
            len = mat(info, -j, cind, cval);
            xassert(0 <= len && len <= m);
            goto drop;
         }
         /* take column j from the singleton list */
         j = list[ns--];
         xassert(flag[j] == 2);
         /* j-th column has at least one row singleton in the active
          * submatrix; choose one having maximal magnitude */
         len = mat(info, -j, cind, cval);
         xassert(0 <= len && len <= m);
         kk = 0;
         for (k = 1; k <= len; k++)
         {  i = cind[k];
            xassert(1 <= i && i <= m);
            if (cnt[i] == 1)
            {  /* a[i,j] is row singleton */
               if (kk == 0 || fabs(cval[kk]) < fabs(cval[k]))
                  kk = k;
            }
         }
         xassert(kk > 0);
         /* check magnitude of the row singleton chosen */
         if (fabs(cval[kk]) < tol * big[j])
         {  /* all row singletons are too small in magnitude; drop j-th
             * column */
            goto drop;
         }
         /* row singleton a[i,j] is ok; add i-th row and j-th column to
          * the triangular part */
         size++;
         rn[size] = cind[kk];
         cn[size] = j;
drop:    /* remove j-th column from the active submatrix */
         xassert(flag[j]);
         flag[j] = 0;
         if (prev[j] == 0)
            head = next[j];
         else
            next[prev[j]] = next[j];
         if (next[j] == 0)
            ;
         else
            prev[next[j]] = prev[j];
         /* decrease row counts */
         for (k = 1; k <= len; k++)
         {  i = cind[k];
            xassert(1 <= i && i <= m);
            xassert(cnt[i] > 0);
            cnt[i]--;
            if (cnt[i] == 1)
            {  /* new singleton appeared in i-th row; determine number
                * of corresponding column (it is the only active column
                * in this row) */
               len2 = mat(info, +i, rind, rval);
               xassert(0 <= len2 && len2 <= n);
               ks = 0;
               for (kk = 1; kk <= len2; kk++)
               {  jj = rind[kk];
                  xassert(1 <= jj && jj <= n);
                  if (flag[jj])
                  {  xassert(ks == 0);
                     ks = kk;
                  }
               }
               xassert(ks > 0);
               /* a[i,jj] is new row singleton */
               jj = rind[ks];
               if (flag[jj] != 2)
               {  /* include jj-th column in the singleton list */
                  flag[jj] = 2;
                  list[++ns] = jj;
               }
            }
         }
      }
      /* now all row counts should be zero */
      for (i = 1; i <= m; i++)
         xassert(cnt[i] == 0);
      /* deallocate working arrays */
      tfree(cind);
      tfree(cval);
      tfree(rind);
      tfree(rval);
      tfree(ptr);
      tfree(list);
      tfree(prev);
      tfree(next);
      tfree(big);
      tfree(flag);
      return size;
}

/* eof */
