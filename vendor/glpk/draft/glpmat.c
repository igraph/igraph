/* glpmat.c */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2000-2013 Free Software Foundation, Inc.
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
#include "glpmat.h"
#include "qmd.h"
#include "amd.h"
#include "colamd.h"

/*----------------------------------------------------------------------
-- check_fvs - check sparse vector in full-vector storage format.
--
-- SYNOPSIS
--
-- #include "glpmat.h"
-- int check_fvs(int n, int nnz, int ind[], double vec[]);
--
-- DESCRIPTION
--
-- The routine check_fvs checks if a given vector of dimension n in
-- full-vector storage format has correct representation.
--
-- RETURNS
--
-- The routine returns one of the following codes:
--
-- 0 - the vector is correct;
-- 1 - the number of elements (n) is negative;
-- 2 - the number of non-zero elements (nnz) is negative;
-- 3 - some element index is out of range;
-- 4 - some element index is duplicate;
-- 5 - some non-zero element is out of pattern. */

int check_fvs(int n, int nnz, int ind[], double vec[])
{     int i, t, ret, *flag = NULL;
      /* check the number of elements */
      if (n < 0)
      {  ret = 1;
         goto done;
      }
      /* check the number of non-zero elements */
      if (nnz < 0)
      {  ret = 2;
         goto done;
      }
      /* check vector indices */
      flag = xcalloc(1+n, sizeof(int));
      for (i = 1; i <= n; i++) flag[i] = 0;
      for (t = 1; t <= nnz; t++)
      {  i = ind[t];
         if (!(1 <= i && i <= n))
         {  ret = 3;
            goto done;
         }
         if (flag[i])
         {  ret = 4;
            goto done;
         }
         flag[i] = 1;
      }
      /* check vector elements */
      for (i = 1; i <= n; i++)
      {  if (!flag[i] && vec[i] != 0.0)
         {  ret = 5;
            goto done;
         }
      }
      /* the vector is ok */
      ret = 0;
done: if (flag != NULL) xfree(flag);
      return ret;
}

/*----------------------------------------------------------------------
-- check_pattern - check pattern of sparse matrix.
--
-- SYNOPSIS
--
-- #include "glpmat.h"
-- int check_pattern(int m, int n, int A_ptr[], int A_ind[]);
--
-- DESCRIPTION
--
-- The routine check_pattern checks the pattern of a given mxn matrix
-- in storage-by-rows format.
--
-- RETURNS
--
-- The routine returns one of the following codes:
--
-- 0 - the pattern is correct;
-- 1 - the number of rows (m) is negative;
-- 2 - the number of columns (n) is negative;
-- 3 - A_ptr[1] is not 1;
-- 4 - some column index is out of range;
-- 5 - some column indices are duplicate. */

int check_pattern(int m, int n, int A_ptr[], int A_ind[])
{     int i, j, ptr, ret, *flag = NULL;
      /* check the number of rows */
      if (m < 0)
      {  ret = 1;
         goto done;
      }
      /* check the number of columns */
      if (n < 0)
      {  ret = 2;
         goto done;
      }
      /* check location A_ptr[1] */
      if (A_ptr[1] != 1)
      {  ret = 3;
         goto done;
      }
      /* check row patterns */
      flag = xcalloc(1+n, sizeof(int));
      for (j = 1; j <= n; j++) flag[j] = 0;
      for (i = 1; i <= m; i++)
      {  /* check pattern of row i */
         for (ptr = A_ptr[i]; ptr < A_ptr[i+1]; ptr++)
         {  j = A_ind[ptr];
            /* check column index */
            if (!(1 <= j && j <= n))
            {  ret = 4;
               goto done;
            }
            /* check for duplication */
            if (flag[j])
            {  ret = 5;
               goto done;
            }
            flag[j] = 1;
         }
         /* clear flags */
         for (ptr = A_ptr[i]; ptr < A_ptr[i+1]; ptr++)
         {  j = A_ind[ptr];
            flag[j] = 0;
         }
      }
      /* the pattern is ok */
      ret = 0;
done: if (flag != NULL) xfree(flag);
      return ret;
}

/*----------------------------------------------------------------------
-- transpose - transpose sparse matrix.
--
-- *Synopsis*
--
-- #include "glpmat.h"
-- void transpose(int m, int n, int A_ptr[], int A_ind[],
--    double A_val[], int AT_ptr[], int AT_ind[], double AT_val[]);
--
-- *Description*
--
-- For a given mxn sparse matrix A the routine transpose builds a nxm
-- sparse matrix A' which is a matrix transposed to A.
--
-- The arrays A_ptr, A_ind, and A_val specify a given mxn matrix A to
-- be transposed in storage-by-rows format. The parameter A_val can be
-- NULL, in which case numeric values are not copied. The arrays A_ptr,
-- A_ind, and A_val are not changed on exit.
--
-- On entry the arrays AT_ptr, AT_ind, and AT_val must be allocated,
-- but their content is ignored. On exit the routine stores a resultant
-- nxm matrix A' in these arrays in storage-by-rows format. Note that
-- if the parameter A_val is NULL, the array AT_val is not used.
--
-- The routine transpose has a side effect that elements in rows of the
-- resultant matrix A' follow in ascending their column indices. */

void transpose(int m, int n, int A_ptr[], int A_ind[], double A_val[],
      int AT_ptr[], int AT_ind[], double AT_val[])
{     int i, j, t, beg, end, pos, len;
      /* determine row lengths of resultant matrix */
      for (j = 1; j <= n; j++) AT_ptr[j] = 0;
      for (i = 1; i <= m; i++)
      {  beg = A_ptr[i], end = A_ptr[i+1];
         for (t = beg; t < end; t++) AT_ptr[A_ind[t]]++;
      }
      /* set up row pointers of resultant matrix */
      pos = 1;
      for (j = 1; j <= n; j++)
         len = AT_ptr[j], pos += len, AT_ptr[j] = pos;
      AT_ptr[n+1] = pos;
      /* build resultant matrix */
      for (i = m; i >= 1; i--)
      {  beg = A_ptr[i], end = A_ptr[i+1];
         for (t = beg; t < end; t++)
         {  pos = --AT_ptr[A_ind[t]];
            AT_ind[pos] = i;
            if (A_val != NULL) AT_val[pos] = A_val[t];
         }
      }
      return;
}

/*----------------------------------------------------------------------
-- adat_symbolic - compute S = P*A*D*A'*P' (symbolic phase).
--
-- *Synopsis*
--
-- #include "glpmat.h"
-- int *adat_symbolic(int m, int n, int P_per[], int A_ptr[],
--    int A_ind[], int S_ptr[]);
--
-- *Description*
--
-- The routine adat_symbolic implements the symbolic phase to compute
-- symmetric matrix S = P*A*D*A'*P', where P is a permutation matrix,
-- A is a given sparse matrix, D is a diagonal matrix, A' is a matrix
-- transposed to A, P' is an inverse of P.
--
-- The parameter m is the number of rows in A and the order of P.
--
-- The parameter n is the number of columns in A and the order of D.
--
-- The array P_per specifies permutation matrix P. It is not changed on
-- exit.
--
-- The arrays A_ptr and A_ind specify the pattern of matrix A. They are
-- not changed on exit.
--
-- On exit the routine stores the pattern of upper triangular part of
-- matrix S without diagonal elements in the arrays S_ptr and S_ind in
-- storage-by-rows format. The array S_ptr should be allocated on entry,
-- however, its content is ignored. The array S_ind is allocated by the
-- routine itself which returns a pointer to it.
--
-- *Returns*
--
-- The routine returns a pointer to the array S_ind. */

int *adat_symbolic(int m, int n, int P_per[], int A_ptr[], int A_ind[],
      int S_ptr[])
{     int i, j, t, ii, jj, tt, k, size, len;
      int *S_ind, *AT_ptr, *AT_ind, *ind, *map, *temp;
      /* build the pattern of A', which is a matrix transposed to A, to
         efficiently access A in column-wise manner */
      AT_ptr = xcalloc(1+n+1, sizeof(int));
      AT_ind = xcalloc(A_ptr[m+1], sizeof(int));
      transpose(m, n, A_ptr, A_ind, NULL, AT_ptr, AT_ind, NULL);
      /* allocate the array S_ind */
      size = A_ptr[m+1] - 1;
      if (size < m) size = m;
      S_ind = xcalloc(1+size, sizeof(int));
      /* allocate and initialize working arrays */
      ind = xcalloc(1+m, sizeof(int));
      map = xcalloc(1+m, sizeof(int));
      for (jj = 1; jj <= m; jj++) map[jj] = 0;
      /* compute pattern of S; note that symbolically S = B*B', where
         B = P*A, B' is matrix transposed to B */
      S_ptr[1] = 1;
      for (ii = 1; ii <= m; ii++)
      {  /* compute pattern of ii-th row of S */
         len = 0;
         i = P_per[ii]; /* i-th row of A = ii-th row of B */
         for (t = A_ptr[i]; t < A_ptr[i+1]; t++)
         {  k = A_ind[t];
            /* walk through k-th column of A */
            for (tt = AT_ptr[k]; tt < AT_ptr[k+1]; tt++)
            {  j = AT_ind[tt];
               jj = P_per[m+j]; /* j-th row of A = jj-th row of B */
               /* a[i,k] != 0 and a[j,k] != 0 ergo s[ii,jj] != 0 */
               if (ii < jj && !map[jj]) ind[++len] = jj, map[jj] = 1;
            }
         }
         /* now (ind) is pattern of ii-th row of S */
         S_ptr[ii+1] = S_ptr[ii] + len;
         /* at least (S_ptr[ii+1] - 1) locations should be available in
            the array S_ind */
         if (S_ptr[ii+1] - 1 > size)
         {  temp = S_ind;
            size += size;
            S_ind = xcalloc(1+size, sizeof(int));
            memcpy(&S_ind[1], &temp[1], (S_ptr[ii] - 1) * sizeof(int));
            xfree(temp);
         }
         xassert(S_ptr[ii+1] - 1 <= size);
         /* (ii-th row of S) := (ind) */
         memcpy(&S_ind[S_ptr[ii]], &ind[1], len * sizeof(int));
         /* clear the row pattern map */
         for (t = 1; t <= len; t++) map[ind[t]] = 0;
      }
      /* free working arrays */
      xfree(AT_ptr);
      xfree(AT_ind);
      xfree(ind);
      xfree(map);
      /* reallocate the array S_ind to free unused locations */
      temp = S_ind;
      size = S_ptr[m+1] - 1;
      S_ind = xcalloc(1+size, sizeof(int));
      memcpy(&S_ind[1], &temp[1], size * sizeof(int));
      xfree(temp);
      return S_ind;
}

/*----------------------------------------------------------------------
-- adat_numeric - compute S = P*A*D*A'*P' (numeric phase).
--
-- *Synopsis*
--
-- #include "glpmat.h"
-- void adat_numeric(int m, int n, int P_per[],
--    int A_ptr[], int A_ind[], double A_val[], double D_diag[],
--    int S_ptr[], int S_ind[], double S_val[], double S_diag[]);
--
-- *Description*
--
-- The routine adat_numeric implements the numeric phase to compute
-- symmetric matrix S = P*A*D*A'*P', where P is a permutation matrix,
-- A is a given sparse matrix, D is a diagonal matrix, A' is a matrix
-- transposed to A, P' is an inverse of P.
--
-- The parameter m is the number of rows in A and the order of P.
--
-- The parameter n is the number of columns in A and the order of D.
--
-- The matrix P is specified in the array P_per, which is not changed
-- on exit.
--
-- The matrix A is specified in the arrays A_ptr, A_ind, and A_val in
-- storage-by-rows format. These arrays are not changed on exit.
--
-- Diagonal elements of the matrix D are specified in the array D_diag,
-- where D_diag[0] is not used, D_diag[i] = d[i,i] for i = 1, ..., n.
-- The array D_diag is not changed on exit.
--
-- The pattern of the upper triangular part of the matrix S without
-- diagonal elements (previously computed by the routine adat_symbolic)
-- is specified in the arrays S_ptr and S_ind, which are not changed on
-- exit. Numeric values of non-diagonal elements of S are stored in
-- corresponding locations of the array S_val, and values of diagonal
-- elements of S are stored in locations S_diag[1], ..., S_diag[n]. */

void adat_numeric(int m, int n, int P_per[],
      int A_ptr[], int A_ind[], double A_val[], double D_diag[],
      int S_ptr[], int S_ind[], double S_val[], double S_diag[])
{     int i, j, t, ii, jj, tt, beg, end, beg1, end1, k;
      double sum, *work;
      work = xcalloc(1+n, sizeof(double));
      for (j = 1; j <= n; j++) work[j] = 0.0;
      /* compute S = B*D*B', where B = P*A, B' is a matrix transposed
         to B */
      for (ii = 1; ii <= m; ii++)
      {  i = P_per[ii]; /* i-th row of A = ii-th row of B */
         /* (work) := (i-th row of A) */
         beg = A_ptr[i], end = A_ptr[i+1];
         for (t = beg; t < end; t++)
            work[A_ind[t]] = A_val[t];
         /* compute ii-th row of S */
         beg = S_ptr[ii], end = S_ptr[ii+1];
         for (t = beg; t < end; t++)
         {  jj = S_ind[t];
            j = P_per[jj]; /* j-th row of A = jj-th row of B */
            /* s[ii,jj] := sum a[i,k] * d[k,k] * a[j,k] */
            sum = 0.0;
            beg1 = A_ptr[j], end1 = A_ptr[j+1];
            for (tt = beg1; tt < end1; tt++)
            {  k = A_ind[tt];
               sum += work[k] * D_diag[k] * A_val[tt];
            }
            S_val[t] = sum;
         }
         /* s[ii,ii] := sum a[i,k] * d[k,k] * a[i,k] */
         sum = 0.0;
         beg = A_ptr[i], end = A_ptr[i+1];
         for (t = beg; t < end; t++)
         {  k = A_ind[t];
            sum += A_val[t] * D_diag[k] * A_val[t];
            work[k] = 0.0;
         }
         S_diag[ii] = sum;
      }
      xfree(work);
      return;
}

/*----------------------------------------------------------------------
-- min_degree - minimum degree ordering.
--
-- *Synopsis*
--
-- #include "glpmat.h"
-- void min_degree(int n, int A_ptr[], int A_ind[], int P_per[]);
--
-- *Description*
--
-- The routine min_degree uses the minimum degree ordering algorithm
-- to find a permutation matrix P for a given sparse symmetric positive
-- matrix A which minimizes the number of non-zeros in upper triangular
-- factor U for Cholesky factorization P*A*P' = U'*U.
--
-- The parameter n is the order of matrices A and P.
--
-- The pattern of the given matrix A is specified on entry in the arrays
-- A_ptr and A_ind in storage-by-rows format. Only the upper triangular
-- part without diagonal elements (which all are assumed to be non-zero)
-- should be specified as if A were upper triangular. The arrays A_ptr
-- and A_ind are not changed on exit.
--
-- The permutation matrix P is stored by the routine in the array P_per
-- on exit.
--
-- *Algorithm*
--
-- The routine min_degree is based on some subroutines from the package
-- SPARSPAK (see comments in the module glpqmd). */

void min_degree(int n, int A_ptr[], int A_ind[], int P_per[])
{     int i, j, ne, t, pos, len;
      int *xadj, *adjncy, *deg, *marker, *rchset, *nbrhd, *qsize,
         *qlink, nofsub;
      /* determine number of non-zeros in complete pattern */
      ne = A_ptr[n+1] - 1;
      ne += ne;
      /* allocate working arrays */
      xadj = xcalloc(1+n+1, sizeof(int));
      adjncy = xcalloc(1+ne, sizeof(int));
      deg = xcalloc(1+n, sizeof(int));
      marker = xcalloc(1+n, sizeof(int));
      rchset = xcalloc(1+n, sizeof(int));
      nbrhd = xcalloc(1+n, sizeof(int));
      qsize = xcalloc(1+n, sizeof(int));
      qlink = xcalloc(1+n, sizeof(int));
      /* determine row lengths in complete pattern */
      for (i = 1; i <= n; i++) xadj[i] = 0;
      for (i = 1; i <= n; i++)
      {  for (t = A_ptr[i]; t < A_ptr[i+1]; t++)
         {  j = A_ind[t];
            xassert(i < j && j <= n);
            xadj[i]++, xadj[j]++;
         }
      }
      /* set up row pointers for complete pattern */
      pos = 1;
      for (i = 1; i <= n; i++)
         len = xadj[i], pos += len, xadj[i] = pos;
      xadj[n+1] = pos;
      xassert(pos - 1 == ne);
      /* construct complete pattern */
      for (i = 1; i <= n; i++)
      {  for (t = A_ptr[i]; t < A_ptr[i+1]; t++)
         {  j = A_ind[t];
            adjncy[--xadj[i]] = j, adjncy[--xadj[j]] = i;
         }
      }
      /* call the main minimimum degree ordering routine */
      genqmd(&n, xadj, adjncy, P_per, P_per + n, deg, marker, rchset,
         nbrhd, qsize, qlink, &nofsub);
      /* make sure that permutation matrix P is correct */
      for (i = 1; i <= n; i++)
      {  j = P_per[i];
         xassert(1 <= j && j <= n);
         xassert(P_per[n+j] == i);
      }
      /* free working arrays */
      xfree(xadj);
      xfree(adjncy);
      xfree(deg);
      xfree(marker);
      xfree(rchset);
      xfree(nbrhd);
      xfree(qsize);
      xfree(qlink);
      return;
}

/**********************************************************************/

void amd_order1(int n, int A_ptr[], int A_ind[], int P_per[])
{     /* approximate minimum degree ordering (AMD) */
      int k, ret;
      double Control[AMD_CONTROL], Info[AMD_INFO];
      /* get the default parameters */
      amd_defaults(Control);
#if 0
      /* and print them */
      amd_control(Control);
#endif
      /* make all indices 0-based */
      for (k = 1; k < A_ptr[n+1]; k++) A_ind[k]--;
      for (k = 1; k <= n+1; k++) A_ptr[k]--;
      /* call the ordering routine */
      ret = amd_order(n, &A_ptr[1], &A_ind[1], &P_per[1], Control, Info)
         ;
#if 0
      amd_info(Info);
#endif
      xassert(ret == AMD_OK || ret == AMD_OK_BUT_JUMBLED);
      /* retsore 1-based indices */
      for (k = 1; k <= n+1; k++) A_ptr[k]++;
      for (k = 1; k < A_ptr[n+1]; k++) A_ind[k]++;
      /* patch up permutation matrix */
      memset(&P_per[n+1], 0, n * sizeof(int));
      for (k = 1; k <= n; k++)
      {  P_per[k]++;
         xassert(1 <= P_per[k] && P_per[k] <= n);
         xassert(P_per[n+P_per[k]] == 0);
         P_per[n+P_per[k]] = k;
      }
      return;
}

/**********************************************************************/

static void *allocate(size_t n, size_t size)
{     void *ptr;
      ptr = xcalloc(n, size);
      memset(ptr, 0, n * size);
      return ptr;
}

static void release(void *ptr)
{     xfree(ptr);
      return;
}

void symamd_ord(int n, int A_ptr[], int A_ind[], int P_per[])
{     /* approximate minimum degree ordering (SYMAMD) */
      int k, ok;
      int stats[COLAMD_STATS];
      /* make all indices 0-based */
      for (k = 1; k < A_ptr[n+1]; k++) A_ind[k]--;
      for (k = 1; k <= n+1; k++) A_ptr[k]--;
      /* call the ordering routine */
      ok = symamd(n, &A_ind[1], &A_ptr[1], &P_per[1], NULL, stats,
         allocate, release);
#if 0
      symamd_report(stats);
#endif
      xassert(ok);
      /* restore 1-based indices */
      for (k = 1; k <= n+1; k++) A_ptr[k]++;
      for (k = 1; k < A_ptr[n+1]; k++) A_ind[k]++;
      /* patch up permutation matrix */
      memset(&P_per[n+1], 0, n * sizeof(int));
      for (k = 1; k <= n; k++)
      {  P_per[k]++;
         xassert(1 <= P_per[k] && P_per[k] <= n);
         xassert(P_per[n+P_per[k]] == 0);
         P_per[n+P_per[k]] = k;
      }
      return;
}

/*----------------------------------------------------------------------
-- chol_symbolic - compute Cholesky factorization (symbolic phase).
--
-- *Synopsis*
--
-- #include "glpmat.h"
-- int *chol_symbolic(int n, int A_ptr[], int A_ind[], int U_ptr[]);
--
-- *Description*
--
-- The routine chol_symbolic implements the symbolic phase of Cholesky
-- factorization A = U'*U, where A is a given sparse symmetric positive
-- definite matrix, U is a resultant upper triangular factor, U' is a
-- matrix transposed to U.
--
-- The parameter n is the order of matrices A and U.
--
-- The pattern of the given matrix A is specified on entry in the arrays
-- A_ptr and A_ind in storage-by-rows format. Only the upper triangular
-- part without diagonal elements (which all are assumed to be non-zero)
-- should be specified as if A were upper triangular. The arrays A_ptr
-- and A_ind are not changed on exit.
--
-- The pattern of the matrix U without diagonal elements (which all are
-- assumed to be non-zero) is stored on exit from the routine in the
-- arrays U_ptr and U_ind in storage-by-rows format. The array U_ptr
-- should be allocated on entry, however, its content is ignored. The
-- array U_ind is allocated by the routine which returns a pointer to it
-- on exit.
--
-- *Returns*
--
-- The routine returns a pointer to the array U_ind.
--
-- *Method*
--
-- The routine chol_symbolic computes the pattern of the matrix U in a
-- row-wise manner. No pivoting is used.
--
-- It is known that to compute the pattern of row k of the matrix U we
-- need to merge the pattern of row k of the matrix A and the patterns
-- of each row i of U, where u[i,k] is non-zero (these rows are already
-- computed and placed above row k).
--
-- However, to reduce the number of rows to be merged the routine uses
-- an advanced algorithm proposed in:
--
-- D.J.Rose, R.E.Tarjan, and G.S.Lueker. Algorithmic aspects of vertex
-- elimination on graphs. SIAM J. Comput. 5, 1976, 266-83.
--
-- The authors of the cited paper show that we have the same result if
-- we merge row k of the matrix A and such rows of the matrix U (among
-- rows 1, ..., k-1) whose leftmost non-diagonal non-zero element is
-- placed in k-th column. This feature signficantly reduces the number
-- of rows to be merged, especially on the final steps, where rows of
-- the matrix U become quite dense.
--
-- To determine rows, which should be merged on k-th step, for a fixed
-- time the routine uses linked lists of row numbers of the matrix U.
-- Location head[k] contains the number of a first row, whose leftmost
-- non-diagonal non-zero element is placed in column k, and location
-- next[i] contains the number of a next row with the same property as
-- row i. */

int *chol_symbolic(int n, int A_ptr[], int A_ind[], int U_ptr[])
{     int i, j, k, t, len, size, beg, end, min_j, *U_ind, *head, *next,
         *ind, *map, *temp;
      /* initially we assume that on computing the pattern of U fill-in
         will double the number of non-zeros in A */
      size = A_ptr[n+1] - 1;
      if (size < n) size = n;
      size += size;
      U_ind = xcalloc(1+size, sizeof(int));
      /* allocate and initialize working arrays */
      head = xcalloc(1+n, sizeof(int));
      for (i = 1; i <= n; i++) head[i] = 0;
      next = xcalloc(1+n, sizeof(int));
      ind = xcalloc(1+n, sizeof(int));
      map = xcalloc(1+n, sizeof(int));
      for (j = 1; j <= n; j++) map[j] = 0;
      /* compute the pattern of matrix U */
      U_ptr[1] = 1;
      for (k = 1; k <= n; k++)
      {  /* compute the pattern of k-th row of U, which is the union of
            k-th row of A and those rows of U (among 1, ..., k-1) whose
            leftmost non-diagonal non-zero is placed in k-th column */
         /* (ind) := (k-th row of A) */
         len = A_ptr[k+1] - A_ptr[k];
         memcpy(&ind[1], &A_ind[A_ptr[k]], len * sizeof(int));
         for (t = 1; t <= len; t++)
         {  j = ind[t];
            xassert(k < j && j <= n);
            map[j] = 1;
         }
         /* walk through rows of U whose leftmost non-diagonal non-zero
            is placed in k-th column */
         for (i = head[k]; i != 0; i = next[i])
         {  /* (ind) := (ind) union (i-th row of U) */
            beg = U_ptr[i], end = U_ptr[i+1];
            for (t = beg; t < end; t++)
            {  j = U_ind[t];
               if (j > k && !map[j]) ind[++len] = j, map[j] = 1;
            }
         }
         /* now (ind) is the pattern of k-th row of U */
         U_ptr[k+1] = U_ptr[k] + len;
         /* at least (U_ptr[k+1] - 1) locations should be available in
            the array U_ind */
         if (U_ptr[k+1] - 1 > size)
         {  temp = U_ind;
            size += size;
            U_ind = xcalloc(1+size, sizeof(int));
            memcpy(&U_ind[1], &temp[1], (U_ptr[k] - 1) * sizeof(int));
            xfree(temp);
         }
         xassert(U_ptr[k+1] - 1 <= size);
         /* (k-th row of U) := (ind) */
         memcpy(&U_ind[U_ptr[k]], &ind[1], len * sizeof(int));
         /* determine column index of leftmost non-diagonal non-zero in
            k-th row of U and clear the row pattern map */
         min_j = n + 1;
         for (t = 1; t <= len; t++)
         {  j = ind[t], map[j] = 0;
            if (min_j > j) min_j = j;
         }
         /* include k-th row into corresponding linked list */
         if (min_j <= n) next[k] = head[min_j], head[min_j] = k;
      }
      /* free working arrays */
      xfree(head);
      xfree(next);
      xfree(ind);
      xfree(map);
      /* reallocate the array U_ind to free unused locations */
      temp = U_ind;
      size = U_ptr[n+1] - 1;
      U_ind = xcalloc(1+size, sizeof(int));
      memcpy(&U_ind[1], &temp[1], size * sizeof(int));
      xfree(temp);
      return U_ind;
}

/*----------------------------------------------------------------------
-- chol_numeric - compute Cholesky factorization (numeric phase).
--
-- *Synopsis*
--
-- #include "glpmat.h"
-- int chol_numeric(int n,
--    int A_ptr[], int A_ind[], double A_val[], double A_diag[],
--    int U_ptr[], int U_ind[], double U_val[], double U_diag[]);
--
-- *Description*
--
-- The routine chol_symbolic implements the numeric phase of Cholesky
-- factorization A = U'*U, where A is a given sparse symmetric positive
-- definite matrix, U is a resultant upper triangular factor, U' is a
-- matrix transposed to U.
--
-- The parameter n is the order of matrices A and U.
--
-- Upper triangular part of the matrix A without diagonal elements is
-- specified in the arrays A_ptr, A_ind, and A_val in storage-by-rows
-- format. Diagonal elements of A are specified in the array A_diag,
-- where A_diag[0] is not used, A_diag[i] = a[i,i] for i = 1, ..., n.
-- The arrays A_ptr, A_ind, A_val, and A_diag are not changed on exit.
--
-- The pattern of the matrix U without diagonal elements (previously
-- computed with the routine chol_symbolic) is specified in the arrays
-- U_ptr and U_ind, which are not changed on exit. Numeric values of
-- non-diagonal elements of U are stored in corresponding locations of
-- the array U_val, and values of diagonal elements of U are stored in
-- locations U_diag[1], ..., U_diag[n].
--
-- *Returns*
--
-- The routine returns the number of non-positive diagonal elements of
-- the matrix U which have been replaced by a huge positive number (see
-- the method description below). Zero return code means the matrix A
-- has been successfully factorized.
--
-- *Method*
--
-- The routine chol_numeric computes the matrix U in a row-wise manner
-- using standard gaussian elimination technique. No pivoting is used.
--
-- Initially the routine sets U = A, and before k-th elimination step
-- the matrix U is the following:
--
--       1       k         n
--    1  x x x x x x x x x x
--       . x x x x x x x x x
--       . . x x x x x x x x
--       . . . x x x x x x x
--    k  . . . . * * * * * *
--       . . . . * * * * * *
--       . . . . * * * * * *
--       . . . . * * * * * *
--       . . . . * * * * * *
--    n  . . . . * * * * * *
--
-- where 'x' are elements of already computed rows, '*' are elements of
-- the active submatrix. (Note that the lower triangular part of the
-- active submatrix being symmetric is not stored and diagonal elements
-- are stored separately in the array U_diag.)
--
-- The matrix A is assumed to be positive definite. However, if it is
-- close to semi-definite, on some elimination step a pivot u[k,k] may
-- happen to be non-positive due to round-off errors. In this case the
-- routine uses a technique proposed in:
--
-- S.J.Wright. The Cholesky factorization in interior-point and barrier
-- methods. Preprint MCS-P600-0596, Mathematics and Computer Science
-- Division, Argonne National Laboratory, Argonne, Ill., May 1996.
--
-- The routine just replaces non-positive u[k,k] by a huge positive
-- number. This involves non-diagonal elements in k-th row of U to be
-- close to zero that, in turn, involves k-th component of a solution
-- vector to be close to zero. Note, however, that this technique works
-- only if the system A*x = b is consistent. */

int chol_numeric(int n,
      int A_ptr[], int A_ind[], double A_val[], double A_diag[],
      int U_ptr[], int U_ind[], double U_val[], double U_diag[])
{     int i, j, k, t, t1, beg, end, beg1, end1, count = 0;
      double ukk, uki, *work;
      work = xcalloc(1+n, sizeof(double));
      for (j = 1; j <= n; j++) work[j] = 0.0;
      /* U := (upper triangle of A) */
      /* note that the upper traingle of A is a subset of U */
      for (i = 1; i <= n; i++)
      {  beg = A_ptr[i], end = A_ptr[i+1];
         for (t = beg; t < end; t++)
            j = A_ind[t], work[j] = A_val[t];
         beg = U_ptr[i], end = U_ptr[i+1];
         for (t = beg; t < end; t++)
            j = U_ind[t], U_val[t] = work[j], work[j] = 0.0;
         U_diag[i] = A_diag[i];
      }
      /* main elimination loop */
      for (k = 1; k <= n; k++)
      {  /* transform k-th row of U */
         ukk = U_diag[k];
         if (ukk > 0.0)
            U_diag[k] = ukk = sqrt(ukk);
         else
            U_diag[k] = ukk = DBL_MAX, count++;
         /* (work) := (transformed k-th row) */
         beg = U_ptr[k], end = U_ptr[k+1];
         for (t = beg; t < end; t++)
            work[U_ind[t]] = (U_val[t] /= ukk);
         /* transform other rows of U */
         for (t = beg; t < end; t++)
         {  i = U_ind[t];
            xassert(i > k);
            /* (i-th row) := (i-th row) - u[k,i] * (k-th row) */
            uki = work[i];
            beg1 = U_ptr[i], end1 = U_ptr[i+1];
            for (t1 = beg1; t1 < end1; t1++)
               U_val[t1] -= uki * work[U_ind[t1]];
            U_diag[i] -= uki * uki;
         }
         /* (work) := 0 */
         for (t = beg; t < end; t++)
            work[U_ind[t]] = 0.0;
      }
      xfree(work);
      return count;
}

/*----------------------------------------------------------------------
-- u_solve - solve upper triangular system U*x = b.
--
-- *Synopsis*
--
-- #include "glpmat.h"
-- void u_solve(int n, int U_ptr[], int U_ind[], double U_val[],
--    double U_diag[], double x[]);
--
-- *Description*
--
-- The routine u_solve solves an linear system U*x = b, where U is an
-- upper triangular matrix.
--
-- The parameter n is the order of matrix U.
--
-- The matrix U without diagonal elements is specified in the arrays
-- U_ptr, U_ind, and U_val in storage-by-rows format. Diagonal elements
-- of U are specified in the array U_diag, where U_diag[0] is not used,
-- U_diag[i] = u[i,i] for i = 1, ..., n. All these four arrays are not
-- changed on exit.
--
-- The right-hand side vector b is specified on entry in the array x,
-- where x[0] is not used, and x[i] = b[i] for i = 1, ..., n. On exit
-- the routine stores computed components of the vector of unknowns x
-- in the array x in the same manner. */

void u_solve(int n, int U_ptr[], int U_ind[], double U_val[],
      double U_diag[], double x[])
{     int i, t, beg, end;
      double temp;
      for (i = n; i >= 1; i--)
      {  temp = x[i];
         beg = U_ptr[i], end = U_ptr[i+1];
         for (t = beg; t < end; t++)
            temp -= U_val[t] * x[U_ind[t]];
         xassert(U_diag[i] != 0.0);
         x[i] = temp / U_diag[i];
      }
      return;
}

/*----------------------------------------------------------------------
-- ut_solve - solve lower triangular system U'*x = b.
--
-- *Synopsis*
--
-- #include "glpmat.h"
-- void ut_solve(int n, int U_ptr[], int U_ind[], double U_val[],
--    double U_diag[], double x[]);
--
-- *Description*
--
-- The routine ut_solve solves an linear system U'*x = b, where U is a
-- matrix transposed to an upper triangular matrix.
--
-- The parameter n is the order of matrix U.
--
-- The matrix U without diagonal elements is specified in the arrays
-- U_ptr, U_ind, and U_val in storage-by-rows format. Diagonal elements
-- of U are specified in the array U_diag, where U_diag[0] is not used,
-- U_diag[i] = u[i,i] for i = 1, ..., n. All these four arrays are not
-- changed on exit.
--
-- The right-hand side vector b is specified on entry in the array x,
-- where x[0] is not used, and x[i] = b[i] for i = 1, ..., n. On exit
-- the routine stores computed components of the vector of unknowns x
-- in the array x in the same manner. */

void ut_solve(int n, int U_ptr[], int U_ind[], double U_val[],
      double U_diag[], double x[])
{     int i, t, beg, end;
      double temp;
      for (i = 1; i <= n; i++)
      {  xassert(U_diag[i] != 0.0);
         temp = (x[i] /= U_diag[i]);
         if (temp == 0.0) continue;
         beg = U_ptr[i], end = U_ptr[i+1];
         for (t = beg; t < end; t++)
            x[U_ind[t]] -= U_val[t] * temp;
      }
      return;
}

/* eof */
