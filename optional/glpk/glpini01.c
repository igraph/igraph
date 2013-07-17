/* glpini01.c */

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
#pragma clang diagnostic ignored "-Wself-assign"
#endif

#include "glpapi.h"

/*----------------------------------------------------------------------
-- triang - find maximal triangular part of a rectangular matrix.
--
-- *Synopsis*
--
-- int triang(int m, int n,
--    void *info, int (*mat)(void *info, int k, int ndx[]),
--    int rn[], int cn[]);
--
-- *Description*
--
-- For a given rectangular (sparse) matrix A with m rows and n columns
-- the routine triang tries to find such permutation matrices P and Q
-- that the first rows and columns of the matrix B = P*A*Q form a lower
-- triangular submatrix of as greatest size as possible:
--
--                   1                       n
--                1  * . . . . . . x x x x x x
--                   * * . . . . . x x x x x x
--                   * * * . . . . x x x x x x
--                   * * * * . . . x x x x x x
--    B = P*A*Q =    * * * * * . . x x x x x x
--                   * * * * * * . x x x x x x
--                   * * * * * * * x x x x x x
--                   x x x x x x x x x x x x x
--                   x x x x x x x x x x x x x
--                m  x x x x x x x x x x x x x
--
-- where: '*' - elements of the lower triangular part, '.' - structural
-- zeros, 'x' - other (either non-zero or zero) elements.
--
-- The parameter info is a transit pointer passed to the formal routine
-- mat (see below).
--
-- The formal routine mat specifies the given matrix A in both row- and
-- column-wise formats. In order to obtain an i-th row of the matrix A
-- the routine triang calls the routine mat with the parameter k = +i,
-- 1 <= i <= m. In response the routine mat should store column indices
-- of (non-zero) elements of the i-th row to the locations ndx[1], ...,
-- ndx[len], where len is number of non-zeros in the i-th row returned
-- on exit. Analogously, in order to obtain a j-th column of the matrix
-- A, the routine mat is called with the parameter k = -j, 1 <= j <= n,
-- and should return pattern of the j-th column in the same way as for
-- row patterns. Note that the routine mat may be called more than once
-- for the same rows and columns.
--
-- On exit the routine computes two resultant arrays rn and cn, which
-- define the permutation matrices P and Q, respectively. The array rn
-- should have at least 1+m locations, where rn[i] = i' (1 <= i <= m)
-- means that i-th row of the original matrix A corresponds to i'-th row
-- of the matrix B = P*A*Q. Similarly, the array cn should have at least
-- 1+n locations, where cn[j] = j' (1 <= j <= n) means that j-th column
-- of the matrix A corresponds to j'-th column of the matrix B.
--
-- *Returns*
--
-- The routine triang returns the size of the lower tringular part of
-- the matrix B = P*A*Q (see the figure above).
--
-- *Complexity*
--
-- The time complexity of the routine triang is O(nnz), where nnz is
-- number of non-zeros in the given matrix A.
--
-- *Algorithm*
--
-- The routine triang starts from the matrix B = P*Q*A, where P and Q
-- are unity matrices, so initially B = A.
--
-- Before the next iteration B = (B1 | B2 | B3), where B1 is partially
-- built a lower triangular submatrix, B2 is the active submatrix, and
-- B3 is a submatrix that contains rejected columns. Thus, the current
-- matrix B looks like follows (initially k1 = 1 and k2 = n):
--
--       1         k1         k2         n
--    1  x . . . . . . . . . . . . . # # #
--       x x . . . . . . . . . . . . # # #
--       x x x . . . . . . . . . . # # # #
--       x x x x . . . . . . . . . # # # #
--       x x x x x . . . . . . . # # # # #
--    k1 x x x x x * * * * * * * # # # # #
--       x x x x x * * * * * * * # # # # #
--       x x x x x * * * * * * * # # # # #
--       x x x x x * * * * * * * # # # # #
--    m  x x x x x * * * * * * * # # # # #
--       <--B1---> <----B2-----> <---B3-->
--
-- On each iteartion the routine looks for a singleton row, i.e. some
-- row that has the only non-zero in the active submatrix B2. If such
-- row exists and the corresponding non-zero is b[i,j], where (by the
-- definition) k1 <= i <= m and k1 <= j <= k2, the routine permutes
-- k1-th and i-th rows and k1-th and j-th columns of the matrix B (in
-- order to place the element in the position b[k1,k1]), removes the
-- k1-th column from the active submatrix B2, and adds this column to
-- the submatrix B1. If no row singletons exist, but B2 is not empty
-- yet, the routine chooses a j-th column, which has maximal number of
-- non-zeros among other columns of B2, removes this column from B2 and
-- adds it to the submatrix B3 in the hope that new row singletons will
-- appear in the active submatrix. */

static int triang(int m, int n,
      void *info, int (*mat)(void *info, int k, int ndx[]),
      int rn[], int cn[])
{     int *ndx; /* int ndx[1+max(m,n)]; */
      /* this array is used for querying row and column patterns of the
         given matrix A (the third parameter to the routine mat) */
      int *rs_len; /* int rs_len[1+m]; */
      /* rs_len[0] is not used;
         rs_len[i], 1 <= i <= m, is number of non-zeros in the i-th row
         of the matrix A, which (non-zeros) belong to the current active
         submatrix */
      int *rs_head; /* int rs_head[1+n]; */
      /* rs_head[len], 0 <= len <= n, is the number i of the first row
         of the matrix A, for which rs_len[i] = len */
      int *rs_prev; /* int rs_prev[1+m]; */
      /* rs_prev[0] is not used;
         rs_prev[i], 1 <= i <= m, is a number i' of the previous row of
         the matrix A, for which rs_len[i] = rs_len[i'] (zero marks the
         end of this linked list) */
      int *rs_next; /* int rs_next[1+m]; */
      /* rs_next[0] is not used;
         rs_next[i], 1 <= i <= m, is a number i' of the next row of the
         matrix A, for which rs_len[i] = rs_len[i'] (zero marks the end
         this linked list) */
      int cs_head;
      /* is a number j of the first column of the matrix A, which has
         maximal number of non-zeros among other columns */
      int *cs_prev; /* cs_prev[1+n]; */
      /* cs_prev[0] is not used;
         cs_prev[j], 1 <= j <= n, is a number of the previous column of
         the matrix A with the same or greater number of non-zeros than
         in the j-th column (zero marks the end of this linked list) */
      int *cs_next; /* cs_next[1+n]; */
      /* cs_next[0] is not used;
         cs_next[j], 1 <= j <= n, is a number of the next column of
         the matrix A with the same or lesser number of non-zeros than
         in the j-th column (zero marks the end of this linked list) */
      int i, j, ii, jj, k1, k2, len, t, size = 0;
      int *head, *rn_inv, *cn_inv;
      if (!(m > 0 && n > 0))
         xerror("triang: m = %d; n = %d; invalid dimension\n", m, n);
      /* allocate working arrays */
      ndx = xcalloc(1+(m >= n ? m : n), sizeof(int));
      rs_len = xcalloc(1+m, sizeof(int));
      rs_head = xcalloc(1+n, sizeof(int));
      rs_prev = xcalloc(1+m, sizeof(int));
      rs_next = xcalloc(1+m, sizeof(int));
      cs_prev = xcalloc(1+n, sizeof(int));
      cs_next = xcalloc(1+n, sizeof(int));
      /* build linked lists of columns of the matrix A with the same
         number of non-zeros */
      head = rs_len; /* currently rs_len is used as working array */
      for (len = 0; len <= m; len ++) head[len] = 0;
      for (j = 1; j <= n; j++)
      {  /* obtain length of the j-th column */
         len = mat(info, -j, ndx);
         xassert(0 <= len && len <= m);
         /* include the j-th column in the corresponding linked list */
         cs_prev[j] = head[len];
         head[len] = j;
      }
      /* merge all linked lists of columns in one linked list, where
         columns are ordered by descending of their lengths */
      cs_head = 0;
      for (len = 0; len <= m; len++)
      {  for (j = head[len]; j != 0; j = cs_prev[j])
         {  cs_next[j] = cs_head;
            cs_head = j;
         }
      }
      jj = 0;
      for (j = cs_head; j != 0; j = cs_next[j])
      {  cs_prev[j] = jj;
         jj = j;
      }
      /* build initial doubly linked lists of rows of the matrix A with
         the same number of non-zeros */
      for (len = 0; len <= n; len++) rs_head[len] = 0;
      for (i = 1; i <= m; i++)
      {  /* obtain length of the i-th row */
         rs_len[i] = len = mat(info, +i, ndx);
         xassert(0 <= len && len <= n);
         /* include the i-th row in the correspondng linked list */
         rs_prev[i] = 0;
         rs_next[i] = rs_head[len];
         if (rs_next[i] != 0) rs_prev[rs_next[i]] = i;
         rs_head[len] = i;
      }
      /* initially all rows and columns of the matrix A are active */
      for (i = 1; i <= m; i++) rn[i] = 0;
      for (j = 1; j <= n; j++) cn[j] = 0;
      /* set initial bounds of the active submatrix */
      k1 = 1, k2 = n;
      /* main loop starts here */
      while (k1 <= k2)
      {  i = rs_head[1];
         if (i != 0)
         {  /* the i-th row of the matrix A is a row singleton, since
               it has the only non-zero in the active submatrix */
            xassert(rs_len[i] == 1);
            /* determine the number j of an active column of the matrix
               A, in which this non-zero is placed */
            j = 0;
            t = mat(info, +i, ndx);
            xassert(0 <= t && t <= n);
            for (t = t; t >= 1; t--)
            {  jj = ndx[t];
               xassert(1 <= jj && jj <= n);
               if (cn[jj] == 0)
               {  xassert(j == 0);
                  j = jj;
               }
            }
            xassert(j != 0);
            /* the singleton is a[i,j]; move a[i,j] to the position
               b[k1,k1] of the matrix B */
            rn[i] = cn[j] = k1;
            /* shift the left bound of the active submatrix */
            k1++;
            /* increase the size of the lower triangular part */
            size++;
         }
         else
         {  /* the current active submatrix has no row singletons */
            /* remove an active column with maximal number of non-zeros
               from the active submatrix */
            j = cs_head;
            xassert(j != 0);
            cn[j] = k2;
            /* shift the right bound of the active submatrix */
            k2--;
         }
         /* the j-th column of the matrix A has been removed from the
            active submatrix */
         /* remove the j-th column from the linked list */
         if (cs_prev[j] == 0)
            cs_head = cs_next[j];
         else
            cs_next[cs_prev[j]] = cs_next[j];
         if (cs_next[j] == 0)
            /* nop */;
         else
            cs_prev[cs_next[j]] = cs_prev[j];
         /* go through non-zeros of the j-th columns and update active
            lengths of the corresponding rows */
         t = mat(info, -j, ndx);
         xassert(0 <= t && t <= m);
         for (t = t; t >= 1; t--)
         {  i = ndx[t];
            xassert(1 <= i && i <= m);
            /* the non-zero a[i,j] has left the active submatrix */
            len = rs_len[i];
            xassert(len >= 1);
            /* remove the i-th row from the linked list of rows with
               active length len */
            if (rs_prev[i] == 0)
               rs_head[len] = rs_next[i];
            else
               rs_next[rs_prev[i]] = rs_next[i];
            if (rs_next[i] == 0)
               /* nop */;
            else
               rs_prev[rs_next[i]] = rs_prev[i];
            /* decrease the active length of the i-th row */
            rs_len[i] = --len;
            /* return the i-th row to the corresponding linked list */
            rs_prev[i] = 0;
            rs_next[i] = rs_head[len];
            if (rs_next[i] != 0) rs_prev[rs_next[i]] = i;
            rs_head[len] = i;
         }
      }
      /* other rows of the matrix A, which are still active, correspond
         to rows k1, ..., m of the matrix B (in arbitrary order) */
      for (i = 1; i <= m; i++) if (rn[i] == 0) rn[i] = k1++;
      /* but for columns this is not needed, because now the submatrix
         B2 has no columns */
      for (j = 1; j <= n; j++) xassert(cn[j] != 0);
      /* perform some optional checks */
      /* make sure that rn is a permutation of {1, ..., m} and cn is a
         permutation of {1, ..., n} */
      rn_inv = rs_len; /* used as working array */
      for (ii = 1; ii <= m; ii++) rn_inv[ii] = 0;
      for (i = 1; i <= m; i++)
      {  ii = rn[i];
         xassert(1 <= ii && ii <= m);
         xassert(rn_inv[ii] == 0);
         rn_inv[ii] = i;
      }
      cn_inv = rs_head; /* used as working array */
      for (jj = 1; jj <= n; jj++) cn_inv[jj] = 0;
      for (j = 1; j <= n; j++)
      {  jj = cn[j];
         xassert(1 <= jj && jj <= n);
         xassert(cn_inv[jj] == 0);
         cn_inv[jj] = j;
      }
      /* make sure that the matrix B = P*A*Q really has the form, which
         was declared */
      for (ii = 1; ii <= size; ii++)
      {  int diag = 0;
         i = rn_inv[ii];
         t = mat(info, +i, ndx);
         xassert(0 <= t && t <= n);
         for (t = t; t >= 1; t--)
         {  j = ndx[t];
            xassert(1 <= j && j <= n);
            jj = cn[j];
            if (jj <= size) xassert(jj <= ii);
            if (jj == ii)
            {  xassert(!diag);
               diag = 1;
            }
         }
         xassert(diag);
      }
      /* free working arrays */
      xfree(ndx);
      xfree(rs_len);
      xfree(rs_head);
      xfree(rs_prev);
      xfree(rs_next);
      xfree(cs_prev);
      xfree(cs_next);
      /* return to the calling program */
      return size;
}

/*----------------------------------------------------------------------
-- adv_basis - construct advanced initial LP basis.
--
-- *Synopsis*
--
-- #include "glpini.h"
-- void adv_basis(glp_prob *lp);
--
-- *Description*
--
-- The routine adv_basis constructs an advanced initial basis for an LP
-- problem object, which the parameter lp points to.
--
-- In order to build the initial basis the routine does the following:
--
-- 1) includes in the basis all non-fixed auxiliary variables;
--
-- 2) includes in the basis as many as possible non-fixed structural
--    variables preserving triangular form of the basis matrix;
--
-- 3) includes in the basis appropriate (fixed) auxiliary variables
--    in order to complete the basis.
--
-- As a result the initial basis has minimum of fixed variables and the
-- corresponding basis matrix is triangular. */

static int mat(void *info, int k, int ndx[])
{     /* this auxiliary routine returns the pattern of a given row or
         a given column of the augmented constraint matrix A~ = (I|-A),
         in which columns of fixed variables are implicitly cleared */
      LPX *lp = info;
      int m = lpx_get_num_rows(lp);
      int n = lpx_get_num_cols(lp);
      int typx, i, j, lll, len = 0;
      if (k > 0)
      {  /* the pattern of the i-th row is required */
         i = +k;
         xassert(1 <= i && i <= m);
#if 0 /* 22/XII-2003 */
         /* if the auxiliary variable x[i] is non-fixed, include its
            element (placed in the i-th column) in the pattern */
         lpx_get_row_bnds(lp, i, &typx, NULL, NULL);
         if (typx != LPX_FX) ndx[++len] = i;
         /* include in the pattern elements placed in columns, which
            correspond to non-fixed structural varables */
         i_beg = aa_ptr[i];
         i_end = i_beg + aa_len[i] - 1;
         for (i_ptr = i_beg; i_ptr <= i_end; i_ptr++)
         {  j = m + sv_ndx[i_ptr];
            lpx_get_col_bnds(lp, j-m, &typx, NULL, NULL);
            if (typx != LPX_FX) ndx[++len] = j;
         }
#else
         lll = lpx_get_mat_row(lp, i, ndx, NULL);
         for (k = 1; k <= lll; k++)
         {  lpx_get_col_bnds(lp, ndx[k], &typx, NULL, NULL);
            if (typx != LPX_FX) ndx[++len] = m + ndx[k];
         }
         lpx_get_row_bnds(lp, i, &typx, NULL, NULL);
         if (typx != LPX_FX) ndx[++len] = i;
#endif
      }
      else
      {  /* the pattern of the j-th column is required */
         j = -k;
         xassert(1 <= j && j <= m+n);
         /* if the (auxiliary or structural) variable x[j] is fixed,
            the pattern of its column is empty */
         if (j <= m)
            lpx_get_row_bnds(lp, j, &typx, NULL, NULL);
         else
            lpx_get_col_bnds(lp, j-m, &typx, NULL, NULL);
         if (typx != LPX_FX)
         {  if (j <= m)
            {  /* x[j] is non-fixed auxiliary variable */
               ndx[++len] = j;
            }
            else
            {  /* x[j] is non-fixed structural variables */
#if 0 /* 22/XII-2003 */
               j_beg = aa_ptr[j];
               j_end = j_beg + aa_len[j] - 1;
               for (j_ptr = j_beg; j_ptr <= j_end; j_ptr++)
                  ndx[++len] = sv_ndx[j_ptr];
#else
               len = lpx_get_mat_col(lp, j-m, ndx, NULL);
#endif
            }
         }
      }
      /* return the length of the row/column pattern */
      return len;
}

static void adv_basis(glp_prob *lp)
{     int m = lpx_get_num_rows(lp);
      int n = lpx_get_num_cols(lp);
      int i, j, jj, k, size;
      int *rn, *cn, *rn_inv, *cn_inv;
      int typx, *tagx = xcalloc(1+m+n, sizeof(int));
      double lb, ub;
      xprintf("Constructing initial basis...\n");
#if 0 /* 13/V-2009 */
      if (m == 0)
         xerror("glp_adv_basis: problem has no rows\n");
      if (n == 0)
         xerror("glp_adv_basis: problem has no columns\n");
#else
      if (m == 0 || n == 0)
      {  glp_std_basis(lp);
         return;
      }
#endif
      /* use the routine triang (see above) to find maximal triangular
         part of the augmented constraint matrix A~ = (I|-A); in order
         to prevent columns of fixed variables to be included in the
         triangular part, such columns are implictly removed from the
         matrix A~ by the routine adv_mat */
      rn = xcalloc(1+m, sizeof(int));
      cn = xcalloc(1+m+n, sizeof(int));
      size = triang(m, m+n, lp, mat, rn, cn);
      if (lpx_get_int_parm(lp, LPX_K_MSGLEV) >= 3)
         xprintf("Size of triangular part = %d\n", size);
      /* the first size rows and columns of the matrix P*A~*Q (where
         P and Q are permutation matrices defined by the arrays rn and
         cn) form a lower triangular matrix; build the arrays (rn_inv
         and cn_inv), which define the matrices inv(P) and inv(Q) */
      rn_inv = xcalloc(1+m, sizeof(int));
      cn_inv = xcalloc(1+m+n, sizeof(int));
      for (i = 1; i <= m; i++) rn_inv[rn[i]] = i;
      for (j = 1; j <= m+n; j++) cn_inv[cn[j]] = j;
      /* include the columns of the matrix A~, which correspond to the
         first size columns of the matrix P*A~*Q, in the basis */
      for (k = 1; k <= m+n; k++) tagx[k] = -1;
      for (jj = 1; jj <= size; jj++)
      {  j = cn_inv[jj];
         /* the j-th column of A~ is the jj-th column of P*A~*Q */
         tagx[j] = LPX_BS;
      }
      /* if size < m, we need to add appropriate columns of auxiliary
         variables to the basis */
      for (jj = size + 1; jj <= m; jj++)
      {  /* the jj-th column of P*A~*Q should be replaced by the column
            of the auxiliary variable, for which the only unity element
            is placed in the position [jj,jj] */
         i = rn_inv[jj];
         /* the jj-th row of P*A~*Q is the i-th row of A~, but in the
            i-th row of A~ the unity element belongs to the i-th column
            of A~; therefore the disired column corresponds to the i-th
            auxiliary variable (note that this column doesn't belong to
            the triangular part found by the routine triang) */
         xassert(1 <= i && i <= m);
         xassert(cn[i] > size);
         tagx[i] = LPX_BS;
      }
      /* free working arrays */
      xfree(rn);
      xfree(cn);
      xfree(rn_inv);
      xfree(cn_inv);
      /* build tags of non-basic variables */
      for (k = 1; k <= m+n; k++)
      {  if (tagx[k] != LPX_BS)
         {  if (k <= m)
               lpx_get_row_bnds(lp, k, &typx, &lb, &ub);
            else
               lpx_get_col_bnds(lp, k-m, &typx, &lb, &ub);
            switch (typx)
            {  case LPX_FR:
                  tagx[k] = LPX_NF; break;
               case LPX_LO:
                  tagx[k] = LPX_NL; break;
               case LPX_UP:
                  tagx[k] = LPX_NU; break;
               case LPX_DB:
                  tagx[k] =
                     (fabs(lb) <= fabs(ub) ? LPX_NL : LPX_NU);
                  break;
               case LPX_FX:
                  tagx[k] = LPX_NS; break;
               default:
                  xassert(typx != typx);
            }
         }
      }
      for (k = 1; k <= m+n; k++)
      {  if (k <= m)
            lpx_set_row_stat(lp, k, tagx[k]);
         else
            lpx_set_col_stat(lp, k-m, tagx[k]);
      }
      xfree(tagx);
      return;
}

/***********************************************************************
*  NAME
*
*  glp_adv_basis - construct advanced initial LP basis
*
*  SYNOPSIS
*
*  void glp_adv_basis(glp_prob *lp, int flags);
*
*  DESCRIPTION
*
*  The routine glp_adv_basis constructs an advanced initial basis for
*  the specified problem object.
*
*  The parameter flags is reserved for use in the future and must be
*  specified as zero. */

void glp_adv_basis(glp_prob *lp, int flags)
{     if (flags != 0)
         xerror("glp_adv_basis: flags = %d; invalid flags\n", flags);
      if (lp->m == 0 || lp->n == 0)
         glp_std_basis(lp);
      else
         adv_basis(lp);
      return;
}

/* eof */
