/* glpluf.c (LU-factorization) */

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
#pragma clang diagnostic ignored "-Wsign-conversion"
#endif

#include "glpenv.h"
#include "glpluf.h"
#define xfault xerror

/* CAUTION: DO NOT CHANGE THE LIMIT BELOW */

#define N_MAX 100000000 /* = 100*10^6 */
/* maximal order of the original matrix */

/***********************************************************************
*  NAME
*
*  luf_create_it - create LU-factorization
*
*  SYNOPSIS
*
*  #include "glpluf.h"
*  LUF *luf_create_it(void);
*
*  DESCRIPTION
*
*  The routine luf_create_it creates a program object, which represents
*  LU-factorization of a square matrix.
*
*  RETURNS
*
*  The routine luf_create_it returns a pointer to the object created. */

LUF *luf_create_it(void)
{     LUF *luf;
      luf = xmalloc(sizeof(LUF));
      luf->n_max = luf->n = 0;
      luf->valid = 0;
      luf->fr_ptr = luf->fr_len = NULL;
      luf->fc_ptr = luf->fc_len = NULL;
      luf->vr_ptr = luf->vr_len = luf->vr_cap = NULL;
      luf->vr_piv = NULL;
      luf->vc_ptr = luf->vc_len = luf->vc_cap = NULL;
      luf->pp_row = luf->pp_col = NULL;
      luf->qq_row = luf->qq_col = NULL;
      luf->sv_size = 0;
      luf->sv_beg = luf->sv_end = 0;
      luf->sv_ind = NULL;
      luf->sv_val = NULL;
      luf->sv_head = luf->sv_tail = 0;
      luf->sv_prev = luf->sv_next = NULL;
      luf->vr_max = NULL;
      luf->rs_head = luf->rs_prev = luf->rs_next = NULL;
      luf->cs_head = luf->cs_prev = luf->cs_next = NULL;
      luf->flag = NULL;
      luf->work = NULL;
      luf->new_sva = 0;
      luf->piv_tol = 0.10;
      luf->piv_lim = 4;
      luf->suhl = 1;
      luf->eps_tol = 1e-15;
      luf->max_gro = 1e+10;
      luf->nnz_a = luf->nnz_f = luf->nnz_v = 0;
      luf->max_a = luf->big_v = 0.0;
      luf->rank = 0;
      return luf;
}

/***********************************************************************
*  NAME
*
*  luf_defrag_sva - defragment the sparse vector area
*
*  SYNOPSIS
*
*  #include "glpluf.h"
*  void luf_defrag_sva(LUF *luf);
*
*  DESCRIPTION
*
*  The routine luf_defrag_sva defragments the sparse vector area (SVA)
*  gathering all unused locations in one continuous extent. In order to
*  do that the routine moves all unused locations from the left part of
*  SVA (which contains rows and columns of the matrix V) to the middle
*  part (which contains free locations). This is attained by relocating
*  elements of rows and columns of the matrix V toward the beginning of
*  the left part.
*
*  NOTE that this "garbage collection" involves changing row and column
*  pointers of the matrix V. */

void luf_defrag_sva(LUF *luf)
{     int n = luf->n;
      int *vr_ptr = luf->vr_ptr;
      int *vr_len = luf->vr_len;
      int *vr_cap = luf->vr_cap;
      int *vc_ptr = luf->vc_ptr;
      int *vc_len = luf->vc_len;
      int *vc_cap = luf->vc_cap;
      int *sv_ind = luf->sv_ind;
      double *sv_val = luf->sv_val;
      int *sv_next = luf->sv_next;
      int sv_beg = 1;
      int i, j, k;
      /* skip rows and columns, which do not need to be relocated */
      for (k = luf->sv_head; k != 0; k = sv_next[k])
      {  if (k <= n)
         {  /* i-th row of the matrix V */
            i = k;
            if (vr_ptr[i] != sv_beg) break;
            vr_cap[i] = vr_len[i];
            sv_beg += vr_cap[i];
         }
         else
         {  /* j-th column of the matrix V */
            j = k - n;
            if (vc_ptr[j] != sv_beg) break;
            vc_cap[j] = vc_len[j];
            sv_beg += vc_cap[j];
         }
      }
      /* relocate other rows and columns in order to gather all unused
         locations in one continuous extent */
      for (k = k; k != 0; k = sv_next[k])
      {  if (k <= n)
         {  /* i-th row of the matrix V */
            i = k;
            memmove(&sv_ind[sv_beg], &sv_ind[vr_ptr[i]],
               vr_len[i] * sizeof(int));
            memmove(&sv_val[sv_beg], &sv_val[vr_ptr[i]],
               vr_len[i] * sizeof(double));
            vr_ptr[i] = sv_beg;
            vr_cap[i] = vr_len[i];
            sv_beg += vr_cap[i];
         }
         else
         {  /* j-th column of the matrix V */
            j = k - n;
            memmove(&sv_ind[sv_beg], &sv_ind[vc_ptr[j]],
               vc_len[j] * sizeof(int));
            memmove(&sv_val[sv_beg], &sv_val[vc_ptr[j]],
               vc_len[j] * sizeof(double));
            vc_ptr[j] = sv_beg;
            vc_cap[j] = vc_len[j];
            sv_beg += vc_cap[j];
         }
      }
      /* set new pointer to the beginning of the free part */
      luf->sv_beg = sv_beg;
      return;
}

/***********************************************************************
*  NAME
*
*  luf_enlarge_row - enlarge row capacity
*
*  SYNOPSIS
*
*  #include "glpluf.h"
*  int luf_enlarge_row(LUF *luf, int i, int cap);
*
*  DESCRIPTION
*
*  The routine luf_enlarge_row enlarges capacity of the i-th row of the
*  matrix V to cap locations (assuming that its current capacity is less
*  than cap). In order to do that the routine relocates elements of the
*  i-th row to the end of the left part of SVA (which contains rows and
*  columns of the matrix V) and then expands the left part by allocating
*  cap free locations from the free part. If there are less than cap
*  free locations, the routine defragments the sparse vector area.
*
*  Due to "garbage collection" this operation may change row and column
*  pointers of the matrix V.
*
*  RETURNS
*
*  If no error occured, the routine returns zero. Otherwise, in case of
*  overflow of the sparse vector area, the routine returns non-zero. */

int luf_enlarge_row(LUF *luf, int i, int cap)
{     int n = luf->n;
      int *vr_ptr = luf->vr_ptr;
      int *vr_len = luf->vr_len;
      int *vr_cap = luf->vr_cap;
      int *vc_cap = luf->vc_cap;
      int *sv_ind = luf->sv_ind;
      double *sv_val = luf->sv_val;
      int *sv_prev = luf->sv_prev;
      int *sv_next = luf->sv_next;
      int ret = 0;
      int cur, k, kk;
      xassert(1 <= i && i <= n);
      xassert(vr_cap[i] < cap);
      /* if there are less than cap free locations, defragment SVA */
      if (luf->sv_end - luf->sv_beg < cap)
      {  luf_defrag_sva(luf);
         if (luf->sv_end - luf->sv_beg < cap)
         {  ret = 1;
            goto done;
         }
      }
      /* save current capacity of the i-th row */
      cur = vr_cap[i];
      /* copy existing elements to the beginning of the free part */
      memmove(&sv_ind[luf->sv_beg], &sv_ind[vr_ptr[i]],
         vr_len[i] * sizeof(int));
      memmove(&sv_val[luf->sv_beg], &sv_val[vr_ptr[i]],
         vr_len[i] * sizeof(double));
      /* set new pointer and new capacity of the i-th row */
      vr_ptr[i] = luf->sv_beg;
      vr_cap[i] = cap;
      /* set new pointer to the beginning of the free part */
      luf->sv_beg += cap;
      /* now the i-th row starts in the rightmost location among other
         rows and columns of the matrix V, so its node should be moved
         to the end of the row/column linked list */
      k = i;
      /* remove the i-th row node from the linked list */
      if (sv_prev[k] == 0)
         luf->sv_head = sv_next[k];
      else
      {  /* capacity of the previous row/column can be increased at the
            expense of old locations of the i-th row */
         kk = sv_prev[k];
         if (kk <= n) vr_cap[kk] += cur; else vc_cap[kk-n] += cur;
         sv_next[sv_prev[k]] = sv_next[k];
      }
      if (sv_next[k] == 0)
         luf->sv_tail = sv_prev[k];
      else
         sv_prev[sv_next[k]] = sv_prev[k];
      /* insert the i-th row node to the end of the linked list */
      sv_prev[k] = luf->sv_tail;
      sv_next[k] = 0;
      if (sv_prev[k] == 0)
         luf->sv_head = k;
      else
         sv_next[sv_prev[k]] = k;
      luf->sv_tail = k;
done: return ret;
}

/***********************************************************************
*  NAME
*
*  luf_enlarge_col - enlarge column capacity
*
*  SYNOPSIS
*
*  #include "glpluf.h"
*  int luf_enlarge_col(LUF *luf, int j, int cap);
*
*  DESCRIPTION
*
*  The routine luf_enlarge_col enlarges capacity of the j-th column of
*  the matrix V to cap locations (assuming that its current capacity is
*  less than cap). In order to do that the routine relocates elements
*  of the j-th column to the end of the left part of SVA (which contains
*  rows and columns of the matrix V) and then expands the left part by
*  allocating cap free locations from the free part. If there are less
*  than cap free locations, the routine defragments the sparse vector
*  area.
*
*  Due to "garbage collection" this operation may change row and column
*  pointers of the matrix V.
*
*  RETURNS
*
*  If no error occured, the routine returns zero. Otherwise, in case of
*  overflow of the sparse vector area, the routine returns non-zero. */

int luf_enlarge_col(LUF *luf, int j, int cap)
{     int n = luf->n;
      int *vr_cap = luf->vr_cap;
      int *vc_ptr = luf->vc_ptr;
      int *vc_len = luf->vc_len;
      int *vc_cap = luf->vc_cap;
      int *sv_ind = luf->sv_ind;
      double *sv_val = luf->sv_val;
      int *sv_prev = luf->sv_prev;
      int *sv_next = luf->sv_next;
      int ret = 0;
      int cur, k, kk;
      xassert(1 <= j && j <= n);
      xassert(vc_cap[j] < cap);
      /* if there are less than cap free locations, defragment SVA */
      if (luf->sv_end - luf->sv_beg < cap)
      {  luf_defrag_sva(luf);
         if (luf->sv_end - luf->sv_beg < cap)
         {  ret = 1;
            goto done;
         }
      }
      /* save current capacity of the j-th column */
      cur = vc_cap[j];
      /* copy existing elements to the beginning of the free part */
      memmove(&sv_ind[luf->sv_beg], &sv_ind[vc_ptr[j]],
         vc_len[j] * sizeof(int));
      memmove(&sv_val[luf->sv_beg], &sv_val[vc_ptr[j]],
         vc_len[j] * sizeof(double));
      /* set new pointer and new capacity of the j-th column */
      vc_ptr[j] = luf->sv_beg;
      vc_cap[j] = cap;
      /* set new pointer to the beginning of the free part */
      luf->sv_beg += cap;
      /* now the j-th column starts in the rightmost location among
         other rows and columns of the matrix V, so its node should be
         moved to the end of the row/column linked list */
      k = n + j;
      /* remove the j-th column node from the linked list */
      if (sv_prev[k] == 0)
         luf->sv_head = sv_next[k];
      else
      {  /* capacity of the previous row/column can be increased at the
            expense of old locations of the j-th column */
         kk = sv_prev[k];
         if (kk <= n) vr_cap[kk] += cur; else vc_cap[kk-n] += cur;
         sv_next[sv_prev[k]] = sv_next[k];
      }
      if (sv_next[k] == 0)
         luf->sv_tail = sv_prev[k];
      else
         sv_prev[sv_next[k]] = sv_prev[k];
      /* insert the j-th column node to the end of the linked list */
      sv_prev[k] = luf->sv_tail;
      sv_next[k] = 0;
      if (sv_prev[k] == 0)
         luf->sv_head = k;
      else
         sv_next[sv_prev[k]] = k;
      luf->sv_tail = k;
done: return ret;
}

/***********************************************************************
*  reallocate - reallocate LU-factorization arrays
*
*  This routine reallocates arrays, whose size depends of n, the order
*  of the matrix A to be factorized. */

static void reallocate(LUF *luf, int n)
{     int n_max = luf->n_max;
      luf->n = n;
      if (n <= n_max) goto done;
      if (luf->fr_ptr != NULL) xfree(luf->fr_ptr);
      if (luf->fr_len != NULL) xfree(luf->fr_len);
      if (luf->fc_ptr != NULL) xfree(luf->fc_ptr);
      if (luf->fc_len != NULL) xfree(luf->fc_len);
      if (luf->vr_ptr != NULL) xfree(luf->vr_ptr);
      if (luf->vr_len != NULL) xfree(luf->vr_len);
      if (luf->vr_cap != NULL) xfree(luf->vr_cap);
      if (luf->vr_piv != NULL) xfree(luf->vr_piv);
      if (luf->vc_ptr != NULL) xfree(luf->vc_ptr);
      if (luf->vc_len != NULL) xfree(luf->vc_len);
      if (luf->vc_cap != NULL) xfree(luf->vc_cap);
      if (luf->pp_row != NULL) xfree(luf->pp_row);
      if (luf->pp_col != NULL) xfree(luf->pp_col);
      if (luf->qq_row != NULL) xfree(luf->qq_row);
      if (luf->qq_col != NULL) xfree(luf->qq_col);
      if (luf->sv_prev != NULL) xfree(luf->sv_prev);
      if (luf->sv_next != NULL) xfree(luf->sv_next);
      if (luf->vr_max != NULL) xfree(luf->vr_max);
      if (luf->rs_head != NULL) xfree(luf->rs_head);
      if (luf->rs_prev != NULL) xfree(luf->rs_prev);
      if (luf->rs_next != NULL) xfree(luf->rs_next);
      if (luf->cs_head != NULL) xfree(luf->cs_head);
      if (luf->cs_prev != NULL) xfree(luf->cs_prev);
      if (luf->cs_next != NULL) xfree(luf->cs_next);
      if (luf->flag != NULL) xfree(luf->flag);
      if (luf->work != NULL) xfree(luf->work);
      luf->n_max = n_max = n + 100;
      luf->fr_ptr = xcalloc(1+n_max, sizeof(int));
      luf->fr_len = xcalloc(1+n_max, sizeof(int));
      luf->fc_ptr = xcalloc(1+n_max, sizeof(int));
      luf->fc_len = xcalloc(1+n_max, sizeof(int));
      luf->vr_ptr = xcalloc(1+n_max, sizeof(int));
      luf->vr_len = xcalloc(1+n_max, sizeof(int));
      luf->vr_cap = xcalloc(1+n_max, sizeof(int));
      luf->vr_piv = xcalloc(1+n_max, sizeof(double));
      luf->vc_ptr = xcalloc(1+n_max, sizeof(int));
      luf->vc_len = xcalloc(1+n_max, sizeof(int));
      luf->vc_cap = xcalloc(1+n_max, sizeof(int));
      luf->pp_row = xcalloc(1+n_max, sizeof(int));
      luf->pp_col = xcalloc(1+n_max, sizeof(int));
      luf->qq_row = xcalloc(1+n_max, sizeof(int));
      luf->qq_col = xcalloc(1+n_max, sizeof(int));
      luf->sv_prev = xcalloc(1+n_max+n_max, sizeof(int));
      luf->sv_next = xcalloc(1+n_max+n_max, sizeof(int));
      luf->vr_max = xcalloc(1+n_max, sizeof(double));
      luf->rs_head = xcalloc(1+n_max, sizeof(int));
      luf->rs_prev = xcalloc(1+n_max, sizeof(int));
      luf->rs_next = xcalloc(1+n_max, sizeof(int));
      luf->cs_head = xcalloc(1+n_max, sizeof(int));
      luf->cs_prev = xcalloc(1+n_max, sizeof(int));
      luf->cs_next = xcalloc(1+n_max, sizeof(int));
      luf->flag = xcalloc(1+n_max, sizeof(int));
      luf->work = xcalloc(1+n_max, sizeof(double));
done: return;
}

/***********************************************************************
*  initialize - initialize LU-factorization data structures
*
*  This routine initializes data structures for subsequent computing
*  the LU-factorization of a given matrix A, which is specified by the
*  formal routine col. On exit V = A and F = P = Q = I, where I is the
*  unity matrix. (Row-wise representation of the matrix F is not used
*  at the factorization stage and therefore is not initialized.)
*
*  If no error occured, the routine returns zero. Otherwise, in case of
*  overflow of the sparse vector area, the routine returns non-zero. */

static int initialize(LUF *luf, int (*col)(void *info, int j, int rn[],
      double aj[]), void *info)
{     int n = luf->n;
      int *fc_ptr = luf->fc_ptr;
      int *fc_len = luf->fc_len;
      int *vr_ptr = luf->vr_ptr;
      int *vr_len = luf->vr_len;
      int *vr_cap = luf->vr_cap;
      int *vc_ptr = luf->vc_ptr;
      int *vc_len = luf->vc_len;
      int *vc_cap = luf->vc_cap;
      int *pp_row = luf->pp_row;
      int *pp_col = luf->pp_col;
      int *qq_row = luf->qq_row;
      int *qq_col = luf->qq_col;
      int *sv_ind = luf->sv_ind;
      double *sv_val = luf->sv_val;
      int *sv_prev = luf->sv_prev;
      int *sv_next = luf->sv_next;
      double *vr_max = luf->vr_max;
      int *rs_head = luf->rs_head;
      int *rs_prev = luf->rs_prev;
      int *rs_next = luf->rs_next;
      int *cs_head = luf->cs_head;
      int *cs_prev = luf->cs_prev;
      int *cs_next = luf->cs_next;
      int *flag = luf->flag;
      double *work = luf->work;
      int ret = 0;
      int i, i_ptr, j, j_beg, j_end, k, len, nnz, sv_beg, sv_end, ptr;
      double big, val;
      /* free all locations of the sparse vector area */
      sv_beg = 1;
      sv_end = luf->sv_size + 1;
      /* (row-wise representation of the matrix F is not initialized,
         because it is not used at the factorization stage) */
      /* build the matrix F in column-wise format (initially F = I) */
      for (j = 1; j <= n; j++)
      {  fc_ptr[j] = sv_end;
         fc_len[j] = 0;
      }
      /* clear rows of the matrix V; clear the flag array */
      for (i = 1; i <= n; i++)
         vr_len[i] = vr_cap[i] = 0, flag[i] = 0;
      /* build the matrix V in column-wise format (initially V = A);
         count non-zeros in rows of this matrix; count total number of
         non-zeros; compute largest of absolute values of elements */
      nnz = 0;
      big = 0.0;
      for (j = 1; j <= n; j++)
      {  int *rn = pp_row;
         double *aj = work;
         /* obtain j-th column of the matrix A */
         len = col(info, j, rn, aj);
         if (!(0 <= len && len <= n))
            xfault("luf_factorize: j = %d; len = %d; invalid column len"
               "gth\n", j, len);
         /* check for free locations */
         if (sv_end - sv_beg < len)
         {  /* overflow of the sparse vector area */
            ret = 1;
            goto done;
         }
         /* set pointer to the j-th column */
         vc_ptr[j] = sv_beg;
         /* set length of the j-th column */
         vc_len[j] = vc_cap[j] = len;
         /* count total number of non-zeros */
         nnz += len;
         /* walk through elements of the j-th column */
         for (ptr = 1; ptr <= len; ptr++)
         {  /* get row index and numerical value of a[i,j] */
            i = rn[ptr];
            val = aj[ptr];
            if (!(1 <= i && i <= n))
               xfault("luf_factorize: i = %d; j = %d; invalid row index"
                  "\n", i, j);
            if (flag[i])
               xfault("luf_factorize: i = %d; j = %d; duplicate element"
                  " not allowed\n", i, j);
            if (val == 0.0)
               xfault("luf_factorize: i = %d; j = %d; zero element not "
                  "allowed\n", i, j);
            /* add new element v[i,j] = a[i,j] to j-th column */
            sv_ind[sv_beg] = i;
            sv_val[sv_beg] = val;
            sv_beg++;
            /* big := max(big, |a[i,j]|) */
            if (val < 0.0) val = - val;
            if (big < val) big = val;
            /* mark non-zero in the i-th position of the j-th column */
            flag[i] = 1;
            /* increase length of the i-th row */
            vr_cap[i]++;
         }
         /* reset all non-zero marks */
         for (ptr = 1; ptr <= len; ptr++) flag[rn[ptr]] = 0;
      }
      /* allocate rows of the matrix V */
      for (i = 1; i <= n; i++)
      {  /* get length of the i-th row */
         len = vr_cap[i];
         /* check for free locations */
         if (sv_end - sv_beg < len)
         {  /* overflow of the sparse vector area */
            ret = 1;
            goto done;
         }
         /* set pointer to the i-th row */
         vr_ptr[i] = sv_beg;
         /* reserve locations for the i-th row */
         sv_beg += len;
      }
      /* build the matrix V in row-wise format using representation of
         this matrix in column-wise format */
      for (j = 1; j <= n; j++)
      {  /* walk through elements of the j-th column */
         j_beg = vc_ptr[j];
         j_end = j_beg + vc_len[j] - 1;
         for (k = j_beg; k <= j_end; k++)
         {  /* get row index and numerical value of v[i,j] */
            i = sv_ind[k];
            val = sv_val[k];
            /* store element in the i-th row */
            i_ptr = vr_ptr[i] + vr_len[i];
            sv_ind[i_ptr] = j;
            sv_val[i_ptr] = val;
            /* increase count of the i-th row */
            vr_len[i]++;
         }
      }
      /* initialize the matrices P and Q (initially P = Q = I) */
      for (k = 1; k <= n; k++)
         pp_row[k] = pp_col[k] = qq_row[k] = qq_col[k] = k;
      /* set sva partitioning pointers */
      luf->sv_beg = sv_beg;
      luf->sv_end = sv_end;
      /* the initial physical order of rows and columns of the matrix V
         is n+1, ..., n+n, 1, ..., n (firstly columns, then rows) */
      luf->sv_head = n+1;
      luf->sv_tail = n;
      for (i = 1; i <= n; i++)
      {  sv_prev[i] = i-1;
         sv_next[i] = i+1;
      }
      sv_prev[1] = n+n;
      sv_next[n] = 0;
      for (j = 1; j <= n; j++)
      {  sv_prev[n+j] = n+j-1;
         sv_next[n+j] = n+j+1;
      }
      sv_prev[n+1] = 0;
      sv_next[n+n] = 1;
      /* clear working arrays */
      for (k = 1; k <= n; k++)
      {  flag[k] = 0;
         work[k] = 0.0;
      }
      /* initialize some statistics */
      luf->nnz_a = nnz;
      luf->nnz_f = 0;
      luf->nnz_v = nnz;
      luf->max_a = big;
      luf->big_v = big;
      luf->rank = -1;
      /* initially the active submatrix is the entire matrix V */
      /* largest of absolute values of elements in each active row is
         unknown yet */
      for (i = 1; i <= n; i++) vr_max[i] = -1.0;
      /* build linked lists of active rows */
      for (len = 0; len <= n; len++) rs_head[len] = 0;
      for (i = 1; i <= n; i++)
      {  len = vr_len[i];
         rs_prev[i] = 0;
         rs_next[i] = rs_head[len];
         if (rs_next[i] != 0) rs_prev[rs_next[i]] = i;
         rs_head[len] = i;
      }
      /* build linked lists of active columns */
      for (len = 0; len <= n; len++) cs_head[len] = 0;
      for (j = 1; j <= n; j++)
      {  len = vc_len[j];
         cs_prev[j] = 0;
         cs_next[j] = cs_head[len];
         if (cs_next[j] != 0) cs_prev[cs_next[j]] = j;
         cs_head[len] = j;
      }
done: /* return to the factorizing routine */
      return ret;
}

/***********************************************************************
*  find_pivot - choose a pivot element
*
*  This routine chooses a pivot element in the active submatrix of the
*  matrix U = P*V*Q.
*
*  It is assumed that on entry the matrix U has the following partially
*  triangularized form:
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
*  (elements of the active submatrix are marked by '*').
*
*  Since the matrix U = P*V*Q is not stored, the routine works with the
*  matrix V. It is assumed that the row-wise representation corresponds
*  to the matrix V, but the column-wise representation corresponds to
*  the active submatrix of the matrix V, i.e. elements of the matrix V,
*  which doesn't belong to the active submatrix, are missing from the
*  column linked lists. It is also assumed that each active row of the
*  matrix V is in the set R[len], where len is number of non-zeros in
*  the row, and each active column of the matrix V is in the set C[len],
*  where len is number of non-zeros in the column (in the latter case
*  only elements of the active submatrix are counted; such elements are
*  marked by '*' on the figure above).
* 
*  For the reason of numerical stability the routine applies so called
*  threshold pivoting proposed by J.Reid. It is assumed that an element
*  v[i,j] can be selected as a pivot candidate if it is not very small
*  (in absolute value) among other elements in the same row, i.e. if it
*  satisfies to the stability condition |v[i,j]| >= tol * max|v[i,*]|,
*  where 0 < tol < 1 is a given tolerance.
* 
*  In order to keep sparsity of the matrix V the routine uses Markowitz
*  strategy, trying to choose such element v[p,q], which satisfies to
*  the stability condition (see above) and has smallest Markowitz cost
*  (nr[p]-1) * (nc[q]-1), where nr[p] and nc[q] are numbers of non-zero
*  elements, respectively, in the p-th row and in the q-th column of the
*  active submatrix.
* 
*  In order to reduce the search, i.e. not to walk through all elements
*  of the active submatrix, the routine exploits a technique proposed by
*  I.Duff. This technique is based on using the sets R[len] and C[len]
*  of active rows and columns.
* 
*  If the pivot element v[p,q] has been chosen, the routine stores its
*  indices to the locations *p and *q and returns zero. Otherwise, if
*  the active submatrix is empty and therefore the pivot element can't
*  be chosen, the routine returns non-zero. */

static int find_pivot(LUF *luf, int *_p, int *_q)
{     int n = luf->n;
      int *vr_ptr = luf->vr_ptr;
      int *vr_len = luf->vr_len;
      int *vc_ptr = luf->vc_ptr;
      int *vc_len = luf->vc_len;
      int *sv_ind = luf->sv_ind;
      double *sv_val = luf->sv_val;
      double *vr_max = luf->vr_max;
      int *rs_head = luf->rs_head;
      int *rs_next = luf->rs_next;
      int *cs_head = luf->cs_head;
      int *cs_prev = luf->cs_prev;
      int *cs_next = luf->cs_next;
      double piv_tol = luf->piv_tol;
      int piv_lim = luf->piv_lim;
      int suhl = luf->suhl;
      int p, q, len, i, i_beg, i_end, i_ptr, j, j_beg, j_end, j_ptr,
         ncand, next_j, min_p, min_q, min_len;
      double best, cost, big, temp;
      /* initially no pivot candidates have been found so far */
      p = q = 0, best = DBL_MAX, ncand = 0;
      /* if in the active submatrix there is a column that has the only
         non-zero (column singleton), choose it as pivot */
      j = cs_head[1];
      if (j != 0)
      {  xassert(vc_len[j] == 1);
         p = sv_ind[vc_ptr[j]], q = j;
         goto done;
      }
      /* if in the active submatrix there is a row that has the only
         non-zero (row singleton), choose it as pivot */
      i = rs_head[1];
      if (i != 0)
      {  xassert(vr_len[i] == 1);
         p = i, q = sv_ind[vr_ptr[i]];
         goto done;
      }
      /* there are no singletons in the active submatrix; walk through
         other non-empty rows and columns */
      for (len = 2; len <= n; len++)
      {  /* consider active columns that have len non-zeros */
         for (j = cs_head[len]; j != 0; j = next_j)
         {  /* the j-th column has len non-zeros */
            j_beg = vc_ptr[j];
            j_end = j_beg + vc_len[j] - 1;
            /* save pointer to the next column with the same length */
            next_j = cs_next[j];
            /* find an element in the j-th column, which is placed in a
               row with minimal number of non-zeros and satisfies to the
               stability condition (such element may not exist) */
            min_p = min_q = 0, min_len = INT_MAX;
            for (j_ptr = j_beg; j_ptr <= j_end; j_ptr++)
            {  /* get row index of v[i,j] */
               i = sv_ind[j_ptr];
               i_beg = vr_ptr[i];
               i_end = i_beg + vr_len[i] - 1;
               /* if the i-th row is not shorter than that one, where
                  minimal element is currently placed, skip v[i,j] */
               if (vr_len[i] >= min_len) continue;
               /* determine the largest of absolute values of elements
                  in the i-th row */
               big = vr_max[i];
               if (big < 0.0)
               {  /* the largest value is unknown yet; compute it */
                  for (i_ptr = i_beg; i_ptr <= i_end; i_ptr++)
                  {  temp = sv_val[i_ptr];
                     if (temp < 0.0) temp = - temp;
                     if (big < temp) big = temp;
                  }
                  vr_max[i] = big;
               }
               /* find v[i,j] in the i-th row */
               for (i_ptr = vr_ptr[i]; sv_ind[i_ptr] != j; i_ptr++);
               xassert(i_ptr <= i_end);
               /* if v[i,j] doesn't satisfy to the stability condition,
                  skip it */
               temp = sv_val[i_ptr];
               if (temp < 0.0) temp = - temp;
               if (temp < piv_tol * big) continue;
               /* v[i,j] is better than the current minimal element */
               min_p = i, min_q = j, min_len = vr_len[i];
               /* if Markowitz cost of the current minimal element is
                  not greater than (len-1)**2, it can be chosen right
                  now; this heuristic reduces the search and works well
                  in many cases */
               if (min_len <= len)
               {  p = min_p, q = min_q;
                  goto done;
               }
            }
            /* the j-th column has been scanned */
            if (min_p != 0)
            {  /* the minimal element is a next pivot candidate */
               ncand++;
               /* compute its Markowitz cost */
               cost = (double)(min_len - 1) * (double)(len - 1);
               /* choose between the minimal element and the current
                  candidate */
               if (cost < best) p = min_p, q = min_q, best = cost;
               /* if piv_lim candidates have been considered, there are
                  doubts that a much better candidate exists; therefore
                  it's time to terminate the search */
               if (ncand == piv_lim) goto done;
            }
            else
            {  /* the j-th column has no elements, which satisfy to the
                  stability condition; Uwe Suhl suggests to exclude such
                  column from the further consideration until it becomes
                  a column singleton; in hard cases this significantly
                  reduces a time needed for pivot searching */
               if (suhl)
               {  /* remove the j-th column from the active set */
                  if (cs_prev[j] == 0)
                     cs_head[len] = cs_next[j];
                  else
                     cs_next[cs_prev[j]] = cs_next[j];
                  if (cs_next[j] == 0)
                     /* nop */;
                  else
                     cs_prev[cs_next[j]] = cs_prev[j];
                  /* the following assignment is used to avoid an error
                     when the routine eliminate (see below) will try to
                     remove the j-th column from the active set */
                  cs_prev[j] = cs_next[j] = j;
               }
            }
         }
         /* consider active rows that have len non-zeros */
         for (i = rs_head[len]; i != 0; i = rs_next[i])
         {  /* the i-th row has len non-zeros */
            i_beg = vr_ptr[i];
            i_end = i_beg + vr_len[i] - 1;
            /* determine the largest of absolute values of elements in
               the i-th row */
            big = vr_max[i];
            if (big < 0.0)
            {  /* the largest value is unknown yet; compute it */
               for (i_ptr = i_beg; i_ptr <= i_end; i_ptr++)
               {  temp = sv_val[i_ptr];
                  if (temp < 0.0) temp = - temp;
                  if (big < temp) big = temp;
               }
               vr_max[i] = big;
            }
            /* find an element in the i-th row, which is placed in a
               column with minimal number of non-zeros and satisfies to
               the stability condition (such element always exists) */
            min_p = min_q = 0, min_len = INT_MAX;
            for (i_ptr = i_beg; i_ptr <= i_end; i_ptr++)
            {  /* get column index of v[i,j] */
               j = sv_ind[i_ptr];
               /* if the j-th column is not shorter than that one, where
                  minimal element is currently placed, skip v[i,j] */
               if (vc_len[j] >= min_len) continue;
               /* if v[i,j] doesn't satisfy to the stability condition,
                  skip it */
               temp = sv_val[i_ptr];
               if (temp < 0.0) temp = - temp;
               if (temp < piv_tol * big) continue;
               /* v[i,j] is better than the current minimal element */
               min_p = i, min_q = j, min_len = vc_len[j];
               /* if Markowitz cost of the current minimal element is
                  not greater than (len-1)**2, it can be chosen right
                  now; this heuristic reduces the search and works well
                  in many cases */
               if (min_len <= len)
               {  p = min_p, q = min_q;
                  goto done;
               }
            }
            /* the i-th row has been scanned */
            if (min_p != 0)
            {  /* the minimal element is a next pivot candidate */
               ncand++;
               /* compute its Markowitz cost */
               cost = (double)(len - 1) * (double)(min_len - 1);
               /* choose between the minimal element and the current
                  candidate */
               if (cost < best) p = min_p, q = min_q, best = cost;
               /* if piv_lim candidates have been considered, there are
                  doubts that a much better candidate exists; therefore
                  it's time to terminate the search */
               if (ncand == piv_lim) goto done;
            }
            else
            {  /* this can't be because this can never be */
               xassert(min_p != min_p);
            }
         }
      }
done: /* bring the pivot to the factorizing routine */
      *_p = p, *_q = q;
      return (p == 0);
}

/***********************************************************************
*  eliminate - perform gaussian elimination.
* 
*  This routine performs elementary gaussian transformations in order
*  to eliminate subdiagonal elements in the k-th column of the matrix
*  U = P*V*Q using the pivot element u[k,k], where k is the number of
*  the current elimination step.
* 
*  The parameters p and q are, respectively, row and column indices of
*  the element v[p,q], which corresponds to the element u[k,k].
* 
*  Each time when the routine applies the elementary transformation to
*  a non-pivot row of the matrix V, it stores the corresponding element
*  to the matrix F in order to keep the main equality A = F*V.
* 
*  The routine assumes that on entry the matrices L = P*F*inv(P) and
*  U = P*V*Q are the following:
* 
*        1       k                  1       k         n
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
*  where rows and columns of the matrix U with numbers k, k+1, ..., n
*  form the active submatrix (eliminated elements are marked by '#' and
*  other elements of the active submatrix are marked by '*'). Note that
*  each eliminated non-zero element u[i,k] of the matrix U gives the
*  corresponding element l[i,k] of the matrix L (marked by '_').
* 
*  Actually all operations are performed on the matrix V. Should note
*  that the row-wise representation corresponds to the matrix V, but the
*  column-wise representation corresponds to the active submatrix of the
*  matrix V, i.e. elements of the matrix V, which doesn't belong to the
*  active submatrix, are missing from the column linked lists.
* 
*  Let u[k,k] = v[p,q] be the pivot. In order to eliminate subdiagonal
*  elements u[i',k] = v[i,q], i' = k+1, k+2, ..., n, the routine applies
*  the following elementary gaussian transformations:
* 
*     (i-th row of V) := (i-th row of V) - f[i,p] * (p-th row of V),
* 
*  where f[i,p] = v[i,q] / v[p,q] is a gaussian multiplier.
*
*  Additionally, in order to keep the main equality A = F*V, each time
*  when the routine applies the transformation to i-th row of the matrix
*  V, it also adds f[i,p] as a new element to the matrix F.
*
*  IMPORTANT: On entry the working arrays flag and work should contain
*  zeros. This status is provided by the routine on exit.
*
*  If no error occured, the routine returns zero. Otherwise, in case of
*  overflow of the sparse vector area, the routine returns non-zero. */

static int eliminate(LUF *luf, int p, int q)
{     int n = luf->n;
      int *fc_ptr = luf->fc_ptr;
      int *fc_len = luf->fc_len;
      int *vr_ptr = luf->vr_ptr;
      int *vr_len = luf->vr_len;
      int *vr_cap = luf->vr_cap;
      double *vr_piv = luf->vr_piv;
      int *vc_ptr = luf->vc_ptr;
      int *vc_len = luf->vc_len;
      int *vc_cap = luf->vc_cap;
      int *sv_ind = luf->sv_ind;
      double *sv_val = luf->sv_val;
      int *sv_prev = luf->sv_prev;
      int *sv_next = luf->sv_next;
      double *vr_max = luf->vr_max;
      int *rs_head = luf->rs_head;
      int *rs_prev = luf->rs_prev;
      int *rs_next = luf->rs_next;
      int *cs_head = luf->cs_head;
      int *cs_prev = luf->cs_prev;
      int *cs_next = luf->cs_next;
      int *flag = luf->flag;
      double *work = luf->work;
      double eps_tol = luf->eps_tol;
      /* at this stage the row-wise representation of the matrix F is
         not used, so fr_len can be used as a working array */
      int *ndx = luf->fr_len;
      int ret = 0;
      int len, fill, i, i_beg, i_end, i_ptr, j, j_beg, j_end, j_ptr, k,
         p_beg, p_end, p_ptr, q_beg, q_end, q_ptr;
      double fip, val, vpq, temp;
      xassert(1 <= p && p <= n);
      xassert(1 <= q && q <= n);
      /* remove the p-th (pivot) row from the active set; this row will
         never return there */
      if (rs_prev[p] == 0)
         rs_head[vr_len[p]] = rs_next[p];
      else
         rs_next[rs_prev[p]] = rs_next[p];
      if (rs_next[p] == 0)
         ;
      else
         rs_prev[rs_next[p]] = rs_prev[p];
      /* remove the q-th (pivot) column from the active set; this column
         will never return there */
      if (cs_prev[q] == 0)
         cs_head[vc_len[q]] = cs_next[q];
      else
         cs_next[cs_prev[q]] = cs_next[q];
      if (cs_next[q] == 0)
         ;
      else
         cs_prev[cs_next[q]] = cs_prev[q];
      /* find the pivot v[p,q] = u[k,k] in the p-th row */
      p_beg = vr_ptr[p];
      p_end = p_beg + vr_len[p] - 1;
      for (p_ptr = p_beg; sv_ind[p_ptr] != q; p_ptr++) /* nop */;
      xassert(p_ptr <= p_end);
      /* store value of the pivot */
      vpq = (vr_piv[p] = sv_val[p_ptr]);
      /* remove the pivot from the p-th row */
      sv_ind[p_ptr] = sv_ind[p_end];
      sv_val[p_ptr] = sv_val[p_end];
      vr_len[p]--;
      p_end--;
      /* find the pivot v[p,q] = u[k,k] in the q-th column */
      q_beg = vc_ptr[q];
      q_end = q_beg + vc_len[q] - 1;
      for (q_ptr = q_beg; sv_ind[q_ptr] != p; q_ptr++) /* nop */;
      xassert(q_ptr <= q_end);
      /* remove the pivot from the q-th column */
      sv_ind[q_ptr] = sv_ind[q_end];
      vc_len[q]--;
      q_end--;
      /* walk through the p-th (pivot) row, which doesn't contain the
         pivot v[p,q] already, and do the following... */
      for (p_ptr = p_beg; p_ptr <= p_end; p_ptr++)
      {  /* get column index of v[p,j] */
         j = sv_ind[p_ptr];
         /* store v[p,j] to the working array */
         flag[j] = 1;
         work[j] = sv_val[p_ptr];
         /* remove the j-th column from the active set; this column will
            return there later with new length */
         if (cs_prev[j] == 0)
            cs_head[vc_len[j]] = cs_next[j];
         else
            cs_next[cs_prev[j]] = cs_next[j];
         if (cs_next[j] == 0)
            ;
         else
            cs_prev[cs_next[j]] = cs_prev[j];
         /* find v[p,j] in the j-th column */
         j_beg = vc_ptr[j];
         j_end = j_beg + vc_len[j] - 1;
         for (j_ptr = j_beg; sv_ind[j_ptr] != p; j_ptr++) /* nop */;
         xassert(j_ptr <= j_end);
         /* since v[p,j] leaves the active submatrix, remove it from the
            j-th column; however, v[p,j] is kept in the p-th row */
         sv_ind[j_ptr] = sv_ind[j_end];
         vc_len[j]--;
      }
      /* walk through the q-th (pivot) column, which doesn't contain the
         pivot v[p,q] already, and perform gaussian elimination */
      while (q_beg <= q_end)
      {  /* element v[i,q] should be eliminated */
         /* get row index of v[i,q] */
         i = sv_ind[q_beg];
         /* remove the i-th row from the active set; later this row will
            return there with new length */
         if (rs_prev[i] == 0)
            rs_head[vr_len[i]] = rs_next[i];
         else
            rs_next[rs_prev[i]] = rs_next[i];
         if (rs_next[i] == 0)
            ;
         else
            rs_prev[rs_next[i]] = rs_prev[i];
         /* find v[i,q] in the i-th row */
         i_beg = vr_ptr[i];
         i_end = i_beg + vr_len[i] - 1;
         for (i_ptr = i_beg; sv_ind[i_ptr] != q; i_ptr++) /* nop */;
         xassert(i_ptr <= i_end);
         /* compute gaussian multiplier f[i,p] = v[i,q] / v[p,q] */
         fip = sv_val[i_ptr] / vpq;
         /* since v[i,q] should be eliminated, remove it from the i-th
            row */
         sv_ind[i_ptr] = sv_ind[i_end];
         sv_val[i_ptr] = sv_val[i_end];
         vr_len[i]--;
         i_end--;
         /* and from the q-th column */
         sv_ind[q_beg] = sv_ind[q_end];
         vc_len[q]--;
         q_end--;
         /* perform gaussian transformation:
            (i-th row) := (i-th row) - f[i,p] * (p-th row)
            note that now the p-th row, which is in the working array,
            doesn't contain the pivot v[p,q], and the i-th row doesn't
            contain the eliminated element v[i,q] */
         /* walk through the i-th row and transform existing non-zero
            elements */
         fill = vr_len[p];
         for (i_ptr = i_beg; i_ptr <= i_end; i_ptr++)
         {  /* get column index of v[i,j] */
            j = sv_ind[i_ptr];
            /* v[i,j] := v[i,j] - f[i,p] * v[p,j] */
            if (flag[j])
            {  /* v[p,j] != 0 */
               temp = (sv_val[i_ptr] -= fip * work[j]);
               if (temp < 0.0) temp = - temp;
               flag[j] = 0;
               fill--; /* since both v[i,j] and v[p,j] exist */
               if (temp == 0.0 || temp < eps_tol)
               {  /* new v[i,j] is closer to zero; replace it by exact
                     zero, i.e. remove it from the active submatrix */
                  /* remove v[i,j] from the i-th row */
                  sv_ind[i_ptr] = sv_ind[i_end];
                  sv_val[i_ptr] = sv_val[i_end];
                  vr_len[i]--;
                  i_ptr--;
                  i_end--;
                  /* find v[i,j] in the j-th column */
                  j_beg = vc_ptr[j];
                  j_end = j_beg + vc_len[j] - 1;
                  for (j_ptr = j_beg; sv_ind[j_ptr] != i; j_ptr++);
                  xassert(j_ptr <= j_end);
                  /* remove v[i,j] from the j-th column */
                  sv_ind[j_ptr] = sv_ind[j_end];
                  vc_len[j]--;
               }
               else
               {  /* v_big := max(v_big, |v[i,j]|) */
                  if (luf->big_v < temp) luf->big_v = temp;
               }
            }
         }
         /* now flag is the pattern of the set v[p,*] \ v[i,*], and fill
            is number of non-zeros in this set; therefore up to fill new
            non-zeros may appear in the i-th row */
         if (vr_len[i] + fill > vr_cap[i])
         {  /* enlarge the i-th row */
            if (luf_enlarge_row(luf, i, vr_len[i] + fill))
            {  /* overflow of the sparse vector area */
               ret = 1;
               goto done;
            }
            /* defragmentation may change row and column pointers of the
               matrix V */
            p_beg = vr_ptr[p];
            p_end = p_beg + vr_len[p] - 1;
            q_beg = vc_ptr[q];
            q_end = q_beg + vc_len[q] - 1;
         }
         /* walk through the p-th (pivot) row and create new elements
            of the i-th row that appear due to fill-in; column indices
            of these new elements are accumulated in the array ndx */
         len = 0;
         for (p_ptr = p_beg; p_ptr <= p_end; p_ptr++)
         {  /* get column index of v[p,j], which may cause fill-in */
            j = sv_ind[p_ptr];
            if (flag[j])
            {  /* compute new non-zero v[i,j] = 0 - f[i,p] * v[p,j] */
               temp = (val = - fip * work[j]);
               if (temp < 0.0) temp = - temp;
               if (temp == 0.0 || temp < eps_tol)
                  /* if v[i,j] is closer to zero; just ignore it */;
               else
               {  /* add v[i,j] to the i-th row */
                  i_ptr = vr_ptr[i] + vr_len[i];
                  sv_ind[i_ptr] = j;
                  sv_val[i_ptr] = val;
                  vr_len[i]++;
                  /* remember column index of v[i,j] */
                  ndx[++len] = j;
                  /* big_v := max(big_v, |v[i,j]|) */
                  if (luf->big_v < temp) luf->big_v = temp;
               }
            }
            else
            {  /* there is no fill-in, because v[i,j] already exists in
                  the i-th row; restore the flag of the element v[p,j],
                  which was reset before */
               flag[j] = 1;
            }
         }
         /* add new non-zeros v[i,j] to the corresponding columns */
         for (k = 1; k <= len; k++)
         {  /* get column index of new non-zero v[i,j] */
            j = ndx[k];
            /* one free location is needed in the j-th column */
            if (vc_len[j] + 1 > vc_cap[j])
            {  /* enlarge the j-th column */
               if (luf_enlarge_col(luf, j, vc_len[j] + 10))
               {  /* overflow of the sparse vector area */
                  ret = 1;
                  goto done;
               }
               /* defragmentation may change row and column pointers of
                  the matrix V */
               p_beg = vr_ptr[p];
               p_end = p_beg + vr_len[p] - 1;
               q_beg = vc_ptr[q];
               q_end = q_beg + vc_len[q] - 1;
            }
            /* add new non-zero v[i,j] to the j-th column */
            j_ptr = vc_ptr[j] + vc_len[j];
            sv_ind[j_ptr] = i;
            vc_len[j]++;
         }
         /* now the i-th row has been completely transformed, therefore
            it can return to the active set with new length */
         rs_prev[i] = 0;
         rs_next[i] = rs_head[vr_len[i]];
         if (rs_next[i] != 0) rs_prev[rs_next[i]] = i;
         rs_head[vr_len[i]] = i;
         /* the largest of absolute values of elements in the i-th row
            is currently unknown */
         vr_max[i] = -1.0;
         /* at least one free location is needed to store the gaussian
            multiplier */
         if (luf->sv_end - luf->sv_beg < 1)
         {  /* there are no free locations at all; defragment SVA */
            luf_defrag_sva(luf);
            if (luf->sv_end - luf->sv_beg < 1)
            {  /* overflow of the sparse vector area */
               ret = 1;
               goto done;
            }
            /* defragmentation may change row and column pointers of the
               matrix V */
            p_beg = vr_ptr[p];
            p_end = p_beg + vr_len[p] - 1;
            q_beg = vc_ptr[q];
            q_end = q_beg + vc_len[q] - 1;
         }
         /* add the element f[i,p], which is the gaussian multiplier,
            to the matrix F */
         luf->sv_end--;
         sv_ind[luf->sv_end] = i;
         sv_val[luf->sv_end] = fip;
         fc_len[p]++;
         /* end of elimination loop */
      }
      /* at this point the q-th (pivot) column should be empty */
      xassert(vc_len[q] == 0);
      /* reset capacity of the q-th column */
      vc_cap[q] = 0;
      /* remove node of the q-th column from the addressing list */
      k = n + q;
      if (sv_prev[k] == 0)
         luf->sv_head = sv_next[k];
      else
         sv_next[sv_prev[k]] = sv_next[k];
      if (sv_next[k] == 0)
         luf->sv_tail = sv_prev[k];
      else
         sv_prev[sv_next[k]] = sv_prev[k];
      /* the p-th column of the matrix F has been completely built; set
         its pointer */
      fc_ptr[p] = luf->sv_end;
      /* walk through the p-th (pivot) row and do the following... */
      for (p_ptr = p_beg; p_ptr <= p_end; p_ptr++)
      {  /* get column index of v[p,j] */
         j = sv_ind[p_ptr];
         /* erase v[p,j] from the working array */
         flag[j] = 0;
         work[j] = 0.0;
         /* the j-th column has been completely transformed, therefore
            it can return to the active set with new length; however
            the special case c_prev[j] = c_next[j] = j means that the
            routine find_pivot excluded the j-th column from the active
            set due to Uwe Suhl's rule, and therefore in this case the
            column can return to the active set only if it is a column
            singleton */
         if (!(vc_len[j] != 1 && cs_prev[j] == j && cs_next[j] == j))
         {  cs_prev[j] = 0;
            cs_next[j] = cs_head[vc_len[j]];
            if (cs_next[j] != 0) cs_prev[cs_next[j]] = j;
            cs_head[vc_len[j]] = j;
         }
      }
done: /* return to the factorizing routine */
      return ret;
}

/***********************************************************************
*  build_v_cols - build the matrix V in column-wise format
*
*  This routine builds the column-wise representation of the matrix V
*  using its row-wise representation.
*
*  If no error occured, the routine returns zero. Otherwise, in case of
*  overflow of the sparse vector area, the routine returns non-zero. */

static int build_v_cols(LUF *luf)
{     int n = luf->n;
      int *vr_ptr = luf->vr_ptr;
      int *vr_len = luf->vr_len;
      int *vc_ptr = luf->vc_ptr;
      int *vc_len = luf->vc_len;
      int *vc_cap = luf->vc_cap;
      int *sv_ind = luf->sv_ind;
      double *sv_val = luf->sv_val;
      int *sv_prev = luf->sv_prev;
      int *sv_next = luf->sv_next;
      int ret = 0;
      int i, i_beg, i_end, i_ptr, j, j_ptr, k, nnz;
      /* it is assumed that on entry all columns of the matrix V are
         empty, i.e. vc_len[j] = vc_cap[j] = 0 for all j = 1, ..., n,
         and have been removed from the addressing list */
      /* count non-zeros in columns of the matrix V; count total number
         of non-zeros in this matrix */
      nnz = 0;
      for (i = 1; i <= n; i++)
      {  /* walk through elements of the i-th row and count non-zeros
            in the corresponding columns */
         i_beg = vr_ptr[i];
         i_end = i_beg + vr_len[i] - 1;
         for (i_ptr = i_beg; i_ptr <= i_end; i_ptr++)
            vc_cap[sv_ind[i_ptr]]++;
         /* count total number of non-zeros */
         nnz += vr_len[i];
      }
      /* store total number of non-zeros */
      luf->nnz_v = nnz;
      /* check for free locations */
      if (luf->sv_end - luf->sv_beg < nnz)
      {  /* overflow of the sparse vector area */
         ret = 1;
         goto done;
      }
      /* allocate columns of the matrix V */
      for (j = 1; j <= n; j++)
      {  /* set pointer to the j-th column */
         vc_ptr[j] = luf->sv_beg;
         /* reserve locations for the j-th column */
         luf->sv_beg += vc_cap[j];
      }
      /* build the matrix V in column-wise format using this matrix in
         row-wise format */
      for (i = 1; i <= n; i++)
      {  /* walk through elements of the i-th row */
         i_beg = vr_ptr[i];
         i_end = i_beg + vr_len[i] - 1;
         for (i_ptr = i_beg; i_ptr <= i_end; i_ptr++)
         {  /* get column index */
            j = sv_ind[i_ptr];
            /* store element in the j-th column */
            j_ptr = vc_ptr[j] + vc_len[j];
            sv_ind[j_ptr] = i;
            sv_val[j_ptr] = sv_val[i_ptr];
            /* increase length of the j-th column */
            vc_len[j]++;
         }
      }
      /* now columns are placed in the sparse vector area behind rows
         in the order n+1, n+2, ..., n+n; so insert column nodes in the
         addressing list using this order */
      for (k = n+1; k <= n+n; k++)
      {  sv_prev[k] = k-1;
         sv_next[k] = k+1;
      }
      sv_prev[n+1] = luf->sv_tail;
      sv_next[luf->sv_tail] = n+1;
      sv_next[n+n] = 0;
      luf->sv_tail = n+n;
done: /* return to the factorizing routine */
      return ret;
}

/***********************************************************************
*  build_f_rows - build the matrix F in row-wise format
*
*  This routine builds the row-wise representation of the matrix F using
*  its column-wise representation.
*
*  If no error occured, the routine returns zero. Otherwise, in case of
*  overflow of the sparse vector area, the routine returns non-zero. */

static int build_f_rows(LUF *luf)
{     int n = luf->n;
      int *fr_ptr = luf->fr_ptr;
      int *fr_len = luf->fr_len;
      int *fc_ptr = luf->fc_ptr;
      int *fc_len = luf->fc_len;
      int *sv_ind = luf->sv_ind;
      double *sv_val = luf->sv_val;
      int ret = 0;
      int i, j, j_beg, j_end, j_ptr, ptr, nnz;
      /* clear rows of the matrix F */
      for (i = 1; i <= n; i++) fr_len[i] = 0;
      /* count non-zeros in rows of the matrix F; count total number of
         non-zeros in this matrix */
      nnz = 0;
      for (j = 1; j <= n; j++)
      {  /* walk through elements of the j-th column and count non-zeros
            in the corresponding rows */
         j_beg = fc_ptr[j];
         j_end = j_beg + fc_len[j] - 1;
         for (j_ptr = j_beg; j_ptr <= j_end; j_ptr++)
            fr_len[sv_ind[j_ptr]]++;
         /* increase total number of non-zeros */
         nnz += fc_len[j];
      }
      /* store total number of non-zeros */
      luf->nnz_f = nnz;
      /* check for free locations */
      if (luf->sv_end - luf->sv_beg < nnz)
      {  /* overflow of the sparse vector area */
         ret = 1;
         goto done;
      }
      /* allocate rows of the matrix F */
      for (i = 1; i <= n; i++)
      {  /* set pointer to the end of the i-th row; later this pointer
            will be set to the beginning of the i-th row */
         fr_ptr[i] = luf->sv_end;
         /* reserve locations for the i-th row */
         luf->sv_end -= fr_len[i];
      }
      /* build the matrix F in row-wise format using this matrix in
         column-wise format */
      for (j = 1; j <= n; j++)
      {  /* walk through elements of the j-th column */
         j_beg = fc_ptr[j];
         j_end = j_beg + fc_len[j] - 1;
         for (j_ptr = j_beg; j_ptr <= j_end; j_ptr++)
         {  /* get row index */
            i = sv_ind[j_ptr];
            /* store element in the i-th row */
            ptr = --fr_ptr[i];
            sv_ind[ptr] = j;
            sv_val[ptr] = sv_val[j_ptr];
         }
      }
done: /* return to the factorizing routine */
      return ret;
}

/***********************************************************************
*  NAME
*
*  luf_factorize - compute LU-factorization
*
*  SYNOPSIS
*
*  #include "glpluf.h"
*  int luf_factorize(LUF *luf, int n, int (*col)(void *info, int j,
*     int ind[], double val[]), void *info);
*
*  DESCRIPTION
*
*  The routine luf_factorize computes LU-factorization of a specified
*  square matrix A.
*
*  The parameter luf specifies LU-factorization program object created
*  by the routine luf_create_it.
*
*  The parameter n specifies the order of A, n > 0.
*
*  The formal routine col specifies the matrix A to be factorized. To
*  obtain j-th column of A the routine luf_factorize calls the routine
*  col with the parameter j (1 <= j <= n). In response the routine col
*  should store row indices and numerical values of non-zero elements
*  of j-th column of A to locations ind[1,...,len] and val[1,...,len],
*  respectively, where len is the number of non-zeros in j-th column
*  returned on exit. Neither zero nor duplicate elements are allowed.
*
*  The parameter info is a transit pointer passed to the routine col.
*
*  RETURNS
*
*  0  LU-factorization has been successfully computed.
*
*  LUF_ESING
*     The specified matrix is singular within the working precision.
*     (On some elimination step the active submatrix is exactly zero,
*     so no pivot can be chosen.)
*
*  LUF_ECOND
*     The specified matrix is ill-conditioned.
*     (On some elimination step too intensive growth of elements of the
*     active submatix has been detected.)
*
*  If matrix A is well scaled, the return code LUF_ECOND may also mean
*  that the threshold pivoting tolerance piv_tol should be increased.
*
*  In case of non-zero return code the factorization becomes invalid.
*  It should not be used in other operations until the cause of failure
*  has been eliminated and the factorization has been recomputed again
*  with the routine luf_factorize.
*
*  REPAIRING SINGULAR MATRIX
*
*  If the routine luf_factorize returns non-zero code, it provides all
*  necessary information that can be used for "repairing" the matrix A,
*  where "repairing" means replacing linearly dependent columns of the
*  matrix A by appropriate columns of the unity matrix. This feature is
*  needed when this routine is used for factorizing the basis matrix
*  within the simplex method procedure.
*
*  On exit linearly dependent columns of the (partially transformed)
*  matrix U have numbers rank+1, rank+2, ..., n, where rank is estimated
*  rank of the matrix A stored by the routine to the member luf->rank.
*  The correspondence between columns of A and U is the same as between
*  columns of V and U. Thus, linearly dependent columns of the matrix A
*  have numbers qq_col[rank+1], qq_col[rank+2], ..., qq_col[n], where
*  qq_col is the column-like representation of the permutation matrix Q.
*  It is understood that each j-th linearly dependent column of the
*  matrix U should be replaced by the unity vector, where all elements
*  are zero except the unity diagonal element u[j,j]. On the other hand
*  j-th row of the matrix U corresponds to the row of the matrix V (and
*  therefore of the matrix A) with the number pp_row[j], where pp_row is
*  the row-like representation of the permutation matrix P. Thus, each
*  j-th linearly dependent column of the matrix U should be replaced by
*  column of the unity matrix with the number pp_row[j].
*
*  The code that repairs the matrix A may look like follows:
*
*     for (j = rank+1; j <= n; j++)
*     {  replace the column qq_col[j] of the matrix A by the column
*        pp_row[j] of the unity matrix;
*     }
*
*  where rank, pp_row, and qq_col are members of the structure LUF. */

int luf_factorize(LUF *luf, int n, int (*col)(void *info, int j,
      int ind[], double val[]), void *info)
{     int *pp_row, *pp_col, *qq_row, *qq_col;
      double max_gro = luf->max_gro;
      int i, j, k, p, q, t, ret;
      if (n < 1)
         xfault("luf_factorize: n = %d; invalid parameter\n", n);
      if (n > N_MAX)
         xfault("luf_factorize: n = %d; matrix too big\n", n);
      /* invalidate the factorization */
      luf->valid = 0;
      /* reallocate arrays, if necessary */
      reallocate(luf, n);
      pp_row = luf->pp_row;
      pp_col = luf->pp_col;
      qq_row = luf->qq_row;
      qq_col = luf->qq_col;
      /* estimate initial size of the SVA, if not specified */
      if (luf->sv_size == 0 && luf->new_sva == 0)
         luf->new_sva = 5 * (n + 10);
more: /* reallocate the sparse vector area, if required */
      if (luf->new_sva > 0)
      {  if (luf->sv_ind != NULL) xfree(luf->sv_ind);
         if (luf->sv_val != NULL) xfree(luf->sv_val);
         luf->sv_size = luf->new_sva;
         luf->sv_ind = xcalloc(1+luf->sv_size, sizeof(int));
         luf->sv_val = xcalloc(1+luf->sv_size, sizeof(double));
         luf->new_sva = 0;
      }
      /* initialize LU-factorization data structures */
      if (initialize(luf, col, info))
      {  /* overflow of the sparse vector area */
         luf->new_sva = luf->sv_size + luf->sv_size;
         xassert(luf->new_sva > luf->sv_size);
         goto more;
      }
      /* main elimination loop */
      for (k = 1; k <= n; k++)
      {  /* choose a pivot element v[p,q] */
         if (find_pivot(luf, &p, &q))
         {  /* no pivot can be chosen, because the active submatrix is
               exactly zero */
            luf->rank = k - 1;
            ret = LUF_ESING;
            goto done;
         }
         /* let v[p,q] correspond to u[i',j']; permute k-th and i'-th
            rows and k-th and j'-th columns of the matrix U = P*V*Q to
            move the element u[i',j'] to the position u[k,k] */
         i = pp_col[p], j = qq_row[q];
         xassert(k <= i && i <= n && k <= j && j <= n);
         /* permute k-th and i-th rows of the matrix U */
         t = pp_row[k];
         pp_row[i] = t, pp_col[t] = i;
         pp_row[k] = p, pp_col[p] = k;
         /* permute k-th and j-th columns of the matrix U */
         t = qq_col[k];
         qq_col[j] = t, qq_row[t] = j;
         qq_col[k] = q, qq_row[q] = k;
         /* eliminate subdiagonal elements of k-th column of the matrix
            U = P*V*Q using the pivot element u[k,k] = v[p,q] */
         if (eliminate(luf, p, q))
         {  /* overflow of the sparse vector area */
            luf->new_sva = luf->sv_size + luf->sv_size;
            xassert(luf->new_sva > luf->sv_size);
            goto more;
         }
         /* check relative growth of elements of the matrix V */
         if (luf->big_v > max_gro * luf->max_a)
         {  /* the growth is too intensive, therefore most probably the
               matrix A is ill-conditioned */
            luf->rank = k - 1;
            ret = LUF_ECOND;
            goto done;
         }
      }
      /* now the matrix U = P*V*Q is upper triangular, the matrix V has
         been built in row-wise format, and the matrix F has been built
         in column-wise format */
      /* defragment the sparse vector area in order to merge all free
         locations in one continuous extent */
      luf_defrag_sva(luf);
      /* build the matrix V in column-wise format */
      if (build_v_cols(luf))
      {  /* overflow of the sparse vector area */
         luf->new_sva = luf->sv_size + luf->sv_size;
         xassert(luf->new_sva > luf->sv_size);
         goto more;
      }
      /* build the matrix F in row-wise format */
      if (build_f_rows(luf))
      {  /* overflow of the sparse vector area */
         luf->new_sva = luf->sv_size + luf->sv_size;
         xassert(luf->new_sva > luf->sv_size);
         goto more;
      }
      /* the LU-factorization has been successfully computed */
      luf->valid = 1;
      luf->rank = n;
      ret = 0;
      /* if there are few free locations in the sparse vector area, try
         increasing its size in the future */
      t = 3 * (n + luf->nnz_v) + 2 * luf->nnz_f;
      if (luf->sv_size < t)
      {  luf->new_sva = luf->sv_size;
         while (luf->new_sva < t)
         {  k = luf->new_sva;
            luf->new_sva = k + k;
            xassert(luf->new_sva > k);
         }
      }
done: /* return to the calling program */
      return ret;
}

/***********************************************************************
*  NAME
*
*  luf_f_solve - solve system F*x = b or F'*x = b
*
*  SYNOPSIS
*
*  #include "glpluf.h"
*  void luf_f_solve(LUF *luf, int tr, double x[]);
*
*  DESCRIPTION
*
*  The routine luf_f_solve solves either the system F*x = b (if the
*  flag tr is zero) or the system F'*x = b (if the flag tr is non-zero),
*  where the matrix F is a component of LU-factorization specified by
*  the parameter luf, F' is a matrix transposed to F.
*
*  On entry the array x should contain elements of the right-hand side
*  vector b in locations x[1], ..., x[n], where n is the order of the
*  matrix F. On exit this array will contain elements of the solution
*  vector x in the same locations. */

void luf_f_solve(LUF *luf, int tr, double x[])
{     int n = luf->n;
      int *fr_ptr = luf->fr_ptr;
      int *fr_len = luf->fr_len;
      int *fc_ptr = luf->fc_ptr;
      int *fc_len = luf->fc_len;
      int *pp_row = luf->pp_row;
      int *sv_ind = luf->sv_ind;
      double *sv_val = luf->sv_val;
      int i, j, k, beg, end, ptr;
      double xk;
      if (!luf->valid)
         xfault("luf_f_solve: LU-factorization is not valid\n");
      if (!tr)
      {  /* solve the system F*x = b */
         for (j = 1; j <= n; j++)
         {  k = pp_row[j];
            xk = x[k];
            if (xk != 0.0)
            {  beg = fc_ptr[k];
               end = beg + fc_len[k] - 1;
               for (ptr = beg; ptr <= end; ptr++)
                  x[sv_ind[ptr]] -= sv_val[ptr] * xk;
            }
         }
      }
      else
      {  /* solve the system F'*x = b */
         for (i = n; i >= 1; i--)
         {  k = pp_row[i];
            xk = x[k];
            if (xk != 0.0)
            {  beg = fr_ptr[k];
               end = beg + fr_len[k] - 1;
               for (ptr = beg; ptr <= end; ptr++)
                  x[sv_ind[ptr]] -= sv_val[ptr] * xk;
            }
         }
      }
      return;
}

/***********************************************************************
*  NAME
*
*  luf_v_solve - solve system V*x = b or V'*x = b
*
*  SYNOPSIS
*
*  #include "glpluf.h"
*  void luf_v_solve(LUF *luf, int tr, double x[]);
*
*  DESCRIPTION
*
*  The routine luf_v_solve solves either the system V*x = b (if the
*  flag tr is zero) or the system V'*x = b (if the flag tr is non-zero),
*  where the matrix V is a component of LU-factorization specified by
*  the parameter luf, V' is a matrix transposed to V.
*
*  On entry the array x should contain elements of the right-hand side
*  vector b in locations x[1], ..., x[n], where n is the order of the
*  matrix V. On exit this array will contain elements of the solution
*  vector x in the same locations. */

void luf_v_solve(LUF *luf, int tr, double x[])
{     int n = luf->n;
      int *vr_ptr = luf->vr_ptr;
      int *vr_len = luf->vr_len;
      double *vr_piv = luf->vr_piv;
      int *vc_ptr = luf->vc_ptr;
      int *vc_len = luf->vc_len;
      int *pp_row = luf->pp_row;
      int *qq_col = luf->qq_col;
      int *sv_ind = luf->sv_ind;
      double *sv_val = luf->sv_val;
      double *b = luf->work;
      int i, j, k, beg, end, ptr;
      double temp;
      if (!luf->valid)
         xfault("luf_v_solve: LU-factorization is not valid\n");
      for (k = 1; k <= n; k++) b[k] = x[k], x[k] = 0.0;
      if (!tr)
      {  /* solve the system V*x = b */
         for (k = n; k >= 1; k--)
         {  i = pp_row[k], j = qq_col[k];
            temp = b[i];
            if (temp != 0.0)
            {  x[j] = (temp /= vr_piv[i]);
               beg = vc_ptr[j];
               end = beg + vc_len[j] - 1;
               for (ptr = beg; ptr <= end; ptr++)
                  b[sv_ind[ptr]] -= sv_val[ptr] * temp;
            }
         }
      }
      else
      {  /* solve the system V'*x = b */
         for (k = 1; k <= n; k++)
         {  i = pp_row[k], j = qq_col[k];
            temp = b[j];
            if (temp != 0.0)
            {  x[i] = (temp /= vr_piv[i]);
               beg = vr_ptr[i];
               end = beg + vr_len[i] - 1;
               for (ptr = beg; ptr <= end; ptr++)
                  b[sv_ind[ptr]] -= sv_val[ptr] * temp;
            }
         }
      }
      return;
}

/***********************************************************************
*  NAME
*
*  luf_a_solve - solve system A*x = b or A'*x = b
*
*  SYNOPSIS
*
*  #include "glpluf.h"
*  void luf_a_solve(LUF *luf, int tr, double x[]);
*
*  DESCRIPTION
*
*  The routine luf_a_solve solves either the system A*x = b (if the
*  flag tr is zero) or the system A'*x = b (if the flag tr is non-zero),
*  where the parameter luf specifies LU-factorization of the matrix A,
*  A' is a matrix transposed to A.
*
*  On entry the array x should contain elements of the right-hand side
*  vector b in locations x[1], ..., x[n], where n is the order of the
*  matrix A. On exit this array will contain elements of the solution
*  vector x in the same locations. */

void luf_a_solve(LUF *luf, int tr, double x[])
{     if (!luf->valid)
         xfault("luf_a_solve: LU-factorization is not valid\n");
      if (!tr)
      {  /* A = F*V, therefore inv(A) = inv(V)*inv(F) */
         luf_f_solve(luf, 0, x);
         luf_v_solve(luf, 0, x);
      }
      else
      {  /* A' = V'*F', therefore inv(A') = inv(F')*inv(V') */
         luf_v_solve(luf, 1, x);
         luf_f_solve(luf, 1, x);
      }
      return;
}

/***********************************************************************
*  NAME
*
*  luf_delete_it - delete LU-factorization
*
*  SYNOPSIS
*
*  #include "glpluf.h"
*  void luf_delete_it(LUF *luf);
*
*  DESCRIPTION
*
*  The routine luf_delete deletes LU-factorization specified by the
*  parameter luf and frees all the memory allocated to this program
*  object. */

void luf_delete_it(LUF *luf)
{     if (luf->fr_ptr != NULL) xfree(luf->fr_ptr);
      if (luf->fr_len != NULL) xfree(luf->fr_len);
      if (luf->fc_ptr != NULL) xfree(luf->fc_ptr);
      if (luf->fc_len != NULL) xfree(luf->fc_len);
      if (luf->vr_ptr != NULL) xfree(luf->vr_ptr);
      if (luf->vr_len != NULL) xfree(luf->vr_len);
      if (luf->vr_cap != NULL) xfree(luf->vr_cap);
      if (luf->vr_piv != NULL) xfree(luf->vr_piv);
      if (luf->vc_ptr != NULL) xfree(luf->vc_ptr);
      if (luf->vc_len != NULL) xfree(luf->vc_len);
      if (luf->vc_cap != NULL) xfree(luf->vc_cap);
      if (luf->pp_row != NULL) xfree(luf->pp_row);
      if (luf->pp_col != NULL) xfree(luf->pp_col);
      if (luf->qq_row != NULL) xfree(luf->qq_row);
      if (luf->qq_col != NULL) xfree(luf->qq_col);
      if (luf->sv_ind != NULL) xfree(luf->sv_ind);
      if (luf->sv_val != NULL) xfree(luf->sv_val);
      if (luf->sv_prev != NULL) xfree(luf->sv_prev);
      if (luf->sv_next != NULL) xfree(luf->sv_next);
      if (luf->vr_max != NULL) xfree(luf->vr_max);
      if (luf->rs_head != NULL) xfree(luf->rs_head);
      if (luf->rs_prev != NULL) xfree(luf->rs_prev);
      if (luf->rs_next != NULL) xfree(luf->rs_next);
      if (luf->cs_head != NULL) xfree(luf->cs_head);
      if (luf->cs_prev != NULL) xfree(luf->cs_prev);
      if (luf->cs_next != NULL) xfree(luf->cs_next);
      if (luf->flag != NULL) xfree(luf->flag);
      if (luf->work != NULL) xfree(luf->work);
      xfree(luf);
      return;
}

/* eof */
