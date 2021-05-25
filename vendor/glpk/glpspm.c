/* glpspm.c */

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

#include "glphbm.h"
#include "glprgr.h"
#include "glpspm.h"

/***********************************************************************
*  NAME
*
*  spm_create_mat - create general sparse matrix
*
*  SYNOPSIS
*
*  #include "glpspm.h"
*  SPM *spm_create_mat(int m, int n);
*
*  DESCRIPTION
*
*  The routine spm_create_mat creates a general sparse matrix having
*  m rows and n columns. Being created the matrix is zero (empty), i.e.
*  has no elements.
*
*  RETURNS
*
*  The routine returns a pointer to the matrix created. */

SPM *spm_create_mat(int m, int n)
{     SPM *A;
      xassert(0 <= m && m < INT_MAX);
      xassert(0 <= n && n < INT_MAX);
      A = xmalloc(sizeof(SPM));
      A->m = m;
      A->n = n;
      if (m == 0 || n == 0)
      {  A->pool = NULL;
         A->row = NULL;
         A->col = NULL;
      }
      else
      {  int i, j;
         A->pool = dmp_create_pool();
         A->row = xcalloc(1+m, sizeof(SPME *));
         for (i = 1; i <= m; i++) A->row[i] = NULL;
         A->col = xcalloc(1+n, sizeof(SPME *));
         for (j = 1; j <= n; j++) A->col[j] = NULL;
      }
      return A;
}

/***********************************************************************
*  NAME
*
*  spm_new_elem - add new element to sparse matrix
*
*  SYNOPSIS
*
*  #include "glpspm.h"
*  SPME *spm_new_elem(SPM *A, int i, int j, double val);
*
*  DESCRIPTION
*
*  The routine spm_new_elem adds a new element to the specified sparse
*  matrix. Parameters i, j, and val specify the row number, the column
*  number, and a numerical value of the element, respectively.
*
*  RETURNS
*
*  The routine returns a pointer to the new element added. */

SPME *spm_new_elem(SPM *A, int i, int j, double val)
{     SPME *e;
      xassert(1 <= i && i <= A->m);
      xassert(1 <= j && j <= A->n);
      e = dmp_get_atom(A->pool, sizeof(SPME));
      e->i = i;
      e->j = j;
      e->val = val;
      e->r_prev = NULL;
      e->r_next = A->row[i];
      if (e->r_next != NULL) e->r_next->r_prev = e;
      e->c_prev = NULL;
      e->c_next = A->col[j];
      if (e->c_next != NULL) e->c_next->c_prev = e;
      A->row[i] = A->col[j] = e;
      return e;
}

/***********************************************************************
*  NAME
*
*  spm_delete_mat - delete general sparse matrix
*
*  SYNOPSIS
*
*  #include "glpspm.h"
*  void spm_delete_mat(SPM *A);
*
*  DESCRIPTION
*
*  The routine deletes the specified general sparse matrix freeing all
*  the memory allocated to this object. */

void spm_delete_mat(SPM *A)
{     /* delete sparse matrix */
      if (A->pool != NULL) dmp_delete_pool(A->pool);
      if (A->row != NULL) xfree(A->row);
      if (A->col != NULL) xfree(A->col);
      xfree(A);
      return;
}

/***********************************************************************
*  NAME
*
*  spm_test_mat_e - create test sparse matrix of E(n,c) class
*
*  SYNOPSIS
*
*  #include "glpspm.h"
*  SPM *spm_test_mat_e(int n, int c);
*
*  DESCRIPTION
*
*  The routine spm_test_mat_e creates a test sparse matrix of E(n,c)
*  class as described in the book: Ole 0sterby, Zahari Zlatev. Direct
*  Methods for Sparse Matrices. Springer-Verlag, 1983.
*
*  Matrix of E(n,c) class is a symmetric positive definite matrix of
*  the order n. It has the number 4 on its main diagonal and the number
*  -1 on its four co-diagonals, two of which are neighbour to the main
*  diagonal and two others are shifted from the main diagonal on the
*  distance c.
*
*  It is necessary that n >= 3 and 2 <= c <= n-1.
*
*  RETURNS
*
*  The routine returns a pointer to the matrix created. */

SPM *spm_test_mat_e(int n, int c)
{     SPM *A;
      int i;
      xassert(n >= 3 && 2 <= c && c <= n-1);
      A = spm_create_mat(n, n);
      for (i = 1; i <= n; i++)
         spm_new_elem(A, i, i, 4.0);
      for (i = 1; i <= n-1; i++)
      {  spm_new_elem(A, i, i+1, -1.0);
         spm_new_elem(A, i+1, i, -1.0);
      }
      for (i = 1; i <= n-c; i++)
      {  spm_new_elem(A, i, i+c, -1.0);
         spm_new_elem(A, i+c, i, -1.0);
      }
      return A;
}

/***********************************************************************
*  NAME
*
*  spm_test_mat_d - create test sparse matrix of D(n,c) class
*
*  SYNOPSIS
*
*  #include "glpspm.h"
*  SPM *spm_test_mat_d(int n, int c);
*
*  DESCRIPTION
*
*  The routine spm_test_mat_d creates a test sparse matrix of D(n,c)
*  class as described in the book: Ole 0sterby, Zahari Zlatev. Direct
*  Methods for Sparse Matrices. Springer-Verlag, 1983.
*
*  Matrix of D(n,c) class is a non-singular matrix of the order n. It
*  has unity main diagonal, three co-diagonals above the main diagonal
*  on the distance c, which are cyclically continued below the main
*  diagonal, and a triangle block of the size 10x10 in the upper right
*  corner.
*
*  It is necessary that n >= 14 and 1 <= c <= n-13.
*
*  RETURNS
*
*  The routine returns a pointer to the matrix created. */

SPM *spm_test_mat_d(int n, int c)
{     SPM *A;
      int i, j;
      xassert(n >= 14 && 1 <= c && c <= n-13);
      A = spm_create_mat(n, n);
      for (i = 1; i <= n; i++)
         spm_new_elem(A, i, i, 1.0);
      for (i = 1; i <= n-c; i++)
         spm_new_elem(A, i, i+c, (double)(i+1));
      for (i = n-c+1; i <= n; i++)
         spm_new_elem(A, i, i-n+c, (double)(i+1));
      for (i = 1; i <= n-c-1; i++)
         spm_new_elem(A, i, i+c+1, (double)(-i));
      for (i = n-c; i <= n; i++)
         spm_new_elem(A, i, i-n+c+1, (double)(-i));
      for (i = 1; i <= n-c-2; i++)
         spm_new_elem(A, i, i+c+2, 16.0);
      for (i = n-c-1; i <= n; i++)
         spm_new_elem(A, i, i-n+c+2, 16.0);
      for (j = 1; j <= 10; j++)
         for (i = 1; i <= 11-j; i++)
            spm_new_elem(A, i, n-11+i+j, 100.0 * (double)j);
      return A;
}

/***********************************************************************
*  NAME
*
*  spm_show_mat - write sparse matrix pattern in BMP file format
*
*  SYNOPSIS
*
*  #include "glpspm.h"
*  int spm_show_mat(const SPM *A, const char *fname);
*
*  DESCRIPTION
*
*  The routine spm_show_mat writes pattern of the specified sparse
*  matrix in uncompressed BMP file format (Windows bitmap) to a binary
*  file whose name is specified by the character string fname.
*
*  Each pixel corresponds to one matrix element. The pixel colors have
*  the following meaning:
*
*  Black    structurally zero element
*  White    positive element
*  Cyan     negative element
*  Green    zero element
*  Red      duplicate element
*
*  RETURNS
*
*  If no error occured, the routine returns zero. Otherwise, it prints
*  an appropriate error message and returns non-zero. */

int spm_show_mat(const SPM *A, const char *fname)
{     int m = A->m;
      int n = A->n;
      int i, j, k, ret;
      char *map;
      xprintf("spm_show_mat: writing matrix pattern to `%s'...\n",
         fname);
      xassert(1 <= m && m <= 32767);
      xassert(1 <= n && n <= 32767);
      map = xmalloc(m * n);
      memset(map, 0x08, m * n);
      for (i = 1; i <= m; i++)
      {  SPME *e;
         for (e = A->row[i]; e != NULL; e = e->r_next)
         {  j = e->j;
            xassert(1 <= j && j <= n);
            k = n * (i - 1) + (j - 1);
            if (map[k] != 0x08)
               map[k] = 0x0C;
            else if (e->val > 0.0)
               map[k] = 0x0F;
            else if (e->val < 0.0)
               map[k] = 0x0B;
            else
               map[k] = 0x0A;
         }
      }
      ret = rgr_write_bmp16(fname, m, n, map);
      xfree(map);
      return ret;
}

/***********************************************************************
*  NAME
*
*  spm_read_hbm - read sparse matrix in Harwell-Boeing format
*
*  SYNOPSIS
*
*  #include "glpspm.h"
*  SPM *spm_read_hbm(const char *fname);
*
*  DESCRIPTION
*
*  The routine spm_read_hbm reads a sparse matrix in the Harwell-Boeing
*  format from a text file whose name is the character string fname.
*
*  Detailed description of the Harwell-Boeing format recognised by this
*  routine can be found in the following report:
*
*  I.S.Duff, R.G.Grimes, J.G.Lewis. User's Guide for the Harwell-Boeing
*  Sparse Matrix Collection (Release I), TR/PA/92/86, October 1992.
*
*  NOTE
*
*  The routine spm_read_hbm reads the matrix "as is", due to which zero
*  and/or duplicate elements can appear in the matrix.
*
*  RETURNS
*
*  If no error occured, the routine returns a pointer to the matrix
*  created. Otherwise, the routine prints an appropriate error message
*  and returns NULL. */

SPM *spm_read_hbm(const char *fname)
{     SPM *A = NULL;
      HBM *hbm;
      int nrow, ncol, nnzero, i, j, beg, end, ptr, *colptr, *rowind;
      double val, *values;
      char *mxtype;
      hbm = hbm_read_mat(fname);
      if (hbm == NULL)
      {  xprintf("spm_read_hbm: unable to read matrix\n");
         goto fini;
      }
      mxtype = hbm->mxtype;
      nrow = hbm->nrow;
      ncol = hbm->ncol;
      nnzero = hbm->nnzero;
      colptr = hbm->colptr;
      rowind = hbm->rowind;
      values = hbm->values;
      if (!(strcmp(mxtype, "RSA") == 0 || strcmp(mxtype, "PSA") == 0 ||
            strcmp(mxtype, "RUA") == 0 || strcmp(mxtype, "PUA") == 0 ||
            strcmp(mxtype, "RRA") == 0 || strcmp(mxtype, "PRA") == 0))
      {  xprintf("spm_read_hbm: matrix type `%s' not supported\n",
            mxtype);
         goto fini;
      }
      A = spm_create_mat(nrow, ncol);
      if (mxtype[1] == 'S' || mxtype[1] == 'U')
         xassert(nrow == ncol);
      for (j = 1; j <= ncol; j++)
      {  beg = colptr[j];
         end = colptr[j+1];
         xassert(1 <= beg && beg <= end && end <= nnzero + 1);
         for (ptr = beg; ptr < end; ptr++)
         {  i = rowind[ptr];
            xassert(1 <= i && i <= nrow);
            if (mxtype[0] == 'R')
               val = values[ptr];
            else
               val = 1.0;
            spm_new_elem(A, i, j, val);
            if (mxtype[1] == 'S' && i != j)
               spm_new_elem(A, j, i, val);
         }
      }
fini: if (hbm != NULL) hbm_free_mat(hbm);
      return A;
}

/***********************************************************************
*  NAME
*
*  spm_count_nnz - determine number of non-zeros in sparse matrix
*
*  SYNOPSIS
*
*  #include "glpspm.h"
*  int spm_count_nnz(const SPM *A);
*
*  RETURNS
*
*  The routine spm_count_nnz returns the number of structural non-zero
*  elements in the specified sparse matrix. */

int spm_count_nnz(const SPM *A)
{     SPME *e;
      int i, nnz = 0;
      for (i = 1; i <= A->m; i++)
         for (e = A->row[i]; e != NULL; e = e->r_next) nnz++;
      return nnz;
}

/***********************************************************************
*  NAME
*
*  spm_drop_zeros - remove zero elements from sparse matrix
*
*  SYNOPSIS
*
*  #include "glpspm.h"
*  int spm_drop_zeros(SPM *A, double eps);
*
*  DESCRIPTION
*
*  The routine spm_drop_zeros removes all elements from the specified
*  sparse matrix, whose absolute value is less than eps.
*
*  If the parameter eps is 0, only zero elements are removed from the
*  matrix.
*
*  RETURNS
*
*  The routine returns the number of elements removed. */

int spm_drop_zeros(SPM *A, double eps)
{     SPME *e, *next;
      int i, count = 0;
      for (i = 1; i <= A->m; i++)
      {  for (e = A->row[i]; e != NULL; e = next)
         {  next = e->r_next;
            if (e->val == 0.0 || fabs(e->val) < eps)
            {  /* remove element from the row list */
               if (e->r_prev == NULL)
                  A->row[e->i] = e->r_next;
               else
                  e->r_prev->r_next = e->r_next;
               if (e->r_next == NULL)
                  ;
               else
                  e->r_next->r_prev = e->r_prev;
               /* remove element from the column list */
               if (e->c_prev == NULL)
                  A->col[e->j] = e->c_next;
               else
                  e->c_prev->c_next = e->c_next;
               if (e->c_next == NULL)
                  ;
               else
                  e->c_next->c_prev = e->c_prev;
               /* return element to the memory pool */
               dmp_free_atom(A->pool, e, sizeof(SPME));
               count++;
            }
         }
      }
      return count;
}

/***********************************************************************
*  NAME
*
*  spm_read_mat - read sparse matrix from text file
*
*  SYNOPSIS
*
*  #include "glpspm.h"
*  SPM *spm_read_mat(const char *fname);
*
*  DESCRIPTION
*
*  The routine reads a sparse matrix from a text file whose name is
*  specified by the parameter fname.
*
*  For the file format see description of the routine spm_write_mat.
*
*  RETURNS
*
*  On success the routine returns a pointer to the matrix created,
*  otherwise NULL. */

#if 1
SPM *spm_read_mat(const char *fname)
{     xassert(fname != fname);
      return NULL;
}
#else
SPM *spm_read_mat(const char *fname)
{     SPM *A = NULL;
      PDS *pds;
      jmp_buf jump;
      int i, j, k, m, n, nnz, fail = 0;
      double val;
      xprintf("spm_read_mat: reading matrix from `%s'...\n", fname);
      pds = pds_open_file(fname);
      if (pds == NULL)
      {  xprintf("spm_read_mat: unable to open `%s' - %s\n", fname,
            strerror(errno));
         fail = 1;
         goto done;
      }
      if (setjmp(jump))
      {  fail = 1;
         goto done;
      }
      pds_set_jump(pds, jump);
      /* number of rows, number of columns, number of non-zeros */
      m = pds_scan_int(pds);
      if (m < 0)
         pds_error(pds, "invalid number of rows\n");
      n = pds_scan_int(pds);
      if (n < 0)
         pds_error(pds, "invalid number of columns\n");
      nnz = pds_scan_int(pds);
      if (nnz < 0)
         pds_error(pds, "invalid number of non-zeros\n");
      /* create matrix */
      xprintf("spm_read_mat: %d rows, %d columns, %d non-zeros\n",
         m, n, nnz);
      A = spm_create_mat(m, n);
      /* read matrix elements */
      for (k = 1; k <= nnz; k++)
      {  /* row index, column index, element value */
         i = pds_scan_int(pds);
         if (!(1 <= i && i <= m))
            pds_error(pds, "row index out of range\n");
         j = pds_scan_int(pds);
         if (!(1 <= j && j <= n))
            pds_error(pds, "column index out of range\n");
         val = pds_scan_num(pds);
         /* add new element to the matrix */
         spm_new_elem(A, i, j, val);
      }
      xprintf("spm_read_mat: %d lines were read\n", pds->count);
done: if (pds != NULL) pds_close_file(pds);
      if (fail && A != NULL) spm_delete_mat(A), A = NULL;
      return A;
}
#endif

/***********************************************************************
*  NAME
*
*  spm_write_mat - write sparse matrix to text file
*
*  SYNOPSIS
*
*  #include "glpspm.h"
*  int spm_write_mat(const SPM *A, const char *fname);
*
*  DESCRIPTION
*
*  The routine spm_write_mat writes the specified sparse matrix to a
*  text file whose name is specified by the parameter fname. This file
*  can be read back with the routine spm_read_mat.
*
*  RETURNS
*
*  On success the routine returns zero, otherwise non-zero.
*
*  FILE FORMAT
*
*  The file created by the routine spm_write_mat is a plain text file,
*  which contains the following information:
*
*     m n nnz
*     row[1] col[1] val[1]
*     row[2] col[2] val[2]
*     . . .
*     row[nnz] col[nnz] val[nnz]
*
*  where:
*  m is the number of rows;
*  n is the number of columns;
*  nnz is the number of non-zeros;
*  row[k], k = 1,...,nnz, are row indices;
*  col[k], k = 1,...,nnz, are column indices;
*  val[k], k = 1,...,nnz, are element values. */

#if 1
int spm_write_mat(const SPM *A, const char *fname)
{     xassert(A != A);
      xassert(fname != fname);
      return 0;
}
#else
int spm_write_mat(const SPM *A, const char *fname)
{     FILE *fp;
      int i, nnz, ret = 0;
      xprintf("spm_write_mat: writing matrix to `%s'...\n", fname);
      fp = fopen(fname, "w");
      if (fp == NULL)
      {  xprintf("spm_write_mat: unable to create `%s' - %s\n", fname,
            strerror(errno));
         ret = 1;
         goto done;
      }
      /* number of rows, number of columns, number of non-zeros */
      nnz = spm_count_nnz(A);
      fprintf(fp, "%d %d %d\n", A->m, A->n, nnz);
      /* walk through rows of the matrix */
      for (i = 1; i <= A->m; i++)
      {  SPME *e;
         /* walk through elements of i-th row */
         for (e = A->row[i]; e != NULL; e = e->r_next)
         {  /* row index, column index, element value */
            fprintf(fp, "%d %d %.*g\n", e->i, e->j, DBL_DIG, e->val);
         }
      }
      fflush(fp);
      if (ferror(fp))
      {  xprintf("spm_write_mat: writing error on `%s' - %s\n", fname,
            strerror(errno));
         ret = 1;
         goto done;
      }
      xprintf("spm_write_mat: %d lines were written\n", 1 + nnz);
done: if (fp != NULL) fclose(fp);
      return ret;
}
#endif

/***********************************************************************
*  NAME
*
*  spm_transpose - transpose sparse matrix
*
*  SYNOPSIS
*
*  #include "glpspm.h"
*  SPM *spm_transpose(const SPM *A);
*
*  RETURNS
*
*  The routine computes and returns sparse matrix B, which is a matrix
*  transposed to sparse matrix A. */

SPM *spm_transpose(const SPM *A)
{     SPM *B;
      int i;
      B = spm_create_mat(A->n, A->m);
      for (i = 1; i <= A->m; i++)
      {  SPME *e;
         for (e = A->row[i]; e != NULL; e = e->r_next)
            spm_new_elem(B, e->j, i, e->val);
      }
      return B;
}

SPM *spm_add_sym(const SPM *A, const SPM *B)
{     /* add two sparse matrices (symbolic phase) */
      SPM *C;
      int i, j, *flag;
      xassert(A->m == B->m);
      xassert(A->n == B->n);
      /* create resultant matrix */
      C = spm_create_mat(A->m, A->n);
      /* allocate and clear the flag array */
      flag = xcalloc(1+C->n, sizeof(int));
      for (j = 1; j <= C->n; j++)
         flag[j] = 0;
      /* compute pattern of C = A + B */
      for (i = 1; i <= C->m; i++)
      {  SPME *e;
         /* at the beginning i-th row of C is empty */
         /* (i-th row of C) := (i-th row of C) union (i-th row of A) */
         for (e = A->row[i]; e != NULL; e = e->r_next)
         {  /* (note that i-th row of A may have duplicate elements) */
            j = e->j;
            if (!flag[j])
            {  spm_new_elem(C, i, j, 0.0);
               flag[j] = 1;
            }
         }
         /* (i-th row of C) := (i-th row of C) union (i-th row of B) */
         for (e = B->row[i]; e != NULL; e = e->r_next)
         {  /* (note that i-th row of B may have duplicate elements) */
            j = e->j;
            if (!flag[j])
            {  spm_new_elem(C, i, j, 0.0);
               flag[j] = 1;
            }
         }
         /* reset the flag array */
         for (e = C->row[i]; e != NULL; e = e->r_next)
            flag[e->j] = 0;
      }
      /* check and deallocate the flag array */
      for (j = 1; j <= C->n; j++)
         xassert(!flag[j]);
      xfree(flag);
      return C;
}

void spm_add_num(SPM *C, double alfa, const SPM *A, double beta,
      const SPM *B)
{     /* add two sparse matrices (numeric phase) */
      int i, j;
      double *work;
      /* allocate and clear the working array */
      work = xcalloc(1+C->n, sizeof(double));
      for (j = 1; j <= C->n; j++)
         work[j] = 0.0;
      /* compute matrix C = alfa * A + beta * B */
      for (i = 1; i <= C->n; i++)
      {  SPME *e;
         /* work := alfa * (i-th row of A) + beta * (i-th row of B) */
         /* (note that A and/or B may have duplicate elements) */
         for (e = A->row[i]; e != NULL; e = e->r_next)
            work[e->j] += alfa * e->val;
         for (e = B->row[i]; e != NULL; e = e->r_next)
            work[e->j] += beta * e->val;
         /* (i-th row of C) := work, work := 0 */
         for (e = C->row[i]; e != NULL; e = e->r_next)
         {  j = e->j;
            e->val = work[j];
            work[j] = 0.0;
         }
      }
      /* check and deallocate the working array */
      for (j = 1; j <= C->n; j++)
         xassert(work[j] == 0.0);
      xfree(work);
      return;
}

SPM *spm_add_mat(double alfa, const SPM *A, double beta, const SPM *B)
{     /* add two sparse matrices (driver routine) */
      SPM *C;
      C = spm_add_sym(A, B);
      spm_add_num(C, alfa, A, beta, B);
      return C;
}

SPM *spm_mul_sym(const SPM *A, const SPM *B)
{     /* multiply two sparse matrices (symbolic phase) */
      int i, j, k, *flag;
      SPM *C;
      xassert(A->n == B->m);
      /* create resultant matrix */
      C = spm_create_mat(A->m, B->n);
      /* allocate and clear the flag array */
      flag = xcalloc(1+C->n, sizeof(int));
      for (j = 1; j <= C->n; j++)
         flag[j] = 0;
      /* compute pattern of C = A * B */
      for (i = 1; i <= C->m; i++)
      {  SPME *e, *ee;
         /* compute pattern of i-th row of C */
         for (e = A->row[i]; e != NULL; e = e->r_next)
         {  k = e->j;
            for (ee = B->row[k]; ee != NULL; ee = ee->r_next)
            {  j = ee->j;
               /* if a[i,k] != 0 and b[k,j] != 0 then c[i,j] != 0 */
               if (!flag[j])
               {  /* c[i,j] does not exist, so create it */
                  spm_new_elem(C, i, j, 0.0);
                  flag[j] = 1;
               }
            }
         }
         /* reset the flag array */
         for (e = C->row[i]; e != NULL; e = e->r_next)
            flag[e->j] = 0;
      }
      /* check and deallocate the flag array */
      for (j = 1; j <= C->n; j++)
         xassert(!flag[j]);
      xfree(flag);
      return C;
}

void spm_mul_num(SPM *C, const SPM *A, const SPM *B)
{     /* multiply two sparse matrices (numeric phase) */
      int i, j;
      double *work;
      /* allocate and clear the working array */
      work = xcalloc(1+A->n, sizeof(double));
      for (j = 1; j <= A->n; j++)
         work[j] = 0.0;
      /* compute matrix C = A * B */
      for (i = 1; i <= C->m; i++)
      {  SPME *e, *ee;
         double temp;
         /* work := (i-th row of A) */
         /* (note that A may have duplicate elements) */
         for (e = A->row[i]; e != NULL; e = e->r_next)
            work[e->j] += e->val;
         /* compute i-th row of C */
         for (e = C->row[i]; e != NULL; e = e->r_next)
         {  j = e->j;
            /* c[i,j] := work * (j-th column of B) */
            temp = 0.0;
            for (ee = B->col[j]; ee != NULL; ee = ee->c_next)
               temp += work[ee->i] * ee->val;
            e->val = temp;
         }
         /* reset the working array */
         for (e = A->row[i]; e != NULL; e = e->r_next)
            work[e->j] = 0.0;
      }
      /* check and deallocate the working array */
      for (j = 1; j <= A->n; j++)
         xassert(work[j] == 0.0);
      xfree(work);
      return;
}

SPM *spm_mul_mat(const SPM *A, const SPM *B)
{     /* multiply two sparse matrices (driver routine) */
      SPM *C;
      C = spm_mul_sym(A, B);
      spm_mul_num(C, A, B);
      return C;
}

PER *spm_create_per(int n)
{     /* create permutation matrix */
      PER *P;
      int k;
      xassert(n >= 0);
      P = xmalloc(sizeof(PER));
      P->n = n;
      P->row = xcalloc(1+n, sizeof(int));
      P->col = xcalloc(1+n, sizeof(int));
      /* initially it is identity matrix */
      for (k = 1; k <= n; k++)
         P->row[k] = P->col[k] = k;
      return P;
}

void spm_check_per(PER *P)
{     /* check permutation matrix for correctness */
      int i, j;
      xassert(P->n >= 0);
      for (i = 1; i <= P->n; i++)
      {  j = P->row[i];
         xassert(1 <= j && j <= P->n);
         xassert(P->col[j] == i);
      }
      return;
}

void spm_delete_per(PER *P)
{     /* delete permutation matrix */
      xfree(P->row);
      xfree(P->col);
      xfree(P);
      return;
}

/* eof */
