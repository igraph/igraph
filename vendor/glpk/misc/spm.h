/* spm.h (general sparse matrices) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2004-2018 Free Software Foundation, Inc.
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

#ifndef SPM_H
#define SPM_H

#include "dmp.h"

typedef struct SPM SPM;
typedef struct SPME SPME;

struct SPM
{     /* general sparse matrix */
      int m;
      /* number of rows, m >= 0 */
      int n;
      /* number of columns, n >= 0 */
      DMP *pool;
      /* memory pool to store matrix elements */
      SPME **row; /* SPME *row[1+m]; */
      /* row[i], 1 <= i <= m, is a pointer to i-th row list */
      SPME **col; /* SPME *col[1+n]; */
      /* col[j], 1 <= j <= n, is a pointer to j-th column list */
};

struct SPME
{     /* sparse matrix element */
      int i;
      /* row number */
      int j;
      /* column number */
      double val;
      /* element value */
      SPME *r_prev;
      /* pointer to previous element in the same row */
      SPME *r_next;
      /* pointer to next element in the same row */
      SPME *c_prev;
      /* pointer to previous element in the same column */
      SPME *c_next;
      /* pointer to next element in the same column */
};

typedef struct PER PER;

struct PER
{     /* permutation matrix */
      int n;
      /* matrix order, n >= 0 */
      int *row; /* int row[1+n]; */
      /* row[i] = j means p[i,j] = 1 */
      int *col; /* int col[1+n]; */
      /* col[j] = i means p[i,j] = 1 */
};

#define spm_create_mat _glp_spm_create_mat
SPM *spm_create_mat(int m, int n);
/* create general sparse matrix */

#define spm_new_elem _glp_spm_new_elem
SPME *spm_new_elem(SPM *A, int i, int j, double val);
/* add new element to sparse matrix */

#define spm_delete_mat _glp_spm_delete_mat
void spm_delete_mat(SPM *A);
/* delete general sparse matrix */

#define spm_test_mat_e _glp_spm_test_mat_e
SPM *spm_test_mat_e(int n, int c);
/* create test sparse matrix of E(n,c) class */

#define spm_test_mat_d _glp_spm_test_mat_d
SPM *spm_test_mat_d(int n, int c);
/* create test sparse matrix of D(n,c) class */

#define spm_show_mat _glp_spm_show_mat
int spm_show_mat(const SPM *A, const char *fname);
/* write sparse matrix pattern in BMP file format */

#define spm_read_hbm _glp_spm_read_hbm
SPM *spm_read_hbm(const char *fname);
/* read sparse matrix in Harwell-Boeing format */

#define spm_count_nnz _glp_spm_count_nnz
int spm_count_nnz(const SPM *A);
/* determine number of non-zeros in sparse matrix */

#define spm_drop_zeros _glp_spm_drop_zeros
int spm_drop_zeros(SPM *A, double eps);
/* remove zero elements from sparse matrix */

#define spm_read_mat _glp_spm_read_mat
SPM *spm_read_mat(const char *fname);
/* read sparse matrix from text file */

#define spm_write_mat _glp_spm_write_mat
int spm_write_mat(const SPM *A, const char *fname);
/* write sparse matrix to text file */

#define spm_transpose _glp_spm_transpose
SPM *spm_transpose(const SPM *A);
/* transpose sparse matrix */

#define spm_add_sym _glp_spm_add_sym
SPM *spm_add_sym(const SPM *A, const SPM *B);
/* add two sparse matrices (symbolic phase) */

#define spm_add_num _glp_spm_add_num
void spm_add_num(SPM *C, double alfa, const SPM *A, double beta,
      const SPM *B);
/* add two sparse matrices (numeric phase) */

#define spm_add_mat _glp_spm_add_mat
SPM *spm_add_mat(double alfa, const SPM *A, double beta,
      const SPM *B);
/* add two sparse matrices (driver routine) */

#define spm_mul_sym _glp_spm_mul_sym
SPM *spm_mul_sym(const SPM *A, const SPM *B);
/* multiply two sparse matrices (symbolic phase) */

#define spm_mul_num _glp_spm_mul_num
void spm_mul_num(SPM *C, const SPM *A, const SPM *B);
/* multiply two sparse matrices (numeric phase) */

#define spm_mul_mat _glp_spm_mul_mat
SPM *spm_mul_mat(const SPM *A, const SPM *B);
/* multiply two sparse matrices (driver routine) */

#define spm_create_per _glp_spm_create_per
PER *spm_create_per(int n);
/* create permutation matrix */

#define spm_check_per _glp_spm_check_per
void spm_check_per(PER *P);
/* check permutation matrix for correctness */

#define spm_delete_per _glp_spm_delete_per
void spm_delete_per(PER *P);
/* delete permutation matrix */

#endif

/* eof */
