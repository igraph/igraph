/* glplux.h (LU-factorization, bignum arithmetic) */

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

#ifndef GLPLUX_H
#define GLPLUX_H

#include "glpdmp.h"
#include "glpgmp.h"

/*----------------------------------------------------------------------
// The structure LUX defines LU-factorization of a square matrix A,
// which is the following quartet:
//
//    [A] = (F, V, P, Q),                                            (1)
//
// where F and V are such matrices that
//
//    A = F * V,                                                     (2)
//
// and P and Q are such permutation matrices that the matrix
//
//    L = P * F * inv(P)                                             (3)
//
// is lower triangular with unity diagonal, and the matrix
//
//    U = P * V * Q                                                  (4)
//
// is upper triangular. All the matrices have the order n.
//
// The matrices F and V are stored in row/column-wise sparse format as
// row and column linked lists of non-zero elements. Unity elements on
// the main diagonal of the matrix F are not stored. Pivot elements of
// the matrix V (that correspond to diagonal elements of the matrix U)
// are also missing from the row and column lists and stored separately
// in an ordinary array.
//
// The permutation matrices P and Q are stored as ordinary arrays using
// both row- and column-like formats.
//
// The matrices L and U being completely defined by the matrices F, V,
// P, and Q are not stored explicitly.
//
// It is easy to show that the factorization (1)-(3) is some version of
// LU-factorization. Indeed, from (3) and (4) it follows that:
//
//    F = inv(P) * L * P,
//
//    V = inv(P) * U * inv(Q),
//
// and substitution into (2) gives:
//
//    A = F * V = inv(P) * L * U * inv(Q).
//
// For more details see the program documentation. */

typedef struct LUX LUX;
typedef struct LUXELM LUXELM;
typedef struct LUXWKA LUXWKA;

struct LUX
{     /* LU-factorization of a square matrix */
      int n;
      /* the order of matrices A, F, V, P, Q */
      DMP *pool;
      /* memory pool for elements of matrices F and V */
      LUXELM **F_row; /* LUXELM *F_row[1+n]; */
      /* F_row[0] is not used;
         F_row[i], 1 <= i <= n, is a pointer to the list of elements in
         i-th row of matrix F (diagonal elements are not stored) */
      LUXELM **F_col; /* LUXELM *F_col[1+n]; */
      /* F_col[0] is not used;
         F_col[j], 1 <= j <= n, is a pointer to the list of elements in
         j-th column of matrix F (diagonal elements are not stored) */
      mpq_t *V_piv; /* mpq_t V_piv[1+n]; */
      /* V_piv[0] is not used;
         V_piv[p], 1 <= p <= n, is a pivot element v[p,q] corresponding
         to a diagonal element u[k,k] of matrix U = P*V*Q (used on k-th
         elimination step, k = 1, 2, ..., n) */
      LUXELM **V_row; /* LUXELM *V_row[1+n]; */
      /* V_row[0] is not used;
         V_row[i], 1 <= i <= n, is a pointer to the list of elements in
         i-th row of matrix V (except pivot elements) */
      LUXELM **V_col; /* LUXELM *V_col[1+n]; */
      /* V_col[0] is not used;
         V_col[j], 1 <= j <= n, is a pointer to the list of elements in
         j-th column of matrix V (except pivot elements) */
      int *P_row; /* int P_row[1+n]; */
      /* P_row[0] is not used;
         P_row[i] = j means that p[i,j] = 1, where p[i,j] is an element
         of permutation matrix P */
      int *P_col; /* int P_col[1+n]; */
      /* P_col[0] is not used;
         P_col[j] = i means that p[i,j] = 1, where p[i,j] is an element
         of permutation matrix P */
      /* if i-th row or column of matrix F is i'-th row or column of
         matrix L = P*F*inv(P), or if i-th row of matrix V is i'-th row
         of matrix U = P*V*Q, then P_row[i'] = i and P_col[i] = i' */
      int *Q_row; /* int Q_row[1+n]; */
      /* Q_row[0] is not used;
         Q_row[i] = j means that q[i,j] = 1, where q[i,j] is an element
         of permutation matrix Q */
      int *Q_col; /* int Q_col[1+n]; */
      /* Q_col[0] is not used;
         Q_col[j] = i means that q[i,j] = 1, where q[i,j] is an element
         of permutation matrix Q */
      /* if j-th column of matrix V is j'-th column of matrix U = P*V*Q,
         then Q_row[j] = j' and Q_col[j'] = j */
      int rank;
      /* the (exact) rank of matrices A and V */
};

struct LUXELM
{     /* element of matrix F or V */
      int i;
      /* row index, 1 <= i <= m */
      int j;
      /* column index, 1 <= j <= n */
      mpq_t val;
      /* numeric (non-zero) element value */
      LUXELM *r_prev;
      /* pointer to previous element in the same row */
      LUXELM *r_next;
      /* pointer to next element in the same row */
      LUXELM *c_prev;
      /* pointer to previous element in the same column */
      LUXELM *c_next;
      /* pointer to next element in the same column */
};

struct LUXWKA
{     /* working area (used only during factorization) */
      /* in order to efficiently implement Markowitz strategy and Duff
         search technique there are two families {R[0], R[1], ..., R[n]}
         and {C[0], C[1], ..., C[n]}; member R[k] is a set of active
         rows of matrix V having k non-zeros, and member C[k] is a set
         of active columns of matrix V having k non-zeros (in the active
         submatrix); each set R[k] and C[k] is implemented as a separate
         doubly linked list */
      int *R_len; /* int R_len[1+n]; */
      /* R_len[0] is not used;
         R_len[i], 1 <= i <= n, is the number of non-zero elements in
         i-th row of matrix V (that is the length of i-th row) */
      int *R_head; /* int R_head[1+n]; */
      /* R_head[k], 0 <= k <= n, is the number of a first row, which is
         active and whose length is k */
      int *R_prev; /* int R_prev[1+n]; */
      /* R_prev[0] is not used;
         R_prev[i], 1 <= i <= n, is the number of a previous row, which
         is active and has the same length as i-th row */
      int *R_next; /* int R_next[1+n]; */
      /* R_prev[0] is not used;
         R_prev[i], 1 <= i <= n, is the number of a next row, which is
         active and has the same length as i-th row */
      int *C_len; /* int C_len[1+n]; */
      /* C_len[0] is not used;
         C_len[j], 1 <= j <= n, is the number of non-zero elements in
         j-th column of the active submatrix of matrix V (that is the
         length of j-th column in the active submatrix) */
      int *C_head; /* int C_head[1+n]; */
      /* C_head[k], 0 <= k <= n, is the number of a first column, which
         is active and whose length is k */
      int *C_prev; /* int C_prev[1+n]; */
      /* C_prev[0] is not used;
         C_prev[j], 1 <= j <= n, is the number of a previous column,
         which is active and has the same length as j-th column */
      int *C_next; /* int C_next[1+n]; */
      /* C_next[0] is not used;
         C_next[j], 1 <= j <= n, is the number of a next column, which
         is active and has the same length as j-th column */
};

#define lux_create            _glp_lux_create
#define lux_decomp            _glp_lux_decomp
#define lux_f_solve           _glp_lux_f_solve
#define lux_v_solve           _glp_lux_v_solve
#define lux_solve             _glp_lux_solve
#define lux_delete            _glp_lux_delete

LUX *lux_create(int n);
/* create LU-factorization */

int lux_decomp(LUX *lux, int (*col)(void *info, int j, int ind[],
      mpq_t val[]), void *info);
/* compute LU-factorization */

void lux_f_solve(LUX *lux, int tr, mpq_t x[]);
/* solve system F*x = b or F'*x = b */

void lux_v_solve(LUX *lux, int tr, mpq_t x[]);
/* solve system V*x = b or V'*x = b */

void lux_solve(LUX *lux, int tr, mpq_t x[]);
/* solve system A*x = b or A'*x = b */

void lux_delete(LUX *lux);
/* delete LU-factorization */

#endif

/* eof */
