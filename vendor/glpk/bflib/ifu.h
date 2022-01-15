/* ifu.h (dense updatable IFU-factorization) */

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

#ifndef IFU_H
#define IFU_H

/***********************************************************************
*  The structure IFU describes dense updatable IFU-factorization.
*
*  The IFU-factorization has the following format:
*
*     A = inv(F) * U,                                                (1)
*
*  where A is a given (unsymmetric) nxn square matrix, F is a square
*  matrix, U is an upper triangular matrix. Obviously, the equality (1)
*  is equivalent to the following equality:
*
*     F * A = U.                                                     (2)
*
*  It is assumed that matrix A is small and dense, so matrices F and U
*  are stored by rows in dense format as follows:
*
*        1         n       n_max      1         n       n_max
*      1 * * * * * * x x x x        1 * * * * * * x x x x
*        * * * * * * x x x x          ? * * * * * x x x x
*        * * * * * * x x x x          ? ? * * * * x x x x
*        * * * * * * x x x x          ? ? ? * * * x x x x
*        * * * * * * x x x x          ? ? ? ? * * x x x x
*      n * * * * * * x x x x        n ? ? ? ? ? * x x x x
*        x x x x x x x x x x          x x x x x x x x x x
*        x x x x x x x x x x          x x x x x x x x x x
*        x x x x x x x x x x          x x x x x x x x x x
*  n_max x x x x x x x x x x    n_max x x x x x x x x x x
*
*             matrix F                     matrix U
*
*  where '*' are matrix elements, '?' are unused locations, 'x' are
*  reserved locations. */

typedef struct IFU IFU;

struct IFU
{     /* IFU-factorization */
      int n_max;
      /* maximal order of matrices A, F, U; n_max >= 1 */
      int n;
      /* current order of matrices A, F, U; 0 <= n <= n_max */
      double *f; /* double f[n_max*n_max]; */
      /* matrix F stored by rows */
      double *u; /* double u[n_max*n_max]; */
      /* matrix U stored by rows */
};

#define ifu_expand _glp_ifu_expand
void ifu_expand(IFU *ifu, double c[/*1+n*/], double r[/*1+n*/],
      double d);
/* expand IFU-factorization */

#define ifu_bg_update _glp_ifu_bg_update
int ifu_bg_update(IFU *ifu, double c[/*1+n*/], double r[/*1+n*/],
      double d);
/* update IFU-factorization (Bartels-Golub) */

#define ifu_gr_update _glp_ifu_gr_update
int ifu_gr_update(IFU *ifu, double c[/*1+n*/], double r[/*1+n*/],
      double d);
/* update IFU-factorization (Givens rotations) */

#define ifu_a_solve _glp_ifu_a_solve
void ifu_a_solve(IFU *ifu, double x[/*1+n*/], double w[/*1+n*/]);
/* solve system A * x = b */

#define ifu_at_solve _glp_ifu_at_solve
void ifu_at_solve(IFU *ifu, double x[/*1+n*/], double w[/*1+n*/]);
/* solve system A'* x = b */

#endif

/* eof */
