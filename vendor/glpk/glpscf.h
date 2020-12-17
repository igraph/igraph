/* glpscf.h (Schur complement factorization) */

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

#ifndef GLPSCF_H
#define GLPSCF_H

/***********************************************************************
*  The structure SCF defines the following factorization of a square
*  nxn matrix C (which is the Schur complement):
*
*     F * C = U * P,
*
*  where F is a square transforming matrix, U is an upper triangular
*  matrix, P is a permutation matrix.
*
*  It is assumed that matrix C is small and dense, so matrices F and U
*  are stored in the dense format by rows as follows:
*
*        1         n       n_max    1         n       n_max
*      1 * * * * * * x x x x      1 * * * * * * x x x x
*        * * * * * * x x x x        . * * * * * x x x x
*        * * * * * * x x x x        . . * * * * x x x x
*        * * * * * * x x x x        . . . * * * x x x x
*        * * * * * * x x x x        . . . . * * x x x x
*      n * * * * * * x x x x      n . . . . . * x x x x
*        x x x x x x x x x x        . . . . . . x x x x
*        x x x x x x x x x x        . . . . . . . x x x
*        x x x x x x x x x x        . . . . . . . . x x
*  n_max x x x x x x x x x x  n_max . . . . . . . . . x
*
*             matrix F                   matrix U
*
*  where '*' are matrix elements, 'x' are reserved locations.
*
*  Permutation matrix P is stored in row-like format.
*
*  Matrix C normally is not stored.
*
*  REFERENCES
*
*  1. M.A.Saunders, "LUSOL: A basis package for constrained optimiza-
*     tion," SCCM, Stanford University, 2006.
*
*  2. M.A.Saunders, "Notes 5: Basis Updates," CME 318, Stanford Univer-
*     sity, Spring 2006.
*
*  3. M.A.Saunders, "Notes 6: LUSOL---a Basis Factorization Package,"
*     ibid. */

typedef struct SCF SCF;

struct SCF
{     /* Schur complement factorization */
      int n_max;
      /* maximal order of matrices C, F, U, P; n_max >= 1 */
      int n;
      /* current order of matrices C, F, U, P; n >= 0 */
      double *f; /* double f[1+n_max*n_max]; */
      /* matrix F stored by rows */
      double *u; /* double u[1+n_max*(n_max+1)/2]; */
      /* upper triangle of matrix U stored by rows */
      int *p; /* int p[1+n_max]; */
      /* matrix P; p[i] = j means that P[i,j] = 1 */
      int t_opt;
      /* type of transformation used to restore triangular structure of
         matrix U: */
#define SCF_TBG      1  /* Bartels-Golub elimination */
#define SCF_TGR      2  /* Givens plane rotation */
      int rank;
      /* estimated rank of matrices C and U */
      double *c; /* double c[1+n_max*n_max]; */
      /* matrix C stored in the same format as matrix F and used only
         for debugging; normally this array is not allocated */
      double *w; /* double w[1+n_max]; */
      /* working array */
};

/* return codes: */
#define SCF_ESING    1  /* singular matrix */
#define SCF_ELIMIT   2  /* update limit reached */

#define scf_create_it _glp_scf_create_it
SCF *scf_create_it(int n_max);
/* create Schur complement factorization */

#define scf_update_exp _glp_scf_update_exp
int scf_update_exp(SCF *scf, const double x[], const double y[],
      double z);
/* update factorization on expanding C */

#define scf_solve_it _glp_scf_solve_it
void scf_solve_it(SCF *scf, int tr, double x[]);
/* solve either system C * x = b or C' * x = b */

#define scf_reset_it _glp_scf_reset_it
void scf_reset_it(SCF *scf);
/* reset factorization for empty matrix C */

#define scf_delete_it _glp_scf_delete_it
void scf_delete_it(SCF *scf);
/* delete Schur complement factorization */

#endif

/* eof */
