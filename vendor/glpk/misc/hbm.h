/* hbm.h (Harwell-Boeing sparse matrix format) */

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

#ifndef HBM_H
#define HBM_H

typedef struct HBM HBM;

struct HBM
{     /* sparse matrix in Harwell-Boeing format; for details see the
         report: I.S.Duff, R.G.Grimes, J.G.Lewis. User's Guide for the
         Harwell-Boeing Sparse Matrix Collection (Release I), 1992 */
      char title[72+1];
      /* matrix title (informative) */
      char key[8+1];
      /* matrix key (informative) */
      char mxtype[3+1];
      /* matrix type:
         R.. real matrix
         C.. complex matrix
         P.. pattern only (no numerical values supplied)
         .S. symmetric (lower triangle + main diagonal)
         .U. unsymmetric
         .H. hermitian (lower triangle + main diagonal)
         .Z. skew symmetric (lower triangle only)
         .R. rectangular
         ..A assembled
         ..E elemental (unassembled) */
      char rhstyp[3+1];
      /* optional types:
         F.. right-hand sides in dense format
         M.. right-hand sides in same format as matrix
         .G. starting vector(s) (guess) is supplied
         ..X exact solution vector(s) is supplied */
      char ptrfmt[16+1];
      /* format for pointers */
      char indfmt[16+1];
      /* format for row (or variable) indices */
      char valfmt[20+1];
      /* format for numerical values of coefficient matrix */
      char rhsfmt[20+1];
      /* format for numerical values of right-hand sides */
      int totcrd;
      /* total number of cards excluding header */
      int ptrcrd;
      /* number of cards for ponters */
      int indcrd;
      /* number of cards for row (or variable) indices */
      int valcrd;
      /* number of cards for numerical values */
      int rhscrd;
      /* number of lines for right-hand sides;
         including starting guesses and solution vectors if present;
         zero indicates no right-hand side data is present */
      int nrow;
      /* number of rows (or variables) */
      int ncol;
      /* number of columns (or elements) */
      int nnzero;
      /* number of row (or variable) indices;
         equal to number of entries for assembled matrix */
      int neltvl;
      /* number of elemental matrix entries;
         zero in case of assembled matrix */
      int nrhs;
      /* number of right-hand sides */
      int nrhsix;
      /* number of row indices;
         ignored in case of unassembled matrix */
      int nrhsvl;
      /* total number of entries in all right-hand sides */
      int nguess;
      /* total number of entries in all starting guesses */
      int nexact;
      /* total number of entries in all solution vectors */
      int *colptr; /* alias: eltptr */
      /* column pointers (in case of assembled matrix);
         elemental matrix pointers (in case of unassembled matrix) */
      int *rowind; /* alias: varind */
      /* row indices (in case of assembled matrix);
         variable indices (in case of unassembled matrix) */
      int *rhsptr;
      /* right-hand side pointers */
      int *rhsind;
      /* right-hand side indices */
      double *values;
      /* matrix values */
      double *rhsval;
      /* right-hand side values */
      double *sguess;
      /* starting guess values */
      double *xexact;
      /* solution vector values */
};

#define hbm_read_mat _glp_hbm_read_mat
HBM *hbm_read_mat(const char *fname);
/* read sparse matrix in Harwell-Boeing format */

#define hbm_free_mat _glp_hbm_free_mat
void hbm_free_mat(HBM *hbm);
/* free sparse matrix in Harwell-Boeing format */

#endif

/* eof */
