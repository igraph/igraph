/* lufint.h (interface to LU-factorization) */

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

#ifndef LUFINT_H
#define LUFINT_H

#include "sgf.h"

typedef struct LUFINT LUFINT;

struct LUFINT
{     /* interface to LU-factorization */
      int n_max;
      /* maximal value of n (increased automatically) */
      int valid;
      /* factorization is valid only if this flag is set */
      SVA *sva;
      /* sparse vector area (SVA) */
      LUF *luf;
      /* sparse LU-factorization */
      SGF *sgf;
      /* sparse Gaussian factorizer workspace */
      /*--------------------------------------------------------------*/
      /* control parameters */
      int sva_n_max, sva_size;
      /* parameters passed to sva_create_area */
      int delta_n0, delta_n;
      /* if n_max = 0, set n_max = n + delta_n0
       * if n_max < n, set n_max = n + delta_n */
      int sgf_updat;
      double sgf_piv_tol;
      int sgf_piv_lim;
      int sgf_suhl;
      double sgf_eps_tol;
      /* factorizer control parameters */
};

#define lufint_create _glp_lufint_create
LUFINT *lufint_create(void);
/* create interface to LU-factorization */

#define lufint_factorize _glp_lufint_factorize
int lufint_factorize(LUFINT *fi, int n, int (*col)(void *info, int j,
      int ind[], double val[]), void *info);
/* compute LU-factorization of specified matrix A */

#define lufint_delete _glp_lufint_delete
void lufint_delete(LUFINT *fi);
/* delete interface to LU-factorization */

#endif

/* eof */
