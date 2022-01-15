/* spxchuzc.h */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2015 Free Software Foundation, Inc.
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

#ifndef SPXCHUZC_H
#define SPXCHUZC_H

#include "spxlp.h"

#define spx_chuzc_sel _glp_spx_chuzc_sel
int spx_chuzc_sel(SPXLP *lp, const double d[/*1+n-m*/], double tol,
      double tol1, int list[/*1+n-m*/]);
/* select eligible non-basic variables */

#define spx_chuzc_std _glp_spx_chuzc_std
int spx_chuzc_std(SPXLP *lp, const double d[/*1+n-m*/], int num,
      const int list[]);
/* choose non-basic variable (Dantzig's rule) */

typedef struct SPXSE SPXSE;

struct SPXSE
{     /* projected steepest edge and Devex pricing data block */
      int valid;
      /* content validity flag */
      char *refsp; /* char refsp[1+n]; */
      /* refsp[0] is not used;
       * refsp[k], 1 <= k <= n, is the flag meaning that variable x[k]
       * is in the reference space */
      double *gamma; /* double gamma[1+n-m]; */
      /* gamma[0] is not used;
       * gamma[j], 1 <= j <= n-m, is the weight for reduced cost d[j]
       * of non-basic variable xN[j] in the current basis */
      double *work; /* double work[1+m]; */
      /* working array */
};

#define spx_alloc_se _glp_spx_alloc_se
void spx_alloc_se(SPXLP *lp, SPXSE *se);
/* allocate pricing data block */

#define spx_reset_refsp _glp_spx_reset_refsp
void spx_reset_refsp(SPXLP *lp, SPXSE *se);
/* reset reference space */

#define spx_eval_gamma_j _glp_spx_eval_gamma_j
double spx_eval_gamma_j(SPXLP *lp, SPXSE *se, int j);
/* compute projeted steepest edge weight directly */

#define spx_chuzc_pse _glp_spx_chuzc_pse
int spx_chuzc_pse(SPXLP *lp, SPXSE *se, const double d[/*1+n-m*/],
      int num, const int list[]);
/* choose non-basic variable (projected steepest edge) */

#define spx_update_gamma _glp_spx_update_gamma
double spx_update_gamma(SPXLP *lp, SPXSE *se, int p, int q,
      const double trow[/*1+n-m*/], const double tcol[/*1+m*/]);
/* update projected steepest edge weights exactly */

#define spx_free_se _glp_spx_free_se
void spx_free_se(SPXLP *lp, SPXSE *se);
/* deallocate pricing data block */

#endif

/* eof */
