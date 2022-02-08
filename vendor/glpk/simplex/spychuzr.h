/* spychuzr.h */

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

#ifndef SPYCHUZR_H
#define SPYCHUZR_H

#include "spxlp.h"

#define spy_chuzr_sel _glp_spy_chuzr_sel
int spy_chuzr_sel(SPXLP *lp, const double beta[/*1+m*/], double tol,
      double tol1, int list[/*1+m*/]);
/* select eligible basic variables */

#define spy_chuzr_std _glp_spy_chuzr_std
int spy_chuzr_std(SPXLP *lp, const double beta[/*1+m*/], int num,
      const int list[]);
/* choose basic variable (dual Dantzig's rule) */

typedef struct SPYSE SPYSE;

struct SPYSE
{     /* dual projected steepest edge and Devex pricing data block */
      int valid;
      /* content validity flag */
      char *refsp; /* char refsp[1+n]; */
      /* refsp[0] is not used;
       * refsp[k], 1 <= k <= n, is the flag meaning that dual variable
       * lambda[k] is in the dual reference space */
      double *gamma; /* double gamma[1+m]; */
      /* gamma[0] is not used;
       * gamma[i], 1 <= i <= m, is the weight for reduced cost r[i]
       * of dual non-basic variable lambdaB[j] in the current basis
       * (r[i] is bound violation for basic variable xB[i]) */
      double *work; /* double work[1+m]; */
      /* working array */
#if 1 /* 30/III-2016 */
      FVS u; /* FVS u[1:m]; */
      /* working vector */
#endif
};

#define spy_alloc_se _glp_spy_alloc_se
void spy_alloc_se(SPXLP *lp, SPYSE *se);
/* allocate dual pricing data block */

#define spy_reset_refsp _glp_spy_reset_refsp
void spy_reset_refsp(SPXLP *lp, SPYSE *se);
/* reset dual reference space */

#define spy_eval_gamma_i _glp_spy_eval_gamma_i
double spy_eval_gamma_i(SPXLP *lp, SPYSE *se, int i);
/* compute dual projected steepest edge weight directly */

#define spy_chuzr_pse _glp_spy_chuzr_pse
int spy_chuzr_pse(SPXLP *lp, SPYSE *se, const double beta[/*1+m*/],
      int num, const int list[]);
/* choose basic variable (dual projected steepest edge) */

#define spy_update_gamma _glp_spy_update_gamma
double spy_update_gamma(SPXLP *lp, SPYSE *se, int p, int q,
      const double trow[/*1+n-m*/], const double tcol[/*1+m*/]);
/* update dual projected steepest edge weights exactly */

#if 1 /* 30/III-2016 */
#define spy_update_gamma_s _glp_spy_update_gamma_s
double spy_update_gamma_s(SPXLP *lp, SPYSE *se, int p, int q,
      const FVS *trow, const FVS *tcol);
/* sparse version of spy_update_gamma */
#endif

#define spy_free_se _glp_spy_free_se
void spy_free_se(SPXLP *lp, SPYSE *se);
/* deallocate dual pricing data block */

#endif

/* eof */
