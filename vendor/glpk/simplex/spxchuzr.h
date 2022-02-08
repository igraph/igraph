/* spxchuzr.h */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2015-2017 Free Software Foundation, Inc.
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

#ifndef SPXCHUZR_H
#define SPXCHUZR_H

#include "spxlp.h"

#define spx_chuzr_std _glp_spx_chuzr_std
int spx_chuzr_std(SPXLP *lp, int phase, const double beta[/*1+m*/],
      int q, double s, const double tcol[/*1+m*/], int *p_flag,
      double tol_piv, double tol, double tol1);
/* choose basic variable (textbook ratio test) */

#define spx_chuzr_harris _glp_spx_chuzr_harris
int spx_chuzr_harris(SPXLP *lp, int phase, const double beta[/*1+m*/],
      int q, double s, const double tcol[/*1+m*/], int *p_flag,
      double tol_piv, double tol, double tol1);
/* choose basic variable (Harris' ratio test) */

#if 1 /* 22/VI-2017 */
typedef struct SPXBP SPXBP;

struct SPXBP
{     /* penalty function (sum of infeasibilities) break point */
      int i;
      /* basic variable xB[i], 1 <= i <= m, that intersects its bound
       * at this break point
       * i > 0 if xB[i] intersects its lower bound (or fixed value)
       * i < 0 if xB[i] intersects its upper bound
       * i = 0 if xN[q] intersects its opposite bound */
      double teta;
      /* ray parameter value, teta >= 0, at this break point */
      double dc;
      /* increment of the penalty function coefficient cB[i] at this
       * break point */
      double dz;
      /* increment, z[t] - z[0], of the penalty function at this break
       * point */
};

#define spx_ls_eval_bp _glp_spx_ls_eval_bp
int spx_ls_eval_bp(SPXLP *lp, const double beta[/*1+m*/],
      int q, double dq, const double tcol[/*1+m*/], double tol_piv,
      SPXBP bp[/*1+2*m+1*/]);
/* determine penalty function break points */

#define spx_ls_select_bp _glp_spx_ls_select_bp
int spx_ls_select_bp(SPXLP *lp, const double tcol[/*1+m*/],
      int nbp, SPXBP bp[/*1+m+m+1*/], int num, double *slope, double
      teta_lim);
/* select and process penalty function break points */
#endif

#endif

/* eof */
