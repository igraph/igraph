/* spychuzc.h */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2015-2016 Free Software Foundation, Inc.
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

#ifndef SPYCHUZC_H
#define SPYCHUZC_H

#include "spxlp.h"

#define spy_chuzc_std _glp_spy_chuzc_std
int spy_chuzc_std(SPXLP *lp, const double d[/*1+n-m*/],
#if 0 /* 14/III-2016 */
      double s, const double trow[/*1+n-m*/], double tol_piv,
#else
      double r, const double trow[/*1+n-m*/], double tol_piv,
#endif
      double tol, double tol1);
/* choose non-basic variable (dual textbook ratio test) */

#define spy_chuzc_harris _glp_spy_chuzc_harris
int spy_chuzc_harris(SPXLP *lp, const double d[/*1+n-m*/],
#if 0 /* 14/III-2016 */
      double s, const double trow[/*1+n-m*/], double tol_piv,
#else
      double r, const double trow[/*1+n-m*/], double tol_piv,
#endif
      double tol, double tol1);
/* choose non-basic variable (dual Harris' ratio test) */

typedef struct SPYBP SPYBP;

struct SPYBP
{     /* dual objective function break point */
      int j;
      /* dual basic variable lambdaN[j], 1 <= j <= n-m, that intersects
       * zero at this break point */
      double teta;
      /* ray parameter value, teta[j] >= 0, at this break point */
      double dz;
      /* increment, zeta[j] - zeta[0], of the dual objective function
       * at this break point */
};

#if 0 /* 23/III-2016 */
#define spy_eval_bp _glp_spy_eval_bp
int spy_eval_bp(SPXLP *lp, const double d[/*1+n-m*/],
      double r, const double trow[/*1+n-m*/], double tol_piv,
      SPYBP bp[/*1+n-m*/]);
/* determine dual objective function break-points */
#endif

#define spy_ls_eval_bp _glp_spy_ls_eval_bp
int spy_ls_eval_bp(SPXLP *lp, const double d[/*1+n-m*/],
      double r, const double trow[/*1+n-m*/], double tol_piv,
      SPYBP bp[/*1+n-m*/]);
/* determine dual objective function break-points */

#define spy_ls_select_bp _glp_spy_ls_select_bp
int spy_ls_select_bp(SPXLP *lp, const double trow[/*1+n-m*/],
      int nbp, SPYBP bp[/*1+n-m*/], int num, double *slope, double
      teta_lim);
/* select and process dual objective break-points */

#endif

/* eof */
