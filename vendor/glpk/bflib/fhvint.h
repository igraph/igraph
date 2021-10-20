/* fhvint.h (interface to FHV-factorization) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2012-2014 Free Software Foundation, Inc.
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

#ifndef FHVINT_H
#define FHVINT_H

#include "fhv.h"
#include "lufint.h"

typedef struct FHVINT FHVINT;

struct FHVINT
{     /* interface to FHV-factorization */
      int valid;
      /* factorization is valid only if this flag is set */
      FHV fhv;
      /* FHV-factorization */
      LUFINT *lufi;
      /* interface to underlying LU-factorization */
      /*--------------------------------------------------------------*/
      /* control parameters */
      int nfs_max;
      /* required maximal number of row-like factors */
};

#define fhvint_create _glp_fhvint_create
FHVINT *fhvint_create(void);
/* create interface to FHV-factorization */

#define fhvint_factorize _glp_fhvint_factorize
int fhvint_factorize(FHVINT *fi, int n, int (*col)(void *info, int j,
      int ind[], double val[]), void *info);
/* compute FHV-factorization of specified matrix A */

#define fhvint_update _glp_fhvint_update
int fhvint_update(FHVINT *fi, int j, int len, const int ind[],
      const double val[]);
/* update FHV-factorization after replacing j-th column of A */

#define fhvint_ftran _glp_fhvint_ftran
void fhvint_ftran(FHVINT *fi, double x[]);
/* solve system A * x = b */

#define fhvint_btran _glp_fhvint_btran
void fhvint_btran(FHVINT *fi, double x[]);
/* solve system A'* x = b */

#define fhvint_estimate _glp_fhvint_estimate
double fhvint_estimate(FHVINT *fi);
/* estimate 1-norm of inv(A) */

#define fhvint_delete _glp_fhvint_delete
void fhvint_delete(FHVINT *fi);
/* delete interface to FHV-factorization */

#endif

/* eof */
