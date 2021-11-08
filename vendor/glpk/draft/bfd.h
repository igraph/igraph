/* bfd.h (LP basis factorization driver) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2007-2014 Free Software Foundation, Inc.
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

#ifndef BFD_H
#define BFD_H

#if 1 /* 30/III-2016 */
#include "fvs.h"
#endif

typedef struct BFD BFD;

/* return codes: */
#define BFD_ESING    1  /* singular matrix */
#define BFD_ECOND    2  /* ill-conditioned matrix */
#define BFD_ECHECK   3  /* insufficient accuracy */
#define BFD_ELIMIT   4  /* update limit reached */
#if 0 /* 05/III-2014 */
#define BFD_EROOM    5  /* SVA overflow */
#endif

#define bfd_create_it _glp_bfd_create_it
BFD *bfd_create_it(void);
/* create LP basis factorization */

#if 0 /* 08/III-2014 */
#define bfd_set_parm _glp_bfd_set_parm
void bfd_set_parm(BFD *bfd, const void *parm);
/* change LP basis factorization control parameters */
#endif

#define bfd_get_bfcp _glp_bfd_get_bfcp
void bfd_get_bfcp(BFD *bfd, void /* glp_bfcp */ *parm);
/* retrieve LP basis factorization control parameters */

#define bfd_set_bfcp _glp_bfd_set_bfcp
void bfd_set_bfcp(BFD *bfd, const void /* glp_bfcp */ *parm);
/* change LP basis factorization control parameters */

#define bfd_factorize _glp_bfd_factorize
int bfd_factorize(BFD *bfd, int m, /*const int bh[],*/ int (*col)
      (void *info, int j, int ind[], double val[]), void *info);
/* compute LP basis factorization */

#if 1 /* 21/IV-2014 */
#define bfd_condest _glp_bfd_condest
double bfd_condest(BFD *bfd);
/* estimate condition of B */
#endif

#define bfd_ftran _glp_bfd_ftran
void bfd_ftran(BFD *bfd, double x[]);
/* perform forward transformation (solve system B*x = b) */

#if 1 /* 30/III-2016 */
#define bfd_ftran_s _glp_bfd_ftran_s
void bfd_ftran_s(BFD *bfd, FVS *x);
/* sparse version of bfd_ftran */
#endif

#define bfd_btran _glp_bfd_btran
void bfd_btran(BFD *bfd, double x[]);
/* perform backward transformation (solve system B'*x = b) */

#if 1 /* 30/III-2016 */
#define bfd_btran_s _glp_bfd_btran_s
void bfd_btran_s(BFD *bfd, FVS *x);
/* sparse version of bfd_btran */
#endif

#define bfd_update _glp_bfd_update
int bfd_update(BFD *bfd, int j, int len, const int ind[], const double
      val[]);
/* update LP basis factorization */

#define bfd_get_count _glp_bfd_get_count
int bfd_get_count(BFD *bfd);
/* determine factorization update count */

#define bfd_delete_it _glp_bfd_delete_it
void bfd_delete_it(BFD *bfd);
/* delete LP basis factorization */

#endif

/* eof */
