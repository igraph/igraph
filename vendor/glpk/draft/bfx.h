/* bfx.h (LP basis factorization driver, rational arithmetic) */

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

#ifndef BFX_H
#define BFX_H

#include "mygmp.h"

typedef struct BFX BFX;

#define bfx_create_binv _glp_bfx_create_binv
BFX *bfx_create_binv(void);
/* create factorization of the basis matrix */

#define bfx_is_valid _glp_bfx_is_valid
int bfx_is_valid(BFX *binv);
/* check if factorization is valid */

#define bfx_invalidate _glp_bfx_invalidate
void bfx_invalidate(BFX *binv);
/* invalidate factorization of the basis matrix */

#define bfx_factorize _glp_bfx_factorize
int bfx_factorize(BFX *binv, int m, int (*col)(void *info, int j,
      int ind[], mpq_t val[]), void *info);
/* compute factorization of the basis matrix */

#define bfx_ftran _glp_bfx_ftran
void bfx_ftran(BFX *binv, mpq_t x[], int save);
/* perform forward transformation (FTRAN) */

#define bfx_btran _glp_bfx_btran
void bfx_btran(BFX *binv, mpq_t x[]);
/* perform backward transformation (BTRAN) */

#define bfx_update _glp_bfx_update
int bfx_update(BFX *binv, int j);
/* update factorization of the basis matrix */

#define bfx_delete_binv _glp_bfx_delete_binv
void bfx_delete_binv(BFX *binv);
/* delete factorization of the basis matrix */

#endif

/* eof */
