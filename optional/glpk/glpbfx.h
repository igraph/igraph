/* glpbfx.h (basis factorization interface, bignum arithmetic) */

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

#ifndef GLPBFX_H
#define GLPBFX_H

#include "glpgmp.h"

#ifndef GLPBFX_DEFINED
#define GLPBFX_DEFINED
typedef struct BFX BFX;
#endif

#define bfx_create_binv       _glp_bfx_create_binv
#define bfx_is_valid          _glp_bfx_is_valid
#define bfx_invalidate        _glp_bfx_invalidate
#define bfx_factorize         _glp_bfx_factorize
#define bfx_ftran             _glp_bfx_ftran
#define bfx_btran             _glp_bfx_btran
#define bfx_update            _glp_bfx_update
#define bfx_delete_binv       _glp_bfx_delete_binv

BFX *bfx_create_binv(void);
/* create factorization of the basis matrix */

int bfx_is_valid(BFX *binv);
/* check if factorization is valid */

void bfx_invalidate(BFX *binv);
/* invalidate factorization of the basis matrix */

int bfx_factorize(BFX *binv, int m, int (*col)(void *info, int j,
      int ind[], mpq_t val[]), void *info);
/* compute factorization of the basis matrix */

void bfx_ftran(BFX *binv, mpq_t x[], int save);
/* perform forward transformation (FTRAN) */

void bfx_btran(BFX *binv, mpq_t x[]);
/* perform backward transformation (BTRAN) */

int bfx_update(BFX *binv, int j);
/* update factorization of the basis matrix */

void bfx_delete_binv(BFX *binv);
/* delete factorization of the basis matrix */

#endif

/* eof */
