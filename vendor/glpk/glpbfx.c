/* glpbfx.c */

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

typedef struct BFX BFX;
#define GLPBFX_DEFINED
#include "glpbfx.h"
#include "glpenv.h"
#include "glplux.h"

struct BFX
{     int valid;
      LUX *lux;
};

BFX *bfx_create_binv(void)
{     /* create factorization of the basis matrix */
      BFX *bfx;
      bfx = xmalloc(sizeof(BFX));
      bfx->valid = 0;
      bfx->lux = NULL;
      return bfx;
}

int bfx_factorize(BFX *binv, int m, int (*col)(void *info, int j,
      int ind[], mpq_t val[]), void *info)
{     /* compute factorization of the basis matrix */
      int ret;
      xassert(m > 0);
      if (binv->lux != NULL && binv->lux->n != m)
      {  lux_delete(binv->lux);
         binv->lux = NULL;
      }
      if (binv->lux == NULL)
         binv->lux = lux_create(m);
      ret = lux_decomp(binv->lux, col, info);
      binv->valid = (ret == 0);
      return ret;
}

void bfx_ftran(BFX *binv, mpq_t x[], int save)
{     /* perform forward transformation (FTRAN) */
      xassert(binv->valid);
      lux_solve(binv->lux, 0, x);
      xassert(save == save);
      return;
}

void bfx_btran(BFX *binv, mpq_t x[])
{     /* perform backward transformation (BTRAN) */
      xassert(binv->valid);
      lux_solve(binv->lux, 1, x);
      return;
}

int bfx_update(BFX *binv, int j)
{     /* update factorization of the basis matrix */
      xassert(binv->valid);
      xassert(1 <= j && j <= binv->lux->n);
      return 1;
}

void bfx_delete_binv(BFX *binv)
{     /* delete factorization of the basis matrix */
      if (binv->lux != NULL)
         lux_delete(binv->lux);
      xfree(binv);
      return;
}

/* eof */
