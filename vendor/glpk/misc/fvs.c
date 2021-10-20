/* fvs.c (sparse vector in FVS format) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2016 Free Software Foundation, Inc.
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

#include "env.h"
#include "fvs.h"

void fvs_alloc_vec(FVS *x, int n)
{     /* allocate sparse vector */
      int j;
      xassert(n >= 0);
      x->n = n;
      x->nnz = 0;
      x->ind = talloc(1+n, int);
      x->vec = talloc(1+n, double);
      for (j = 1; j <= n; j++)
         x->vec[j] = 0.0;
      return;
}

void fvs_check_vec(const FVS *x)
{     /* check sparse vector */
      /* NOTE: for testing/debugging only */
      int n = x->n;
      int nnz = x->nnz;
      int *ind = x->ind;
      double *vec = x->vec;
      char *map;
      int j, k;
      xassert(n >= 0);
      xassert(0 <= nnz && nnz <= n);
      map = talloc(1+n, char);
      for (j = 1; j <= n; j++)
         map[j] = (vec[j] != 0.0);
      for (k = 1; k <= nnz; k++)
      {  j = ind[k];
         xassert(1 <= j && j <= n);
         xassert(map[j]);
         map[j] = 0;
      }
      for (j = 1; j <= n; j++)
         xassert(!map[j]);
      tfree(map);
      return;
}

void fvs_gather_vec(FVS *x, double eps)
{     /* gather sparse vector */
      int n = x->n;
      int *ind = x->ind;
      double *vec = x->vec;
      int j, nnz = 0;
      for (j = n; j >= 1; j--)
      {  if (-eps < vec[j] && vec[j] < +eps)
            vec[j] = 0.0;
         else
            ind[++nnz] = j;
      }
      x->nnz = nnz;
      return;
}

void fvs_clear_vec(FVS *x)
{     /* clear sparse vector */
      int *ind = x->ind;
      double *vec = x->vec;
      int k;
      for (k = x->nnz; k >= 1; k--)
         vec[ind[k]] = 0.0;
      x->nnz = 0;
      return;
}

void fvs_copy_vec(FVS *x, const FVS *y)
{     /* copy sparse vector */
      int *x_ind = x->ind;
      double *x_vec = x->vec;
      int *y_ind = y->ind;
      double *y_vec = y->vec;
      int j, k;
      xassert(x != y);
      xassert(x->n == y->n);
      fvs_clear_vec(x);
      for (k = x->nnz = y->nnz; k >= 1; k--)
      {  j = x_ind[k] = y_ind[k];
         x_vec[j] = y_vec[j];
      }
      return;
}

void fvs_adjust_vec(FVS *x, double eps)
{     /* replace tiny vector elements by exact zeros */
      int nnz = x->nnz;
      int *ind = x->ind;
      double *vec = x->vec;
      int j, k, cnt = 0;
      for (k = 1; k <= nnz; k++)
      {  j = ind[k];
         if (-eps < vec[j] && vec[j] < +eps)
            vec[j] = 0.0;
         else
            ind[++cnt] = j;
      }
      x->nnz = cnt;
      return;
}

void fvs_free_vec(FVS *x)
{     /* deallocate sparse vector */
      tfree(x->ind);
      tfree(x->vec);
      x->n = x->nnz = -1;
      x->ind = NULL;
      x->vec = NULL;
      return;
}

/* eof */
