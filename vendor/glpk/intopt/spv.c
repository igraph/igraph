/* spv.c (operations on sparse vectors) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2007-2017 Free Software Foundation, Inc.
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
#include "spv.h"

/***********************************************************************
*  NAME
*
*  spv_create_vec - create sparse vector
*
*  SYNOPSIS
*
*  #include "glpios.h"
*  SPV *spv_create_vec(int n);
*
*  DESCRIPTION
*
*  The routine spv_create_vec creates a sparse vector of dimension n,
*  which initially is a null vector.
*
*  RETURNS
*
*  The routine returns a pointer to the vector created. */

SPV *spv_create_vec(int n)
{     SPV *v;
      xassert(n >= 0);
      v = xmalloc(sizeof(SPV));
      v->n = n;
      v->nnz = 0;
      v->pos = xcalloc(1+n, sizeof(int));
      memset(&v->pos[1], 0, n * sizeof(int));
      v->ind = xcalloc(1+n, sizeof(int));
      v->val = xcalloc(1+n, sizeof(double));
      return v;
}

/***********************************************************************
*  NAME
*
*  spv_check_vec - check that sparse vector has correct representation
*
*  SYNOPSIS
*
*  #include "glpios.h"
*  void spv_check_vec(SPV *v);
*
*  DESCRIPTION
*
*  The routine spv_check_vec checks that a sparse vector specified by
*  the parameter v has correct representation.
*
*  NOTE
*
*  Complexity of this operation is O(n). */

void spv_check_vec(SPV *v)
{     int j, k, nnz;
      xassert(v->n >= 0);
      nnz = 0;
      for (j = v->n; j >= 1; j--)
      {  k = v->pos[j];
         xassert(0 <= k && k <= v->nnz);
         if (k != 0)
         {  xassert(v->ind[k] == j);
            nnz++;
         }
      }
      xassert(v->nnz == nnz);
      return;
}

/***********************************************************************
*  NAME
*
*  spv_get_vj - retrieve component of sparse vector
*
*  SYNOPSIS
*
*  #include "glpios.h"
*  double spv_get_vj(SPV *v, int j);
*
*  RETURNS
*
*  The routine spv_get_vj returns j-th component of a sparse vector
*  specified by the parameter v. */

double spv_get_vj(SPV *v, int j)
{     int k;
      xassert(1 <= j && j <= v->n);
      k = v->pos[j];
      xassert(0 <= k && k <= v->nnz);
      return (k == 0 ? 0.0 : v->val[k]);
}

/***********************************************************************
*  NAME
*
*  spv_set_vj - set/change component of sparse vector
*
*  SYNOPSIS
*
*  #include "glpios.h"
*  void spv_set_vj(SPV *v, int j, double val);
*
*  DESCRIPTION
*
*  The routine spv_set_vj assigns val to j-th component of a sparse
*  vector specified by the parameter v. */

void spv_set_vj(SPV *v, int j, double val)
{     int k;
      xassert(1 <= j && j <= v->n);
      k = v->pos[j];
      if (val == 0.0)
      {  if (k != 0)
         {  /* remove j-th component */
            v->pos[j] = 0;
            if (k < v->nnz)
            {  v->pos[v->ind[v->nnz]] = k;
               v->ind[k] = v->ind[v->nnz];
               v->val[k] = v->val[v->nnz];
            }
            v->nnz--;
         }
      }
      else
      {  if (k == 0)
         {  /* create j-th component */
            k = ++(v->nnz);
            v->pos[j] = k;
            v->ind[k] = j;
         }
         v->val[k] = val;
      }
      return;
}

/***********************************************************************
*  NAME
*
*  spv_clear_vec - set all components of sparse vector to zero
*
*  SYNOPSIS
*
*  #include "glpios.h"
*  void spv_clear_vec(SPV *v);
*
*  DESCRIPTION
*
*  The routine spv_clear_vec sets all components of a sparse vector
*  specified by the parameter v to zero. */

void spv_clear_vec(SPV *v)
{     int k;
      for (k = 1; k <= v->nnz; k++)
         v->pos[v->ind[k]] = 0;
      v->nnz = 0;
      return;
}

/***********************************************************************
*  NAME
*
*  spv_clean_vec - remove zero or small components from sparse vector
*
*  SYNOPSIS
*
*  #include "glpios.h"
*  void spv_clean_vec(SPV *v, double eps);
*
*  DESCRIPTION
*
*  The routine spv_clean_vec removes zero components and components
*  whose magnitude is less than eps from a sparse vector specified by
*  the parameter v. If eps is 0.0, only zero components are removed. */

void spv_clean_vec(SPV *v, double eps)
{     int k, nnz;
      nnz = 0;
      for (k = 1; k <= v->nnz; k++)
      {  if (fabs(v->val[k]) == 0.0 || fabs(v->val[k]) < eps)
         {  /* remove component */
            v->pos[v->ind[k]] = 0;
         }
         else
         {  /* keep component */
            nnz++;
            v->pos[v->ind[k]] = nnz;
            v->ind[nnz] = v->ind[k];
            v->val[nnz] = v->val[k];
         }
      }
      v->nnz = nnz;
      return;
}

/***********************************************************************
*  NAME
*
*  spv_copy_vec - copy sparse vector (x := y)
*
*  SYNOPSIS
*
*  #include "glpios.h"
*  void spv_copy_vec(SPV *x, SPV *y);
*
*  DESCRIPTION
*
*  The routine spv_copy_vec copies a sparse vector specified by the
*  parameter y to a sparse vector specified by the parameter x. */

void spv_copy_vec(SPV *x, SPV *y)
{     int j;
      xassert(x != y);
      xassert(x->n == y->n);
      spv_clear_vec(x);
      x->nnz = y->nnz;
      memcpy(&x->ind[1], &y->ind[1], x->nnz * sizeof(int));
      memcpy(&x->val[1], &y->val[1], x->nnz * sizeof(double));
      for (j = 1; j <= x->nnz; j++)
         x->pos[x->ind[j]] = j;
      return;
}

/***********************************************************************
*  NAME
*
*  spv_linear_comb - compute linear combination (x := x + a * y)
*
*  SYNOPSIS
*
*  #include "glpios.h"
*  void spv_linear_comb(SPV *x, double a, SPV *y);
*
*  DESCRIPTION
*
*  The routine spv_linear_comb computes the linear combination
*
*     x := x + a * y,
*
*  where x and y are sparse vectors, a is a scalar. */

void spv_linear_comb(SPV *x, double a, SPV *y)
{     int j, k;
      double xj, yj;
      xassert(x != y);
      xassert(x->n == y->n);
      for (k = 1; k <= y->nnz; k++)
      {  j = y->ind[k];
         xj = spv_get_vj(x, j);
         yj = y->val[k];
         spv_set_vj(x, j, xj + a * yj);
      }
      return;
}

/***********************************************************************
*  NAME
*
*  spv_delete_vec - delete sparse vector
*
*  SYNOPSIS
*
*  #include "glpios.h"
*  void spv_delete_vec(SPV *v);
*
*  DESCRIPTION
*
*  The routine spv_delete_vec deletes a sparse vector specified by the
*  parameter v freeing all the memory allocated to this object. */

void spv_delete_vec(SPV *v)
{     /* delete sparse vector */
      xfree(v->pos);
      xfree(v->ind);
      xfree(v->val);
      xfree(v);
      return;
}

/* eof */
