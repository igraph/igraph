/* glpios04.c (operations on sparse vectors) */

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

#ifdef __clang__
#pragma clang diagnostic ignored "-Wsign-conversion"
#endif

#include "glpios.h"

/***********************************************************************
*  NAME
*
*  ios_create_vec - create sparse vector
*
*  SYNOPSIS
*
*  #include "glpios.h"
*  IOSVEC *ios_create_vec(int n);
*
*  DESCRIPTION
*
*  The routine ios_create_vec creates a sparse vector of dimension n,
*  which initially is a null vector.
*
*  RETURNS
*
*  The routine returns a pointer to the vector created. */

IOSVEC *ios_create_vec(int n)
{     IOSVEC *v;
      xassert(n >= 0);
      v = xmalloc(sizeof(IOSVEC));
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
*  ios_check_vec - check that sparse vector has correct representation
*
*  SYNOPSIS
*
*  #include "glpios.h"
*  void ios_check_vec(IOSVEC *v);
*
*  DESCRIPTION
*
*  The routine ios_check_vec checks that a sparse vector specified by
*  the parameter v has correct representation.
*
*  NOTE
*
*  Complexity of this operation is O(n). */

void ios_check_vec(IOSVEC *v)
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
*  ios_get_vj - retrieve component of sparse vector
*
*  SYNOPSIS
*
*  #include "glpios.h"
*  double ios_get_vj(IOSVEC *v, int j);
*
*  RETURNS
*
*  The routine ios_get_vj returns j-th component of a sparse vector
*  specified by the parameter v. */

double ios_get_vj(IOSVEC *v, int j)
{     int k;
      xassert(1 <= j && j <= v->n);
      k = v->pos[j];
      xassert(0 <= k && k <= v->nnz);
      return (k == 0 ? 0.0 : v->val[k]);
}

/***********************************************************************
*  NAME
*
*  ios_set_vj - set/change component of sparse vector
*
*  SYNOPSIS
*
*  #include "glpios.h"
*  void ios_set_vj(IOSVEC *v, int j, double val);
*
*  DESCRIPTION
*
*  The routine ios_set_vj assigns val to j-th component of a sparse
*  vector specified by the parameter v. */

void ios_set_vj(IOSVEC *v, int j, double val)
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
*  ios_clear_vec - set all components of sparse vector to zero
*
*  SYNOPSIS
*
*  #include "glpios.h"
*  void ios_clear_vec(IOSVEC *v);
*
*  DESCRIPTION
*
*  The routine ios_clear_vec sets all components of a sparse vector
*  specified by the parameter v to zero. */

void ios_clear_vec(IOSVEC *v)
{     int k;
      for (k = 1; k <= v->nnz; k++)
         v->pos[v->ind[k]] = 0;
      v->nnz = 0;
      return;
}

/***********************************************************************
*  NAME
*
*  ios_clean_vec - remove zero or small components from sparse vector
*
*  SYNOPSIS
*
*  #include "glpios.h"
*  void ios_clean_vec(IOSVEC *v, double eps);
*
*  DESCRIPTION
*
*  The routine ios_clean_vec removes zero components and components
*  whose magnitude is less than eps from a sparse vector specified by
*  the parameter v. If eps is 0.0, only zero components are removed. */

void ios_clean_vec(IOSVEC *v, double eps)
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
*  ios_copy_vec - copy sparse vector (x := y)
*
*  SYNOPSIS
*
*  #include "glpios.h"
*  void ios_copy_vec(IOSVEC *x, IOSVEC *y);
*
*  DESCRIPTION
*
*  The routine ios_copy_vec copies a sparse vector specified by the
*  parameter y to a sparse vector specified by the parameter x. */

void ios_copy_vec(IOSVEC *x, IOSVEC *y)
{     int j;
      xassert(x != y);
      xassert(x->n == y->n);
      ios_clear_vec(x);
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
*  ios_linear_comb - compute linear combination (x := x + a * y)
*
*  SYNOPSIS
*
*  #include "glpios.h"
*  void ios_linear_comb(IOSVEC *x, double a, IOSVEC *y);
*
*  DESCRIPTION
*
*  The routine ios_linear_comb computes the linear combination
*
*     x := x + a * y,
*
*  where x and y are sparse vectors, a is a scalar. */

void ios_linear_comb(IOSVEC *x, double a, IOSVEC *y)
{     int j, k;
      double xj, yj;
      xassert(x != y);
      xassert(x->n == y->n);
      for (k = 1; k <= y->nnz; k++)
      {  j = y->ind[k];
         xj = ios_get_vj(x, j);
         yj = y->val[k];
         ios_set_vj(x, j, xj + a * yj);
      }
      return;
}

/***********************************************************************
*  NAME
*
*  ios_delete_vec - delete sparse vector
*
*  SYNOPSIS
*
*  #include "glpios.h"
*  void ios_delete_vec(IOSVEC *v);
*
*  DESCRIPTION
*
*  The routine ios_delete_vec deletes a sparse vector specified by the
*  parameter v freeing all the memory allocated to this object. */

void ios_delete_vec(IOSVEC *v)
{     /* delete sparse vector */
      xfree(v->pos);
      xfree(v->ind);
      xfree(v->val);
      xfree(v);
      return;
}

/* eof */
