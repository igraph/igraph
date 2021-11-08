/* ks.c (0-1 knapsack problem) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2017-2018 Free Software Foundation, Inc.
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
#include "ks.h"
#include "mt1.h"

/***********************************************************************
*  0-1 knapsack problem has the following formulation:
*
*     maximize z = sum{j in 1..n} c[j]x[j]                           (1)
*
*         s.t. sum{j in 1..n} a[j]x[j] <= b                          (2)
*
*              x[j] in {0, 1} for all j in 1..n                      (3)
*
*  In general case it is assumed that the instance is non-normalized,
*  i.e. parameters a, b, and c may have any sign.
***********************************************************************/

/***********************************************************************
*  ks_enum - solve 0-1 knapsack problem by complete enumeration
*
*  This routine finds optimal solution to 0-1 knapsack problem (1)-(3)
*  by complete enumeration. It is intended mainly for testing purposes.
*
*  The instance to be solved is specified by parameters n, a, b, and c.
*  Note that these parameters can have any sign, i.e. normalization is
*  not needed.
*
*  On exit the routine stores the optimal point found in locations
*  x[1], ..., x[n] and returns the optimal objective value. However, if
*  the instance is infeasible, the routine returns INT_MIN.
*
*  Since the complete enumeration is inefficient, this routine can be
*  used only for small instances (n <= 20-30). */

#define N_MAX 40

int ks_enum(int n, const int a[/*1+n*/], int b, const int c[/*1+n*/],
      char x[/*1+n*/])
{     int j, s, z, z_best;
      char x_best[1+N_MAX];
      xassert(0 <= n && n <= N_MAX);
      /* initialization */
      memset(&x[1], 0, n * sizeof(char));
      z_best = INT_MIN;
loop: /* compute constraint and objective at current x */
      s = z = 0;
      for (j = 1; j <= n; j++)
      {  if (x[j])
            s += a[j], z += c[j];
      }
      /* check constraint violation */
      if (s > b)
         goto next;
      /* check objective function */
      if (z_best < z)
      {  /* better solution has been found */
         memcpy(&x_best[1], &x[1], n * sizeof(char));
         z_best = z;
      }
next: /* generate next x */
      for (j = 1; j <= n; j++)
      {  if (!x[j])
         {  x[j] = 1;
            goto loop;
         }
         x[j] = 0;
      }
      /* report best (optimal) solution */
      memcpy(&x[1], &x_best[1], n * sizeof(char));
      return z_best;
}

/***********************************************************************
*  reduce - prepare reduced instance of 0-1 knapsack
*
*  Given original instance of 0-1 knapsack (1)-(3) specified by the
*  parameters n, a, b, and c this routine transforms it to equivalent
*  reduced instance in the same format. The reduced instance is
*  normalized, i.e. the following additional conditions are met:
*
*     n >= 2                                                         (4)
*
*     1 <= a[j] <= b for all j in 1..n                               (5)
*
*     sum{j in 1..n} a[j] >= b+1                                     (6)
*
*     c[j] >= 1      for all j in 1..n                               (7)
*
*  The routine creates the structure ks and stores there parameters n,
*  a, b, and c of the reduced instance as well as template of solution
*  to original instance.
*
*  Normally the routine returns a pointer to the structure ks created.
*  However, if the original instance is infeasible, the routine returns
*  a null pointer. */

struct ks
{     int orig_n;
      /* original problem dimension */
      int n;
      /* reduced problem dimension */
      int *a; /* int a[1+orig_n]; */
      /* a{j in 1..n} are constraint coefficients (2) */
      int b;
      /* b is constraint right-hand side (2) */
      int *c; /* int c[1+orig_n]; */
      /* c{j in 1..n} are objective coefficients (1) */
      int c0;
      /* c0 is objective constant term */
      char *x; /* char x[1+orig_n]; */
      /* x{j in 1..orig_n} is solution template to original instance:
       * x[j] = 0       x[j] is fixed at 0
       * x[j] = 1       x[j] is fixed at 1
       * x[j] = 0x10    x[j] = x[j']
       * x[j] = 0x11    x[j] = 1 - x[j']
       * where x[j'] is corresponding solution to reduced instance */
};

static void free_ks(struct ks *ks);

static struct ks *reduce(const int n, const int a[/*1+n*/], int b,
      const int c[/*1+n*/])
{     struct ks *ks;
      int j, s;
      xassert(n >= 0);
      /* initially reduced instance is the same as original one */
      ks = talloc(1, struct ks);
      ks->orig_n = n;
      ks->n = 0;
      ks->a = talloc(1+n, int);
      memcpy(&ks->a[1], &a[1], n * sizeof(int));
      ks->b = b;
      ks->c = talloc(1+n, int);
      memcpy(&ks->c[1], &c[1], n * sizeof(int));
      ks->c0 = 0;
      ks->x = talloc(1+n, char);
      /* make all a[j] non-negative */
      for (j = 1; j <= n; j++)
      {  if (a[j] >= 0)
         {  /* keep original x[j] */
            ks->x[j] = 0x10;
         }
         else /* a[j] < 0 */
         {  /* substitute x[j] = 1 - x'[j] */
            ks->x[j] = 0x11;
            /* ... + a[j]x[j]        + ... <= b
             * ... + a[j](1 - x'[j]) + ... <= b
             * ... - a[j]x'[j]       + ... <= b - a[j] */
            ks->a[j] = - ks->a[j];
            ks->b += ks->a[j];
            /* z = ... + c[j]x[j]        + ... + c0 =
             *   = ... + c[j](1 - x'[j]) + ... + c0 =
             *   = ... - c[j]x'[j]       + ... + (c0 + c[j]) */
            ks->c0 += ks->c[j];
            ks->c[j] = - ks->c[j];
         }
      }
      /* now a[j] >= 0 for all j in 1..n */
      if (ks->b < 0)
      {  /* instance is infeasible */
         free_ks(ks);
         return NULL;
      }
      /* build reduced instance */
      for (j = 1; j <= n; j++)
      {  if (ks->a[j] == 0)
         {  if (ks->c[j] <= 0)
            {  /* fix x[j] at 0 */
               ks->x[j] ^= 0x10;
            }
            else
            {  /* fix x[j] at 1 */
               ks->x[j] ^= 0x11;
               ks->c0 += ks->c[j];
            }
         }
         else if (ks->a[j] > ks->b || ks->c[j] <= 0)
         {  /* fix x[j] at 0 */
            ks->x[j] ^= 0x10;
         }
         else
         {  /* include x[j] in reduced instance */
            ks->n++;
            ks->a[ks->n] = ks->a[j];
            ks->c[ks->n] = ks->c[j];
         }
      }
      /* now conditions (5) and (7) are met */
      /* check condition (6) */
      s = 0;
      for (j = 1; j <= ks->n; j++)
      {  xassert(1 <= ks->a[j] && ks->a[j] <= ks->b);
         xassert(ks->c[j] >= 1);
         s += ks->a[j];
      }
      if (s <= ks->b)
      {  /* sum{j in 1..n} a[j] <= b */
         /* fix all remaining x[j] at 1 to obtain trivial solution */
         for (j = 1; j <= n; j++)
         {  if (ks->x[j] & 0x10)
               ks->x[j] ^= 0x11;
         }
         for (j = 1; j <= ks->n; j++)
            ks->c0 += ks->c[j];
         /* reduced instance is empty */
         ks->n = 0;
      }
      /* here n = 0 or n >= 2 due to condition (6) */
      xassert(ks->n == 0 || ks->n >= 2);
      return ks;
}

/***********************************************************************
*  restore - restore solution to original 0-1 knapsack instance
*
*  Given optimal solution x{j in 1..ks->n} to the reduced 0-1 knapsack
*  instance (previously prepared by the routine reduce) this routine
*  constructs optimal solution to the original instance and stores it
*  in the array ks->x{j in 1..ks->orig_n}.
*
*  On exit the routine returns optimal objective value for the original
*  instance.
*
*  NOTE: This operation should be performed only once. */

static int restore(struct ks *ks, char x[])
{     int j, k, z;
      z = ks->c0;
      for (j = 1, k = 0; j <= ks->orig_n; j++)
      {  if (ks->x[j] & 0x10)
         {  k++;
            xassert(k <= ks->n);
            xassert(x[k] == 0 || x[k] == 1);
            if (ks->x[j] & 1)
               ks->x[j] = 1 - x[k];
            else
               ks->x[j] = x[k];
            if (x[k])
               z += ks->c[k];
         }
      }
      xassert(k == ks->n);
      return z;
}

/***********************************************************************
*  free_ks - deallocate structure ks
*
*  This routine frees memory previously allocated to the structure ks
*  and all its components. */

static void free_ks(struct ks *ks)
{     xassert(ks != NULL);
      tfree(ks->a);
      tfree(ks->c);
      tfree(ks->x);
      tfree(ks);
}

/***********************************************************************
*  ks_mt1 - solve 0-1 knapsack problem with Martello & Toth algorithm
*
*  This routine finds optimal solution to 0-1 knapsack problem (1)-(3)
*  with Martello & Toth algorithm MT1.
*
*  The instance to be solved is specified by parameters n, a, b, and c.
*  Note that these parameters can have any sign, i.e. normalization is
*  not needed.
*
*  On exit the routine stores the optimal point found in locations
*  x[1], ..., x[n] and returns the optimal objective value. However, if
*  the instance is infeasible, the routine returns INT_MIN.
*
*  REFERENCES
*
*  S.Martello, P.Toth. Knapsack Problems: Algorithms and Computer Imp-
*  lementations. John Wiley & Sons, 1990. */

struct mt
{     int j;
      float r; /* r[j] = c[j] / a[j] */
};

static int CDECL fcmp(const void *p1, const void *p2)
{     if (((struct mt *)p1)->r > ((struct mt *)p2)->r)
         return -1;
      else if (((struct mt *)p1)->r < ((struct mt *)p2)->r)
         return +1;
      else
         return 0;
}

static int mt1a(int n, const int a[], int b, const int c[], char x[])
{     /* interface routine to MT1 */
      struct mt *mt;
      int j, z, *p, *w, *x1, *xx, *min, *psign, *wsign, *zsign;
      xassert(n >= 2);
      /* allocate working arrays */
      mt = talloc(1+n, struct mt);
      p = talloc(1+n+1, int);
      w = talloc(1+n+1, int);
      x1 = talloc(1+n+1, int);
      xx = talloc(1+n+1, int);
      min = talloc(1+n+1, int);
      psign = talloc(1+n+1, int);
      wsign = talloc(1+n+1, int);
      zsign = talloc(1+n+1, int);
      /* reorder items to provide c[j] / a[j] >= a[j+1] / a[j+1] */
      for (j = 1; j <= n; j++)
      {  mt[j].j = j;
         mt[j].r = (float)c[j] / (float)a[j];
      }
      qsort(&mt[1], n, sizeof(struct mt), fcmp);
      /* load instance parameters */
      for (j = 1; j <= n; j++)
      {  p[j] = c[mt[j].j];
         w[j] = a[mt[j].j];
      }
      /* find optimal solution */
      z = mt1(n, p, w, b, x1, 1, xx, min, psign, wsign, zsign);
      xassert(z >= 0);
      /* store optimal point found */
      for (j = 1; j <= n; j++)
      {  xassert(x1[j] == 0 || x1[j] == 1);
         x[mt[j].j] = x1[j];
      }
      /* free working arrays */
      tfree(mt);
      tfree(p);
      tfree(w);
      tfree(x1);
      tfree(xx);
      tfree(min);
      tfree(psign);
      tfree(wsign);
      tfree(zsign);
      return z;
}

int ks_mt1(int n, const int a[/*1+n*/], int b, const int c[/*1+n*/],
      char x[/*1+n*/])
{     struct ks *ks;
      int j, s1, s2, z;
      xassert(n >= 0);
      /* prepare reduced instance */
      ks = reduce(n, a, b, c);
      if (ks == NULL)
      {  /* original instance is infeasible */
         return INT_MIN;
      }
      /* find optimal solution to reduced instance */
      if (ks->n > 0)
         mt1a(ks->n, ks->a, ks->b, ks->c, x);
      /* restore solution to original instance */
      z = restore(ks, x);
      memcpy(&x[1], &ks->x[1], n * sizeof(char));
      free_ks(ks);
      /* check solution found */
      s1 = s2 = 0;
      for (j = 1; j <= n; j++)
      {  xassert(x[j] == 0 || x[j] == 1);
         if (x[j])
            s1 += a[j], s2 += c[j];
      }
      xassert(s1 <= b);
      xassert(s2 == z);
      return z;
}

/***********************************************************************
*  ks_greedy - solve 0-1 knapsack problem with greedy heuristic
*
*  This routine finds (sub)optimal solution to 0-1 knapsack problem
*  (1)-(3) with greedy heuristic.
*
*  The instance to be solved is specified by parameters n, a, b, and c.
*  Note that these parameters can have any sign, i.e. normalization is
*  not needed.
*
*  On exit the routine stores the optimal point found in locations
*  x[1], ..., x[n] and returns the optimal objective value. However, if
*  the instance is infeasible, the routine returns INT_MIN. */

static int greedy(int n, const int a[], int b, const int c[], char x[])
{     /* core routine for normalized 0-1 knapsack instance */
      struct mt *mt;
      int j, s, z;
      xassert(n >= 2);
      /* reorder items to provide c[j] / a[j] >= a[j+1] / a[j+1] */
      mt = talloc(1+n, struct mt);
      for (j = 1; j <= n; j++)
      {  mt[j].j = j;
         mt[j].r = (float)c[j] / (float)a[j];
      }
      qsort(&mt[1], n, sizeof(struct mt), fcmp);
      /* take items starting from most valuable ones until the knapsack
       * is full */
      s = z = 0;
      for (j = 1; j <= n; j++)
      {  if (s + a[mt[j].j] > b)
            break;
         x[mt[j].j] = 1;
         s += a[mt[j].j];
         z += c[mt[j].j];
      }
      /* don't take remaining items */
      for (j = j; j <= n; j++)
         x[mt[j].j] = 0;
      tfree(mt);
      return z;
}

int ks_greedy(int n, const int a[/*1+n*/], int b, const int c[/*1+n*/],
      char x[/*1+n*/])
{     struct ks *ks;
      int j, s1, s2, z;
      xassert(n >= 0);
      /* prepare reduced instance */
      ks = reduce(n, a, b, c);
      if (ks == NULL)
      {  /* original instance is infeasible */
         return INT_MIN;
      }
      /* find suboptimal solution to reduced instance */
      if (ks->n > 0)
         greedy(ks->n, ks->a, ks->b, ks->c, x);
      /* restore solution to original instance */
      z = restore(ks, x);
      memcpy(&x[1], &ks->x[1], n * sizeof(char));
      free_ks(ks);
      /* check solution found */
      s1 = s2 = 0;
      for (j = 1; j <= n; j++)
      {  xassert(x[j] == 0 || x[j] == 1);
         if (x[j])
            s1 += a[j], s2 += c[j];
      }
      xassert(s1 <= b);
      xassert(s2 == z);
      return z;
}

/* eof */
