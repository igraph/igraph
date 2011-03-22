/* glpnpp04.c */

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

#include "glpnpp.h"

/***********************************************************************
*  NAME
*
*  npp_binarize_prob - binarize MIP problem
*
*  SYNOPSIS
*
*  #include "glpnpp.h"
*  int npp_binarize_prob(NPP *npp);
*
*  DESCRIPTION
*
*  The routine npp_binarize_prob replaces in the original MIP problem
*  every integer variable:
*
*     l[q] <= x[q] <= u[q],                                          (1)
*
*  where l[q] < u[q], by an equivalent sum of binary variables.
*
*  RETURNS
*
*  The routine returns the number of integer variables for which the
*  transformation failed, because u[q] - l[q] > d_max.
*
*  PROBLEM TRANSFORMATION
*
*  If variable x[q] has non-zero lower bound, it is first processed
*  with the routine npp_lbnd_col. Thus, we can assume that:
*
*     0 <= x[q] <= u[q].                                             (2)
*
*  If u[q] = 1, variable x[q] is already binary, so further processing
*  is not needed. Let, therefore, that 2 <= u[q] <= d_max, and n be a
*  smallest integer such that u[q] <= 2^n - 1 (n >= 2, since u[q] >= 2).
*  Then variable x[q] can be replaced by the following sum:
*
*            n-1
*     x[q] = sum 2^k x[k],                                           (3)
*            k=0
*
*  where x[k] are binary columns (variables). If u[q] < 2^n - 1, the
*  following additional inequality constraint must be also included in
*  the transformed problem:
*
*     n-1
*     sum 2^k x[k] <= u[q].                                          (4)
*     k=0
*
*  Note: Assuming that in the transformed problem x[q] becomes binary
*  variable x[0], this transformation causes new n-1 binary variables
*  to appear.
*
*  Substituting x[q] from (3) to the objective row gives:
*
*     z = sum c[j] x[j] + c[0] =
*          j
*
*       = sum c[j] x[j] + c[q] x[q] + c[0] =
*         j!=q
*                              n-1
*       = sum c[j] x[j] + c[q] sum 2^k x[k] + c[0] =
*         j!=q                 k=0
*                         n-1
*       = sum c[j] x[j] + sum c[k] x[k] + c[0],
*         j!=q            k=0
*
*  where:
*
*     c[k] = 2^k c[q],  k = 0, ..., n-1.                             (5)
*
*  And substituting x[q] from (3) to i-th constraint row i gives:
*
*     L[i] <= sum a[i,j] x[j] <= U[i]  ==>
*              j
*
*     L[i] <= sum a[i,j] x[j] + a[i,q] x[q] <= U[i]  ==>
*             j!=q
*                                      n-1
*     L[i] <= sum a[i,j] x[j] + a[i,q] sum 2^k x[k] <= U[i]  ==>
*             j!=q                     k=0
*                               n-1
*     L[i] <= sum a[i,j] x[j] + sum a[i,k] x[k] <= U[i],
*             j!=q              k=0
*
*  where:
*
*     a[i,k] = 2^k a[i,q],  k = 0, ..., n-1.                         (6)
*
*  RECOVERING SOLUTION
*
*  Value of variable x[q] is computed with formula (3). */

struct binarize
{     int q;
      /* column reference number for x[q] = x[0] */
      int j;
      /* column reference number for x[1]; x[2] has reference number
         j+1, x[3] - j+2, etc. */
      int n;
      /* total number of binary variables, n >= 2 */
};

static int rcv_binarize_prob(NPP *npp, void *info);

int npp_binarize_prob(NPP *npp)
{     /* binarize MIP problem */
      struct binarize *info;
      NPPROW *row;
      NPPCOL *col, *bin;
      NPPAIJ *aij;
      int u, n, k, temp, nfails, nvars, nbins, nrows;
      /* new variables will be added to the end of the column list, so
         we go from the end to beginning of the column list */
      nfails = nvars = nbins = nrows = 0;
      for (col = npp->c_tail; col != NULL; col = col->prev)
      {  /* skip continuous variable */
         if (!col->is_int) continue;
         /* skip fixed variable */
         if (col->lb == col->ub) continue;
         /* skip binary variable */
         if (col->lb == 0.0 && col->ub == 1.0) continue;
         /* check if the transformation is applicable */
         if (col->lb < -1e6 || col->ub > +1e6 ||
             col->ub - col->lb > 4095.0)
         {  /* unfortunately, not */
            nfails++;
            continue;
         }
         /* process integer non-binary variable x[q] */
         nvars++;
         /* make x[q] non-negative, if its lower bound is non-zero */
         if (col->lb != 0.0)
            npp_lbnd_col(npp, col);
         /* now 0 <= x[q] <= u[q] */
         xassert(col->lb == 0.0);
         u = (int)col->ub;
         xassert(col->ub == (double)u);
         /* if x[q] is binary, further processing is not needed */
         if (u == 1) continue;
         /* determine smallest n such that u <= 2^n - 1 (thus, n is the
            number of binary variables needed) */
         n = 2, temp = 4;
         while (u >= temp)
            n++, temp += temp;
         nbins += n;
         /* create transformation stack entry */
         info = npp_push_tse(npp,
            rcv_binarize_prob, sizeof(struct binarize));
         info->q = col->j;
         info->j = 0; /* will be set below */
         info->n = n;
         /* if u < 2^n - 1, we need one additional row for (4) */
         if (u < temp - 1)
         {  row = npp_add_row(npp), nrows++;
            row->lb = -DBL_MAX, row->ub = u;
         }
         else
            row = NULL;
         /* in the transformed problem variable x[q] becomes binary
            variable x[0], so its objective and constraint coefficients
            are not changed */
         col->ub = 1.0;
         /* include x[0] into constraint (4) */
         if (row != NULL)
            npp_add_aij(npp, row, col, 1.0);
         /* add other binary variables x[1], ..., x[n-1] */
         for (k = 1, temp = 2; k < n; k++, temp += temp)
         {  /* add new binary variable x[k] */
            bin = npp_add_col(npp);
            bin->is_int = 1;
            bin->lb = 0.0, bin->ub = 1.0;
            bin->coef = (double)temp * col->coef;
            /* store column reference number for x[1] */
            if (info->j == 0)
               info->j = bin->j;
            else
               xassert(info->j + (k-1) == bin->j);
            /* duplicate constraint coefficients for x[k]; this also
               automatically includes x[k] into constraint (4) */
            for (aij = col->ptr; aij != NULL; aij = aij->c_next)
               npp_add_aij(npp, aij->row, bin, (double)temp * aij->val);
         }
      }
      if (nvars > 0)
         xprintf("%d integer variable(s) were replaced by %d binary one"
            "s\n", nvars, nbins);
      if (nrows > 0)
         xprintf("%d row(s) were added due to binarization\n", nrows);
      if (nfails > 0)
         xprintf("Binarization failed for %d integer variable(s)\n",
            nfails);
      return nfails;
}

static int rcv_binarize_prob(NPP *npp, void *_info)
{     /* recovery binarized variable */
      struct binarize *info = _info;
      int k, temp;
      double sum;
      /* compute value of x[q]; see formula (3) */
      sum = npp->c_value[info->q];
      for (k = 1, temp = 2; k < info->n; k++, temp += temp)
         sum += (double)temp * npp->c_value[info->j + (k-1)];
      npp->c_value[info->q] = sum;
      return 0;
}

/**********************************************************************/

struct elem
{     /* linear form element a[j] x[j] */
      double aj;
      /* non-zero coefficient value */
      NPPCOL *xj;
      /* pointer to variable (column) */
      struct elem *next;
      /* pointer to another term */
};

static struct elem *copy_form(NPP *npp, NPPROW *row, double s)
{     /* copy linear form */
      NPPAIJ *aij;
      struct elem *ptr, *e;
      ptr = NULL;
      for (aij = row->ptr; aij != NULL; aij = aij->r_next)
      {  e = dmp_get_atom(npp->pool, sizeof(struct elem));
         e->aj = s * aij->val;
         e->xj = aij->col;
         e->next = ptr;
         ptr = e;
      }
      return ptr;
}

static void drop_form(NPP *npp, struct elem *ptr)
{     /* drop linear form */
      struct elem *e;
      while (ptr != NULL)
      {  e = ptr;
         ptr = e->next;
         dmp_free_atom(npp->pool, e, sizeof(struct elem));
      }
      return;
}

/***********************************************************************
*  NAME
*
*  npp_is_packing - test if constraint is packing inequality
*
*  SYNOPSIS
*
*  #include "glpnpp.h"
*  int npp_is_packing(NPP *npp, NPPROW *row);
*
*  RETURNS
*
*  If the specified row (constraint) is packing inequality (see below),
*  the routine npp_is_packing returns non-zero. Otherwise, it returns
*  zero.
*
*  PACKING INEQUALITIES
*
*  In canonical format the packing inequality is the following:
*
*     sum  x[j] <= 1,                                                (1)
*    j in J
*
*  where all variables x[j] are binary. This inequality expresses the
*  condition that in any integer feasible solution at most one variable
*  from set J can take non-zero (unity) value while other variables
*  must be equal to zero. W.l.o.g. it is assumed that |J| >= 2, because
*  if J is empty or |J| = 1, the inequality (1) is redundant.
*
*  In general case the packing inequality may include original variables
*  x[j] as well as their complements x~[j]:
*
*     sum   x[j] + sum   x~[j] <= 1,                                 (2)
*    j in Jp      j in Jn
*
*  where Jp and Jn are not intersected. Therefore, using substitution
*  x~[j] = 1 - x[j] gives the packing inequality in generalized format:
*
*     sum   x[j] - sum   x[j] <= 1 - |Jn|.                           (3)
*    j in Jp      j in Jn */

int npp_is_packing(NPP *npp, NPPROW *row)
{     /* test if constraint is packing inequality */
      NPPCOL *col;
      NPPAIJ *aij;
      int b;
      xassert(npp == npp);
      if (!(row->lb == -DBL_MAX && row->ub != +DBL_MAX))
         return 0;
      b = 1;
      for (aij = row->ptr; aij != NULL; aij = aij->r_next)
      {  col = aij->col;
         if (!(col->is_int && col->lb == 0.0 && col->ub == 1.0))
            return 0;
         if (aij->val == +1.0)
            ;
         else if (aij->val == -1.0)
            b--;
         else
            return 0;
      }
      if (row->ub != (double)b) return 0;
      return 1;
}

/***********************************************************************
*  NAME
*
*  npp_hidden_packing - identify hidden packing inequality
*
*  SYNOPSIS
*
*  #include "glpnpp.h"
*  int npp_hidden_packing(NPP *npp, NPPROW *row);
*
*  DESCRIPTION
*
*  The routine npp_hidden_packing processes specified inequality
*  constraint, which includes only binary variables, and the number of
*  the variables is not less than two. If the original inequality is
*  equivalent to a packing inequality, the routine replaces it by this
*  equivalent inequality. If the original constraint is double-sided
*  inequality, it is replaced by a pair of single-sided inequalities,
*  if necessary.
*
*  RETURNS
*
*  If the original inequality constraint was replaced by equivalent
*  packing inequality, the routine npp_hidden_packing returns non-zero.
*  Otherwise, it returns zero.
*
*  PROBLEM TRANSFORMATION
*
*  Consider an inequality constraint:
*
*     sum  a[j] x[j] <= b,                                           (1)
*    j in J
*
*  where all variables x[j] are binary, and |J| >= 2. (In case of '>='
*  inequality it can be transformed to '<=' format by multiplying both
*  its sides by -1.)
*
*  Let Jp = {j: a[j] > 0}, Jn = {j: a[j] < 0}. Performing substitution
*  x[j] = 1 - x~[j] for all j in Jn, we have:
*
*     sum   a[j] x[j] <= b  ==>
*    j in J
*
*     sum   a[j] x[j] + sum   a[j] x[j] <= b  ==>
*    j in Jp           j in Jn
*
*     sum   a[j] x[j] + sum   a[j] (1 - x~[j]) <= b  ==>
*    j in Jp           j in Jn
*
*     sum   a[j] x[j] - sum   a[j] x~[j] <= b - sum   a[j].
*    j in Jp           j in Jn                 j in Jn
*
*  Thus, meaning the transformation above, we can assume that in
*  inequality (1) all coefficients a[j] are positive. Moreover, we can
*  assume that a[j] <= b. In fact, let a[j] > b; then the following
*  three cases are possible:
*
*  1) b < 0. In this case inequality (1) is infeasible, so the problem
*     has no feasible solution (see the routine npp_analyze_row);
*
*  2) b = 0. In this case inequality (1) is a forcing inequality on its
*     upper bound (see the routine npp_forcing row), from which it
*     follows that all variables x[j] should be fixed at zero;
*
*  3) b > 0. In this case inequality (1) defines an implied zero upper
*     bound for variable x[j] (see the routine npp_implied_bounds), from
*     which it follows that x[j] should be fixed at zero.
*
*  It is assumed that all three cases listed above have been recognized
*  by the routine npp_process_prob, which performs basic MIP processing
*  prior to a call the routine npp_hidden_packing. So, if one of these
*  cases occurs, we should just skip processing such constraint.
*
*  Thus, let 0 < a[j] <= b. Then it is obvious that constraint (1) is
*  equivalent to packing inquality only if:
*
*     a[j] + a[k] > b + eps                                          (2)
*
*  for all j, k in J, j != k, where eps is an absolute tolerance for
*  row (linear form) value. Checking the condition (2) for all j and k,
*  j != k, requires time O(|J|^2). However, this time can be reduced to
*  O(|J|), if use minimal a[j] and a[k], in which case it is sufficient
*  to check the condition (2) only once.
*
*  Once the original inequality (1) is replaced by equivalent packing
*  inequality, we need to perform back substitution x~[j] = 1 - x[j] for
*  all j in Jn (see above).
*
*  RECOVERING SOLUTION
*
*  None needed. */

static int hidden_packing(NPP *npp, struct elem *ptr, double *_b)
{     /* process inequality constraint: sum a[j] x[j] <= b;
         0 - specified row is NOT hidden packing inequality;
         1 - specified row is packing inequality;
         2 - specified row is hidden packing inequality. */
      struct elem *e, *ej, *ek;
      int neg;
      double b = *_b, eps;
      xassert(npp == npp);
      /* a[j] must be non-zero, x[j] must be binary, for all j in J */
      for (e = ptr; e != NULL; e = e->next)
      {  xassert(e->aj != 0.0);
         xassert(e->xj->is_int);
         xassert(e->xj->lb == 0.0 && e->xj->ub == 1.0);
      }
      /* check if the specified inequality constraint already has the
         form of packing inequality */
      neg = 0; /* neg is |Jn| */
      for (e = ptr; e != NULL; e = e->next)
      {  if (e->aj == +1.0)
            ;
         else if (e->aj == -1.0)
            neg++;
         else
            break;
      }
      if (e == NULL)
      {  /* all coefficients a[j] are +1 or -1; check rhs b */
         if (b == (double)(1 - neg))
         {  /* it is packing inequality; no processing is needed */
            return 1;
         }
      }
      /* substitute x[j] = 1 - x~[j] for all j in Jn to make all a[j]
         positive; the result is a~[j] = |a[j]| and new rhs b */
      for (e = ptr; e != NULL; e = e->next)
         if (e->aj < 0) b -= e->aj;
      /* now a[j] > 0 for all j in J (actually |a[j]| are used) */
      /* if a[j] > b, skip processing--this case must not appear */
      for (e = ptr; e != NULL; e = e->next)
         if (fabs(e->aj) > b) return 0;
      /* now 0 < a[j] <= b for all j in J */
      /* find two minimal coefficients a[j] and a[k], j != k */
      ej = NULL;
      for (e = ptr; e != NULL; e = e->next)
         if (ej == NULL || fabs(ej->aj) > fabs(e->aj)) ej = e;
      xassert(ej != NULL);
      ek = NULL;
      for (e = ptr; e != NULL; e = e->next)
         if (e != ej)
            if (ek == NULL || fabs(ek->aj) > fabs(e->aj)) ek = e;
      xassert(ek != NULL);
      /* the specified constraint is equivalent to packing inequality
         iff a[j] + a[k] > b + eps */
      eps = 1e-3 + 1e-6 * fabs(b);
      if (fabs(ej->aj) + fabs(ek->aj) <= b + eps) return 0;
      /* perform back substitution x~[j] = 1 - x[j] and construct the
         final equivalent packing inequality in generalized format */
      b = 1.0;
      for (e = ptr; e != NULL; e = e->next)
      {  if (e->aj > 0.0)
            e->aj = +1.0;
         else /* e->aj < 0.0 */
            e->aj = -1.0, b -= 1.0;
      }
      *_b = b;
      return 2;
}

int npp_hidden_packing(NPP *npp, NPPROW *row)
{     /* identify hidden packing inequality */
      NPPROW *copy;
      NPPAIJ *aij;
      struct elem *ptr, *e;
      int kase, ret, count = 0;
      double b;
      /* the row must be inequality constraint */
      xassert(row->lb < row->ub);
      for (kase = 0; kase <= 1; kase++)
      {  if (kase == 0)
         {  /* process row upper bound */
            if (row->ub == +DBL_MAX) continue;
            ptr = copy_form(npp, row, +1.0);
            b = + row->ub;
         }
         else
         {  /* process row lower bound */
            if (row->lb == -DBL_MAX) continue;
            ptr = copy_form(npp, row, -1.0);
            b = - row->lb;
         }
         /* now the inequality has the form "sum a[j] x[j] <= b" */
         ret = hidden_packing(npp, ptr, &b);
         xassert(0 <= ret && ret <= 2);
         if (kase == 1 && ret == 1 || ret == 2)
         {  /* the original inequality has been identified as hidden
               packing inequality */
            count++;
#ifdef GLP_DEBUG
            xprintf("Original constraint:\n");
            for (aij = row->ptr; aij != NULL; aij = aij->r_next)
               xprintf(" %+g x%d", aij->val, aij->col->j);
            if (row->lb != -DBL_MAX) xprintf(", >= %g", row->lb);
            if (row->ub != +DBL_MAX) xprintf(", <= %g", row->ub);
            xprintf("\n");
            xprintf("Equivalent packing inequality:\n");
            for (e = ptr; e != NULL; e = e->next)
               xprintf(" %sx%d", e->aj > 0.0 ? "+" : "-", e->xj->j);
            xprintf(", <= %g\n", b);
#endif
            if (row->lb == -DBL_MAX || row->ub == +DBL_MAX)
            {  /* the original row is single-sided inequality; no copy
                  is needed */
               copy = NULL;
            }
            else
            {  /* the original row is double-sided inequality; we need
                  to create its copy for other bound before replacing it
                  with the equivalent inequality */
               copy = npp_add_row(npp);
               if (kase == 0)
               {  /* the copy is for lower bound */
                  copy->lb = row->lb, copy->ub = +DBL_MAX;
               }
               else
               {  /* the copy is for upper bound */
                  copy->lb = -DBL_MAX, copy->ub = row->ub;
               }
               /* copy original row coefficients */
               for (aij = row->ptr; aij != NULL; aij = aij->r_next)
                  npp_add_aij(npp, copy, aij->col, aij->val);
            }
            /* replace the original inequality by equivalent one */
            npp_erase_row(npp, row);
            row->lb = -DBL_MAX, row->ub = b;
            for (e = ptr; e != NULL; e = e->next)
               npp_add_aij(npp, row, e->xj, e->aj);
            /* continue processing lower bound for the copy */
            if (copy != NULL) row = copy;
         }
         drop_form(npp, ptr);
      }
      return count;
}

/***********************************************************************
*  NAME
*
*  npp_implied_packing - identify implied packing inequality
*
*  SYNOPSIS
*
*  #include "glpnpp.h"
*  int npp_implied_packing(NPP *npp, NPPROW *row, int which,
*     NPPCOL *var[], char set[]);
*
*  DESCRIPTION
*
*  The routine npp_implied_packing processes specified row (constraint)
*  of general format:
*
*     L <= sum a[j] x[j] <= U.                                       (1)
*           j
*
*  If which = 0, only lower bound L, which must exist, is considered,
*  while upper bound U is ignored. Similarly, if which = 1, only upper
*  bound U, which must exist, is considered, while lower bound L is
*  ignored. Thus, if the specified row is a double-sided inequality or
*  equality constraint, this routine should be called twice for both
*  lower and upper bounds.
*
*  The routine npp_implied_packing attempts to find a non-trivial (i.e.
*  having not less than two binary variables) packing inequality:
*
*     sum   x[j] - sum   x[j] <= 1 - |Jn|,                           (2)
*    j in Jp      j in Jn
*
*  which is relaxation of the constraint (1) in the sense that any
*  solution satisfying to that constraint also satisfies to the packing
*  inequality (2). If such relaxation exists, the routine stores
*  pointers to descriptors of corresponding binary variables and their
*  flags, resp., to locations var[1], var[2], ..., var[len] and set[1],
*  set[2], ..., set[len], where set[j] = 0 means that j in Jp and
*  set[j] = 1 means that j in Jn.
*
*  RETURNS
*
*  The routine npp_implied_packing returns len, which is the total
*  number of binary variables in the packing inequality found, len >= 2.
*  However, if the relaxation does not exist, the routine returns zero.
*
*  ALGORITHM
*
*  If which = 0, the constraint coefficients (1) are multiplied by -1
*  and b is assigned -L; if which = 1, the constraint coefficients (1)
*  are not changed and b is assigned +U. In both cases the specified
*  constraint gets the following format:
*
*     sum a[j] x[j] <= b.                                            (3)
*      j
*
*  (Note that (3) is a relaxation of (1), because one of bounds L or U
*  is ignored.)
*
*  Let J be set of binary variables, Kp be set of non-binary (integer
*  or continuous) variables with a[j] > 0, and Kn be set of non-binary
*  variables with a[j] < 0. Then the inequality (3) can be written as
*  follows:
*
*     sum  a[j] x[j] <= b - sum   a[j] x[j] - sum   a[j] x[j].       (4)
*    j in J                j in Kp           j in Kn
*
*  To get rid of non-binary variables we can replace the inequality (4)
*  by the following relaxed inequality:
*
*     sum  a[j] x[j] <= b~,                                          (5)
*    j in J
*
*  where:
*
*     b~ = sup(b - sum   a[j] x[j] - sum   a[j] x[j]) =
*                 j in Kp           j in Kn
*
*        = b - inf sum   a[j] x[j] - inf sum   a[j] x[j] =           (6)
*                 j in Kp               j in Kn
*
*        = b - sum   a[j] l[j] - sum   a[j] u[j].
*             j in Kp           j in Kn
*
*  Note that if lower bound l[j] (if j in Kp) or upper bound u[j]
*  (if j in Kn) of some non-binary variable x[j] does not exist, then
*  formally b = +oo, in which case further analysis is not performed.
*
*  Let Bp = {j in J: a[j] > 0}, Bn = {j in J: a[j] < 0}. To make all
*  the inequality coefficients in (5) positive, we replace all x[j] in
*  Bn by their complementaries, substituting x[j] = 1 - x~[j] for all
*  j in Bn, that gives:
*
*     sum   a[j] x[j] - sum   a[j] x~[j] <= b~ - sum   a[j].         (7)
*    j in Bp           j in Bn                  j in Bn
*
*  This inequality is a relaxation of the original constraint (1), and
*  it is a binary knapsack inequality. Writing it in the standard format
*  we have:
*
*     sum  alfa[j] z[j] <= beta,                                     (8)
*    j in J
*
*  where:
*               ( + a[j],   if j in Bp,
*     alfa[j] = <                                                    (9)
*               ( - a[j],   if j in Bn,
*
*               ( x[j],     if j in Bp,
*        z[j] = <                                                   (10)
*               ( 1 - x[j], if j in Bn,
*
*        beta = b~ - sum   a[j].                                    (11)
*                   j in Bn
*
*  In the inequality (8) all coefficients are positive, therefore, the
*  packing relaxation to be found for this inequality is the following:
*
*     sum  z[j] <= 1.                                               (12)
*    j in P
*
*  It is obvious that set P within J, which we would like to find, must
*  satisfy to the following condition:
*
*     alfa[j] + alfa[k] > beta + eps  for all j, k in P, j != k,    (13)
*
*  where eps is an absolute tolerance for value of the linear form.
*  Thus, it is natural to take P = {j: alpha[j] > (beta + eps) / 2}.
*  Moreover, if in the equality (8) there exist coefficients alfa[k],
*  for which alfa[k] <= (beta + eps) / 2, but which, nevertheless,
*  satisfies to the condition (13) for all j in P, *one* corresponding
*  variable z[k] (having, for example, maximal coefficient alfa[k]) can
*  be included in set P, that allows increasing the number of binary
*  variables in (12) by one.
*
*  Once the set P has been built, for the inequality (12) we need to
*  perform back substitution according to (10) in order to express it
*  through the original binary variables. As the result of such back
*  substitution the relaxed packing inequality get its final format (2),
*  where Jp = J intersect Bp, and Jn = J intersect Bn. */

int npp_implied_packing(NPP *npp, NPPROW *row, int which,
      NPPCOL *var[], char set[])
{     struct elem *ptr, *e, *i, *k;
      int len = 0;
      double b, eps;
      /* build inequality (3) */
      if (which == 0)
      {  ptr = copy_form(npp, row, -1.0);
         xassert(row->lb != -DBL_MAX);
         b = - row->lb;
      }
      else if (which == 1)
      {  ptr = copy_form(npp, row, +1.0);
         xassert(row->ub != +DBL_MAX);
         b = + row->ub;
      }
      /* remove non-binary variables to build relaxed inequality (5);
         compute its right-hand side b~ with formula (6) */
      for (e = ptr; e != NULL; e = e->next)
      {  if (!(e->xj->is_int && e->xj->lb == 0.0 && e->xj->ub == 1.0))
         {  /* x[j] is non-binary variable */
            if (e->aj > 0.0)
            {  if (e->xj->lb == -DBL_MAX) goto done;
               b -= e->aj * e->xj->lb;
            }
            else /* e->aj < 0.0 */
            {  if (e->xj->ub == +DBL_MAX) goto done;
               b -= e->aj * e->xj->ub;
            }
            /* a[j] = 0 means that variable x[j] is removed */
            e->aj = 0.0;
         }
      }
      /* substitute x[j] = 1 - x~[j] to build knapsack inequality (8);
         compute its right-hand side beta with formula (11) */
      for (e = ptr; e != NULL; e = e->next)
         if (e->aj < 0.0) b -= e->aj;
      /* if beta is close to zero, the knapsack inequality is either
         infeasible or forcing inequality; this must never happen, so
         we skip further analysis */
      if (b < 1e-3) goto done;
      /* build set P as well as sets Jp and Jn, and determine x[k] as
         explained above in comments to the routine */
      eps = 1e-3 + 1e-6 * b;
      i = k = NULL;
      for (e = ptr; e != NULL; e = e->next)
      {  /* note that alfa[j] = |a[j]| */
         if (fabs(e->aj) > 0.5 * (b + eps))
         {  /* alfa[j] > (b + eps) / 2; include x[j] in set P, i.e. in
               set Jp or Jn */
            var[++len] = e->xj;
            set[len] = (char)(e->aj > 0.0 ? 0 : 1);
            /* alfa[i] = min alfa[j] over all j included in set P */
            if (i == NULL || fabs(i->aj) > fabs(e->aj)) i = e;
         }
         else if (fabs(e->aj) >= 1e-3)
         {  /* alfa[k] = max alfa[j] over all j not included in set P;
               we skip coefficient a[j] if it is close to zero to avoid
               numerically unreliable results */
            if (k == NULL || fabs(k->aj) < fabs(e->aj)) k = e;
         }
      }
      /* if alfa[k] satisfies to condition (13) for all j in P, include
         x[k] in P */
      if (i != NULL && k != NULL && fabs(i->aj) + fabs(k->aj) > b + eps)
      {  var[++len] = k->xj;
         set[len] = (char)(k->aj > 0.0 ? 0 : 1);
      }
      /* trivial packing inequality being redundant must never appear,
         so we just ignore it */
      if (len < 2) len = 0;
done: drop_form(npp, ptr);
      return len;
}

/***********************************************************************
*  NAME
*
*  npp_is_covering - test if constraint is covering inequality
*
*  SYNOPSIS
*
*  #include "glpnpp.h"
*  int npp_is_covering(NPP *npp, NPPROW *row);
*
*  RETURNS
*
*  If the specified row (constraint) is covering inequality (see below),
*  the routine npp_is_covering returns non-zero. Otherwise, it returns
*  zero.
*
*  COVERING INEQUALITIES
*
*  In canonical format the covering inequality is the following:
*
*     sum  x[j] >= 1,                                                (1)
*    j in J
*
*  where all variables x[j] are binary. This inequality expresses the
*  condition that in any integer feasible solution variables in set J
*  cannot be all equal to zero at the same time, i.e. at least one
*  variable must take non-zero (unity) value. W.l.o.g. it is assumed
*  that |J| >= 2, because if J is empty, the inequality (1) is
*  infeasible, and if |J| = 1, the inequality (1) is a forcing row.
*
*  In general case the covering inequality may include original
*  variables x[j] as well as their complements x~[j]:
*
*     sum   x[j] + sum   x~[j] >= 1,                                 (2)
*    j in Jp      j in Jn
*
*  where Jp and Jn are not intersected. Therefore, using substitution
*  x~[j] = 1 - x[j] gives the packing inequality in generalized format:
*
*     sum   x[j] - sum   x[j] >= 1 - |Jn|.                           (3)
*    j in Jp      j in Jn
*
*  (May note that the inequality (3) cuts off infeasible solutions,
*  where x[j] = 0 for all j in Jp and x[j] = 1 for all j in Jn.)
*
*  NOTE: If |J| = 2, the inequality (3) is equivalent to packing
*        inequality (see the routine npp_is_packing). */

int npp_is_covering(NPP *npp, NPPROW *row)
{     /* test if constraint is covering inequality */
      NPPCOL *col;
      NPPAIJ *aij;
      int b;
      xassert(npp == npp);
      if (!(row->lb != -DBL_MAX && row->ub == +DBL_MAX))
         return 0;
      b = 1;
      for (aij = row->ptr; aij != NULL; aij = aij->r_next)
      {  col = aij->col;
         if (!(col->is_int && col->lb == 0.0 && col->ub == 1.0))
            return 0;
         if (aij->val == +1.0)
            ;
         else if (aij->val == -1.0)
            b--;
         else
            return 0;
      }
      if (row->lb != (double)b) return 0;
      return 1;
}

/***********************************************************************
*  NAME
*
*  npp_hidden_covering - identify hidden covering inequality
*
*  SYNOPSIS
*
*  #include "glpnpp.h"
*  int npp_hidden_covering(NPP *npp, NPPROW *row);
*
*  DESCRIPTION
*
*  The routine npp_hidden_covering processes specified inequality
*  constraint, which includes only binary variables, and the number of
*  the variables is not less than three. If the original inequality is
*  equivalent to a covering inequality (see below), the routine
*  replaces it by the equivalent inequality. If the original constraint
*  is double-sided inequality, it is replaced by a pair of single-sided
*  inequalities, if necessary.
*
*  RETURNS
*
*  If the original inequality constraint was replaced by equivalent
*  covering inequality, the routine npp_hidden_covering returns
*  non-zero. Otherwise, it returns zero.
*
*  PROBLEM TRANSFORMATION
*
*  Consider an inequality constraint:
*
*     sum  a[j] x[j] >= b,                                           (1)
*    j in J
*
*  where all variables x[j] are binary, and |J| >= 3. (In case of '<='
*  inequality it can be transformed to '>=' format by multiplying both
*  its sides by -1.)
*
*  Let Jp = {j: a[j] > 0}, Jn = {j: a[j] < 0}. Performing substitution
*  x[j] = 1 - x~[j] for all j in Jn, we have:
*
*     sum   a[j] x[j] >= b  ==>
*    j in J
*
*     sum   a[j] x[j] + sum   a[j] x[j] >= b  ==>
*    j in Jp           j in Jn
*
*     sum   a[j] x[j] + sum   a[j] (1 - x~[j]) >= b  ==>
*    j in Jp           j in Jn
*
*     sum  m   a[j] x[j] - sum   a[j] x~[j] >= b - sum   a[j].
*    j in Jp              j in Jn                 j in Jn
*
*  Thus, meaning the transformation above, we can assume that in
*  inequality (1) all coefficients a[j] are positive. Moreover, we can
*  assume that b > 0, because otherwise the inequality (1) would be
*  redundant (see the routine npp_analyze_row). It is then obvious that
*  constraint (1) is equivalent to covering inequality only if:
*
*     a[j] >= b,                                                     (2)
*
*  for all j in J.
*
*  Once the original inequality (1) is replaced by equivalent covering
*  inequality, we need to perform back substitution x~[j] = 1 - x[j] for
*  all j in Jn (see above).
*
*  RECOVERING SOLUTION
*
*  None needed. */

static int hidden_covering(NPP *npp, struct elem *ptr, double *_b)
{     /* process inequality constraint: sum a[j] x[j] >= b;
         0 - specified row is NOT hidden covering inequality;
         1 - specified row is covering inequality;
         2 - specified row is hidden covering inequality. */
      struct elem *e;
      int neg;
      double b = *_b, eps;
      xassert(npp == npp);
      /* a[j] must be non-zero, x[j] must be binary, for all j in J */
      for (e = ptr; e != NULL; e = e->next)
      {  xassert(e->aj != 0.0);
         xassert(e->xj->is_int);
         xassert(e->xj->lb == 0.0 && e->xj->ub == 1.0);
      }
      /* check if the specified inequality constraint already has the
         form of covering inequality */
      neg = 0; /* neg is |Jn| */
      for (e = ptr; e != NULL; e = e->next)
      {  if (e->aj == +1.0)
            ;
         else if (e->aj == -1.0)
            neg++;
         else
            break;
      }
      if (e == NULL)
      {  /* all coefficients a[j] are +1 or -1; check rhs b */
         if (b == (double)(1 - neg))
         {  /* it is covering inequality; no processing is needed */
            return 1;
         }
      }
      /* substitute x[j] = 1 - x~[j] for all j in Jn to make all a[j]
         positive; the result is a~[j] = |a[j]| and new rhs b */
      for (e = ptr; e != NULL; e = e->next)
         if (e->aj < 0) b -= e->aj;
      /* now a[j] > 0 for all j in J (actually |a[j]| are used) */
      /* if b <= 0, skip processing--this case must not appear */
      if (b < 1e-3) return 0;
      /* now a[j] > 0 for all j in J, and b > 0 */
      /* the specified constraint is equivalent to covering inequality
         iff a[j] >= b for all j in J */
      eps = 1e-9 + 1e-12 * fabs(b);
      for (e = ptr; e != NULL; e = e->next)
         if (fabs(e->aj) < b - eps) return 0;
      /* perform back substitution x~[j] = 1 - x[j] and construct the
         final equivalent covering inequality in generalized format */
      b = 1.0;
      for (e = ptr; e != NULL; e = e->next)
      {  if (e->aj > 0.0)
            e->aj = +1.0;
         else /* e->aj < 0.0 */
            e->aj = -1.0, b -= 1.0;
      }
      *_b = b;
      return 2;
}

int npp_hidden_covering(NPP *npp, NPPROW *row)
{     /* identify hidden covering inequality */
      NPPROW *copy;
      NPPAIJ *aij;
      struct elem *ptr, *e;
      int kase, ret, count = 0;
      double b;
      /* the row must be inequality constraint */
      xassert(row->lb < row->ub);
      for (kase = 0; kase <= 1; kase++)
      {  if (kase == 0)
         {  /* process row lower bound */
            if (row->lb == -DBL_MAX) continue;
            ptr = copy_form(npp, row, +1.0);
            b = + row->lb;
         }
         else
         {  /* process row upper bound */
            if (row->ub == +DBL_MAX) continue;
            ptr = copy_form(npp, row, -1.0);
            b = - row->ub;
         }
         /* now the inequality has the form "sum a[j] x[j] >= b" */
         ret = hidden_covering(npp, ptr, &b);
         xassert(0 <= ret && ret <= 2);
         if (kase == 1 && ret == 1 || ret == 2)
         {  /* the original inequality has been identified as hidden
               covering inequality */
            count++;
#ifdef GLP_DEBUG
            xprintf("Original constraint:\n");
            for (aij = row->ptr; aij != NULL; aij = aij->r_next)
               xprintf(" %+g x%d", aij->val, aij->col->j);
            if (row->lb != -DBL_MAX) xprintf(", >= %g", row->lb);
            if (row->ub != +DBL_MAX) xprintf(", <= %g", row->ub);
            xprintf("\n");
            xprintf("Equivalent covering inequality:\n");
            for (e = ptr; e != NULL; e = e->next)
               xprintf(" %sx%d", e->aj > 0.0 ? "+" : "-", e->xj->j);
            xprintf(", >= %g\n", b);
#endif
            if (row->lb == -DBL_MAX || row->ub == +DBL_MAX)
            {  /* the original row is single-sided inequality; no copy
                  is needed */
               copy = NULL;
            }
            else
            {  /* the original row is double-sided inequality; we need
                  to create its copy for other bound before replacing it
                  with the equivalent inequality */
               copy = npp_add_row(npp);
               if (kase == 0)
               {  /* the copy is for upper bound */
                  copy->lb = -DBL_MAX, copy->ub = row->ub;
               }
               else
               {  /* the copy is for lower bound */
                  copy->lb = row->lb, copy->ub = +DBL_MAX;
               }
               /* copy original row coefficients */
               for (aij = row->ptr; aij != NULL; aij = aij->r_next)
                  npp_add_aij(npp, copy, aij->col, aij->val);
            }
            /* replace the original inequality by equivalent one */
            npp_erase_row(npp, row);
            row->lb = b, row->ub = +DBL_MAX;
            for (e = ptr; e != NULL; e = e->next)
               npp_add_aij(npp, row, e->xj, e->aj);
            /* continue processing upper bound for the copy */
            if (copy != NULL) row = copy;
         }
         drop_form(npp, ptr);
      }
      return count;
}

/***********************************************************************
*  NAME
*
*  npp_is_partitioning - test if constraint is partitioning equality
*
*  SYNOPSIS
*
*  #include "glpnpp.h"
*  int npp_is_partitioning(NPP *npp, NPPROW *row);
*
*  RETURNS
*
*  If the specified row (constraint) is partitioning equality (see
*  below), the routine npp_is_partitioning returns non-zero. Otherwise,
*  it returns zero.
*
*  PARTITIONING EQUALITIES
*
*  In canonical format the partitioning equality is the following:
*
*     sum  x[j] = 1,                                                 (1)
*    j in J
*
*  where all variables x[j] are binary. This equality expresses the
*  condition that in any integer feasible solution exactly one variable
*  in set J must take non-zero (unity) value while other variables must
*  be equal to zero. W.l.o.g. it is assumed that |J| >= 2, because if
*  J is empty, the inequality (1) is infeasible, and if |J| = 1, the
*  inequality (1) is a fixing row.
*
*  In general case the partitioning equality may include original
*  variables x[j] as well as their complements x~[j]:
*
*     sum   x[j] + sum   x~[j] = 1,                                  (2)
*    j in Jp      j in Jn
*
*  where Jp and Jn are not intersected. Therefore, using substitution
*  x~[j] = 1 - x[j] leads to the partitioning equality in generalized
*  format:
*
*     sum   x[j] - sum   x[j] = 1 - |Jn|.                            (3)
*    j in Jp      j in Jn */

int npp_is_partitioning(NPP *npp, NPPROW *row)
{     /* test if constraint is partitioning equality */
      NPPCOL *col;
      NPPAIJ *aij;
      int b;
      xassert(npp == npp);
      if (row->lb != row->ub) return 0;
      b = 1;
      for (aij = row->ptr; aij != NULL; aij = aij->r_next)
      {  col = aij->col;
         if (!(col->is_int && col->lb == 0.0 && col->ub == 1.0))
            return 0;
         if (aij->val == +1.0)
            ;
         else if (aij->val == -1.0)
            b--;
         else
            return 0;
      }
      if (row->lb != (double)b) return 0;
      return 1;
}

/***********************************************************************
*  NAME
*
*  npp_reduce_ineq_coef - reduce inequality constraint coefficients
*
*  SYNOPSIS
*
*  #include "glpnpp.h"
*  int npp_reduce_ineq_coef(NPP *npp, NPPROW *row);
*
*  DESCRIPTION
*
*  The routine npp_reduce_ineq_coef processes specified inequality
*  constraint attempting to replace it by an equivalent constraint,
*  where magnitude of coefficients at binary variables is smaller than
*  in the original constraint. If the inequality is double-sided, it is
*  replaced by a pair of single-sided inequalities, if necessary.
*
*  RETURNS
*
*  The routine npp_reduce_ineq_coef returns the number of coefficients
*  reduced.
*
*  BACKGROUND
*
*  Consider an inequality constraint:
*
*     sum  a[j] x[j] >= b.                                           (1)
*    j in J
*
*  (In case of '<=' inequality it can be transformed to '>=' format by
*  multiplying both its sides by -1.) Let x[k] be a binary variable;
*  other variables can be integer as well as continuous. We can write
*  constraint (1) as follows:
*
*     a[k] x[k] + t[k] >= b,                                         (2)
*
*  where:
*
*     t[k] = sum      a[j] x[j].                                     (3)
*           j in J\{k}
*
*  Since x[k] is binary, constraint (2) is equivalent to disjunction of
*  the following two constraints:
*
*     x[k] = 0,  t[k] >= b                                           (4)
*
*        OR
*
*     x[k] = 1,  t[k] >= b - a[k].                                   (5)
*
*  Let also that for the partial sum t[k] be known some its implied
*  lower bound inf t[k].
*
*  Case a[k] > 0. Let inf t[k] < b, since otherwise both constraints
*  (4) and (5) and therefore constraint (2) are redundant.
*  If inf t[k] > b - a[k], only constraint (5) is redundant, in which
*  case it can be replaced with the following redundant and therefore
*  equivalent constraint:
*
*     t[k] >= b - a'[k] = inf t[k],                                  (6)
*
*  where:
*
*     a'[k] = b - inf t[k].                                          (7)
*
*  Thus, the original constraint (2) is equivalent to the following
*  constraint with coefficient at variable x[k] changed:
*
*     a'[k] x[k] + t[k] >= b.                                        (8)
*
*  From inf t[k] < b it follows that a'[k] > 0, i.e. the coefficient
*  at x[k] keeps its sign. And from inf t[k] > b - a[k] it follows that
*  a'[k] < a[k], i.e. the coefficient reduces in magnitude.
*
*  Case a[k] < 0. Let inf t[k] < b - a[k], since otherwise both
*  constraints (4) and (5) and therefore constraint (2) are redundant.
*  If inf t[k] > b, only constraint (4) is redundant, in which case it
*  can be replaced with the following redundant and therefore equivalent
*  constraint:
*
*     t[k] >= b' = inf t[k].                                         (9)
*
*  Rewriting constraint (5) as follows:
*
*     t[k] >= b - a[k] = b' - a'[k],                                (10)
*
*  where:
*
*     a'[k] = a[k] + b' - b = a[k] + inf t[k] - b,                  (11)
*
*  we can see that disjunction of constraint (9) and (10) is equivalent
*  to disjunction of constraint (4) and (5), from which it follows that
*  the original constraint (2) is equivalent to the following constraint
*  with both coefficient at variable x[k] and right-hand side changed:
*
*     a'[k] x[k] + t[k] >= b'.                                      (12)
*
*  From inf t[k] < b - a[k] it follows that a'[k] < 0, i.e. the
*  coefficient at x[k] keeps its sign. And from inf t[k] > b it follows
*  that a'[k] > a[k], i.e. the coefficient reduces in magnitude.
*
*  PROBLEM TRANSFORMATION
*
*  In the routine npp_reduce_ineq_coef the following implied lower
*  bound of the partial sum (3) is used:
*
*     inf t[k] = sum       a[j] l[j] + sum       a[j] u[j],         (13)
*               j in Jp\{k}           k in Jn\{k}
*
*  where Jp = {j : a[j] > 0}, Jn = {j : a[j] < 0}, l[j] and u[j] are
*  lower and upper bounds, resp., of variable x[j].
*
*  In order to compute inf t[k] more efficiently, the following formula,
*  which is equivalent to (13), is actually used:
*
*                ( h - a[k] l[k] = h,        if a[k] > 0,
*     inf t[k] = <                                                  (14)
*                ( h - a[k] u[k] = h - a[k], if a[k] < 0,
*
*  where:
*
*     h = sum   a[j] l[j] + sum   a[j] u[j]                         (15)
*        j in Jp           j in Jn
*
*  is the implied lower bound of row (1).
*
*  Reduction of positive coefficient (a[k] > 0) does not change value
*  of h, since l[k] = 0. In case of reduction of negative coefficient
*  (a[k] < 0) from (11) it follows that:
*
*     delta a[k] = a'[k] - a[k] = inf t[k] - b  (> 0),              (16)
*
*  so new value of h (accounting that u[k] = 1) can be computed as
*  follows:
*
*     h := h + delta a[k] = h + (inf t[k] - b).                     (17)
*
*  RECOVERING SOLUTION
*
*  None needed. */

static int reduce_ineq_coef(NPP *npp, struct elem *ptr, double *_b)
{     /* process inequality constraint: sum a[j] x[j] >= b */
      /* returns: the number of coefficients reduced */
      struct elem *e;
      int count = 0;
      double h, inf_t, new_a, b = *_b;
      xassert(npp == npp);
      /* compute h; see (15) */
      h = 0.0;
      for (e = ptr; e != NULL; e = e->next)
      {  if (e->aj > 0.0)
         {  if (e->xj->lb == -DBL_MAX) goto done;
            h += e->aj * e->xj->lb;
         }
         else /* e->aj < 0.0 */
         {  if (e->xj->ub == +DBL_MAX) goto done;
            h += e->aj * e->xj->ub;
         }
      }
      /* perform reduction of coefficients at binary variables */
      for (e = ptr; e != NULL; e = e->next)
      {  /* skip non-binary variable */
         if (!(e->xj->is_int && e->xj->lb == 0.0 && e->xj->ub == 1.0))
            continue;
         if (e->aj > 0.0)
         {  /* compute inf t[k]; see (14) */
            inf_t = h;
            if (b - e->aj < inf_t && inf_t < b)
            {  /* compute reduced coefficient a'[k]; see (7) */
               new_a = b - inf_t;
               if (new_a >= +1e-3 &&
                   e->aj - new_a >= 0.01 * (1.0 + e->aj))
               {  /* accept a'[k] */
#ifdef GLP_DEBUG
                  xprintf("+");
#endif
                  e->aj = new_a;
                  count++;
               }
            }
         }
         else /* e->aj < 0.0 */
         {  /* compute inf t[k]; see (14) */
            inf_t = h - e->aj;
            if (b < inf_t && inf_t < b - e->aj)
            {  /* compute reduced coefficient a'[k]; see (11) */
               new_a = e->aj + (inf_t - b);
               if (new_a <= -1e-3 &&
                   new_a - e->aj >= 0.01 * (1.0 - e->aj))
               {  /* accept a'[k] */
#ifdef GLP_DEBUG
                  xprintf("-");
#endif
                  e->aj = new_a;
                  /* update h; see (17) */
                  h += (inf_t - b);
                  /* compute b'; see (9) */
                  b = inf_t;
                  count++;
               }
            }
         }
      }
      *_b = b;
done: return count;
}

int npp_reduce_ineq_coef(NPP *npp, NPPROW *row)
{     /* reduce inequality constraint coefficients */
      NPPROW *copy;
      NPPAIJ *aij;
      struct elem *ptr, *e;
      int kase, count[2];
      double b;
      /* the row must be inequality constraint */
      xassert(row->lb < row->ub);
      count[0] = count[1] = 0;
      for (kase = 0; kase <= 1; kase++)
      {  if (kase == 0)
         {  /* process row lower bound */
            if (row->lb == -DBL_MAX) continue;
#ifdef GLP_DEBUG
            xprintf("L");
#endif
            ptr = copy_form(npp, row, +1.0);
            b = + row->lb;
         }
         else
         {  /* process row upper bound */
            if (row->ub == +DBL_MAX) continue;
#ifdef GLP_DEBUG
            xprintf("U");
#endif
            ptr = copy_form(npp, row, -1.0);
            b = - row->ub;
         }
         /* now the inequality has the form "sum a[j] x[j] >= b" */
         count[kase] = reduce_ineq_coef(npp, ptr, &b);
         if (count[kase] > 0)
         {  /* the original inequality has been replaced by equivalent
               one with coefficients reduced */
            if (row->lb == -DBL_MAX || row->ub == +DBL_MAX)
            {  /* the original row is single-sided inequality; no copy
                  is needed */
               copy = NULL;
            }
            else
            {  /* the original row is double-sided inequality; we need
                  to create its copy for other bound before replacing it
                  with the equivalent inequality */
#ifdef GLP_DEBUG
               xprintf("*");
#endif
               copy = npp_add_row(npp);
               if (kase == 0)
               {  /* the copy is for upper bound */
                  copy->lb = -DBL_MAX, copy->ub = row->ub;
               }
               else
               {  /* the copy is for lower bound */
                  copy->lb = row->lb, copy->ub = +DBL_MAX;
               }
               /* copy original row coefficients */
               for (aij = row->ptr; aij != NULL; aij = aij->r_next)
                  npp_add_aij(npp, copy, aij->col, aij->val);
            }
            /* replace the original inequality by equivalent one */
            npp_erase_row(npp, row);
            row->lb = b, row->ub = +DBL_MAX;
            for (e = ptr; e != NULL; e = e->next)
               npp_add_aij(npp, row, e->xj, e->aj);
            /* continue processing upper bound for the copy */
            if (copy != NULL) row = copy;
         }
         drop_form(npp, ptr);
      }
      return count[0] + count[1];
}

/* eof */
