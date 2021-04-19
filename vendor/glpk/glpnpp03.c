/* glpnpp03.c */

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
#pragma clang diagnostic ignored "-Wlogical-op-parentheses"
#endif

#include "glpnpp.h"

/***********************************************************************
*  NAME
*
*  npp_empty_row - process empty row
*
*  SYNOPSIS
*
*  #include "glpnpp.h"
*  int npp_empty_row(NPP *npp, NPPROW *p);
*
*  DESCRIPTION
*
*  The routine npp_empty_row processes row p, which is empty, i.e.
*  coefficients at all columns in this row are zero:
*
*     L[p] <= sum 0 x[j] <= U[p],                                    (1)
*
*  where L[p] <= U[p].
*
*  RETURNS
*
*  0 - success;
*
*  1 - problem has no primal feasible solution.
*
*  PROBLEM TRANSFORMATION
*
*  If the following conditions hold:
*
*     L[p] <= +eps,  U[p] >= -eps,                                   (2)
*
*  where eps is an absolute tolerance for row value, the row p is
*  redundant. In this case it can be replaced by equivalent redundant
*  row, which is free (unbounded), and then removed from the problem.
*  Otherwise, the row p is infeasible and, thus, the problem has no
*  primal feasible solution.
*
*  RECOVERING BASIC SOLUTION
*
*  See the routine npp_free_row.
*
*  RECOVERING INTERIOR-POINT SOLUTION
*
*  See the routine npp_free_row.
*
*  RECOVERING MIP SOLUTION
*
*  None needed. */

int npp_empty_row(NPP *npp, NPPROW *p)
{     /* process empty row */
      double eps = 1e-3;
      /* the row must be empty */
      xassert(p->ptr == NULL);
      /* check primal feasibility */
      if (p->lb > +eps || p->ub < -eps)
         return 1;
      /* replace the row by equivalent free (unbounded) row */
      p->lb = -DBL_MAX, p->ub = +DBL_MAX;
      /* and process it */
      npp_free_row(npp, p);
      return 0;
}

/***********************************************************************
*  NAME
*
*  npp_empty_col - process empty column
*
*  SYNOPSIS
*
*  #include "glpnpp.h"
*  int npp_empty_col(NPP *npp, NPPCOL *q);
*
*  DESCRIPTION
*
*  The routine npp_empty_col processes column q:
*
*     l[q] <= x[q] <= u[q],                                          (1)
*
*  where l[q] <= u[q], which is empty, i.e. has zero coefficients in
*  all constraint rows.
*
*  RETURNS
*
*  0 - success;
*
*  1 - problem has no dual feasible solution.
*
*  PROBLEM TRANSFORMATION
*
*  The row of the dual system corresponding to the empty column is the
*  following:
*
*     sum 0 pi[i] + lambda[q] = c[q],                                (2)
*      i
*
*  from which it follows that:
*
*     lambda[q] = c[q].                                              (3)
*
*  If the following condition holds:
*
*     c[q] < - eps,                                                  (4)
*
*  where eps is an absolute tolerance for column multiplier, the lower
*  column bound l[q] must be active to provide dual feasibility (note
*  that being preprocessed the problem is always minimization). In this
*  case the column can be fixed on its lower bound and removed from the
*  problem (if the column is integral, its bounds are also assumed to
*  be integral). And if the column has no lower bound (l[q] = -oo), the
*  problem has no dual feasible solution.
*
*  If the following condition holds:
*
*     c[q] > + eps,                                                  (5)
*
*  the upper column bound u[q] must be active to provide dual
*  feasibility. In this case the column can be fixed on its upper bound
*  and removed from the problem. And if the column has no upper bound
*  (u[q] = +oo), the problem has no dual feasible solution.
*
*  Finally, if the following condition holds:
*
*     - eps <= c[q] <= +eps,                                         (6)
*
*  dual feasibility does not depend on a particular value of column q.
*  In this case the column can be fixed either on its lower bound (if
*  l[q] > -oo) or on its upper bound (if u[q] < +oo) or at zero (if the
*  column is unbounded) and then removed from the problem.
*
*  RECOVERING BASIC SOLUTION
*
*  See the routine npp_fixed_col. Having been recovered the column
*  is assigned status GLP_NS. However, if actually it is not fixed
*  (l[q] < u[q]), its status should be changed to GLP_NL, GLP_NU, or
*  GLP_NF depending on which bound it was fixed on transformation stage.
*
*  RECOVERING INTERIOR-POINT SOLUTION
*
*  See the routine npp_fixed_col.
*
*  RECOVERING MIP SOLUTION
*
*  See the routine npp_fixed_col. */

struct empty_col
{     /* empty column */
      int q;
      /* column reference number */
      char stat;
      /* status in basic solution */
};

static int rcv_empty_col(NPP *npp, void *info);

int npp_empty_col(NPP *npp, NPPCOL *q)
{     /* process empty column */
      struct empty_col *info;
      double eps = 1e-3;
      /* the column must be empty */
      xassert(q->ptr == NULL);
      /* check dual feasibility */
      if (q->coef > +eps && q->lb == -DBL_MAX)
         return 1;
      if (q->coef < -eps && q->ub == +DBL_MAX)
         return 1;
      /* create transformation stack entry */
      info = npp_push_tse(npp,
         rcv_empty_col, sizeof(struct empty_col));
      info->q = q->j;
      /* fix the column */
      if (q->lb == -DBL_MAX && q->ub == +DBL_MAX)
      {  /* free column */
         info->stat = GLP_NF;
         q->lb = q->ub = 0.0;
      }
      else if (q->ub == +DBL_MAX)
lo:   {  /* column with lower bound */
         info->stat = GLP_NL;
         q->ub = q->lb;
      }
      else if (q->lb == -DBL_MAX)
up:   {  /* column with upper bound */
         info->stat = GLP_NU;
         q->lb = q->ub;
      }
      else if (q->lb != q->ub)
      {  /* double-bounded column */
         if (q->coef >= +DBL_EPSILON) goto lo;
         if (q->coef <= -DBL_EPSILON) goto up;
         if (fabs(q->lb) <= fabs(q->ub)) goto lo; else goto up;
      }
      else
      {  /* fixed column */
         info->stat = GLP_NS;
      }
      /* process fixed column */
      npp_fixed_col(npp, q);
      return 0;
}

static int rcv_empty_col(NPP *npp, void *_info)
{     /* recover empty column */
      struct empty_col *info = _info;
      if (npp->sol == GLP_SOL)
         npp->c_stat[info->q] = info->stat;
      return 0;
}

/***********************************************************************
*  NAME
*
*  npp_implied_value - process implied column value
*
*  SYNOPSIS
*
*  #include "glpnpp.h"
*  int npp_implied_value(NPP *npp, NPPCOL *q, double s);
*
*  DESCRIPTION
*
*  For column q:
*
*     l[q] <= x[q] <= u[q],                                          (1)
*
*  where l[q] < u[q], the routine npp_implied_value processes its
*  implied value s[q]. If this implied value satisfies to the current
*  column bounds and integrality condition, the routine fixes column q
*  at the given point. Note that the column is kept in the problem in
*  any case.
*
*  RETURNS
*
*  0 - column has been fixed;
*
*  1 - implied value violates to current column bounds;
*
*  2 - implied value violates integrality condition.
*
*  ALGORITHM
*
*  Implied column value s[q] satisfies to the current column bounds if
*  the following condition holds:
*
*     l[q] - eps <= s[q] <= u[q] + eps,                              (2)
*
*  where eps is an absolute tolerance for column value. If the column
*  is integral, the following condition also must hold:
*
*     |s[q] - floor(s[q]+0.5)| <= eps,                               (3)
*
*  where floor(s[q]+0.5) is the nearest integer to s[q].
*
*  If both condition (2) and (3) are satisfied, the column can be fixed
*  at the value s[q], or, if it is integral, at floor(s[q]+0.5).
*  Otherwise, if s[q] violates (2) or (3), the problem has no feasible
*  solution.
*
*  Note: If s[q] is close to l[q] or u[q], it seems to be reasonable to
*  fix the column at its lower or upper bound, resp. rather than at the
*  implied value. */

int npp_implied_value(NPP *npp, NPPCOL *q, double s)
{     /* process implied column value */
      double eps, nint;
      xassert(npp == npp);
      /* column must not be fixed */
      xassert(q->lb < q->ub);
      /* check integrality */
      if (q->is_int)
      {  nint = floor(s + 0.5);
         if (fabs(s - nint) <= 1e-5)
            s = nint;
         else
            return 2;
      }
      /* check current column lower bound */
      if (q->lb != -DBL_MAX)
      {  eps = (q->is_int ? 1e-5 : 1e-5 + 1e-8 * fabs(q->lb));
         if (s < q->lb - eps) return 1;
         /* if s[q] is close to l[q], fix column at its lower bound
            rather than at the implied value */
         if (s < q->lb + 1e-3 * eps)
         {  q->ub = q->lb;
            return 0;
         }
      }
      /* check current column upper bound */
      if (q->ub != +DBL_MAX)
      {  eps = (q->is_int ? 1e-5 : 1e-5 + 1e-8 * fabs(q->ub));
         if (s > q->ub + eps) return 1;
         /* if s[q] is close to u[q], fix column at its upper bound
            rather than at the implied value */
         if (s > q->ub - 1e-3 * eps)
         {  q->lb = q->ub;
            return 0;
         }
      }
      /* fix column at the implied value */
      q->lb = q->ub = s;
      return 0;
}

/***********************************************************************
*  NAME
*
*  npp_eq_singlet - process row singleton (equality constraint)
*
*  SYNOPSIS
*
*  #include "glpnpp.h"
*  int npp_eq_singlet(NPP *npp, NPPROW *p);
*
*  DESCRIPTION
*
*  The routine npp_eq_singlet processes row p, which is equiality
*  constraint having the only non-zero coefficient:
*
*     a[p,q] x[q] = b.                                               (1)
*
*  RETURNS
*
*  0 - success;
*
*  1 - problem has no primal feasible solution;
*
*  2 - problem has no integer feasible solution.
*
*  PROBLEM TRANSFORMATION
*
*  The equality constraint defines implied value of column q:
*
*     x[q] = s[q] = b / a[p,q].                                      (2)
*
*  If the implied value s[q] satisfies to the column bounds (see the
*  routine npp_implied_value), the column can be fixed at s[q] and
*  removed from the problem. In this case row p becomes redundant, so
*  it can be replaced by equivalent free row and also removed from the
*  problem.
*
*  Note that the routine removes from the problem only row p. Column q
*  becomes fixed, however, it is kept in the problem.
*
*  RECOVERING BASIC SOLUTION
*
*  In solution to the original problem row p is assigned status GLP_NS
*  (active equality constraint), and column q is assigned status GLP_BS
*  (basic column).
*
*  Multiplier for row p can be computed as follows. In the dual system
*  of the original problem column q corresponds to the following row:
*
*     sum a[i,q] pi[i] + lambda[q] = c[q]  ==>
*      i
*
*     sum a[i,q] pi[i] + a[p,q] pi[p] + lambda[q] = c[q].
*     i!=p
*
*  Therefore:
*
*               1
*     pi[p] = ------ (c[q] - lambda[q] - sum a[i,q] pi[i]),          (3)
*             a[p,q]                     i!=q
*
*  where lambda[q] = 0 (since column[q] is basic), and pi[i] for all
*  i != p are known in solution to the transformed problem.
*
*  Value of column q in solution to the original problem is assigned
*  its implied value s[q].
*
*  RECOVERING INTERIOR-POINT SOLUTION
*
*  Multiplier for row p is computed with formula (3). Value of column
*  q is assigned its implied value s[q].
*
*  RECOVERING MIP SOLUTION
*
*  Value of column q is assigned its implied value s[q]. */

struct eq_singlet
{     /* row singleton (equality constraint) */
      int p;
      /* row reference number */
      int q;
      /* column reference number */
      double apq;
      /* constraint coefficient a[p,q] */
      double c;
      /* objective coefficient at x[q] */
      NPPLFE *ptr;
      /* list of non-zero coefficients a[i,q], i != p */
};

static int rcv_eq_singlet(NPP *npp, void *info);

int npp_eq_singlet(NPP *npp, NPPROW *p)
{     /* process row singleton (equality constraint) */
      struct eq_singlet *info;
      NPPCOL *q;
      NPPAIJ *aij;
      NPPLFE *lfe;
      int ret;
      double s;
      /* the row must be singleton equality constraint */
      xassert(p->lb == p->ub);
      xassert(p->ptr != NULL && p->ptr->r_next == NULL);
      /* compute and process implied column value */
      aij = p->ptr;
      q = aij->col;
      s = p->lb / aij->val;
      ret = npp_implied_value(npp, q, s);
      xassert(0 <= ret && ret <= 2);
      if (ret != 0) return ret;
      /* create transformation stack entry */
      info = npp_push_tse(npp,
         rcv_eq_singlet, sizeof(struct eq_singlet));
      info->p = p->i;
      info->q = q->j;
      info->apq = aij->val;
      info->c = q->coef;
      info->ptr = NULL;
      /* save column coefficients a[i,q], i != p (not needed for MIP
         solution) */
      if (npp->sol != GLP_MIP)
      {  for (aij = q->ptr; aij != NULL; aij = aij->c_next)
         {  if (aij->row == p) continue; /* skip a[p,q] */
            lfe = dmp_get_atom(npp->stack, sizeof(NPPLFE));
            lfe->ref = aij->row->i;
            lfe->val = aij->val;
            lfe->next = info->ptr;
            info->ptr = lfe;
         }
      }
      /* remove the row from the problem */
      npp_del_row(npp, p);
      return 0;
}

static int rcv_eq_singlet(NPP *npp, void *_info)
{     /* recover row singleton (equality constraint) */
      struct eq_singlet *info = _info;
      NPPLFE *lfe;
      double temp;
      if (npp->sol == GLP_SOL)
      {  /* column q must be already recovered as GLP_NS */
         if (npp->c_stat[info->q] != GLP_NS)
         {  npp_error();
            return 1;
         }
         npp->r_stat[info->p] = GLP_NS;
         npp->c_stat[info->q] = GLP_BS;
      }
      if (npp->sol != GLP_MIP)
      {  /* compute multiplier for row p with formula (3) */
         temp = info->c;
         for (lfe = info->ptr; lfe != NULL; lfe = lfe->next)
            temp -= lfe->val * npp->r_pi[lfe->ref];
         npp->r_pi[info->p] = temp / info->apq;
      }
      return 0;
}

/***********************************************************************
*  NAME
*
*  npp_implied_lower - process implied column lower bound
*
*  SYNOPSIS
*
*  #include "glpnpp.h"
*  int npp_implied_lower(NPP *npp, NPPCOL *q, double l);
*
*  DESCRIPTION
*
*  For column q:
*
*     l[q] <= x[q] <= u[q],                                          (1)
*
*  where l[q] < u[q], the routine npp_implied_lower processes its
*  implied lower bound l'[q]. As the result the current column lower
*  bound may increase. Note that the column is kept in the problem in
*  any case.
*
*  RETURNS
*
*  0 - current column lower bound has not changed;
*
*  1 - current column lower bound has changed, but not significantly;
*
*  2 - current column lower bound has significantly changed;
*
*  3 - column has been fixed on its upper bound;
*
*  4 - implied lower bound violates current column upper bound.
*
*  ALGORITHM
*
*  If column q is integral, before processing its implied lower bound
*  should be rounded up:
*
*              ( floor(l'[q]+0.5), if |l'[q] - floor(l'[q]+0.5)| <= eps
*     l'[q] := <                                                     (2)
*              ( ceil(l'[q]),      otherwise
*
*  where floor(l'[q]+0.5) is the nearest integer to l'[q], ceil(l'[q])
*  is smallest integer not less than l'[q], and eps is an absolute
*  tolerance for column value.
*
*  Processing implied column lower bound l'[q] includes the following
*  cases:
*
*  1) if l'[q] < l[q] + eps, implied lower bound is redundant;
*
*  2) if l[q] + eps <= l[q] <= u[q] + eps, current column lower bound
*     l[q] can be strengthened by replacing it with l'[q]. If in this
*     case new column lower bound becomes close to current column upper
*     bound u[q], the column can be fixed on its upper bound;
*
*  3) if l'[q] > u[q] + eps, implied lower bound violates current
*     column upper bound u[q], in which case the problem has no primal
*     feasible solution. */

int npp_implied_lower(NPP *npp, NPPCOL *q, double l)
{     /* process implied column lower bound */
      int ret;
      double eps, nint;
      xassert(npp == npp);
      /* column must not be fixed */
      xassert(q->lb < q->ub);
      /* implied lower bound must be finite */
      xassert(l != -DBL_MAX);
      /* if column is integral, round up l'[q] */
      if (q->is_int)
      {  nint = floor(l + 0.5);
         if (fabs(l - nint) <= 1e-5)
            l = nint;
         else
            l = ceil(l);
      }
      /* check current column lower bound */
      if (q->lb != -DBL_MAX)
      {  eps = (q->is_int ? 1e-3 : 1e-3 + 1e-6 * fabs(q->lb));
         if (l < q->lb + eps)
         {  ret = 0; /* redundant */
            goto done;
         }
      }
      /* check current column upper bound */
      if (q->ub != +DBL_MAX)
      {  eps = (q->is_int ? 1e-5 : 1e-5 + 1e-8 * fabs(q->ub));
         if (l > q->ub + eps)
         {  ret = 4; /* infeasible */
            goto done;
         }
         /* if l'[q] is close to u[q], fix column at its upper bound */
         if (l > q->ub - 1e-3 * eps)
         {  q->lb = q->ub;
            ret = 3; /* fixed */
            goto done;
         }
      }
      /* check if column lower bound changes significantly */
      if (q->lb == -DBL_MAX)
         ret = 2; /* significantly */
      else if (q->is_int && l > q->lb + 0.5)
         ret = 2; /* significantly */
      else if (l > q->lb + 0.30 * (1.0 + fabs(q->lb)))
         ret = 2; /* significantly */
      else
         ret = 1; /* not significantly */
      /* set new column lower bound */
      q->lb = l;
done: return ret;
}

/***********************************************************************
*  NAME
*
*  npp_implied_upper - process implied column upper bound
*
*  SYNOPSIS
*
*  #include "glpnpp.h"
*  int npp_implied_upper(NPP *npp, NPPCOL *q, double u);
*
*  DESCRIPTION
*
*  For column q:
*
*     l[q] <= x[q] <= u[q],                                          (1)
*
*  where l[q] < u[q], the routine npp_implied_upper processes its
*  implied upper bound u'[q]. As the result the current column upper
*  bound may decrease. Note that the column is kept in the problem in
*  any case.
*
*  RETURNS
*
*  0 - current column upper bound has not changed;
*
*  1 - current column upper bound has changed, but not significantly;
*
*  2 - current column upper bound has significantly changed;
*
*  3 - column has been fixed on its lower bound;
*
*  4 - implied upper bound violates current column lower bound.
*
*  ALGORITHM
*
*  If column q is integral, before processing its implied upper bound
*  should be rounded down:
*
*              ( floor(u'[q]+0.5), if |u'[q] - floor(l'[q]+0.5)| <= eps
*     u'[q] := <                                                     (2)
*              ( floor(l'[q]),     otherwise
*
*  where floor(u'[q]+0.5) is the nearest integer to u'[q],
*  floor(u'[q]) is largest integer not greater than u'[q], and eps is
*  an absolute tolerance for column value.
*
*  Processing implied column upper bound u'[q] includes the following
*  cases:
*
*  1) if u'[q] > u[q] - eps, implied upper bound is redundant;
*
*  2) if l[q] - eps <= u[q] <= u[q] - eps, current column upper bound
*     u[q] can be strengthened by replacing it with u'[q]. If in this
*     case new column upper bound becomes close to current column lower
*     bound, the column can be fixed on its lower bound;
*
*  3) if u'[q] < l[q] - eps, implied upper bound violates current
*     column lower bound l[q], in which case the problem has no primal
*     feasible solution. */

int npp_implied_upper(NPP *npp, NPPCOL *q, double u)
{     int ret;
      double eps, nint;
      xassert(npp == npp);
      /* column must not be fixed */
      xassert(q->lb < q->ub);
      /* implied upper bound must be finite */
      xassert(u != +DBL_MAX);
      /* if column is integral, round down u'[q] */
      if (q->is_int)
      {  nint = floor(u + 0.5);
         if (fabs(u - nint) <= 1e-5)
            u = nint;
         else
            u = floor(u);
      }
      /* check current column upper bound */
      if (q->ub != +DBL_MAX)
      {  eps = (q->is_int ? 1e-3 : 1e-3 + 1e-6 * fabs(q->ub));
         if (u > q->ub - eps)
         {  ret = 0; /* redundant */
            goto done;
         }
      }
      /* check current column lower bound */
      if (q->lb != -DBL_MAX)
      {  eps = (q->is_int ? 1e-5 : 1e-5 + 1e-8 * fabs(q->lb));
         if (u < q->lb - eps)
         {  ret = 4; /* infeasible */
            goto done;
         }
         /* if u'[q] is close to l[q], fix column at its lower bound */
         if (u < q->lb + 1e-3 * eps)
         {  q->ub = q->lb;
            ret = 3; /* fixed */
            goto done;
         }
      }
      /* check if column upper bound changes significantly */
      if (q->ub == +DBL_MAX)
         ret = 2; /* significantly */
      else if (q->is_int && u < q->ub - 0.5)
         ret = 2; /* significantly */
      else if (u < q->ub - 0.30 * (1.0 + fabs(q->ub)))
         ret = 2; /* significantly */
      else
         ret = 1; /* not significantly */
      /* set new column upper bound */
      q->ub = u;
done: return ret;
}

/***********************************************************************
*  NAME
*
*  npp_ineq_singlet - process row singleton (inequality constraint)
*
*  SYNOPSIS
*
*  #include "glpnpp.h"
*  int npp_ineq_singlet(NPP *npp, NPPROW *p);
*
*  DESCRIPTION
*
*  The routine npp_ineq_singlet processes row p, which is inequality
*  constraint having the only non-zero coefficient:
*
*     L[p] <= a[p,q] * x[q] <= U[p],                                 (1)
*
*  where L[p] < U[p], L[p] > -oo and/or U[p] < +oo.
*
*  RETURNS
*
*  0 - current column bounds have not changed;
*
*  1 - current column bounds have changed, but not significantly;
*
*  2 - current column bounds have significantly changed;
*
*  3 - column has been fixed on its lower or upper bound;
*
*  4 - problem has no primal feasible solution.
*
*  PROBLEM TRANSFORMATION
*
*  Inequality constraint (1) defines implied bounds of column q:
*
*             (  L[p] / a[p,q],  if a[p,q] > 0
*     l'[q] = <                                                      (2)
*             (  U[p] / a[p,q],  if a[p,q] < 0
*
*             (  U[p] / a[p,q],  if a[p,q] > 0
*     u'[q] = <                                                      (3)
*             (  L[p] / a[p,q],  if a[p,q] < 0
*
*  If these implied bounds do not violate current bounds of column q:
*
*     l[q] <= x[q] <= u[q],                                          (4)
*
*  they can be used to strengthen the current column bounds:
*
*     l[q] := max(l[q], l'[q]),                                      (5)
*
*     u[q] := min(u[q], u'[q]).                                      (6)
*
*  (See the routines npp_implied_lower and npp_implied_upper.)
*
*  Once bounds of row p (1) have been carried over column q, the row
*  becomes redundant, so it can be replaced by equivalent free row and
*  removed from the problem.
*
*  Note that the routine removes from the problem only row p. Column q,
*  even it has been fixed, is kept in the problem.
*
*  RECOVERING BASIC SOLUTION
*
*  Note that the row in the dual system corresponding to column q is
*  the following:
*
*     sum a[i,q] pi[i] + lambda[q] = c[q]  ==>
*      i
*                                                                    (7)
*     sum a[i,q] pi[i] + a[p,q] pi[p] + lambda[q] = c[q],
*     i!=p
*
*  where pi[i] for all i != p are known in solution to the transformed
*  problem. Row p does not exist in the transformed problem, so it has
*  zero multiplier there. This allows computing multiplier for column q
*  in solution to the transformed problem:
*
*     lambda~[q] = c[q] - sum a[i,q] pi[i].                          (8)
*                         i!=p
*
*  Let in solution to the transformed problem column q be non-basic
*  with lower bound active (GLP_NL, lambda~[q] >= 0), and this lower
*  bound be implied one l'[q]. From the original problem's standpoint
*  this then means that actually the original column lower bound l[q]
*  is inactive, and active is that row bound L[p] or U[p] that defines
*  the implied bound l'[q] (2). In this case in solution to the
*  original problem column q is assigned status GLP_BS while row p is
*  assigned status GLP_NL (if a[p,q] > 0) or GLP_NU (if a[p,q] < 0).
*  Since now column q is basic, its multiplier lambda[q] is zero. This
*  allows using (7) and (8) to find multiplier for row p in solution to
*  the original problem:
*
*               1
*     pi[p] = ------ (c[q] - sum a[i,q] pi[i]) = lambda~[q] / a[p,q] (9)
*             a[p,q]         i!=p
*
*  Now let in solution to the transformed problem column q be non-basic
*  with upper bound active (GLP_NU, lambda~[q] <= 0), and this upper
*  bound be implied one u'[q]. As in the previous case this then means
*  that from the original problem's standpoint actually the original
*  column upper bound u[q] is inactive, and active is that row bound
*  L[p] or U[p] that defines the implied bound u'[q] (3). In this case
*  in solution to the original problem column q is assigned status
*  GLP_BS, row p is assigned status GLP_NU (if a[p,q] > 0) or GLP_NL
*  (if a[p,q] < 0), and its multiplier is computed with formula (9).
*
*  Strengthening bounds of column q according to (5) and (6) may make
*  it fixed. Thus, if in solution to the transformed problem column q is
*  non-basic and fixed (GLP_NS), we can suppose that if lambda~[q] > 0,
*  column q has active lower bound (GLP_NL), and if lambda~[q] < 0,
*  column q has active upper bound (GLP_NU), reducing this case to two
*  previous ones. If, however, lambda~[q] is close to zero or
*  corresponding bound of row p does not exist (this may happen if
*  lambda~[q] has wrong sign due to round-off errors, in which case it
*  is expected to be close to zero, since solution is assumed to be dual
*  feasible), column q can be assigned status GLP_BS (basic), and row p
*  can be made active on its existing bound. In the latter case row
*  multiplier pi[p] computed with formula (9) will be also close to
*  zero, and dual feasibility will be kept.
*
*  In all other cases, namely, if in solution to the transformed
*  problem column q is basic (GLP_BS), or non-basic with original lower
*  bound l[q] active (GLP_NL), or non-basic with original upper bound
*  u[q] active (GLP_NU), constraint (1) is inactive. So in solution to
*  the original problem status of column q remains unchanged, row p is
*  assigned status GLP_BS, and its multiplier pi[p] is assigned zero
*  value.
*
*  RECOVERING INTERIOR-POINT SOLUTION
*
*  First, value of multiplier for column q in solution to the original
*  problem is computed with formula (8). If lambda~[q] > 0 and column q
*  has implied lower bound, or if lambda~[q] < 0 and column q has
*  implied upper bound, this means that from the original problem's
*  standpoint actually row p has corresponding active bound, in which
*  case its multiplier pi[p] is computed with formula (9). In other
*  cases, when the sign of lambda~[q] corresponds to original bound of
*  column q, or when lambda~[q] =~ 0, value of row multiplier pi[p] is
*  assigned zero value.
*
*  RECOVERING MIP SOLUTION
*
*  None needed. */

struct ineq_singlet
{     /* row singleton (inequality constraint) */
      int p;
      /* row reference number */
      int q;
      /* column reference number */
      double apq;
      /* constraint coefficient a[p,q] */
      double c;
      /* objective coefficient at x[q] */
      double lb;
      /* row lower bound */
      double ub;
      /* row upper bound */
      char lb_changed;
      /* this flag is set if column lower bound was changed */
      char ub_changed;
      /* this flag is set if column upper bound was changed */
      NPPLFE *ptr;
      /* list of non-zero coefficients a[i,q], i != p */
};

static int rcv_ineq_singlet(NPP *npp, void *info);

int npp_ineq_singlet(NPP *npp, NPPROW *p)
{     /* process row singleton (inequality constraint) */
      struct ineq_singlet *info;
      NPPCOL *q;
      NPPAIJ *apq, *aij;
      NPPLFE *lfe;
      int lb_changed, ub_changed;
      double ll, uu;
      /* the row must be singleton inequality constraint */
      xassert(p->lb != -DBL_MAX || p->ub != +DBL_MAX);
      xassert(p->lb < p->ub);
      xassert(p->ptr != NULL && p->ptr->r_next == NULL);
      /* compute implied column bounds */
      apq = p->ptr;
      q = apq->col;
      xassert(q->lb < q->ub);
      if (apq->val > 0.0)
      {  ll = (p->lb == -DBL_MAX ? -DBL_MAX : p->lb / apq->val);
         uu = (p->ub == +DBL_MAX ? +DBL_MAX : p->ub / apq->val);
      }
      else
      {  ll = (p->ub == +DBL_MAX ? -DBL_MAX : p->ub / apq->val);
         uu = (p->lb == -DBL_MAX ? +DBL_MAX : p->lb / apq->val);
      }
      /* process implied column lower bound */
      if (ll == -DBL_MAX)
         lb_changed = 0;
      else
      {  lb_changed = npp_implied_lower(npp, q, ll);
         xassert(0 <= lb_changed && lb_changed <= 4);
         if (lb_changed == 4) return 4; /* infeasible */
      }
      /* process implied column upper bound */
      if (uu == +DBL_MAX)
         ub_changed = 0;
      else if (lb_changed == 3)
      {  /* column was fixed on its upper bound due to l'[q] = u[q] */
         /* note that L[p] < U[p], so l'[q] = u[q] < u'[q] */
         ub_changed = 0;
      }
      else
      {  ub_changed = npp_implied_upper(npp, q, uu);
         xassert(0 <= ub_changed && ub_changed <= 4);
         if (ub_changed == 4) return 4; /* infeasible */
      }
      /* if neither lower nor upper column bound was changed, the row
         is originally redundant and can be replaced by free row */
      if (!lb_changed && !ub_changed)
      {  p->lb = -DBL_MAX, p->ub = +DBL_MAX;
         npp_free_row(npp, p);
         return 0;
      }
      /* create transformation stack entry */
      info = npp_push_tse(npp,
         rcv_ineq_singlet, sizeof(struct ineq_singlet));
      info->p = p->i;
      info->q = q->j;
      info->apq = apq->val;
      info->c = q->coef;
      info->lb = p->lb;
      info->ub = p->ub;
      info->lb_changed = (char)lb_changed;
      info->ub_changed = (char)ub_changed;
      info->ptr = NULL;
      /* save column coefficients a[i,q], i != p (not needed for MIP
         solution) */
      if (npp->sol != GLP_MIP)
      {  for (aij = q->ptr; aij != NULL; aij = aij->c_next)
         {  if (aij == apq) continue; /* skip a[p,q] */
            lfe = dmp_get_atom(npp->stack, sizeof(NPPLFE));
            lfe->ref = aij->row->i;
            lfe->val = aij->val;
            lfe->next = info->ptr;
            info->ptr = lfe;
         }
      }
      /* remove the row from the problem */
      npp_del_row(npp, p);
      return lb_changed >= ub_changed ? lb_changed : ub_changed;
}

static int rcv_ineq_singlet(NPP *npp, void *_info)
{     /* recover row singleton (inequality constraint) */
      struct ineq_singlet *info = _info;
      NPPLFE *lfe;
      double lambda;
      if (npp->sol == GLP_MIP) goto done;
      /* compute lambda~[q] in solution to the transformed problem
         with formula (8) */
      lambda = info->c;
      for (lfe = info->ptr; lfe != NULL; lfe = lfe->next)
         lambda -= lfe->val * npp->r_pi[lfe->ref];
      if (npp->sol == GLP_SOL)
      {  /* recover basic solution */
         if (npp->c_stat[info->q] == GLP_BS)
         {  /* column q is basic, so row p is inactive */
            npp->r_stat[info->p] = GLP_BS;
            npp->r_pi[info->p] = 0.0;
         }
         else if (npp->c_stat[info->q] == GLP_NL)
nl:      {  /* column q is non-basic with lower bound active */
            if (info->lb_changed)
            {  /* it is implied bound, so actually row p is active
                  while column q is basic */
               npp->r_stat[info->p] =
                  (char)(info->apq > 0.0 ? GLP_NL : GLP_NU);
               npp->c_stat[info->q] = GLP_BS;
               npp->r_pi[info->p] = lambda / info->apq;
            }
            else
            {  /* it is original bound, so row p is inactive */
               npp->r_stat[info->p] = GLP_BS;
               npp->r_pi[info->p] = 0.0;
            }
         }
         else if (npp->c_stat[info->q] == GLP_NU)
nu:      {  /* column q is non-basic with upper bound active */
            if (info->ub_changed)
            {  /* it is implied bound, so actually row p is active
                  while column q is basic */
               npp->r_stat[info->p] =
                  (char)(info->apq > 0.0 ? GLP_NU : GLP_NL);
               npp->c_stat[info->q] = GLP_BS;
               npp->r_pi[info->p] = lambda / info->apq;
            }
            else
            {  /* it is original bound, so row p is inactive */
               npp->r_stat[info->p] = GLP_BS;
               npp->r_pi[info->p] = 0.0;
            }
         }
         else if (npp->c_stat[info->q] == GLP_NS)
         {  /* column q is non-basic and fixed; note, however, that in
               in the original problem it is non-fixed */
            if (lambda > +1e-7)
            {  if (info->apq > 0.0 && info->lb != -DBL_MAX ||
                   info->apq < 0.0 && info->ub != +DBL_MAX ||
                  !info->lb_changed)
               {  /* either corresponding bound of row p exists or
                     column q remains non-basic with its original lower
                     bound active */
                  npp->c_stat[info->q] = GLP_NL;
                  goto nl;
               }
            }
            if (lambda < -1e-7)
            {  if (info->apq > 0.0 && info->ub != +DBL_MAX ||
                   info->apq < 0.0 && info->lb != -DBL_MAX ||
                  !info->ub_changed)
               {  /* either corresponding bound of row p exists or
                     column q remains non-basic with its original upper
                     bound active */
                  npp->c_stat[info->q] = GLP_NU;
                  goto nu;
               }
            }
            /* either lambda~[q] is close to zero, or corresponding
               bound of row p does not exist, because lambda~[q] has
               wrong sign due to round-off errors; in the latter case
               lambda~[q] is also assumed to be close to zero; so, we
               can make row p active on its existing bound and column q
               basic; pi[p] will have wrong sign, but it also will be
               close to zero (rarus casus of dual degeneracy) */
            if (info->lb != -DBL_MAX && info->ub == +DBL_MAX)
            {  /* row lower bound exists, but upper bound doesn't */
               npp->r_stat[info->p] = GLP_NL;
            }
            else if (info->lb == -DBL_MAX && info->ub != +DBL_MAX)
            {  /* row upper bound exists, but lower bound doesn't */
               npp->r_stat[info->p] = GLP_NU;
            }
            else if (info->lb != -DBL_MAX && info->ub != +DBL_MAX)
            {  /* both row lower and upper bounds exist */
               /* to choose proper active row bound we should not use
                  lambda~[q], because its value being close to zero is
                  unreliable; so we choose that bound which provides
                  primal feasibility for original constraint (1) */
               if (info->apq * npp->c_value[info->q] <=
                   0.5 * (info->lb + info->ub))
                  npp->r_stat[info->p] = GLP_NL;
               else
                  npp->r_stat[info->p] = GLP_NU;
            }
            else
            {  npp_error();
               return 1;
            }
            npp->c_stat[info->q] = GLP_BS;
            npp->r_pi[info->p] = lambda / info->apq;
         }
         else
         {  npp_error();
            return 1;
         }
      }
      if (npp->sol == GLP_IPT)
      {  /* recover interior-point solution */
         if (lambda > +DBL_EPSILON && info->lb_changed ||
             lambda < -DBL_EPSILON && info->ub_changed)
         {  /* actually row p has corresponding active bound */
            npp->r_pi[info->p] = lambda / info->apq;
         }
         else
         {  /* either bounds of column q are both inactive or its
               original bound is active */
            npp->r_pi[info->p] = 0.0;
         }
      }
done: return 0;
}

/***********************************************************************
*  NAME
*
*  npp_implied_slack - process column singleton (implied slack variable)
*
*  SYNOPSIS
*
*  #include "glpnpp.h"
*  void npp_implied_slack(NPP *npp, NPPCOL *q);
*
*  DESCRIPTION
*
*  The routine npp_implied_slack processes column q:
*
*     l[q] <= x[q] <= u[q],                                          (1)
*
*  where l[q] < u[q], having the only non-zero coefficient in row p,
*  which is equality constraint:
*
*     sum a[p,j] x[j] + a[p,q] x[q] = b.                             (2)
*     j!=q
*
*  PROBLEM TRANSFORMATION
*
*  (If x[q] is integral, this transformation must not be used.)
*
*  The term a[p,q] x[q] in constraint (2) can be considered as a slack
*  variable that allows to carry bounds of column q over row p and then
*  remove column q from the problem.
*
*  Constraint (2) can be written as follows:
*
*     sum a[p,j] x[j] = b - a[p,q] x[q].                             (3)
*     j!=q
*
*  According to (1) constraint (3) is equivalent to the following
*  inequality constraint:
*
*     L[p] <= sum a[p,j] x[j] <= U[p],                               (4)
*             j!=q
*
*  where
*
*            ( b - a[p,q] u[q],  if a[p,q] > 0
*     L[p] = <                                                       (5)
*            ( b - a[p,q] l[q],  if a[p,q] < 0
*
*            ( b - a[p,q] l[q],  if a[p,q] > 0
*     U[p] = <                                                       (6)
*            ( b - a[p,q] u[q],  if a[p,q] < 0
*
*  From (2) it follows that:
*
*              1
*     x[q] = ------ (b - sum a[p,j] x[j]).                           (7)
*            a[p,q]      j!=q
*
*  In order to eliminate x[q] from the objective row we substitute it
*  from (6) to that row:
*
*     z = sum c[j] x[j] + c[q] x[q] + c[0] =
*         j!=q
*                                 1
*       = sum c[j] x[j] + c[q] [------ (b - sum a[p,j] x[j])] + c0 =
*         j!=q                  a[p,q]      j!=q
*
*       = sum c~[j] x[j] + c~[0],
*         j!=q
*                         a[p,j]                     b
*     c~[j] = c[j] - c[q] ------,  c~0 = c0 - c[q] ------            (8)
*                         a[p,q]                   a[p,q]
*
*  are values of objective coefficients and constant term, resp., in
*  the transformed problem.
*
*  Note that column q is column singleton, so in the dual system of the
*  original problem it corresponds to the following row singleton:
*
*     a[p,q] pi[p] + lambda[q] = c[q].                               (9)
*
*  In the transformed problem row (9) would be the following:
*
*     a[p,q] pi~[p] + lambda[q] = c~[q] = 0.                        (10)
*
*  Subtracting (10) from (9) we have:
*
*     a[p,q] (pi[p] - pi~[p]) = c[q]
*
*  that gives the following formula to compute multiplier for row p in
*  solution to the original problem using its value in solution to the
*  transformed problem:
*
*     pi[p] = pi~[p] + c[q] / a[p,q].                               (11)
*
*  RECOVERING BASIC SOLUTION
*
*  Status of column q in solution to the original problem is defined
*  by status of row p in solution to the transformed problem and the
*  sign of coefficient a[p,q] in the original inequality constraint (2)
*  as follows:
*
*     +-----------------------+---------+--------------------+
*     |    Status of row p    | Sign of | Status of column q |
*     | (transformed problem) | a[p,q]  | (original problem) |
*     +-----------------------+---------+--------------------+
*     |        GLP_BS         |  + / -  |       GLP_BS       |
*     |        GLP_NL         |    +    |       GLP_NU       |
*     |        GLP_NL         |    -    |       GLP_NL       |
*     |        GLP_NU         |    +    |       GLP_NL       |
*     |        GLP_NU         |    -    |       GLP_NU       |
*     |        GLP_NF         |  + / -  |       GLP_NF       |
*     +-----------------------+---------+--------------------+
*
*  Value of column q is computed with formula (7). Since originally row
*  p is equality constraint, its status is assigned GLP_NS, and value of
*  its multiplier pi[p] is computed with formula (11).
*
*  RECOVERING INTERIOR-POINT SOLUTION
*
*  Value of column q is computed with formula (7). Row multiplier value
*  pi[p] is computed with formula (11).
*
*  RECOVERING MIP SOLUTION
*
*  Value of column q is computed with formula (7). */

struct implied_slack
{     /* column singleton (implied slack variable) */
      int p;
      /* row reference number */
      int q;
      /* column reference number */
      double apq;
      /* constraint coefficient a[p,q] */
      double b;
      /* right-hand side of original equality constraint */
      double c;
      /* original objective coefficient at x[q] */
      NPPLFE *ptr;
      /* list of non-zero coefficients a[p,j], j != q */
};

static int rcv_implied_slack(NPP *npp, void *info);

void npp_implied_slack(NPP *npp, NPPCOL *q)
{     /* process column singleton (implied slack variable) */
      struct implied_slack *info;
      NPPROW *p;
      NPPAIJ *aij;
      NPPLFE *lfe;
      /* the column must be non-integral non-fixed singleton */
      xassert(!q->is_int);
      xassert(q->lb < q->ub);
      xassert(q->ptr != NULL && q->ptr->c_next == NULL);
      /* corresponding row must be equality constraint */
      aij = q->ptr;
      p = aij->row;
      xassert(p->lb == p->ub);
      /* create transformation stack entry */
      info = npp_push_tse(npp,
         rcv_implied_slack, sizeof(struct implied_slack));
      info->p = p->i;
      info->q = q->j;
      info->apq = aij->val;
      info->b = p->lb;
      info->c = q->coef;
      info->ptr = NULL;
      /* save row coefficients a[p,j], j != q, and substitute x[q]
         into the objective row */
      for (aij = p->ptr; aij != NULL; aij = aij->r_next)
      {  if (aij->col == q) continue; /* skip a[p,q] */
         lfe = dmp_get_atom(npp->stack, sizeof(NPPLFE));
         lfe->ref = aij->col->j;
         lfe->val = aij->val;
         lfe->next = info->ptr;
         info->ptr = lfe;
         aij->col->coef -= info->c * (aij->val / info->apq);
      }
      npp->c0 += info->c * (info->b / info->apq);
      /* compute new row bounds */
      if (info->apq > 0.0)
      {  p->lb = (q->ub == +DBL_MAX ?
            -DBL_MAX : info->b - info->apq * q->ub);
         p->ub = (q->lb == -DBL_MAX ?
            +DBL_MAX : info->b - info->apq * q->lb);
      }
      else
      {  p->lb = (q->lb == -DBL_MAX ?
            -DBL_MAX : info->b - info->apq * q->lb);
         p->ub = (q->ub == +DBL_MAX ?
            +DBL_MAX : info->b - info->apq * q->ub);
      }
      /* remove the column from the problem */
      npp_del_col(npp, q);
      return;
}

static int rcv_implied_slack(NPP *npp, void *_info)
{     /* recover column singleton (implied slack variable) */
      struct implied_slack *info = _info;
      NPPLFE *lfe;
      double temp;
      if (npp->sol == GLP_SOL)
      {  /* assign statuses to row p and column q */
         if (npp->r_stat[info->p] == GLP_BS ||
             npp->r_stat[info->p] == GLP_NF)
            npp->c_stat[info->q] = npp->r_stat[info->p];
         else if (npp->r_stat[info->p] == GLP_NL)
            npp->c_stat[info->q] =
               (char)(info->apq > 0.0 ? GLP_NU : GLP_NL);
         else if (npp->r_stat[info->p] == GLP_NU)
            npp->c_stat[info->q] =
               (char)(info->apq > 0.0 ? GLP_NL : GLP_NU);
         else
         {  npp_error();
            return 1;
         }
         npp->r_stat[info->p] = GLP_NS;
      }
      if (npp->sol != GLP_MIP)
      {  /* compute multiplier for row p */
         npp->r_pi[info->p] += info->c / info->apq;
      }
      /* compute value of column q */
      temp = info->b;
      for (lfe = info->ptr; lfe != NULL; lfe = lfe->next)
         temp -= lfe->val * npp->c_value[lfe->ref];
      npp->c_value[info->q] = temp / info->apq;
      return 0;
}

/***********************************************************************
*  NAME
*
*  npp_implied_free - process column singleton (implied free variable)
*
*  SYNOPSIS
*
*  #include "glpnpp.h"
*  int npp_implied_free(NPP *npp, NPPCOL *q);
*
*  DESCRIPTION
*
*  The routine npp_implied_free processes column q:
*
*     l[q] <= x[q] <= u[q],                                          (1)
*
*  having non-zero coefficient in the only row p, which is inequality
*  constraint:
*
*     L[p] <= sum a[p,j] x[j] + a[p,q] x[q] <= U[p],                 (2)
*             j!=q
*
*  where l[q] < u[q], L[p] < U[p], L[p] > -oo and/or U[p] < +oo.
*
*  RETURNS
*
*  0 - success;
*
*  1 - column lower and/or upper bound(s) can be active;
*
*  2 - problem has no dual feasible solution.
*
*  PROBLEM TRANSFORMATION
*
*  Constraint (2) can be written as follows:
*
*     L[p] - sum a[p,j] x[j] <= a[p,q] x[q] <= U[p] - sum a[p,j] x[j],
*            j!=q                                     j!=q
*
*  from which it follows that:
*
*     alfa <= a[p,q] x[q] <= beta,                                   (3)
*
*  where
*
*     alfa = inf(L[p] - sum a[p,j] x[j]) =
*                       j!=q
*
*          = L[p] - sup sum a[p,j] x[j] =                            (4)
*                       j!=q
*
*          = L[p] -  sum  a[p,j] u[j] -  sum  a[p,j] l[j],
*                  j in Jp             j in Jn
*
*     beta = sup(L[p] - sum a[p,j] x[j]) =
*                       j!=q
*
*          = L[p] - inf sum a[p,j] x[j] =                            (5)
*                       j!=q
*
*          = L[p] -  sum  a[p,j] l[j] -  sum  a[p,j] u[j],
*                  j in Jp             j in Jn
*
*     Jp = {j != q: a[p,j] > 0},  Jn = {j != q: a[p,j] < 0}.         (6)
*
*  Inequality (3) defines implied bounds of variable x[q]:
*
*     l'[q] <= x[q] <= u'[q],                                        (7)
*
*  where
*
*             ( alfa / a[p,q], if a[p,q] > 0
*     l'[q] = <                                                     (8a)
*             ( beta / a[p,q], if a[p,q] < 0
*
*             ( beta / a[p,q], if a[p,q] > 0
*     u'[q] = <                                                     (8b)
*             ( alfa / a[p,q], if a[p,q] < 0
*
*  Thus, if l'[q] > l[q] - eps and u'[q] < u[q] + eps, where eps is
*  an absolute tolerance for column value, column bounds (1) cannot be
*  active, in which case column q can be replaced by equivalent free
*  (unbounded) column.
*
*  Note that column q is column singleton, so in the dual system of the
*  original problem it corresponds to the following row singleton:
*
*     a[p,q] pi[p] + lambda[q] = c[q],                               (9)
*
*  from which it follows that:
*
*     pi[p] = (c[q] - lambda[q]) / a[p,q].                          (10)
*
*  Let x[q] be implied free (unbounded) variable. Then column q can be
*  only basic, so its multiplier lambda[q] is equal to zero, and from
*  (10) we have:
*
*     pi[p] = c[q] / a[p,q].                                        (11)
*
*  There are possible three cases:
*
*  1) pi[p] < -eps, where eps is an absolute tolerance for row
*     multiplier. In this case, to provide dual feasibility of the
*     original problem, row p must be active on its lower bound, and
*     if its lower bound does not exist (L[p] = -oo), the problem has
*     no dual feasible solution;
*
*  2) pi[p] > +eps. In this case row p must be active on its upper
*     bound, and if its upper bound does not exist (U[p] = +oo), the
*     problem has no dual feasible solution;
*
*  3) -eps <= pi[p] <= +eps. In this case any (either lower or upper)
*     bound of row p can be active, because this does not affect dual
*     feasibility.
*
*  Thus, in all three cases original inequality constraint (2) can be
*  replaced by equality constraint, where the right-hand side is either
*  lower or upper bound of row p, and bounds of column q can be removed
*  that makes it free (unbounded). (May note that this transformation
*  can be followed by transformation "Column singleton (implied slack
*  variable)" performed by the routine npp_implied_slack.)
*
*  RECOVERING BASIC SOLUTION
*
*  Status of row p in solution to the original problem is determined
*  by its status in solution to the transformed problem and its bound,
*  which was choosen to be active:
*
*     +-----------------------+--------+--------------------+
*     |    Status of row p    | Active | Status of row p    |
*     | (transformed problem) | bound  | (original problem) |
*     +-----------------------+--------+--------------------+
*     |        GLP_BS         |  L[p]  |       GLP_BS       |
*     |        GLP_BS         |  U[p]  |       GLP_BS       |
*     |        GLP_NS         |  L[p]  |       GLP_NL       |
*     |        GLP_NS         |  U[p]  |       GLP_NU       |
*     +-----------------------+--------+--------------------+
*
*  Value of row multiplier pi[p] (as well as value of column q) in
*  solution to the original problem is the same as in solution to the
*  transformed problem.
*
*  RECOVERING INTERIOR-POINT SOLUTION
*
*  Value of row multiplier pi[p] in solution to the original problem is
*  the same as in solution to the transformed problem.
*
*  RECOVERING MIP SOLUTION
*
*  None needed. */

struct implied_free
{     /* column singleton (implied free variable) */
      int p;
      /* row reference number */
      char stat;
      /* row status:
         GLP_NL - active constraint on lower bound
         GLP_NU - active constraint on upper bound */
};

static int rcv_implied_free(NPP *npp, void *info);

int npp_implied_free(NPP *npp, NPPCOL *q)
{     /* process column singleton (implied free variable) */
      struct implied_free *info;
      NPPROW *p;
      NPPAIJ *apq, *aij;
      double alfa, beta, l, u, pi, eps;
      /* the column must be non-fixed singleton */
      xassert(q->lb < q->ub);
      xassert(q->ptr != NULL && q->ptr->c_next == NULL);
      /* corresponding row must be inequality constraint */
      apq = q->ptr;
      p = apq->row;
      xassert(p->lb != -DBL_MAX || p->ub != +DBL_MAX);
      xassert(p->lb < p->ub);
      /* compute alfa */
      alfa = p->lb;
      if (alfa != -DBL_MAX)
      {  for (aij = p->ptr; aij != NULL; aij = aij->r_next)
         {  if (aij == apq) continue; /* skip a[p,q] */
            if (aij->val > 0.0)
            {  if (aij->col->ub == +DBL_MAX)
               {  alfa = -DBL_MAX;
                  break;
               }
               alfa -= aij->val * aij->col->ub;
            }
            else /* < 0.0 */
            {  if (aij->col->lb == -DBL_MAX)
               {  alfa = -DBL_MAX;
                  break;
               }
               alfa -= aij->val * aij->col->lb;
            }
         }
      }
      /* compute beta */
      beta = p->ub;
      if (beta != +DBL_MAX)
      {  for (aij = p->ptr; aij != NULL; aij = aij->r_next)
         {  if (aij == apq) continue; /* skip a[p,q] */
            if (aij->val > 0.0)
            {  if (aij->col->lb == -DBL_MAX)
               {  beta = +DBL_MAX;
                  break;
               }
               beta -= aij->val * aij->col->lb;
            }
            else /* < 0.0 */
            {  if (aij->col->ub == +DBL_MAX)
               {  beta = +DBL_MAX;
                  break;
               }
               beta -= aij->val * aij->col->ub;
            }
         }
      }
      /* compute implied column lower bound l'[q] */
      if (apq->val > 0.0)
         l = (alfa == -DBL_MAX ? -DBL_MAX : alfa / apq->val);
      else /* < 0.0 */
         l = (beta == +DBL_MAX ? -DBL_MAX : beta / apq->val);
      /* compute implied column upper bound u'[q] */
      if (apq->val > 0.0)
         u = (beta == +DBL_MAX ? +DBL_MAX : beta / apq->val);
      else
         u = (alfa == -DBL_MAX ? +DBL_MAX : alfa / apq->val);
      /* check if column lower bound l[q] can be active */
      if (q->lb != -DBL_MAX)
      {  eps = 1e-9 + 1e-12 * fabs(q->lb);
         if (l < q->lb - eps) return 1; /* yes, it can */
      }
      /* check if column upper bound u[q] can be active */
      if (q->ub != +DBL_MAX)
      {  eps = 1e-9 + 1e-12 * fabs(q->ub);
         if (u > q->ub + eps) return 1; /* yes, it can */
      }
      /* okay; make column q free (unbounded) */
      q->lb = -DBL_MAX, q->ub = +DBL_MAX;
      /* create transformation stack entry */
      info = npp_push_tse(npp,
         rcv_implied_free, sizeof(struct implied_free));
      info->p = p->i;
      info->stat = -1;
      /* compute row multiplier pi[p] */
      pi = q->coef / apq->val;
      /* check dual feasibility for row p */
      if (pi > +DBL_EPSILON)
      {  /* lower bound L[p] must be active */
         if (p->lb != -DBL_MAX)
nl:      {  info->stat = GLP_NL;
            p->ub = p->lb;
         }
         else
         {  if (pi > +1e-5) return 2; /* dual infeasibility */
            /* take a chance on U[p] */
            xassert(p->ub != +DBL_MAX);
            goto nu;
         }
      }
      else if (pi < -DBL_EPSILON)
      {  /* upper bound U[p] must be active */
         if (p->ub != +DBL_MAX)
nu:      {  info->stat = GLP_NU;
            p->lb = p->ub;
         }
         else
         {  if (pi < -1e-5) return 2; /* dual infeasibility */
            /* take a chance on L[p] */
            xassert(p->lb != -DBL_MAX);
            goto nl;
         }
      }
      else
      {  /* any bound (either L[p] or U[p]) can be made active  */
         if (p->ub == +DBL_MAX)
         {  xassert(p->lb != -DBL_MAX);
            goto nl;
         }
         if (p->lb == -DBL_MAX)
         {  xassert(p->ub != +DBL_MAX);
            goto nu;
         }
         if (fabs(p->lb) <= fabs(p->ub)) goto nl; else goto nu;
      }
      return 0;
}

static int rcv_implied_free(NPP *npp, void *_info)
{     /* recover column singleton (implied free variable) */
      struct implied_free *info = _info;
      if (npp->sol == GLP_SOL)
      {  if (npp->r_stat[info->p] == GLP_BS)
            npp->r_stat[info->p] = GLP_BS;
         else if (npp->r_stat[info->p] == GLP_NS)
         {  xassert(info->stat == GLP_NL || info->stat == GLP_NU);
            npp->r_stat[info->p] = info->stat;
         }
         else
         {  npp_error();
            return 1;
         }
      }
      return 0;
}

/***********************************************************************
*  NAME
*
*  npp_eq_doublet - process row doubleton (equality constraint)
*
*  SYNOPSIS
*
*  #include "glpnpp.h"
*  NPPCOL *npp_eq_doublet(NPP *npp, NPPROW *p);
*
*  DESCRIPTION
*
*  The routine npp_eq_doublet processes row p, which is equality
*  constraint having exactly two non-zero coefficients:
*
*     a[p,q] x[q] + a[p,r] x[r] = b.                                 (1)
*
*  As the result of processing one of columns q or r is eliminated from
*  all other rows and, thus, becomes column singleton of type "implied
*  slack variable". Row p is not changed and along with column q and r
*  remains in the problem.
*
*  RETURNS
*
*  The routine npp_eq_doublet returns pointer to the descriptor of that
*  column q or r which has been eliminated. If, due to some reason, the
*  elimination was not performed, the routine returns NULL.
*
*  PROBLEM TRANSFORMATION
*
*  First, we decide which column q or r will be eliminated. Let it be
*  column q. Consider i-th constraint row, where column q has non-zero
*  coefficient a[i,q] != 0:
*
*     L[i] <= sum a[i,j] x[j] <= U[i].                               (2)
*              j
*
*  In order to eliminate column q from row (2) we subtract from it row
*  (1) multiplied by gamma[i] = a[i,q] / a[p,q], i.e. we replace in the
*  transformed problem row (2) by its linear combination with row (1).
*  This transformation changes only coefficients in columns q and r,
*  and bounds of row i as follows:
*
*     a~[i,q] = a[i,q] - gamma[i] a[p,q] = 0,                        (3)
*
*     a~[i,r] = a[i,r] - gamma[i] a[p,r],                            (4)
*
*       L~[i] = L[i] - gamma[i] b,                                   (5)
*
*       U~[i] = U[i] - gamma[i] b.                                   (6)
*
*  RECOVERING BASIC SOLUTION
*
*  The transformation of the primal system of the original problem:
*
*     L <= A x <= U                                                  (7)
*
*  is equivalent to multiplying from the left a transformation matrix F
*  by components of this primal system, which in the transformed problem
*  becomes the following:
*
*     F L <= F A x <= F U  ==>  L~ <= A~x <= U~.                     (8)
*
*  The matrix F has the following structure:
*
*         ( 1           -gamma[1]            )
*         (                                  )
*         (    1        -gamma[2]            )
*         (                                  )
*         (      ...       ...               )
*         (                                  )
*     F = (          1  -gamma[p-1]          )                       (9)
*         (                                  )
*         (                 1                )
*         (                                  )
*         (             -gamma[p+1]  1       )
*         (                                  )
*         (                ...          ...  )
*
*  where its column containing elements -gamma[i] corresponds to row p
*  of the primal system.
*
*  From (8) it follows that the dual system of the original problem:
*
*     A'pi + lambda = c,                                            (10)
*
*  in the transformed problem becomes the following:
*
*     A'F'inv(F')pi + lambda = c  ==>  (A~)'pi~ + lambda = c,       (11)
*
*  where:
*
*     pi~ = inv(F')pi                                               (12)
*
*  is the vector of row multipliers in the transformed problem. Thus:
*
*     pi = F'pi~.                                                   (13)
*
*  Therefore, as it follows from (13), value of multiplier for row p in
*  solution to the original problem can be computed as follows:
*
*     pi[p] = pi~[p] - sum gamma[i] pi~[i],                         (14)
*                       i
*
*  where pi~[i] = pi[i] is multiplier for row i (i != p).
*
*  Note that the statuses of all rows and columns are not changed.
*
*  RECOVERING INTERIOR-POINT SOLUTION
*
*  Multiplier for row p in solution to the original problem is computed
*  with formula (14).
*
*  RECOVERING MIP SOLUTION
*
*  None needed. */

struct eq_doublet
{     /* row doubleton (equality constraint) */
      int p;
      /* row reference number */
      double apq;
      /* constraint coefficient a[p,q] */
      NPPLFE *ptr;
      /* list of non-zero coefficients a[i,q], i != p */
};

static int rcv_eq_doublet(NPP *npp, void *info);

NPPCOL *npp_eq_doublet(NPP *npp, NPPROW *p)
{     /* process row doubleton (equality constraint) */
      struct eq_doublet *info;
      NPPROW *i;
      NPPCOL *q, *r;
      NPPAIJ *apq, *apr, *aiq, *air, *next;
      NPPLFE *lfe;
      double gamma;
      /* the row must be doubleton equality constraint */
      xassert(p->lb == p->ub);
      xassert(p->ptr != NULL && p->ptr->r_next != NULL &&
              p->ptr->r_next->r_next == NULL);
      /* choose column to be eliminated */
      {  NPPAIJ *a1, *a2;
         a1 = p->ptr, a2 = a1->r_next;
         if (fabs(a2->val) < 0.001 * fabs(a1->val))
         {  /* only first column can be eliminated, because second one
               has too small constraint coefficient */
            apq = a1, apr = a2;
         }
         else if (fabs(a1->val) < 0.001 * fabs(a2->val))
         {  /* only second column can be eliminated, because first one
               has too small constraint coefficient */
            apq = a2, apr = a1;
         }
         else
         {  /* both columns are appropriate; choose that one which is
               shorter to minimize fill-in */
            if (npp_col_nnz(npp, a1->col) <= npp_col_nnz(npp, a2->col))
            {  /* first column is shorter */
               apq = a1, apr = a2;
            }
            else
            {  /* second column is shorter */
               apq = a2, apr = a1;
            }
         }
      }
      /* now columns q and r have been chosen */
      q = apq->col, r = apr->col;
      /* create transformation stack entry */
      info = npp_push_tse(npp,
         rcv_eq_doublet, sizeof(struct eq_doublet));
      info->p = p->i;
      info->apq = apq->val;
      info->ptr = NULL;
      /* transform each row i (i != p), where a[i,q] != 0, to eliminate
         column q */
      for (aiq = q->ptr; aiq != NULL; aiq = next)
      {  next = aiq->c_next;
         if (aiq == apq) continue; /* skip row p */
         i = aiq->row; /* row i to be transformed */
         /* save constraint coefficient a[i,q] */
         if (npp->sol != GLP_MIP)
         {  lfe = dmp_get_atom(npp->stack, sizeof(NPPLFE));
            lfe->ref = i->i;
            lfe->val = aiq->val;
            lfe->next = info->ptr;
            info->ptr = lfe;
         }
         /* find coefficient a[i,r] in row i */
         for (air = i->ptr; air != NULL; air = air->r_next)
            if (air->col == r) break;
         /* if a[i,r] does not exist, create a[i,r] = 0 */
         if (air == NULL)
            air = npp_add_aij(npp, i, r, 0.0);
         /* compute gamma[i] = a[i,q] / a[p,q] */
         gamma = aiq->val / apq->val;
         /* (row i) := (row i) - gamma[i] * (row p); see (3)-(6) */
         /* new a[i,q] is exact zero due to elimnation; remove it from
            row i */
         npp_del_aij(npp, aiq);
         /* compute new a[i,r] */
         air->val -= gamma * apr->val;
         /* if new a[i,r] is close to zero due to numeric cancelation,
            remove it from row i */
         if (fabs(air->val) <= 1e-10)
            npp_del_aij(npp, air);
         /* compute new lower and upper bounds of row i */
         if (i->lb == i->ub)
            i->lb = i->ub = (i->lb - gamma * p->lb);
         else
         {  if (i->lb != -DBL_MAX)
               i->lb -= gamma * p->lb;
            if (i->ub != +DBL_MAX)
               i->ub -= gamma * p->lb;
         }
      }
      return q;
}

static int rcv_eq_doublet(NPP *npp, void *_info)
{     /* recover row doubleton (equality constraint) */
      struct eq_doublet *info = _info;
      NPPLFE *lfe;
      double gamma, temp;
      /* we assume that processing row p is followed by processing
         column q as singleton of type "implied slack variable", in
         which case row p must always be active equality constraint */
      if (npp->sol == GLP_SOL)
      {  if (npp->r_stat[info->p] != GLP_NS)
         {  npp_error();
            return 1;
         }
      }
      if (npp->sol != GLP_MIP)
      {  /* compute value of multiplier for row p; see (14) */
         temp = npp->r_pi[info->p];
         for (lfe = info->ptr; lfe != NULL; lfe = lfe->next)
         {  gamma = lfe->val / info->apq; /* a[i,q] / a[p,q] */
            temp -= gamma * npp->r_pi[lfe->ref];
         }
         npp->r_pi[info->p] = temp;
      }
      return 0;
}

/***********************************************************************
*  NAME
*
*  npp_forcing_row - process forcing row
*
*  SYNOPSIS
*
*  #include "glpnpp.h"
*  int npp_forcing_row(NPP *npp, NPPROW *p, int at);
*
*  DESCRIPTION
*
*  The routine npp_forcing row processes row p of general format:
*
*     L[p] <= sum a[p,j] x[j] <= U[p],                               (1)
*              j
*
*     l[j] <= x[j] <= u[j],                                          (2)
*
*  where L[p] <= U[p] and l[j] < u[j] for all a[p,j] != 0. It is also
*  assumed that:
*
*  1) if at = 0 then |L[p] - U'[p]| <= eps, where U'[p] is implied
*     row upper bound (see below), eps is an absolute tolerance for row
*     value;
*
*  2) if at = 1 then |U[p] - L'[p]| <= eps, where L'[p] is implied
*     row lower bound (see below).
*
*  RETURNS
*
*  0 - success;
*
*  1 - cannot fix columns due to too small constraint coefficients.
*
*  PROBLEM TRANSFORMATION
*
*  Implied lower and upper bounds of row (1) are determined by bounds
*  of corresponding columns (variables) as follows:
*
*     L'[p] = inf sum a[p,j] x[j] =
*                  j
*                                                                    (3)
*           =  sum  a[p,j] l[j] +  sum  a[p,j] u[j],
*            j in Jp             j in Jn
*
*     U'[p] = sup sum a[p,j] x[j] =
*                                                                    (4)
*           =  sum  a[p,j] u[j] +  sum  a[p,j] l[j],
*            j in Jp             j in Jn
*
*     Jp = {j: a[p,j] > 0},  Jn = {j: a[p,j] < 0}.                   (5)
*
*  If L[p] =~ U'[p] (at = 0), solution can be primal feasible only when
*  all variables take their boundary values as defined by (4):
*
*            ( u[j], if j in Jp
*     x[j] = <                                                       (6)
*            ( l[j], if j in Jn
*
*  Similarly, if U[p] =~ L'[p] (at = 1), solution can be primal feasible
*  only when all variables take their boundary values as defined by (3):
*
*            ( l[j], if j in Jp
*     x[j] = <                                                       (7)
*            ( u[j], if j in Jn
*
*  Condition (6) or (7) allows fixing all columns (variables x[j])
*  in row (1) on their bounds and then removing them from the problem
*  (see the routine npp_fixed_col). Due to this row p becomes redundant,
*  so it can be replaced by equivalent free (unbounded) row and also
*  removed from the problem (see the routine npp_free_row).
*
*  1. To apply this transformation row (1) should not have coefficients
*     whose magnitude is too small, i.e. all a[p,j] should satisfy to
*     the following condition:
*
*        |a[p,j]| >= eps * max(1, |a[p,k]|),                         (8)
*                           k
*     where eps is a relative tolerance for constraint coefficients.
*     Otherwise, fixing columns may be numerically unreliable and may
*     lead to wrong solution.
*
*  2. The routine fixes columns and remove bounds of row p, however,
*     it does not remove the row and columns from the problem.
*
*  RECOVERING BASIC SOLUTION
*
*  In the transformed problem row p being inactive constraint is
*  assigned status GLP_BS (as the result of transformation of free
*  row), and all columns in this row are assigned status GLP_NS (as the
*  result of transformation of fixed columns).
*
*  Note that in the dual system of the transformed (as well as original)
*  problem every column j in row p corresponds to the following row:
*
*     sum  a[i,j] pi[i] + a[p,j] pi[p] + lambda[j] = c[j],           (9)
*     i!=p
*
*  from which it follows that:
*
*     lambda[j] = c[j] - sum a[i,j] pi[i] - a[p,j] pi[p].           (10)
*                        i!=p
*
*  In the transformed problem values of all multipliers pi[i] are known
*  (including pi[i], whose value is zero, since row p is inactive).
*  Thus, using formula (10) it is possible to compute values of
*  multipliers lambda[j] for all columns in row p.
*
*  Note also that in the original problem all columns in row p are
*  bounded, not fixed. So status GLP_NS assigned to every such column
*  must be changed to GLP_NL or GLP_NU depending on which bound the
*  corresponding column has been fixed. This status change may lead to
*  dual feasibility violation for solution of the original problem,
*  because now column multipliers must satisfy to the following
*  condition:
*
*               ( >= 0, if status of column j is GLP_NL,
*     lambda[j] <                                                   (11)
*               ( <= 0, if status of column j is GLP_NU.
*
*  If this condition holds, solution to the original problem is the
*  same as to the transformed problem. Otherwise, we have to perform
*  one degenerate pivoting step of the primal simplex method to obtain
*  dual feasible (hence, optimal) solution to the original problem as
*  follows. If, on problem transformation, row p was made active on its
*  lower bound (case at = 0), we change its status to GLP_NL (or GLP_NS)
*  and start increasing its multiplier pi[p]. Otherwise, if row p was
*  made active on its upper bound (case at = 1), we change its status
*  to GLP_NU (or GLP_NS) and start decreasing pi[p]. From (10) it
*  follows that:
*
*     delta lambda[j] = - a[p,j] * delta pi[p] = - a[p,j] pi[p].    (12)
*
*  Simple analysis of formulae (3)-(5) shows that changing pi[p] in the
*  specified direction causes increasing lambda[j] for every column j
*  assigned status GLP_NL (delta lambda[j] > 0) and decreasing lambda[j]
*  for every column j assigned status GLP_NU (delta lambda[j] < 0). It
*  is understood that once the last lambda[q], which violates condition
*  (11), has reached zero, multipliers lambda[j] for all columns get
*  valid signs. Such column q can be determined as follows. Let d[j] be
*  initial value of lambda[j] (i.e. reduced cost of column j) in the
*  transformed problem computed with formula (10) when pi[p] = 0. Then
*  lambda[j] = d[j] + delta lambda[j], and from (12) it follows that
*  lambda[j] becomes zero if:
*
*     delta lambda[j] = - a[p,j] pi[p] = - d[j]  ==>
*                                                                   (13)
*     pi[p] = d[j] / a[p,j].
*
*  Therefore, the last column q, for which lambda[q] becomes zero, can
*  be determined from the following condition:
*
*     |d[q] / a[p,q]| = max  |pi[p]| = max  |d[j] / a[p,j]|,        (14)
*                      j in D         j in D
*
*  where D is a set of columns j whose, reduced costs d[j] have invalid
*  signs, i.e. violate condition (11). (Thus, if D is empty, solution
*  to the original problem is the same as solution to the transformed
*  problem, and no correction is needed as was noticed above.) In
*  solution to the original problem column q is assigned status GLP_BS,
*  since it replaces column of auxiliary variable of row p (becoming
*  active) in the basis, and multiplier for row p is assigned its new
*  value, which is pi[p] = d[q] / a[p,q]. Note that due to primal
*  degeneracy values of all columns having non-zero coefficients in row
*  p remain unchanged.
*
*  RECOVERING INTERIOR-POINT SOLUTION
*
*  Value of multiplier pi[p] in solution to the original problem is
*  corrected in the same way as for basic solution. Values of all
*  columns having non-zero coefficients in row p remain unchanged.
*
*  RECOVERING MIP SOLUTION
*
*  None needed. */

struct forcing_col
{     /* column fixed on its bound by forcing row */
      int j;
      /* column reference number */
      char stat;
      /* original column status:
         GLP_NL - fixed on lower bound
         GLP_NU - fixed on upper bound */
      double a;
      /* constraint coefficient a[p,j] */
      double c;
      /* objective coefficient c[j] */
      NPPLFE *ptr;
      /* list of non-zero coefficients a[i,j], i != p */
      struct forcing_col *next;
      /* pointer to another column fixed by forcing row */
};

struct forcing_row
{     /* forcing row */
      int p;
      /* row reference number */
      char stat;
      /* status assigned to the row if it becomes active:
         GLP_NS - active equality constraint
         GLP_NL - inequality constraint with lower bound active
         GLP_NU - inequality constraint with upper bound active */
      struct forcing_col *ptr;
      /* list of all columns having non-zero constraint coefficient
         a[p,j] in the forcing row */
};

static int rcv_forcing_row(NPP *npp, void *info);

int npp_forcing_row(NPP *npp, NPPROW *p, int at)
{     /* process forcing row */
      struct forcing_row *info;
      struct forcing_col *col = NULL;
      NPPCOL *j;
      NPPAIJ *apj, *aij;
      NPPLFE *lfe;
      double big;
      xassert(at == 0 || at == 1);
      /* determine maximal magnitude of the row coefficients */
      big = 1.0;
      for (apj = p->ptr; apj != NULL; apj = apj->r_next)
         if (big < fabs(apj->val)) big = fabs(apj->val);
      /* if there are too small coefficients in the row, transformation
         should not be applied */
      for (apj = p->ptr; apj != NULL; apj = apj->r_next)
         if (fabs(apj->val) < 1e-7 * big) return 1;
      /* create transformation stack entry */
      info = npp_push_tse(npp,
         rcv_forcing_row, sizeof(struct forcing_row));
      info->p = p->i;
      if (p->lb == p->ub)
      {  /* equality constraint */
         info->stat = GLP_NS;
      }
      else if (at == 0)
      {  /* inequality constraint; case L[p] = U'[p] */
         info->stat = GLP_NL;
         xassert(p->lb != -DBL_MAX);
      }
      else /* at == 1 */
      {  /* inequality constraint; case U[p] = L'[p] */
         info->stat = GLP_NU;
         xassert(p->ub != +DBL_MAX);
      }
      info->ptr = NULL;
      /* scan the forcing row, fix columns at corresponding bounds, and
         save column information (the latter is not needed for MIP) */
      for (apj = p->ptr; apj != NULL; apj = apj->r_next)
      {  /* column j has non-zero coefficient in the forcing row */
         j = apj->col;
         /* it must be non-fixed */
         xassert(j->lb < j->ub);
         /* allocate stack entry to save column information */
         if (npp->sol != GLP_MIP)
         {  col = dmp_get_atom(npp->stack, sizeof(struct forcing_col));
            col->j = j->j;
            col->stat = -1; /* will be set below */
            col->a = apj->val;
            col->c = j->coef;
            col->ptr = NULL;
            col->next = info->ptr;
            info->ptr = col;
         }
         /* fix column j */
         if (at == 0 && apj->val < 0.0 || at != 0 && apj->val > 0.0)
         {  /* at its lower bound */
            if (npp->sol != GLP_MIP)
               col->stat = GLP_NL;
            xassert(j->lb != -DBL_MAX);
            j->ub = j->lb;
         }
         else
         {  /* at its upper bound */
            if (npp->sol != GLP_MIP)
               col->stat = GLP_NU;
            xassert(j->ub != +DBL_MAX);
            j->lb = j->ub;
         }
         /* save column coefficients a[i,j], i != p */
         if (npp->sol != GLP_MIP)
         {  for (aij = j->ptr; aij != NULL; aij = aij->c_next)
            {  if (aij == apj) continue; /* skip a[p,j] */
               lfe = dmp_get_atom(npp->stack, sizeof(NPPLFE));
               lfe->ref = aij->row->i;
               lfe->val = aij->val;
               lfe->next = col->ptr;
               col->ptr = lfe;
            }
         }
      }
      /* make the row free (unbounded) */
      p->lb = -DBL_MAX, p->ub = +DBL_MAX;
      return 0;
}

static int rcv_forcing_row(NPP *npp, void *_info)
{     /* recover forcing row */
      struct forcing_row *info = _info;
      struct forcing_col *col, *piv;
      NPPLFE *lfe;
      double d, big, temp;
      if (npp->sol == GLP_MIP) goto done;
      /* initially solution to the original problem is the same as
         to the transformed problem, where row p is inactive constraint
         with pi[p] = 0, and all columns are non-basic */
      if (npp->sol == GLP_SOL)
      {  if (npp->r_stat[info->p] != GLP_BS)
         {  npp_error();
            return 1;
         }
         for (col = info->ptr; col != NULL; col = col->next)
         {  if (npp->c_stat[col->j] != GLP_NS)
            {  npp_error();
               return 1;
            }
            npp->c_stat[col->j] = col->stat; /* original status */
         }
      }
      /* compute reduced costs d[j] for all columns with formula (10)
         and store them in col.c instead objective coefficients */
      for (col = info->ptr; col != NULL; col = col->next)
      {  d = col->c;
         for (lfe = col->ptr; lfe != NULL; lfe = lfe->next)
            d -= lfe->val * npp->r_pi[lfe->ref];
         col->c = d;
      }
      /* consider columns j, whose multipliers lambda[j] has wrong
         sign in solution to the transformed problem (where lambda[j] =
         d[j]), and choose column q, whose multipler lambda[q] reaches
         zero last on changing row multiplier pi[p]; see (14) */
      piv = NULL, big = 0.0;
      for (col = info->ptr; col != NULL; col = col->next)
      {  d = col->c; /* d[j] */
         temp = fabs(d / col->a);
         if (col->stat == GLP_NL)
         {  /* column j has active lower bound */
            if (d < 0.0 && big < temp)
               piv = col, big = temp;
         }
         else if (col->stat == GLP_NU)
         {  /* column j has active upper bound */
            if (d > 0.0 && big < temp)
               piv = col, big = temp;
         }
         else
         {  npp_error();
            return 1;
         }
      }
      /* if column q does not exist, no correction is needed */
      if (piv != NULL)
      {  /* correct solution; row p becomes active constraint while
            column q becomes basic */
         if (npp->sol == GLP_SOL)
         {  npp->r_stat[info->p] = info->stat;
            npp->c_stat[piv->j] = GLP_BS;
         }
         /* assign new value to row multiplier pi[p] = d[p] / a[p,q] */
         npp->r_pi[info->p] = piv->c / piv->a;
      }
done: return 0;
}

/***********************************************************************
*  NAME
*
*  npp_analyze_row - perform general row analysis
*
*  SYNOPSIS
*
*  #include "glpnpp.h"
*  int npp_analyze_row(NPP *npp, NPPROW *p);
*
*  DESCRIPTION
*
*  The routine npp_analyze_row performs analysis of row p of general
*  format:
*
*     L[p] <= sum a[p,j] x[j] <= U[p],                               (1)
*              j
*
*     l[j] <= x[j] <= u[j],                                          (2)
*
*  where L[p] <= U[p] and l[j] <= u[j] for all a[p,j] != 0.
*
*  RETURNS
*
*  0x?0 - row lower bound does not exist or is redundant;
*
*  0x?1 - row lower bound can be active;
*
*  0x?2 - row lower bound is a forcing bound;
*
*  0x0? - row upper bound does not exist or is redundant;
*
*  0x1? - row upper bound can be active;
*
*  0x2? - row upper bound is a forcing bound;
*
*  0x33 - row bounds are inconsistent with column bounds.
*
*  ALGORITHM
*
*  Analysis of row (1) is based on analysis of its implied lower and
*  upper bounds, which are determined by bounds of corresponding columns
*  (variables) as follows:
*
*     L'[p] = inf sum a[p,j] x[j] =
*                  j
*                                                                    (3)
*           =  sum  a[p,j] l[j] +  sum  a[p,j] u[j],
*            j in Jp             j in Jn
*
*     U'[p] = sup sum a[p,j] x[j] =
*                                                                    (4)
*           =  sum  a[p,j] u[j] +  sum  a[p,j] l[j],
*            j in Jp             j in Jn
*
*     Jp = {j: a[p,j] > 0},  Jn = {j: a[p,j] < 0}.                   (5)
*
*  (Note that bounds of all columns in row p are assumed to be correct,
*  so L'[p] <= U'[p].)
*
*  Analysis of row lower bound L[p] includes the following cases:
*
*  1) if L[p] > U'[p] + eps, where eps is an absolute tolerance for row
*     value, row lower bound L[p] and implied row upper bound U'[p] are
*     inconsistent, ergo, the problem has no primal feasible solution;
*
*  2) if U'[p] - eps <= L[p] <= U'[p] + eps, i.e. if L[p] =~ U'[p],
*     the row is a forcing row on its lower bound (see description of
*     the routine npp_forcing_row);
*
*  3) if L[p] > L'[p] + eps, row lower bound L[p] can be active (this
*     conclusion does not account other rows in the problem);
*
*  4) if L[p] <= L'[p] + eps, row lower bound L[p] cannot be active, so
*     it is redundant and can be removed (replaced by -oo).
*
*  Analysis of row upper bound U[p] is performed in a similar way and
*  includes the following cases:
*
*  1) if U[p] < L'[p] - eps, row upper bound U[p] and implied row lower
*     bound L'[p] are inconsistent, ergo the problem has no primal
*     feasible solution;
*
*  2) if L'[p] - eps <= U[p] <= L'[p] + eps, i.e. if U[p] =~ L'[p],
*     the row is a forcing row on its upper bound (see description of
*     the routine npp_forcing_row);
*
*  3) if U[p] < U'[p] - eps, row upper bound U[p] can be active (this
*     conclusion does not account other rows in the problem);
*
*  4) if U[p] >= U'[p] - eps, row upper bound U[p] cannot be active, so
*     it is redundant and can be removed (replaced by +oo). */

int npp_analyze_row(NPP *npp, NPPROW *p)
{     /* perform general row analysis */
      NPPAIJ *aij;
      int ret = 0x00;
      double l, u, eps;
      xassert(npp == npp);
      /* compute implied lower bound L'[p]; see (3) */
      l = 0.0;
      for (aij = p->ptr; aij != NULL; aij = aij->r_next)
      {  if (aij->val > 0.0)
         {  if (aij->col->lb == -DBL_MAX)
            {  l = -DBL_MAX;
               break;
            }
            l += aij->val * aij->col->lb;
         }
         else /* aij->val < 0.0 */
         {  if (aij->col->ub == +DBL_MAX)
            {  l = -DBL_MAX;
               break;
            }
            l += aij->val * aij->col->ub;
         }
      }
      /* compute implied upper bound U'[p]; see (4) */
      u = 0.0;
      for (aij = p->ptr; aij != NULL; aij = aij->r_next)
      {  if (aij->val > 0.0)
         {  if (aij->col->ub == +DBL_MAX)
            {  u = +DBL_MAX;
               break;
            }
            u += aij->val * aij->col->ub;
         }
         else /* aij->val < 0.0 */
         {  if (aij->col->lb == -DBL_MAX)
            {  u = +DBL_MAX;
               break;
            }
            u += aij->val * aij->col->lb;
         }
      }
      /* column bounds are assumed correct, so L'[p] <= U'[p] */
      /* check if row lower bound is consistent */
      if (p->lb != -DBL_MAX)
      {  eps = 1e-3 + 1e-6 * fabs(p->lb);
         if (p->lb - eps > u)
         {  ret = 0x33;
            goto done;
         }
      }
      /* check if row upper bound is consistent */
      if (p->ub != +DBL_MAX)
      {  eps = 1e-3 + 1e-6 * fabs(p->ub);
         if (p->ub + eps < l)
         {  ret = 0x33;
            goto done;
         }
      }
      /* check if row lower bound can be active/forcing */
      if (p->lb != -DBL_MAX)
      {  eps = 1e-9 + 1e-12 * fabs(p->lb);
         if (p->lb - eps > l)
         {  if (p->lb + eps <= u)
               ret |= 0x01;
            else
               ret |= 0x02;
         }
      }
      /* check if row upper bound can be active/forcing */
      if (p->ub != +DBL_MAX)
      {  eps = 1e-9 + 1e-12 * fabs(p->ub);
         if (p->ub + eps < u)
         {  /* check if the upper bound is forcing */
            if (p->ub - eps >= l)
               ret |= 0x10;
            else
               ret |= 0x20;
         }
      }
done: return ret;
}

/***********************************************************************
*  NAME
*
*  npp_inactive_bound - remove row lower/upper inactive bound
*
*  SYNOPSIS
*
*  #include "glpnpp.h"
*  void npp_inactive_bound(NPP *npp, NPPROW *p, int which);
*
*  DESCRIPTION
*
*  The routine npp_inactive_bound removes lower (if which = 0) or upper
*  (if which = 1) bound of row p:
*
*     L[p] <= sum a[p,j] x[j] <= U[p],
*
*  which (bound) is assumed to be redundant.
*
*  PROBLEM TRANSFORMATION
*
*  If which = 0, current lower bound L[p] of row p is assigned -oo.
*  If which = 1, current upper bound U[p] of row p is assigned +oo.
*
*  RECOVERING BASIC SOLUTION
*
*  If in solution to the transformed problem row p is inactive
*  constraint (GLP_BS), its status is not changed in solution to the
*  original problem. Otherwise, status of row p in solution to the
*  original problem is defined by its type before transformation and
*  its status in solution to the transformed problem as follows:
*
*     +---------------------+-------+---------------+---------------+
*     |        Row          | Flag  | Row status in | Row status in |
*     |        type         | which | transfmd soln | original soln |
*     +---------------------+-------+---------------+---------------+
*     |     sum >= L[p]     |   0   |    GLP_NF     |    GLP_NL     |
*     |     sum <= U[p]     |   1   |    GLP_NF     |    GLP_NU     |
*     | L[p] <= sum <= U[p] |   0   |    GLP_NU     |    GLP_NU     |
*     | L[p] <= sum <= U[p] |   1   |    GLP_NL     |    GLP_NL     |
*     |  sum = L[p] = U[p]  |   0   |    GLP_NU     |    GLP_NS     |
*     |  sum = L[p] = U[p]  |   1   |    GLP_NL     |    GLP_NS     |
*     +---------------------+-------+---------------+---------------+
*
*  RECOVERING INTERIOR-POINT SOLUTION
*
*  None needed.
*
*  RECOVERING MIP SOLUTION
*
*  None needed. */

struct inactive_bound
{     /* row inactive bound */
      int p;
      /* row reference number */
      char stat;
      /* row status (if active constraint) */
};

static int rcv_inactive_bound(NPP *npp, void *info);

void npp_inactive_bound(NPP *npp, NPPROW *p, int which)
{     /* remove row lower/upper inactive bound */
      struct inactive_bound *info;
      if (npp->sol == GLP_SOL)
      {  /* create transformation stack entry */
         info = npp_push_tse(npp,
            rcv_inactive_bound, sizeof(struct inactive_bound));
         info->p = p->i;
         if (p->ub == +DBL_MAX)
            info->stat = GLP_NL;
         else if (p->lb == -DBL_MAX)
            info->stat = GLP_NU;
         else if (p->lb != p->ub)
            info->stat = (char)(which == 0 ? GLP_NU : GLP_NL);
         else
            info->stat = GLP_NS;
      }
      /* remove row inactive bound */
      if (which == 0)
      {  xassert(p->lb != -DBL_MAX);
         p->lb = -DBL_MAX;
      }
      else if (which == 1)
      {  xassert(p->ub != +DBL_MAX);
         p->ub = +DBL_MAX;
      }
      else
         xassert(which != which);
      return;
}

static int rcv_inactive_bound(NPP *npp, void *_info)
{     /* recover row status */
      struct inactive_bound *info = _info;
      if (npp->sol != GLP_SOL)
      {  npp_error();
         return 1;
      }
      if (npp->r_stat[info->p] == GLP_BS)
         npp->r_stat[info->p] = GLP_BS;
      else
         npp->r_stat[info->p] = info->stat;
      return 0;
}

/***********************************************************************
*  NAME
*
*  npp_implied_bounds - determine implied column bounds
*
*  SYNOPSIS
*
*  #include "glpnpp.h"
*  void npp_implied_bounds(NPP *npp, NPPROW *p);
*
*  DESCRIPTION
*
*  The routine npp_implied_bounds inspects general row (constraint) p:
*
*     L[p] <= sum a[p,j] x[j] <= U[p],                               (1)
*
*     l[j] <= x[j] <= u[j],                                          (2)
*
*  where L[p] <= U[p] and l[j] <= u[j] for all a[p,j] != 0, to compute
*  implied bounds of columns (variables x[j]) in this row.
*
*  The routine stores implied column bounds l'[j] and u'[j] in column
*  descriptors (NPPCOL); it does not change current column bounds l[j]
*  and u[j]. (Implied column bounds can be then used to strengthen the
*  current column bounds; see the routines npp_implied_lower and
*  npp_implied_upper).
*
*  ALGORITHM
*
*  Current column bounds (2) define implied lower and upper bounds of
*  row (1) as follows:
*
*     L'[p] = inf sum a[p,j] x[j] =
*                  j
*                                                                    (3)
*           =  sum  a[p,j] l[j] +  sum  a[p,j] u[j],
*            j in Jp             j in Jn
*
*     U'[p] = sup sum a[p,j] x[j] =
*                                                                    (4)
*           =  sum  a[p,j] u[j] +  sum  a[p,j] l[j],
*            j in Jp             j in Jn
*
*     Jp = {j: a[p,j] > 0},  Jn = {j: a[p,j] < 0}.                   (5)
*
*  (Note that bounds of all columns in row p are assumed to be correct,
*  so L'[p] <= U'[p].)
*
*  If L[p] > L'[p] and/or U[p] < U'[p], the lower and/or upper bound of
*  row (1) can be active, in which case such row defines implied bounds
*  of its variables.
*
*  Let x[k] be some variable having in row (1) coefficient a[p,k] != 0.
*  Consider a case when row lower bound can be active (L[p] > L'[p]):
*
*     sum a[p,j] x[j] >= L[p]  ==>
*      j
*
*     sum a[p,j] x[j] + a[p,k] x[k] >= L[p]  ==>
*     j!=k
*                                                                    (6)
*     a[p,k] x[k] >= L[p] - sum a[p,j] x[j]  ==>
*                           j!=k
*
*     a[p,k] x[k] >= L[p,k],
*
*  where
*
*     L[p,k] = inf(L[p] - sum a[p,j] x[j]) =
*                         j!=k
*
*            = L[p] - sup sum a[p,j] x[j] =                          (7)
*                         j!=k
*
*            = L[p] - sum a[p,j] u[j] - sum a[p,j] l[j].
*                    j in Jp\{k}       j in Jn\{k}
*
*  Thus:
*
*     x[k] >= l'[k] = L[p,k] / a[p,k],  if a[p,k] > 0,               (8)
*
*     x[k] <= u'[k] = L[p,k] / a[p,k],  if a[p,k] < 0.               (9)
*
*  where l'[k] and u'[k] are implied lower and upper bounds of variable
*  x[k], resp.
*
*  Now consider a similar case when row upper bound can be active
*  (U[p] < U'[p]):
*
*     sum a[p,j] x[j] <= U[p]  ==>
*      j
*
*     sum a[p,j] x[j] + a[p,k] x[k] <= U[p]  ==>
*     j!=k
*                                                                   (10)
*     a[p,k] x[k] <= U[p] - sum a[p,j] x[j]  ==>
*                           j!=k
*
*     a[p,k] x[k] <= U[p,k],
*
*  where:
*
*     U[p,k] = sup(U[p] - sum a[p,j] x[j]) =
*                         j!=k
*
*            = U[p] - inf sum a[p,j] x[j] =                         (11)
*                         j!=k
*
*            = U[p] - sum a[p,j] l[j] - sum a[p,j] u[j].
*                    j in Jp\{k}       j in Jn\{k}
*
*  Thus:
*
*     x[k] <= u'[k] = U[p,k] / a[p,k],  if a[p,k] > 0,              (12)
*
*     x[k] >= l'[k] = U[p,k] / a[p,k],  if a[p,k] < 0.              (13)
*
*  Note that in formulae (8), (9), (12), and (13) coefficient a[p,k]
*  must not be too small in magnitude relatively to other non-zero
*  coefficients in row (1), i.e. the following condition must hold:
*
*     |a[p,k]| >= eps * max(1, |a[p,j]|),                           (14)
*                        j
*
*  where eps is a relative tolerance for constraint coefficients.
*  Otherwise the implied column bounds can be numerical inreliable. For
*  example, using formula (8) for the following inequality constraint:
*
*     1e-12 x1 - x2 - x3 >= 0,
*
*  where x1 >= -1, x2, x3, >= 0, may lead to numerically unreliable
*  conclusion that x1 >= 0.
*
*  Using formulae (8), (9), (12), and (13) to compute implied bounds
*  for one variable requires |J| operations, where J = {j: a[p,j] != 0},
*  because this needs computing L[p,k] and U[p,k]. Thus, computing
*  implied bounds for all variables in row (1) would require |J|^2
*  operations, that is not a good technique. However, the total number
*  of operations can be reduced to |J| as follows.
*
*  Let a[p,k] > 0. Then from (7) and (11) we have:
*
*     L[p,k] = L[p] - (U'[p] - a[p,k] u[k]) =
*
*            = L[p] - U'[p] + a[p,k] u[k],
*
*     U[p,k] = U[p] - (L'[p] - a[p,k] l[k]) =
*
*            = U[p] - L'[p] + a[p,k] l[k],
*
*  where L'[p] and U'[p] are implied row lower and upper bounds defined
*  by formulae (3) and (4). Substituting these expressions into (8) and
*  (12) gives:
*
*     l'[k] = L[p,k] / a[p,k] = u[k] + (L[p] - U'[p]) / a[p,k],     (15)
*
*     u'[k] = U[p,k] / a[p,k] = l[k] + (U[p] - L'[p]) / a[p,k].     (16)
*
*  Similarly, if a[p,k] < 0, according to (7) and (11) we have:
*
*     L[p,k] = L[p] - (U'[p] - a[p,k] l[k]) =
*
*            = L[p] - U'[p] + a[p,k] l[k],
*
*     U[p,k] = U[p] - (L'[p] - a[p,k] u[k]) =
*
*            = U[p] - L'[p] + a[p,k] u[k],
*
*  and substituting these expressions into (8) and (12) gives:
*
*     l'[k] = U[p,k] / a[p,k] = u[k] + (U[p] - L'[p]) / a[p,k],     (17)
*
*     u'[k] = L[p,k] / a[p,k] = l[k] + (L[p] - U'[p]) / a[p,k].     (18)
*
*  Note that formulae (15)-(18) can be used only if L'[p] and U'[p]
*  exist. However, if for some variable x[j] it happens that l[j] = -oo
*  and/or u[j] = +oo, values of L'[p] (if a[p,j] > 0) and/or U'[p] (if
*  a[p,j] < 0) are undefined. Consider, therefore, the most general
*  situation, when some column bounds (2) may not exist.
*
*  Let:
*
*     J' = {j : (a[p,j] > 0 and l[j] = -oo) or
*                                                                   (19)
*               (a[p,j] < 0 and u[j] = +oo)}.
*
*  Then (assuming that row upper bound U[p] can be active) the following
*  three cases are possible:
*
*  1) |J'| = 0. In this case L'[p] exists, thus, for all variables x[j]
*     in row (1) we can use formulae (16) and (17);
*
*  2) J' = {k}. In this case L'[p] = -oo, however, U[p,k] (11) exists,
*     so for variable x[k] we can use formulae (12) and (13). Note that
*     for all other variables x[j] (j != k) l'[j] = -oo (if a[p,j] < 0)
*     or u'[j] = +oo (if a[p,j] > 0);
*
*  3) |J'| > 1. In this case for all variables x[j] in row [1] we have
*     l'[j] = -oo (if a[p,j] < 0) or u'[j] = +oo (if a[p,j] > 0).
*
*  Similarly, let:
*
*     J'' = {j : (a[p,j] > 0 and u[j] = +oo) or
*                                                                   (20)
*                (a[p,j] < 0 and l[j] = -oo)}.
*
*  Then (assuming that row lower bound L[p] can be active) the following
*  three cases are possible:
*
*  1) |J''| = 0. In this case U'[p] exists, thus, for all variables x[j]
*     in row (1) we can use formulae (15) and (18);
*
*  2) J'' = {k}. In this case U'[p] = +oo, however, L[p,k] (7) exists,
*     so for variable x[k] we can use formulae (8) and (9). Note that
*     for all other variables x[j] (j != k) l'[j] = -oo (if a[p,j] > 0)
*     or u'[j] = +oo (if a[p,j] < 0);
*
*  3) |J''| > 1. In this case for all variables x[j] in row (1) we have
*     l'[j] = -oo (if a[p,j] > 0) or u'[j] = +oo (if a[p,j] < 0). */

void npp_implied_bounds(NPP *npp, NPPROW *p)
{     NPPAIJ *apj, *apk;
      double big, eps, temp;
      xassert(npp == npp);
      /* initialize implied bounds for all variables and determine
         maximal magnitude of row coefficients a[p,j] */
      big = 1.0;
      for (apj = p->ptr; apj != NULL; apj = apj->r_next)
      {  apj->col->ll.ll = -DBL_MAX, apj->col->uu.uu = +DBL_MAX;
         if (big < fabs(apj->val)) big = fabs(apj->val);
      }
      eps = 1e-6 * big;
      /* process row lower bound (assuming that it can be active) */
      if (p->lb != -DBL_MAX)
      {  apk = NULL;
         for (apj = p->ptr; apj != NULL; apj = apj->r_next)
         {  if (apj->val > 0.0 && apj->col->ub == +DBL_MAX ||
                apj->val < 0.0 && apj->col->lb == -DBL_MAX)
            {  if (apk == NULL)
                  apk = apj;
               else
                  goto skip1;
            }
         }
         /* if a[p,k] = NULL then |J'| = 0 else J' = { k } */
         temp = p->lb;
         for (apj = p->ptr; apj != NULL; apj = apj->r_next)
         {  if (apj == apk)
               /* skip a[p,k] */;
            else if (apj->val > 0.0)
               temp -= apj->val * apj->col->ub;
            else /* apj->val < 0.0 */
               temp -= apj->val * apj->col->lb;
         }
         /* compute column implied bounds */
         if (apk == NULL)
         {  /* temp = L[p] - U'[p] */
            for (apj = p->ptr; apj != NULL; apj = apj->r_next)
            {  if (apj->val >= +eps)
               {  /* l'[j] := u[j] + (L[p] - U'[p]) / a[p,j] */
                  apj->col->ll.ll = apj->col->ub + temp / apj->val;
               }
               else if (apj->val <= -eps)
               {  /* u'[j] := l[j] + (L[p] - U'[p]) / a[p,j] */
                  apj->col->uu.uu = apj->col->lb + temp / apj->val;
               }
            }
         }
         else
         {  /* temp = L[p,k] */
            if (apk->val >= +eps)
            {  /* l'[k] := L[p,k] / a[p,k] */
               apk->col->ll.ll = temp / apk->val;
            }
            else if (apk->val <= -eps)
            {  /* u'[k] := L[p,k] / a[p,k] */
               apk->col->uu.uu = temp / apk->val;
            }
         }
skip1:   ;
      }
      /* process row upper bound (assuming that it can be active) */
      if (p->ub != +DBL_MAX)
      {  apk = NULL;
         for (apj = p->ptr; apj != NULL; apj = apj->r_next)
         {  if (apj->val > 0.0 && apj->col->lb == -DBL_MAX ||
                apj->val < 0.0 && apj->col->ub == +DBL_MAX)
            {  if (apk == NULL)
                  apk = apj;
               else
                  goto skip2;
            }
         }
         /* if a[p,k] = NULL then |J''| = 0 else J'' = { k } */
         temp = p->ub;
         for (apj = p->ptr; apj != NULL; apj = apj->r_next)
         {  if (apj == apk)
               /* skip a[p,k] */;
            else if (apj->val > 0.0)
               temp -= apj->val * apj->col->lb;
            else /* apj->val < 0.0 */
               temp -= apj->val * apj->col->ub;
         }
         /* compute column implied bounds */
         if (apk == NULL)
         {  /* temp = U[p] - L'[p] */
            for (apj = p->ptr; apj != NULL; apj = apj->r_next)
            {  if (apj->val >= +eps)
               {  /* u'[j] := l[j] + (U[p] - L'[p]) / a[p,j] */
                  apj->col->uu.uu = apj->col->lb + temp / apj->val;
               }
               else if (apj->val <= -eps)
               {  /* l'[j] := u[j] + (U[p] - L'[p]) / a[p,j] */
                  apj->col->ll.ll = apj->col->ub + temp / apj->val;
               }
            }
         }
         else
         {  /* temp = U[p,k] */
            if (apk->val >= +eps)
            {  /* u'[k] := U[p,k] / a[p,k] */
               apk->col->uu.uu = temp / apk->val;
            }
            else if (apk->val <= -eps)
            {  /* l'[k] := U[p,k] / a[p,k] */
               apk->col->ll.ll = temp / apk->val;
            }
         }
skip2:   ;
      }
      return;
}

/* eof */
