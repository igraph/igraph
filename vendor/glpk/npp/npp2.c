/* npp2.c */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2009-2017 Free Software Foundation, Inc.
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
#include "npp.h"

/***********************************************************************
*  NAME
*
*  npp_free_row - process free (unbounded) row
*
*  SYNOPSIS
*
*  #include "glpnpp.h"
*  void npp_free_row(NPP *npp, NPPROW *p);
*
*  DESCRIPTION
*
*  The routine npp_free_row processes row p, which is free (i.e. has
*  no finite bounds):
*
*     -inf < sum a[p,j] x[j] < +inf.                                 (1)
*             j
*
*  PROBLEM TRANSFORMATION
*
*  Constraint (1) cannot be active, so it is redundant and can be
*  removed from the original problem.
*
*  Removing row p leads to removing a column of multiplier pi[p] for
*  this row in the dual system. Since row p has no bounds, pi[p] = 0,
*  so removing the column does not affect the dual solution.
*
*  RECOVERING BASIC SOLUTION
*
*  In solution to the original problem row p is inactive constraint,
*  so it is assigned status GLP_BS, and multiplier pi[p] is assigned
*  zero value.
*
*  RECOVERING INTERIOR-POINT SOLUTION
*
*  In solution to the original problem row p is inactive constraint,
*  so its multiplier pi[p] is assigned zero value.
*
*  RECOVERING MIP SOLUTION
*
*  None needed. */

struct free_row
{     /* free (unbounded) row */
      int p;
      /* row reference number */
};

static int rcv_free_row(NPP *npp, void *info);

void npp_free_row(NPP *npp, NPPROW *p)
{     /* process free (unbounded) row */
      struct free_row *info;
      /* the row must be free */
      xassert(p->lb == -DBL_MAX && p->ub == +DBL_MAX);
      /* create transformation stack entry */
      info = npp_push_tse(npp,
         rcv_free_row, sizeof(struct free_row));
      info->p = p->i;
      /* remove the row from the problem */
      npp_del_row(npp, p);
      return;
}

static int rcv_free_row(NPP *npp, void *_info)
{     /* recover free (unbounded) row */
      struct free_row *info = _info;
      if (npp->sol == GLP_SOL)
         npp->r_stat[info->p] = GLP_BS;
      if (npp->sol != GLP_MIP)
         npp->r_pi[info->p] = 0.0;
      return 0;
}

/***********************************************************************
*  NAME
*
*  npp_geq_row - process row of 'not less than' type
*
*  SYNOPSIS
*
*  #include "glpnpp.h"
*  void npp_geq_row(NPP *npp, NPPROW *p);
*
*  DESCRIPTION
*
*  The routine npp_geq_row processes row p, which is 'not less than'
*  inequality constraint:
*
*     L[p] <= sum a[p,j] x[j] (<= U[p]),                             (1)
*              j
*
*  where L[p] < U[p], and upper bound may not exist (U[p] = +oo).
*
*  PROBLEM TRANSFORMATION
*
*  Constraint (1) can be replaced by equality constraint:
*
*     sum a[p,j] x[j] - s = L[p],                                    (2)
*      j
*
*  where
*
*     0 <= s (<= U[p] - L[p])                                        (3)
*
*  is a non-negative surplus variable.
*
*  Since in the primal system there appears column s having the only
*  non-zero coefficient in row p, in the dual system there appears a
*  new row:
*
*     (-1) pi[p] + lambda = 0,                                       (4)
*
*  where (-1) is coefficient of column s in row p, pi[p] is multiplier
*  of row p, lambda is multiplier of column q, 0 is coefficient of
*  column s in the objective row.
*
*  RECOVERING BASIC SOLUTION
*
*  Status of row p in solution to the original problem is determined
*  by its status and status of column q in solution to the transformed
*  problem as follows:
*
*     +--------------------------------------+------------------+
*     |         Transformed problem          | Original problem |
*     +-----------------+--------------------+------------------+
*     | Status of row p | Status of column s | Status of row p  |
*     +-----------------+--------------------+------------------+
*     |     GLP_BS      |       GLP_BS       |       N/A        |
*     |     GLP_BS      |       GLP_NL       |      GLP_BS      |
*     |     GLP_BS      |       GLP_NU       |      GLP_BS      |
*     |     GLP_NS      |       GLP_BS       |      GLP_BS      |
*     |     GLP_NS      |       GLP_NL       |      GLP_NL      |
*     |     GLP_NS      |       GLP_NU       |      GLP_NU      |
*     +-----------------+--------------------+------------------+
*
*  Value of row multiplier pi[p] in solution to the original problem
*  is the same as in solution to the transformed problem.
*
*  1. In solution to the transformed problem row p and column q cannot
*     be basic at the same time; otherwise the basis matrix would have
*     two linear dependent columns: unity column of auxiliary variable
*     of row p and unity column of variable s.
*
*  2. Though in the transformed problem row p is equality constraint,
*     it may be basic due to primal degenerate solution.
*
*  RECOVERING INTERIOR-POINT SOLUTION
*
*  Value of row multiplier pi[p] in solution to the original problem
*  is the same as in solution to the transformed problem.
*
*  RECOVERING MIP SOLUTION
*
*  None needed. */

struct ineq_row
{     /* inequality constraint row */
      int p;
      /* row reference number */
      int s;
      /* column reference number for slack/surplus variable */
};

static int rcv_geq_row(NPP *npp, void *info);

void npp_geq_row(NPP *npp, NPPROW *p)
{     /* process row of 'not less than' type */
      struct ineq_row *info;
      NPPCOL *s;
      /* the row must have lower bound */
      xassert(p->lb != -DBL_MAX);
      xassert(p->lb < p->ub);
      /* create column for surplus variable */
      s = npp_add_col(npp);
      s->lb = 0.0;
      s->ub = (p->ub == +DBL_MAX ? +DBL_MAX : p->ub - p->lb);
      /* and add it to the transformed problem */
      npp_add_aij(npp, p, s, -1.0);
      /* create transformation stack entry */
      info = npp_push_tse(npp,
         rcv_geq_row, sizeof(struct ineq_row));
      info->p = p->i;
      info->s = s->j;
      /* replace the row by equality constraint */
      p->ub = p->lb;
      return;
}

static int rcv_geq_row(NPP *npp, void *_info)
{     /* recover row of 'not less than' type */
      struct ineq_row *info = _info;
      if (npp->sol == GLP_SOL)
      {  if (npp->r_stat[info->p] == GLP_BS)
         {  if (npp->c_stat[info->s] == GLP_BS)
            {  npp_error();
               return 1;
            }
            else if (npp->c_stat[info->s] == GLP_NL ||
                     npp->c_stat[info->s] == GLP_NU)
               npp->r_stat[info->p] = GLP_BS;
            else
            {  npp_error();
               return 1;
            }
         }
         else if (npp->r_stat[info->p] == GLP_NS)
         {  if (npp->c_stat[info->s] == GLP_BS)
               npp->r_stat[info->p] = GLP_BS;
            else if (npp->c_stat[info->s] == GLP_NL)
               npp->r_stat[info->p] = GLP_NL;
            else if (npp->c_stat[info->s] == GLP_NU)
               npp->r_stat[info->p] = GLP_NU;
            else
            {  npp_error();
               return 1;
            }
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
*  npp_leq_row - process row of 'not greater than' type
*
*  SYNOPSIS
*
*  #include "glpnpp.h"
*  void npp_leq_row(NPP *npp, NPPROW *p);
*
*  DESCRIPTION
*
*  The routine npp_leq_row processes row p, which is 'not greater than'
*  inequality constraint:
*
*     (L[p] <=) sum a[p,j] x[j] <= U[p],                             (1)
*                j
*
*  where L[p] < U[p], and lower bound may not exist (L[p] = +oo).
*
*  PROBLEM TRANSFORMATION
*
*  Constraint (1) can be replaced by equality constraint:
*
*     sum a[p,j] x[j] + s = L[p],                                    (2)
*      j
*
*  where
*
*     0 <= s (<= U[p] - L[p])                                        (3)
*
*  is a non-negative slack variable.
*
*  Since in the primal system there appears column s having the only
*  non-zero coefficient in row p, in the dual system there appears a
*  new row:
*
*     (+1) pi[p] + lambda = 0,                                       (4)
*
*  where (+1) is coefficient of column s in row p, pi[p] is multiplier
*  of row p, lambda is multiplier of column q, 0 is coefficient of
*  column s in the objective row.
*
*  RECOVERING BASIC SOLUTION
*
*  Status of row p in solution to the original problem is determined
*  by its status and status of column q in solution to the transformed
*  problem as follows:
*
*     +--------------------------------------+------------------+
*     |         Transformed problem          | Original problem |
*     +-----------------+--------------------+------------------+
*     | Status of row p | Status of column s | Status of row p  |
*     +-----------------+--------------------+------------------+
*     |     GLP_BS      |       GLP_BS       |       N/A        |
*     |     GLP_BS      |       GLP_NL       |      GLP_BS      |
*     |     GLP_BS      |       GLP_NU       |      GLP_BS      |
*     |     GLP_NS      |       GLP_BS       |      GLP_BS      |
*     |     GLP_NS      |       GLP_NL       |      GLP_NU      |
*     |     GLP_NS      |       GLP_NU       |      GLP_NL      |
*     +-----------------+--------------------+------------------+
*
*  Value of row multiplier pi[p] in solution to the original problem
*  is the same as in solution to the transformed problem.
*
*  1. In solution to the transformed problem row p and column q cannot
*     be basic at the same time; otherwise the basis matrix would have
*     two linear dependent columns: unity column of auxiliary variable
*     of row p and unity column of variable s.
*
*  2. Though in the transformed problem row p is equality constraint,
*     it may be basic due to primal degeneracy.
*
*  RECOVERING INTERIOR-POINT SOLUTION
*
*  Value of row multiplier pi[p] in solution to the original problem
*  is the same as in solution to the transformed problem.
*
*  RECOVERING MIP SOLUTION
*
*  None needed. */

static int rcv_leq_row(NPP *npp, void *info);

void npp_leq_row(NPP *npp, NPPROW *p)
{     /* process row of 'not greater than' type */
      struct ineq_row *info;
      NPPCOL *s;
      /* the row must have upper bound */
      xassert(p->ub != +DBL_MAX);
      xassert(p->lb < p->ub);
      /* create column for slack variable */
      s = npp_add_col(npp);
      s->lb = 0.0;
      s->ub = (p->lb == -DBL_MAX ? +DBL_MAX : p->ub - p->lb);
      /* and add it to the transformed problem */
      npp_add_aij(npp, p, s, +1.0);
      /* create transformation stack entry */
      info = npp_push_tse(npp,
         rcv_leq_row, sizeof(struct ineq_row));
      info->p = p->i;
      info->s = s->j;
      /* replace the row by equality constraint */
      p->lb = p->ub;
      return;
}

static int rcv_leq_row(NPP *npp, void *_info)
{     /* recover row of 'not greater than' type */
      struct ineq_row *info = _info;
      if (npp->sol == GLP_SOL)
      {  if (npp->r_stat[info->p] == GLP_BS)
         {  if (npp->c_stat[info->s] == GLP_BS)
            {  npp_error();
               return 1;
            }
            else if (npp->c_stat[info->s] == GLP_NL ||
                     npp->c_stat[info->s] == GLP_NU)
               npp->r_stat[info->p] = GLP_BS;
            else
            {  npp_error();
               return 1;
            }
         }
         else if (npp->r_stat[info->p] == GLP_NS)
         {  if (npp->c_stat[info->s] == GLP_BS)
               npp->r_stat[info->p] = GLP_BS;
            else if (npp->c_stat[info->s] == GLP_NL)
               npp->r_stat[info->p] = GLP_NU;
            else if (npp->c_stat[info->s] == GLP_NU)
               npp->r_stat[info->p] = GLP_NL;
            else
            {  npp_error();
               return 1;
            }
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
*  npp_free_col - process free (unbounded) column
*
*  SYNOPSIS
*
*  #include "glpnpp.h"
*  void npp_free_col(NPP *npp, NPPCOL *q);
*
*  DESCRIPTION
*
*  The routine npp_free_col processes column q, which is free (i.e. has
*  no finite bounds):
*
*     -oo < x[q] < +oo.                                              (1)
*
*  PROBLEM TRANSFORMATION
*
*  Free (unbounded) variable can be replaced by the difference of two
*  non-negative variables:
*
*     x[q] = s' - s'',   s', s'' >= 0.                               (2)
*
*  Assuming that in the transformed problem x[q] becomes s',
*  transformation (2) causes new column s'' to appear, which differs
*  from column s' only in the sign of coefficients in constraint and
*  objective rows. Thus, if in the dual system the following row
*  corresponds to column s':
*
*     sum a[i,q] pi[i] + lambda' = c[q],                             (3)
*      i
*
*  the row which corresponds to column s'' is the following:
*
*     sum (-a[i,q]) pi[i] + lambda'' = -c[q].                        (4)
*      i
*
*  Then from (3) and (4) it follows that:
*
*     lambda' + lambda'' = 0   =>   lambda' = lmabda'' = 0,          (5)
*
*  where lambda' and lambda'' are multipliers for columns s' and s'',
*  resp.
*
*  RECOVERING BASIC SOLUTION
*
*  With respect to (5) status of column q in solution to the original
*  problem is determined by statuses of columns s' and s'' in solution
*  to the transformed problem as follows:
*
*     +--------------------------------------+------------------+
*     |         Transformed problem          | Original problem |
*     +------------------+-------------------+------------------+
*     | Status of col s' | Status of col s'' | Status of col q  |
*     +------------------+-------------------+------------------+
*     |      GLP_BS      |      GLP_BS       |       N/A        |
*     |      GLP_BS      |      GLP_NL       |      GLP_BS      |
*     |      GLP_NL      |      GLP_BS       |      GLP_BS      |
*     |      GLP_NL      |      GLP_NL       |      GLP_NF      |
*     +------------------+-------------------+------------------+
*
*  Value of column q is computed with formula (2).
*
*  1. In solution to the transformed problem columns s' and s'' cannot
*     be basic at the same time, because they differ only in the sign,
*     hence, are linear dependent.
*
*  2. Though column q is free, it can be non-basic due to dual
*     degeneracy.
*
*  3. If column q is integral, columns s' and s'' are also integral.
*
*  RECOVERING INTERIOR-POINT SOLUTION
*
*  Value of column q is computed with formula (2).
*
*  RECOVERING MIP SOLUTION
*
*  Value of column q is computed with formula (2). */

struct free_col
{     /* free (unbounded) column */
      int q;
      /* column reference number for variables x[q] and s' */
      int s;
      /* column reference number for variable s'' */
};

static int rcv_free_col(NPP *npp, void *info);

void npp_free_col(NPP *npp, NPPCOL *q)
{     /* process free (unbounded) column */
      struct free_col *info;
      NPPCOL *s;
      NPPAIJ *aij;
      /* the column must be free */
      xassert(q->lb == -DBL_MAX && q->ub == +DBL_MAX);
      /* variable x[q] becomes s' */
      q->lb = 0.0, q->ub = +DBL_MAX;
      /* create variable s'' */
      s = npp_add_col(npp);
      s->is_int = q->is_int;
      s->lb = 0.0, s->ub = +DBL_MAX;
      /* duplicate objective coefficient */
      s->coef = -q->coef;
      /* duplicate column of the constraint matrix */
      for (aij = q->ptr; aij != NULL; aij = aij->c_next)
         npp_add_aij(npp, aij->row, s, -aij->val);
      /* create transformation stack entry */
      info = npp_push_tse(npp,
         rcv_free_col, sizeof(struct free_col));
      info->q = q->j;
      info->s = s->j;
      return;
}

static int rcv_free_col(NPP *npp, void *_info)
{     /* recover free (unbounded) column */
      struct free_col *info = _info;
      if (npp->sol == GLP_SOL)
      {  if (npp->c_stat[info->q] == GLP_BS)
         {  if (npp->c_stat[info->s] == GLP_BS)
            {  npp_error();
               return 1;
            }
            else if (npp->c_stat[info->s] == GLP_NL)
               npp->c_stat[info->q] = GLP_BS;
            else
            {  npp_error();
               return -1;
            }
         }
         else if (npp->c_stat[info->q] == GLP_NL)
         {  if (npp->c_stat[info->s] == GLP_BS)
               npp->c_stat[info->q] = GLP_BS;
            else if (npp->c_stat[info->s] == GLP_NL)
               npp->c_stat[info->q] = GLP_NF;
            else
            {  npp_error();
               return -1;
            }
         }
         else
         {  npp_error();
            return -1;
         }
      }
      /* compute value of x[q] with formula (2) */
      npp->c_value[info->q] -= npp->c_value[info->s];
      return 0;
}

/***********************************************************************
*  NAME
*
*  npp_lbnd_col - process column with (non-zero) lower bound
*
*  SYNOPSIS
*
*  #include "glpnpp.h"
*  void npp_lbnd_col(NPP *npp, NPPCOL *q);
*
*  DESCRIPTION
*
*  The routine npp_lbnd_col processes column q, which has (non-zero)
*  lower bound:
*
*     l[q] <= x[q] (<= u[q]),                                        (1)
*
*  where l[q] < u[q], and upper bound may not exist (u[q] = +oo).
*
*  PROBLEM TRANSFORMATION
*
*  Column q can be replaced as follows:
*
*     x[q] = l[q] + s,                                               (2)
*
*  where
*
*     0 <= s (<= u[q] - l[q])                                        (3)
*
*  is a non-negative variable.
*
*  Substituting x[q] from (2) into the objective row, we have:
*
*     z = sum c[j] x[j] + c0 =
*          j
*
*       = sum c[j] x[j] + c[q] x[q] + c0 =
*         j!=q
*
*       = sum c[j] x[j] + c[q] (l[q] + s) + c0 =
*         j!=q
*
*       = sum c[j] x[j] + c[q] s + c~0,
*
*  where
*
*     c~0 = c0 + c[q] l[q]                                           (4)
*
*  is the constant term of the objective in the transformed problem.
*  Similarly, substituting x[q] into constraint row i, we have:
*
*     L[i] <= sum a[i,j] x[j] <= U[i]  ==>
*              j
*
*     L[i] <= sum a[i,j] x[j] + a[i,q] x[q] <= U[i]  ==>
*             j!=q
*
*     L[i] <= sum a[i,j] x[j] + a[i,q] (l[q] + s) <= U[i]  ==>
*             j!=q
*
*     L~[i] <= sum a[i,j] x[j] + a[i,q] s <= U~[i],
*              j!=q
*
*  where
*
*     L~[i] = L[i] - a[i,q] l[q],  U~[i] = U[i] - a[i,q] l[q]        (5)
*
*  are lower and upper bounds of row i in the transformed problem,
*  resp.
*
*  Transformation (2) does not affect the dual system.
*
*  RECOVERING BASIC SOLUTION
*
*  Status of column q in solution to the original problem is the same
*  as in solution to the transformed problem (GLP_BS, GLP_NL or GLP_NU).
*  Value of column q is computed with formula (2).
*
*  RECOVERING INTERIOR-POINT SOLUTION
*
*  Value of column q is computed with formula (2).
*
*  RECOVERING MIP SOLUTION
*
*  Value of column q is computed with formula (2). */

struct bnd_col
{     /* bounded column */
      int q;
      /* column reference number for variables x[q] and s */
      double bnd;
      /* lower/upper bound l[q] or u[q] */
};

static int rcv_lbnd_col(NPP *npp, void *info);

void npp_lbnd_col(NPP *npp, NPPCOL *q)
{     /* process column with (non-zero) lower bound */
      struct bnd_col *info;
      NPPROW *i;
      NPPAIJ *aij;
      /* the column must have non-zero lower bound */
      xassert(q->lb != 0.0);
      xassert(q->lb != -DBL_MAX);
      xassert(q->lb < q->ub);
      /* create transformation stack entry */
      info = npp_push_tse(npp,
         rcv_lbnd_col, sizeof(struct bnd_col));
      info->q = q->j;
      info->bnd = q->lb;
      /* substitute x[q] into objective row */
      npp->c0 += q->coef * q->lb;
      /* substitute x[q] into constraint rows */
      for (aij = q->ptr; aij != NULL; aij = aij->c_next)
      {  i = aij->row;
         if (i->lb == i->ub)
            i->ub = (i->lb -= aij->val * q->lb);
         else
         {  if (i->lb != -DBL_MAX)
               i->lb -= aij->val * q->lb;
            if (i->ub != +DBL_MAX)
               i->ub -= aij->val * q->lb;
         }
      }
      /* column x[q] becomes column s */
      if (q->ub != +DBL_MAX)
         q->ub -= q->lb;
      q->lb = 0.0;
      return;
}

static int rcv_lbnd_col(NPP *npp, void *_info)
{     /* recover column with (non-zero) lower bound */
      struct bnd_col *info = _info;
      if (npp->sol == GLP_SOL)
      {  if (npp->c_stat[info->q] == GLP_BS ||
             npp->c_stat[info->q] == GLP_NL ||
             npp->c_stat[info->q] == GLP_NU)
            npp->c_stat[info->q] = npp->c_stat[info->q];
         else
         {  npp_error();
            return 1;
         }
      }
      /* compute value of x[q] with formula (2) */
      npp->c_value[info->q] = info->bnd + npp->c_value[info->q];
      return 0;
}

/***********************************************************************
*  NAME
*
*  npp_ubnd_col - process column with upper bound
*
*  SYNOPSIS
*
*  #include "glpnpp.h"
*  void npp_ubnd_col(NPP *npp, NPPCOL *q);
*
*  DESCRIPTION
*
*  The routine npp_ubnd_col processes column q, which has upper bound:
*
*     (l[q] <=) x[q] <= u[q],                                        (1)
*
*  where l[q] < u[q], and lower bound may not exist (l[q] = -oo).
*
*  PROBLEM TRANSFORMATION
*
*  Column q can be replaced as follows:
*
*     x[q] = u[q] - s,                                               (2)
*
*  where
*
*     0 <= s (<= u[q] - l[q])                                        (3)
*
*  is a non-negative variable.
*
*  Substituting x[q] from (2) into the objective row, we have:
*
*     z = sum c[j] x[j] + c0 =
*          j
*
*       = sum c[j] x[j] + c[q] x[q] + c0 =
*         j!=q
*
*       = sum c[j] x[j] + c[q] (u[q] - s) + c0 =
*         j!=q
*
*       = sum c[j] x[j] - c[q] s + c~0,
*
*  where
*
*     c~0 = c0 + c[q] u[q]                                           (4)
*
*  is the constant term of the objective in the transformed problem.
*  Similarly, substituting x[q] into constraint row i, we have:
*
*     L[i] <= sum a[i,j] x[j] <= U[i]  ==>
*              j
*
*     L[i] <= sum a[i,j] x[j] + a[i,q] x[q] <= U[i]  ==>
*             j!=q
*
*     L[i] <= sum a[i,j] x[j] + a[i,q] (u[q] - s) <= U[i]  ==>
*             j!=q
*
*     L~[i] <= sum a[i,j] x[j] - a[i,q] s <= U~[i],
*              j!=q
*
*  where
*
*     L~[i] = L[i] - a[i,q] u[q],  U~[i] = U[i] - a[i,q] u[q]        (5)
*
*  are lower and upper bounds of row i in the transformed problem,
*  resp.
*
*  Note that in the transformed problem coefficients c[q] and a[i,q]
*  change their sign. Thus, the row of the dual system corresponding to
*  column q:
*
*     sum a[i,q] pi[i] + lambda[q] = c[q]                            (6)
*      i
*
*  in the transformed problem becomes the following:
*
*     sum (-a[i,q]) pi[i] + lambda[s] = -c[q].                       (7)
*      i
*
*  Therefore:
*
*     lambda[q] = - lambda[s],                                       (8)
*
*  where lambda[q] is multiplier for column q, lambda[s] is multiplier
*  for column s.
*
*  RECOVERING BASIC SOLUTION
*
*  With respect to (8) status of column q in solution to the original
*  problem is determined by status of column s in solution to the
*  transformed problem as follows:
*
*     +-----------------------+--------------------+
*     |  Status of column s   | Status of column q |
*     | (transformed problem) | (original problem) |
*     +-----------------------+--------------------+
*     |        GLP_BS         |       GLP_BS       |
*     |        GLP_NL         |       GLP_NU       |
*     |        GLP_NU         |       GLP_NL       |
*     +-----------------------+--------------------+
*
*  Value of column q is computed with formula (2).
*
*  RECOVERING INTERIOR-POINT SOLUTION
*
*  Value of column q is computed with formula (2).
*
*  RECOVERING MIP SOLUTION
*
*  Value of column q is computed with formula (2). */

static int rcv_ubnd_col(NPP *npp, void *info);

void npp_ubnd_col(NPP *npp, NPPCOL *q)
{     /* process column with upper bound */
      struct bnd_col *info;
      NPPROW *i;
      NPPAIJ *aij;
      /* the column must have upper bound */
      xassert(q->ub != +DBL_MAX);
      xassert(q->lb < q->ub);
      /* create transformation stack entry */
      info = npp_push_tse(npp,
         rcv_ubnd_col, sizeof(struct bnd_col));
      info->q = q->j;
      info->bnd = q->ub;
      /* substitute x[q] into objective row */
      npp->c0 += q->coef * q->ub;
      q->coef = -q->coef;
      /* substitute x[q] into constraint rows */
      for (aij = q->ptr; aij != NULL; aij = aij->c_next)
      {  i = aij->row;
         if (i->lb == i->ub)
            i->ub = (i->lb -= aij->val * q->ub);
         else
         {  if (i->lb != -DBL_MAX)
               i->lb -= aij->val * q->ub;
            if (i->ub != +DBL_MAX)
               i->ub -= aij->val * q->ub;
         }
         aij->val = -aij->val;
      }
      /* column x[q] becomes column s */
      if (q->lb != -DBL_MAX)
         q->ub -= q->lb;
      else
         q->ub = +DBL_MAX;
      q->lb = 0.0;
      return;
}

static int rcv_ubnd_col(NPP *npp, void *_info)
{     /* recover column with upper bound */
      struct bnd_col *info = _info;
      if (npp->sol == GLP_BS)
      {  if (npp->c_stat[info->q] == GLP_BS)
            npp->c_stat[info->q] = GLP_BS;
         else if (npp->c_stat[info->q] == GLP_NL)
            npp->c_stat[info->q] = GLP_NU;
         else if (npp->c_stat[info->q] == GLP_NU)
            npp->c_stat[info->q] = GLP_NL;
         else
         {  npp_error();
            return 1;
         }
      }
      /* compute value of x[q] with formula (2) */
      npp->c_value[info->q] = info->bnd - npp->c_value[info->q];
      return 0;
}

/***********************************************************************
*  NAME
*
*  npp_dbnd_col - process non-negative column with upper bound
*
*  SYNOPSIS
*
*  #include "glpnpp.h"
*  void npp_dbnd_col(NPP *npp, NPPCOL *q);
*
*  DESCRIPTION
*
*  The routine npp_dbnd_col processes column q, which is non-negative
*  and has upper bound:
*
*     0 <= x[q] <= u[q],                                             (1)
*
*  where u[q] > 0.
*
*  PROBLEM TRANSFORMATION
*
*  Upper bound of column q can be replaced by the following equality
*  constraint:
*
*     x[q] + s = u[q],                                               (2)
*
*  where s >= 0 is a non-negative complement variable.
*
*  Since in the primal system along with new row (2) there appears a
*  new column s having the only non-zero coefficient in this row, in
*  the dual system there appears a new row:
*
*     (+1)pi + lambda[s] = 0,                                        (3)
*
*  where (+1) is coefficient at column s in row (2), pi is multiplier
*  for row (2), lambda[s] is multiplier for column s, 0 is coefficient
*  at column s in the objective row.
*
*  RECOVERING BASIC SOLUTION
*
*  Status of column q in solution to the original problem is determined
*  by its status and status of column s in solution to the transformed
*  problem as follows:
*
*     +-----------------------------------+------------------+
*     |         Transformed problem       | Original problem |
*     +-----------------+-----------------+------------------+
*     | Status of col q | Status of col s | Status of col q  |
*     +-----------------+-----------------+------------------+
*     |     GLP_BS      |     GLP_BS      |      GLP_BS      |
*     |     GLP_BS      |     GLP_NL      |      GLP_NU      |
*     |     GLP_NL      |     GLP_BS      |      GLP_NL      |
*     |     GLP_NL      |     GLP_NL      |      GLP_NL (*)  |
*     +-----------------+-----------------+------------------+
*
*  Value of column q in solution to the original problem is the same as
*  in solution to the transformed problem.
*
*  1. Formally, in solution to the transformed problem columns q and s
*     cannot be non-basic at the same time, since the constraint (2)
*     would be violated. However, if u[q] is close to zero, violation
*     may be less than a working precision even if both columns q and s
*     are non-basic. In this degenerate case row (2) can be only basic,
*     i.e. non-active constraint (otherwise corresponding row of the
*     basis matrix would be zero). This allows to pivot out auxiliary
*     variable and pivot in column s, in which case the row becomes
*     active while column s becomes basic.
*
*  2. If column q is integral, column s is also integral.
*
*  RECOVERING INTERIOR-POINT SOLUTION
*
*  Value of column q in solution to the original problem is the same as
*  in solution to the transformed problem.
*
*  RECOVERING MIP SOLUTION
*
*  Value of column q in solution to the original problem is the same as
*  in solution to the transformed problem. */

struct dbnd_col
{     /* double-bounded column */
      int q;
      /* column reference number for variable x[q] */
      int s;
      /* column reference number for complement variable s */
};

static int rcv_dbnd_col(NPP *npp, void *info);

void npp_dbnd_col(NPP *npp, NPPCOL *q)
{     /* process non-negative column with upper bound */
      struct dbnd_col *info;
      NPPROW *p;
      NPPCOL *s;
      /* the column must be non-negative with upper bound */
      xassert(q->lb == 0.0);
      xassert(q->ub > 0.0);
      xassert(q->ub != +DBL_MAX);
      /* create variable s */
      s = npp_add_col(npp);
      s->is_int = q->is_int;
      s->lb = 0.0, s->ub = +DBL_MAX;
      /* create equality constraint (2) */
      p = npp_add_row(npp);
      p->lb = p->ub = q->ub;
      npp_add_aij(npp, p, q, +1.0);
      npp_add_aij(npp, p, s, +1.0);
      /* create transformation stack entry */
      info = npp_push_tse(npp,
         rcv_dbnd_col, sizeof(struct dbnd_col));
      info->q = q->j;
      info->s = s->j;
      /* remove upper bound of x[q] */
      q->ub = +DBL_MAX;
      return;
}

static int rcv_dbnd_col(NPP *npp, void *_info)
{     /* recover non-negative column with upper bound */
      struct dbnd_col *info = _info;
      if (npp->sol == GLP_BS)
      {  if (npp->c_stat[info->q] == GLP_BS)
         {  if (npp->c_stat[info->s] == GLP_BS)
               npp->c_stat[info->q] = GLP_BS;
            else if (npp->c_stat[info->s] == GLP_NL)
               npp->c_stat[info->q] = GLP_NU;
            else
            {  npp_error();
               return 1;
            }
         }
         else if (npp->c_stat[info->q] == GLP_NL)
         {  if (npp->c_stat[info->s] == GLP_BS ||
                npp->c_stat[info->s] == GLP_NL)
               npp->c_stat[info->q] = GLP_NL;
            else
            {  npp_error();
               return 1;
            }
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
*  npp_fixed_col - process fixed column
*
*  SYNOPSIS
*
*  #include "glpnpp.h"
*  void npp_fixed_col(NPP *npp, NPPCOL *q);
*
*  DESCRIPTION
*
*  The routine npp_fixed_col processes column q, which is fixed:
*
*     x[q] = s[q],                                                   (1)
*
*  where s[q] is a fixed column value.
*
*  PROBLEM TRANSFORMATION
*
*  The value of a fixed column can be substituted into the objective
*  and constraint rows that allows removing the column from the problem.
*
*  Substituting x[q] = s[q] into the objective row, we have:
*
*     z = sum c[j] x[j] + c0 =
*          j
*
*       = sum c[j] x[j] + c[q] x[q] + c0 =
*         j!=q
*
*       = sum c[j] x[j] + c[q] s[q] + c0 =
*         j!=q
*
*       = sum c[j] x[j] + c~0,
*         j!=q
*
*  where
*
*     c~0 = c0 + c[q] s[q]                                           (2)
*
*  is the constant term of the objective in the transformed problem.
*  Similarly, substituting x[q] = s[q] into constraint row i, we have:
*
*     L[i] <= sum a[i,j] x[j] <= U[i]  ==>
*              j
*
*     L[i] <= sum a[i,j] x[j] + a[i,q] x[q] <= U[i]  ==>
*             j!=q
*
*     L[i] <= sum a[i,j] x[j] + a[i,q] s[q] <= U[i]  ==>
*             j!=q
*
*     L~[i] <= sum a[i,j] x[j] + a[i,q] s <= U~[i],
*              j!=q
*
*  where
*
*     L~[i] = L[i] - a[i,q] s[q],  U~[i] = U[i] - a[i,q] s[q]        (3)
*
*  are lower and upper bounds of row i in the transformed problem,
*  resp.
*
*  RECOVERING BASIC SOLUTION
*
*  Column q is assigned status GLP_NS and its value is assigned s[q].
*
*  RECOVERING INTERIOR-POINT SOLUTION
*
*  Value of column q is assigned s[q].
*
*  RECOVERING MIP SOLUTION
*
*  Value of column q is assigned s[q]. */

struct fixed_col
{     /* fixed column */
      int q;
      /* column reference number for variable x[q] */
      double s;
      /* value, at which x[q] is fixed */
};

static int rcv_fixed_col(NPP *npp, void *info);

void npp_fixed_col(NPP *npp, NPPCOL *q)
{     /* process fixed column */
      struct fixed_col *info;
      NPPROW *i;
      NPPAIJ *aij;
      /* the column must be fixed */
      xassert(q->lb == q->ub);
      /* create transformation stack entry */
      info = npp_push_tse(npp,
         rcv_fixed_col, sizeof(struct fixed_col));
      info->q = q->j;
      info->s = q->lb;
      /* substitute x[q] = s[q] into objective row */
      npp->c0 += q->coef * q->lb;
      /* substitute x[q] = s[q] into constraint rows */
      for (aij = q->ptr; aij != NULL; aij = aij->c_next)
      {  i = aij->row;
         if (i->lb == i->ub)
            i->ub = (i->lb -= aij->val * q->lb);
         else
         {  if (i->lb != -DBL_MAX)
               i->lb -= aij->val * q->lb;
            if (i->ub != +DBL_MAX)
               i->ub -= aij->val * q->lb;
         }
      }
      /* remove the column from the problem */
      npp_del_col(npp, q);
      return;
}

static int rcv_fixed_col(NPP *npp, void *_info)
{     /* recover fixed column */
      struct fixed_col *info = _info;
      if (npp->sol == GLP_SOL)
         npp->c_stat[info->q] = GLP_NS;
      npp->c_value[info->q] = info->s;
      return 0;
}

/***********************************************************************
*  NAME
*
*  npp_make_equality - process row with almost identical bounds
*
*  SYNOPSIS
*
*  #include "glpnpp.h"
*  int npp_make_equality(NPP *npp, NPPROW *p);
*
*  DESCRIPTION
*
*  The routine npp_make_equality processes row p:
*
*     L[p] <= sum a[p,j] x[j] <= U[p],                               (1)
*              j
*
*  where -oo < L[p] < U[p] < +oo, i.e. which is double-sided inequality
*  constraint.
*
*  RETURNS
*
*  0 - row bounds have not been changed;
*
*  1 - row has been replaced by equality constraint.
*
*  PROBLEM TRANSFORMATION
*
*  If bounds of row (1) are very close to each other:
*
*     U[p] - L[p] <= eps,                                            (2)
*
*  where eps is an absolute tolerance for row value, the row can be
*  replaced by the following almost equivalent equiality constraint:
*
*     sum a[p,j] x[j] = b,                                           (3)
*      j
*
*  where b = (L[p] + U[p]) / 2. If the right-hand side in (3) happens
*  to be very close to its nearest integer:
*
*     |b - floor(b + 0.5)| <= eps,                                   (4)
*
*  it is reasonable to use this nearest integer as the right-hand side.
*
*  RECOVERING BASIC SOLUTION
*
*  Status of row p in solution to the original problem is determined
*  by its status and the sign of its multiplier pi[p] in solution to
*  the transformed problem as follows:
*
*     +-----------------------+---------+--------------------+
*     |    Status of row p    | Sign of |  Status of row p   |
*     | (transformed problem) |  pi[p]  | (original problem) |
*     +-----------------------+---------+--------------------+
*     |        GLP_BS         |  + / -  |       GLP_BS       |
*     |        GLP_NS         |    +    |       GLP_NL       |
*     |        GLP_NS         |    -    |       GLP_NU       |
*     +-----------------------+---------+--------------------+
*
*  Value of row multiplier pi[p] in solution to the original problem is
*  the same as in solution to the transformed problem.
*
*  RECOVERING INTERIOR POINT SOLUTION
*
*  Value of row multiplier pi[p] in solution to the original problem is
*  the same as in solution to the transformed problem.
*
*  RECOVERING MIP SOLUTION
*
*  None needed. */

struct make_equality
{     /* row with almost identical bounds */
      int p;
      /* row reference number */
};

static int rcv_make_equality(NPP *npp, void *info);

int npp_make_equality(NPP *npp, NPPROW *p)
{     /* process row with almost identical bounds */
      struct make_equality *info;
      double b, eps, nint;
      /* the row must be double-sided inequality */
      xassert(p->lb != -DBL_MAX);
      xassert(p->ub != +DBL_MAX);
      xassert(p->lb < p->ub);
      /* check row bounds */
      eps = 1e-9 + 1e-12 * fabs(p->lb);
      if (p->ub - p->lb > eps) return 0;
      /* row bounds are very close to each other */
      /* create transformation stack entry */
      info = npp_push_tse(npp,
         rcv_make_equality, sizeof(struct make_equality));
      info->p = p->i;
      /* compute right-hand side */
      b = 0.5 * (p->ub + p->lb);
      nint = floor(b + 0.5);
      if (fabs(b - nint) <= eps) b = nint;
      /* replace row p by almost equivalent equality constraint */
      p->lb = p->ub = b;
      return 1;
}

int rcv_make_equality(NPP *npp, void *_info)
{     /* recover row with almost identical bounds */
      struct make_equality *info = _info;
      if (npp->sol == GLP_SOL)
      {  if (npp->r_stat[info->p] == GLP_BS)
            npp->r_stat[info->p] = GLP_BS;
         else if (npp->r_stat[info->p] == GLP_NS)
         {  if (npp->r_pi[info->p] >= 0.0)
               npp->r_stat[info->p] = GLP_NL;
            else
               npp->r_stat[info->p] = GLP_NU;
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
*  npp_make_fixed - process column with almost identical bounds
*
*  SYNOPSIS
*
*  #include "glpnpp.h"
*  int npp_make_fixed(NPP *npp, NPPCOL *q);
*
*  DESCRIPTION
*
*  The routine npp_make_fixed processes column q:
*
*     l[q] <= x[q] <= u[q],                                          (1)
*
*  where -oo < l[q] < u[q] < +oo, i.e. which has both lower and upper
*  bounds.
*
*  RETURNS
*
*  0 - column bounds have not been changed;
*
*  1 - column has been fixed.
*
*  PROBLEM TRANSFORMATION
*
*  If bounds of column (1) are very close to each other:
*
*     u[q] - l[q] <= eps,                                            (2)
*
*  where eps is an absolute tolerance for column value, the column can
*  be fixed:
*
*     x[q] = s[q],                                                   (3)
*
*  where s[q] = (l[q] + u[q]) / 2. And if the fixed column value s[q]
*  happens to be very close to its nearest integer:
*
*     |s[q] - floor(s[q] + 0.5)| <= eps,                             (4)
*
*  it is reasonable to use this nearest integer as the fixed value.
*
*  RECOVERING BASIC SOLUTION
*
*  In the dual system of the original (as well as transformed) problem
*  column q corresponds to the following row:
*
*     sum a[i,q] pi[i] + lambda[q] = c[q].                           (5)
*      i
*
*  Since multipliers pi[i] are known for all rows from solution to the
*  transformed problem, formula (5) allows computing value of multiplier
*  (reduced cost) for column q:
*
*     lambda[q] = c[q] - sum a[i,q] pi[i].                           (6)
*                         i
*
*  Status of column q in solution to the original problem is determined
*  by its status and the sign of its multiplier lambda[q] in solution to
*  the transformed problem as follows:
*
*     +-----------------------+-----------+--------------------+
*     |  Status of column q   |  Sign of  | Status of column q |
*     | (transformed problem) | lambda[q] | (original problem) |
*     +-----------------------+-----------+--------------------+
*     |        GLP_BS         |   + / -   |       GLP_BS       |
*     |        GLP_NS         |     +     |       GLP_NL       |
*     |        GLP_NS         |     -     |       GLP_NU       |
*     +-----------------------+-----------+--------------------+
*
*  Value of column q in solution to the original problem is the same as
*  in solution to the transformed problem.
*
*  RECOVERING INTERIOR POINT SOLUTION
*
*  Value of column q in solution to the original problem is the same as
*  in solution to the transformed problem.
*
*  RECOVERING MIP SOLUTION
*
*  None needed. */

struct make_fixed
{     /* column with almost identical bounds */
      int q;
      /* column reference number */
      double c;
      /* objective coefficient at x[q] */
      NPPLFE *ptr;
      /* list of non-zero coefficients a[i,q] */
};

static int rcv_make_fixed(NPP *npp, void *info);

int npp_make_fixed(NPP *npp, NPPCOL *q)
{     /* process column with almost identical bounds */
      struct make_fixed *info;
      NPPAIJ *aij;
      NPPLFE *lfe;
      double s, eps, nint;
      /* the column must be double-bounded */
      xassert(q->lb != -DBL_MAX);
      xassert(q->ub != +DBL_MAX);
      xassert(q->lb < q->ub);
      /* check column bounds */
      eps = 1e-9 + 1e-12 * fabs(q->lb);
      if (q->ub - q->lb > eps) return 0;
      /* column bounds are very close to each other */
      /* create transformation stack entry */
      info = npp_push_tse(npp,
         rcv_make_fixed, sizeof(struct make_fixed));
      info->q = q->j;
      info->c = q->coef;
      info->ptr = NULL;
      /* save column coefficients a[i,q] (needed for basic solution
         only) */
      if (npp->sol == GLP_SOL)
      {  for (aij = q->ptr; aij != NULL; aij = aij->c_next)
         {  lfe = dmp_get_atom(npp->stack, sizeof(NPPLFE));
            lfe->ref = aij->row->i;
            lfe->val = aij->val;
            lfe->next = info->ptr;
            info->ptr = lfe;
         }
      }
      /* compute column fixed value */
      s = 0.5 * (q->ub + q->lb);
      nint = floor(s + 0.5);
      if (fabs(s - nint) <= eps) s = nint;
      /* make column q fixed */
      q->lb = q->ub = s;
      return 1;
}

static int rcv_make_fixed(NPP *npp, void *_info)
{     /* recover column with almost identical bounds */
      struct make_fixed *info = _info;
      NPPLFE *lfe;
      double lambda;
      if (npp->sol == GLP_SOL)
      {  if (npp->c_stat[info->q] == GLP_BS)
            npp->c_stat[info->q] = GLP_BS;
         else if (npp->c_stat[info->q] == GLP_NS)
         {  /* compute multiplier for column q with formula (6) */
            lambda = info->c;
            for (lfe = info->ptr; lfe != NULL; lfe = lfe->next)
               lambda -= lfe->val * npp->r_pi[lfe->ref];
            /* assign status to non-basic column */
            if (lambda >= 0.0)
               npp->c_stat[info->q] = GLP_NL;
            else
               npp->c_stat[info->q] = GLP_NU;
         }
         else
         {  npp_error();
            return 1;
         }
      }
      return 0;
}

/* eof */
