/* covgen.c */

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
#include "fvs.h"
#include "ks.h"
#include "prob.h"

struct glp_cov
{     /* cover cut generator working area */
      int n;
      /* number of columns (variables) */
      glp_prob *set;
      /* set of globally valid 0-1 knapsack inequalities chosen from
       * the root problem; each inequality is either original row or
       * its relaxation (surrogate 0-1 knapsack) which is constructed
       * by substitution of lower/upper single/variable bounds for
       * continuous and general integer (non-binary) variables */
};

struct bnd
{     /* simple or variable bound */
      /* if z = 0, it is a simple bound x >= or <= b; if b = -DBL_MAX
       * (b = +DBL_MAX), x has no lower (upper) bound; otherwise, if
       * z != 0, it is a variable bound x >= or <= a * z + b */
      int z;
      /* number of binary variable or 0 */
      double a, b;
      /* bound parameters */
};

struct csa
{     /* common storage area */
      glp_prob *P;
      /* original (root) MIP */
      struct bnd *l; /* struct bnd l[1+P->n]; */
      /* lower simple/variable bounds of variables */
      struct bnd *u; /* struct bnd u[1+P->n]; */
      /* upper simple/variable bounds of variables */
      glp_prob *set;
      /* see struct glp_cov above */
};

/***********************************************************************
*  init_bounds - initialize bounds of variables with simple bounds
*
*  This routine initializes lower and upper bounds of all variables
*  with simple bounds specified in the original mip. */

static void init_bounds(struct csa *csa)
{     glp_prob *P = csa->P;
      struct bnd *l = csa->l, *u = csa->u;
      int j;
      for (j = 1; j <= P->n; j++)
      {  l[j].z = u[j].z = 0;
         l[j].a = u[j].a = 0;
         l[j].b = glp_get_col_lb(P, j);
         u[j].b = glp_get_col_ub(P, j);
      }
      return;
}

/***********************************************************************
*  check_vb - check variable bound
*
*  This routine checks if the specified i-th row has the form
*
*     a1 * x + a2 * z >= or <= rhs,                                  (1)
*
*  where x is a non-fixed continuous or general integer variable, and
*  z is a binary variable. If it is, the routine converts the row to
*  the following variable lower/upper bound (VLB/VUB) of x:
*
*     x >= or <= a * z + b,                                          (2)
*
*  where a = - a2 / a1, b = rhs / a1. Note that the inequality type is
*  changed to opposite one when a1 < 0.
*
*  If the row is identified as a variable bound, the routine returns
*  GLP_LO for VLB or GLP_UP for VUB and provides the reference numbers
*  of variables x and z and values of a and b. Otherwise, the routine
*  returns zero. */

static int check_vb(struct csa *csa, int i, int *x, int *z, double *a,
      double *b)
{     glp_prob *P = csa->P;
      GLPROW *row;
      GLPAIJ *a1, *a2;
      int type;
      double rhs;
      xassert(1 <= i && i <= P->m);
      row = P->row[i];
      /* check row type */
      switch (row->type)
      {  case GLP_LO:
         case GLP_UP:
            break;
         default:
            return 0;
      }
      /* take first term of the row */
      a1 = row->ptr;
      if (a1 == NULL)
         return 0;
      /* take second term of the row */
      a2 = a1->r_next;
      if (a2 == NULL)
         return 0;
      /* there should be exactly two terms in the row */
      if (a2->r_next != NULL)
         return 0;
      /* if first term is a binary variable, swap the terms */
      if (glp_get_col_kind(P, a1->col->j) == GLP_BV)
      {  GLPAIJ *a;
         a = a1, a1 = a2, a2 = a;
      }
      /* now first term should be a non-fixed continuous or general
       * integer variable */
      if (a1->col->type == GLP_FX)
         return 0;
      if (glp_get_col_kind(P, a1->col->j) == GLP_BV)
         return 0;
      /* and second term should be a binary variable */
      if (glp_get_col_kind(P, a2->col->j) != GLP_BV)
         return 0;
      /* VLB/VUB row has been identified */
      switch (row->type)
      {  case GLP_LO:
            type = a1->val > 0 ? GLP_LO : GLP_UP;
            rhs = row->lb;
            break;
         case GLP_UP:
            type = a1->val > 0 ? GLP_UP : GLP_LO;
            rhs = row->ub;
            break;
         default:
            xassert(type != type);
      }
      *x = a1->col->j;
      *z = a2->col->j;
      *a = - a2->val / a1->val;
      *b = rhs / a1->val;
      return type;
}

/***********************************************************************
*  set_vb - set variable bound
*
*  This routine sets lower or upper variable bound specified as
*
*     x >= a * z + b    (type = GLP_LO)
*
*     x <= a * z + b    (type = GLP_UP) */

static void set_vb(struct csa *csa, int type, int x, int z, double a,
      double b)
{     glp_prob *P = csa->P;
      struct bnd *l = csa->l, *u = csa->u;
      xassert(glp_get_col_type(P, x) != GLP_FX);
      xassert(glp_get_col_kind(P, x) != GLP_BV);
      xassert(glp_get_col_kind(P, z) == GLP_BV);
      xassert(a != 0);
      switch (type)
      {  case GLP_LO:
            /* FIXME: check existing simple lower bound? */
            l[x].z = z, l[x].a = a, l[x].b = b;
            break;
         case GLP_UP:
            /* FIXME: check existing simple upper bound? */
            u[x].z = z, u[x].a = a, u[x].b = b;
            break;
         default:
            xassert(type != type);
      }
      return;
}

/***********************************************************************
*  obtain_vbs - obtain and set variable bounds
*
*  This routine walks thru all rows of the original mip, identifies
*  rows specifying variable lower/upper bounds, and sets these bounds
*  for corresponding (non-binary) variables. */

static void obtain_vbs(struct csa *csa)
{     glp_prob *P = csa->P;
      int i, x, z, type, save;
      double a, b;
      for (i = 1; i <= P->m; i++)
      {  switch (P->row[i]->type)
         {  case GLP_FR:
               break;
            case GLP_LO:
            case GLP_UP:
               type = check_vb(csa, i, &x, &z, &a, &b);
               if (type)
                  set_vb(csa, type, x, z, a, b);
               break;
            case GLP_DB:
            case GLP_FX:
               /* double-side inequality l <= ... <= u and equality
                * ... = l = u are considered as two single inequalities
                * ... >= l and ... <= u */
               save = P->row[i]->type;
               P->row[i]->type = GLP_LO;
               type = check_vb(csa, i, &x, &z, &a, &b);
               if (type)
                  set_vb(csa, type, x, z, a, b);
               P->row[i]->type = GLP_UP;
               type = check_vb(csa, i, &x, &z, &a, &b);
               if (type)
                  set_vb(csa, type, x, z, a, b);
               P->row[i]->type = save;
               break;
            default:
               xassert(P != P);
         }
      }
      return;
}

/***********************************************************************
*  add_term - add term to sparse vector
*
*  This routine computes the following linear combination:
*
*     v := v + a * e[j],
*
*  where v is a sparse vector in full storage format, a is a non-zero
*  scalar, e[j] is j-th column of unity matrix. */

static void add_term(FVS *v, int j, double a)
{     xassert(1 <= j && j <= v->n);
      xassert(a != 0);
      if (v->vec[j] == 0)
      {  /* create j-th component */
         v->nnz++;
         xassert(v->nnz <= v->n);
         v->ind[v->nnz] = j;
      }
      /* perform addition */
      v->vec[j] += a;
      if (fabs(v->vec[j]) < 1e-9 * (1 + fabs(a)))
      {  /* remove j-th component */
         v->vec[j] = DBL_MIN;
      }
      return;
}

/***********************************************************************
*  build_ks - build "0-1 knapsack" inequality
*
*  Given an inequality of "not greater" type:
*
*     sum{j in 1..n} a[j]*x[j] <= b,                                 (1)
*
*  this routine attempts to transform it to equivalent or relaxed "0-1
*  knapsack" inequality that contains only binary variables.
*
*  If x[j] is a binary variable, the term a[j]*x[j] is not changed.
*  Otherwise, if x[j] is a continuous or integer non-binary variable,
*  it is replaced by its lower (if a[j] > 0) or upper (if a[j] < 0)
*  single or variable bound. In the latter case, if x[j] is a non-fixed
*  variable, this results in a relaxation of original inequality known
*  as "surrogate knapsack". Thus, if the specified inequality is valid
*  for the original mip, the resulting inequality is also valid.
*
*  Note that in both source and resulting inequalities coefficients
*  a[j] can have any sign.
*
*  On entry to the routine the source inequality is specified by the
*  parameters n, ind (contains original numbers of x[j]), a, and b. The
*  parameter v is a working sparse vector whose components are assumed
*  to be zero.
*
*  On exit the routine stores the resulting "0-1 knapsack" inequality
*  in the parameters ind, a, and b, and returns n which is the number
*  of terms in the resulting inequality. Zero content of the vector v
*  is restored before exit.
*
*  If the resulting inequality cannot be constructed due to missing
*  lower/upper bounds of some variable, the routine returns a negative
*  value. */

static int build_ks(struct csa *csa, int n, int ind[], double a[],
      double *b, FVS *v)
{     glp_prob *P = csa->P;
      struct bnd *l = csa->l, *u = csa->u;
      int j, k;
      /* check that v = 0 */
#ifdef GLP_DEBUG
      fvs_check_vec(v);
#endif
      xassert(v->nnz == 0);
      /* walk thru terms of original inequality */
      for (j = 1; j <= n; j++)
      {  /* process term a[j]*x[j] */
         k = ind[j]; /* original number of x[j] in mip */
         if (glp_get_col_kind(P, k) == GLP_BV)
         {  /* x[j] is a binary variable */
            /* include its term into resulting inequality */
            add_term(v, k, a[j]);
         }
         else if (a[j] > 0)
         {  /* substitute x[j] by its lower bound */
            if (l[k].b == -DBL_MAX)
            {  /* x[j] has no lower bound */
               n = -1;
               goto skip;
            }
            else if (l[k].z == 0)
            {  /* x[j] has simple lower bound */
               *b -= a[j] * l[k].b;
            }
            else
            {  /* x[j] has variable lower bound (a * z + b) */
               add_term(v, l[k].z, a[j] * l[k].a);
               *b -= a[j] * l[k].b;
            }
         }
         else /* a[j] < 0 */
         {  /* substitute x[j] by its upper bound */
            if (u[k].b == +DBL_MAX)
            {  /* x[j] has no upper bound */
               n = -1;
               goto skip;
            }
            else if (u[k].z == 0)
            {  /* x[j] has simple upper bound */
               *b -= a[j] * u[k].b;
            }
            else
            {  /* x[j] has variable upper bound (a * z + b) */
               add_term(v, u[k].z, a[j] * u[k].a);
               *b -= a[j] * u[k].b;
            }
         }
      }
      /* replace tiny coefficients by exact zeros (see add_term) */
      fvs_adjust_vec(v, 2 * DBL_MIN);
      /* copy terms of resulting inequality */
      xassert(v->nnz <= n);
      n = v->nnz;
      for (j = 1; j <= n; j++)
      {  ind[j] = v->ind[j];
         a[j] = v->vec[ind[j]];
      }
skip: /* restore zero content of v */
      fvs_clear_vec(v);
      return n;
}

/***********************************************************************
*  can_be_active - check if inequality can be active
*
*  This routine checks if the specified "0-1 knapsack" inequality
*
*     sum{j in 1..n} a[j]*x[j] <= b
*
*  can be active. If so, the routine returns true, otherwise false. */

static int can_be_active(int n, const double a[], double b)
{     int j;
      double s;
      s = 0;
      for (j = 1; j <= n; j++)
      {  if (a[j] > 0)
            s += a[j];
      }
      return s > b + .001 * (1 + fabs(b));
}

/***********************************************************************
*  is_sos_ineq - check if inequality is packing (SOS) constraint
*
*  This routine checks if the specified "0-1 knapsack" inequality
*
*     sum{j in 1..n} a[j]*x[j] <= b                                  (1)
*
*  is equivalent to packing inequality (Padberg calls such inequalities
*  special ordered set or SOS constraints)
*
*     sum{j in J'} x[j] - sum{j in J"} x[j] <= 1 - |J"|.             (2)
*
*  If so, the routine returns true, otherwise false.
*
*  Note that if X is a set of feasible binary points satisfying to (2),
*  its convex hull conv(X) equals to the set of feasible points of LP
*  relaxation of (2), which is a n-dimensional simplex, so inequalities
*  (2) are useless for generating cover cuts (due to unimodularity).
*
*  ALGORITHM
*
*  First, we make all a[j] positive by complementing x[j] = 1 - x'[j]
*  in (1). This is performed implicitly (i.e. actually the array a is
*  not changed), but b is replaced by b - sum{j : a[j] < 0}.
*
*  Then we find two smallest coefficients a[p] = min{j in 1..n} a[j]
*  and a[q] = min{j in 1..n : j != p} a[j]. It is obvious that if
*  a[p] + a[q] > b, then a[i] + a[j] > b for all i != j, from which it
*  follows that x[i] + x[j] <= 1 for all i != j. But the latter means
*  that the original inequality (with all a[j] > 0) is equivalent to
*  packing inequality
*
*     sum{j in 1..n} x[j] <= 1.                                      (3)
*
*  Returning to original (uncomplemented) variables x'[j] = 1 - x[j]
*  we have that the original inequality is equivalent to (2), where
*  J' = {j : a[j] > 0} and J" = {j : a[j] < 0}. */

static int is_sos_ineq(int n, const double a[], double b)
{     int j, p, q;
      xassert(n >= 2);
      /* compute b := b - sum{j : a[j] < 0} */
      for (j = 1; j <= n; j++)
      {  if (a[j] < 0)
            b -= a[j];
      }
      /* find a[p] = min{j in 1..n} a[j] */
      p = 1;
      for (j = 2; j <= n; j++)
      {  if (fabs(a[p]) > fabs(a[j]))
            p = j;
      }
      /* find a[q] = min{j in 1..n : j != p} a[j] */
      q = 0;
      for (j = 1; j <= n; j++)
      {  if (j != p)
         {  if (q == 0 || fabs(a[q]) > fabs(a[j]))
               q = j;
         }
      }
      xassert(q != 0);
      /* check condition a[p] + a[q] > b */
      return fabs(a[p]) + fabs(a[q]) > b + .001 * (1 + fabs(b));
}

/***********************************************************************
*  process_ineq - basic inequality processing
*
*  This routine performs basic processing of an inequality of "not
*  greater" type
*
*     sum{j in 1..n} a[j]*x[j] <= b
*
*  specified by the parameters, n, ind, a, and b.
*
*  If the inequality can be transformed to "0-1 knapsack" ineqiality
*  suitable for generating cover cuts, the routine adds it to the set
*  of "0-1 knapsack" inequalities.
*
*  Note that the arrays ind and a are not saved on exit. */

static void process_ineq(struct csa *csa, int n, int ind[], double a[],
      double b, FVS *v)
{     int i;
      /* attempt to transform the specified inequality to equivalent or
       * relaxed "0-1 knapsack" inequality */
      n = build_ks(csa, n, ind, a, &b, v);
      if (n <= 1)
      {  /* uninteresting inequality (in principle, such inequalities
          * should be removed by the preprocessor) */
         goto done;
      }
      if (!can_be_active(n, a, b))
      {  /* inequality is redundant (i.e. cannot be active) */
         goto done;
      }
      if (is_sos_ineq(n, a, b))
      {  /* packing (SOS) inequality is useless for generating cover
          * cuts; currently such inequalities are just ignored */
         goto done;
      }
      /* add resulting "0-1 knapsack" inequality to the set */
      i = glp_add_rows(csa->set, 1);
      glp_set_mat_row(csa->set, i, n, ind, a);
      glp_set_row_bnds(csa->set, i, GLP_UP, b, b);
done: return;
}

/**********************************************************************/

glp_cov *glp_cov_init(glp_prob *P)
{     /* create and initialize cover cut generator */
      glp_cov *cov;
      struct csa csa;
      int i, k, len, *ind;
      double rhs, *val;
      FVS fvs;
      csa.P = P;
      csa.l = talloc(1+P->n, struct bnd);
      csa.u = talloc(1+P->n, struct bnd);
      csa.set = glp_create_prob();
      glp_add_cols(csa.set, P->n);
      /* initialize bounds of variables with simple bounds */
      init_bounds(&csa);
      /* obtain and set variable bounds */
      obtain_vbs(&csa);
      /* allocate working arrays */
      ind = talloc(1+P->n, int);
      val = talloc(1+P->n, double);
      fvs_alloc_vec(&fvs, P->n);
      /* process all rows of the root mip */
      for (i = 1; i <= P->m; i++)
      {  switch (P->row[i]->type)
         {  case GLP_FR:
               break;
            case GLP_LO:
               /* obtain row of ">=" type */
               len = glp_get_mat_row(P, i, ind, val);
               rhs = P->row[i]->lb;
               /* transforms it to row of "<=" type */
               for (k = 1; k <= len; k++)
                  val[k] = - val[k];
               rhs = - rhs;
               /* process the row */
               process_ineq(&csa, len, ind, val, rhs, &fvs);
               break;
            case GLP_UP:
               /* obtain row of "<=" type */
               len = glp_get_mat_row(P, i, ind, val);
               rhs = P->row[i]->ub;
               /* and process it */
               process_ineq(&csa, len, ind, val, rhs, &fvs);
               break;
            case GLP_DB:
            case GLP_FX:
               /* double-sided inequalitiy and equality constraints are
                * processed as two separate inequalities */
               /* obtain row as if it were of ">=" type */
               len = glp_get_mat_row(P, i, ind, val);
               rhs = P->row[i]->lb;
               /* transforms it to row of "<=" type */
               for (k = 1; k <= len; k++)
                  val[k] = - val[k];
               rhs = - rhs;
               /* and process it */
               process_ineq(&csa, len, ind, val, rhs, &fvs);
               /* obtain the same row as if it were of "<=" type */
               len = glp_get_mat_row(P, i, ind, val);
               rhs = P->row[i]->ub;
               /* and process it */
               process_ineq(&csa, len, ind, val, rhs, &fvs);
               break;
            default:
               xassert(P != P);
         }
      }
      /* free working arrays */
      tfree(ind);
      tfree(val);
      fvs_check_vec(&fvs);
      fvs_free_vec(&fvs);
      /* the set of "0-1 knapsack" inequalities has been built */
      if (csa.set->m == 0)
      {  /* the set is empty */
         xprintf("No 0-1 knapsack inequalities detected\n");
         cov = NULL;
         glp_delete_prob(csa.set);
      }
      else
      {  /* create the cover cut generator working area */
         xprintf("Number of 0-1 knapsack inequalities = %d\n",
            csa.set->m);
         cov = talloc(1, glp_cov);
         cov->n = P->n;
         cov->set = csa.set;
#if 0
         glp_write_lp(cov->set, 0, "set.lp");
#endif
      }
      tfree(csa.l);
      tfree(csa.u);
      return cov;
}

/***********************************************************************
*  solve_ks - solve 0-1 knapsack problem
*
*  This routine finds (sub)optimal solution to 0-1 knapsack problem:
*
*     maximize z = sum{j in 1..n} c[j]x[j]                           (1)
*
*         s.t. sum{j in 1..n} a[j]x[j] <= b                          (2)
*
*              x[j] in {0, 1} for all j in 1..n                      (3)
*
*  It is assumed that the instance is non-normalized, i.e. parameters
*  a, b, and c may have any sign.
*
*  On exit the routine stores the (sub)optimal point found in locations
*  x[1], ..., x[n] and returns the optimal objective value. However, if
*  the instance is infeasible, the routine returns INT_MIN. */

static int solve_ks(int n, const int a[], int b, const int c[],
      char x[])
{     int z;
      /* surprisingly, even for some small instances (n = 50-100)
       * MT1 routine takes too much time, so it is used only for tiny
       * instances */
      if (n <= 16)
#if 0
         z = ks_enum(n, a, b, c, x);
#else
         z = ks_mt1(n, a, b, c, x);
#endif
      else
         z = ks_greedy(n, a, b, c, x);
      return z;
}

/***********************************************************************
*  simple_cover - find simple cover cut
*
*  Given a 0-1 knapsack inequality (which may be globally as well as
*  locally valid)
*
*     sum{j in 1..n} a[j]x[j] <= b,                                  (1)
*
*  where all x[j] are binary variables and all a[j] are positive, and
*  a fractional point x~{j in 1..n}, which is feasible to LP relaxation
*  of (1), this routine attempts to find a simple cover inequality
*
*     sum{j in C} (1 - x[j]) >= 1,                                   (2)
*
*  which is valid for (1) and violated at x~.
*
*  Actually, the routine finds a cover C, i.e. a subset of {1, ..., n}
*  such that
*
*     sum{j in C} a[j] > b,                                          (3)
*
*  and which minimizes the left-hand side of (2) at x~
*
*     zeta = sum{j in C} (1 - x~[j]).                                (4)
*
*  On exit the routine stores the characteritic vector z{j in 1..n}
*  of the cover found (i.e. z[j] = 1 means j in C, and z[j] = 0 means
*  j not in C), and returns corresponding minimal value of zeta (4).
*  However, if no cover is found, the routine returns DBL_MAX.
*
*  ALGORITHM
*
*  The separation problem (3)-(4) is converted to 0-1 knapsack problem
*  as follows.
*
*  First, note that the constraint (3) is equivalent to
*
*     sum{j in 1..n} a[j]z[j] >= b + eps,                            (5)
*
*  where eps > 0 is a sufficiently small number (in case of integral
*  a and b we may take eps = 1). Multiplying both sides of (5) by (-1)
*  gives
*
*     sum{j in 1..n} (-a[j])z[j] <= - b - eps.                       (6)
*
*  To make all coefficients in (6) positive, z[j] is complemented by
*  substitution z[j] = 1 - z'[j] that finally gives
*
*     sum{j in 1..n} a[j]z'[j] <= sum{j in 1..n} a[j] - b - eps.     (7)
*
*  Minimization of zeta (4) is equivalent to maximization of
*
*     -zeta = sum{j in 1..n} (x~[j] - 1)z[j].                        (8)
*
*  Substitution z[j] = 1 - z'[j] gives
*
*     -zeta = sum{j in 1..n} (1 - x~[j])z'[j] - zeta0,               (9)
*
*  where zeta0 = sum{j in 1..n} (1 - x~[j]) is a constant term.
*
*  Thus, the 0-1 knapsack problem to be solved is the following:
*
*     maximize
*
*        -zeta = sum{j in 1..n} (1 - x~[j])z'[j] - zeta0            (10)
*
*     subject to
*
*        sum{j in 1..n} a[j]z'[j] <= sum{j in 1..n} a[j] - b - eps  (11)
*
*        z'[j] in {0,1} for all j = 1,...,n                         (12)
*
*  (The constant term zeta0 doesn't affect the solution, so it can be
*  dropped.) */

static double simple_cover(int n, const double a[], double b, const
      double x[], char z[])
{     int j, *aa, bb, *cc;
      double max_aj, min_aj, s, eps;
      xassert(n >= 3);
      /* allocate working arrays */
      aa = talloc(1+n, int);
      cc = talloc(1+n, int);
      /* compute max{j in 1..n} a[j] and min{j in 1..n} a[j] */
      max_aj = 0, min_aj = DBL_MAX;
      for (j = 1; j <= n; j++)
      {  xassert(a[j] > 0);
         if (max_aj < a[j])
            max_aj = a[j];
         if (min_aj > a[j])
            min_aj = a[j];
      }
      /* scale and round constraint parameters to make them integral;
       * note that we make the resulting inequality stronger than (11),
       * so a[j]'s are rounded up while rhs is rounded down */
      s = 0;
      for (j = 1; j <= n; j++)
      {  s += a[j];
         aa[j] = ceil(a[j] / max_aj * 1000);
      }
      bb = floor((s - b) / max_aj * 1000) - 1;
      /* scale and round obj. coefficients to make them integral;
       * again we make the objective function stronger than (10), so
       * the coefficients are rounded down */
      for (j = 1; j <= n; j++)
      {  xassert(0 <= x[j] && x[j] <= 1);
         cc[j] = floor((1 - x[j]) * 1000);
      }
      /* solve separation problem */
      if (solve_ks(n, aa, bb, cc, z) == INT_MIN)
      {  /* no cover exists */
         s = DBL_MAX;
         goto skip;
      }
      /* determine z[j] = 1 - z'[j] */
      for (j = 1; j <= n; j++)
      {  xassert(z[j] == 0 || z[j] == 1);
         z[j] ^= 1;
      }
      /* check condition (11) for original (non-scaled) parameters */
      s = 0;
      for (j = 1; j <= n; j++)
      {  if (z[j])
            s += a[j];
      }
      eps = 0.01 * (min_aj >= 1 ? min_aj : 1);
      if (!(s >= b + eps))
      {  /* no cover found within a precision req'd */
         s = DBL_MAX;
         goto skip;
      }
      /* compute corresponding zeta (4) for cover found */
      s = 0;
      for (j = 1; j <= n; j++)
      {  if (z[j])
            s += 1 - x[j];
      }
skip: /* free working arrays */
      tfree(aa);
      tfree(cc);
      return s;
}

/**********************************************************************/

void glp_cov_gen1(glp_prob *P, glp_cov *cov, glp_prob *pool)
{     /* generate locally valid simple cover cuts */
      int i, k, len, new_len, *ind;
      double *val, rhs, *x, zeta;
      char *z;
      xassert(P->n == cov->n && P->n == cov->set->n);
      xassert(glp_get_status(P) == GLP_OPT);
      /* allocate working arrays */
      ind = talloc(1+P->n, int);
      val = talloc(1+P->n, double);
      x = talloc(1+P->n, double);
      z = talloc(1+P->n, char);
      /* walk thru 0-1 knapsack inequalities */
      for (i = 1; i <= cov->set->m; i++)
      {  /* retrieve 0-1 knapsack inequality */
         len = glp_get_mat_row(cov->set, i, ind, val);
         rhs = glp_get_row_ub(cov->set, i);
         xassert(rhs != +DBL_MAX);
         /* FIXME: skip, if slack is too large? */
         /* substitute and eliminate binary variables which have been
          * fixed in the current subproblem (this makes the inequality
          * only locally valid) */
         new_len = 0;
         for (k = 1; k <= len; k++)
         {  if (glp_get_col_type(P, ind[k]) == GLP_FX)
               rhs -= val[k] * glp_get_col_prim(P, ind[k]);
            else
            {  new_len++;
               ind[new_len] = ind[k];
               val[new_len] = val[k];
            }
         }
         len = new_len;
         /* we need at least 3 binary variables in the inequality */
         if (len <= 2)
            continue;
         /* obtain values of binary variables from optimal solution to
          * LP relaxation of current subproblem */
         for (k = 1; k <= len; k++)
         {  xassert(glp_get_col_kind(P, ind[k]) == GLP_BV);
            x[k] = glp_get_col_prim(P, ind[k]);
            if (x[k] < 0.00001)
               x[k] = 0;
            else if (x[k] > 0.99999)
               x[k] = 1;
            /* if val[k] < 0, perform substitution x[k] = 1 - x'[k] to
             * make all coefficients positive */
            if (val[k] < 0)
            {  ind[k] = - ind[k]; /* x[k] is complemented */
               val[k] = - val[k];
               rhs += val[k];
               x[k] = 1 - x[k];
            }
         }
         /* find locally valid simple cover cut */
         zeta = simple_cover(len, val, rhs, x, z);
         if (zeta > 0.95)
         {  /* no violation or insufficient violation; see (2) */
            continue;
         }
         /* construct cover inequality (2) for the cover found, which
          * for original binary variables x[k] is equivalent to:
          *    sum{k in C'} x[k] + sum{k in C"} x'[k] <= |C| - 1
          * or
          *    sum{k in C'} x[k] + sum{k in C"} (1 - x[k]) <= |C| - 1
          * or
          *    sum{k in C'} x[k] - sum{k in C"} x[k] <= |C'| - 1
          * since |C| - |C"| = |C'| */
         new_len = 0;
         rhs = -1;
         for (k = 1; k <= len; k++)
         {  if (z[k])
            {  new_len++;
               if (ind[k] > 0)
               {  ind[new_len] = +ind[k];
                  val[new_len] = +1;
                  rhs++;
               }
               else /* ind[k] < 0 */
               {  ind[new_len] = -ind[k];
                  val[new_len] = -1;
               }
            }
         }
         len = new_len;
         /* add the cover inequality to the local cut pool */
         k = glp_add_rows(pool, 1);
         glp_set_mat_row(pool, k, len, ind, val);
         glp_set_row_bnds(pool, k, GLP_UP, rhs, rhs);
      }
      /* free working arrays */
      tfree(ind);
      tfree(val);
      tfree(x);
      tfree(z);
      return;
}

/**********************************************************************/

void glp_cov_free(glp_cov *cov)
{     /* delete cover cut generator workspace */
      xassert(cov != NULL);
      glp_delete_prob(cov->set);
      tfree(cov);
      return;
}

/* eof */
