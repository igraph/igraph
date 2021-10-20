/* spydual.c */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2015-2017 Free Software Foundation, Inc.
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

#if 1 /* 18/VII-2017 */
#define SCALE_Z 1
#endif

#include "env.h"
#include "simplex.h"
#include "spxat.h"
#include "spxnt.h"
#include "spxprob.h"
#include "spychuzc.h"
#include "spychuzr.h"
#if 0 /* 11/VI-2017 */
#if 1 /* 29/III-2016 */
#include "fvs.h"
#endif
#endif

#define CHECK_ACCURACY 0
/* (for debugging) */

struct csa
{     /* common storage area */
      SPXLP *lp;
      /* LP problem data and its (current) basis; this LP has m rows
       * and n columns */
      int dir;
      /* original optimization direction:
       * +1 - minimization
       * -1 - maximization */
#if SCALE_Z
      double fz;
      /* factor used to scale original objective */
#endif
      double *orig_b; /* double orig_b[1+m]; */
      /* copy of original right-hand sides */
      double *orig_c; /* double orig_c[1+n]; */
      /* copy of original objective coefficients */
      double *orig_l; /* double orig_l[1+n]; */
      /* copy of original lower bounds */
      double *orig_u; /* double orig_u[1+n]; */
      /* copy of original upper bounds */
      SPXAT *at;
      /* mxn-matrix A of constraint coefficients, in sparse row-wise
       * format (NULL if not used) */
      SPXNT *nt;
      /* mx(n-m)-matrix N composed of non-basic columns of constraint
       * matrix A, in sparse row-wise format (NULL if not used) */
      int phase;
      /* search phase:
       * 0 - not determined yet
       * 1 - searching for dual feasible solution
       * 2 - searching for optimal solution */
      double *beta; /* double beta[1+m]; */
      /* beta[i] is primal value of basic variable xB[i] */
      int beta_st;
      /* status of the vector beta:
       * 0 - undefined
       * 1 - just computed
       * 2 - updated */
      double *d; /* double d[1+n-m]; */
      /* d[j] is reduced cost of non-basic variable xN[j] */
      int d_st;
      /* status of the vector d:
       * 0 - undefined
       * 1 - just computed
       * 2 - updated */
      SPYSE *se;
      /* dual projected steepest edge and Devex pricing data block
       * (NULL if not used) */
#if 0 /* 30/III-2016 */
      int num;
      /* number of eligible basic variables */
      int *list; /* int list[1+m]; */
      /* list[1], ..., list[num] are indices i of eligible basic
       * variables xB[i] */
#else
      FVS r; /* FVS r[1:m]; */
      /* vector of primal infeasibilities */
      /* r->nnz = num; r->ind = list */
      /* vector r has the same status as vector beta (see above) */
#endif
      int p;
      /* xB[p] is a basic variable chosen to leave the basis */
#if 0 /* 29/III-2016 */
      double *trow; /* double trow[1+n-m]; */
#else
      FVS trow; /* FVS trow[1:n-m]; */
#endif
      /* p-th (pivot) row of the simplex table */
#if 1 /* 16/III-2016 */
      SPYBP *bp; /* SPYBP bp[1+n-m]; */
      /* dual objective break-points */
#endif
      int q;
      /* xN[q] is a non-basic variable chosen to enter the basis */
#if 0 /* 29/III-2016 */
      double *tcol; /* double tcol[1+m]; */
#else
      FVS tcol; /* FVS tcol[1:m]; */
#endif
      /* q-th (pivot) column of the simplex table */
      double *work; /* double work[1+m]; */
      /* working array */
      double *work1; /* double work1[1+n-m]; */
      /* another working array */
#if 0 /* 11/VI-2017 */
#if 1 /* 31/III-2016 */
      FVS wrow; /* FVS wrow[1:n-m]; */
      FVS wcol; /* FVS wcol[1:m]; */
      /* working sparse vectors */
#endif
#endif
      int p_stat, d_stat;
      /* primal and dual solution statuses */
      /*--------------------------------------------------------------*/
      /* control parameters (see struct glp_smcp) */
      int msg_lev;
      /* message level */
      int dualp;
      /* if this flag is set, report failure in case of instability */
#if 0 /* 16/III-2016 */
      int harris;
      /* dual ratio test technique:
       * 0 - textbook ratio test
       * 1 - Harris' two pass ratio test */
#else
      int r_test;
      /* dual ratio test technique:
       * GLP_RT_STD  - textbook ratio test
       * GLP_RT_HAR  - Harris' two pass ratio test
       * GLP_RT_FLIP - long-step (flip-flop) ratio test */
#endif
      double tol_bnd, tol_bnd1;
      /* primal feasibility tolerances */
      double tol_dj, tol_dj1;
      /* dual feasibility tolerances */
      double tol_piv;
      /* pivot tolerance */
      double obj_lim;
      /* objective limit */
      int it_lim;
      /* iteration limit */
      int tm_lim;
      /* time limit, milliseconds */
      int out_frq;
#if 0 /* 15/VII-2017 */
      /* display output frequency, iterations */
#else
      /* display output frequency, milliseconds */
#endif
      int out_dly;
      /* display output delay, milliseconds */
      /*--------------------------------------------------------------*/
      /* working parameters */
      double tm_beg;
      /* time value at the beginning of the search */
      int it_beg;
      /* simplex iteration count at the beginning of the search */
      int it_cnt;
      /* simplex iteration count; it increases by one every time the
       * basis changes */
      int it_dpy;
      /* simplex iteration count at most recent display output */
#if 1 /* 15/VII-2017 */
      double tm_dpy;
      /* time value at most recent display output */
#endif
      int inv_cnt;
      /* basis factorization count since most recent display output */
#if 1 /* 11/VII-2017 */
      int degen;
      /* count of successive degenerate iterations; this count is used
       * to detect stalling */
#endif
#if 1 /* 23/III-2016 */
      int ns_cnt, ls_cnt;
      /* normal and long-step iteration count */
#endif
};

/***********************************************************************
*  check_flags - check correctness of active bound flags
*
*  This routine checks that flags specifying active bounds of all
*  non-basic variables are correct.
*
*  NOTE: It is important to note that if bounds of variables have been
*  changed, active bound flags should be corrected accordingly. */

static void check_flags(struct csa *csa)
{     SPXLP *lp = csa->lp;
      int m = lp->m;
      int n = lp->n;
      double *l = lp->l;
      double *u = lp->u;
      int *head = lp->head;
      char *flag = lp->flag;
      int j, k;
      for (j = 1; j <= n-m; j++)
      {  k = head[m+j]; /* x[k] = xN[j] */
         if (l[k] == -DBL_MAX && u[k] == +DBL_MAX)
            xassert(!flag[j]);
         else if (l[k] != -DBL_MAX && u[k] == +DBL_MAX)
            xassert(!flag[j]);
         else if (l[k] == -DBL_MAX && u[k] != +DBL_MAX)
            xassert(flag[j]);
         else if (l[k] == u[k])
            xassert(!flag[j]);
      }
      return;
}

/***********************************************************************
*  set_art_bounds - set artificial right-hand sides and bounds
*
*  This routine sets artificial right-hand sides and artificial bounds
*  for all variables to minimize the sum of dual infeasibilities on
*  phase I. Given current reduced costs d = (d[j]) this routine also
*  sets active artificial bounds of non-basic variables to provide dual
*  feasibility (this is always possible because all variables have both
*  lower and upper artificial bounds). */

static void set_art_bounds(struct csa *csa)
{     SPXLP *lp = csa->lp;
      int m = lp->m;
      int n = lp->n;
      double *b = lp->b;
      double *l = lp->l;
      double *u = lp->u;
      int *head = lp->head;
      char *flag = lp->flag;
      double *d = csa->d;
      int i, j, k;
#if 1 /* 31/III-2016: FIXME */
      /* set artificial right-hand sides */
      for (i = 1; i <= m; i++)
         b[i] = 0.0;
      /* set artificial bounds depending on types of variables */
      for (k = 1; k <= n; k++)
      {  if (csa->orig_l[k] == -DBL_MAX && csa->orig_u[k] == +DBL_MAX)
         {  /* force free variables to enter the basis */
            l[k] = -1e3, u[k] = +1e3;
         }
      else if (csa->orig_l[k] != -DBL_MAX && csa->orig_u[k] == +DBL_MAX)
            l[k] = 0.0, u[k] = +1.0;
      else if (csa->orig_l[k] == -DBL_MAX && csa->orig_u[k] != +DBL_MAX)
            l[k] = -1.0, u[k] = 0.0;
         else
            l[k] = u[k] = 0.0;
      }
#endif
      /* set active artificial bounds for non-basic variables */
      xassert(csa->d_st == 1);
      for (j = 1; j <= n-m; j++)
      {  k = head[m+j]; /* x[k] = xN[j] */
         flag[j] = (l[k] != u[k] && d[j] < 0.0);
      }
      /* invalidate values of basic variables, since active bounds of
       * non-basic variables have been changed */
      csa->beta_st = 0;
      return;
}

/***********************************************************************
*  set_orig_bounds - restore original right-hand sides and bounds
*
*  This routine restores original right-hand sides and original bounds
*  for all variables. This routine also sets active original bounds for
*  non-basic variables; for double-bounded non-basic variables current
*  reduced costs d = (d[j]) are used to decide which bound (lower or
*  upper) should be made active. */

static void set_orig_bounds(struct csa *csa)
{     SPXLP *lp = csa->lp;
      int m = lp->m;
      int n = lp->n;
      double *b = lp->b;
      double *l = lp->l;
      double *u = lp->u;
      int *head = lp->head;
      char *flag = lp->flag;
      double *d = csa->d;
      int j, k;
      /* restore original right-hand sides */
      memcpy(b, csa->orig_b, (1+m) * sizeof(double));
      /* restore original bounds of all variables */
      memcpy(l, csa->orig_l, (1+n) * sizeof(double));
      memcpy(u, csa->orig_u, (1+n) * sizeof(double));
      /* set active original bounds for non-basic variables */
      xassert(csa->d_st == 1);
      for (j = 1; j <= n-m; j++)
      {  k = head[m+j]; /* x[k] = xN[j] */
         if (l[k] == -DBL_MAX && u[k] == +DBL_MAX)
            flag[j] = 0;
         else if (l[k] != -DBL_MAX && u[k] == +DBL_MAX)
            flag[j] = 0;
         else if (l[k] == -DBL_MAX && u[k] != +DBL_MAX)
            flag[j] = 1;
         else if (l[k] != u[k])
            flag[j] = (d[j] < 0.0);
         else
            flag[j] = 0;
      }
      /* invalidate values of basic variables, since active bounds of
       * non-basic variables have been changed */
      csa->beta_st = 0;
      return;
}

/***********************************************************************
*  check_feas - check dual feasibility of basic solution
*
*  This routine checks that reduced costs of all non-basic variables
*  d = (d[j]) have correct signs.
*
*  Reduced cost d[j] is considered as having correct sign within the
*  specified tolerance depending on status of non-basic variable xN[j]
*  if one of the following conditions is met:
*
*     xN[j] is free                       -eps <= d[j] <= +eps
*
*     xN[j] has its lower bound active    d[j] >= -eps
*
*     xN[j] has its upper bound active    d[j] <= +eps
*
*     xN[j] is fixed                      d[j] has any value
*
*  where eps = tol + tol1 * |cN[j]|, cN[j] is the objective coefficient
*  at xN[j]. (See also the routine spx_chuzc_sel.)
*
*  The flag recov allows the routine to recover dual feasibility by
*  changing active bounds of non-basic variables. (For example, if
*  xN[j] has its lower bound active and d[j] < -eps, the feasibility
*  can be recovered by making xN[j] active on its upper bound.)
*
*  If the basic solution is dual feasible, the routine returns zero.
*  If the basic solution is dual infeasible, but its dual feasibility
*  can be recovered (or has been recovered, if the flag recov is set),
*  the routine returns a negative value. Otherwise, the routine returns
*  the number j of some non-basic variable xN[j], whose reduced cost
*  d[j] is dual infeasible and cannot be recovered. */

static int check_feas(struct csa *csa, double tol, double tol1,
      int recov)
{     SPXLP *lp = csa->lp;
      int m = lp->m;
      int n = lp->n;
      double *c = lp->c;
      double *l = lp->l;
      double *u = lp->u;
      int *head = lp->head;
      char *flag = lp->flag;
      double *d = csa->d;
      int j, k, ret = 0;
      double eps;
      /* reduced costs should be just computed */
      xassert(csa->d_st == 1);
      /* walk thru list of non-basic variables */
      for (j = 1; j <= n-m; j++)
      {  k = head[m+j]; /* x[k] = xN[j] */
         if (l[k] == u[k])
         {  /* xN[j] is fixed variable; skip it */
            continue;
         }
         /* determine absolute tolerance eps[j] */
         eps = tol + tol1 * (c[k] >= 0.0 ? +c[k] : -c[k]);
         /* check dual feasibility of xN[j] */
         if (d[j] > +eps)
         {  /* xN[j] should have its lower bound active */
            if (l[k] == -DBL_MAX || flag[j])
            {  /* but it either has no lower bound or its lower bound
                * is inactive */
               if (l[k] == -DBL_MAX)
               {  /* cannot recover, since xN[j] has no lower bound */
                  ret = j;
                  break;
               }
               /* recovering is possible */
               if (recov)
                  flag[j] = 0;
               ret = -1;
            }
         }
         else if (d[j] < -eps)
         {  /* xN[j] should have its upper bound active */
            if (!flag[j])
            {  /* but it either has no upper bound or its upper bound
                * is inactive */
               if (u[k] == +DBL_MAX)
               {  /* cannot recover, since xN[j] has no upper bound */
                  ret = j;
                  break;
               }
               /* recovering is possible */
               if (recov)
                  flag[j] = 1;
               ret = -1;
            }
         }
      }
      if (recov && ret)
      {  /* invalidate values of basic variables, since active bounds
          * of non-basic variables have been changed */
         csa->beta_st = 0;
      }
      return ret;
}

#if CHECK_ACCURACY
/***********************************************************************
*  err_in_vec - compute maximal relative error between two vectors
*
*  This routine computes and returns maximal relative error between
*  n-vectors x and y:
*
*     err_max = max |x[i] - y[i]| / (1 + |x[i]|).
*
*  NOTE: This routine is intended only for debugging purposes. */

static double err_in_vec(int n, const double x[], const double y[])
{     int i;
      double err, err_max;
      err_max = 0.0;
      for (i = 1; i <= n; i++)
      {  err = fabs(x[i] - y[i]) / (1.0 + fabs(x[i]));
         if (err_max < err)
            err_max = err;
      }
      return err_max;
}
#endif

#if CHECK_ACCURACY
/***********************************************************************
*  err_in_beta - compute maximal relative error in vector beta
*
*  This routine computes and returns maximal relative error in vector
*  of values of basic variables beta = (beta[i]).
*
*  NOTE: This routine is intended only for debugging purposes. */

static double err_in_beta(struct csa *csa)
{     SPXLP *lp = csa->lp;
      int m = lp->m;
      double err, *beta;
      beta = talloc(1+m, double);
      spx_eval_beta(lp, beta);
      err = err_in_vec(m, beta, csa->beta);
      tfree(beta);
      return err;
}
#endif

#if CHECK_ACCURACY
static double err_in_r(struct csa *csa)
{     SPXLP *lp = csa->lp;
      int m = lp->m;
      int i, k;
      double err, *r;
      r = talloc(1+m, double);
      for (i = 1; i <= m; i++)
      {  k = lp->head[i];
         if (csa->beta[i] < lp->l[k])
            r[i] = lp->l[k] - csa->beta[i];
         else if (csa->beta[i] > lp->u[k])
            r[i] = lp->u[k] - csa->beta[i];
         else
            r[i] = 0.0;

if (fabs(r[i] - csa->r.vec[i]) > 1e-6)
printf("i = %d; r = %g; csa->r = %g\n", i, r[i], csa->r.vec[i]);


      }
      err = err_in_vec(m, r, csa->r.vec);
      tfree(r);
      return err;
}
#endif

#if CHECK_ACCURACY
/***********************************************************************
*  err_in_d - compute maximal relative error in vector d
*
*  This routine computes and returns maximal relative error in vector
*  of reduced costs of non-basic variables d = (d[j]).
*
*  NOTE: This routine is intended only for debugging purposes. */

static double err_in_d(struct csa *csa)
{     SPXLP *lp = csa->lp;
      int m = lp->m;
      int n = lp->n;
      int j;
      double err, *pi, *d;
      pi = talloc(1+m, double);
      d = talloc(1+n-m, double);
      spx_eval_pi(lp, pi);
      for (j = 1; j <= n-m; j++)
         d[j] = spx_eval_dj(lp, pi, j);
      err = err_in_vec(n-m, d, csa->d);
      tfree(pi);
      tfree(d);
      return err;
}
#endif

#if CHECK_ACCURACY
/***********************************************************************
*  err_in_gamma - compute maximal relative error in vector gamma
*
*  This routine computes and returns maximal relative error in vector
*  of projected steepest edge weights gamma = (gamma[j]).
*
*  NOTE: This routine is intended only for debugging purposes. */

static double err_in_gamma(struct csa *csa)
{     SPXLP *lp = csa->lp;
      int m = lp->m;
      int n = lp->n;
      SPYSE *se = csa->se;
      int i;
      double err, *gamma;
      xassert(se != NULL);
gamma = talloc(1+m, double);
      for (i = 1; i <= m; i++)
         gamma[i] = spy_eval_gamma_i(lp, se, i);
      err = err_in_vec(m, gamma, se->gamma);
      tfree(gamma);
      return err;
}
#endif

#if CHECK_ACCURACY
/***********************************************************************
*  check_accuracy - check accuracy of basic solution components
*
*  This routine checks accuracy of current basic solution components.
*
*  NOTE: This routine is intended only for debugging purposes. */

static void check_accuracy(struct csa *csa)
{     double e_beta, e_r, e_d, e_gamma;
      e_beta = err_in_beta(csa);
      e_r = err_in_r(csa);
      e_d = err_in_d(csa);
      if (csa->se == NULL)
         e_gamma = 0.;
      else
         e_gamma = err_in_gamma(csa);
      xprintf("e_beta = %10.3e; e_r = %10.3e; e_d = %10.3e; e_gamma = %"
         "10.3e\n", e_beta, e_r, e_d, e_gamma);
      xassert(e_beta <= 1e-5 && e_d <= 1e-5 && e_gamma <= 1e-3);
      return;
}
#endif

#if 1 /* 30/III-2016 */
static
void spy_eval_r(SPXLP *lp, const double beta[/*1+m*/], double tol,
      double tol1, FVS *r)
{     /* this routine computes the vector of primal infeasibilities:
       *
       *        ( lB[i] - beta[i] > 0, if beta[i] < lb[i]
       * r[i] = { 0,                   if lb[i] <= beta[i] <= ub[i]
       *        ( ub[i] - beta[i] < 0, if beta[i] > ub[i]
       *
       * (this routine replaces spy_chuzr_sel) */
      int m = lp->m;
      double *l = lp->l;
      double *u = lp->u;
      int *head = lp->head;
      int *ind = r->ind;
      double *vec = r->vec;
      int i, k, nnz = 0;
      double lk, uk, eps;
      xassert(r->n == m);
      /* walk thru the list of basic variables */
      for (i = 1; i <= m; i++)
      {  vec[i] = 0.0;
         k = head[i]; /* x[k] = xB[i] */
         lk = l[k], uk = u[k];
         /* check primal feasibility */
         if (beta[i] < lk)
         {  /* determine absolute tolerance eps1[i] */
            eps = tol + tol1 * (lk >= 0.0 ? +lk : -lk);
            if (beta[i] < lk - eps)
            {  /* lower bound is violated */
               ind[++nnz] = i;
               vec[i] = lk - beta[i];
            }
         }
         else if (beta[i] > uk)
         {  /* determine absolute tolerance eps2[i] */
            eps = tol + tol1 * (uk >= 0.0 ? +uk : -uk);
            if (beta[i] > uk + eps)
            {  /* upper bound is violated */
               ind[++nnz] = i;
               vec[i] = uk - beta[i];
            }
         }
      }
      r->nnz = nnz;
      return;
}
#endif

/***********************************************************************
*  choose_pivot - choose xB[p] and xN[q]
*
*  Given the list of eligible basic variables this routine first
*  chooses basic variable xB[p]. This choice is always possible,
*  because the list is assumed to be non-empty. Then the routine
*  computes p-th row T[p,*] of the simplex table T[i,j] and chooses
*  non-basic variable xN[q]. If the pivot T[p,q] is small in magnitude,
*  the routine attempts to choose another xB[p] and xN[q] in order to
*  avoid badly conditioned adjacent bases.
*
*  If the normal choice was made, the routine returns zero. Otherwise,
*  if the long-step choice was made, the routine returns non-zero. */

#ifdef TIMING /* 31/III-2016 */

#include "choose_pivot.c"

#else

#define MIN_RATIO 0.0001

static int choose_pivot(struct csa *csa)
{     SPXLP *lp = csa->lp;
      int m = lp->m;
      int n = lp->n;
      double *l = lp->l;
      double *u = lp->u;
      int *head = lp->head;
      SPXAT *at = csa->at;
      SPXNT *nt = csa->nt;
      double *beta = csa->beta;
      double *d = csa->d;
      SPYSE *se = csa->se;
#if 0 /* 30/III-2016 */
      int *list = csa->list;
#else
      int *list = csa->r.ind;
#endif
      double *rho = csa->work;
      double *trow = csa->work1;
      SPYBP *bp = csa->bp;
      double tol_piv = csa->tol_piv;
      int try, nnn, j, k, p, q, t, t_best, nbp, ret;
      double big, temp, r, best_ratio, dz_best;
      xassert(csa->beta_st);
      xassert(csa->d_st);
more: /* initial number of eligible basic variables */
#if 0 /* 30/III-2016 */
      nnn = csa->num;
#else
      nnn = csa->r.nnz;
#endif
      /* nothing has been chosen so far */
      csa->p = 0;
      best_ratio = 0.0;
      try = ret = 0;
try:  /* choose basic variable xB[p] */
      xassert(nnn > 0);
      try++;
      if (se == NULL)
      {  /* dual Dantzig's rule */
         p = spy_chuzr_std(lp, beta, nnn, list);
      }
      else
      {  /* dual projected steepest edge */
         p = spy_chuzr_pse(lp, se, beta, nnn, list);
      }
      xassert(1 <= p && p <= m);
      /* compute p-th row of inv(B) */
      spx_eval_rho(lp, p, rho);
      /* compute p-th row of the simplex table */
      if (at != NULL)
         spx_eval_trow1(lp, at, rho, trow);
      else
         spx_nt_prod(lp, nt, trow, 1, -1.0, rho);
#if 1 /* 23/III-2016 */
      /* big := max(1, |trow[1]|, ..., |trow[n-m]|) */
      big = 1.0;
      for (j = 1; j <= n-m; j++)
      {  temp = trow[j];
         if (temp < 0.0)
            temp = - temp;
         if (big < temp)
            big = temp;
      }
#else
      /* this still puzzles me */
      big = 1.0;
#endif
      /* choose non-basic variable xN[q] */
      k = head[p]; /* x[k] = xB[p] */
      xassert(beta[p] < l[k] || beta[p] > u[k]);
      r = beta[p] < l[k] ? l[k] - beta[p] : u[k] - beta[p];
      if (csa->r_test == GLP_RT_FLIP && try <= 2)
      {  /* long-step ratio test */
#if 0 /* 23/III-2016 */
         /* determine dual objective break-points */
         nbp = spy_eval_bp(lp, d, r, trow, tol_piv, bp);
         if (nbp <= 1)
            goto skip;
         /* choose appropriate break-point */
         t_best = 0, dz_best = -DBL_MAX;
         for (t = 1; t <= nbp; t++)
         {  if (fabs(trow[bp[t].j]) / big >= MIN_RATIO)
            {  if (dz_best < bp[t].dz)
                  t_best = t, dz_best = bp[t].dz;
            }
         }
         if (t_best == 0)
            goto skip;
#else
         int t, num, num1;
         double slope, teta_lim;
         /* determine dual objective break-points */
         nbp = spy_ls_eval_bp(lp, d, r, trow, tol_piv, bp);
         if (nbp < 2)
            goto skip;
         /* set initial slope */
         slope = fabs(r);
         /* estimate initial teta_lim */
         teta_lim = DBL_MAX;
         for (t = 1; t <= nbp; t++)
         {  if (teta_lim > bp[t].teta)
               teta_lim = bp[t].teta;
         }
         xassert(teta_lim >= 0.0);
         if (teta_lim < 1e-6)
            teta_lim = 1e-6;
         /* nothing has been chosen so far */
         t_best = 0, dz_best = 0.0, num = 0;
         /* choose appropriate break-point */
         while (num < nbp)
         {  /* select and process a new portion of break-points */
            num1 = spy_ls_select_bp(lp, trow, nbp, bp, num, &slope,
               teta_lim);
            for (t = num+1; t <= num1; t++)
            {  if (fabs(trow[bp[t].j]) / big >= MIN_RATIO)
               {  if (dz_best < bp[t].dz)
                     t_best = t, dz_best = bp[t].dz;
               }
            }
            if (slope < 0.0)
            {  /* the dual objective starts decreasing */
               break;
            }
            /* the dual objective continues increasing */
            num = num1;
            teta_lim += teta_lim;
         }
         if (dz_best == 0.0)
            goto skip;
         xassert(1 <= t_best && t_best <= num1);
#endif
         /* the choice has been made */
         csa->p = p;
#if 0 /* 29/III-2016 */
         memcpy(&csa->trow[1], &trow[1], (n-m) * sizeof(double));
#else
         memcpy(&csa->trow.vec[1], &trow[1], (n-m) * sizeof(double));
         fvs_gather_vec(&csa->trow, DBL_EPSILON);
#endif
         csa->q = bp[t_best].j;
         best_ratio = fabs(trow[bp[t_best].j]) / big;
#if 0
         xprintf("num = %d; t_best = %d; dz = %g\n", num, t_best,
            bp[t_best].dz);
#endif
         ret = 1;
         goto done;
skip:    ;
      }
      if (csa->r_test == GLP_RT_STD)
      {  /* textbook dual ratio test */
         q = spy_chuzc_std(lp, d, r, trow, tol_piv,
            .30 * csa->tol_dj, .30 * csa->tol_dj1);
      }
      else
      {  /* Harris' two-pass dual ratio test */
         q = spy_chuzc_harris(lp, d, r, trow, tol_piv,
            .35 * csa->tol_dj, .35 * csa->tol_dj1);
      }
      if (q == 0)
      {  /* dual unboundedness */
         csa->p = p;
#if 0 /* 29/III-2016 */
         memcpy(&csa->trow[1], &trow[1], (n-m) * sizeof(double));
#else
         memcpy(&csa->trow.vec[1], &trow[1], (n-m) * sizeof(double));
         fvs_gather_vec(&csa->trow, DBL_EPSILON);
#endif
         csa->q = q;
         best_ratio = 1.0;
         goto done;
      }
      /* either keep previous choice or accept new choice depending on
       * which one is better */
      if (best_ratio < fabs(trow[q]) / big)
      {  csa->p = p;
#if 0 /* 29/III-2016 */
         memcpy(&csa->trow[1], &trow[1], (n-m) * sizeof(double));
#else
         memcpy(&csa->trow.vec[1], &trow[1], (n-m) * sizeof(double));
         fvs_gather_vec(&csa->trow, DBL_EPSILON);
#endif
         csa->q = q;
         best_ratio = fabs(trow[q]) / big;
      }
      /* check if the current choice is acceptable */
      if (best_ratio >= MIN_RATIO || nnn == 1 || try == 5)
         goto done;
      /* try to choose other xB[p] and xN[q] */
      /* find xB[p] in the list */
      for (t = 1; t <= nnn; t++)
         if (list[t] == p) break;
      xassert(t <= nnn);
      /* move xB[p] to the end of the list */
      list[t] = list[nnn], list[nnn] = p;
      /* and exclude it from consideration */
      nnn--;
      /* repeat the choice */
      goto try;
done: /* the choice has been made */
#if 1 /* FIXME: currently just to avoid badly conditioned basis */
      if (best_ratio < .001 * MIN_RATIO)
      {  /* looks like this helps */
         if (bfd_get_count(lp->bfd) > 0)
            return -1;
         /* didn't help; last chance to improve the choice */
         if (tol_piv == csa->tol_piv)
         {  tol_piv *= 1000.;
            goto more;
         }
      }
#endif
#if 1 /* FIXME */
      if (ret)
      {  /* invalidate basic solution components */
#if 0 /* 28/III-2016 */
         csa->beta_st = csa->d_st = 0;
#else
         /* dual solution remains valid */
         csa->beta_st = 0;
#endif
         /* set double-bounded non-basic variables to opposite bounds
          * for all break-points preceding the chosen one */
         for (t = 1; t < t_best; t++)
         {  k = head[m + bp[t].j];
            xassert(-DBL_MAX < l[k] && l[k] < u[k] && u[k] < +DBL_MAX);
            lp->flag[bp[t].j] = !(lp->flag[bp[t].j]);
         }
      }
#endif
      return ret;
}

#endif

/***********************************************************************
*  play_coef - play objective coefficients
*
*  This routine is called after the reduced costs d[j] was updated and
*  the basis was changed to the adjacent one.
*
*  It is assumed that before updating all the reduced costs d[j] were
*  strongly feasible, so in the adjacent basis d[j] remain feasible
*  within a tolerance, i.e. if some d[j] violates its zero bound, the
*  violation is insignificant.
*
*  If some d[j] violates its zero bound, the routine changes (perturbs)
*  objective coefficient cN[j] to provide d[j] = 0, i.e. to make all
*  d[j] strongly feasible. Otherwise, if d[j] has a feasible value, the
*  routine attempts to reduce (or remove) perturbation in cN[j] by
*  shifting d[j] to its zero bound keeping strong feasibility. */

static void play_coef(struct csa *csa, int all)
{     SPXLP *lp = csa->lp;
      int m = lp->m;
      int n = lp->n;
      double *c = lp->c;
      double *l = lp->l;
      double *u = lp->u;
      int *head = lp->head;
      char *flag = lp->flag;
      double *orig_c = csa->orig_c;
      double *d = csa->d;
      const double *trow = csa->trow.vec;
      /* this vector was used to update d = (d[j]) */
      int j, k;
      static const double eps = 1e-9;
      /* reduced costs d = (d[j]) should be valid */
      xassert(csa->d_st);
      /* walk thru the list of non-basic variables xN = (xN[j]) */
      for (j = 1; j <= n-m; j++)
      {  if (all || trow[j] != 0.0)
         {  /* d[j] has changed in the adjacent basis */
            k = head[m+j]; /* x[k] = xN[j] */
            if (l[k] == u[k])
            {  /* xN[j] is fixed variable */
               /* d[j] may have any sign */
            }
            else if (l[k] == -DBL_MAX && u[k] == +DBL_MAX)
            {  /* xN[j] is free (unbounded) variable */
               /* strong feasibility means d[j] = 0 */
               c[k] -= d[j], d[j] = 0.0;
               /* in this case dual degeneracy is not critical, since
                * if xN[j] enters the basis, it never leaves it */
            }
            else if (!flag[j])
            {  /* xN[j] has its lower bound active */
               xassert(l[k] != -DBL_MAX);
               /* first, we remove current perturbation to provide
                * c[k] = orig_c[k] */
               d[j] -= c[k] - orig_c[k], c[k] = orig_c[k];
               /* strong feasibility means d[j] >= 0, but we provide
                * d[j] >= +eps to prevent dual degeneracy */
               if (d[j] < +eps)
                  c[k] -= d[j] - eps, d[j] = +eps;
            }
            else
            {  /* xN[j] has its upper bound active */
               xassert(u[k] != +DBL_MAX);
               /* similarly, we remove current perturbation to provide
                * c[k] = orig_c[k] */
               d[j] -= c[k] - orig_c[k], c[k] = orig_c[k];
               /* strong feasibility means d[j] <= 0, but we provide
                * d[j] <= -eps to prevent dual degeneracy */
               if (d[j] > -eps)
                  c[k] -= d[j] + eps, d[j] = -eps;
            }
         }
      }
      return;
}

#if 1 /* 11/VII-2017 */
static void remove_perturb(struct csa *csa)
{     /* remove perturbation */
      SPXLP *lp = csa->lp;
      int n = lp->n;
      double *c = lp->c;
      double *orig_c = csa->orig_c;
      memcpy(c, orig_c, (1+n) * sizeof(double));
      /* removing perturbation changes dual solution components */
      csa->phase = csa->d_st = 0;
#if 1
      if (csa->msg_lev >= GLP_MSG_ALL)
         xprintf("Removing LP perturbation [%d]...\n",
            csa->it_cnt);
#endif
      return;
}
#endif

/***********************************************************************
*  display - display search progress
*
*  This routine displays some information about the search progress
*  that includes:
*
*  search phase;
*
*  number of simplex iterations performed by the solver;
*
*  original objective value (only on phase II);
*
*  sum of (scaled) dual infeasibilities for original bounds;
*
*  number of dual infeasibilities (phase I) or primal infeasibilities
*  (phase II);
*
*  number of basic factorizations since last display output. */

static void display(struct csa *csa, int spec)
{     SPXLP *lp = csa->lp;
      int m = lp->m;
      int n = lp->n;
      int *head = lp->head;
      char *flag = lp->flag;
      double *l = csa->orig_l; /* original lower bounds */
      double *u = csa->orig_u; /* original upper bounds */
      double *beta = csa->beta;
      double *d = csa->d;
      int j, k, nnn;
      double sum;
#if 1 /* 15/VII-2017 */
      double tm_cur;
#endif
      /* check if the display output should be skipped */
      if (csa->msg_lev < GLP_MSG_ON) goto skip;
#if 1 /* 15/VII-2017 */
      tm_cur = xtime();
#endif
      if (csa->out_dly > 0 &&
#if 0 /* 15/VII-2017 */
         1000.0 * xdifftime(xtime(), csa->tm_beg) < csa->out_dly)
#else
         1000.0 * xdifftime(tm_cur, csa->tm_beg) < csa->out_dly)
#endif
         goto skip;
      if (csa->it_cnt == csa->it_dpy) goto skip;
#if 0 /* 15/VII-2017 */
      if (!spec && csa->it_cnt % csa->out_frq != 0) goto skip;
#else
      if (!spec &&
         1000.0 * xdifftime(tm_cur, csa->tm_dpy) < csa->out_frq)
         goto skip;
#endif
      /* display search progress depending on search phase */
      switch (csa->phase)
      {  case 1:
            /* compute sum and number of (scaled) dual infeasibilities
             * for original bounds */
            sum = 0.0, nnn = 0;
            for (j = 1; j <= n-m; j++)
            {  k = head[m+j]; /* x[k] = xN[j] */
               if (d[j] > 0.0)
               {  /* xN[j] should have lower bound */
                  if (l[k] == -DBL_MAX)
                  {  sum += d[j];
                     if (d[j] > +1e-7)
                        nnn++;
                  }
               }
               else if (d[j] < 0.0)
               {  /* xN[j] should have upper bound */
                  if (u[k] == +DBL_MAX)
                  {  sum -= d[j];
                     if (d[j] < -1e-7)
                        nnn++;
                  }
               }
            }
            /* on phase I variables have artificial bounds which are
             * meaningless for original LP, so corresponding objective
             * function value is also meaningless */
#if 0 /* 27/III-2016 */
            xprintf(" %6d: %23s inf = %11.3e (%d)",
               csa->it_cnt, "", sum, nnn);
#else
            xprintf(" %6d: sum = %17.9e inf = %11.3e (%d)",
               csa->it_cnt, lp->c[0] - spx_eval_obj(lp, beta),
               sum, nnn);
#endif
            break;
         case 2:
            /* compute sum of (scaled) dual infeasibilities */
            sum = 0.0, nnn = 0;
            for (j = 1; j <= n-m; j++)
            {  k = head[m+j]; /* x[k] = xN[j] */
               if (d[j] > 0.0)
               {  /* xN[j] should have its lower bound active */
                  if (l[k] == -DBL_MAX || flag[j])
                     sum += d[j];
               }
               else if (d[j] < 0.0)
               {  /* xN[j] should have its upper bound active */
                  if (l[k] != u[k] && !flag[j])
                     sum -= d[j];
               }
            }
            /* compute number of primal infeasibilities */
            nnn = spy_chuzr_sel(lp, beta, csa->tol_bnd, csa->tol_bnd1,
               NULL);
            xprintf("#%6d: obj = %17.9e inf = %11.3e (%d)",
#if SCALE_Z
               csa->it_cnt,
               (double)csa->dir * csa->fz * spx_eval_obj(lp, beta),
#else
               csa->it_cnt, (double)csa->dir * spx_eval_obj(lp, beta),
#endif
               sum, nnn);
            break;
         default:
            xassert(csa != csa);
      }
      if (csa->inv_cnt)
      {  /* number of basis factorizations performed */
         xprintf(" %d", csa->inv_cnt);
         csa->inv_cnt = 0;
      }
#if 1 /* 23/III-2016 */
      if (csa->r_test == GLP_RT_FLIP)
      {  /*xprintf("   %d,%d", csa->ns_cnt, csa->ls_cnt);*/
         if (csa->ns_cnt + csa->ls_cnt)
            xprintf(" %d%%",
               (100 * csa->ls_cnt) / (csa->ns_cnt + csa->ls_cnt));
         csa->ns_cnt = csa->ls_cnt = 0;
      }
#endif
      xprintf("\n");
      csa->it_dpy = csa->it_cnt;
#if 1 /* 15/VII-2017 */
      csa->tm_dpy = tm_cur;
#endif
skip: return;
}

#if 1 /* 31/III-2016 */
static
void spy_update_r(SPXLP *lp, int p, int q, const double beta[/*1+m*/],
      const FVS *tcol, double tol, double tol1, FVS *r)
{     /* update vector r of primal infeasibilities */
      /* it is assumed that xB[p] leaves the basis, xN[q] enters the
       * basis, and beta corresponds to the adjacent basis (i.e. this
       * routine should be called after spx_update_beta) */
      int m = lp->m;
      int n = lp->n;
      double *l = lp->l;
      double *u = lp->u;
      int *head = lp->head;
      int *tcol_ind = tcol->ind;
      int *ind = r->ind;
      double *vec = r->vec;
      int i, k, t, nnz;
      double lk, uk, ri, eps;
      xassert(1 <= p && p <= m);
      xassert(1 <= q && q <= n-m);
      nnz = r->nnz;
      for (t = tcol->nnz; t >= 1; t--)
      {  i = tcol_ind[t];
         /* xB[i] changes in the adjacent basis to beta[i], so only
          * r[i] should be updated */
         if (i == p)
            k = head[m+q]; /* x[k] = new xB[p] = old xN[q] */
         else
            k = head[i];   /* x[k] = new xB[i] = old xB[i] */
         lk = l[k], uk = u[k];
         /* determine new value of r[i]; see spy_eval_r */
         ri = 0.0;
         if (beta[i] < lk)
         {  /* determine absolute tolerance eps1[i] */
            eps = tol + tol1 * (lk >= 0.0 ? +lk : -lk);
            if (beta[i] < lk - eps)
            {  /* lower bound is violated */
               ri = lk - beta[i];
            }
         }
         else if (beta[i] > uk)
         {  /* determine absolute tolerance eps2[i] */
            eps = tol + tol1 * (uk >= 0.0 ? +uk : -uk);
            if (beta[i] > uk + eps)
            {  /* upper bound is violated */
               ri = uk - beta[i];
            }
         }
         if (ri == 0.0)
         {  if (vec[i] != 0.0)
               vec[i] = DBL_MIN; /* will be removed */
         }
         else
         {  if (vec[i] == 0.0)
               ind[++nnz] = i;
            vec[i] = ri;
         }

      }
      r->nnz = nnz;
      /* remove zero elements */
      fvs_adjust_vec(r, DBL_MIN + DBL_MIN);
      return;
}
#endif

/***********************************************************************
*  spy_dual - driver to the dual simplex method
*
*  This routine is a driver to the two-phase dual simplex method.
*
*  On exit this routine returns one of the following codes:
*
*  0  LP instance has been successfully solved.
*
*  GLP_EOBJLL
*     Objective lower limit has been reached (maximization).
*
*  GLP_EOBJUL
*     Objective upper limit has been reached (minimization).
*
*  GLP_EITLIM
*     Iteration limit has been exhausted.
*
*  GLP_ETMLIM
*     Time limit has been exhausted.
*
*  GLP_EFAIL
*     The solver failed to solve LP instance. */

static int dual_simplex(struct csa *csa)
{     /* dual simplex method main logic routine */
      SPXLP *lp = csa->lp;
      int m = lp->m;
      int n = lp->n;
      double *l = lp->l;
      double *u = lp->u;
      int *head = lp->head;
      SPXNT *nt = csa->nt;
      double *beta = csa->beta;
      double *d = csa->d;
      SPYSE *se = csa->se;
#if 0 /* 30/III-2016 */
      int *list = csa->list;
#endif
#if 0 /* 31/III-2016 */
      double *trow = csa->trow;
      double *tcol = csa->tcol;
#endif
      double *pi = csa->work;
      int msg_lev = csa->msg_lev;
      double tol_bnd = csa->tol_bnd;
      double tol_bnd1 = csa->tol_bnd1;
      double tol_dj = csa->tol_dj;
      double tol_dj1 = csa->tol_dj1;
      int j, k, p_flag, refct, ret;
      int perturb = -1;
      /* -1 = perturbation is not used, but enabled
       *  0 = perturbation is not used and disabled
       * +1 = perturbation is being used */
#if 1 /* 27/III-2016 */
      int instab = 0; /* instability count */
#endif
#ifdef TIMING
      double t_total  = timer(); /* total time */
      double t_fact   = 0.0;     /* computing factorization */
      double t_rtest  = 0.0;     /* performing ratio test */
      double t_pivcol = 0.0;     /* computing pivot column */
      double t_upd1   = 0.0;     /* updating primal values */
      double t_upd2   = 0.0;     /* updating dual values */
      double t_upd3   = 0.0;     /* updating se weights */
      double t_upd4   = 0.0;     /* updating matrix N */
      double t_upd5   = 0.0;     /* updating factorization */
      double t_start;
#endif
      check_flags(csa);
loop: /* main loop starts here */
      /* compute factorization of the basis matrix */
      if (!lp->valid)
      {  double cond;
#ifdef TIMING
         t_start = timer();
#endif
         ret = spx_factorize(lp);
#ifdef TIMING
         t_fact += timer() - t_start;
#endif
         csa->inv_cnt++;
         if (ret != 0)
         {  if (msg_lev >= GLP_MSG_ERR)
               xprintf("Error: unable to factorize the basis matrix (%d"
                  ")\n", ret);
            csa->p_stat = csa->d_stat = GLP_UNDEF;
            ret = GLP_EFAIL;
            goto fini;
         }
         /* check condition of the basis matrix */
         cond = bfd_condest(lp->bfd);
         if (cond > 1.0 / DBL_EPSILON)
         {  if (msg_lev >= GLP_MSG_ERR)
               xprintf("Error: basis matrix is singular to working prec"
                  "ision (cond = %.3g)\n", cond);
            csa->p_stat = csa->d_stat = GLP_UNDEF;
            ret = GLP_EFAIL;
            goto fini;
         }
         if (cond > 0.001 / DBL_EPSILON)
         {  if (msg_lev >= GLP_MSG_ERR)
               xprintf("Warning: basis matrix is ill-conditioned (cond "
                  "= %.3g)\n", cond);
         }
         /* invalidate basic solution components */
         csa->beta_st = csa->d_st = 0;
      }
      /* compute reduced costs of non-basic variables d = (d[j]) */
      if (!csa->d_st)
      {  spx_eval_pi(lp, pi);
         for (j = 1; j <= n-m; j++)
            d[j] = spx_eval_dj(lp, pi, j);
         csa->d_st = 1; /* just computed */
         /* determine the search phase, if not determined yet (this is
          * performed only once at the beginning of the search for the
          * original bounds) */
         if (!csa->phase)
         {  j = check_feas(csa, 0.97 * tol_dj, 0.97 * tol_dj1, 1);
            if (j > 0)
            {  /* initial basic solution is dual infeasible and cannot
                * be recovered */
               /* start to search for dual feasible solution */
               set_art_bounds(csa);
               csa->phase = 1;
            }
            else
            {  /* initial basic solution is either dual feasible or its
                * dual feasibility has been recovered */
               /* start to search for optimal solution */
               csa->phase = 2;
            }
         }
         /* make sure that current basic solution is dual feasible */
#if 1 /* 11/VII-2017 */
         if (perturb <= 0)
         {  if (check_feas(csa, tol_dj, tol_dj1, 0))
            {  /* dual feasibility is broken due to excessive round-off
                * errors */
               if (perturb < 0)
               {  if (msg_lev >= GLP_MSG_ALL)
                     xprintf("Perturbing LP to avoid instability [%d].."
                        ".\n", csa->it_cnt);
                  perturb = 1;
                  goto loop;
               }
               if (msg_lev >= GLP_MSG_ERR)
                  xprintf("Warning: numerical instability (dual simplex"
                     ", phase %s)\n", csa->phase == 1 ? "I" : "II");
               instab++;
               if (csa->dualp && instab >= 10)
               {  /* do not continue the search; report failure */
                  if (msg_lev >= GLP_MSG_ERR)
                     xprintf("Warning: dual simplex failed due to exces"
                        "sive numerical instability\n");
                  csa->p_stat = csa->d_stat = GLP_UNDEF;
                  ret = -1; /* special case of GLP_EFAIL */
                  goto fini;
               }
               /* try to recover dual feasibility */
               j = check_feas(csa, 0.97 * tol_dj, 0.97 * tol_dj1, 1);
               if (j > 0)
               {  /* dual feasibility cannot be recovered (this may
                   * happen only on phase II) */
                  xassert(csa->phase == 2);
                  /* restart to search for dual feasible solution */
                  set_art_bounds(csa);
                  csa->phase = 1;
               }
            }
         }
         else
         {  /* FIXME */
            play_coef(csa, 1);
         }
      }
#endif
      /* at this point the search phase is determined */
      xassert(csa->phase == 1 || csa->phase == 2);
      /* compute values of basic variables beta = (beta[i]) */
      if (!csa->beta_st)
      {  spx_eval_beta(lp, beta);
#if 1 /* 31/III-2016 */
         /* also compute vector r of primal infeasibilities */
         switch (csa->phase)
         {  case 1:
               spy_eval_r(lp, beta, 1e-8, 0.0, &csa->r);
               break;
            case 2:
               spy_eval_r(lp, beta, tol_bnd, tol_bnd1, &csa->r);
               break;
            default:
               xassert(csa != csa);
         }
#endif
         csa->beta_st = 1; /* just computed */
      }
      /* reset the dual reference space, if necessary */
      if (se != NULL && !se->valid)
         spy_reset_refsp(lp, se), refct = 1000;
      /* at this point the basis factorization and all basic solution
       * components are valid */
      xassert(lp->valid && csa->beta_st && csa->d_st);
#ifdef GLP_DEBUG
      check_flags(csa);
#endif
#if CHECK_ACCURACY
      /* check accuracy of current basic solution components (only for
       * debugging) */
      check_accuracy(csa);
#endif
      /* check if the objective limit has been reached */
      if (csa->phase == 2 && csa->obj_lim != DBL_MAX
         && spx_eval_obj(lp, beta) >= csa->obj_lim)
      {
#if 1 /* 26/V-2017 by mao */
         if (perturb > 0)
         {  /* remove perturbation */
            /* [Should note that perturbing of objective coefficients
             * implemented in play_coef is equivalent to *relaxing* of
             * (zero) bounds of dual variables, so the perturbed
             * objective is always better (*greater*) that the original
             * one at the same basic point.] */
            remove_perturb(csa);
            perturb = 0;
         }
#endif
         if (csa->beta_st != 1)
            csa->beta_st = 0;
         if (csa->d_st != 1)
            csa->d_st = 0;
         if (!(csa->beta_st && csa->d_st))
            goto loop;
         display(csa, 1);
         if (msg_lev >= GLP_MSG_ALL)
            xprintf("OBJECTIVE %s LIMIT REACHED; SEARCH TERMINATED\n",
               csa->dir > 0 ? "UPPER" : "LOWER");
#if 0 /* 30/III-2016 */
         csa->num = spy_chuzr_sel(lp, beta, tol_bnd, tol_bnd1, list);
         csa->p_stat = (csa->num == 0 ? GLP_FEAS : GLP_INFEAS);
#else
         spy_eval_r(lp, beta, tol_bnd, tol_bnd1, &csa->r);
         csa->p_stat = (csa->r.nnz == 0 ? GLP_FEAS : GLP_INFEAS);
#endif
         csa->d_stat = GLP_FEAS;
         ret = (csa->dir > 0 ? GLP_EOBJUL : GLP_EOBJLL);
         goto fini;
      }
      /* check if the iteration limit has been exhausted */
      if (csa->it_cnt - csa->it_beg >= csa->it_lim)
      {  if (perturb > 0)
         {  /* remove perturbation */
            remove_perturb(csa);
            perturb = 0;
         }
         if (csa->beta_st != 1)
            csa->beta_st = 0;
         if (csa->d_st != 1)
            csa->d_st = 0;
         if (!(csa->beta_st && csa->d_st))
            goto loop;
         display(csa, 1);
         if (msg_lev >= GLP_MSG_ALL)
            xprintf("ITERATION LIMIT EXCEEDED; SEARCH TERMINATED\n");
         if (csa->phase == 1)
         {  set_orig_bounds(csa);
            check_flags(csa);
            spx_eval_beta(lp, beta);
         }
#if 0 /* 30/III-2016 */
         csa->num = spy_chuzr_sel(lp, beta, tol_bnd, tol_bnd1, list);
         csa->p_stat = (csa->num == 0 ? GLP_FEAS : GLP_INFEAS);
#else
         spy_eval_r(lp, beta, tol_bnd, tol_bnd1, &csa->r);
         csa->p_stat = (csa->r.nnz == 0 ? GLP_FEAS : GLP_INFEAS);
#endif
         csa->d_stat = (csa->phase == 1 ? GLP_INFEAS : GLP_FEAS);
         ret = GLP_EITLIM;
         goto fini;
      }
      /* check if the time limit has been exhausted */
      if (1000.0 * xdifftime(xtime(), csa->tm_beg) >= csa->tm_lim)
      {  if (perturb > 0)
         {  /* remove perturbation */
            remove_perturb(csa);
            perturb = 0;
         }
         if (csa->beta_st != 1)
            csa->beta_st = 0;
         if (csa->d_st != 1)
            csa->d_st = 0;
         if (!(csa->beta_st && csa->d_st))
            goto loop;
         display(csa, 1);
         if (msg_lev >= GLP_MSG_ALL)
            xprintf("TIME LIMIT EXCEEDED; SEARCH TERMINATED\n");
         if (csa->phase == 1)
         {  set_orig_bounds(csa);
            check_flags(csa);
            spx_eval_beta(lp, beta);
         }
#if 0 /* 30/III-2016 */
         csa->num = spy_chuzr_sel(lp, beta, tol_bnd, tol_bnd1, list);
         csa->p_stat = (csa->num == 0 ? GLP_FEAS : GLP_INFEAS);
#else
         spy_eval_r(lp, beta, tol_bnd, tol_bnd1, &csa->r);
         csa->p_stat = (csa->r.nnz == 0 ? GLP_FEAS : GLP_INFEAS);
#endif
         csa->d_stat = (csa->phase == 1 ? GLP_INFEAS : GLP_FEAS);
         ret = GLP_ETMLIM;
         goto fini;
      }
      /* display the search progress */
      display(csa, 0);
      /* select eligible basic variables */
#if 0 /* 31/III-2016; not needed because r is valid */
      switch (csa->phase)
      {  case 1:
#if 0 /* 30/III-2016 */
            csa->num = spy_chuzr_sel(lp, beta, 1e-8, 0.0, list);
#else
            spy_eval_r(lp, beta, 1e-8, 0.0, &csa->r);
#endif
            break;
         case 2:
#if 0 /* 30/III-2016 */
            csa->num = spy_chuzr_sel(lp, beta, tol_bnd, tol_bnd1, list);
#else
            spy_eval_r(lp, beta, tol_bnd, tol_bnd1, &csa->r);
#endif
            break;
         default:
            xassert(csa != csa);
      }
#endif
      /* check for optimality */
#if 0 /* 30/III-2016 */
      if (csa->num == 0)
#else
      if (csa->r.nnz == 0)
#endif
      {  if (perturb > 0 && csa->phase == 2)
         {  /* remove perturbation */
            remove_perturb(csa);
            perturb = 0;
         }
         if (csa->beta_st != 1)
            csa->beta_st = 0;
         if (csa->d_st != 1)
            csa->d_st = 0;
         if (!(csa->beta_st && csa->d_st))
            goto loop;
         /* current basis is optimal */
         display(csa, 1);
         switch (csa->phase)
         {  case 1:
               /* check for dual feasibility */
               set_orig_bounds(csa);
               check_flags(csa);
               if (check_feas(csa, tol_dj, tol_dj1, 0) == 0)
               {  /* dual feasible solution found; switch to phase II */
                  csa->phase = 2;
                  xassert(!csa->beta_st);
                  goto loop;
               }
#if 1 /* 26/V-2017 by cmatraki */
               if (perturb > 0)
               {  /* remove perturbation */
                  remove_perturb(csa);
                  perturb = 0;
                  goto loop;
               }
#endif
               /* no dual feasible solution exists */
               if (msg_lev >= GLP_MSG_ALL)
                  xprintf("LP HAS NO DUAL FEASIBLE SOLUTION\n");
               spx_eval_beta(lp, beta);
#if 0 /* 30/III-2016 */
               csa->num = spy_chuzr_sel(lp, beta, tol_bnd, tol_bnd1,
                  list);
               csa->p_stat = (csa->num == 0 ? GLP_FEAS : GLP_INFEAS);
#else
               spy_eval_r(lp, beta, tol_bnd, tol_bnd1, &csa->r);
               csa->p_stat = (csa->r.nnz == 0 ? GLP_FEAS : GLP_INFEAS);
#endif
               csa->d_stat = GLP_NOFEAS;
               ret = 0;
               goto fini;
            case 2:
               /* optimal solution found */
               if (msg_lev >= GLP_MSG_ALL)
                  xprintf("OPTIMAL LP SOLUTION FOUND\n");
               csa->p_stat = csa->d_stat = GLP_FEAS;
               ret = 0;
               goto fini;
            default:
               xassert(csa != csa);
         }
      }
      /* choose xB[p] and xN[q] */
#if 0 /* 23/III-2016 */
      choose_pivot(csa);
#else
#ifdef TIMING
      t_start = timer();
#endif
#if 1 /* 31/III-2016 */
      ret = choose_pivot(csa);
#endif
#ifdef TIMING
      t_rtest += timer() - t_start;
#endif
      if (ret < 0)
      {  lp->valid = 0;
         goto loop;
      }
      if (ret == 0)
         csa->ns_cnt++;
      else
         csa->ls_cnt++;
#endif
      /* check for dual unboundedness */
      if (csa->q == 0)
      {  if (perturb > 0)
         {  /* remove perturbation */
            remove_perturb(csa);
            perturb = 0;
         }
         if (csa->beta_st != 1)
            csa->beta_st = 0;
         if (csa->d_st != 1)
            csa->d_st = 0;
         if (!(csa->beta_st && csa->d_st))
            goto loop;
         display(csa, 1);
         switch (csa->phase)
         {  case 1:
               /* this should never happen */
               if (msg_lev >= GLP_MSG_ERR)
                  xprintf("Error: dual simplex failed\n");
               csa->p_stat = csa->d_stat = GLP_UNDEF;
               ret = GLP_EFAIL;
               goto fini;
            case 2:
               /* dual unboundedness detected */
               if (msg_lev >= GLP_MSG_ALL)
                  xprintf("LP HAS NO PRIMAL FEASIBLE SOLUTION\n");
               csa->p_stat = GLP_NOFEAS;
               csa->d_stat = GLP_FEAS;
               ret = 0;
               goto fini;
            default:
               xassert(csa != csa);
         }
      }
      /* compute q-th column of the simplex table */
#ifdef TIMING
      t_start = timer();
#endif
#if 0 /* 31/III-2016 */
      spx_eval_tcol(lp, csa->q, tcol);
#else
      spx_eval_tcol(lp, csa->q, csa->tcol.vec);
      fvs_gather_vec(&csa->tcol, DBL_EPSILON);
#endif
#ifdef TIMING
      t_pivcol += timer() - t_start;
#endif
      /* FIXME: tcol[p] and trow[q] should be close to each other */
#if 0 /* 26/V-2017 by cmatraki */
      xassert(csa->tcol.vec[csa->p] != 0.0);
#else
      if (csa->tcol.vec[csa->p] == 0.0)
      {  if (msg_lev >= GLP_MSG_ERR)
            xprintf("Error: tcol[p] = 0.0\n");
         csa->p_stat = csa->d_stat = GLP_UNDEF;
         ret = GLP_EFAIL;
         goto fini;
      }
#endif
      /* update values of basic variables for adjacent basis */
      k = head[csa->p]; /* x[k] = xB[p] */
      p_flag = (l[k] != u[k] && beta[csa->p] > u[k]);
#if 0 /* 16/III-2016 */
      spx_update_beta(lp, beta, csa->p, p_flag, csa->q, tcol);
      csa->beta_st = 2;
#else
      /* primal solution may be invalidated due to long step */
#ifdef TIMING
      t_start = timer();
#endif
      if (csa->beta_st)
#if 0 /* 30/III-2016 */
      {  spx_update_beta(lp, beta, csa->p, p_flag, csa->q, tcol);
#else
      {  spx_update_beta_s(lp, beta, csa->p, p_flag, csa->q,
            &csa->tcol);
         /* also update vector r of primal infeasibilities */
         /*fvs_check_vec(&csa->r);*/
         switch (csa->phase)
         {  case 1:
               spy_update_r(lp, csa->p, csa->q, beta, &csa->tcol,
                  1e-8, 0.0, &csa->r);
               break;
            case 2:
               spy_update_r(lp, csa->p, csa->q, beta, &csa->tcol,
                  tol_bnd, tol_bnd1, &csa->r);
               break;
            default:
               xassert(csa != csa);
         }
         /*fvs_check_vec(&csa->r);*/
#endif
         csa->beta_st = 2;
      }
#ifdef TIMING
      t_upd1 += timer() - t_start;
#endif
#endif
#if 1 /* 11/VII-2017 */
      /* check for stalling */
      {  int k;
         xassert(1 <= csa->p && csa->p <= m);
         xassert(1 <= csa->q && csa->q <= n-m);
         /* FIXME: recompute d[q]; see spx_update_d */
         k = head[m+csa->q]; /* x[k] = xN[q] */
         if (!(lp->l[k] == -DBL_MAX && lp->u[k] == +DBL_MAX))
         {  if (fabs(d[csa->q]) >= 1e-6)
            {  csa->degen = 0;
               goto skip1;
            }
            /* degenerate iteration has been detected */
            csa->degen++;
            if (perturb < 0 && csa->degen >= 200)
            {  if (msg_lev >= GLP_MSG_ALL)
                  xprintf("Perturbing LP to avoid stalling [%d]...\n",
                     csa->it_cnt);
               perturb = 1;
            }
skip1:      ;
         }
      }
#endif
      /* update reduced costs of non-basic variables for adjacent
       * basis */
#if 1 /* 28/III-2016 */
      xassert(csa->d_st);
#endif
#ifdef TIMING
      t_start = timer();
#endif
#if 0 /* 30/III-2016 */
      if (spx_update_d(lp, d, csa->p, csa->q, trow, tcol) <= 1e-9)
#else
      if (spx_update_d_s(lp, d, csa->p, csa->q, &csa->trow, &csa->tcol)
         <= 1e-9)
#endif
      {  /* successful updating */
         csa->d_st = 2;
      }
      else
      {  /* new reduced costs are inaccurate */
         csa->d_st = 0;
      }
#ifdef TIMING
      t_upd2 += timer() - t_start;
#endif
      /* update steepest edge weights for adjacent basis, if used */
#ifdef TIMING
      t_start = timer();
#endif
      if (se != NULL)
      {  if (refct > 0)
#if 0 /* 30/III-2016 */
         {  if (spy_update_gamma(lp, se, csa->p, csa->q, trow, tcol)
               <= 1e-3)
#else
         {  if (spy_update_gamma_s(lp, se, csa->p, csa->q, &csa->trow,
               &csa->tcol) <= 1e-3)
#endif
            {  /* successful updating */
               refct--;
            }
            else
            {  /* new weights are inaccurate; reset reference space */
               se->valid = 0;
            }
         }
         else
         {  /* too many updates; reset reference space */
            se->valid = 0;
         }
      }
#ifdef TIMING
      t_upd3 += timer() - t_start;
#endif
#ifdef TIMING
      t_start = timer();
#endif
      /* update matrix N for adjacent basis, if used */
      if (nt != NULL)
         spx_update_nt(lp, nt, csa->p, csa->q);
#ifdef TIMING
      t_upd4 += timer() - t_start;
#endif
      /* change current basis header to adjacent one */
      spx_change_basis(lp, csa->p, p_flag, csa->q);
      /* and update factorization of the basis matrix */
#ifdef TIMING
      t_start = timer();
#endif
#if 0 /* 16/III-2016 */
      if (csa->p > 0)
#endif
         spx_update_invb(lp, csa->p, head[csa->p]);
#ifdef TIMING
      t_upd5 += timer() - t_start;
#endif
      if (perturb > 0 && csa->d_st)
         play_coef(csa, 0);
      /* dual simplex iteration complete */
      csa->it_cnt++;
      goto loop;
fini:
#ifdef TIMING
      t_total = timer() - t_total;
      xprintf("Total time      = %10.3f\n", t_total);
      xprintf("Factorization   = %10.3f\n", t_fact);
      xprintf("Ratio test      = %10.3f\n", t_rtest);
      xprintf("Pivot column    = %10.3f\n", t_pivcol);
      xprintf("Updating beta   = %10.3f\n", t_upd1);
      xprintf("Updating d      = %10.3f\n", t_upd2);
      xprintf("Updating gamma  = %10.3f\n", t_upd3);
      xprintf("Updating N      = %10.3f\n", t_upd4);
      xprintf("Updating inv(B) = %10.3f\n", t_upd5);
#endif
      return ret;
}

int spy_dual(glp_prob *P, const glp_smcp *parm)
{     /* driver to the dual simplex method */
      struct csa csa_, *csa = &csa_;
      SPXLP lp;
      SPXAT at;
      SPXNT nt;
      SPYSE se;
      int ret, *map, *daeh;
#if SCALE_Z
      int i, j, k;
#endif
      /* build working LP and its initial basis */
      memset(csa, 0, sizeof(struct csa));
      csa->lp = &lp;
      spx_init_lp(csa->lp, P, parm->excl);
      spx_alloc_lp(csa->lp);
      map = talloc(1+P->m+P->n, int);
      spx_build_lp(csa->lp, P, parm->excl, parm->shift, map);
      spx_build_basis(csa->lp, P, map);
      switch (P->dir)
      {  case GLP_MIN:
            csa->dir = +1;
            break;
         case GLP_MAX:
            csa->dir = -1;
            break;
         default:
            xassert(P != P);
      }
#if SCALE_Z
      csa->fz = 0.0;
      for (k = 1; k <= csa->lp->n; k++)
      {  double t = fabs(csa->lp->c[k]);
         if (csa->fz < t)
            csa->fz = t;
      }
      if (csa->fz <= 1000.0)
         csa->fz = 1.0;
      else
         csa->fz /= 1000.0;
      /*xprintf("csa->fz = %g\n", csa->fz);*/
      for (k = 0; k <= csa->lp->n; k++)
         csa->lp->c[k] /= csa->fz;
#endif
      csa->orig_b = talloc(1+csa->lp->m, double);
      memcpy(csa->orig_b, csa->lp->b, (1+csa->lp->m) * sizeof(double));
      csa->orig_c = talloc(1+csa->lp->n, double);
      memcpy(csa->orig_c, csa->lp->c, (1+csa->lp->n) * sizeof(double));
      csa->orig_l = talloc(1+csa->lp->n, double);
      memcpy(csa->orig_l, csa->lp->l, (1+csa->lp->n) * sizeof(double));
      csa->orig_u = talloc(1+csa->lp->n, double);
      memcpy(csa->orig_u, csa->lp->u, (1+csa->lp->n) * sizeof(double));
      switch (parm->aorn)
      {  case GLP_USE_AT:
            /* build matrix A in row-wise format */
            csa->at = &at;
            csa->nt = NULL;
            spx_alloc_at(csa->lp, csa->at);
            spx_build_at(csa->lp, csa->at);
            break;
         case GLP_USE_NT:
            /* build matrix N in row-wise format for initial basis */
            csa->at = NULL;
            csa->nt = &nt;
            spx_alloc_nt(csa->lp, csa->nt);
            spx_init_nt(csa->lp, csa->nt);
            spx_build_nt(csa->lp, csa->nt);
            break;
         default:
            xassert(parm != parm);
      }
      /* allocate and initialize working components */
      csa->phase = 0;
      csa->beta = talloc(1+csa->lp->m, double);
      csa->beta_st = 0;
      csa->d = talloc(1+csa->lp->n-csa->lp->m, double);
      csa->d_st = 0;
      switch (parm->pricing)
      {  case GLP_PT_STD:
            csa->se = NULL;
            break;
         case GLP_PT_PSE:
            csa->se = &se;
            spy_alloc_se(csa->lp, csa->se);
            break;
         default:
            xassert(parm != parm);
      }
#if 0 /* 30/III-2016 */
      csa->list = talloc(1+csa->lp->m, int);
      csa->trow = talloc(1+csa->lp->n-csa->lp->m, double);
      csa->tcol = talloc(1+csa->lp->m, double);
#else
      fvs_alloc_vec(&csa->r, csa->lp->m);
      fvs_alloc_vec(&csa->trow, csa->lp->n-csa->lp->m);
      fvs_alloc_vec(&csa->tcol, csa->lp->m);
#endif
#if 1 /* 16/III-2016 */
      csa->bp = NULL;
#endif
      csa->work = talloc(1+csa->lp->m, double);
      csa->work1 = talloc(1+csa->lp->n-csa->lp->m, double);
#if 0 /* 11/VI-2017 */
#if 1 /* 31/III-2016 */
      fvs_alloc_vec(&csa->wrow, csa->lp->n-csa->lp->m);
      fvs_alloc_vec(&csa->wcol, csa->lp->m);
#endif
#endif
      /* initialize control parameters */
      csa->msg_lev = parm->msg_lev;
      csa->dualp = (parm->meth == GLP_DUALP);
#if 0 /* 16/III-2016 */
      switch (parm->r_test)
      {  case GLP_RT_STD:
            csa->harris = 0;
            break;
         case GLP_RT_HAR:
            csa->harris = 1;
            break;
         default:
            xassert(parm != parm);
      }
#else
      switch (parm->r_test)
      {  case GLP_RT_STD:
         case GLP_RT_HAR:
            break;
         case GLP_RT_FLIP:
            csa->bp = talloc(1+csa->lp->n-csa->lp->m, SPYBP);
            break;
         default:
            xassert(parm != parm);
      }
      csa->r_test = parm->r_test;
#endif
      csa->tol_bnd = parm->tol_bnd;
      csa->tol_bnd1 = .001 * parm->tol_bnd;
      csa->tol_dj = parm->tol_dj;
      csa->tol_dj1 = .001 * parm->tol_dj;
#if 0
      csa->tol_dj1 = 1e-9 * parm->tol_dj;
#endif
      csa->tol_piv = parm->tol_piv;
      switch (P->dir)
      {  case GLP_MIN:
            csa->obj_lim = + parm->obj_ul;
            break;
         case GLP_MAX:
            csa->obj_lim = - parm->obj_ll;
            break;
         default:
            xassert(parm != parm);
      }
#if SCALE_Z
      if (csa->obj_lim != DBL_MAX)
         csa->obj_lim /= csa->fz;
#endif
      csa->it_lim = parm->it_lim;
      csa->tm_lim = parm->tm_lim;
      csa->out_frq = parm->out_frq;
      csa->out_dly = parm->out_dly;
      /* initialize working parameters */
      csa->tm_beg = xtime();
      csa->it_beg = csa->it_cnt = P->it_cnt;
      csa->it_dpy = -1;
#if 1 /* 15/VII-2017 */
      csa->tm_dpy = 0.0;
#endif
      csa->inv_cnt = 0;
#if 1 /* 11/VII-2017 */
      csa->degen = 0;
#endif
#if 1 /* 23/III-2016 */
      csa->ns_cnt = csa->ls_cnt = 0;
#endif
      /* try to solve working LP */
      ret = dual_simplex(csa);
      /* return basis factorization back to problem object */
      P->valid = csa->lp->valid;
      P->bfd = csa->lp->bfd;
      /* set solution status */
      P->pbs_stat = csa->p_stat;
      P->dbs_stat = csa->d_stat;
      /* if the solver failed, do not store basis header and basic
       * solution components to problem object */
      if (ret == GLP_EFAIL)
         goto skip;
      /* convert working LP basis to original LP basis and store it to
       * problem object */
      daeh = talloc(1+csa->lp->n, int);
      spx_store_basis(csa->lp, P, map, daeh);
      /* compute simplex multipliers for final basic solution found by
       * the solver */
      spx_eval_pi(csa->lp, csa->work);
      /* convert working LP solution to original LP solution and store
       * it to problem object */
#if SCALE_Z
      for (i = 1; i <= csa->lp->m; i++)
         csa->work[i] *= csa->fz;
      for (j = 1; j <= csa->lp->n-csa->lp->m; j++)
         csa->d[j] *= csa->fz;
#endif
      spx_store_sol(csa->lp, P, parm->shift, map, daeh, csa->beta,
         csa->work, csa->d);
      tfree(daeh);
      /* save simplex iteration count */
      P->it_cnt = csa->it_cnt;
      /* report auxiliary/structural variable causing unboundedness */
      P->some = 0;
      if (csa->p_stat == GLP_NOFEAS && csa->d_stat == GLP_FEAS)
      {  int k, kk;
         /* xB[p] = x[k] causes dual unboundedness */
         xassert(1 <= csa->p && csa->p <= csa->lp->m);
         k = csa->lp->head[csa->p];
         xassert(1 <= k && k <= csa->lp->n);
         /* convert to number of original variable */
         for (kk = 1; kk <= P->m + P->n; kk++)
         {  if (abs(map[kk]) == k)
            {  P->some = kk;
               break;
            }
         }
         xassert(P->some != 0);
      }
skip: /* deallocate working objects and arrays */
      spx_free_lp(csa->lp);
      tfree(map);
      tfree(csa->orig_b);
      tfree(csa->orig_c);
      tfree(csa->orig_l);
      tfree(csa->orig_u);
      if (csa->at != NULL)
         spx_free_at(csa->lp, csa->at);
      if (csa->nt != NULL)
         spx_free_nt(csa->lp, csa->nt);
      tfree(csa->beta);
      tfree(csa->d);
      if (csa->se != NULL)
         spy_free_se(csa->lp, csa->se);
#if 0 /* 30/III-2016 */
      tfree(csa->list);
      tfree(csa->trow);
#else
      fvs_free_vec(&csa->r);
      fvs_free_vec(&csa->trow);
#endif
#if 1 /* 16/III-2016 */
      if (csa->bp != NULL)
         tfree(csa->bp);
#endif
#if 0 /* 29/III-2016 */
      tfree(csa->tcol);
#else
      fvs_free_vec(&csa->tcol);
#endif
      tfree(csa->work);
      tfree(csa->work1);
#if 0 /* 11/VI-2017 */
#if 1 /* 31/III-2016 */
      fvs_free_vec(&csa->wrow);
      fvs_free_vec(&csa->wcol);
#endif
#endif
      /* return to calling program */
      return ret >= 0 ? ret : GLP_EFAIL;
}

/* eof */
