/* glpapi07.c (exact simplex solver) */

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

#include "draft.h"
#include "glpssx.h"
#include "misc.h"
#include "prob.h"

/***********************************************************************
*  NAME
*
*  glp_exact - solve LP problem in exact arithmetic
*
*  SYNOPSIS
*
*  int glp_exact(glp_prob *lp, const glp_smcp *parm);
*
*  DESCRIPTION
*
*  The routine glp_exact is a tentative implementation of the primal
*  two-phase simplex method based on exact (rational) arithmetic. It is
*  similar to the routine glp_simplex, however, for all internal
*  computations it uses arithmetic of rational numbers, which is exact
*  in mathematical sense, i.e. free of round-off errors unlike floating
*  point arithmetic.
*
*  Note that the routine glp_exact uses inly two control parameters
*  passed in the structure glp_smcp, namely, it_lim and tm_lim.
*
*  RETURNS
*
*  0  The LP problem instance has been successfully solved. This code
*     does not necessarily mean that the solver has found optimal
*     solution. It only means that the solution process was successful.
*
*  GLP_EBADB
*     Unable to start the search, because the initial basis specified
*     in the problem object is invalid--the number of basic (auxiliary
*     and structural) variables is not the same as the number of rows in
*     the problem object.
*
*  GLP_ESING
*     Unable to start the search, because the basis matrix correspodning
*     to the initial basis is exactly singular.
*
*  GLP_EBOUND
*     Unable to start the search, because some double-bounded variables
*     have incorrect bounds.
*
*  GLP_EFAIL
*     The problem has no rows/columns.
*
*  GLP_EITLIM
*     The search was prematurely terminated, because the simplex
*     iteration limit has been exceeded.
*
*  GLP_ETMLIM
*     The search was prematurely terminated, because the time limit has
*     been exceeded. */

static void set_d_eps(mpq_t x, double val)
{     /* convert double val to rational x obtaining a more adequate
         fraction than provided by mpq_set_d due to allowing a small
         approximation error specified by a given relative tolerance;
         for example, mpq_set_d would give the following
         1/3 ~= 0.333333333333333314829616256247391... ->
             -> 6004799503160661/18014398509481984
         while this routine gives exactly 1/3 */
      int s, n, j;
      double f, p, q, eps = 1e-9;
      mpq_t temp;
      xassert(-DBL_MAX <= val && val <= +DBL_MAX);
#if 1 /* 30/VII-2008 */
      if (val == floor(val))
      {  /* if val is integral, do not approximate */
         mpq_set_d(x, val);
         goto done;
      }
#endif
      if (val > 0.0)
         s = +1;
      else if (val < 0.0)
         s = -1;
      else
      {  mpq_set_si(x, 0, 1);
         goto done;
      }
      f = frexp(fabs(val), &n);
      /* |val| = f * 2^n, where 0.5 <= f < 1.0 */
      fp2rat(f, 0.1 * eps, &p, &q);
      /* f ~= p / q, where p and q are integers */
      mpq_init(temp);
      mpq_set_d(x, p);
      mpq_set_d(temp, q);
      mpq_div(x, x, temp);
      mpq_set_si(temp, 1, 1);
      for (j = 1; j <= abs(n); j++)
         mpq_add(temp, temp, temp);
      if (n > 0)
         mpq_mul(x, x, temp);
      else if (n < 0)
         mpq_div(x, x, temp);
      mpq_clear(temp);
      if (s < 0) mpq_neg(x, x);
      /* check that the desired tolerance has been attained */
      xassert(fabs(val - mpq_get_d(x)) <= eps * (1.0 + fabs(val)));
done: return;
}

static void load_data(SSX *ssx, glp_prob *lp)
{     /* load LP problem data into simplex solver workspace */
      int m = ssx->m;
      int n = ssx->n;
      int nnz = ssx->A_ptr[n+1]-1;
      int j, k, type, loc, len, *ind;
      double lb, ub, coef, *val;
      xassert(lp->m == m);
      xassert(lp->n == n);
      xassert(lp->nnz == nnz);
      /* types and bounds of rows and columns */
      for (k = 1; k <= m+n; k++)
      {  if (k <= m)
         {  type = lp->row[k]->type;
            lb = lp->row[k]->lb;
            ub = lp->row[k]->ub;
         }
         else
         {  type = lp->col[k-m]->type;
            lb = lp->col[k-m]->lb;
            ub = lp->col[k-m]->ub;
         }
         switch (type)
         {  case GLP_FR: type = SSX_FR; break;
            case GLP_LO: type = SSX_LO; break;
            case GLP_UP: type = SSX_UP; break;
            case GLP_DB: type = SSX_DB; break;
            case GLP_FX: type = SSX_FX; break;
            default: xassert(type != type);
         }
         ssx->type[k] = type;
         set_d_eps(ssx->lb[k], lb);
         set_d_eps(ssx->ub[k], ub);
      }
      /* optimization direction */
      switch (lp->dir)
      {  case GLP_MIN: ssx->dir = SSX_MIN; break;
         case GLP_MAX: ssx->dir = SSX_MAX; break;
         default: xassert(lp != lp);
      }
      /* objective coefficients */
      for (k = 0; k <= m+n; k++)
      {  if (k == 0)
            coef = lp->c0;
         else if (k <= m)
            coef = 0.0;
         else
            coef = lp->col[k-m]->coef;
         set_d_eps(ssx->coef[k], coef);
      }
      /* constraint coefficients */
      ind = xcalloc(1+m, sizeof(int));
      val = xcalloc(1+m, sizeof(double));
      loc = 0;
      for (j = 1; j <= n; j++)
      {  ssx->A_ptr[j] = loc+1;
         len = glp_get_mat_col(lp, j, ind, val);
         for (k = 1; k <= len; k++)
         {  loc++;
            ssx->A_ind[loc] = ind[k];
            set_d_eps(ssx->A_val[loc], val[k]);
         }
      }
      xassert(loc == nnz);
      xfree(ind);
      xfree(val);
      return;
}

static int load_basis(SSX *ssx, glp_prob *lp)
{     /* load current LP basis into simplex solver workspace */
      int m = ssx->m;
      int n = ssx->n;
      int *type = ssx->type;
      int *stat = ssx->stat;
      int *Q_row = ssx->Q_row;
      int *Q_col = ssx->Q_col;
      int i, j, k;
      xassert(lp->m == m);
      xassert(lp->n == n);
      /* statuses of rows and columns */
      for (k = 1; k <= m+n; k++)
      {  if (k <= m)
            stat[k] = lp->row[k]->stat;
         else
            stat[k] = lp->col[k-m]->stat;
         switch (stat[k])
         {  case GLP_BS:
               stat[k] = SSX_BS;
               break;
            case GLP_NL:
               stat[k] = SSX_NL;
               xassert(type[k] == SSX_LO || type[k] == SSX_DB);
               break;
            case GLP_NU:
               stat[k] = SSX_NU;
               xassert(type[k] == SSX_UP || type[k] == SSX_DB);
               break;
            case GLP_NF:
               stat[k] = SSX_NF;
               xassert(type[k] == SSX_FR);
               break;
            case GLP_NS:
               stat[k] = SSX_NS;
               xassert(type[k] == SSX_FX);
               break;
            default:
               xassert(stat != stat);
         }
      }
      /* build permutation matix Q */
      i = j = 0;
      for (k = 1; k <= m+n; k++)
      {  if (stat[k] == SSX_BS)
         {  i++;
            if (i > m) return 1;
            Q_row[k] = i, Q_col[i] = k;
         }
         else
         {  j++;
            if (j > n) return 1;
            Q_row[k] = m+j, Q_col[m+j] = k;
         }
      }
      xassert(i == m && j == n);
      return 0;
}

int glp_exact(glp_prob *lp, const glp_smcp *parm)
{     glp_smcp _parm;
      SSX *ssx;
      int m = lp->m;
      int n = lp->n;
      int nnz = lp->nnz;
      int i, j, k, type, pst, dst, ret, stat;
      double lb, ub, prim, dual, sum;
      if (parm == NULL)
         parm = &_parm, glp_init_smcp((glp_smcp *)parm);
      /* check control parameters */
#if 1 /* 25/XI-2017 */
      switch (parm->msg_lev)
      {  case GLP_MSG_OFF:
         case GLP_MSG_ERR:
         case GLP_MSG_ON:
         case GLP_MSG_ALL:
         case GLP_MSG_DBG:
            break;
         default:
            xerror("glp_exact: msg_lev = %d; invalid parameter\n",
               parm->msg_lev);
      }
#endif
      if (parm->it_lim < 0)
         xerror("glp_exact: it_lim = %d; invalid parameter\n",
            parm->it_lim);
      if (parm->tm_lim < 0)
         xerror("glp_exact: tm_lim = %d; invalid parameter\n",
            parm->tm_lim);
      /* the problem must have at least one row and one column */
      if (!(m > 0 && n > 0))
#if 0 /* 25/XI-2017 */
      {  xprintf("glp_exact: problem has no rows/columns\n");
#else
      {  if (parm->msg_lev >= GLP_MSG_ERR)
            xprintf("glp_exact: problem has no rows/columns\n");
#endif
         return GLP_EFAIL;
      }
#if 1
      /* basic solution is currently undefined */
      lp->pbs_stat = lp->dbs_stat = GLP_UNDEF;
      lp->obj_val = 0.0;
      lp->some = 0;
#endif
      /* check that all double-bounded variables have correct bounds */
      for (k = 1; k <= m+n; k++)
      {  if (k <= m)
         {  type = lp->row[k]->type;
            lb = lp->row[k]->lb;
            ub = lp->row[k]->ub;
         }
         else
         {  type = lp->col[k-m]->type;
            lb = lp->col[k-m]->lb;
            ub = lp->col[k-m]->ub;
         }
         if (type == GLP_DB && lb >= ub)
#if 0 /* 25/XI-2017 */
         {  xprintf("glp_exact: %s %d has invalid bounds\n",
               k <= m ? "row" : "column", k <= m ? k : k-m);
#else
         {  if (parm->msg_lev >= GLP_MSG_ERR)
               xprintf("glp_exact: %s %d has invalid bounds\n",
                  k <= m ? "row" : "column", k <= m ? k : k-m);
#endif
            return GLP_EBOUND;
         }
      }
      /* create the simplex solver workspace */
#if 1 /* 25/XI-2017 */
      if (parm->msg_lev >= GLP_MSG_ALL)
      {
#endif
      xprintf("glp_exact: %d rows, %d columns, %d non-zeros\n",
         m, n, nnz);
#ifdef HAVE_GMP
      xprintf("GNU MP bignum library is being used\n");
#else
      xprintf("GLPK bignum module is being used\n");
      xprintf("(Consider installing GNU MP to attain a much better perf"
         "ormance.)\n");
#endif
#if 1 /* 25/XI-2017 */
      }
#endif
      ssx = ssx_create(m, n, nnz);
      /* load LP problem data into the workspace */
      load_data(ssx, lp);
      /* load current LP basis into the workspace */
      if (load_basis(ssx, lp))
#if 0 /* 25/XI-2017 */
      {  xprintf("glp_exact: initial LP basis is invalid\n");
#else
      {  if (parm->msg_lev >= GLP_MSG_ERR)
            xprintf("glp_exact: initial LP basis is invalid\n");
#endif
         ret = GLP_EBADB;
         goto done;
      }
#if 0
      /* inherit some control parameters from the LP object */
      ssx->it_lim = lpx_get_int_parm(lp, LPX_K_ITLIM);
      ssx->it_cnt = lpx_get_int_parm(lp, LPX_K_ITCNT);
      ssx->tm_lim = lpx_get_real_parm(lp, LPX_K_TMLIM);
#else
#if 1 /* 25/XI-2017 */
      ssx->msg_lev = parm->msg_lev;
#endif
      ssx->it_lim = parm->it_lim;
      ssx->it_cnt = lp->it_cnt;
      ssx->tm_lim = (double)parm->tm_lim / 1000.0;
#endif
      ssx->out_frq = 5.0;
      ssx->tm_beg = xtime();
#if 0 /* 10/VI-2013 */
      ssx->tm_lag = xlset(0);
#else
      ssx->tm_lag = 0.0;
#endif
      /* solve LP */
      ret = ssx_driver(ssx);
#if 0
      /* copy back some statistics to the LP object */
      lpx_set_int_parm(lp, LPX_K_ITLIM, ssx->it_lim);
      lpx_set_int_parm(lp, LPX_K_ITCNT, ssx->it_cnt);
      lpx_set_real_parm(lp, LPX_K_TMLIM, ssx->tm_lim);
#else
      lp->it_cnt = ssx->it_cnt;
#endif
      /* analyze the return code */
      switch (ret)
      {  case 0:
            /* optimal solution found */
            ret = 0;
            pst = dst = GLP_FEAS;
            break;
         case 1:
            /* problem has no feasible solution */
            ret = 0;
            pst = GLP_NOFEAS, dst = GLP_INFEAS;
            break;
         case 2:
            /* problem has unbounded solution */
            ret = 0;
            pst = GLP_FEAS, dst = GLP_NOFEAS;
#if 1
            xassert(1 <= ssx->q && ssx->q <= n);
            lp->some = ssx->Q_col[m + ssx->q];
            xassert(1 <= lp->some && lp->some <= m+n);
#endif
            break;
         case 3:
            /* iteration limit exceeded (phase I) */
            ret = GLP_EITLIM;
            pst = dst = GLP_INFEAS;
            break;
         case 4:
            /* iteration limit exceeded (phase II) */
            ret = GLP_EITLIM;
            pst = GLP_FEAS, dst = GLP_INFEAS;
            break;
         case 5:
            /* time limit exceeded (phase I) */
            ret = GLP_ETMLIM;
            pst = dst = GLP_INFEAS;
            break;
         case 6:
            /* time limit exceeded (phase II) */
            ret = GLP_ETMLIM;
            pst = GLP_FEAS, dst = GLP_INFEAS;
            break;
         case 7:
            /* initial basis matrix is singular */
            ret = GLP_ESING;
            goto done;
         default:
            xassert(ret != ret);
      }
      /* store final basic solution components into LP object */
      lp->pbs_stat = pst;
      lp->dbs_stat = dst;
      sum = lp->c0;
      for (k = 1; k <= m+n; k++)
      {  if (ssx->stat[k] == SSX_BS)
         {  i = ssx->Q_row[k]; /* x[k] = xB[i] */
            xassert(1 <= i && i <= m);
            stat = GLP_BS;
            prim = mpq_get_d(ssx->bbar[i]);
            dual = 0.0;
         }
         else
         {  j = ssx->Q_row[k] - m; /* x[k] = xN[j] */
            xassert(1 <= j && j <= n);
            switch (ssx->stat[k])
            {  case SSX_NF:
                  stat = GLP_NF;
                  prim = 0.0;
                  break;
               case SSX_NL:
                  stat = GLP_NL;
                  prim = mpq_get_d(ssx->lb[k]);
                  break;
               case SSX_NU:
                  stat = GLP_NU;
                  prim = mpq_get_d(ssx->ub[k]);
                  break;
               case SSX_NS:
                  stat = GLP_NS;
                  prim = mpq_get_d(ssx->lb[k]);
                  break;
               default:
                  xassert(ssx != ssx);
            }
            dual = mpq_get_d(ssx->cbar[j]);
         }
         if (k <= m)
         {  glp_set_row_stat(lp, k, stat);
            lp->row[k]->prim = prim;
            lp->row[k]->dual = dual;
         }
         else
         {  glp_set_col_stat(lp, k-m, stat);
            lp->col[k-m]->prim = prim;
            lp->col[k-m]->dual = dual;
            sum += lp->col[k-m]->coef * prim;
         }
      }
      lp->obj_val = sum;
done: /* delete the simplex solver workspace */
      ssx_delete(ssx);
#if 1 /* 23/XI-2015 */
      xassert(gmp_pool_count() == 0);
      gmp_free_mem();
#endif
      /* return to the application program */
      return ret;
}

/* eof */
