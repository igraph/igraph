/* glpapi07.c (exact simplex solver) */

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

#include "glpapi.h"
#include "glpssx.h"

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

static void load_data(SSX *ssx, LPX *lp)
{     /* load LP problem data into simplex solver workspace */
      int m = ssx->m;
      int n = ssx->n;
      int nnz = ssx->A_ptr[n+1]-1;
      int j, k, type, loc, len, *ind;
      double lb, ub, coef, *val;
      xassert(lpx_get_num_rows(lp) == m);
      xassert(lpx_get_num_cols(lp) == n);
      xassert(lpx_get_num_nz(lp) == nnz);
      /* types and bounds of rows and columns */
      for (k = 1; k <= m+n; k++)
      {  if (k <= m)
         {  type = lpx_get_row_type(lp, k);
            lb = lpx_get_row_lb(lp, k);
            ub = lpx_get_row_ub(lp, k);
         }
         else
         {  type = lpx_get_col_type(lp, k-m);
            lb = lpx_get_col_lb(lp, k-m);
            ub = lpx_get_col_ub(lp, k-m);
         }
         switch (type)
         {  case LPX_FR: type = SSX_FR; break;
            case LPX_LO: type = SSX_LO; break;
            case LPX_UP: type = SSX_UP; break;
            case LPX_DB: type = SSX_DB; break;
            case LPX_FX: type = SSX_FX; break;
            default: xassert(type != type);
         }
         ssx->type[k] = type;
         set_d_eps(ssx->lb[k], lb);
         set_d_eps(ssx->ub[k], ub);
      }
      /* optimization direction */
      switch (lpx_get_obj_dir(lp))
      {  case LPX_MIN: ssx->dir = SSX_MIN; break;
         case LPX_MAX: ssx->dir = SSX_MAX; break;
         default: xassert(lp != lp);
      }
      /* objective coefficients */
      for (k = 0; k <= m+n; k++)
      {  if (k == 0)
            coef = lpx_get_obj_coef(lp, 0);
         else if (k <= m)
            coef = 0.0;
         else
            coef = lpx_get_obj_coef(lp, k-m);
         set_d_eps(ssx->coef[k], coef);
      }
      /* constraint coefficients */
      ind = xcalloc(1+m, sizeof(int));
      val = xcalloc(1+m, sizeof(double));
      loc = 0;
      for (j = 1; j <= n; j++)
      {  ssx->A_ptr[j] = loc+1;
         len = lpx_get_mat_col(lp, j, ind, val);
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

static int load_basis(SSX *ssx, LPX *lp)
{     /* load current LP basis into simplex solver workspace */
      int m = ssx->m;
      int n = ssx->n;
      int *type = ssx->type;
      int *stat = ssx->stat;
      int *Q_row = ssx->Q_row;
      int *Q_col = ssx->Q_col;
      int i, j, k;
      xassert(lpx_get_num_rows(lp) == m);
      xassert(lpx_get_num_cols(lp) == n);
      /* statuses of rows and columns */
      for (k = 1; k <= m+n; k++)
      {  if (k <= m)
            stat[k] = lpx_get_row_stat(lp, k);
         else
            stat[k] = lpx_get_col_stat(lp, k-m);
         switch (stat[k])
         {  case LPX_BS:
               stat[k] = SSX_BS;
               break;
            case LPX_NL:
               stat[k] = SSX_NL;
               xassert(type[k] == SSX_LO || type[k] == SSX_DB);
               break;
            case LPX_NU:
               stat[k] = SSX_NU;
               xassert(type[k] == SSX_UP || type[k] == SSX_DB);
               break;
            case LPX_NF:
               stat[k] = SSX_NF;
               xassert(type[k] == SSX_FR);
               break;
            case LPX_NS:
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
      int m = lpx_get_num_rows(lp);
      int n = lpx_get_num_cols(lp);
      int nnz = lpx_get_num_nz(lp);
      int i, j, k, type, pst, dst, ret, *stat;
      double lb, ub, *prim, *dual, sum;
      if (parm == NULL)
         parm = &_parm, glp_init_smcp((glp_smcp *)parm);
      /* check control parameters */
      if (parm->it_lim < 0)
         xerror("glp_exact: it_lim = %d; invalid parameter\n",
            parm->it_lim);
      if (parm->tm_lim < 0)
         xerror("glp_exact: tm_lim = %d; invalid parameter\n",
            parm->tm_lim);
      /* the problem must have at least one row and one column */
      if (!(m > 0 && n > 0))
      {  xprintf("glp_exact: problem has no rows/columns\n");
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
         {  type = lpx_get_row_type(lp, k);
            lb = lpx_get_row_lb(lp, k);
            ub = lpx_get_row_ub(lp, k);
         }
         else
         {  type = lpx_get_col_type(lp, k-m);
            lb = lpx_get_col_lb(lp, k-m);
            ub = lpx_get_col_ub(lp, k-m);
         }
         if (type == LPX_DB && lb >= ub)
         {  xprintf("glp_exact: %s %d has invalid bounds\n",
               k <= m ? "row" : "column", k <= m ? k : k-m);
            return GLP_EBOUND;
         }
      }
      /* create the simplex solver workspace */
      xprintf("glp_exact: %d rows, %d columns, %d non-zeros\n",
         m, n, nnz);
#ifdef HAVE_GMP
      xprintf("GNU MP bignum library is being used\n");
#else
      xprintf("GLPK bignum module is being used\n");
      xprintf("(Consider installing GNU MP to attain a much better perf"
         "ormance.)\n");
#endif
      ssx = ssx_create(m, n, nnz);
      /* load LP problem data into the workspace */
      load_data(ssx, lp);
      /* load current LP basis into the workspace */
      if (load_basis(ssx, lp))
      {  xprintf("glp_exact: initial LP basis is invalid\n");
         ret = GLP_EBADB;
         goto done;
      }
      /* inherit some control parameters from the LP object */
#if 0
      ssx->it_lim = lpx_get_int_parm(lp, LPX_K_ITLIM);
      ssx->it_cnt = lpx_get_int_parm(lp, LPX_K_ITCNT);
      ssx->tm_lim = lpx_get_real_parm(lp, LPX_K_TMLIM);
#else
      ssx->it_lim = parm->it_lim;
      ssx->it_cnt = lp->it_cnt;
      ssx->tm_lim = (double)parm->tm_lim / 1000.0;
#endif
      ssx->out_frq = 5.0;
      ssx->tm_beg = xtime();
      ssx->tm_lag = xlset(0);
      /* solve LP */
      ret = ssx_driver(ssx);
      /* copy back some statistics to the LP object */
#if 0
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
            pst = LPX_P_FEAS, dst = LPX_D_FEAS;
            break;
         case 1:
            /* problem has no feasible solution */
            ret = 0;
            pst = LPX_P_NOFEAS, dst = LPX_D_INFEAS;
            break;
         case 2:
            /* problem has unbounded solution */
            ret = 0;
            pst = LPX_P_FEAS, dst = LPX_D_NOFEAS;
#if 1
            xassert(1 <= ssx->q && ssx->q <= n);
            lp->some = ssx->Q_col[m + ssx->q];
            xassert(1 <= lp->some && lp->some <= m+n);
#endif
            break;
         case 3:
            /* iteration limit exceeded (phase I) */
            ret = GLP_EITLIM;
            pst = LPX_P_INFEAS, dst = LPX_D_INFEAS;
            break;
         case 4:
            /* iteration limit exceeded (phase II) */
            ret = GLP_EITLIM;
            pst = LPX_P_FEAS, dst = LPX_D_INFEAS;
            break;
         case 5:
            /* time limit exceeded (phase I) */
            ret = GLP_ETMLIM;
            pst = LPX_P_INFEAS, dst = LPX_D_INFEAS;
            break;
         case 6:
            /* time limit exceeded (phase II) */
            ret = GLP_ETMLIM;
            pst = LPX_P_FEAS, dst = LPX_D_INFEAS;
            break;
         case 7:
            /* initial basis matrix is singular */
            ret = GLP_ESING;
            goto done;
         default:
            xassert(ret != ret);
      }
      /* obtain final basic solution components */
      stat = xcalloc(1+m+n, sizeof(int));
      prim = xcalloc(1+m+n, sizeof(double));
      dual = xcalloc(1+m+n, sizeof(double));
      for (k = 1; k <= m+n; k++)
      {  if (ssx->stat[k] == SSX_BS)
         {  i = ssx->Q_row[k]; /* x[k] = xB[i] */
            xassert(1 <= i && i <= m);
            stat[k] = LPX_BS;
            prim[k] = mpq_get_d(ssx->bbar[i]);
            dual[k] = 0.0;
         }
         else
         {  j = ssx->Q_row[k] - m; /* x[k] = xN[j] */
            xassert(1 <= j && j <= n);
            switch (ssx->stat[k])
            {  case SSX_NF:
                  stat[k] = LPX_NF;
                  prim[k] = 0.0;
                  break;
               case SSX_NL:
                  stat[k] = LPX_NL;
                  prim[k] = mpq_get_d(ssx->lb[k]);
                  break;
               case SSX_NU:
                  stat[k] = LPX_NU;
                  prim[k] = mpq_get_d(ssx->ub[k]);
                  break;
               case SSX_NS:
                  stat[k] = LPX_NS;
                  prim[k] = mpq_get_d(ssx->lb[k]);
                  break;
               default:
                  xassert(ssx != ssx);
            }
            dual[k] = mpq_get_d(ssx->cbar[j]);
         }
      }
      /* and store them into the LP object */
      pst = pst - LPX_P_UNDEF + GLP_UNDEF;
      dst = dst - LPX_D_UNDEF + GLP_UNDEF;
      for (k = 1; k <= m+n; k++)
         stat[k] = stat[k] - LPX_BS + GLP_BS;
      sum = lpx_get_obj_coef(lp, 0);
      for (j = 1; j <= n; j++)
         sum += lpx_get_obj_coef(lp, j) * prim[m+j];
      lpx_put_solution(lp, 1, &pst, &dst, &sum,
         &stat[0], &prim[0], &dual[0], &stat[m], &prim[m], &dual[m]);
      xfree(stat);
      xfree(prim);
      xfree(dual);
done: /* delete the simplex solver workspace */
      ssx_delete(ssx);
      /* return to the application program */
      return ret;
}

/* eof */
