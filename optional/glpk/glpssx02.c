/* glpssx02.c */

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

#include "glpenv.h"
#include "glpssx.h"

static void show_progress(SSX *ssx, int phase)
{     /* this auxiliary routine displays information about progress of
         the search */
      int i, def = 0;
      for (i = 1; i <= ssx->m; i++)
         if (ssx->type[ssx->Q_col[i]] == SSX_FX) def++;
      xprintf("%s%6d:   %s = %22.15g   (%d)\n", phase == 1 ? " " : "*",
         ssx->it_cnt, phase == 1 ? "infsum" : "objval",
         mpq_get_d(ssx->bbar[0]), def);
#if 0
      ssx->tm_lag = utime();
#else
      ssx->tm_lag = xtime();
#endif
      return;
}

/*----------------------------------------------------------------------
// ssx_phase_I - find primal feasible solution.
//
// This routine implements phase I of the primal simplex method.
//
// On exit the routine returns one of the following codes:
//
// 0 - feasible solution found;
// 1 - problem has no feasible solution;
// 2 - iterations limit exceeded;
// 3 - time limit exceeded.
----------------------------------------------------------------------*/

int ssx_phase_I(SSX *ssx)
{     int m = ssx->m;
      int n = ssx->n;
      int *type = ssx->type;
      mpq_t *lb = ssx->lb;
      mpq_t *ub = ssx->ub;
      mpq_t *coef = ssx->coef;
      int *A_ptr = ssx->A_ptr;
      int *A_ind = ssx->A_ind;
      mpq_t *A_val = ssx->A_val;
      int *Q_col = ssx->Q_col;
      mpq_t *bbar = ssx->bbar;
      mpq_t *pi = ssx->pi;
      mpq_t *cbar = ssx->cbar;
      int *orig_type, orig_dir;
      mpq_t *orig_lb, *orig_ub, *orig_coef;
      int i, k, ret;
      /* save components of the original LP problem, which are changed
         by the routine */
      orig_type = xcalloc(1+m+n, sizeof(int));
      orig_lb = xcalloc(1+m+n, sizeof(mpq_t));
      orig_ub = xcalloc(1+m+n, sizeof(mpq_t));
      orig_coef = xcalloc(1+m+n, sizeof(mpq_t));
      for (k = 1; k <= m+n; k++)
      {  orig_type[k] = type[k];
         mpq_init(orig_lb[k]);
         mpq_set(orig_lb[k], lb[k]);
         mpq_init(orig_ub[k]);
         mpq_set(orig_ub[k], ub[k]);
      }
      orig_dir = ssx->dir;
      for (k = 0; k <= m+n; k++)
      {  mpq_init(orig_coef[k]);
         mpq_set(orig_coef[k], coef[k]);
      }
      /* build an artificial basic solution, which is primal feasible,
         and also build an auxiliary objective function to minimize the
         sum of infeasibilities for the original problem */
      ssx->dir = SSX_MIN;
      for (k = 0; k <= m+n; k++) mpq_set_si(coef[k], 0, 1);
      mpq_set_si(bbar[0], 0, 1);
      for (i = 1; i <= m; i++)
      {  int t;
         k = Q_col[i]; /* x[k] = xB[i] */
         t = type[k];
         if (t == SSX_LO || t == SSX_DB || t == SSX_FX)
         {  /* in the original problem x[k] has lower bound */
            if (mpq_cmp(bbar[i], lb[k]) < 0)
            {  /* which is violated */
               type[k] = SSX_UP;
               mpq_set(ub[k], lb[k]);
               mpq_set_si(lb[k], 0, 1);
               mpq_set_si(coef[k], -1, 1);
               mpq_add(bbar[0], bbar[0], ub[k]);
               mpq_sub(bbar[0], bbar[0], bbar[i]);
            }
         }
         if (t == SSX_UP || t == SSX_DB || t == SSX_FX)
         {  /* in the original problem x[k] has upper bound */
            if (mpq_cmp(bbar[i], ub[k]) > 0)
            {  /* which is violated */
               type[k] = SSX_LO;
               mpq_set(lb[k], ub[k]);
               mpq_set_si(ub[k], 0, 1);
               mpq_set_si(coef[k], +1, 1);
               mpq_add(bbar[0], bbar[0], bbar[i]);
               mpq_sub(bbar[0], bbar[0], lb[k]);
            }
         }
      }
      /* now the initial basic solution should be primal feasible due
         to changes of bounds of some basic variables, which turned to
         implicit artifical variables */
      /* compute simplex multipliers and reduced costs */
      ssx_eval_pi(ssx);
      ssx_eval_cbar(ssx);
      /* display initial progress of the search */
      show_progress(ssx, 1);
      /* main loop starts here */
      for (;;)
      {  /* display current progress of the search */
#if 0
         if (utime() - ssx->tm_lag >= ssx->out_frq - 0.001)
#else
         if (xdifftime(xtime(), ssx->tm_lag) >= ssx->out_frq - 0.001)
#endif
            show_progress(ssx, 1);
         /* we do not need to wait until all artificial variables have
            left the basis */
         if (mpq_sgn(bbar[0]) == 0)
         {  /* the sum of infeasibilities is zero, therefore the current
               solution is primal feasible for the original problem */
            ret = 0;
            break;
         }
         /* check if the iterations limit has been exhausted */
         if (ssx->it_lim == 0)
         {  ret = 2;
            break;
         }
         /* check if the time limit has been exhausted */
#if 0
         if (ssx->tm_lim >= 0.0 && ssx->tm_lim <= utime() - ssx->tm_beg)
#else
         if (ssx->tm_lim >= 0.0 &&
             ssx->tm_lim <= xdifftime(xtime(), ssx->tm_beg))
#endif
         {  ret = 3;
            break;
         }
         /* choose non-basic variable xN[q] */
         ssx_chuzc(ssx);
         /* if xN[q] cannot be chosen, the sum of infeasibilities is
            minimal but non-zero; therefore the original problem has no
            primal feasible solution */
         if (ssx->q == 0)
         {  ret = 1;
            break;
         }
         /* compute q-th column of the simplex table */
         ssx_eval_col(ssx);
         /* choose basic variable xB[p] */
         ssx_chuzr(ssx);
         /* the sum of infeasibilities cannot be negative, therefore
            the auxiliary lp problem cannot have unbounded solution */
         xassert(ssx->p != 0);
         /* update values of basic variables */
         ssx_update_bbar(ssx);
         if (ssx->p > 0)
         {  /* compute p-th row of the inverse inv(B) */
            ssx_eval_rho(ssx);
            /* compute p-th row of the simplex table */
            ssx_eval_row(ssx);
            xassert(mpq_cmp(ssx->aq[ssx->p], ssx->ap[ssx->q]) == 0);
            /* update simplex multipliers */
            ssx_update_pi(ssx);
            /* update reduced costs of non-basic variables */
            ssx_update_cbar(ssx);
         }
         /* xB[p] is leaving the basis; if it is implicit artificial
            variable, the corresponding residual vanishes; therefore
            bounds of this variable should be restored to the original
            values */
         if (ssx->p > 0)
         {  k = Q_col[ssx->p]; /* x[k] = xB[p] */
            if (type[k] != orig_type[k])
            {  /* x[k] is implicit artificial variable */
               type[k] = orig_type[k];
               mpq_set(lb[k], orig_lb[k]);
               mpq_set(ub[k], orig_ub[k]);
               xassert(ssx->p_stat == SSX_NL || ssx->p_stat == SSX_NU);
               ssx->p_stat = (ssx->p_stat == SSX_NL ? SSX_NU : SSX_NL);
               if (type[k] == SSX_FX) ssx->p_stat = SSX_NS;
               /* nullify the objective coefficient at x[k] */
               mpq_set_si(coef[k], 0, 1);
               /* since coef[k] has been changed, we need to compute
                  new reduced cost of x[k], which it will have in the
                  adjacent basis */
               /* the formula d[j] = cN[j] - pi' * N[j] is used (note
                  that the vector pi is not changed, because it depends
                  on objective coefficients at basic variables, but in
                  the adjacent basis, for which the vector pi has been
                  just recomputed, x[k] is non-basic) */
               if (k <= m)
               {  /* x[k] is auxiliary variable */
                  mpq_neg(cbar[ssx->q], pi[k]);
               }
               else
               {  /* x[k] is structural variable */
                  int ptr;
                  mpq_t temp;
                  mpq_init(temp);
                  mpq_set_si(cbar[ssx->q], 0, 1);
                  for (ptr = A_ptr[k-m]; ptr < A_ptr[k-m+1]; ptr++)
                  {  mpq_mul(temp, pi[A_ind[ptr]], A_val[ptr]);
                     mpq_add(cbar[ssx->q], cbar[ssx->q], temp);
                  }
                  mpq_clear(temp);
               }
            }
         }
         /* jump to the adjacent vertex of the polyhedron */
         ssx_change_basis(ssx);
         /* one simplex iteration has been performed */
         if (ssx->it_lim > 0) ssx->it_lim--;
         ssx->it_cnt++;
      }
      /* display final progress of the search */
      show_progress(ssx, 1);
      /* restore components of the original problem, which were changed
         by the routine */
      for (k = 1; k <= m+n; k++)
      {  type[k] = orig_type[k];
         mpq_set(lb[k], orig_lb[k]);
         mpq_clear(orig_lb[k]);
         mpq_set(ub[k], orig_ub[k]);
         mpq_clear(orig_ub[k]);
      }
      ssx->dir = orig_dir;
      for (k = 0; k <= m+n; k++)
      {  mpq_set(coef[k], orig_coef[k]);
         mpq_clear(orig_coef[k]);
      }
      xfree(orig_type);
      xfree(orig_lb);
      xfree(orig_ub);
      xfree(orig_coef);
      /* return to the calling program */
      return ret;
}

/*----------------------------------------------------------------------
// ssx_phase_II - find optimal solution.
//
// This routine implements phase II of the primal simplex method.
//
// On exit the routine returns one of the following codes:
//
// 0 - optimal solution found;
// 1 - problem has unbounded solution;
// 2 - iterations limit exceeded;
// 3 - time limit exceeded.
----------------------------------------------------------------------*/

int ssx_phase_II(SSX *ssx)
{     int ret;
      /* display initial progress of the search */
      show_progress(ssx, 2);
      /* main loop starts here */
      for (;;)
      {  /* display current progress of the search */
#if 0
         if (utime() - ssx->tm_lag >= ssx->out_frq - 0.001)
#else
         if (xdifftime(xtime(), ssx->tm_lag) >= ssx->out_frq - 0.001)
#endif
            show_progress(ssx, 2);
         /* check if the iterations limit has been exhausted */
         if (ssx->it_lim == 0)
         {  ret = 2;
            break;
         }
         /* check if the time limit has been exhausted */
#if 0
         if (ssx->tm_lim >= 0.0 && ssx->tm_lim <= utime() - ssx->tm_beg)
#else
         if (ssx->tm_lim >= 0.0 &&
             ssx->tm_lim <= xdifftime(xtime(), ssx->tm_beg))
#endif
         {  ret = 3;
            break;
         }
         /* choose non-basic variable xN[q] */
         ssx_chuzc(ssx);
         /* if xN[q] cannot be chosen, the current basic solution is
            dual feasible and therefore optimal */
         if (ssx->q == 0)
         {  ret = 0;
            break;
         }
         /* compute q-th column of the simplex table */
         ssx_eval_col(ssx);
         /* choose basic variable xB[p] */
         ssx_chuzr(ssx);
         /* if xB[p] cannot be chosen, the problem has no dual feasible
            solution (i.e. unbounded) */
         if (ssx->p == 0)
         {  ret = 1;
            break;
         }
         /* update values of basic variables */
         ssx_update_bbar(ssx);
         if (ssx->p > 0)
         {  /* compute p-th row of the inverse inv(B) */
            ssx_eval_rho(ssx);
            /* compute p-th row of the simplex table */
            ssx_eval_row(ssx);
            xassert(mpq_cmp(ssx->aq[ssx->p], ssx->ap[ssx->q]) == 0);
#if 0
            /* update simplex multipliers */
            ssx_update_pi(ssx);
#endif
            /* update reduced costs of non-basic variables */
            ssx_update_cbar(ssx);
         }
         /* jump to the adjacent vertex of the polyhedron */
         ssx_change_basis(ssx);
         /* one simplex iteration has been performed */
         if (ssx->it_lim > 0) ssx->it_lim--;
         ssx->it_cnt++;
      }
      /* display final progress of the search */
      show_progress(ssx, 2);
      /* return to the calling program */
      return ret;
}

/*----------------------------------------------------------------------
// ssx_driver - base driver to exact simplex method.
//
// This routine is a base driver to a version of the primal simplex
// method using exact (bignum) arithmetic.
//
// On exit the routine returns one of the following codes:
//
// 0 - optimal solution found;
// 1 - problem has no feasible solution;
// 2 - problem has unbounded solution;
// 3 - iterations limit exceeded (phase I);
// 4 - iterations limit exceeded (phase II);
// 5 - time limit exceeded (phase I);
// 6 - time limit exceeded (phase II);
// 7 - initial basis matrix is exactly singular.
----------------------------------------------------------------------*/

int ssx_driver(SSX *ssx)
{     int m = ssx->m;
      int *type = ssx->type;
      mpq_t *lb = ssx->lb;
      mpq_t *ub = ssx->ub;
      int *Q_col = ssx->Q_col;
      mpq_t *bbar = ssx->bbar;
      int i, k, ret;
      ssx->tm_beg = xtime();
      /* factorize the initial basis matrix */
      if (ssx_factorize(ssx))
      {  xprintf("Initial basis matrix is singular\n");
         ret = 7;
         goto done;
      }
      /* compute values of basic variables */
      ssx_eval_bbar(ssx);
      /* check if the initial basic solution is primal feasible */
      for (i = 1; i <= m; i++)
      {  int t;
         k = Q_col[i]; /* x[k] = xB[i] */
         t = type[k];
         if (t == SSX_LO || t == SSX_DB || t == SSX_FX)
         {  /* x[k] has lower bound */
            if (mpq_cmp(bbar[i], lb[k]) < 0)
            {  /* which is violated */
               break;
            }
         }
         if (t == SSX_UP || t == SSX_DB || t == SSX_FX)
         {  /* x[k] has upper bound */
            if (mpq_cmp(bbar[i], ub[k]) > 0)
            {  /* which is violated */
               break;
            }
         }
      }
      if (i > m)
      {  /* no basic variable violates its bounds */
         ret = 0;
         goto skip;
      }
      /* phase I: find primal feasible solution */
      ret = ssx_phase_I(ssx);
      switch (ret)
      {  case 0:
            ret = 0;
            break;
         case 1:
            xprintf("PROBLEM HAS NO FEASIBLE SOLUTION\n");
            ret = 1;
            break;
         case 2:
            xprintf("ITERATIONS LIMIT EXCEEDED; SEARCH TERMINATED\n");
            ret = 3;
            break;
         case 3:
            xprintf("TIME LIMIT EXCEEDED; SEARCH TERMINATED\n");
            ret = 5;
            break;
         default:
            xassert(ret != ret);
      }
      /* compute values of basic variables (actually only the objective
         value needs to be computed) */
      ssx_eval_bbar(ssx);
skip: /* compute simplex multipliers */
      ssx_eval_pi(ssx);
      /* compute reduced costs of non-basic variables */
      ssx_eval_cbar(ssx);
      /* if phase I failed, do not start phase II */
      if (ret != 0) goto done;
      /* phase II: find optimal solution */
      ret = ssx_phase_II(ssx);
      switch (ret)
      {  case 0:
            xprintf("OPTIMAL SOLUTION FOUND\n");
            ret = 0;
            break;
         case 1:
            xprintf("PROBLEM HAS UNBOUNDED SOLUTION\n");
            ret = 2;
            break;
         case 2:
            xprintf("ITERATIONS LIMIT EXCEEDED; SEARCH TERMINATED\n");
            ret = 4;
            break;
         case 3:
            xprintf("TIME LIMIT EXCEEDED; SEARCH TERMINATED\n");
            ret = 6;
            break;
         default:
            xassert(ret != ret);
      }
done: /* decrease the time limit by the spent amount of time */
      if (ssx->tm_lim >= 0.0)
#if 0
      {  ssx->tm_lim -= utime() - ssx->tm_beg;
#else
      {  ssx->tm_lim -= xdifftime(xtime(), ssx->tm_beg);
#endif
         if (ssx->tm_lim < 0.0) ssx->tm_lim = 0.0;
      }
      return ret;
}

/* eof */
