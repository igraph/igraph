/* glpapi08.c (interior-point method routines) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2000-2013 Free Software Foundation, Inc.
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
#include "glpipm.h"
#include "npp.h"

/***********************************************************************
*  NAME
*
*  glp_interior - solve LP problem with the interior-point method
*
*  SYNOPSIS
*
*  int glp_interior(glp_prob *P, const glp_iptcp *parm);
*
*  The routine glp_interior is a driver to the LP solver based on the
*  interior-point method.
*
*  The interior-point solver has a set of control parameters. Values of
*  the control parameters can be passed in a structure glp_iptcp, which
*  the parameter parm points to.
*
*  Currently this routine implements an easy variant of the primal-dual
*  interior-point method based on Mehrotra's technique.
*
*  This routine transforms the original LP problem to an equivalent LP
*  problem in the standard formulation (all constraints are equalities,
*  all variables are non-negative), calls the routine ipm_main to solve
*  the transformed problem, and then transforms an obtained solution to
*  the solution of the original problem.
*
*  RETURNS
*
*  0  The LP problem instance has been successfully solved. This code
*     does not necessarily mean that the solver has found optimal
*     solution. It only means that the solution process was successful.
*
*  GLP_EFAIL
*     The problem has no rows/columns.
*
*  GLP_ENOCVG
*     Very slow convergence or divergence.
*
*  GLP_EITLIM
*     Iteration limit exceeded.
*
*  GLP_EINSTAB
*     Numerical instability on solving Newtonian system. */

static void transform(NPP *npp)
{     /* transform LP to the standard formulation */
      NPPROW *row, *prev_row;
      NPPCOL *col, *prev_col;
      for (row = npp->r_tail; row != NULL; row = prev_row)
      {  prev_row = row->prev;
         if (row->lb == -DBL_MAX && row->ub == +DBL_MAX)
            npp_free_row(npp, row);
         else if (row->lb == -DBL_MAX)
            npp_leq_row(npp, row);
         else if (row->ub == +DBL_MAX)
            npp_geq_row(npp, row);
         else if (row->lb != row->ub)
         {  if (fabs(row->lb) < fabs(row->ub))
               npp_geq_row(npp, row);
            else
               npp_leq_row(npp, row);
         }
      }
      for (col = npp->c_tail; col != NULL; col = prev_col)
      {  prev_col = col->prev;
         if (col->lb == -DBL_MAX && col->ub == +DBL_MAX)
            npp_free_col(npp, col);
         else if (col->lb == -DBL_MAX)
            npp_ubnd_col(npp, col);
         else if (col->ub == +DBL_MAX)
         {  if (col->lb != 0.0)
               npp_lbnd_col(npp, col);
         }
         else if (col->lb != col->ub)
         {  if (fabs(col->lb) < fabs(col->ub))
            {  if (col->lb != 0.0)
                  npp_lbnd_col(npp, col);
            }
            else
               npp_ubnd_col(npp, col);
            npp_dbnd_col(npp, col);
         }
         else
            npp_fixed_col(npp, col);
      }
      for (row = npp->r_head; row != NULL; row = row->next)
         xassert(row->lb == row->ub);
      for (col = npp->c_head; col != NULL; col = col->next)
         xassert(col->lb == 0.0 && col->ub == +DBL_MAX);
      return;
}

int glp_interior(glp_prob *P, const glp_iptcp *parm)
{     glp_iptcp _parm;
      GLPROW *row;
      GLPCOL *col;
      NPP *npp = NULL;
      glp_prob *prob = NULL;
      int i, j, ret;
      /* check control parameters */
      if (parm == NULL)
         glp_init_iptcp(&_parm), parm = &_parm;
      if (!(parm->msg_lev == GLP_MSG_OFF ||
            parm->msg_lev == GLP_MSG_ERR ||
            parm->msg_lev == GLP_MSG_ON  ||
            parm->msg_lev == GLP_MSG_ALL))
         xerror("glp_interior: msg_lev = %d; invalid parameter\n",
            parm->msg_lev);
      if (!(parm->ord_alg == GLP_ORD_NONE ||
            parm->ord_alg == GLP_ORD_QMD ||
            parm->ord_alg == GLP_ORD_AMD ||
            parm->ord_alg == GLP_ORD_SYMAMD))
         xerror("glp_interior: ord_alg = %d; invalid parameter\n",
            parm->ord_alg);
      /* interior-point solution is currently undefined */
      P->ipt_stat = GLP_UNDEF;
      P->ipt_obj = 0.0;
      /* check bounds of double-bounded variables */
      for (i = 1; i <= P->m; i++)
      {  row = P->row[i];
         if (row->type == GLP_DB && row->lb >= row->ub)
         {  if (parm->msg_lev >= GLP_MSG_ERR)
               xprintf("glp_interior: row %d: lb = %g, ub = %g; incorre"
                  "ct bounds\n", i, row->lb, row->ub);
            ret = GLP_EBOUND;
            goto done;
         }
      }
      for (j = 1; j <= P->n; j++)
      {  col = P->col[j];
         if (col->type == GLP_DB && col->lb >= col->ub)
         {  if (parm->msg_lev >= GLP_MSG_ERR)
               xprintf("glp_interior: column %d: lb = %g, ub = %g; inco"
                  "rrect bounds\n", j, col->lb, col->ub);
            ret = GLP_EBOUND;
            goto done;
         }
      }
      /* transform LP to the standard formulation */
      if (parm->msg_lev >= GLP_MSG_ALL)
         xprintf("Original LP has %d row(s), %d column(s), and %d non-z"
            "ero(s)\n", P->m, P->n, P->nnz);
      npp = npp_create_wksp();
      npp_load_prob(npp, P, GLP_OFF, GLP_IPT, GLP_ON);
      transform(npp);
      prob = glp_create_prob();
      npp_build_prob(npp, prob);
      if (parm->msg_lev >= GLP_MSG_ALL)
         xprintf("Working LP has %d row(s), %d column(s), and %d non-ze"
            "ro(s)\n", prob->m, prob->n, prob->nnz);
#if 1
      /* currently empty problem cannot be solved */
      if (!(prob->m > 0 && prob->n > 0))
      {  if (parm->msg_lev >= GLP_MSG_ERR)
            xprintf("glp_interior: unable to solve empty problem\n");
         ret = GLP_EFAIL;
         goto done;
      }
#endif
      /* scale the resultant LP */
      {  ENV *env = get_env_ptr();
         int term_out = env->term_out;
         env->term_out = GLP_OFF;
         glp_scale_prob(prob, GLP_SF_EQ);
         env->term_out = term_out;
      }
      /* warn about dense columns */
      if (parm->msg_lev >= GLP_MSG_ON && prob->m >= 200)
      {  int len, cnt = 0;
         for (j = 1; j <= prob->n; j++)
         {  len = glp_get_mat_col(prob, j, NULL, NULL);
            if ((double)len >= 0.20 * (double)prob->m) cnt++;
         }
         if (cnt == 1)
            xprintf("WARNING: PROBLEM HAS ONE DENSE COLUMN\n");
         else if (cnt > 0)
            xprintf("WARNING: PROBLEM HAS %d DENSE COLUMNS\n", cnt);
      }
      /* solve the transformed LP */
      ret = ipm_solve(prob, parm);
      /* postprocess solution from the transformed LP */
      npp_postprocess(npp, prob);
      /* and store solution to the original LP */
      npp_unload_sol(npp, P);
done: /* free working program objects */
      if (npp != NULL) npp_delete_wksp(npp);
      if (prob != NULL) glp_delete_prob(prob);
      /* return to the application program */
      return ret;
}

/***********************************************************************
*  NAME
*
*  glp_init_iptcp - initialize interior-point solver control parameters
*
*  SYNOPSIS
*
*  void glp_init_iptcp(glp_iptcp *parm);
*
*  DESCRIPTION
*
*  The routine glp_init_iptcp initializes control parameters, which are
*  used by the interior-point solver, with default values.
*
*  Default values of the control parameters are stored in the glp_iptcp
*  structure, which the parameter parm points to. */

void glp_init_iptcp(glp_iptcp *parm)
{     parm->msg_lev = GLP_MSG_ALL;
      parm->ord_alg = GLP_ORD_AMD;
      return;
}

/***********************************************************************
*  NAME
*
*  glp_ipt_status - retrieve status of interior-point solution
*
*  SYNOPSIS
*
*  int glp_ipt_status(glp_prob *lp);
*
*  RETURNS
*
*  The routine glp_ipt_status reports the status of solution found by
*  the interior-point solver as follows:
*
*  GLP_UNDEF  - interior-point solution is undefined;
*  GLP_OPT    - interior-point solution is optimal;
*  GLP_INFEAS - interior-point solution is infeasible;
*  GLP_NOFEAS - no feasible solution exists. */

int glp_ipt_status(glp_prob *lp)
{     int ipt_stat = lp->ipt_stat;
      return ipt_stat;
}

/***********************************************************************
*  NAME
*
*  glp_ipt_obj_val - retrieve objective value (interior point)
*
*  SYNOPSIS
*
*  double glp_ipt_obj_val(glp_prob *lp);
*
*  RETURNS
*
*  The routine glp_ipt_obj_val returns value of the objective function
*  for interior-point solution. */

double glp_ipt_obj_val(glp_prob *lp)
{     /*struct LPXCPS *cps = lp->cps;*/
      double z;
      z = lp->ipt_obj;
      /*if (cps->round && fabs(z) < 1e-9) z = 0.0;*/
      return z;
}

/***********************************************************************
*  NAME
*
*  glp_ipt_row_prim - retrieve row primal value (interior point)
*
*  SYNOPSIS
*
*  double glp_ipt_row_prim(glp_prob *lp, int i);
*
*  RETURNS
*
*  The routine glp_ipt_row_prim returns primal value of the auxiliary
*  variable associated with i-th row. */

double glp_ipt_row_prim(glp_prob *lp, int i)
{     /*struct LPXCPS *cps = lp->cps;*/
      double pval;
      if (!(1 <= i && i <= lp->m))
         xerror("glp_ipt_row_prim: i = %d; row number out of range\n",
            i);
      pval = lp->row[i]->pval;
      /*if (cps->round && fabs(pval) < 1e-9) pval = 0.0;*/
      return pval;
}

/***********************************************************************
*  NAME
*
*  glp_ipt_row_dual - retrieve row dual value (interior point)
*
*  SYNOPSIS
*
*  double glp_ipt_row_dual(glp_prob *lp, int i);
*
*  RETURNS
*
*  The routine glp_ipt_row_dual returns dual value (i.e. reduced cost)
*  of the auxiliary variable associated with i-th row. */

double glp_ipt_row_dual(glp_prob *lp, int i)
{     /*struct LPXCPS *cps = lp->cps;*/
      double dval;
      if (!(1 <= i && i <= lp->m))
         xerror("glp_ipt_row_dual: i = %d; row number out of range\n",
            i);
      dval = lp->row[i]->dval;
      /*if (cps->round && fabs(dval) < 1e-9) dval = 0.0;*/
      return dval;
}

/***********************************************************************
*  NAME
*
*  glp_ipt_col_prim - retrieve column primal value (interior point)
*
*  SYNOPSIS
*
*  double glp_ipt_col_prim(glp_prob *lp, int j);
*
*  RETURNS
*
*  The routine glp_ipt_col_prim returns primal value of the structural
*  variable associated with j-th column. */

double glp_ipt_col_prim(glp_prob *lp, int j)
{     /*struct LPXCPS *cps = lp->cps;*/
      double pval;
      if (!(1 <= j && j <= lp->n))
         xerror("glp_ipt_col_prim: j = %d; column number out of range\n"
            , j);
      pval = lp->col[j]->pval;
      /*if (cps->round && fabs(pval) < 1e-9) pval = 0.0;*/
      return pval;
}

/***********************************************************************
*  NAME
*
*  glp_ipt_col_dual - retrieve column dual value (interior point)
*
*  SYNOPSIS
*
*  double glp_ipt_col_dual(glp_prob *lp, int j);
*
*  RETURNS
*
*  The routine glp_ipt_col_dual returns dual value (i.e. reduced cost)
*  of the structural variable associated with j-th column. */

double glp_ipt_col_dual(glp_prob *lp, int j)
{     /*struct LPXCPS *cps = lp->cps;*/
      double dval;
      if (!(1 <= j && j <= lp->n))
         xerror("glp_ipt_col_dual: j = %d; column number out of range\n"
            , j);
      dval = lp->col[j]->dval;
      /*if (cps->round && fabs(dval) < 1e-9) dval = 0.0;*/
      return dval;
}

/* eof */
