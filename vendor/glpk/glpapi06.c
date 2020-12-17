/* glpapi06.c (simplex method routines) */

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
#pragma clang diagnostic ignored "-Wself-assign"
#endif

#include "glpios.h"
#include "glpnpp.h"
#include "glpspx.h"

/***********************************************************************
*  NAME
*
*  glp_simplex - solve LP problem with the simplex method
*
*  SYNOPSIS
*
*  int glp_simplex(glp_prob *P, const glp_smcp *parm);
*
*  DESCRIPTION
*
*  The routine glp_simplex is a driver to the LP solver based on the
*  simplex method. This routine retrieves problem data from the
*  specified problem object, calls the solver to solve the problem
*  instance, and stores results of computations back into the problem
*  object.
*
*  The simplex solver has a set of control parameters. Values of the
*  control parameters can be passed in a structure glp_smcp, which the
*  parameter parm points to.
*
*  The parameter parm can be specified as NULL, in which case the LP
*  solver uses default settings.
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
*     to the initial basis is singular within the working precision.
*
*  GLP_ECOND
*     Unable to start the search, because the basis matrix correspodning
*     to the initial basis is ill-conditioned, i.e. its condition number
*     is too large.
*
*  GLP_EBOUND
*     Unable to start the search, because some double-bounded variables
*     have incorrect bounds.
*
*  GLP_EFAIL
*     The search was prematurely terminated due to the solver failure.
*
*  GLP_EOBJLL
*     The search was prematurely terminated, because the objective
*     function being maximized has reached its lower limit and continues
*     decreasing (dual simplex only).
*
*  GLP_EOBJUL
*     The search was prematurely terminated, because the objective
*     function being minimized has reached its upper limit and continues
*     increasing (dual simplex only).
*
*  GLP_EITLIM
*     The search was prematurely terminated, because the simplex
*     iteration limit has been exceeded.
*
*  GLP_ETMLIM
*     The search was prematurely terminated, because the time limit has
*     been exceeded.
*
*  GLP_ENOPFS
*     The LP problem instance has no primal feasible solution (only if
*     the LP presolver is used).
*
*  GLP_ENODFS
*     The LP problem instance has no dual feasible solution (only if the
*     LP presolver is used). */

static void trivial_lp(glp_prob *P, const glp_smcp *parm)
{     /* solve trivial LP which has empty constraint matrix */
      GLPROW *row;
      GLPCOL *col;
      int i, j;
      double p_infeas, d_infeas, zeta;
      P->valid = 0;
      P->pbs_stat = P->dbs_stat = GLP_FEAS;
      P->obj_val = P->c0;
      P->some = 0;
      p_infeas = d_infeas = 0.0;
      /* make all auxiliary variables basic */
      for (i = 1; i <= P->m; i++)
      {  row = P->row[i];
         row->stat = GLP_BS;
         row->prim = row->dual = 0.0;
         /* check primal feasibility */
         if (row->type == GLP_LO || row->type == GLP_DB ||
             row->type == GLP_FX)
         {  /* row has lower bound */
            if (row->lb > + parm->tol_bnd)
            {  P->pbs_stat = GLP_NOFEAS;
               if (P->some == 0 && parm->meth != GLP_PRIMAL)
                  P->some = i;
            }
            if (p_infeas < + row->lb)
               p_infeas = + row->lb;
         }
         if (row->type == GLP_UP || row->type == GLP_DB ||
             row->type == GLP_FX)
         {  /* row has upper bound */
            if (row->ub < - parm->tol_bnd)
            {  P->pbs_stat = GLP_NOFEAS;
               if (P->some == 0 && parm->meth != GLP_PRIMAL)
                  P->some = i;
            }
            if (p_infeas < - row->ub)
               p_infeas = - row->ub;
         }
      }
      /* determine scale factor for the objective row */
      zeta = 1.0;
      for (j = 1; j <= P->n; j++)
      {  col = P->col[j];
         if (zeta < fabs(col->coef)) zeta = fabs(col->coef);
      }
      zeta = (P->dir == GLP_MIN ? +1.0 : -1.0) / zeta;
      /* make all structural variables non-basic */
      for (j = 1; j <= P->n; j++)
      {  col = P->col[j];
         if (col->type == GLP_FR)
            col->stat = GLP_NF, col->prim = 0.0;
         else if (col->type == GLP_LO)
lo:         col->stat = GLP_NL, col->prim = col->lb;
         else if (col->type == GLP_UP)
up:         col->stat = GLP_NU, col->prim = col->ub;
         else if (col->type == GLP_DB)
         {  if (zeta * col->coef > 0.0)
               goto lo;
            else if (zeta * col->coef < 0.0)
               goto up;
            else if (fabs(col->lb) <= fabs(col->ub))
               goto lo;
            else
               goto up;
         }
         else if (col->type == GLP_FX)
            col->stat = GLP_NS, col->prim = col->lb;
         col->dual = col->coef;
         P->obj_val += col->coef * col->prim;
         /* check dual feasibility */
         if (col->type == GLP_FR || col->type == GLP_LO)
         {  /* column has no upper bound */
            if (zeta * col->dual < - parm->tol_dj)
            {  P->dbs_stat = GLP_NOFEAS;
               if (P->some == 0 && parm->meth == GLP_PRIMAL)
                  P->some = P->m + j;
            }
            if (d_infeas < - zeta * col->dual)
               d_infeas = - zeta * col->dual;
         }
         if (col->type == GLP_FR || col->type == GLP_UP)
         {  /* column has no lower bound */
            if (zeta * col->dual > + parm->tol_dj)
            {  P->dbs_stat = GLP_NOFEAS;
               if (P->some == 0 && parm->meth == GLP_PRIMAL)
                  P->some = P->m + j;
            }
            if (d_infeas < + zeta * col->dual)
               d_infeas = + zeta * col->dual;
         }
      }
      /* simulate the simplex solver output */
      if (parm->msg_lev >= GLP_MSG_ON && parm->out_dly == 0)
      {  xprintf("~%6d: obj = %17.9e  infeas = %10.3e\n", P->it_cnt,
            P->obj_val, parm->meth == GLP_PRIMAL ? p_infeas : d_infeas);
      }
      if (parm->msg_lev >= GLP_MSG_ALL && parm->out_dly == 0)
      {  if (P->pbs_stat == GLP_FEAS && P->dbs_stat == GLP_FEAS)
            xprintf("OPTIMAL SOLUTION FOUND\n");
         else if (P->pbs_stat == GLP_NOFEAS)
            xprintf("PROBLEM HAS NO FEASIBLE SOLUTION\n");
         else if (parm->meth == GLP_PRIMAL)
            xprintf("PROBLEM HAS UNBOUNDED SOLUTION\n");
         else
            xprintf("PROBLEM HAS NO DUAL FEASIBLE SOLUTION\n");
      }
      return;
}

static int solve_lp(glp_prob *P, const glp_smcp *parm)
{     /* solve LP directly without using the preprocessor */
      int ret;
      if (!glp_bf_exists(P))
      {  ret = glp_factorize(P);
         if (ret == 0)
            ;
         else if (ret == GLP_EBADB)
         {  if (parm->msg_lev >= GLP_MSG_ERR)
               xprintf("glp_simplex: initial basis is invalid\n");
         }
         else if (ret == GLP_ESING)
         {  if (parm->msg_lev >= GLP_MSG_ERR)
               xprintf("glp_simplex: initial basis is singular\n");
         }
         else if (ret == GLP_ECOND)
         {  if (parm->msg_lev >= GLP_MSG_ERR)
               xprintf(
                  "glp_simplex: initial basis is ill-conditioned\n");
         }
         else
            xassert(ret != ret);
         if (ret != 0) goto done;
      }
      if (parm->meth == GLP_PRIMAL)
         ret = spx_primal(P, parm);
      else if (parm->meth == GLP_DUALP)
      {  ret = spx_dual(P, parm);
         if (ret == GLP_EFAIL && P->valid)
            ret = spx_primal(P, parm);
      }
      else if (parm->meth == GLP_DUAL)
         ret = spx_dual(P, parm);
      else
         xassert(parm != parm);
done: return ret;
}

static int preprocess_and_solve_lp(glp_prob *P, const glp_smcp *parm)
{     /* solve LP using the preprocessor */
      NPP *npp;
      glp_prob *lp = NULL;
      glp_bfcp bfcp;
      int ret;
      if (parm->msg_lev >= GLP_MSG_ALL)
         xprintf("Preprocessing...\n");
      /* create preprocessor workspace */
      npp = npp_create_wksp();
      /* load original problem into the preprocessor workspace */
      npp_load_prob(npp, P, GLP_OFF, GLP_SOL, GLP_OFF);
      /* process LP prior to applying primal/dual simplex method */
      ret = npp_simplex(npp, parm);
      if (ret == 0)
         ;
      else if (ret == GLP_ENOPFS)
      {  if (parm->msg_lev >= GLP_MSG_ALL)
            xprintf("PROBLEM HAS NO PRIMAL FEASIBLE SOLUTION\n");
      }
      else if (ret == GLP_ENODFS)
      {  if (parm->msg_lev >= GLP_MSG_ALL)
            xprintf("PROBLEM HAS NO DUAL FEASIBLE SOLUTION\n");
      }
      else
         xassert(ret != ret);
      if (ret != 0) goto done;
      /* build transformed LP */
      lp = glp_create_prob();
      npp_build_prob(npp, lp);
      /* if the transformed LP is empty, it has empty solution, which
         is optimal */
      if (lp->m == 0 && lp->n == 0)
      {  lp->pbs_stat = lp->dbs_stat = GLP_FEAS;
         lp->obj_val = lp->c0;
         if (parm->msg_lev >= GLP_MSG_ON && parm->out_dly == 0)
         {  xprintf("~%6d: obj = %17.9e  infeas = %10.3e\n", P->it_cnt,
               lp->obj_val, 0.0);
         }
         if (parm->msg_lev >= GLP_MSG_ALL)
            xprintf("OPTIMAL SOLUTION FOUND BY LP PREPROCESSOR\n");
         goto post;
      }
      if (parm->msg_lev >= GLP_MSG_ALL)
      {  xprintf("%d row%s, %d column%s, %d non-zero%s\n",
            lp->m, lp->m == 1 ? "" : "s", lp->n, lp->n == 1 ? "" : "s",
            lp->nnz, lp->nnz == 1 ? "" : "s");
      }
      /* inherit basis factorization control parameters */
      glp_get_bfcp(P, &bfcp);
      glp_set_bfcp(lp, &bfcp);
      /* scale the transformed problem */
      {  ENV *env = get_env_ptr();
         int term_out = env->term_out;
         if (!term_out || parm->msg_lev < GLP_MSG_ALL)
            env->term_out = GLP_OFF;
         else
            env->term_out = GLP_ON;
         glp_scale_prob(lp, GLP_SF_AUTO);
         env->term_out = term_out;
      }
      /* build advanced initial basis */
      {  ENV *env = get_env_ptr();
         int term_out = env->term_out;
         if (!term_out || parm->msg_lev < GLP_MSG_ALL)
            env->term_out = GLP_OFF;
         else
            env->term_out = GLP_ON;
         glp_adv_basis(lp, 0);
         env->term_out = term_out;
      }
      /* solve the transformed LP */
      lp->it_cnt = P->it_cnt;
      ret = solve_lp(lp, parm);
      P->it_cnt = lp->it_cnt;
      /* only optimal solution can be postprocessed */
      if (!(ret == 0 && lp->pbs_stat == GLP_FEAS && lp->dbs_stat ==
            GLP_FEAS))
      {  if (parm->msg_lev >= GLP_MSG_ERR)
            xprintf("glp_simplex: unable to recover undefined or non-op"
               "timal solution\n");
         if (ret == 0)
         {  if (lp->pbs_stat == GLP_NOFEAS)
               ret = GLP_ENOPFS;
            else if (lp->dbs_stat == GLP_NOFEAS)
               ret = GLP_ENODFS;
            else
               xassert(lp != lp);
         }
         goto done;
      }
post: /* postprocess solution from the transformed LP */
      npp_postprocess(npp, lp);
      /* the transformed LP is no longer needed */
      glp_delete_prob(lp), lp = NULL;
      /* store solution to the original problem */
      npp_unload_sol(npp, P);
      /* the original LP has been successfully solved */
      ret = 0;
done: /* delete the transformed LP, if it exists */
      if (lp != NULL) glp_delete_prob(lp);
      /* delete preprocessor workspace */
      npp_delete_wksp(npp);
      return ret;
}

int glp_simplex(glp_prob *P, const glp_smcp *parm)
{     /* solve LP problem with the simplex method */
      glp_smcp _parm;
      int i, j, ret;
      /* check problem object */
      if (P == NULL || P->magic != GLP_PROB_MAGIC)
         xerror("glp_simplex: P = %p; invalid problem object\n", P);
      if (P->tree != NULL && P->tree->reason != 0)
         xerror("glp_simplex: operation not allowed\n");
      /* check control parameters */
      if (parm == NULL)
         parm = &_parm, glp_init_smcp((glp_smcp *)parm);
      if (!(parm->msg_lev == GLP_MSG_OFF ||
            parm->msg_lev == GLP_MSG_ERR ||
            parm->msg_lev == GLP_MSG_ON  ||
            parm->msg_lev == GLP_MSG_ALL ||
            parm->msg_lev == GLP_MSG_DBG))
         xerror("glp_simplex: msg_lev = %d; invalid parameter\n",
            parm->msg_lev);
      if (!(parm->meth == GLP_PRIMAL ||
            parm->meth == GLP_DUALP  ||
            parm->meth == GLP_DUAL))
         xerror("glp_simplex: meth = %d; invalid parameter\n",
            parm->meth);
      if (!(parm->pricing == GLP_PT_STD ||
            parm->pricing == GLP_PT_PSE))
         xerror("glp_simplex: pricing = %d; invalid parameter\n",
            parm->pricing);
      if (!(parm->r_test == GLP_RT_STD ||
            parm->r_test == GLP_RT_HAR))
         xerror("glp_simplex: r_test = %d; invalid parameter\n",
            parm->r_test);
      if (!(0.0 < parm->tol_bnd && parm->tol_bnd < 1.0))
         xerror("glp_simplex: tol_bnd = %g; invalid parameter\n",
            parm->tol_bnd);
      if (!(0.0 < parm->tol_dj && parm->tol_dj < 1.0))
         xerror("glp_simplex: tol_dj = %g; invalid parameter\n",
            parm->tol_dj);
      if (!(0.0 < parm->tol_piv && parm->tol_piv < 1.0))
         xerror("glp_simplex: tol_piv = %g; invalid parameter\n",
            parm->tol_piv);
      if (parm->it_lim < 0)
         xerror("glp_simplex: it_lim = %d; invalid parameter\n",
            parm->it_lim);
      if (parm->tm_lim < 0)
         xerror("glp_simplex: tm_lim = %d; invalid parameter\n",
            parm->tm_lim);
      if (parm->out_frq < 1)
         xerror("glp_simplex: out_frq = %d; invalid parameter\n",
            parm->out_frq);
      if (parm->out_dly < 0)
         xerror("glp_simplex: out_dly = %d; invalid parameter\n",
            parm->out_dly);
      if (!(parm->presolve == GLP_ON || parm->presolve == GLP_OFF))
         xerror("glp_simplex: presolve = %d; invalid parameter\n",
            parm->presolve);
      /* basic solution is currently undefined */
      P->pbs_stat = P->dbs_stat = GLP_UNDEF;
      P->obj_val = 0.0;
      P->some = 0;
      /* check bounds of double-bounded variables */
      for (i = 1; i <= P->m; i++)
      {  GLPROW *row = P->row[i];
         if (row->type == GLP_DB && row->lb >= row->ub)
         {  if (parm->msg_lev >= GLP_MSG_ERR)
               xprintf("glp_simplex: row %d: lb = %g, ub = %g; incorrec"
                  "t bounds\n", i, row->lb, row->ub);
            ret = GLP_EBOUND;
            goto done;
         }
      }
      for (j = 1; j <= P->n; j++)
      {  GLPCOL *col = P->col[j];
         if (col->type == GLP_DB && col->lb >= col->ub)
         {  if (parm->msg_lev >= GLP_MSG_ERR)
               xprintf("glp_simplex: column %d: lb = %g, ub = %g; incor"
                  "rect bounds\n", j, col->lb, col->ub);
            ret = GLP_EBOUND;
            goto done;
         }
      }
      /* solve LP problem */
      if (parm->msg_lev >= GLP_MSG_ALL)
      {  xprintf("GLPK Simplex Optimizer, v%s\n", glp_version());
         xprintf("%d row%s, %d column%s, %d non-zero%s\n",
            P->m, P->m == 1 ? "" : "s", P->n, P->n == 1 ? "" : "s",
            P->nnz, P->nnz == 1 ? "" : "s");
      }
      if (P->nnz == 0)
         trivial_lp(P, parm), ret = 0;
      else if (!parm->presolve)
         ret = solve_lp(P, parm);
      else
         ret = preprocess_and_solve_lp(P, parm);
done: /* return to the application program */
      return ret;
}

/***********************************************************************
*  NAME
*
*  glp_init_smcp - initialize simplex method control parameters
*
*  SYNOPSIS
*
*  void glp_init_smcp(glp_smcp *parm);
*
*  DESCRIPTION
*
*  The routine glp_init_smcp initializes control parameters, which are
*  used by the simplex solver, with default values.
*
*  Default values of the control parameters are stored in a glp_smcp
*  structure, which the parameter parm points to. */

void glp_init_smcp(glp_smcp *parm)
{     parm->msg_lev = GLP_MSG_ALL;
      parm->meth = GLP_PRIMAL;
      parm->pricing = GLP_PT_PSE;
      parm->r_test = GLP_RT_HAR;
      parm->tol_bnd = 1e-7;
      parm->tol_dj = 1e-7;
      parm->tol_piv = 1e-10;
      parm->obj_ll = -DBL_MAX;
      parm->obj_ul = +DBL_MAX;
      parm->it_lim = INT_MAX;
      parm->tm_lim = INT_MAX;
      parm->out_frq = 500;
      parm->out_dly = 0;
      parm->presolve = GLP_OFF;
      return;
}

/***********************************************************************
*  NAME
*
*  glp_get_status - retrieve generic status of basic solution
*
*  SYNOPSIS
*
*  int glp_get_status(glp_prob *lp);
*
*  RETURNS
*
*  The routine glp_get_status reports the generic status of the basic
*  solution for the specified problem object as follows:
*
*  GLP_OPT    - solution is optimal;
*  GLP_FEAS   - solution is feasible;
*  GLP_INFEAS - solution is infeasible;
*  GLP_NOFEAS - problem has no feasible solution;
*  GLP_UNBND  - problem has unbounded solution;
*  GLP_UNDEF  - solution is undefined. */

int glp_get_status(glp_prob *lp)
{     int status;
      status = glp_get_prim_stat(lp);
      switch (status)
      {  case GLP_FEAS:
            switch (glp_get_dual_stat(lp))
            {  case GLP_FEAS:
                  status = GLP_OPT;
                  break;
               case GLP_NOFEAS:
                  status = GLP_UNBND;
                  break;
               case GLP_UNDEF:
               case GLP_INFEAS:
                  status = status;
                  break;
               default:
                  xassert(lp != lp);
            }
            break;
         case GLP_UNDEF:
         case GLP_INFEAS:
         case GLP_NOFEAS:
            status = status;
            break;
         default:
            xassert(lp != lp);
      }
      return status;
}

/***********************************************************************
*  NAME
*
*  glp_get_prim_stat - retrieve status of primal basic solution
*
*  SYNOPSIS
*
*  int glp_get_prim_stat(glp_prob *lp);
*
*  RETURNS
*
*  The routine glp_get_prim_stat reports the status of the primal basic
*  solution for the specified problem object as follows:
*
*  GLP_UNDEF  - primal solution is undefined;
*  GLP_FEAS   - primal solution is feasible;
*  GLP_INFEAS - primal solution is infeasible;
*  GLP_NOFEAS - no primal feasible solution exists. */

int glp_get_prim_stat(glp_prob *lp)
{     int pbs_stat = lp->pbs_stat;
      return pbs_stat;
}

/***********************************************************************
*  NAME
*
*  glp_get_dual_stat - retrieve status of dual basic solution
*
*  SYNOPSIS
*
*  int glp_get_dual_stat(glp_prob *lp);
*
*  RETURNS
*
*  The routine glp_get_dual_stat reports the status of the dual basic
*  solution for the specified problem object as follows:
*
*  GLP_UNDEF  - dual solution is undefined;
*  GLP_FEAS   - dual solution is feasible;
*  GLP_INFEAS - dual solution is infeasible;
*  GLP_NOFEAS - no dual feasible solution exists. */

int glp_get_dual_stat(glp_prob *lp)
{     int dbs_stat = lp->dbs_stat;
      return dbs_stat;
}

/***********************************************************************
*  NAME
*
*  glp_get_obj_val - retrieve objective value (basic solution)
*
*  SYNOPSIS
*
*  double glp_get_obj_val(glp_prob *lp);
*
*  RETURNS
*
*  The routine glp_get_obj_val returns value of the objective function
*  for basic solution. */

double glp_get_obj_val(glp_prob *lp)
{     /*struct LPXCPS *cps = lp->cps;*/
      double z;
      z = lp->obj_val;
      /*if (cps->round && fabs(z) < 1e-9) z = 0.0;*/
      return z;
}

/***********************************************************************
*  NAME
*
*  glp_get_row_stat - retrieve row status
*
*  SYNOPSIS
*
*  int glp_get_row_stat(glp_prob *lp, int i);
*
*  RETURNS
*
*  The routine glp_get_row_stat returns current status assigned to the
*  auxiliary variable associated with i-th row as follows:
*
*  GLP_BS - basic variable;
*  GLP_NL - non-basic variable on its lower bound;
*  GLP_NU - non-basic variable on its upper bound;
*  GLP_NF - non-basic free (unbounded) variable;
*  GLP_NS - non-basic fixed variable. */

int glp_get_row_stat(glp_prob *lp, int i)
{     if (!(1 <= i && i <= lp->m))
         xerror("glp_get_row_stat: i = %d; row number out of range\n",
            i);
      return lp->row[i]->stat;
}

/***********************************************************************
*  NAME
*
*  glp_get_row_prim - retrieve row primal value (basic solution)
*
*  SYNOPSIS
*
*  double glp_get_row_prim(glp_prob *lp, int i);
*
*  RETURNS
*
*  The routine glp_get_row_prim returns primal value of the auxiliary
*  variable associated with i-th row. */

double glp_get_row_prim(glp_prob *lp, int i)
{     /*struct LPXCPS *cps = lp->cps;*/
      double prim;
      if (!(1 <= i && i <= lp->m))
         xerror("glp_get_row_prim: i = %d; row number out of range\n",
            i);
      prim = lp->row[i]->prim;
      /*if (cps->round && fabs(prim) < 1e-9) prim = 0.0;*/
      return prim;
}

/***********************************************************************
*  NAME
*
*  glp_get_row_dual - retrieve row dual value (basic solution)
*
*  SYNOPSIS
*
*  double glp_get_row_dual(glp_prob *lp, int i);
*
*  RETURNS
*
*  The routine glp_get_row_dual returns dual value (i.e. reduced cost)
*  of the auxiliary variable associated with i-th row. */

double glp_get_row_dual(glp_prob *lp, int i)
{     /*struct LPXCPS *cps = lp->cps;*/
      double dual;
      if (!(1 <= i && i <= lp->m))
         xerror("glp_get_row_dual: i = %d; row number out of range\n",
            i);
      dual = lp->row[i]->dual;
      /*if (cps->round && fabs(dual) < 1e-9) dual = 0.0;*/
      return dual;
}

/***********************************************************************
*  NAME
*
*  glp_get_col_stat - retrieve column status
*
*  SYNOPSIS
*
*  int glp_get_col_stat(glp_prob *lp, int j);
*
*  RETURNS
*
*  The routine glp_get_col_stat returns current status assigned to the
*  structural variable associated with j-th column as follows:
*
*  GLP_BS - basic variable;
*  GLP_NL - non-basic variable on its lower bound;
*  GLP_NU - non-basic variable on its upper bound;
*  GLP_NF - non-basic free (unbounded) variable;
*  GLP_NS - non-basic fixed variable. */

int glp_get_col_stat(glp_prob *lp, int j)
{     if (!(1 <= j && j <= lp->n))
         xerror("glp_get_col_stat: j = %d; column number out of range\n"
            , j);
      return lp->col[j]->stat;
}

/***********************************************************************
*  NAME
*
*  glp_get_col_prim - retrieve column primal value (basic solution)
*
*  SYNOPSIS
*
*  double glp_get_col_prim(glp_prob *lp, int j);
*
*  RETURNS
*
*  The routine glp_get_col_prim returns primal value of the structural
*  variable associated with j-th column. */

double glp_get_col_prim(glp_prob *lp, int j)
{     /*struct LPXCPS *cps = lp->cps;*/
      double prim;
      if (!(1 <= j && j <= lp->n))
         xerror("glp_get_col_prim: j = %d; column number out of range\n"
            , j);
      prim = lp->col[j]->prim;
      /*if (cps->round && fabs(prim) < 1e-9) prim = 0.0;*/
      return prim;
}

/***********************************************************************
*  NAME
*
*  glp_get_col_dual - retrieve column dual value (basic solution)
*
*  SYNOPSIS
*
*  double glp_get_col_dual(glp_prob *lp, int j);
*
*  RETURNS
*
*  The routine glp_get_col_dual returns dual value (i.e. reduced cost)
*  of the structural variable associated with j-th column. */

double glp_get_col_dual(glp_prob *lp, int j)
{     /*struct LPXCPS *cps = lp->cps;*/
      double dual;
      if (!(1 <= j && j <= lp->n))
         xerror("glp_get_col_dual: j = %d; column number out of range\n"
            , j);
      dual = lp->col[j]->dual;
      /*if (cps->round && fabs(dual) < 1e-9) dual = 0.0;*/
      return dual;
}

/***********************************************************************
*  NAME
*
*  glp_get_unbnd_ray - determine variable causing unboundedness
*
*  SYNOPSIS
*
*  int glp_get_unbnd_ray(glp_prob *lp);
*
*  RETURNS
*
*  The routine glp_get_unbnd_ray returns the number k of a variable,
*  which causes primal or dual unboundedness. If 1 <= k <= m, it is
*  k-th auxiliary variable, and if m+1 <= k <= m+n, it is (k-m)-th
*  structural variable, where m is the number of rows, n is the number
*  of columns in the problem object. If such variable is not defined,
*  the routine returns 0.
*
*  COMMENTS
*
*  If it is not exactly known which version of the simplex solver
*  detected unboundedness, i.e. whether the unboundedness is primal or
*  dual, it is sufficient to check the status of the variable reported
*  with the routine glp_get_row_stat or glp_get_col_stat. If the
*  variable is non-basic, the unboundedness is primal, otherwise, if
*  the variable is basic, the unboundedness is dual (the latter case
*  means that the problem has no primal feasible dolution). */

int glp_get_unbnd_ray(glp_prob *lp)
{     int k;
      k = lp->some;
      xassert(k >= 0);
      if (k > lp->m + lp->n) k = 0;
      return k;
}

/* eof */
