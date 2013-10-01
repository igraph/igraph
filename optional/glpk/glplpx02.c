/* glplpx02.c */

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

#include "glpapi.h"

/***********************************************************************
*  NAME
*
*  lpx_put_solution - store basic solution components
*
*  SYNOPSIS
*
*  void lpx_put_solution(glp_prob *lp, int inval, const int *p_stat,
*     const int *d_stat, const double *obj_val, const int r_stat[],
*     const double r_prim[], const double r_dual[], const int c_stat[],
*     const double c_prim[], const double c_dual[])
*
*  DESCRIPTION
*
*  The routine lpx_put_solution stores basic solution components to the
*  specified problem object.
*
*  The parameter inval is the basis factorization invalidity flag.
*  If this flag is clear, the current status of the basis factorization
*  remains unchanged. If this flag is set, the routine invalidates the
*  basis factorization.
*
*  The parameter p_stat is a pointer to the status of primal basic
*  solution, which should be specified as follows:
*
*  GLP_UNDEF  - primal solution is undefined;
*  GLP_FEAS   - primal solution is feasible;
*  GLP_INFEAS - primal solution is infeasible;
*  GLP_NOFEAS - no primal feasible solution exists.
*
*  If the parameter p_stat is NULL, the current status of primal basic
*  solution remains unchanged.
*
*  The parameter d_stat is a pointer to the status of dual basic
*  solution, which should be specified as follows:
*
*  GLP_UNDEF  - dual solution is undefined;
*  GLP_FEAS   - dual solution is feasible;
*  GLP_INFEAS - dual solution is infeasible;
*  GLP_NOFEAS - no dual feasible solution exists.
*
*  If the parameter d_stat is NULL, the current status of dual basic
*  solution remains unchanged.
*
*  The parameter obj_val is a pointer to the objective function value.
*  If it is NULL, the current value of the objective function remains
*  unchanged.
*
*  The array element r_stat[i], 1 <= i <= m (where m is the number of
*  rows in the problem object), specifies the status of i-th auxiliary
*  variable, which should be specified as follows:
*
*  GLP_BS - basic variable;
*  GLP_NL - non-basic variable on lower bound;
*  GLP_NU - non-basic variable on upper bound;
*  GLP_NF - non-basic free variable;
*  GLP_NS - non-basic fixed variable.
*
*  If the parameter r_stat is NULL, the current statuses of auxiliary
*  variables remain unchanged.
*
*  The array element r_prim[i], 1 <= i <= m (where m is the number of
*  rows in the problem object), specifies a primal value of i-th
*  auxiliary variable. If the parameter r_prim is NULL, the current
*  primal values of auxiliary variables remain unchanged.
*
*  The array element r_dual[i], 1 <= i <= m (where m is the number of
*  rows in the problem object), specifies a dual value (reduced cost)
*  of i-th auxiliary variable. If the parameter r_dual is NULL, the
*  current dual values of auxiliary variables remain unchanged.
*
*  The array element c_stat[j], 1 <= j <= n (where n is the number of
*  columns in the problem object), specifies the status of j-th
*  structural variable, which should be specified as follows:
*
*  GLP_BS - basic variable;
*  GLP_NL - non-basic variable on lower bound;
*  GLP_NU - non-basic variable on upper bound;
*  GLP_NF - non-basic free variable;
*  GLP_NS - non-basic fixed variable.
*
*  If the parameter c_stat is NULL, the current statuses of structural
*  variables remain unchanged.
*
*  The array element c_prim[j], 1 <= j <= n (where n is the number of
*  columns in the problem object), specifies a primal value of j-th
*  structural variable. If the parameter c_prim is NULL, the current
*  primal values of structural variables remain unchanged.
*
*  The array element c_dual[j], 1 <= j <= n (where n is the number of
*  columns in the problem object), specifies a dual value (reduced cost)
*  of j-th structural variable. If the parameter c_dual is NULL, the
*  current dual values of structural variables remain unchanged. */

void lpx_put_solution(glp_prob *lp, int inval, const int *p_stat,
      const int *d_stat, const double *obj_val, const int r_stat[],
      const double r_prim[], const double r_dual[], const int c_stat[],
      const double c_prim[], const double c_dual[])
{     GLPROW *row;
      GLPCOL *col;
      int i, j;
      /* invalidate the basis factorization, if required */
      if (inval) lp->valid = 0;
      /* store primal status */
      if (p_stat != NULL)
      {  if (!(*p_stat == GLP_UNDEF  || *p_stat == GLP_FEAS ||
               *p_stat == GLP_INFEAS || *p_stat == GLP_NOFEAS))
            xerror("lpx_put_solution: p_stat = %d; invalid primal statu"
               "s\n", *p_stat);
         lp->pbs_stat = *p_stat;
      }
      /* store dual status */
      if (d_stat != NULL)
      {  if (!(*d_stat == GLP_UNDEF  || *d_stat == GLP_FEAS ||
               *d_stat == GLP_INFEAS || *d_stat == GLP_NOFEAS))
            xerror("lpx_put_solution: d_stat = %d; invalid dual status "
               "\n", *d_stat);
         lp->dbs_stat = *d_stat;
      }
      /* store objective function value */
      if (obj_val != NULL) lp->obj_val = *obj_val;
      /* store row solution components */
      for (i = 1; i <= lp->m; i++)
      {  row = lp->row[i];
         if (r_stat != NULL)
         {  if (!(r_stat[i] == GLP_BS ||
                  row->type == GLP_FR && r_stat[i] == GLP_NF ||
                  row->type == GLP_LO && r_stat[i] == GLP_NL ||
                  row->type == GLP_UP && r_stat[i] == GLP_NU ||
                  row->type == GLP_DB && r_stat[i] == GLP_NL ||
                  row->type == GLP_DB && r_stat[i] == GLP_NU ||
                  row->type == GLP_FX && r_stat[i] == GLP_NS))
               xerror("lpx_put_solution: r_stat[%d] = %d; invalid row s"
                  "tatus\n", i, r_stat[i]);
            row->stat = r_stat[i];
         }
         if (r_prim != NULL) row->prim = r_prim[i];
         if (r_dual != NULL) row->dual = r_dual[i];
      }
      /* store column solution components */
      for (j = 1; j <= lp->n; j++)
      {  col = lp->col[j];
         if (c_stat != NULL)
         {  if (!(c_stat[j] == GLP_BS ||
                  col->type == GLP_FR && c_stat[j] == GLP_NF ||
                  col->type == GLP_LO && c_stat[j] == GLP_NL ||
                  col->type == GLP_UP && c_stat[j] == GLP_NU ||
                  col->type == GLP_DB && c_stat[j] == GLP_NL ||
                  col->type == GLP_DB && c_stat[j] == GLP_NU ||
                  col->type == GLP_FX && c_stat[j] == GLP_NS))
               xerror("lpx_put_solution: c_stat[%d] = %d; invalid colum"
                  "n status\n", j, c_stat[j]);
            col->stat = c_stat[j];
         }
         if (c_prim != NULL) col->prim = c_prim[j];
         if (c_dual != NULL) col->dual = c_dual[j];
      }
      return;
}

/*----------------------------------------------------------------------
-- lpx_put_mip_soln - store mixed integer solution components.
--
-- *Synopsis*
--
-- #include "glplpx.h"
-- void lpx_put_mip_soln(glp_prob *lp, int i_stat, double row_mipx[],
--    double col_mipx[]);
--
-- *Description*
--
-- The routine lpx_put_mip_soln stores solution components obtained by
-- branch-and-bound solver into the specified problem object.
--
-- NOTE: This routine is intended for internal use only. */

void lpx_put_mip_soln(glp_prob *lp, int i_stat, double row_mipx[],
      double col_mipx[])
{     GLPROW *row;
      GLPCOL *col;
      int i, j;
      double sum;
      /* store mixed integer status */
#if 0
      if (!(i_stat == LPX_I_UNDEF || i_stat == LPX_I_OPT ||
            i_stat == LPX_I_FEAS  || i_stat == LPX_I_NOFEAS))
         fault("lpx_put_mip_soln: i_stat = %d; invalid mixed integer st"
            "atus", i_stat);
      lp->i_stat = i_stat;
#else
      switch (i_stat)
      {  case LPX_I_UNDEF:
            lp->mip_stat = GLP_UNDEF; break;
         case LPX_I_OPT:
            lp->mip_stat = GLP_OPT;  break;
         case LPX_I_FEAS:
            lp->mip_stat = GLP_FEAS; break;
         case LPX_I_NOFEAS:
            lp->mip_stat = GLP_NOFEAS; break;
         default:
            xerror("lpx_put_mip_soln: i_stat = %d; invalid mixed intege"
               "r status\n", i_stat);
      }
#endif
      /* store row solution components */
      if (row_mipx != NULL)
      {  for (i = 1; i <= lp->m; i++)
         {  row = lp->row[i];
            row->mipx = row_mipx[i];
         }
      }
      /* store column solution components */
      if (col_mipx != NULL)
      {  for (j = 1; j <= lp->n; j++)
         {  col = lp->col[j];
            col->mipx = col_mipx[j];
         }
      }
      /* if the solution is claimed to be integer feasible, check it */
      if (lp->mip_stat == GLP_OPT || lp->mip_stat == GLP_FEAS)
      {  for (j = 1; j <= lp->n; j++)
         {  col = lp->col[j];
            if (col->kind == GLP_IV && col->mipx != floor(col->mipx))
               xerror("lpx_put_mip_soln: col_mipx[%d] = %.*g; must be i"
                  "ntegral\n", j, DBL_DIG, col->mipx);
         }
      }
      /* compute the objective function value */
      sum = lp->c0;
      for (j = 1; j <= lp->n; j++)
      {  col = lp->col[j];
         sum += col->coef * col->mipx;
      }
      lp->mip_obj = sum;
      return;
}

/* eof */
