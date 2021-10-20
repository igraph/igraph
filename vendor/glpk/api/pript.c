/* pript.c (write interior-point solution in printable format) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2009-2016 Free Software Foundation, Inc.
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
#include "prob.h"

#define xfprintf glp_format

int glp_print_ipt(glp_prob *P, const char *fname)
{     /* write interior-point solution in printable format */
      glp_file *fp;
      GLPROW *row;
      GLPCOL *col;
      int i, j, t, ae_ind, re_ind, ret;
      double ae_max, re_max;
      xprintf("Writing interior-point solution to '%s'...\n", fname);
      fp = glp_open(fname, "w");
      if (fp == NULL)
      {  xprintf("Unable to create '%s' - %s\n", fname, get_err_msg());
         ret = 1;
         goto done;
      }
      xfprintf(fp, "%-12s%s\n", "Problem:",
         P->name == NULL ? "" : P->name);
      xfprintf(fp, "%-12s%d\n", "Rows:", P->m);
      xfprintf(fp, "%-12s%d\n", "Columns:", P->n);
      xfprintf(fp, "%-12s%d\n", "Non-zeros:", P->nnz);
      t = glp_ipt_status(P);
      xfprintf(fp, "%-12s%s\n", "Status:",
         t == GLP_OPT    ? "OPTIMAL" :
         t == GLP_UNDEF  ? "UNDEFINED" :
         t == GLP_INFEAS ? "INFEASIBLE (INTERMEDIATE)" :
         t == GLP_NOFEAS ? "INFEASIBLE (FINAL)" : "???");
      xfprintf(fp, "%-12s%s%s%.10g (%s)\n", "Objective:",
         P->obj == NULL ? "" : P->obj,
         P->obj == NULL ? "" : " = ", P->ipt_obj,
         P->dir == GLP_MIN ? "MINimum" :
         P->dir == GLP_MAX ? "MAXimum" : "???");
      xfprintf(fp, "\n");
      xfprintf(fp, "   No.   Row name        Activity     Lower bound  "
         " Upper bound    Marginal\n");
      xfprintf(fp, "------ ------------    ------------- ------------- "
         "------------- -------------\n");
      for (i = 1; i <= P->m; i++)
      {  row = P->row[i];
         xfprintf(fp, "%6d ", i);
         if (row->name == NULL || strlen(row->name) <= 12)
            xfprintf(fp, "%-12s ", row->name == NULL ? "" : row->name);
         else
            xfprintf(fp, "%s\n%20s", row->name, "");
         xfprintf(fp, "%3s", "");
         xfprintf(fp, "%13.6g ",
            fabs(row->pval) <= 1e-9 ? 0.0 : row->pval);
         if (row->type == GLP_LO || row->type == GLP_DB ||
             row->type == GLP_FX)
            xfprintf(fp, "%13.6g ", row->lb);
         else
            xfprintf(fp, "%13s ", "");
         if (row->type == GLP_UP || row->type == GLP_DB)
            xfprintf(fp, "%13.6g ", row->ub);
         else
            xfprintf(fp, "%13s ", row->type == GLP_FX ? "=" : "");
         if (fabs(row->dval) <= 1e-9)
            xfprintf(fp, "%13s", "< eps");
         else
            xfprintf(fp, "%13.6g ", row->dval);
         xfprintf(fp, "\n");
      }
      xfprintf(fp, "\n");
      xfprintf(fp, "   No. Column name       Activity     Lower bound  "
         " Upper bound    Marginal\n");
      xfprintf(fp, "------ ------------    ------------- ------------- "
         "------------- -------------\n");
      for (j = 1; j <= P->n; j++)
      {  col = P->col[j];
         xfprintf(fp, "%6d ", j);
         if (col->name == NULL || strlen(col->name) <= 12)
            xfprintf(fp, "%-12s ", col->name == NULL ? "" : col->name);
         else
            xfprintf(fp, "%s\n%20s", col->name, "");
         xfprintf(fp, "%3s", "");
         xfprintf(fp, "%13.6g ",
            fabs(col->pval) <= 1e-9 ? 0.0 : col->pval);
         if (col->type == GLP_LO || col->type == GLP_DB ||
             col->type == GLP_FX)
            xfprintf(fp, "%13.6g ", col->lb);
         else
            xfprintf(fp, "%13s ", "");
         if (col->type == GLP_UP || col->type == GLP_DB)
            xfprintf(fp, "%13.6g ", col->ub);
         else
            xfprintf(fp, "%13s ", col->type == GLP_FX ? "=" : "");
         if (fabs(col->dval) <= 1e-9)
            xfprintf(fp, "%13s", "< eps");
         else
            xfprintf(fp, "%13.6g ", col->dval);
         xfprintf(fp, "\n");
      }
      xfprintf(fp, "\n");
      xfprintf(fp, "Karush-Kuhn-Tucker optimality conditions:\n");
      xfprintf(fp, "\n");
      glp_check_kkt(P, GLP_IPT, GLP_KKT_PE, &ae_max, &ae_ind, &re_max,
         &re_ind);
      xfprintf(fp, "KKT.PE: max.abs.err = %.2e on row %d\n",
         ae_max, ae_ind);
      xfprintf(fp, "        max.rel.err = %.2e on row %d\n",
         re_max, re_ind);
      xfprintf(fp, "%8s%s\n", "",
         re_max <= 1e-9 ? "High quality" :
         re_max <= 1e-6 ? "Medium quality" :
         re_max <= 1e-3 ? "Low quality" : "PRIMAL SOLUTION IS WRONG");
      xfprintf(fp, "\n");
      glp_check_kkt(P, GLP_IPT, GLP_KKT_PB, &ae_max, &ae_ind, &re_max,
         &re_ind);
      xfprintf(fp, "KKT.PB: max.abs.err = %.2e on %s %d\n",
            ae_max, ae_ind <= P->m ? "row" : "column",
            ae_ind <= P->m ? ae_ind : ae_ind - P->m);
      xfprintf(fp, "        max.rel.err = %.2e on %s %d\n",
            re_max, re_ind <= P->m ? "row" : "column",
            re_ind <= P->m ? re_ind : re_ind - P->m);
      xfprintf(fp, "%8s%s\n", "",
         re_max <= 1e-9 ? "High quality" :
         re_max <= 1e-6 ? "Medium quality" :
         re_max <= 1e-3 ? "Low quality" : "PRIMAL SOLUTION IS INFEASIBL"
            "E");
      xfprintf(fp, "\n");
      glp_check_kkt(P, GLP_IPT, GLP_KKT_DE, &ae_max, &ae_ind, &re_max,
         &re_ind);
      xfprintf(fp, "KKT.DE: max.abs.err = %.2e on column %d\n",
         ae_max, ae_ind == 0 ? 0 : ae_ind - P->m);
      xfprintf(fp, "        max.rel.err = %.2e on column %d\n",
         re_max, re_ind == 0 ? 0 : re_ind - P->m);
      xfprintf(fp, "%8s%s\n", "",
         re_max <= 1e-9 ? "High quality" :
         re_max <= 1e-6 ? "Medium quality" :
         re_max <= 1e-3 ? "Low quality" : "DUAL SOLUTION IS WRONG");
      xfprintf(fp, "\n");
      glp_check_kkt(P, GLP_IPT, GLP_KKT_DB, &ae_max, &ae_ind, &re_max,
         &re_ind);
      xfprintf(fp, "KKT.DB: max.abs.err = %.2e on %s %d\n",
            ae_max, ae_ind <= P->m ? "row" : "column",
            ae_ind <= P->m ? ae_ind : ae_ind - P->m);
      xfprintf(fp, "        max.rel.err = %.2e on %s %d\n",
            re_max, re_ind <= P->m ? "row" : "column",
            re_ind <= P->m ? re_ind : re_ind - P->m);
      xfprintf(fp, "%8s%s\n", "",
         re_max <= 1e-9 ? "High quality" :
         re_max <= 1e-6 ? "Medium quality" :
         re_max <= 1e-3 ? "Low quality" : "DUAL SOLUTION IS INFEASIBLE")
            ;
      xfprintf(fp, "\n");
      xfprintf(fp, "End of output\n");
#if 0 /* FIXME */
      xfflush(fp);
#endif
      if (glp_ioerr(fp))
      {  xprintf("Write error on '%s' - %s\n", fname, get_err_msg());
         ret = 1;
         goto done;
      }
      ret = 0;
done: if (fp != NULL) glp_close(fp);
      return ret;
}

/* eof */
