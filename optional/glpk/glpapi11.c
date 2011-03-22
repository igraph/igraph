/* glpapi11.c (utility routines) */

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

int glp_print_sol(glp_prob *P, const char *fname)
{     /* write basic solution in printable format */
      XFILE *fp;
      GLPROW *row;
      GLPCOL *col;
      int i, j, t, ae_ind, re_ind, ret;
      double ae_max, re_max;
      xprintf("Writing basic solution to `%s'...\n", fname);
      fp = xfopen(fname, "w");
      if (fp == NULL)
      {  xprintf("Unable to create `%s' - %s\n", fname, xerrmsg());
         ret = 1;
         goto done;
      }
      xfprintf(fp, "%-12s%s\n", "Problem:",
         P->name == NULL ? "" : P->name);
      xfprintf(fp, "%-12s%d\n", "Rows:", P->m);
      xfprintf(fp, "%-12s%d\n", "Columns:", P->n);
      xfprintf(fp, "%-12s%d\n", "Non-zeros:", P->nnz);
      t = glp_get_status(P);
      xfprintf(fp, "%-12s%s\n", "Status:",
         t == GLP_OPT    ? "OPTIMAL" :
         t == GLP_FEAS   ? "FEASIBLE" :
         t == GLP_INFEAS ? "INFEASIBLE (INTERMEDIATE)" :
         t == GLP_NOFEAS ? "INFEASIBLE (FINAL)" :
         t == GLP_UNBND  ? "UNBOUNDED" :
         t == GLP_UNDEF  ? "UNDEFINED" : "???");
      xfprintf(fp, "%-12s%s%s%.10g (%s)\n", "Objective:",
         P->obj == NULL ? "" : P->obj,
         P->obj == NULL ? "" : " = ", P->obj_val,
         P->dir == GLP_MIN ? "MINimum" :
         P->dir == GLP_MAX ? "MAXimum" : "???");
      xfprintf(fp, "\n");
      xfprintf(fp, "   No.   Row name   St   Activity     Lower bound  "
         " Upper bound    Marginal\n");
      xfprintf(fp, "------ ------------ -- ------------- ------------- "
         "------------- -------------\n");
      for (i = 1; i <= P->m; i++)
      {  row = P->row[i];
         xfprintf(fp, "%6d ", i);
         if (row->name == NULL || strlen(row->name) <= 12)
            xfprintf(fp, "%-12s ", row->name == NULL ? "" : row->name);
         else
            xfprintf(fp, "%s\n%20s", row->name, "");
         xfprintf(fp, "%s ",
            row->stat == GLP_BS ? "B " :
            row->stat == GLP_NL ? "NL" :
            row->stat == GLP_NU ? "NU" :
            row->stat == GLP_NF ? "NF" :
            row->stat == GLP_NS ? "NS" : "??");
         xfprintf(fp, "%13.6g ",
            fabs(row->prim) <= 1e-9 ? 0.0 : row->prim);
         if (row->type == GLP_LO || row->type == GLP_DB ||
             row->type == GLP_FX)
            xfprintf(fp, "%13.6g ", row->lb);
         else
            xfprintf(fp, "%13s ", "");
         if (row->type == GLP_UP || row->type == GLP_DB)
            xfprintf(fp, "%13.6g ", row->ub);
         else
            xfprintf(fp, "%13s ", row->type == GLP_FX ? "=" : "");
         if (row->stat != GLP_BS)
         {  if (fabs(row->dual) <= 1e-9)
               xfprintf(fp, "%13s", "< eps");
            else
               xfprintf(fp, "%13.6g ", row->dual);
         }
         xfprintf(fp, "\n");
      }
      xfprintf(fp, "\n");
      xfprintf(fp, "   No. Column name  St   Activity     Lower bound  "
         " Upper bound    Marginal\n");
      xfprintf(fp, "------ ------------ -- ------------- ------------- "
         "------------- -------------\n");
      for (j = 1; j <= P->n; j++)
      {  col = P->col[j];
         xfprintf(fp, "%6d ", j);
         if (col->name == NULL || strlen(col->name) <= 12)
            xfprintf(fp, "%-12s ", col->name == NULL ? "" : col->name);
         else
            xfprintf(fp, "%s\n%20s", col->name, "");
         xfprintf(fp, "%s ",
            col->stat == GLP_BS ? "B " :
            col->stat == GLP_NL ? "NL" :
            col->stat == GLP_NU ? "NU" :
            col->stat == GLP_NF ? "NF" :
            col->stat == GLP_NS ? "NS" : "??");
         xfprintf(fp, "%13.6g ",
            fabs(col->prim) <= 1e-9 ? 0.0 : col->prim);
         if (col->type == GLP_LO || col->type == GLP_DB ||
             col->type == GLP_FX)
            xfprintf(fp, "%13.6g ", col->lb);
         else
            xfprintf(fp, "%13s ", "");
         if (col->type == GLP_UP || col->type == GLP_DB)
            xfprintf(fp, "%13.6g ", col->ub);
         else
            xfprintf(fp, "%13s ", col->type == GLP_FX ? "=" : "");
         if (col->stat != GLP_BS)
         {  if (fabs(col->dual) <= 1e-9)
               xfprintf(fp, "%13s", "< eps");
            else
               xfprintf(fp, "%13.6g ", col->dual);
         }
         xfprintf(fp, "\n");
      }
      xfprintf(fp, "\n");
      xfprintf(fp, "Karush-Kuhn-Tucker optimality conditions:\n");
      xfprintf(fp, "\n");
      _glp_check_kkt(P, GLP_SOL, GLP_KKT_PE, &ae_max, &ae_ind, &re_max,
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
      _glp_check_kkt(P, GLP_SOL, GLP_KKT_PB, &ae_max, &ae_ind, &re_max,
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
      _glp_check_kkt(P, GLP_SOL, GLP_KKT_DE, &ae_max, &ae_ind, &re_max,
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
      _glp_check_kkt(P, GLP_SOL, GLP_KKT_DB, &ae_max, &ae_ind, &re_max,
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
      xfflush(fp);
      if (xferror(fp))
      {  xprintf("Write error on `%s' - %s\n", fname, xerrmsg());
         ret = 1;
         goto done;
      }
      ret = 0;
done: if (fp != NULL) xfclose(fp);
      return ret;
}

/***********************************************************************
*  NAME
*
*  glp_read_sol - read basic solution from text file
*
*  SYNOPSIS
*
*  int glp_read_sol(glp_prob *lp, const char *fname);
*
*  DESCRIPTION
*
*  The routine glp_read_sol reads basic solution from a text file whose
*  name is specified by the parameter fname into the problem object.
*
*  For the file format see description of the routine glp_write_sol.
*
*  RETURNS
*
*  On success the routine returns zero, otherwise non-zero. */

int glp_read_sol(glp_prob *lp, const char *fname)
{     glp_data *data;
      jmp_buf jump;
      int i, j, k, ret = 0;
      xprintf("Reading basic solution from `%s'...\n", fname);
      data = glp_sdf_open_file(fname);
      if (data == NULL)
      {  ret = 1;
         goto done;
      }
      if (setjmp(jump))
      {  ret = 1;
         goto done;
      }
      glp_sdf_set_jump(data, jump);
      /* number of rows, number of columns */
      k = glp_sdf_read_int(data);
      if (k != lp->m)
         glp_sdf_error(data, "wrong number of rows\n");
      k = glp_sdf_read_int(data);
      if (k != lp->n)
         glp_sdf_error(data, "wrong number of columns\n");
      /* primal status, dual status, objective value */
      k = glp_sdf_read_int(data);
      if (!(k == GLP_UNDEF || k == GLP_FEAS || k == GLP_INFEAS ||
            k == GLP_NOFEAS))
         glp_sdf_error(data, "invalid primal status\n");
      lp->pbs_stat = k;
      k = glp_sdf_read_int(data);
      if (!(k == GLP_UNDEF || k == GLP_FEAS || k == GLP_INFEAS ||
            k == GLP_NOFEAS))
         glp_sdf_error(data, "invalid dual status\n");
      lp->dbs_stat = k;
      lp->obj_val = glp_sdf_read_num(data);
      /* rows (auxiliary variables) */
      for (i = 1; i <= lp->m; i++)
      {  GLPROW *row = lp->row[i];
         /* status, primal value, dual value */
         k = glp_sdf_read_int(data);
         if (!(k == GLP_BS || k == GLP_NL || k == GLP_NU ||
               k == GLP_NF || k == GLP_NS))
            glp_sdf_error(data, "invalid row status\n");
         glp_set_row_stat(lp, i, k);
         row->prim = glp_sdf_read_num(data);
         row->dual = glp_sdf_read_num(data);
      }
      /* columns (structural variables) */
      for (j = 1; j <= lp->n; j++)
      {  GLPCOL *col = lp->col[j];
         /* status, primal value, dual value */
         k = glp_sdf_read_int(data);
         if (!(k == GLP_BS || k == GLP_NL || k == GLP_NU ||
               k == GLP_NF || k == GLP_NS))
            glp_sdf_error(data, "invalid column status\n");
         glp_set_col_stat(lp, j, k);
         col->prim = glp_sdf_read_num(data);
         col->dual = glp_sdf_read_num(data);
      }
      xprintf("%d lines were read\n", glp_sdf_line(data));
done: if (ret) lp->pbs_stat = lp->dbs_stat = GLP_UNDEF;
      if (data != NULL) glp_sdf_close_file(data);
      return ret;
}

/***********************************************************************
*  NAME
*
*  glp_write_sol - write basic solution to text file
*
*  SYNOPSIS
*
*  int glp_write_sol(glp_prob *lp, const char *fname);
*
*  DESCRIPTION
*
*  The routine glp_write_sol writes the current basic solution to a
*  text file whose name is specified by the parameter fname. This file
*  can be read back with the routine glp_read_sol.
*
*  RETURNS
*
*  On success the routine returns zero, otherwise non-zero.
*
*  FILE FORMAT
*
*  The file created by the routine glp_write_sol is a plain text file,
*  which contains the following information:
*
*     m n
*     p_stat d_stat obj_val
*     r_stat[1] r_prim[1] r_dual[1]
*     . . .
*     r_stat[m] r_prim[m] r_dual[m]
*     c_stat[1] c_prim[1] c_dual[1]
*     . . .
*     c_stat[n] c_prim[n] c_dual[n]
*
*  where:
*  m is the number of rows (auxiliary variables);
*  n is the number of columns (structural variables);
*  p_stat is the primal status of the basic solution (GLP_UNDEF = 1,
*     GLP_FEAS = 2, GLP_INFEAS = 3, or GLP_NOFEAS = 4);
*  d_stat is the dual status of the basic solution (GLP_UNDEF = 1,
*     GLP_FEAS = 2, GLP_INFEAS = 3, or GLP_NOFEAS = 4);
*  obj_val is the objective value;
*  r_stat[i], i = 1,...,m, is the status of i-th row (GLP_BS = 1,
*     GLP_NL = 2, GLP_NU = 3, GLP_NF = 4, or GLP_NS = 5);
*  r_prim[i], i = 1,...,m, is the primal value of i-th row;
*  r_dual[i], i = 1,...,m, is the dual value of i-th row;
*  c_stat[j], j = 1,...,n, is the status of j-th column (GLP_BS = 1,
*     GLP_NL = 2, GLP_NU = 3, GLP_NF = 4, or GLP_NS = 5);
*  c_prim[j], j = 1,...,n, is the primal value of j-th column;
*  c_dual[j], j = 1,...,n, is the dual value of j-th column. */

int glp_write_sol(glp_prob *lp, const char *fname)
{     XFILE *fp;
      int i, j, ret = 0;
      xprintf("Writing basic solution to `%s'...\n", fname);
      fp = xfopen(fname, "w");
      if (fp == NULL)
      {  xprintf("Unable to create `%s' - %s\n", fname, xerrmsg());
         ret = 1;
         goto done;
      }
      /* number of rows, number of columns */
      xfprintf(fp, "%d %d\n", lp->m, lp->n);
      /* primal status, dual status, objective value */
      xfprintf(fp, "%d %d %.*g\n", lp->pbs_stat, lp->dbs_stat, DBL_DIG,
         lp->obj_val);
      /* rows (auxiliary variables) */
      for (i = 1; i <= lp->m; i++)
      {  GLPROW *row = lp->row[i];
         /* status, primal value, dual value */
         xfprintf(fp, "%d %.*g %.*g\n", row->stat, DBL_DIG, row->prim,
            DBL_DIG, row->dual);
      }
      /* columns (structural variables) */
      for (j = 1; j <= lp->n; j++)
      {  GLPCOL *col = lp->col[j];
         /* status, primal value, dual value */
         xfprintf(fp, "%d %.*g %.*g\n", col->stat, DBL_DIG, col->prim,
            DBL_DIG, col->dual);
      }
      xfflush(fp);
      if (xferror(fp))
      {  xprintf("Write error on `%s' - %s\n", fname, xerrmsg());
         ret = 1;
         goto done;
      }
      xprintf("%d lines were written\n", 2 + lp->m + lp->n);
done: if (fp != NULL) xfclose(fp);
      return ret;
}

/**********************************************************************/

static char *format(char buf[13+1], double x)
{     /* format floating-point number in MPS/360-like style */
      if (x == -DBL_MAX)
         strcpy(buf, "         -Inf");
      else if (x == +DBL_MAX)
         strcpy(buf, "         +Inf");
      else if (fabs(x) <= 999999.99998)
      {  sprintf(buf, "%13.5f", x);
#if 1
         if (strcmp(buf, "      0.00000") == 0 ||
             strcmp(buf, "     -0.00000") == 0)
            strcpy(buf, "       .     ");
         else if (memcmp(buf, "      0.", 8) == 0)
            memcpy(buf, "       .", 8);
         else if (memcmp(buf, "     -0.", 8) == 0)
            memcpy(buf, "      -.", 8);
#endif
      }
      else
         sprintf(buf, "%13.6g", x);
      return buf;
}

int glp_print_ranges(glp_prob *P, int len, const int list[],
      int flags, const char *fname)
{     /* print sensitivity analysis report */
      XFILE *fp = NULL;
      GLPROW *row;
      GLPCOL *col;
      int m, n, pass, k, t, numb, type, stat, var1, var2, count, page,
         ret;
      double lb, ub, slack, coef, prim, dual, value1, value2, coef1,
         coef2, obj1, obj2;
      const char *name, *limit;
      char buf[13+1];
      /* sanity checks */
      if (P == NULL || P->magic != GLP_PROB_MAGIC)
         xerror("glp_print_ranges: P = %p; invalid problem object\n",
            P);
      m = P->m, n = P->n;
      if (len < 0)
         xerror("glp_print_ranges: len = %d; invalid list length\n",
            len);
      if (len > 0)
      {  if (list == NULL)
            xerror("glp_print_ranges: list = %p: invalid parameter\n",
               list);
         for (t = 1; t <= len; t++)
         {  k = list[t];
            if (!(1 <= k && k <= m+n))
               xerror("glp_print_ranges: list[%d] = %d; row/column numb"
                  "er out of range\n", t, k);
         }
      }
      if (flags != 0)
         xerror("glp_print_ranges: flags = %d; invalid parameter\n",
            flags);
      if (fname == NULL)
         xerror("glp_print_ranges: fname = %p; invalid parameter\n",
            fname);
      if (glp_get_status(P) != GLP_OPT)
      {  xprintf("glp_print_ranges: optimal basic solution required\n");
         ret = 1;
         goto done;
      }
      if (!glp_bf_exists(P))
      {  xprintf("glp_print_ranges: basis factorization required\n");
         ret = 2;
         goto done;
      }
      /* start reporting */
      xprintf("Write sensitivity analysis report to `%s'...\n", fname);
      fp = xfopen(fname, "w");
      if (fp == NULL)
      {  xprintf("Unable to create `%s' - %s\n", fname, xerrmsg());
         ret = 3;
         goto done;
      }
      page = count = 0;
      for (pass = 1; pass <= 2; pass++)
      for (t = 1; t <= (len == 0 ? m+n : len); t++)
      {  if (t == 1) count = 0;
         k = (len == 0 ? t : list[t]);
         if (pass == 1 && k > m || pass == 2 && k <= m)
            continue;
         if (count == 0)
         {  xfprintf(fp, "GLPK %-4s - SENSITIVITY ANALYSIS REPORT%73sPa"
               "ge%4d\n", glp_version(), "", ++page);
            xfprintf(fp, "\n");
            xfprintf(fp, "%-12s%s\n", "Problem:",
               P->name == NULL ? "" : P->name);
            xfprintf(fp, "%-12s%s%s%.10g (%s)\n", "Objective:",
               P->obj == NULL ? "" : P->obj,
               P->obj == NULL ? "" : " = ", P->obj_val,
               P->dir == GLP_MIN ? "MINimum" :
               P->dir == GLP_MAX ? "MAXimum" : "???");
            xfprintf(fp, "\n");
            xfprintf(fp, "%6s %-12s %2s %13s %13s %13s  %13s %13s %13s "
               "%s\n", "No.", pass == 1 ? "Row name" : "Column name",
               "St", "Activity", pass == 1 ? "Slack" : "Obj coef",
               "Lower bound", "Activity", "Obj coef", "Obj value at",
               "Limiting");
            xfprintf(fp, "%6s %-12s %2s %13s %13s %13s  %13s %13s %13s "
               "%s\n", "", "", "", "", "Marginal", "Upper bound",
               "range", "range", "break point", "variable");
            xfprintf(fp, "------ ------------ -- ------------- --------"
               "----- -------------  ------------- ------------- ------"
               "------- ------------\n");
         }
         if (pass == 1)
         {  numb = k;
            xassert(1 <= numb && numb <= m);
            row = P->row[numb];
            name = row->name;
            type = row->type;
            lb = glp_get_row_lb(P, numb);
            ub = glp_get_row_ub(P, numb);
            coef = 0.0;
            stat = row->stat;
            prim = row->prim;
            if (type == GLP_FR)
               slack = - prim;
            else if (type == GLP_LO)
               slack = lb - prim;
            else if (type == GLP_UP || type == GLP_DB || type == GLP_FX)
               slack = ub - prim;
            dual = row->dual;
         }
         else
         {  numb = k - m;
            xassert(1 <= numb && numb <= n);
            col = P->col[numb];
            name = col->name;
            lb = glp_get_col_lb(P, numb);
            ub = glp_get_col_ub(P, numb);
            coef = col->coef;
            stat = col->stat;
            prim = col->prim;
            slack = 0.0;
            dual = col->dual;
         }
         if (stat != GLP_BS)
         {  glp_analyze_bound(P, k, &value1, &var1, &value2, &var2);
            if (stat == GLP_NF)
               coef1 = coef2 = coef;
            else if (stat == GLP_NS)
               coef1 = -DBL_MAX, coef2 = +DBL_MAX;
            else if (stat == GLP_NL && P->dir == GLP_MIN ||
                     stat == GLP_NU && P->dir == GLP_MAX)
               coef1 = coef - dual, coef2 = +DBL_MAX;
            else
               coef1 = -DBL_MAX, coef2 = coef - dual;
            if (value1 == -DBL_MAX)
            {  if (dual < -1e-9)
                  obj1 = +DBL_MAX;
               else if (dual > +1e-9)
                  obj1 = -DBL_MAX;
               else
                  obj1 = P->obj_val;
            }
            else
               obj1 = P->obj_val + dual * (value1 - prim);
            if (value2 == +DBL_MAX)
            {  if (dual < -1e-9)
                  obj2 = -DBL_MAX;
               else if (dual > +1e-9)
                  obj2 = +DBL_MAX;
               else
                  obj2 = P->obj_val;
            }
            else
               obj2 = P->obj_val + dual * (value2 - prim);
         }
         else
         {  glp_analyze_coef(P, k, &coef1, &var1, &value1, &coef2,
               &var2, &value2);
            if (coef1 == -DBL_MAX)
            {  if (prim < -1e-9)
                  obj1 = +DBL_MAX;
               else if (prim > +1e-9)
                  obj1 = -DBL_MAX;
               else
                  obj1 = P->obj_val;
            }
            else
               obj1 = P->obj_val + (coef1 - coef) * prim;
            if (coef2 == +DBL_MAX)
            {  if (prim < -1e-9)
                  obj2 = -DBL_MAX;
               else if (prim > +1e-9)
                  obj2 = +DBL_MAX;
               else
                  obj2 = P->obj_val;
            }
            else
               obj2 = P->obj_val + (coef2 - coef) * prim;
         }
         /*** first line ***/
         /* row/column number */
         xfprintf(fp, "%6d", numb);
         /* row/column name */
         xfprintf(fp, " %-12.12s", name == NULL ? "" : name);
         if (name != NULL && strlen(name) > 12)
            xfprintf(fp, "%s\n%6s %12s", name+12, "", "");
         /* row/column status */
         xfprintf(fp, " %2s",
            stat == GLP_BS ? "BS" : stat == GLP_NL ? "NL" :
            stat == GLP_NU ? "NU" : stat == GLP_NF ? "NF" :
            stat == GLP_NS ? "NS" : "??");
         /* row/column activity */
         xfprintf(fp, " %s", format(buf, prim));
         /* row slack, column objective coefficient */
         xfprintf(fp, " %s", format(buf, k <= m ? slack : coef));
         /* row/column lower bound */
         xfprintf(fp, " %s", format(buf, lb));
         /* row/column activity range */
         xfprintf(fp, "  %s", format(buf, value1));
         /* row/column objective coefficient range */
         xfprintf(fp, " %s", format(buf, coef1));
         /* objective value at break point */
         xfprintf(fp, " %s", format(buf, obj1));
         /* limiting variable name */
         if (var1 != 0)
         {  if (var1 <= m)
               limit = glp_get_row_name(P, var1);
            else
               limit = glp_get_col_name(P, var1 - m);
            if (limit != NULL)
               xfprintf(fp, " %s", limit);
         }
         xfprintf(fp, "\n");
         /*** second line ***/
         xfprintf(fp, "%6s %-12s %2s %13s", "", "", "", "");
         /* row/column reduced cost */
         xfprintf(fp, " %s", format(buf, dual));
         /* row/column upper bound */
         xfprintf(fp, " %s", format(buf, ub));
         /* row/column activity range */
         xfprintf(fp, "  %s", format(buf, value2));
         /* row/column objective coefficient range */
         xfprintf(fp, " %s", format(buf, coef2));
         /* objective value at break point */
         xfprintf(fp, " %s", format(buf, obj2));
         /* limiting variable name */
         if (var2 != 0)
         {  if (var2 <= m)
               limit = glp_get_row_name(P, var2);
            else
               limit = glp_get_col_name(P, var2 - m);
            if (limit != NULL)
               xfprintf(fp, " %s", limit);
         }
         xfprintf(fp, "\n");
         xfprintf(fp, "\n");
         /* print 10 items per page */
         count = (count + 1) % 10;
      }
      xfprintf(fp, "End of report\n");
      xfflush(fp);
      if (xferror(fp))
      {  xprintf("Write error on `%s' - %s\n", fname, xerrmsg());
         ret = 4;
         goto done;
      }
      ret = 0;
done: if (fp != NULL) xfclose(fp);
      return ret;
}

/**********************************************************************/

int glp_print_ipt(glp_prob *P, const char *fname)
{     /* write interior-point solution in printable format */
      XFILE *fp;
      GLPROW *row;
      GLPCOL *col;
      int i, j, t, ae_ind, re_ind, ret;
      double ae_max, re_max;
      xprintf("Writing interior-point solution to `%s'...\n", fname);
      fp = xfopen(fname, "w");
      if (fp == NULL)
      {  xprintf("Unable to create `%s' - %s\n", fname, xerrmsg());
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
      _glp_check_kkt(P, GLP_IPT, GLP_KKT_PE, &ae_max, &ae_ind, &re_max,
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
      _glp_check_kkt(P, GLP_IPT, GLP_KKT_PB, &ae_max, &ae_ind, &re_max,
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
      _glp_check_kkt(P, GLP_IPT, GLP_KKT_DE, &ae_max, &ae_ind, &re_max,
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
      _glp_check_kkt(P, GLP_IPT, GLP_KKT_DB, &ae_max, &ae_ind, &re_max,
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
      xfflush(fp);
      if (xferror(fp))
      {  xprintf("Write error on `%s' - %s\n", fname, xerrmsg());
         ret = 1;
         goto done;
      }
      ret = 0;
done: if (fp != NULL) xfclose(fp);
      return ret;
}

/***********************************************************************
*  NAME
*
*  glp_read_ipt - read interior-point solution from text file
*
*  SYNOPSIS
*
*  int glp_read_ipt(glp_prob *lp, const char *fname);
*
*  DESCRIPTION
*
*  The routine glp_read_ipt reads interior-point solution from a text
*  file whose name is specified by the parameter fname into the problem
*  object.
*
*  For the file format see description of the routine glp_write_ipt.
*
*  RETURNS
*
*  On success the routine returns zero, otherwise non-zero. */

int glp_read_ipt(glp_prob *lp, const char *fname)
{     glp_data *data;
      jmp_buf jump;
      int i, j, k, ret = 0;
      xprintf("Reading interior-point solution from `%s'...\n", fname);
      data = glp_sdf_open_file(fname);
      if (data == NULL)
      {  ret = 1;
         goto done;
      }
      if (setjmp(jump))
      {  ret = 1;
         goto done;
      }
      glp_sdf_set_jump(data, jump);
      /* number of rows, number of columns */
      k = glp_sdf_read_int(data);
      if (k != lp->m)
         glp_sdf_error(data, "wrong number of rows\n");
      k = glp_sdf_read_int(data);
      if (k != lp->n)
         glp_sdf_error(data, "wrong number of columns\n");
      /* solution status, objective value */
      k = glp_sdf_read_int(data);
      if (!(k == GLP_UNDEF || k == GLP_OPT))
         glp_sdf_error(data, "invalid solution status\n");
      lp->ipt_stat = k;
      lp->ipt_obj = glp_sdf_read_num(data);
      /* rows (auxiliary variables) */
      for (i = 1; i <= lp->m; i++)
      {  GLPROW *row = lp->row[i];
         /* primal value, dual value */
         row->pval = glp_sdf_read_num(data);
         row->dval = glp_sdf_read_num(data);
      }
      /* columns (structural variables) */
      for (j = 1; j <= lp->n; j++)
      {  GLPCOL *col = lp->col[j];
         /* primal value, dual value */
         col->pval = glp_sdf_read_num(data);
         col->dval = glp_sdf_read_num(data);
      }
      xprintf("%d lines were read\n", glp_sdf_line(data));
done: if (ret) lp->ipt_stat = GLP_UNDEF;
      if (data != NULL) glp_sdf_close_file(data);
      return ret;
}

/***********************************************************************
*  NAME
*
*  glp_write_ipt - write interior-point solution to text file
*
*  SYNOPSIS
*
*  int glp_write_ipt(glp_prob *lp, const char *fname);
*
*  DESCRIPTION
*
*  The routine glp_write_ipt writes the current interior-point solution
*  to a text file whose name is specified by the parameter fname. This
*  file can be read back with the routine glp_read_ipt.
*
*  RETURNS
*
*  On success the routine returns zero, otherwise non-zero.
*
*  FILE FORMAT
*
*  The file created by the routine glp_write_ipt is a plain text file,
*  which contains the following information:
*
*     m n
*     stat obj_val
*     r_prim[1] r_dual[1]
*     . . .
*     r_prim[m] r_dual[m]
*     c_prim[1] c_dual[1]
*     . . .
*     c_prim[n] c_dual[n]
*
*  where:
*  m is the number of rows (auxiliary variables);
*  n is the number of columns (structural variables);
*  stat is the solution status (GLP_UNDEF = 1 or GLP_OPT = 5);
*  obj_val is the objective value;
*  r_prim[i], i = 1,...,m, is the primal value of i-th row;
*  r_dual[i], i = 1,...,m, is the dual value of i-th row;
*  c_prim[j], j = 1,...,n, is the primal value of j-th column;
*  c_dual[j], j = 1,...,n, is the dual value of j-th column. */

int glp_write_ipt(glp_prob *lp, const char *fname)
{     XFILE *fp;
      int i, j, ret = 0;
      xprintf("Writing interior-point solution to `%s'...\n", fname);
      fp = xfopen(fname, "w");
      if (fp == NULL)
      {  xprintf("Unable to create `%s' - %s\n", fname, xerrmsg());
         ret = 1;
         goto done;
      }
      /* number of rows, number of columns */
      xfprintf(fp, "%d %d\n", lp->m, lp->n);
      /* solution status, objective value */
      xfprintf(fp, "%d %.*g\n", lp->ipt_stat, DBL_DIG, lp->ipt_obj);
      /* rows (auxiliary variables) */
      for (i = 1; i <= lp->m; i++)
      {  GLPROW *row = lp->row[i];
         /* primal value, dual value */
         xfprintf(fp, "%.*g %.*g\n", DBL_DIG, row->pval, DBL_DIG,
            row->dval);
      }
      /* columns (structural variables) */
      for (j = 1; j <= lp->n; j++)
      {  GLPCOL *col = lp->col[j];
         /* primal value, dual value */
         xfprintf(fp, "%.*g %.*g\n", DBL_DIG, col->pval, DBL_DIG,
            col->dval);
      }
      xfflush(fp);
      if (xferror(fp))
      {  xprintf("Write error on `%s' - %s\n", fname, xerrmsg());
         ret = 1;
         goto done;
      }
      xprintf("%d lines were written\n", 2 + lp->m + lp->n);
done: if (fp != NULL) xfclose(fp);
      return ret;
}

/**********************************************************************/

int glp_print_mip(glp_prob *P, const char *fname)
{     /* write MIP solution in printable format */
      XFILE *fp;
      GLPROW *row;
      GLPCOL *col;
      int i, j, t, ae_ind, re_ind, ret;
      double ae_max, re_max;
      xprintf("Writing MIP solution to `%s'...\n", fname);
      fp = xfopen(fname, "w");
      if (fp == NULL)
      {  xprintf("Unable to create `%s' - %s\n", fname, xerrmsg());
         ret = 1;
         goto done;
      }
      xfprintf(fp, "%-12s%s\n", "Problem:",
         P->name == NULL ? "" : P->name);
      xfprintf(fp, "%-12s%d\n", "Rows:", P->m);
      xfprintf(fp, "%-12s%d (%d integer, %d binary)\n", "Columns:",
         P->n, glp_get_num_int(P), glp_get_num_bin(P));
      xfprintf(fp, "%-12s%d\n", "Non-zeros:", P->nnz);
      t = glp_mip_status(P);
      xfprintf(fp, "%-12s%s\n", "Status:",
         t == GLP_OPT    ? "INTEGER OPTIMAL" :
         t == GLP_FEAS   ? "INTEGER NON-OPTIMAL" :
         t == GLP_NOFEAS ? "INTEGER EMPTY" :
         t == GLP_UNDEF  ? "INTEGER UNDEFINED" : "???");
      xfprintf(fp, "%-12s%s%s%.10g (%s)\n", "Objective:",
         P->obj == NULL ? "" : P->obj,
         P->obj == NULL ? "" : " = ", P->mip_obj,
         P->dir == GLP_MIN ? "MINimum" :
         P->dir == GLP_MAX ? "MAXimum" : "???");
      xfprintf(fp, "\n");
      xfprintf(fp, "   No.   Row name        Activity     Lower bound  "
         " Upper bound\n");
      xfprintf(fp, "------ ------------    ------------- ------------- "
         "-------------\n");
      for (i = 1; i <= P->m; i++)
      {  row = P->row[i];
         xfprintf(fp, "%6d ", i);
         if (row->name == NULL || strlen(row->name) <= 12)
            xfprintf(fp, "%-12s ", row->name == NULL ? "" : row->name);
         else
            xfprintf(fp, "%s\n%20s", row->name, "");
         xfprintf(fp, "%3s", "");
         xfprintf(fp, "%13.6g ",
            fabs(row->mipx) <= 1e-9 ? 0.0 : row->mipx);
         if (row->type == GLP_LO || row->type == GLP_DB ||
             row->type == GLP_FX)
            xfprintf(fp, "%13.6g ", row->lb);
         else
            xfprintf(fp, "%13s ", "");
         if (row->type == GLP_UP || row->type == GLP_DB)
            xfprintf(fp, "%13.6g ", row->ub);
         else
            xfprintf(fp, "%13s ", row->type == GLP_FX ? "=" : "");
         xfprintf(fp, "\n");
      }
      xfprintf(fp, "\n");
      xfprintf(fp, "   No. Column name       Activity     Lower bound  "
         " Upper bound\n");
      xfprintf(fp, "------ ------------    ------------- ------------- "
         "-------------\n");
      for (j = 1; j <= P->n; j++)
      {  col = P->col[j];
         xfprintf(fp, "%6d ", j);
         if (col->name == NULL || strlen(col->name) <= 12)
            xfprintf(fp, "%-12s ", col->name == NULL ? "" : col->name);
         else
            xfprintf(fp, "%s\n%20s", col->name, "");
         xfprintf(fp, "%s  ",
            col->kind == GLP_CV ? " " :
            col->kind == GLP_IV ? "*" : "?");
         xfprintf(fp, "%13.6g ",
            fabs(col->mipx) <= 1e-9 ? 0.0 : col->mipx);
         if (col->type == GLP_LO || col->type == GLP_DB ||
             col->type == GLP_FX)
            xfprintf(fp, "%13.6g ", col->lb);
         else
            xfprintf(fp, "%13s ", "");
         if (col->type == GLP_UP || col->type == GLP_DB)
            xfprintf(fp, "%13.6g ", col->ub);
         else
            xfprintf(fp, "%13s ", col->type == GLP_FX ? "=" : "");
         xfprintf(fp, "\n");
      }
      xfprintf(fp, "\n");
      xfprintf(fp, "Integer feasibility conditions:\n");
      xfprintf(fp, "\n");
      _glp_check_kkt(P, GLP_MIP, GLP_KKT_PE, &ae_max, &ae_ind, &re_max,
         &re_ind);
      xfprintf(fp, "KKT.PE: max.abs.err = %.2e on row %d\n",
         ae_max, ae_ind);
      xfprintf(fp, "        max.rel.err = %.2e on row %d\n",
         re_max, re_ind);
      xfprintf(fp, "%8s%s\n", "",
         re_max <= 1e-9 ? "High quality" :
         re_max <= 1e-6 ? "Medium quality" :
         re_max <= 1e-3 ? "Low quality" : "SOLUTION IS WRONG");
      xfprintf(fp, "\n");
      _glp_check_kkt(P, GLP_MIP, GLP_KKT_PB, &ae_max, &ae_ind, &re_max,
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
         re_max <= 1e-3 ? "Low quality" : "SOLUTION IS INFEASIBLE");
      xfprintf(fp, "\n");
      xfprintf(fp, "End of output\n");
      xfflush(fp);
      if (xferror(fp))
      {  xprintf("Write error on `%s' - %s\n", fname, xerrmsg());
         ret = 1;
         goto done;
      }
      ret = 0;
done: if (fp != NULL) xfclose(fp);
      return ret;
}

/***********************************************************************
*  NAME
*
*  glp_read_mip - read MIP solution from text file
*
*  SYNOPSIS
*
*  int glp_read_mip(glp_prob *mip, const char *fname);
*
*  DESCRIPTION
*
*  The routine glp_read_mip reads MIP solution from a text file whose
*  name is specified by the parameter fname into the problem object.
*
*  For the file format see description of the routine glp_write_mip.
*
*  RETURNS
*
*  On success the routine returns zero, otherwise non-zero. */

int glp_read_mip(glp_prob *mip, const char *fname)
{     glp_data *data;
      jmp_buf jump;
      int i, j, k, ret = 0;
      xprintf("Reading MIP solution from `%s'...\n", fname);
      data = glp_sdf_open_file(fname);
      if (data == NULL)
      {  ret = 1;
         goto done;
      }
      if (setjmp(jump))
      {  ret = 1;
         goto done;
      }
      glp_sdf_set_jump(data, jump);
      /* number of rows, number of columns */
      k = glp_sdf_read_int(data);
      if (k != mip->m)
         glp_sdf_error(data, "wrong number of rows\n");
      k = glp_sdf_read_int(data);
      if (k != mip->n)
         glp_sdf_error(data, "wrong number of columns\n");
      /* solution status, objective value */
      k = glp_sdf_read_int(data);
      if (!(k == GLP_UNDEF || k == GLP_OPT || k == GLP_FEAS ||
            k == GLP_NOFEAS))
         glp_sdf_error(data, "invalid solution status\n");
      mip->mip_stat = k;
      mip->mip_obj = glp_sdf_read_num(data);
      /* rows (auxiliary variables) */
      for (i = 1; i <= mip->m; i++)
      {  GLPROW *row = mip->row[i];
         row->mipx = glp_sdf_read_num(data);
      }
      /* columns (structural variables) */
      for (j = 1; j <= mip->n; j++)
      {  GLPCOL *col = mip->col[j];
         col->mipx = glp_sdf_read_num(data);
         if (col->kind == GLP_IV && col->mipx != floor(col->mipx))
            glp_sdf_error(data, "non-integer column value");
      }
      xprintf("%d lines were read\n", glp_sdf_line(data));
done: if (ret) mip->mip_stat = GLP_UNDEF;
      if (data != NULL) glp_sdf_close_file(data);
      return ret;
}

/***********************************************************************
*  NAME
*
*  glp_write_mip - write MIP solution to text file
*
*  SYNOPSIS
*
*  int glp_write_mip(glp_prob *mip, const char *fname);
*
*  DESCRIPTION
*
*  The routine glp_write_mip writes the current MIP solution to a text
*  file whose name is specified by the parameter fname. This file can
*  be read back with the routine glp_read_mip.
*
*  RETURNS
*
*  On success the routine returns zero, otherwise non-zero.
*
*  FILE FORMAT
*
*  The file created by the routine glp_write_sol is a plain text file,
*  which contains the following information:
*
*     m n
*     stat obj_val
*     r_val[1]
*     . . .
*     r_val[m]
*     c_val[1]
*     . . .
*     c_val[n]
*
*  where:
*  m is the number of rows (auxiliary variables);
*  n is the number of columns (structural variables);
*  stat is the solution status (GLP_UNDEF = 1, GLP_FEAS = 2,
*     GLP_NOFEAS = 4, or GLP_OPT = 5);
*  obj_val is the objective value;
*  r_val[i], i = 1,...,m, is the value of i-th row;
*  c_val[j], j = 1,...,n, is the value of j-th column. */

int glp_write_mip(glp_prob *mip, const char *fname)
{     XFILE *fp;
      int i, j, ret = 0;
      xprintf("Writing MIP solution to `%s'...\n", fname);
      fp = xfopen(fname, "w");
      if (fp == NULL)
      {  xprintf("Unable to create `%s' - %s\n", fname, xerrmsg());
         ret = 1;
         goto done;
      }
      /* number of rows, number of columns */
      xfprintf(fp, "%d %d\n", mip->m, mip->n);
      /* solution status, objective value */
      xfprintf(fp, "%d %.*g\n", mip->mip_stat, DBL_DIG, mip->mip_obj);
      /* rows (auxiliary variables) */
      for (i = 1; i <= mip->m; i++)
         xfprintf(fp, "%.*g\n", DBL_DIG, mip->row[i]->mipx);
      /* columns (structural variables) */
      for (j = 1; j <= mip->n; j++)
         xfprintf(fp, "%.*g\n", DBL_DIG, mip->col[j]->mipx);
      xfflush(fp);
      if (xferror(fp))
      {  xprintf("Write error on `%s' - %s\n", fname, xerrmsg());
         ret = 1;
         goto done;
      }
      xprintf("%d lines were written\n", 2 + mip->m + mip->n);
done: if (fp != NULL) xfclose(fp);
      return ret;
}

/* eof */
