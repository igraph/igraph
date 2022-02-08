/* prrngs.c (print sensitivity analysis report) */

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
      glp_file *fp = NULL;
      GLPROW *row;
      GLPCOL *col;
      int m, n, pass, k, t, numb, type, stat, var1, var2, count, page,
         ret;
      double lb, ub, slack, coef, prim, dual, value1, value2, coef1,
         coef2, obj1, obj2;
      const char *name, *limit;
      char buf[13+1];
      /* sanity checks */
#if 0 /* 04/IV-2016 */
      if (P == NULL || P->magic != GLP_PROB_MAGIC)
         xerror("glp_print_ranges: P = %p; invalid problem object\n",
            P);
#endif
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
      xprintf("Write sensitivity analysis report to '%s'...\n", fname);
      fp = glp_open(fname, "w");
      if (fp == NULL)
      {  xprintf("Unable to create '%s' - %s\n", fname, get_err_msg());
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
#if 0 /* FIXME */
      xfflush(fp);
#endif
      if (glp_ioerr(fp))
      {  xprintf("Write error on '%s' - %s\n", fname, get_err_msg());
         ret = 4;
         goto done;
      }
      ret = 0;
done: if (fp != NULL) glp_close(fp);
      return ret;
}

/* eof */
