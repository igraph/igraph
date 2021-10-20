/* wrsol.c (write basic solution in GLPK format) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2010-2016 Free Software Foundation, Inc.
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

/***********************************************************************
*  NAME
*
*  glp_write_sol - write basic solution in GLPK format
*
*  SYNOPSIS
*
*  int glp_write_sol(glp_prob *P, const char *fname);
*
*  DESCRIPTION
*
*  The routine glp_write_sol writes basic solution to a text file in
*  GLPK format.
*
*  RETURNS
*
*  If the operation was successful, the routine returns zero. Otherwise
*  it prints an error message and returns non-zero. */

int glp_write_sol(glp_prob *P, const char *fname)
{     glp_file *fp;
      GLPROW *row;
      GLPCOL *col;
      int i, j, count, ret = 1;
      char *s;
#if 0 /* 04/IV-2016 */
      if (P == NULL || P->magic != GLP_PROB_MAGIC)
         xerror("glp_write_sol: P = %p; invalid problem object\n", P);
#endif
      if (fname == NULL)
         xerror("glp_write_sol: fname = %d; invalid parameter\n", fname)
            ;
      xprintf("Writing basic solution to '%s'...\n", fname);
      fp = glp_open(fname, "w"), count = 0;
      if (fp == NULL)
      {  xprintf("Unable to create '%s' - %s\n", fname, get_err_msg());
         goto done;
      }
      /* write comment lines */
      glp_format(fp, "c %-12s%s\n", "Problem:",
         P->name == NULL ? "" : P->name), count++;
      glp_format(fp, "c %-12s%d\n", "Rows:", P->m), count++;
      glp_format(fp, "c %-12s%d\n", "Columns:", P->n), count++;
      glp_format(fp, "c %-12s%d\n", "Non-zeros:", P->nnz), count++;
      switch (glp_get_status(P))
      {  case GLP_OPT:    s = "OPTIMAL";                   break;
         case GLP_FEAS:   s = "FEASIBLE";                  break;
         case GLP_INFEAS: s = "INFEASIBLE (INTERMEDIATE)"; break;
         case GLP_NOFEAS: s = "INFEASIBLE (FINAL)";        break;
         case GLP_UNBND:  s = "UNBOUNDED";                 break;
         case GLP_UNDEF:  s = "UNDEFINED";                 break;
         default:         s = "???";                       break;
      }
      glp_format(fp, "c %-12s%s\n", "Status:", s), count++;
      switch (P->dir)
      {  case GLP_MIN: s = "MINimum"; break;
         case GLP_MAX: s = "MAXimum"; break;
         default:      s = "???";     break;
      }
      glp_format(fp, "c %-12s%s%s%.10g (%s)\n", "Objective:",
         P->obj == NULL ? "" : P->obj,
         P->obj == NULL ? "" : " = ", P->obj_val, s), count++;
      glp_format(fp, "c\n"), count++;
      /* write solution line */
      glp_format(fp, "s bas %d %d ", P->m, P->n), count++;
      switch (P->pbs_stat)
      {  case GLP_UNDEF:  glp_format(fp, "u"); break;
         case GLP_FEAS:   glp_format(fp, "f"); break;
         case GLP_INFEAS: glp_format(fp, "i"); break;
         case GLP_NOFEAS: glp_format(fp, "n"); break;
         default:         glp_format(fp, "?"); break;
      }
      glp_format(fp, " ");
      switch (P->dbs_stat)
      {  case GLP_UNDEF:  glp_format(fp, "u"); break;
         case GLP_FEAS:   glp_format(fp, "f"); break;
         case GLP_INFEAS: glp_format(fp, "i"); break;
         case GLP_NOFEAS: glp_format(fp, "n"); break;
         default:         glp_format(fp, "?"); break;
      }
      glp_format(fp, " %.*g\n", DBL_DIG, P->obj_val);
      /* write row solution descriptor lines */
      for (i = 1; i <= P->m; i++)
      {  row = P->row[i];
         glp_format(fp, "i %d ", i), count++;
         switch (row->stat)
         {  case GLP_BS:
               glp_format(fp, "b");
               break;
            case GLP_NL:
               glp_format(fp, "l");
               break;
            case GLP_NU:
               glp_format(fp, "u");
               break;
            case GLP_NF:
               glp_format(fp, "f");
               break;
            case GLP_NS:
               glp_format(fp, "s");
               break;
            default:
               xassert(row != row);
         }
         glp_format(fp, " %.*g %.*g\n", DBL_DIG, row->prim, DBL_DIG,
            row->dual);
      }
      /* write column solution descriptor lines */
      for (j = 1; j <= P->n; j++)
      {  col = P->col[j];
         glp_format(fp, "j %d ", j), count++;
         switch (col->stat)
         {  case GLP_BS:
               glp_format(fp, "b");
               break;
            case GLP_NL:
               glp_format(fp, "l");
               break;
            case GLP_NU:
               glp_format(fp, "u");
               break;
            case GLP_NF:
               glp_format(fp, "f");
               break;
            case GLP_NS:
               glp_format(fp, "s");
               break;
            default:
               xassert(col != col);
         }
         glp_format(fp, " %.*g %.*g\n", DBL_DIG, col->prim, DBL_DIG,
            col->dual);
      }
      /* write end line */
      glp_format(fp, "e o f\n"), count++;
      if (glp_ioerr(fp))
      {  xprintf("Write error on '%s' - %s\n", fname, get_err_msg());
         goto done;
      }
      /* basic solution has been successfully written */
      xprintf("%d lines were written\n", count);
      ret = 0;
done: if (fp != NULL)
         glp_close(fp);
      return ret;
}

/* eof */
