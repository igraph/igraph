/* wrprob.c (write problem data in GLPK format) */

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

#define xfprintf        glp_format

/***********************************************************************
*  NAME
*
*  glp_write_prob - write problem data in GLPK format
*
*  SYNOPSIS
*
*  int glp_write_prob(glp_prob *P, int flags, const char *fname);
*
*  The routine glp_write_prob writes problem data in GLPK LP/MIP format
*  to a text file.
*
*  RETURNS
*
*  If the operation was successful, the routine returns zero. Otherwise
*  it prints an error message and returns non-zero. */

int glp_write_prob(glp_prob *P, int flags, const char *fname)
{     glp_file *fp;
      GLPROW *row;
      GLPCOL *col;
      GLPAIJ *aij;
      int mip, i, j, count, ret;
#if 0 /* 04/IV-2016 */
      if (P == NULL || P->magic != GLP_PROB_MAGIC)
         xerror("glp_write_prob: P = %p; invalid problem object\n",
            P);
#endif
      if (flags != 0)
         xerror("glp_write_prob: flags = %d; invalid parameter\n",
            flags);
      if (fname == NULL)
         xerror("glp_write_prob: fname = %d; invalid parameter\n",
            fname);
      xprintf("Writing problem data to '%s'...\n", fname);
      fp = glp_open(fname, "w"), count = 0;
      if (fp == NULL)
      {  xprintf("Unable to create '%s' - %s\n", fname, get_err_msg());
         ret = 1;
         goto done;
      }
      /* write problem line */
      mip = (glp_get_num_int(P) > 0);
      xfprintf(fp, "p %s %s %d %d %d\n", !mip ? "lp" : "mip",
         P->dir == GLP_MIN ? "min" : P->dir == GLP_MAX ? "max" : "???",
         P->m, P->n, P->nnz), count++;
      if (P->name != NULL)
         xfprintf(fp, "n p %s\n", P->name), count++;
      if (P->obj != NULL)
         xfprintf(fp, "n z %s\n", P->obj), count++;
      /* write row descriptors */
      for (i = 1; i <= P->m; i++)
      {  row = P->row[i];
         if (row->type == GLP_FX && row->lb == 0.0)
            goto skip1;
         xfprintf(fp, "i %d ", i), count++;
         if (row->type == GLP_FR)
            xfprintf(fp, "f\n");
         else if (row->type == GLP_LO)
            xfprintf(fp, "l %.*g\n", DBL_DIG, row->lb);
         else if (row->type == GLP_UP)
            xfprintf(fp, "u %.*g\n", DBL_DIG, row->ub);
         else if (row->type == GLP_DB)
            xfprintf(fp, "d %.*g %.*g\n", DBL_DIG, row->lb, DBL_DIG,
                  row->ub);
         else if (row->type == GLP_FX)
            xfprintf(fp, "s %.*g\n", DBL_DIG, row->lb);
         else
            xassert(row != row);
skip1:   if (row->name != NULL)
            xfprintf(fp, "n i %d %s\n", i, row->name), count++;
      }
      /* write column descriptors */
      for (j = 1; j <= P->n; j++)
      {  col = P->col[j];
         if (!mip && col->type == GLP_LO && col->lb == 0.0)
            goto skip2;
         if (mip && col->kind == GLP_IV && col->type == GLP_DB &&
             col->lb == 0.0 && col->ub == 1.0)
            goto skip2;
         xfprintf(fp, "j %d ", j), count++;
         if (mip)
         {  if (col->kind == GLP_CV)
               xfprintf(fp, "c ");
            else if (col->kind == GLP_IV)
               xfprintf(fp, "i ");
            else
               xassert(col != col);
         }
         if (col->type == GLP_FR)
            xfprintf(fp, "f\n");
         else if (col->type == GLP_LO)
            xfprintf(fp, "l %.*g\n", DBL_DIG, col->lb);
         else if (col->type == GLP_UP)
            xfprintf(fp, "u %.*g\n", DBL_DIG, col->ub);
         else if (col->type == GLP_DB)
            xfprintf(fp, "d %.*g %.*g\n", DBL_DIG, col->lb, DBL_DIG,
                  col->ub);
         else if (col->type == GLP_FX)
            xfprintf(fp, "s %.*g\n", DBL_DIG, col->lb);
         else
            xassert(col != col);
skip2:   if (col->name != NULL)
            xfprintf(fp, "n j %d %s\n", j, col->name), count++;
      }
      /* write objective coefficient descriptors */
      if (P->c0 != 0.0)
         xfprintf(fp, "a 0 0 %.*g\n", DBL_DIG, P->c0), count++;
      for (j = 1; j <= P->n; j++)
      {  col = P->col[j];
         if (col->coef != 0.0)
            xfprintf(fp, "a 0 %d %.*g\n", j, DBL_DIG, col->coef),
               count++;
      }
      /* write constraint coefficient descriptors */
      for (i = 1; i <= P->m; i++)
      {  row = P->row[i];
         for (aij = row->ptr; aij != NULL; aij = aij->r_next)
            xfprintf(fp, "a %d %d %.*g\n", i, aij->col->j, DBL_DIG,
               aij->val), count++;
      }
      /* write end line */
      xfprintf(fp, "e o f\n"), count++;
#if 0 /* FIXME */
      xfflush(fp);
#endif
      if (glp_ioerr(fp))
      {  xprintf("Write error on '%s' - %s\n", fname, get_err_msg());
         ret = 1;
         goto done;
      }
      xprintf("%d lines were written\n", count);
      ret = 0;
done: if (fp != NULL) glp_close(fp);
      return ret;
}

/* eof */
