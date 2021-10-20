/* rdprob.c (read problem data in GLPK format) */

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

#include "dimacs.h"
#include "misc.h"
#include "prob.h"

#define xfprintf        glp_format
#define error           dmx_error
#define warning         dmx_warning
#define read_char       dmx_read_char
#define read_designator dmx_read_designator
#define read_field      dmx_read_field
#define end_of_line     dmx_end_of_line
#define check_int       dmx_check_int

/***********************************************************************
*  NAME
*
*  glp_read_prob - read problem data in GLPK format
*
*  SYNOPSIS
*
*  int glp_read_prob(glp_prob *P, int flags, const char *fname);
*
*  The routine glp_read_prob reads problem data in GLPK LP/MIP format
*  from a text file.
*
*  RETURNS
*
*  If the operation was successful, the routine returns zero. Otherwise
*  it prints an error message and returns non-zero. */

int glp_read_prob(glp_prob *P, int flags, const char *fname)
{     DMX _csa, *csa = &_csa;
      int mip, m, n, nnz, ne, i, j, k, type, kind, ret, *ln = NULL,
         *ia = NULL, *ja = NULL;
      double lb, ub, temp, *ar = NULL;
      char *rf = NULL, *cf = NULL;
#if 0 /* 04/IV-2016 */
      if (P == NULL || P->magic != GLP_PROB_MAGIC)
         xerror("glp_read_prob: P = %p; invalid problem object\n",
            P);
#endif
      if (flags != 0)
         xerror("glp_read_prob: flags = %d; invalid parameter\n",
            flags);
      if (fname == NULL)
         xerror("glp_read_prob: fname = %d; invalid parameter\n",
            fname);
      glp_erase_prob(P);
      if (setjmp(csa->jump))
      {  ret = 1;
         goto done;
      }
      csa->fname = fname;
      csa->fp = NULL;
      csa->count = 0;
      csa->c = '\n';
      csa->field[0] = '\0';
      csa->empty = csa->nonint = 0;
      xprintf("Reading problem data from '%s'...\n", fname);
      csa->fp = glp_open(fname, "r");
      if (csa->fp == NULL)
      {  xprintf("Unable to open '%s' - %s\n", fname, get_err_msg());
         longjmp(csa->jump, 1);
      }
      /* read problem line */
      read_designator(csa);
      if (strcmp(csa->field, "p") != 0)
         error(csa, "problem line missing or invalid");
      read_field(csa);
      if (strcmp(csa->field, "lp") == 0)
         mip = 0;
      else if (strcmp(csa->field, "mip") == 0)
         mip = 1;
      else
         error(csa, "wrong problem designator; 'lp' or 'mip' expected");
      read_field(csa);
      if (strcmp(csa->field, "min") == 0)
         glp_set_obj_dir(P, GLP_MIN);
      else if (strcmp(csa->field, "max") == 0)
         glp_set_obj_dir(P, GLP_MAX);
      else
         error(csa, "objective sense missing or invalid");
      read_field(csa);
      if (!(str2int(csa->field, &m) == 0 && m >= 0))
         error(csa, "number of rows missing or invalid");
      read_field(csa);
      if (!(str2int(csa->field, &n) == 0 && n >= 0))
         error(csa, "number of columns missing or invalid");
      read_field(csa);
      if (!(str2int(csa->field, &nnz) == 0 && nnz >= 0))
         error(csa, "number of constraint coefficients missing or inval"
            "id");
      if (m > 0)
      {  glp_add_rows(P, m);
         for (i = 1; i <= m; i++)
            glp_set_row_bnds(P, i, GLP_FX, 0.0, 0.0);
      }
      if (n > 0)
      {  glp_add_cols(P, n);
         for (j = 1; j <= n; j++)
         {  if (!mip)
               glp_set_col_bnds(P, j, GLP_LO, 0.0, 0.0);
            else
               glp_set_col_kind(P, j, GLP_BV);
         }
      }
      end_of_line(csa);
      /* allocate working arrays */
      rf = xcalloc(1+m, sizeof(char));
      memset(rf, 0, 1+m);
      cf = xcalloc(1+n, sizeof(char));
      memset(cf, 0, 1+n);
      ln = xcalloc(1+nnz, sizeof(int));
      ia = xcalloc(1+nnz, sizeof(int));
      ja = xcalloc(1+nnz, sizeof(int));
      ar = xcalloc(1+nnz, sizeof(double));
      /* read descriptor lines */
      ne = 0;
      for (;;)
      {  read_designator(csa);
         if (strcmp(csa->field, "i") == 0)
         {  /* row descriptor */
            read_field(csa);
            if (str2int(csa->field, &i) != 0)
               error(csa, "row number missing or invalid");
            if (!(1 <= i && i <= m))
               error(csa, "row number out of range");
            read_field(csa);
            if (strcmp(csa->field, "f") == 0)
               type = GLP_FR;
            else if (strcmp(csa->field, "l") == 0)
               type = GLP_LO;
            else if (strcmp(csa->field, "u") == 0)
               type = GLP_UP;
            else if (strcmp(csa->field, "d") == 0)
               type = GLP_DB;
            else if (strcmp(csa->field, "s") == 0)
               type = GLP_FX;
            else
               error(csa, "row type missing or invalid");
            if (type == GLP_LO || type == GLP_DB || type == GLP_FX)
            {  read_field(csa);
               if (str2num(csa->field, &lb) != 0)
                  error(csa, "row lower bound/fixed value missing or in"
                     "valid");
            }
            else
               lb = 0.0;
            if (type == GLP_UP || type == GLP_DB)
            {  read_field(csa);
               if (str2num(csa->field, &ub) != 0)
                  error(csa, "row upper bound missing or invalid");
            }
            else
               ub = 0.0;
            if (rf[i] & 0x01)
               error(csa, "duplicate row descriptor");
            glp_set_row_bnds(P, i, type, lb, ub), rf[i] |= 0x01;
         }
         else if (strcmp(csa->field, "j") == 0)
         {  /* column descriptor */
            read_field(csa);
            if (str2int(csa->field, &j) != 0)
               error(csa, "column number missing or invalid");
            if (!(1 <= j && j <= n))
               error(csa, "column number out of range");
            if (!mip)
               kind = GLP_CV;
            else
            {  read_field(csa);
               if (strcmp(csa->field, "c") == 0)
                  kind = GLP_CV;
               else if (strcmp(csa->field, "i") == 0)
                  kind = GLP_IV;
               else if (strcmp(csa->field, "b") == 0)
               {  kind = GLP_IV;
                  type = GLP_DB, lb = 0.0, ub = 1.0;
                  goto skip;
               }
               else
                  error(csa, "column kind missing or invalid");
            }
            read_field(csa);
            if (strcmp(csa->field, "f") == 0)
               type = GLP_FR;
            else if (strcmp(csa->field, "l") == 0)
               type = GLP_LO;
            else if (strcmp(csa->field, "u") == 0)
               type = GLP_UP;
            else if (strcmp(csa->field, "d") == 0)
               type = GLP_DB;
            else if (strcmp(csa->field, "s") == 0)
               type = GLP_FX;
            else
               error(csa, "column type missing or invalid");
            if (type == GLP_LO || type == GLP_DB || type == GLP_FX)
            {  read_field(csa);
               if (str2num(csa->field, &lb) != 0)
                  error(csa, "column lower bound/fixed value missing or"
                     " invalid");
            }
            else
               lb = 0.0;
            if (type == GLP_UP || type == GLP_DB)
            {  read_field(csa);
               if (str2num(csa->field, &ub) != 0)
                  error(csa, "column upper bound missing or invalid");
            }
            else
               ub = 0.0;
skip:       if (cf[j] & 0x01)
               error(csa, "duplicate column descriptor");
            glp_set_col_kind(P, j, kind);
            glp_set_col_bnds(P, j, type, lb, ub), cf[j] |= 0x01;
         }
         else if (strcmp(csa->field, "a") == 0)
         {  /* coefficient descriptor */
            read_field(csa);
            if (str2int(csa->field, &i) != 0)
               error(csa, "row number missing or invalid");
            if (!(0 <= i && i <= m))
               error(csa, "row number out of range");
            read_field(csa);
            if (str2int(csa->field, &j) != 0)
               error(csa, "column number missing or invalid");
            if (!((i == 0 ? 0 : 1) <= j && j <= n))
               error(csa, "column number out of range");
            read_field(csa);
            if (i == 0)
            {  if (str2num(csa->field, &temp) != 0)
                  error(csa, "objective %s missing or invalid",
                     j == 0 ? "constant term" : "coefficient");
               if (cf[j] & 0x10)
                  error(csa, "duplicate objective %s",
                     j == 0 ? "constant term" : "coefficient");
               glp_set_obj_coef(P, j, temp), cf[j] |= 0x10;
            }
            else
            {  if (str2num(csa->field, &temp) != 0)
                  error(csa, "constraint coefficient missing or invalid"
                     );
               if (ne == nnz)
                  error(csa, "too many constraint coefficient descripto"
                     "rs");
               ln[++ne] = csa->count;
               ia[ne] = i, ja[ne] = j, ar[ne] = temp;
            }
         }
         else if (strcmp(csa->field, "n") == 0)
         {  /* symbolic name descriptor */
            read_field(csa);
            if (strcmp(csa->field, "p") == 0)
            {  /* problem name */
               read_field(csa);
               if (P->name != NULL)
                  error(csa, "duplicate problem name");
               glp_set_prob_name(P, csa->field);
            }
            else if (strcmp(csa->field, "z") == 0)
            {  /* objective name */
               read_field(csa);
               if (P->obj != NULL)
                  error(csa, "duplicate objective name");
               glp_set_obj_name(P, csa->field);
            }
            else if (strcmp(csa->field, "i") == 0)
            {  /* row name */
               read_field(csa);
               if (str2int(csa->field, &i) != 0)
                  error(csa, "row number missing or invalid");
               if (!(1 <= i && i <= m))
                  error(csa, "row number out of range");
               read_field(csa);
               if (P->row[i]->name != NULL)
                  error(csa, "duplicate row name");
               glp_set_row_name(P, i, csa->field);
            }
            else if (strcmp(csa->field, "j") == 0)
            {  /* column name */
               read_field(csa);
               if (str2int(csa->field, &j) != 0)
                  error(csa, "column number missing or invalid");
               if (!(1 <= j && j <= n))
                  error(csa, "column number out of range");
               read_field(csa);
               if (P->col[j]->name != NULL)
                  error(csa, "duplicate column name");
               glp_set_col_name(P, j, csa->field);
            }
            else
               error(csa, "object designator missing or invalid");
         }
         else if (strcmp(csa->field, "e") == 0)
            break;
         else
            error(csa, "line designator missing or invalid");
         end_of_line(csa);
      }
      if (ne < nnz)
         error(csa, "too few constraint coefficient descriptors");
      xassert(ne == nnz);
      k = glp_check_dup(m, n, ne, ia, ja);
      xassert(0 <= k && k <= nnz);
      if (k > 0)
      {  csa->count = ln[k];
         error(csa, "duplicate constraint coefficient");
      }
      glp_load_matrix(P, ne, ia, ja, ar);
      /* print some statistics */
      if (P->name != NULL)
         xprintf("Problem: %s\n", P->name);
      if (P->obj != NULL)
         xprintf("Objective: %s\n", P->obj);
      xprintf("%d row%s, %d column%s, %d non-zero%s\n",
         m, m == 1 ? "" : "s", n, n == 1 ? "" : "s", nnz, nnz == 1 ?
         "" : "s");
      if (glp_get_num_int(P) > 0)
      {  int ni = glp_get_num_int(P);
         int nb = glp_get_num_bin(P);
         if (ni == 1)
         {  if (nb == 0)
               xprintf("One variable is integer\n");
            else
               xprintf("One variable is binary\n");
         }
         else
         {  xprintf("%d integer variables, ", ni);
            if (nb == 0)
               xprintf("none");
            else if (nb == 1)
               xprintf("one");
            else if (nb == ni)
               xprintf("all");
            else
               xprintf("%d", nb);
            xprintf(" of which %s binary\n", nb == 1 ? "is" : "are");
         }
      }
      xprintf("%d lines were read\n", csa->count);
      /* problem data has been successfully read */
      glp_sort_matrix(P);
      ret = 0;
done: if (csa->fp != NULL) glp_close(csa->fp);
      if (rf != NULL) xfree(rf);
      if (cf != NULL) xfree(cf);
      if (ln != NULL) xfree(ln);
      if (ia != NULL) xfree(ia);
      if (ja != NULL) xfree(ja);
      if (ar != NULL) xfree(ar);
      if (ret) glp_erase_prob(P);
      return ret;
}

/* eof */
