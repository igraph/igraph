/* rdcnf.c (read CNF-SAT problem data in DIMACS format) */

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

int glp_read_cnfsat(glp_prob *P, const char *fname)
{     /* read CNF-SAT problem data in DIMACS format */
      DMX _csa, *csa = &_csa;
      int m, n, i, j, len, neg, rhs, ret = 0, *ind = NULL;
      double *val = NULL;
      char *map = NULL;
#if 0 /* 04/IV-2016 */
      if (P == NULL || P->magic != GLP_PROB_MAGIC)
         xerror("glp_read_cnfsat: P = %p; invalid problem object\n",
            P);
#endif
      if (fname == NULL)
         xerror("glp_read_cnfsat: fname = %p; invalid parameter\n",
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
      xprintf("Reading CNF-SAT problem data from '%s'...\n", fname);
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
      if (strcmp(csa->field, "cnf") != 0)
         error(csa, "wrong problem designator; 'cnf' expected\n");
      read_field(csa);
      if (!(str2int(csa->field, &n) == 0 && n >= 0))
         error(csa, "number of variables missing or invalid\n");
      read_field(csa);
      if (!(str2int(csa->field, &m) == 0 && m >= 0))
         error(csa, "number of clauses missing or invalid\n");
      xprintf("Instance has %d variable%s and %d clause%s\n",
         n, n == 1 ? "" : "s", m, m == 1 ? "" : "s");
      end_of_line(csa);
      if (m > 0)
         glp_add_rows(P, m);
      if (n > 0)
      {  glp_add_cols(P, n);
         for (j = 1; j <= n; j++)
            glp_set_col_kind(P, j, GLP_BV);
      }
      /* allocate working arrays */
      ind = xcalloc(1+n, sizeof(int));
      val = xcalloc(1+n, sizeof(double));
      map = xcalloc(1+n, sizeof(char));
      for (j = 1; j <= n; j++) map[j] = 0;
      /* read clauses */
      for (i = 1; i <= m; i++)
      {  /* read i-th clause */
         len = 0, rhs = 1;
         for (;;)
         {  /* skip white-space characters */
            while (csa->c == ' ' || csa->c == '\n')
               read_char(csa);
            /* read term */
            read_field(csa);
            if (str2int(csa->field, &j) != 0)
               error(csa, "variable number missing or invalid\n");
            if (j > 0)
               neg = 0;
            else if (j < 0)
               neg = 1, j = -j, rhs--;
            else
               break;
            if (!(1 <= j && j <= n))
               error(csa, "variable number out of range\n");
            if (map[j])
               error(csa, "duplicate variable number\n");
            len++, ind[len] = j, val[len] = (neg ? -1.0 : +1.0);
            map[j] = 1;
         }
         glp_set_row_bnds(P, i, GLP_LO, (double)rhs, 0.0);
         glp_set_mat_row(P, i, len, ind, val);
         while (len > 0) map[ind[len--]] = 0;
      }
      xprintf("%d lines were read\n", csa->count);
      /* problem data has been successfully read */
      glp_sort_matrix(P);
done: if (csa->fp != NULL) glp_close(csa->fp);
      if (ind != NULL) xfree(ind);
      if (val != NULL) xfree(val);
      if (map != NULL) xfree(map);
      if (ret) glp_erase_prob(P);
      return ret;
}

/* eof */
