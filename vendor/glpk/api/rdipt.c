/* rdipt.c (read interior-point solution in GLPK format) */

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
#include "env.h"
#include "misc.h"
#include "prob.h"

/***********************************************************************
*  NAME
*
*  glp_read_ipt - read interior-point solution in GLPK format
*
*  SYNOPSIS
*
*  int glp_read_ipt(glp_prob *P, const char *fname);
*
*  DESCRIPTION
*
*  The routine glp_read_ipt reads interior-point solution from a text
*  file in GLPK format.
*
*  RETURNS
*
*  If the operation was successful, the routine returns zero. Otherwise
*  it prints an error message and returns non-zero. */

int glp_read_ipt(glp_prob *P, const char *fname)
{     DMX dmx_, *dmx = &dmx_;
      int i, j, k, m, n, sst, ret = 1;
      char *stat = NULL;
      double obj, *prim = NULL, *dual = NULL;
#if 0 /* 04/IV-2016 */
      if (P == NULL || P->magic != GLP_PROB_MAGIC)
         xerror("glp_read_ipt: P = %p; invalid problem object\n", P);
#endif
      if (fname == NULL)
         xerror("glp_read_ipt: fname = %d; invalid parameter\n", fname);
      if (setjmp(dmx->jump))
         goto done;
      dmx->fname = fname;
      dmx->fp = NULL;
      dmx->count = 0;
      dmx->c = '\n';
      dmx->field[0] = '\0';
      dmx->empty = dmx->nonint = 0;
      xprintf("Reading interior-point solution from '%s'...\n", fname);
      dmx->fp = glp_open(fname, "r");
      if (dmx->fp == NULL)
      {  xprintf("Unable to open '%s' - %s\n", fname, get_err_msg());
         goto done;
      }
      /* read solution line */
      dmx_read_designator(dmx);
      if (strcmp(dmx->field, "s") != 0)
         dmx_error(dmx, "solution line missing or invalid");
      dmx_read_field(dmx);
      if (strcmp(dmx->field, "ipt") != 0)
         dmx_error(dmx, "wrong solution designator; 'ipt' expected");
      dmx_read_field(dmx);
      if (!(str2int(dmx->field, &m) == 0 && m >= 0))
         dmx_error(dmx, "number of rows missing or invalid");
      if (m != P->m)
         dmx_error(dmx, "number of rows mismatch");
      dmx_read_field(dmx);
      if (!(str2int(dmx->field, &n) == 0 && n >= 0))
         dmx_error(dmx, "number of columns missing or invalid");
      if (n != P->n)
         dmx_error(dmx, "number of columns mismatch");
      dmx_read_field(dmx);
      if (strcmp(dmx->field, "o") == 0)
         sst = GLP_OPT;
      else if (strcmp(dmx->field, "i") == 0)
         sst = GLP_INFEAS;
      else if (strcmp(dmx->field, "n") == 0)
         sst = GLP_NOFEAS;
      else if (strcmp(dmx->field, "u") == 0)
         sst = GLP_UNDEF;
      else
         dmx_error(dmx, "solution status missing or invalid");
      dmx_read_field(dmx);
      if (str2num(dmx->field, &obj) != 0)
         dmx_error(dmx, "objective value missing or invalid");
      dmx_end_of_line(dmx);
      /* allocate working arrays */
      stat = xalloc(1+m+n, sizeof(stat[0]));
      for (k = 1; k <= m+n; k++)
         stat[k] = '?';
      prim = xalloc(1+m+n, sizeof(prim[0]));
      dual = xalloc(1+m+n, sizeof(dual[0]));
      /* read solution descriptor lines */
      for (;;)
      {  dmx_read_designator(dmx);
         if (strcmp(dmx->field, "i") == 0)
         {  /* row solution descriptor */
            dmx_read_field(dmx);
            if (str2int(dmx->field, &i) != 0)
               dmx_error(dmx, "row number missing or invalid");
            if (!(1 <= i && i <= m))
               dmx_error(dmx, "row number out of range");
            if (stat[i] != '?')
               dmx_error(dmx, "duplicate row solution descriptor");
            stat[i] = GLP_BS;
            dmx_read_field(dmx);
            if (str2num(dmx->field, &prim[i]) != 0)
               dmx_error(dmx, "row primal value missing or invalid");
            dmx_read_field(dmx);
            if (str2num(dmx->field, &dual[i]) != 0)
               dmx_error(dmx, "row dual value missing or invalid");
            dmx_end_of_line(dmx);
         }
         else if (strcmp(dmx->field, "j") == 0)
         {  /* column solution descriptor */
            dmx_read_field(dmx);
            if (str2int(dmx->field, &j) != 0)
               dmx_error(dmx, "column number missing or invalid");
            if (!(1 <= j && j <= n))
               dmx_error(dmx, "column number out of range");
            if (stat[m+j] != '?')
               dmx_error(dmx, "duplicate column solution descriptor");
            stat[m+j] = GLP_BS;
            dmx_read_field(dmx);
            if (str2num(dmx->field, &prim[m+j]) != 0)
               dmx_error(dmx, "column primal value missing or invalid");
            dmx_read_field(dmx);
            if (str2num(dmx->field, &dual[m+j]) != 0)
               dmx_error(dmx, "column dual value missing or invalid");
            dmx_end_of_line(dmx);
         }
         else if (strcmp(dmx->field, "e") == 0)
            break;
         else
            dmx_error(dmx, "line designator missing or invalid");
         dmx_end_of_line(dmx);
      }
      /* store solution components into problem object */
      for (k = 1; k <= m+n; k++)
      {  if (stat[k] == '?')
            dmx_error(dmx, "incomplete interior-point solution");
      }
      P->ipt_stat = sst;
      P->ipt_obj = obj;
      for (i = 1; i <= m; i++)
      {  P->row[i]->pval = prim[i];
         P->row[i]->dval = dual[i];
      }
      for (j = 1; j <= n; j++)
      {  P->col[j]->pval = prim[m+j];
         P->col[j]->dval = dual[m+j];
      }
      /* interior-point solution has been successfully read */
      xprintf("%d lines were read\n", dmx->count);
      ret = 0;
done: if (dmx->fp != NULL)
         glp_close(dmx->fp);
      if (stat != NULL)
         xfree(stat);
      if (prim != NULL)
         xfree(prim);
      if (dual != NULL)
         xfree(dual);
      return ret;
}

/* eof */
