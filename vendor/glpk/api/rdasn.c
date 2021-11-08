/* rdasn.c (read assignment problem data in DIMACS format) */

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

#include "dimacs.h"
#include "glpk.h"
#include "misc.h"

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
*  glp_read_asnprob - read assignment problem data in DIMACS format
*
*  SYNOPSIS
*
*  int glp_read_asnprob(glp_graph *G, int v_set, int a_cost,
*     const char *fname);
*
*  DESCRIPTION
*
*  The routine glp_read_asnprob reads assignment problem data in DIMACS
*  format from a text file.
*
*  RETURNS
*
*  If the operation was successful, the routine returns zero. Otherwise
*  it prints an error message and returns non-zero. */

int glp_read_asnprob(glp_graph *G, int v_set, int a_cost, const char
      *fname)
{     DMX _csa, *csa = &_csa;
      glp_vertex *v;
      glp_arc *a;
      int nv, na, n1, i, j, k, ret = 0;
      double cost;
      char *flag = NULL;
      if (v_set >= 0 && v_set > G->v_size - (int)sizeof(int))
         xerror("glp_read_asnprob: v_set = %d; invalid offset\n",
            v_set);
      if (a_cost >= 0 && a_cost > G->a_size - (int)sizeof(double))
         xerror("glp_read_asnprob: a_cost = %d; invalid offset\n",
            a_cost);
      glp_erase_graph(G, G->v_size, G->a_size);
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
      xprintf("Reading assignment problem data from '%s'...\n", fname);
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
      if (strcmp(csa->field, "asn") != 0)
         error(csa, "wrong problem designator; 'asn' expected");
      read_field(csa);
      if (!(str2int(csa->field, &nv) == 0 && nv >= 0))
         error(csa, "number of nodes missing or invalid");
      read_field(csa);
      if (!(str2int(csa->field, &na) == 0 && na >= 0))
         error(csa, "number of arcs missing or invalid");
      if (nv > 0) glp_add_vertices(G, nv);
      end_of_line(csa);
      /* read node descriptor lines */
      flag = xcalloc(1+nv, sizeof(char));
      memset(&flag[1], 0, nv * sizeof(char));
      n1 = 0;
      for (;;)
      {  read_designator(csa);
         if (strcmp(csa->field, "n") != 0) break;
         read_field(csa);
         if (str2int(csa->field, &i) != 0)
            error(csa, "node number missing or invalid");
         if (!(1 <= i && i <= nv))
            error(csa, "node number %d out of range", i);
         if (flag[i])
            error(csa, "duplicate descriptor of node %d", i);
         flag[i] = 1, n1++;
         end_of_line(csa);
      }
      xprintf(
         "Assignment problem has %d + %d = %d node%s and %d arc%s\n",
         n1, nv - n1, nv, nv == 1 ? "" : "s", na, na == 1 ? "" : "s");
      if (v_set >= 0)
      {  for (i = 1; i <= nv; i++)
         {  v = G->v[i];
            k = (flag[i] ? 0 : 1);
            memcpy((char *)v->data + v_set, &k, sizeof(int));
         }
      }
      /* read arc descriptor lines */
      for (k = 1; k <= na; k++)
      {  if (k > 1) read_designator(csa);
         if (strcmp(csa->field, "a") != 0)
            error(csa, "wrong line designator; 'a' expected");
         read_field(csa);
         if (str2int(csa->field, &i) != 0)
            error(csa, "starting node number missing or invalid");
         if (!(1 <= i && i <= nv))
            error(csa, "starting node number %d out of range", i);
         if (!flag[i])
            error(csa, "node %d cannot be a starting node", i);
         read_field(csa);
         if (str2int(csa->field, &j) != 0)
            error(csa, "ending node number missing or invalid");
         if (!(1 <= j && j <= nv))
            error(csa, "ending node number %d out of range", j);
         if (flag[j])
            error(csa, "node %d cannot be an ending node", j);
         read_field(csa);
         if (str2num(csa->field, &cost) != 0)
            error(csa, "arc cost missing or invalid");
         check_int(csa, cost);
         a = glp_add_arc(G, i, j);
         if (a_cost >= 0)
            memcpy((char *)a->data + a_cost, &cost, sizeof(double));
         end_of_line(csa);
      }
      xprintf("%d lines were read\n", csa->count);
done: if (ret) glp_erase_graph(G, G->v_size, G->a_size);
      if (csa->fp != NULL) glp_close(csa->fp);
      if (flag != NULL) xfree(flag);
      return ret;
}

/* eof */
