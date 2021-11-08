/* rdmcf.c (read min-cost flow problem data in DIMACS format) */

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
*  glp_read_mincost - read min-cost flow problem data in DIMACS format
*
*  SYNOPSIS
*
*  int glp_read_mincost(glp_graph *G, int v_rhs, int a_low, int a_cap,
*     int a_cost, const char *fname);
*
*  DESCRIPTION
*
*  The routine glp_read_mincost reads minimum cost flow problem data in
*  DIMACS format from a text file.
*
*  RETURNS
*
*  If the operation was successful, the routine returns zero. Otherwise
*  it prints an error message and returns non-zero. */

int glp_read_mincost(glp_graph *G, int v_rhs, int a_low, int a_cap,
      int a_cost, const char *fname)
{     DMX _csa, *csa = &_csa;
      glp_vertex *v;
      glp_arc *a;
      int i, j, k, nv, na, ret = 0;
      double rhs, low, cap, cost;
      char *flag = NULL;
      if (v_rhs >= 0 && v_rhs > G->v_size - (int)sizeof(double))
         xerror("glp_read_mincost: v_rhs = %d; invalid offset\n",
            v_rhs);
      if (a_low >= 0 && a_low > G->a_size - (int)sizeof(double))
         xerror("glp_read_mincost: a_low = %d; invalid offset\n",
            a_low);
      if (a_cap >= 0 && a_cap > G->a_size - (int)sizeof(double))
         xerror("glp_read_mincost: a_cap = %d; invalid offset\n",
            a_cap);
      if (a_cost >= 0 && a_cost > G->a_size - (int)sizeof(double))
         xerror("glp_read_mincost: a_cost = %d; invalid offset\n",
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
      xprintf("Reading min-cost flow problem data from '%s'...\n",
         fname);
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
      if (strcmp(csa->field, "min") != 0)
         error(csa, "wrong problem designator; 'min' expected");
      read_field(csa);
      if (!(str2int(csa->field, &nv) == 0 && nv >= 0))
         error(csa, "number of nodes missing or invalid");
      read_field(csa);
      if (!(str2int(csa->field, &na) == 0 && na >= 0))
         error(csa, "number of arcs missing or invalid");
      xprintf("Flow network has %d node%s and %d arc%s\n",
         nv, nv == 1 ? "" : "s", na, na == 1 ? "" : "s");
      if (nv > 0) glp_add_vertices(G, nv);
      end_of_line(csa);
      /* read node descriptor lines */
      flag = xcalloc(1+nv, sizeof(char));
      memset(&flag[1], 0, nv * sizeof(char));
      if (v_rhs >= 0)
      {  rhs = 0.0;
         for (i = 1; i <= nv; i++)
         {  v = G->v[i];
            memcpy((char *)v->data + v_rhs, &rhs, sizeof(double));
         }
      }
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
         read_field(csa);
         if (str2num(csa->field, &rhs) != 0)
            error(csa, "node supply/demand missing or invalid");
         check_int(csa, rhs);
         if (v_rhs >= 0)
         {  v = G->v[i];
            memcpy((char *)v->data + v_rhs, &rhs, sizeof(double));
         }
         flag[i] = 1;
         end_of_line(csa);
      }
      xfree(flag), flag = NULL;
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
         read_field(csa);
         if (str2int(csa->field, &j) != 0)
            error(csa, "ending node number missing or invalid");
         if (!(1 <= j && j <= nv))
            error(csa, "ending node number %d out of range", j);
         read_field(csa);
         if (!(str2num(csa->field, &low) == 0 && low >= 0.0))
            error(csa, "lower bound of arc flow missing or invalid");
         check_int(csa, low);
         read_field(csa);
         if (!(str2num(csa->field, &cap) == 0 && cap >= low))
            error(csa, "upper bound of arc flow missing or invalid");
         check_int(csa, cap);
         read_field(csa);
         if (str2num(csa->field, &cost) != 0)
            error(csa, "per-unit cost of arc flow missing or invalid");
         check_int(csa, cost);
         a = glp_add_arc(G, i, j);
         if (a_low >= 0)
            memcpy((char *)a->data + a_low, &low, sizeof(double));
         if (a_cap >= 0)
            memcpy((char *)a->data + a_cap, &cap, sizeof(double));
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
