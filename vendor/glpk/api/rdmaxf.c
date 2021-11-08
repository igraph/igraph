/* rdmaxf.c (read maximum flow problem data in DIMACS format) */

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
*  glp_read_maxflow - read maximum flow problem data in DIMACS format
*
*  SYNOPSIS
*
*  int glp_read_maxflow(glp_graph *G, int *s, int *t, int a_cap,
*     const char *fname);
*
*  DESCRIPTION
*
*  The routine glp_read_maxflow reads maximum flow problem data in
*  DIMACS format from a text file.
*
*  RETURNS
*
*  If the operation was successful, the routine returns zero. Otherwise
*  it prints an error message and returns non-zero. */

int glp_read_maxflow(glp_graph *G, int *_s, int *_t, int a_cap,
      const char *fname)
{     DMX _csa, *csa = &_csa;
      glp_arc *a;
      int i, j, k, s, t, nv, na, ret = 0;
      double cap;
      if (a_cap >= 0 && a_cap > G->a_size - (int)sizeof(double))
         xerror("glp_read_maxflow: a_cap = %d; invalid offset\n",
            a_cap);
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
      xprintf("Reading maximum flow problem data from '%s'...\n",
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
      if (strcmp(csa->field, "max") != 0)
         error(csa, "wrong problem designator; 'max' expected");
      read_field(csa);
      if (!(str2int(csa->field, &nv) == 0 && nv >= 2))
         error(csa, "number of nodes missing or invalid");
      read_field(csa);
      if (!(str2int(csa->field, &na) == 0 && na >= 0))
         error(csa, "number of arcs missing or invalid");
      xprintf("Flow network has %d node%s and %d arc%s\n",
         nv, nv == 1 ? "" : "s", na, na == 1 ? "" : "s");
      if (nv > 0) glp_add_vertices(G, nv);
      end_of_line(csa);
      /* read node descriptor lines */
      s = t = 0;
      for (;;)
      {  read_designator(csa);
         if (strcmp(csa->field, "n") != 0) break;
         read_field(csa);
         if (str2int(csa->field, &i) != 0)
            error(csa, "node number missing or invalid");
         if (!(1 <= i && i <= nv))
            error(csa, "node number %d out of range", i);
         read_field(csa);
         if (strcmp(csa->field, "s") == 0)
         {  if (s > 0)
               error(csa, "only one source node allowed");
            s = i;
         }
         else if (strcmp(csa->field, "t") == 0)
         {  if (t > 0)
               error(csa, "only one sink node allowed");
            t = i;
         }
         else
            error(csa, "wrong node designator; 's' or 't' expected");
         if (s > 0 && s == t)
            error(csa, "source and sink nodes must be distinct");
         end_of_line(csa);
      }
      if (s == 0)
         error(csa, "source node descriptor missing\n");
      if (t == 0)
         error(csa, "sink node descriptor missing\n");
      if (_s != NULL) *_s = s;
      if (_t != NULL) *_t = t;
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
         if (!(str2num(csa->field, &cap) == 0 && cap >= 0.0))
            error(csa, "arc capacity missing or invalid");
         check_int(csa, cap);
         a = glp_add_arc(G, i, j);
         if (a_cap >= 0)
            memcpy((char *)a->data + a_cap, &cap, sizeof(double));
         end_of_line(csa);
      }
      xprintf("%d lines were read\n", csa->count);
done: if (ret) glp_erase_graph(G, G->v_size, G->a_size);
      if (csa->fp != NULL) glp_close(csa->fp);
      return ret;
}

/* eof */
