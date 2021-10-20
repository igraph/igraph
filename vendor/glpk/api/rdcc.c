/* rdcc.c (read graph in DIMACS clique/coloring format) */

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
*  glp_read_ccdata - read graph in DIMACS clique/coloring format
*
*  SYNOPSIS
*
*  int glp_read_ccdata(glp_graph *G, int v_wgt, const char *fname);
*
*  DESCRIPTION
*
*  The routine glp_read_ccdata reads an (undirected) graph in DIMACS
*  clique/coloring format from a text file.
*
*  RETURNS
*
*  If the operation was successful, the routine returns zero. Otherwise
*  it prints an error message and returns non-zero. */

int glp_read_ccdata(glp_graph *G, int v_wgt, const char *fname)
{     DMX _csa, *csa = &_csa;
      glp_vertex *v;
      int i, j, k, nv, ne, ret = 0;
      double w;
      char *flag = NULL;
      if (v_wgt >= 0 && v_wgt > G->v_size - (int)sizeof(double))
         xerror("glp_read_ccdata: v_wgt = %d; invalid offset\n",
            v_wgt);
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
      xprintf("Reading graph from '%s'...\n", fname);
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
      if (strcmp(csa->field, "edge") != 0)
         error(csa, "wrong problem designator; 'edge' expected");
      read_field(csa);
      if (!(str2int(csa->field, &nv) == 0 && nv >= 0))
         error(csa, "number of vertices missing or invalid");
      read_field(csa);
      if (!(str2int(csa->field, &ne) == 0 && ne >= 0))
         error(csa, "number of edges missing or invalid");
      xprintf("Graph has %d vert%s and %d edge%s\n",
         nv, nv == 1 ? "ex" : "ices", ne, ne == 1 ? "" : "s");
      if (nv > 0) glp_add_vertices(G, nv);
      end_of_line(csa);
      /* read node descriptor lines */
      flag = xcalloc(1+nv, sizeof(char));
      memset(&flag[1], 0, nv * sizeof(char));
      if (v_wgt >= 0)
      {  w = 1.0;
         for (i = 1; i <= nv; i++)
         {  v = G->v[i];
            memcpy((char *)v->data + v_wgt, &w, sizeof(double));
         }
      }
      for (;;)
      {  read_designator(csa);
         if (strcmp(csa->field, "n") != 0) break;
         read_field(csa);
         if (str2int(csa->field, &i) != 0)
            error(csa, "vertex number missing or invalid");
         if (!(1 <= i && i <= nv))
            error(csa, "vertex number %d out of range", i);
         if (flag[i])
            error(csa, "duplicate descriptor of vertex %d", i);
         read_field(csa);
         if (str2num(csa->field, &w) != 0)
            error(csa, "vertex weight missing or invalid");
         check_int(csa, w);
         if (v_wgt >= 0)
         {  v = G->v[i];
            memcpy((char *)v->data + v_wgt, &w, sizeof(double));
         }
         flag[i] = 1;
         end_of_line(csa);
      }
      xfree(flag), flag = NULL;
      /* read edge descriptor lines */
      for (k = 1; k <= ne; k++)
      {  if (k > 1) read_designator(csa);
         if (strcmp(csa->field, "e") != 0)
            error(csa, "wrong line designator; 'e' expected");
         read_field(csa);
         if (str2int(csa->field, &i) != 0)
            error(csa, "first vertex number missing or invalid");
         if (!(1 <= i && i <= nv))
            error(csa, "first vertex number %d out of range", i);
         read_field(csa);
         if (str2int(csa->field, &j) != 0)
            error(csa, "second vertex number missing or invalid");
         if (!(1 <= j && j <= nv))
            error(csa, "second vertex number %d out of range", j);
         glp_add_arc(G, i, j);
         end_of_line(csa);
      }
      xprintf("%d lines were read\n", csa->count);
done: if (ret) glp_erase_graph(G, G->v_size, G->a_size);
      if (csa->fp != NULL) glp_close(csa->fp);
      if (flag != NULL) xfree(flag);
      return ret;
}

/**********************************************************************/

int glp_read_graph(glp_graph *G, const char *fname)
{     return
         glp_read_ccdata(G, -1, fname);
}

/* eof */
