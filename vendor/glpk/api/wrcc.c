/* wrcc.c (write graph in DIMACS clique/coloring format) */

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
#include "glpk.h"

#define xfprintf        glp_format

/***********************************************************************
*  NAME
*
*  glp_write_ccdata - write graph in DIMACS clique/coloring format
*
*  SYNOPSIS
*
*  int glp_write_ccdata(glp_graph *G, int v_wgt, const char *fname);
*
*  DESCRIPTION
*
*  The routine glp_write_ccdata writes the specified graph in DIMACS
*  clique/coloring format to a text file.
*
*  RETURNS
*
*  If the operation was successful, the routine returns zero. Otherwise
*  it prints an error message and returns non-zero. */

int glp_write_ccdata(glp_graph *G, int v_wgt, const char *fname)
{     glp_file *fp;
      glp_vertex *v;
      glp_arc *e;
      int i, count = 0, ret;
      double w;
      if (v_wgt >= 0 && v_wgt > G->v_size - (int)sizeof(double))
         xerror("glp_write_ccdata: v_wgt = %d; invalid offset\n",
            v_wgt);
      xprintf("Writing graph to '%s'\n", fname);
      fp = glp_open(fname, "w");
      if (fp == NULL)
      {  xprintf("Unable to create '%s' - %s\n", fname, get_err_msg());
         ret = 1;
         goto done;
      }
      xfprintf(fp, "c %s\n",
         G->name == NULL ? "unknown" : G->name), count++;
      xfprintf(fp, "p edge %d %d\n", G->nv, G->na), count++;
      if (v_wgt >= 0)
      {  for (i = 1; i <= G->nv; i++)
         {  v = G->v[i];
            memcpy(&w, (char *)v->data + v_wgt, sizeof(double));
            if (w != 1.0)
               xfprintf(fp, "n %d %.*g\n", i, DBL_DIG, w), count++;
         }
      }
      for (i = 1; i <= G->nv; i++)
      {  v = G->v[i];
         for (e = v->out; e != NULL; e = e->t_next)
            xfprintf(fp, "e %d %d\n", e->tail->i, e->head->i), count++;
      }
      xfprintf(fp, "c eof\n"), count++;
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

/**********************************************************************/

int glp_write_graph(glp_graph *G, const char *fname)
{     return
         glp_write_ccdata(G, -1, fname);
}

/* eof */
