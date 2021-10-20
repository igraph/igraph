/* wrmaxf.c (write maximum flow problem data in DIMACS format) */

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
*  glp_write_maxflow - write maximum flow problem data in DIMACS format
*
*  SYNOPSIS
*
*  int glp_write_maxflow(glp_graph *G, int s, int t, int a_cap,
*     const char *fname);
*
*  DESCRIPTION
*
*  The routine glp_write_maxflow writes maximum flow problem data in
*  DIMACS format to a text file.
*
*  RETURNS
*
*  If the operation was successful, the routine returns zero. Otherwise
*  it prints an error message and returns non-zero. */

int glp_write_maxflow(glp_graph *G, int s, int t, int a_cap,
      const char *fname)
{     glp_file *fp;
      glp_vertex *v;
      glp_arc *a;
      int i, count = 0, ret;
      double cap;
      if (!(1 <= s && s <= G->nv))
         xerror("glp_write_maxflow: s = %d; source node number out of r"
            "ange\n", s);
      if (!(1 <= t && t <= G->nv))
         xerror("glp_write_maxflow: t = %d: sink node number out of ran"
            "ge\n", t);
      if (a_cap >= 0 && a_cap > G->a_size - (int)sizeof(double))
         xerror("glp_write_mincost: a_cap = %d; invalid offset\n",
            a_cap);
      xprintf("Writing maximum flow problem data to '%s'...\n",
         fname);
      fp = glp_open(fname, "w");
      if (fp == NULL)
      {  xprintf("Unable to create '%s' - %s\n", fname, get_err_msg());
         ret = 1;
         goto done;
      }
      xfprintf(fp, "c %s\n",
         G->name == NULL ? "unknown" : G->name), count++;
      xfprintf(fp, "p max %d %d\n", G->nv, G->na), count++;
      xfprintf(fp, "n %d s\n", s), count++;
      xfprintf(fp, "n %d t\n", t), count++;
      for (i = 1; i <= G->nv; i++)
      {  v = G->v[i];
         for (a = v->out; a != NULL; a = a->t_next)
         {  if (a_cap >= 0)
               memcpy(&cap, (char *)a->data + a_cap, sizeof(double));
            else
               cap = 1.0;
            xfprintf(fp, "a %d %d %.*g\n",
               a->tail->i, a->head->i, DBL_DIG, cap), count++;
         }
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

/* eof */
