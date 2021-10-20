/* wrasn.c (write assignment problem data in DIMACS format) */

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
*  glp_write_asnprob - write assignment problem data in DIMACS format
*
*  SYNOPSIS
*
*  int glp_write_asnprob(glp_graph *G, int v_set, int a_cost,
*     const char *fname);
*
*  DESCRIPTION
*
*  The routine glp_write_asnprob writes assignment problem data in
*  DIMACS format to a text file.
*
*  RETURNS
*
*  If the operation was successful, the routine returns zero. Otherwise
*  it prints an error message and returns non-zero. */

int glp_write_asnprob(glp_graph *G, int v_set, int a_cost, const char
      *fname)
{     glp_file *fp;
      glp_vertex *v;
      glp_arc *a;
      int i, k, count = 0, ret;
      double cost;
      if (v_set >= 0 && v_set > G->v_size - (int)sizeof(int))
         xerror("glp_write_asnprob: v_set = %d; invalid offset\n",
            v_set);
      if (a_cost >= 0 && a_cost > G->a_size - (int)sizeof(double))
         xerror("glp_write_asnprob: a_cost = %d; invalid offset\n",
            a_cost);
      xprintf("Writing assignment problem data to '%s'...\n", fname);
      fp = glp_open(fname, "w");
      if (fp == NULL)
      {  xprintf("Unable to create '%s' - %s\n", fname, get_err_msg());
         ret = 1;
         goto done;
      }
      xfprintf(fp, "c %s\n",
         G->name == NULL ? "unknown" : G->name), count++;
      xfprintf(fp, "p asn %d %d\n", G->nv, G->na), count++;
      for (i = 1; i <= G->nv; i++)
      {  v = G->v[i];
         if (v_set >= 0)
            memcpy(&k, (char *)v->data + v_set, sizeof(int));
         else
            k = (v->out != NULL ? 0 : 1);
         if (k == 0)
            xfprintf(fp, "n %d\n", i), count++;
      }
      for (i = 1; i <= G->nv; i++)
      {  v = G->v[i];
         for (a = v->out; a != NULL; a = a->t_next)
         {  if (a_cost >= 0)
               memcpy(&cost, (char *)a->data + a_cost, sizeof(double));
            else
               cost = 1.0;
            xfprintf(fp, "a %d %d %.*g\n",
               a->tail->i, a->head->i, DBL_DIG, cost), count++;
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
