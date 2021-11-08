/* wrmcf.c (write min-cost flow problem data in DIMACS format) */

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
*  glp_write_mincost - write min-cost flow probl. data in DIMACS format
*
*  SYNOPSIS
*
*  int glp_write_mincost(glp_graph *G, int v_rhs, int a_low, int a_cap,
*     int a_cost, const char *fname);
*
*  DESCRIPTION
*
*  The routine glp_write_mincost writes minimum cost flow problem data
*  in DIMACS format to a text file.
*
*  RETURNS
*
*  If the operation was successful, the routine returns zero. Otherwise
*  it prints an error message and returns non-zero. */

int glp_write_mincost(glp_graph *G, int v_rhs, int a_low, int a_cap,
      int a_cost, const char *fname)
{     glp_file *fp;
      glp_vertex *v;
      glp_arc *a;
      int i, count = 0, ret;
      double rhs, low, cap, cost;
      if (v_rhs >= 0 && v_rhs > G->v_size - (int)sizeof(double))
         xerror("glp_write_mincost: v_rhs = %d; invalid offset\n",
            v_rhs);
      if (a_low >= 0 && a_low > G->a_size - (int)sizeof(double))
         xerror("glp_write_mincost: a_low = %d; invalid offset\n",
            a_low);
      if (a_cap >= 0 && a_cap > G->a_size - (int)sizeof(double))
         xerror("glp_write_mincost: a_cap = %d; invalid offset\n",
            a_cap);
      if (a_cost >= 0 && a_cost > G->a_size - (int)sizeof(double))
         xerror("glp_write_mincost: a_cost = %d; invalid offset\n",
            a_cost);
      xprintf("Writing min-cost flow problem data to '%s'...\n",
         fname);
      fp = glp_open(fname, "w");
      if (fp == NULL)
      {  xprintf("Unable to create '%s' - %s\n", fname, get_err_msg());
         ret = 1;
         goto done;
      }
      xfprintf(fp, "c %s\n",
         G->name == NULL ? "unknown" : G->name), count++;
      xfprintf(fp, "p min %d %d\n", G->nv, G->na), count++;
      if (v_rhs >= 0)
      {  for (i = 1; i <= G->nv; i++)
         {  v = G->v[i];
            memcpy(&rhs, (char *)v->data + v_rhs, sizeof(double));
            if (rhs != 0.0)
               xfprintf(fp, "n %d %.*g\n", i, DBL_DIG, rhs), count++;
         }
      }
      for (i = 1; i <= G->nv; i++)
      {  v = G->v[i];
         for (a = v->out; a != NULL; a = a->t_next)
         {  if (a_low >= 0)
               memcpy(&low, (char *)a->data + a_low, sizeof(double));
            else
               low = 0.0;
            if (a_cap >= 0)
               memcpy(&cap, (char *)a->data + a_cap, sizeof(double));
            else
               cap = 1.0;
            if (a_cost >= 0)
               memcpy(&cost, (char *)a->data + a_cost, sizeof(double));
            else
               cost = 0.0;
            xfprintf(fp, "a %d %d %.*g %.*g %.*g\n",
               a->tail->i, a->head->i, DBL_DIG, low, DBL_DIG, cap,
               DBL_DIG, cost), count++;
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
