/* ckasn.c (check correctness of assignment problem data) */

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

/***********************************************************************
*  NAME
*
*  glp_check_asnprob - check correctness of assignment problem data
*
*  SYNOPSIS
*
*  int glp_check_asnprob(glp_graph *G, int v_set);
*
*  RETURNS
*
*  If the specified assignment problem data are correct, the routine
*  glp_check_asnprob returns zero, otherwise, non-zero. */

int glp_check_asnprob(glp_graph *G, int v_set)
{     glp_vertex *v;
      int i, k, ret = 0;
      if (v_set >= 0 && v_set > G->v_size - (int)sizeof(int))
         xerror("glp_check_asnprob: v_set = %d; invalid offset\n",
            v_set);
      for (i = 1; i <= G->nv; i++)
      {  v = G->v[i];
         if (v_set >= 0)
         {  memcpy(&k, (char *)v->data + v_set, sizeof(int));
            if (k == 0)
            {  if (v->in != NULL)
               {  ret = 1;
                  break;
               }
            }
            else if (k == 1)
            {  if (v->out != NULL)
               {  ret = 2;
                  break;
               }
            }
            else
            {  ret = 3;
               break;
            }
         }
         else
         {  if (v->in != NULL && v->out != NULL)
            {  ret = 4;
               break;
            }
         }
      }
      return ret;
}

/* eof */
