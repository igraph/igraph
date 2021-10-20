/* cfg2.c (conflict graph) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2015-2016 Free Software Foundation, Inc.
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

#include "cfg.h"
#include "env.h"
#include "prob.h"

/***********************************************************************
*  NAME
*
*  glp_cfg_init - create and initialize conflict graph
*
*  SYNOPSIS
*
*  glp_cfg *glp_cfg_init(glp_prob *P);
*
*  DESCRIPTION
*
*  This routine creates and initializes the conflict graph for the
*  specified problem object.
*
*  RETURNS
*
*  The routine returns a pointer to the conflict graph descriptor.
*  However, if the conflict graph is empty (no conflicts have been
*  found), the routine returns NULL. */

glp_cfg *glp_cfg_init(glp_prob *P)
{     glp_cfg *G;
      int j, n1, n2;
      xprintf("Constructing conflict graph...\n");
      G = cfg_build_graph(P);
      n1 = n2 = 0;
      for (j = 1; j <= P->n; j++)
      {  if (G->pos[j])
            n1 ++;
         if (G->neg[j])
            n2++;
      }
      if (n1 == 0 && n2 == 0)
      {  xprintf("No conflicts found\n");
         cfg_delete_graph(G);
         G = NULL;
      }
      else
         xprintf("Conflict graph has %d + %d = %d vertices\n",
            n1, n2, G->nv);
      return G;
}

/***********************************************************************
*  NAME
*
*  glp_cfg_free - delete conflict graph descriptor
*
*  SYNOPSIS
*
*  void glp_cfg_free(glp_cfg *G);
*
*  DESCRIPTION
*
*  This routine deletes the conflict graph descriptor and frees all the
*  memory allocated to it. */

void glp_cfg_free(glp_cfg *G)
{     xassert(G != NULL);
      cfg_delete_graph(G);
      return;
}

/* eof */
