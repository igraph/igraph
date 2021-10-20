/* prob3.c (problem row/column searching routines) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2000-2013 Free Software Foundation, Inc.
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
#include "prob.h"

/***********************************************************************
*  NAME
*
*  glp_create_index - create the name index
*
*  SYNOPSIS
*
*  void glp_create_index(glp_prob *lp);
*
*  DESCRIPTION
*
*  The routine glp_create_index creates the name index for the
*  specified problem object. The name index is an auxiliary data
*  structure, which is intended to quickly (i.e. for logarithmic time)
*  find rows and columns by their names.
*
*  This routine can be called at any time. If the name index already
*  exists, the routine does nothing. */

void glp_create_index(glp_prob *lp)
{     GLPROW *row;
      GLPCOL *col;
      int i, j;
      /* create row name index */
      if (lp->r_tree == NULL)
      {  lp->r_tree = avl_create_tree(avl_strcmp, NULL);
         for (i = 1; i <= lp->m; i++)
         {  row = lp->row[i];
            xassert(row->node == NULL);
            if (row->name != NULL)
            {  row->node = avl_insert_node(lp->r_tree, row->name);
               avl_set_node_link(row->node, row);
            }
         }
      }
      /* create column name index */
      if (lp->c_tree == NULL)
      {  lp->c_tree = avl_create_tree(avl_strcmp, NULL);
         for (j = 1; j <= lp->n; j++)
         {  col = lp->col[j];
            xassert(col->node == NULL);
            if (col->name != NULL)
            {  col->node = avl_insert_node(lp->c_tree, col->name);
               avl_set_node_link(col->node, col);
            }
         }
      }
      return;
}

/***********************************************************************
*  NAME
*
*  glp_find_row - find row by its name
*
*  SYNOPSIS
*
*  int glp_find_row(glp_prob *lp, const char *name);
*
*  RETURNS
*
*  The routine glp_find_row returns the ordinal number of a row,
*  which is assigned (by the routine glp_set_row_name) the specified
*  symbolic name. If no such row exists, the routine returns 0. */

int glp_find_row(glp_prob *lp, const char *name)
{     AVLNODE *node;
      int i = 0;
      if (lp->r_tree == NULL)
         xerror("glp_find_row: row name index does not exist\n");
      if (!(name == NULL || name[0] == '\0' || strlen(name) > 255))
      {  node = avl_find_node(lp->r_tree, name);
         if (node != NULL)
            i = ((GLPROW *)avl_get_node_link(node))->i;
      }
      return i;
}

/***********************************************************************
*  NAME
*
*  glp_find_col - find column by its name
*
*  SYNOPSIS
*
*  int glp_find_col(glp_prob *lp, const char *name);
*
*  RETURNS
*
*  The routine glp_find_col returns the ordinal number of a column,
*  which is assigned (by the routine glp_set_col_name) the specified
*  symbolic name. If no such column exists, the routine returns 0. */

int glp_find_col(glp_prob *lp, const char *name)
{     AVLNODE *node;
      int j = 0;
      if (lp->c_tree == NULL)
         xerror("glp_find_col: column name index does not exist\n");
      if (!(name == NULL || name[0] == '\0' || strlen(name) > 255))
      {  node = avl_find_node(lp->c_tree, name);
         if (node != NULL)
            j = ((GLPCOL *)avl_get_node_link(node))->j;
      }
      return j;
}

/***********************************************************************
*  NAME
*
*  glp_delete_index - delete the name index
*
*  SYNOPSIS
*
*  void glp_delete_index(glp_prob *lp);
*
*  DESCRIPTION
*
*  The routine glp_delete_index deletes the name index previously
*  created by the routine glp_create_index and frees the memory
*  allocated to this auxiliary data structure.
*
*  This routine can be called at any time. If the name index does not
*  exist, the routine does nothing. */

void glp_delete_index(glp_prob *lp)
{     int i, j;
      /* delete row name index */
      if (lp->r_tree != NULL)
      {  for (i = 1; i <= lp->m; i++) lp->row[i]->node = NULL;
         avl_delete_tree(lp->r_tree), lp->r_tree = NULL;
      }
      /* delete column name index */
      if (lp->c_tree != NULL)
      {  for (j = 1; j <= lp->n; j++) lp->col[j]->node = NULL;
         avl_delete_tree(lp->c_tree), lp->c_tree = NULL;
      }
      return;
}

/* eof */
