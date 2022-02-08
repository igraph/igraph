/* advbas.c (construct advanced initial LP basis) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2008-2016 Free Software Foundation, Inc.
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
#include "triang.h"

/***********************************************************************
*  NAME
*
*  glp_adv_basis - construct advanced initial LP basis
*
*  SYNOPSIS
*
*  void glp_adv_basis(glp_prob *P, int flags);
*
*  DESCRIPTION
*
*  The routine glp_adv_basis constructs an advanced initial LP basis
*  for the specified problem object.
*
*  The parameter flag is reserved for use in the future and should be
*  specified as zero.
*
*  NOTE
*
*  The routine glp_adv_basis should be called after the constraint
*  matrix has been scaled (if scaling is used). */

static int mat(void *info, int k, int ind[], double val[])
{     glp_prob *P = info;
      int m = P->m;
      int n = P->n;
      GLPROW **row = P->row;
      GLPCOL **col = P->col;
      GLPAIJ *aij;
      int i, j, len;
      if (k > 0)
      {  /* retrieve scaled row of constraint matrix */
         i = +k;
         xassert(1 <= i && i <= m);
         len = 0;
         if (row[i]->type == GLP_FX)
         {  for (aij = row[i]->ptr; aij != NULL; aij = aij->r_next)
            {  j = aij->col->j;
               if (col[j]->type != GLP_FX)
               {  len++;
                  ind[len] = j;
                  val[len] = aij->row->rii * aij->val * aij->col->sjj;
               }
            }
         }
      }
      else
      {  /* retrieve scaled column of constraint matrix */
         j = -k;
         xassert(1 <= j && j <= n);
         len = 0;
         if (col[j]->type != GLP_FX)
         {  for (aij = col[j]->ptr; aij != NULL; aij = aij->c_next)
            {  i = aij->row->i;
               if (row[i]->type == GLP_FX)
               {  len++;
                  ind[len] = i;
                  val[len] = aij->row->rii * aij->val * aij->col->sjj;
               }
            }
         }
      }
      return len;
}

void glp_adv_basis(glp_prob *P, int flags)
{     int i, j, k, m, n, min_mn, size, *rn, *cn;
      char *flag;
      if (flags != 0)
         xerror("glp_adv_basis: flags = %d; invalid flags\n", flags);
      m = P->m; /* number of rows */
      n = P->n; /* number of columns */
      if (m == 0 || n == 0)
      {  /* trivial case */
         glp_std_basis(P);
         goto done;
      }
      xprintf("Constructing initial basis...\n");
      /* allocate working arrays */
      min_mn = (m < n ? m : n);
      rn = talloc(1+min_mn, int);
      cn = talloc(1+min_mn, int);
      flag = talloc(1+m, char);
      /* make the basis empty */
      for (i = 1; i <= m; i++)
      {  flag[i] = 0;
         glp_set_row_stat(P, i, GLP_NS);
      }
      for (j = 1; j <= n; j++)
         glp_set_col_stat(P, j, GLP_NS);
      /* find maximal triangular part of the constraint matrix;
         to prevent including non-fixed rows and fixed columns in the
         triangular part, such rows and columns are temporarily made
         empty by the routine mat */
#if 1 /* FIXME: tolerance */
      size = triang(m, n, mat, P, 0.001, rn, cn);
#endif
      xassert(0 <= size && size <= min_mn);
      /* include in the basis non-fixed structural variables, whose
         columns constitute the triangular part */
      for (k = 1; k <= size; k++)
      {  i = rn[k];
         xassert(1 <= i && i <= m);
         flag[i] = 1;
         j = cn[k];
         xassert(1 <= j && j <= n);
         glp_set_col_stat(P, j, GLP_BS);
      }
      /* include in the basis appropriate auxiliary variables, whose
         unity columns preserve triangular form of the basis matrix */
      for (i = 1; i <= m; i++)
      {  if (flag[i] == 0)
         {  glp_set_row_stat(P, i, GLP_BS);
            if (P->row[i]->type != GLP_FX)
               size++;
         }
      }
      /* size of triangular part = (number of rows) - (number of basic
         fixed auxiliary variables) */
      xprintf("Size of triangular part is %d\n", size);
      /* deallocate working arrays */
      tfree(rn);
      tfree(cn);
      tfree(flag);
done: return;
}

/* eof */
