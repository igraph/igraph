/* ckcnf.c (check for CNF-SAT problem instance) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2010-2016 Free Software Foundation, Inc.
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

int glp_check_cnfsat(glp_prob *P)
{     /* check for CNF-SAT problem instance */
      int m = P->m;
      int n = P->n;
      GLPROW *row;
      GLPCOL *col;
      GLPAIJ *aij;
      int i, j, neg;
#if 0 /* 04/IV-2016 */
      if (P == NULL || P->magic != GLP_PROB_MAGIC)
         xerror("glp_check_cnfsat: P = %p; invalid problem object\n",
            P);
#endif
      /* check columns */
      for (j = 1; j <= n; j++)
      {  col = P->col[j];
         /* the variable should be binary */
         if (!(col->kind == GLP_IV && col->type == GLP_DB &&
               col->lb == 0.0 && col->ub == 1.0))
            return 1;
      }
      /* objective function should be zero */
      if (P->c0 != 0.0)
         return 2;
      for (j = 1; j <= n; j++)
      {  col = P->col[j];
         if (col->coef != 0.0)
            return 3;
      }
      /* check rows */
      for (i = 1; i <= m; i++)
      {  row = P->row[i];
         /* the row should be of ">=" type */
         if (row->type != GLP_LO)
            return 4;
         /* check constraint coefficients */
         neg = 0;
         for (aij = row->ptr; aij != NULL; aij = aij->r_next)
         {  /* the constraint coefficient should be +1 or -1 */
            if (aij->val == +1.0)
               ;
            else if (aij->val == -1.0)
               neg++;
            else
               return 5;
         }
         /* the right-hand side should be (1 - neg), where neg is the
            number of negative constraint coefficients in the row */
         if (row->lb != (double)(1 - neg))
            return 6;
      }
      /* congratulations; this is CNF-SAT */
      return 0;
}

/* eof */
