/* glpapi04.c (problem scaling routines) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*
*  Copyright (C) 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008,
*  2009, 2010 Andrew Makhorin, Department for Applied Informatics,
*  Moscow Aviation Institute, Moscow, Russia. All rights reserved.
*  E-mail: <mao@gnu.org>.
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

#include "glpapi.h"

/***********************************************************************
*  NAME
*
*  glp_set_rii - set (change) row scale factor
*
*  SYNOPSIS
*
*  void glp_set_rii(glp_prob *lp, int i, double rii);
*
*  DESCRIPTION
*
*  The routine glp_set_rii sets (changes) the scale factor r[i,i] for
*  i-th row of the specified problem object. */

void glp_set_rii(glp_prob *lp, int i, double rii)
{     if (!(1 <= i && i <= lp->m))
         xerror("glp_set_rii: i = %d; row number out of range\n", i);
      if (rii <= 0.0)
         xerror("glp_set_rii: i = %d; rii = %g; invalid scale factor\n",
            i, rii);
      if (lp->valid && lp->row[i]->rii != rii)
      {  GLPAIJ *aij;
         for (aij = lp->row[i]->ptr; aij != NULL; aij = aij->r_next)
         {  if (aij->col->stat == GLP_BS)
            {  /* invalidate the basis factorization */
               lp->valid = 0;
               break;
            }
         }
      }
      lp->row[i]->rii = rii;
      return;
}

/***********************************************************************
*  NAME
*
*  glp_set sjj - set (change) column scale factor
*
*  SYNOPSIS
*
*  void glp_set_sjj(glp_prob *lp, int j, double sjj);
*
*  DESCRIPTION
*
*  The routine glp_set_sjj sets (changes) the scale factor s[j,j] for
*  j-th column of the specified problem object. */

void glp_set_sjj(glp_prob *lp, int j, double sjj)
{     if (!(1 <= j && j <= lp->n))
         xerror("glp_set_sjj: j = %d; column number out of range\n", j);
      if (sjj <= 0.0)
         xerror("glp_set_sjj: j = %d; sjj = %g; invalid scale factor\n",
            j, sjj);
      if (lp->valid && lp->col[j]->sjj != sjj && lp->col[j]->stat ==
         GLP_BS)
      {  /* invalidate the basis factorization */
         lp->valid = 0;
      }
      lp->col[j]->sjj = sjj;
      return;
}

/***********************************************************************
*  NAME
*
*  glp_get_rii - retrieve row scale factor
*
*  SYNOPSIS
*
*  double glp_get_rii(glp_prob *lp, int i);
*
*  RETURNS
*
*  The routine glp_get_rii returns current scale factor r[i,i] for i-th
*  row of the specified problem object. */

double glp_get_rii(glp_prob *lp, int i)
{     if (!(1 <= i && i <= lp->m))
         xerror("glp_get_rii: i = %d; row number out of range\n", i);
      return lp->row[i]->rii;
}

/***********************************************************************
*  NAME
*
*  glp_get_sjj - retrieve column scale factor
*
*  SYNOPSIS
*
*  double glp_get_sjj(glp_prob *lp, int j);
*
*  RETURNS
*
*  The routine glp_get_sjj returns current scale factor s[j,j] for j-th
*  column of the specified problem object. */

double glp_get_sjj(glp_prob *lp, int j)
{     if (!(1 <= j && j <= lp->n))
         xerror("glp_get_sjj: j = %d; column number out of range\n", j);
      return lp->col[j]->sjj;
}

/***********************************************************************
*  NAME
*
*  glp_unscale_prob - unscale problem data
*
*  SYNOPSIS
*
*  void glp_unscale_prob(glp_prob *lp);
*
*  DESCRIPTION
*
*  The routine glp_unscale_prob performs unscaling of problem data for
*  the specified problem object.
*
*  "Unscaling" means replacing the current scaling matrices R and S by
*  unity matrices that cancels the scaling effect. */

void glp_unscale_prob(glp_prob *lp)
{     int m = glp_get_num_rows(lp);
      int n = glp_get_num_cols(lp);
      int i, j;
      for (i = 1; i <= m; i++) glp_set_rii(lp, i, 1.0);
      for (j = 1; j <= n; j++) glp_set_sjj(lp, j, 1.0);
      return;
}

/* eof */
