/* glpapi05.c (LP basis constructing routines) */

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

#ifdef __clang__
#pragma clang diagnostic ignored "-Wlogical-op-parentheses"
#endif

#include "glpapi.h"

/***********************************************************************
*  NAME
*
*  glp_set_row_stat - set (change) row status
*
*  SYNOPSIS
*
*  void glp_set_row_stat(glp_prob *lp, int i, int stat);
*
*  DESCRIPTION
*
*  The routine glp_set_row_stat sets (changes) status of the auxiliary
*  variable associated with i-th row.
*
*  The new status of the auxiliary variable should be specified by the
*  parameter stat as follows:
*
*  GLP_BS - basic variable;
*  GLP_NL - non-basic variable;
*  GLP_NU - non-basic variable on its upper bound; if the variable is
*           not double-bounded, this means the same as GLP_NL (only in
*           case of this routine);
*  GLP_NF - the same as GLP_NL (only in case of this routine);
*  GLP_NS - the same as GLP_NL (only in case of this routine). */

void glp_set_row_stat(glp_prob *lp, int i, int stat)
{     GLPROW *row;
      if (!(1 <= i && i <= lp->m))
         xerror("glp_set_row_stat: i = %d; row number out of range\n",
            i);
      if (!(stat == GLP_BS || stat == GLP_NL || stat == GLP_NU ||
            stat == GLP_NF || stat == GLP_NS))
         xerror("glp_set_row_stat: i = %d; stat = %d; invalid status\n",
            i, stat);
      row = lp->row[i];
      if (stat != GLP_BS)
      {  switch (row->type)
         {  case GLP_FR: stat = GLP_NF; break;
            case GLP_LO: stat = GLP_NL; break;
            case GLP_UP: stat = GLP_NU; break;
            case GLP_DB: if (stat != GLP_NU) stat = GLP_NL; break;
            case GLP_FX: stat = GLP_NS; break;
            default: xassert(row != row);
         }
      }
      if (row->stat == GLP_BS && stat != GLP_BS ||
          row->stat != GLP_BS && stat == GLP_BS)
      {  /* invalidate the basis factorization */
         lp->valid = 0;
      }
      row->stat = stat;
      return;
}

/***********************************************************************
*  NAME
*
*  glp_set_col_stat - set (change) column status
*
*  SYNOPSIS
*
*  void glp_set_col_stat(glp_prob *lp, int j, int stat);
*
*  DESCRIPTION
*
*  The routine glp_set_col_stat sets (changes) status of the structural
*  variable associated with j-th column.
*
*  The new status of the structural variable should be specified by the
*  parameter stat as follows:
*
*  GLP_BS - basic variable;
*  GLP_NL - non-basic variable;
*  GLP_NU - non-basic variable on its upper bound; if the variable is
*           not double-bounded, this means the same as GLP_NL (only in
*           case of this routine);
*  GLP_NF - the same as GLP_NL (only in case of this routine);
*  GLP_NS - the same as GLP_NL (only in case of this routine). */

void glp_set_col_stat(glp_prob *lp, int j, int stat)
{     GLPCOL *col;
      if (!(1 <= j && j <= lp->n))
         xerror("glp_set_col_stat: j = %d; column number out of range\n"
            , j);
      if (!(stat == GLP_BS || stat == GLP_NL || stat == GLP_NU ||
            stat == GLP_NF || stat == GLP_NS))
         xerror("glp_set_col_stat: j = %d; stat = %d; invalid status\n",
            j, stat);
      col = lp->col[j];
      if (stat != GLP_BS)
      {  switch (col->type)
         {  case GLP_FR: stat = GLP_NF; break;
            case GLP_LO: stat = GLP_NL; break;
            case GLP_UP: stat = GLP_NU; break;
            case GLP_DB: if (stat != GLP_NU) stat = GLP_NL; break;
            case GLP_FX: stat = GLP_NS; break;
            default: xassert(col != col);
         }
      }
      if (col->stat == GLP_BS && stat != GLP_BS ||
          col->stat != GLP_BS && stat == GLP_BS)
      {  /* invalidate the basis factorization */
         lp->valid = 0;
      }
      col->stat = stat;
      return;
}

/***********************************************************************
*  NAME
*
*  glp_std_basis - construct standard initial LP basis
*
*  SYNOPSIS
*
*  void glp_std_basis(glp_prob *lp);
*
*  DESCRIPTION
*
*  The routine glp_std_basis builds the "standard" (trivial) initial
*  basis for the specified problem object.
*
*  In the "standard" basis all auxiliary variables are basic, and all
*  structural variables are non-basic. */

void glp_std_basis(glp_prob *lp)
{     int i, j;
      /* make all auxiliary variables basic */
      for (i = 1; i <= lp->m; i++)
         glp_set_row_stat(lp, i, GLP_BS);
      /* make all structural variables non-basic */
      for (j = 1; j <= lp->n; j++)
      {  GLPCOL *col = lp->col[j];
         if (col->type == GLP_DB && fabs(col->lb) > fabs(col->ub))
            glp_set_col_stat(lp, j, GLP_NU);
         else
            glp_set_col_stat(lp, j, GLP_NL);
      }
      return;
}

/* eof */
