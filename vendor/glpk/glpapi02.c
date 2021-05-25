/* glpapi02.c (problem retrieving routines) */

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
#pragma clang diagnostic ignored "-Wsometimes-uninitialized"
#endif

#include "glpapi.h"

/***********************************************************************
*  NAME
*
*  glp_get_prob_name - retrieve problem name
*
*  SYNOPSIS
*
*  const char *glp_get_prob_name(glp_prob *lp);
*
*  RETURNS
*
*  The routine glp_get_prob_name returns a pointer to an internal
*  buffer, which contains symbolic name of the problem. However, if the
*  problem has no assigned name, the routine returns NULL. */

const char *glp_get_prob_name(glp_prob *lp)
{     char *name;
      name = lp->name;
      return name;
}

/***********************************************************************
*  NAME
*
*  glp_get_obj_name - retrieve objective function name
*
*  SYNOPSIS
*
*  const char *glp_get_obj_name(glp_prob *lp);
*
*  RETURNS
*
*  The routine glp_get_obj_name returns a pointer to an internal
*  buffer, which contains a symbolic name of the objective function.
*  However, if the objective function has no assigned name, the routine
*  returns NULL. */

const char *glp_get_obj_name(glp_prob *lp)
{     char *name;
      name = lp->obj;
      return name;
}

/***********************************************************************
*  NAME
*
*  glp_get_obj_dir - retrieve optimization direction flag
*
*  SYNOPSIS
*
*  int glp_get_obj_dir(glp_prob *lp);
*
*  RETURNS
*
*  The routine glp_get_obj_dir returns the optimization direction flag
*  (i.e. "sense" of the objective function):
*
*  GLP_MIN - minimization;
*  GLP_MAX - maximization. */

int glp_get_obj_dir(glp_prob *lp)
{     int dir = lp->dir;
      return dir;
}

/***********************************************************************
*  NAME
*
*  glp_get_num_rows - retrieve number of rows
*
*  SYNOPSIS
*
*  int glp_get_num_rows(glp_prob *lp);
*
*  RETURNS
*
*  The routine glp_get_num_rows returns the current number of rows in
*  the specified problem object. */

int glp_get_num_rows(glp_prob *lp)
{     int m = lp->m;
      return m;
}

/***********************************************************************
*  NAME
*
*  glp_get_num_cols - retrieve number of columns
*
*  SYNOPSIS
*
*  int glp_get_num_cols(glp_prob *lp);
*
*  RETURNS
*
*  The routine glp_get_num_cols returns the current number of columns
*  in the specified problem object. */

int glp_get_num_cols(glp_prob *lp)
{     int n = lp->n;
      return n;
}

/***********************************************************************
*  NAME
*
*  glp_get_row_name - retrieve row name
*
*  SYNOPSIS
*
*  const char *glp_get_row_name(glp_prob *lp, int i);
*
*  RETURNS
*
*  The routine glp_get_row_name returns a pointer to an internal
*  buffer, which contains symbolic name of i-th row. However, if i-th
*  row has no assigned name, the routine returns NULL. */

const char *glp_get_row_name(glp_prob *lp, int i)
{     char *name;
      if (!(1 <= i && i <= lp->m))
         xerror("glp_get_row_name: i = %d; row number out of range\n",
            i);
      name = lp->row[i]->name;
      return name;
}

/***********************************************************************
*  NAME
*
*  glp_get_col_name - retrieve column name
*
*  SYNOPSIS
*
*  const char *glp_get_col_name(glp_prob *lp, int j);
*
*  RETURNS
*
*  The routine glp_get_col_name returns a pointer to an internal
*  buffer, which contains symbolic name of j-th column. However, if j-th
*  column has no assigned name, the routine returns NULL. */

const char *glp_get_col_name(glp_prob *lp, int j)
{     char *name;
      if (!(1 <= j && j <= lp->n))
         xerror("glp_get_col_name: j = %d; column number out of range\n"
            , j);
      name = lp->col[j]->name;
      return name;
}

/***********************************************************************
*  NAME
*
*  glp_get_row_type - retrieve row type
*
*  SYNOPSIS
*
*  int glp_get_row_type(glp_prob *lp, int i);
*
*  RETURNS
*
*  The routine glp_get_row_type returns the type of i-th row, i.e. the
*  type of corresponding auxiliary variable, as follows:
*
*  GLP_FR - free (unbounded) variable;
*  GLP_LO - variable with lower bound;
*  GLP_UP - variable with upper bound;
*  GLP_DB - double-bounded variable;
*  GLP_FX - fixed variable. */

int glp_get_row_type(glp_prob *lp, int i)
{     if (!(1 <= i && i <= lp->m))
         xerror("glp_get_row_type: i = %d; row number out of range\n",
            i);
      return lp->row[i]->type;
}

/***********************************************************************
*  NAME
*
*  glp_get_row_lb - retrieve row lower bound
*
*  SYNOPSIS
*
*  double glp_get_row_lb(glp_prob *lp, int i);
*
*  RETURNS
*
*  The routine glp_get_row_lb returns the lower bound of i-th row, i.e.
*  the lower bound of corresponding auxiliary variable. However, if the
*  row has no lower bound, the routine returns -DBL_MAX. */

double glp_get_row_lb(glp_prob *lp, int i)
{     double lb;
      if (!(1 <= i && i <= lp->m))
         xerror("glp_get_row_lb: i = %d; row number out of range\n", i);
      switch (lp->row[i]->type)
      {  case GLP_FR:
         case GLP_UP:
            lb = -DBL_MAX; break;
         case GLP_LO:
         case GLP_DB:
         case GLP_FX:
            lb = lp->row[i]->lb; break;
         default:
            xassert(lp != lp);
      }
      return lb;
}

/***********************************************************************
*  NAME
*
*  glp_get_row_ub - retrieve row upper bound
*
*  SYNOPSIS
*
*  double glp_get_row_ub(glp_prob *lp, int i);
*
*  RETURNS
*
*  The routine glp_get_row_ub returns the upper bound of i-th row, i.e.
*  the upper bound of corresponding auxiliary variable. However, if the
*  row has no upper bound, the routine returns +DBL_MAX. */

double glp_get_row_ub(glp_prob *lp, int i)
{     double ub;
      if (!(1 <= i && i <= lp->m))
         xerror("glp_get_row_ub: i = %d; row number out of range\n", i);
      switch (lp->row[i]->type)
      {  case GLP_FR:
         case GLP_LO:
            ub = +DBL_MAX; break;
         case GLP_UP:
         case GLP_DB:
         case GLP_FX:
            ub = lp->row[i]->ub; break;
         default:
            xassert(lp != lp);
      }
      return ub;
}

/***********************************************************************
*  NAME
*
*  glp_get_col_type - retrieve column type
*
*  SYNOPSIS
*
*  int glp_get_col_type(glp_prob *lp, int j);
*
*  RETURNS
*
*  The routine glp_get_col_type returns the type of j-th column, i.e.
*  the type of corresponding structural variable, as follows:
*
*  GLP_FR - free (unbounded) variable;
*  GLP_LO - variable with lower bound;
*  GLP_UP - variable with upper bound;
*  GLP_DB - double-bounded variable;
*  GLP_FX - fixed variable. */

int glp_get_col_type(glp_prob *lp, int j)
{     if (!(1 <= j && j <= lp->n))
         xerror("glp_get_col_type: j = %d; column number out of range\n"
            , j);
      return lp->col[j]->type;
}

/***********************************************************************
*  NAME
*
*  glp_get_col_lb - retrieve column lower bound
*
*  SYNOPSIS
*
*  double glp_get_col_lb(glp_prob *lp, int j);
*
*  RETURNS
*
*  The routine glp_get_col_lb returns the lower bound of j-th column,
*  i.e. the lower bound of corresponding structural variable. However,
*  if the column has no lower bound, the routine returns -DBL_MAX. */

double glp_get_col_lb(glp_prob *lp, int j)
{     double lb;
      if (!(1 <= j && j <= lp->n))
         xerror("glp_get_col_lb: j = %d; column number out of range\n",
            j);
      switch (lp->col[j]->type)
      {  case GLP_FR:
         case GLP_UP:
            lb = -DBL_MAX; break;
         case GLP_LO:
         case GLP_DB:
         case GLP_FX:
            lb = lp->col[j]->lb; break;
         default:
            xassert(lp != lp);
      }
      return lb;
}

/***********************************************************************
*  NAME
*
*  glp_get_col_ub - retrieve column upper bound
*
*  SYNOPSIS
*
*  double glp_get_col_ub(glp_prob *lp, int j);
*
*  RETURNS
*
*  The routine glp_get_col_ub returns the upper bound of j-th column,
*  i.e. the upper bound of corresponding structural variable. However,
*  if the column has no upper bound, the routine returns +DBL_MAX. */

double glp_get_col_ub(glp_prob *lp, int j)
{     double ub;
      if (!(1 <= j && j <= lp->n))
         xerror("glp_get_col_ub: j = %d; column number out of range\n",
            j);
      switch (lp->col[j]->type)
      {  case GLP_FR:
         case GLP_LO:
            ub = +DBL_MAX; break;
         case GLP_UP:
         case GLP_DB:
         case GLP_FX:
            ub = lp->col[j]->ub; break;
         default:
            xassert(lp != lp);
      }
      return ub;
}

/***********************************************************************
*  NAME
*
*  glp_get_obj_coef - retrieve obj. coefficient or constant term
*
*  SYNOPSIS
*
*  double glp_get_obj_coef(glp_prob *lp, int j);
*
*  RETURNS
*
*  The routine glp_get_obj_coef returns the objective coefficient at
*  j-th structural variable (column) of the specified problem object.
*
*  If the parameter j is zero, the routine returns the constant term
*  ("shift") of the objective function. */

double glp_get_obj_coef(glp_prob *lp, int j)
{     if (!(0 <= j && j <= lp->n))
         xerror("glp_get_obj_coef: j = %d; column number out of range\n"
            , j);
      return j == 0 ? lp->c0 : lp->col[j]->coef;
}

/***********************************************************************
*  NAME
*
*  glp_get_num_nz - retrieve number of constraint coefficients
*
*  SYNOPSIS
*
*  int glp_get_num_nz(glp_prob *lp);
*
*  RETURNS
*
*  The routine glp_get_num_nz returns the number of (non-zero) elements
*  in the constraint matrix of the specified problem object. */

int glp_get_num_nz(glp_prob *lp)
{     int nnz = lp->nnz;
      return nnz;
}

/***********************************************************************
*  NAME
*
*  glp_get_mat_row - retrieve row of the constraint matrix
*
*  SYNOPSIS
*
*  int glp_get_mat_row(glp_prob *lp, int i, int ind[], double val[]);
*
*  DESCRIPTION
*
*  The routine glp_get_mat_row scans (non-zero) elements of i-th row
*  of the constraint matrix of the specified problem object and stores
*  their column indices and numeric values to locations ind[1], ...,
*  ind[len] and val[1], ..., val[len], respectively, where 0 <= len <= n
*  is the number of elements in i-th row, n is the number of columns.
*
*  The parameter ind and/or val can be specified as NULL, in which case
*  corresponding information is not stored.
*
*  RETURNS
*
*  The routine glp_get_mat_row returns the length len, i.e. the number
*  of (non-zero) elements in i-th row. */

int glp_get_mat_row(glp_prob *lp, int i, int ind[], double val[])
{     GLPAIJ *aij;
      int len;
      if (!(1 <= i && i <= lp->m))
         xerror("glp_get_mat_row: i = %d; row number out of range\n",
            i);
      len = 0;
      for (aij = lp->row[i]->ptr; aij != NULL; aij = aij->r_next)
      {  len++;
         if (ind != NULL) ind[len] = aij->col->j;
         if (val != NULL) val[len] = aij->val;
      }
      xassert(len <= lp->n);
      return len;
}

/***********************************************************************
*  NAME
*
*  glp_get_mat_col - retrieve column of the constraint matrix
*
*  SYNOPSIS
*
*  int glp_get_mat_col(glp_prob *lp, int j, int ind[], double val[]);
*
*  DESCRIPTION
*
*  The routine glp_get_mat_col scans (non-zero) elements of j-th column
*  of the constraint matrix of the specified problem object and stores
*  their row indices and numeric values to locations ind[1], ...,
*  ind[len] and val[1], ..., val[len], respectively, where 0 <= len <= m
*  is the number of elements in j-th column, m is the number of rows.
*
*  The parameter ind or/and val can be specified as NULL, in which case
*  corresponding information is not stored.
*
*  RETURNS
*
*  The routine glp_get_mat_col returns the length len, i.e. the number
*  of (non-zero) elements in j-th column. */

int glp_get_mat_col(glp_prob *lp, int j, int ind[], double val[])
{     GLPAIJ *aij;
      int len;
      if (!(1 <= j && j <= lp->n))
         xerror("glp_get_mat_col: j = %d; column number out of range\n",
            j);
      len = 0;
      for (aij = lp->col[j]->ptr; aij != NULL; aij = aij->c_next)
      {  len++;
         if (ind != NULL) ind[len] = aij->row->i;
         if (val != NULL) val[len] = aij->val;
      }
      xassert(len <= lp->m);
      return len;
}

/* eof */
