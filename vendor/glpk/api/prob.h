/* prob.h (LP/MIP problem object) */

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

#ifndef PROB_H
#define PROB_H

#include "avl.h"
#include "bfd.h"
#include "dmp.h"
#if 1 /* 28/III-2016 */
#define GLP_UNDOC 1
#endif
#include "glpk.h"

typedef struct GLPROW GLPROW;
typedef struct GLPCOL GLPCOL;
typedef struct GLPAIJ GLPAIJ;

#if 0 /* 04/IV-2016 */
#define GLP_PROB_MAGIC 0xD7D9D6C2
#endif

struct glp_prob
{     /* LP/MIP problem object */
#if 0 /* 04/IV-2016 */
      unsigned magic;
      /* magic value used for debugging */
#endif
      DMP *pool;
      /* memory pool to store problem object components */
      glp_tree *tree;
      /* pointer to the search tree; set by the MIP solver when this
         object is used in the tree as a core MIP object */
#if 0 /* 08/III-2014 */
      void *parms;
      /* reserved for backward compatibility */
#endif
      /*--------------------------------------------------------------*/
      /* LP/MIP data */
      char *name;
      /* problem name (1 to 255 chars); NULL means no name is assigned
         to the problem */
      char *obj;
      /* objective function name (1 to 255 chars); NULL means no name
         is assigned to the objective function */
      int dir;
      /* optimization direction flag (objective "sense"):
         GLP_MIN - minimization
         GLP_MAX - maximization */
      double c0;
      /* constant term of the objective function ("shift") */
      int m_max;
      /* length of the array of rows (enlarged automatically) */
      int n_max;
      /* length of the array of columns (enlarged automatically) */
      int m;
      /* number of rows, 0 <= m <= m_max */
      int n;
      /* number of columns, 0 <= n <= n_max */
      int nnz;
      /* number of non-zero constraint coefficients, nnz >= 0 */
      GLPROW **row; /* GLPROW *row[1+m_max]; */
      /* row[i], 1 <= i <= m, is a pointer to i-th row */
      GLPCOL **col; /* GLPCOL *col[1+n_max]; */
      /* col[j], 1 <= j <= n, is a pointer to j-th column */
      AVL *r_tree;
      /* row index to find rows by their names; NULL means this index
         does not exist */
      AVL *c_tree;
      /* column index to find columns by their names; NULL means this
         index does not exist */
      /*--------------------------------------------------------------*/
      /* basis factorization (LP) */
      int valid;
      /* the factorization is valid only if this flag is set */
      int *head; /* int head[1+m_max]; */
      /* basis header (valid only if the factorization is valid);
         head[i] = k is the ordinal number of auxiliary (1 <= k <= m)
         or structural (m+1 <= k <= m+n) variable which corresponds to
         i-th basic variable xB[i], 1 <= i <= m */
#if 0 /* 08/III-2014 */
      glp_bfcp *bfcp;
      /* basis factorization control parameters; may be NULL */
#endif
      BFD *bfd; /* BFD bfd[1:m,1:m]; */
      /* basis factorization driver; may be NULL */
      /*--------------------------------------------------------------*/
      /* basic solution (LP) */
      int pbs_stat;
      /* primal basic solution status:
         GLP_UNDEF  - primal solution is undefined
         GLP_FEAS   - primal solution is feasible
         GLP_INFEAS - primal solution is infeasible
         GLP_NOFEAS - no primal feasible solution exists */
      int dbs_stat;
      /* dual basic solution status:
         GLP_UNDEF  - dual solution is undefined
         GLP_FEAS   - dual solution is feasible
         GLP_INFEAS - dual solution is infeasible
         GLP_NOFEAS - no dual feasible solution exists */
      double obj_val;
      /* objective function value */
      int it_cnt;
      /* simplex method iteration count; increases by one on performing
         one simplex iteration */
      int some;
      /* ordinal number of some auxiliary or structural variable having
         certain property, 0 <= some <= m+n */
      /*--------------------------------------------------------------*/
      /* interior-point solution (LP) */
      int ipt_stat;
      /* interior-point solution status:
         GLP_UNDEF  - interior solution is undefined
         GLP_OPT    - interior solution is optimal
         GLP_INFEAS - interior solution is infeasible
         GLP_NOFEAS - no feasible solution exists */
      double ipt_obj;
      /* objective function value */
      /*--------------------------------------------------------------*/
      /* integer solution (MIP) */
      int mip_stat;
      /* integer solution status:
         GLP_UNDEF  - integer solution is undefined
         GLP_OPT    - integer solution is optimal
         GLP_FEAS   - integer solution is feasible
         GLP_NOFEAS - no integer solution exists */
      double mip_obj;
      /* objective function value */
};

struct GLPROW
{     /* LP/MIP row (auxiliary variable) */
      int i;
      /* ordinal number (1 to m) assigned to this row */
      char *name;
      /* row name (1 to 255 chars); NULL means no name is assigned to
         this row */
      AVLNODE *node;
      /* pointer to corresponding node in the row index; NULL means
         that either the row index does not exist or this row has no
         name assigned */
#if 1 /* 20/IX-2008 */
      int level;
      unsigned char origin;
      unsigned char klass;
#endif
      int type;
      /* type of the auxiliary variable:
         GLP_FR - free variable
         GLP_LO - variable with lower bound
         GLP_UP - variable with upper bound
         GLP_DB - double-bounded variable
         GLP_FX - fixed variable */
      double lb; /* non-scaled */
      /* lower bound; if the row has no lower bound, lb is zero */
      double ub; /* non-scaled */
      /* upper bound; if the row has no upper bound, ub is zero */
      /* if the row type is GLP_FX, ub is equal to lb */
      GLPAIJ *ptr; /* non-scaled */
      /* pointer to doubly linked list of constraint coefficients which
         are placed in this row */
      double rii;
      /* diagonal element r[i,i] of scaling matrix R for this row;
         if the scaling is not used, r[i,i] is 1 */
      int stat;
      /* status of the auxiliary variable:
         GLP_BS - basic variable
         GLP_NL - non-basic variable on lower bound
         GLP_NU - non-basic variable on upper bound
         GLP_NF - non-basic free variable
         GLP_NS - non-basic fixed variable */
      int bind;
      /* if the auxiliary variable is basic, head[bind] refers to this
         row, otherwise, bind is 0; this attribute is valid only if the
         basis factorization is valid */
      double prim; /* non-scaled */
      /* primal value of the auxiliary variable in basic solution */
      double dual; /* non-scaled */
      /* dual value of the auxiliary variable in basic solution */
      double pval; /* non-scaled */
      /* primal value of the auxiliary variable in interior solution */
      double dval; /* non-scaled */
      /* dual value of the auxiliary variable in interior solution */
      double mipx; /* non-scaled */
      /* primal value of the auxiliary variable in integer solution */
};

struct GLPCOL
{     /* LP/MIP column (structural variable) */
      int j;
      /* ordinal number (1 to n) assigned to this column */
      char *name;
      /* column name (1 to 255 chars); NULL means no name is assigned
         to this column */
      AVLNODE *node;
      /* pointer to corresponding node in the column index; NULL means
         that either the column index does not exist or the column has
         no name assigned */
      int kind;
      /* kind of the structural variable:
         GLP_CV - continuous variable
         GLP_IV - integer or binary variable */
      int type;
      /* type of the structural variable:
         GLP_FR - free variable
         GLP_LO - variable with lower bound
         GLP_UP - variable with upper bound
         GLP_DB - double-bounded variable
         GLP_FX - fixed variable */
      double lb; /* non-scaled */
      /* lower bound; if the column has no lower bound, lb is zero */
      double ub; /* non-scaled */
      /* upper bound; if the column has no upper bound, ub is zero */
      /* if the column type is GLP_FX, ub is equal to lb */
      double coef; /* non-scaled */
      /* objective coefficient at the structural variable */
      GLPAIJ *ptr; /* non-scaled */
      /* pointer to doubly linked list of constraint coefficients which
         are placed in this column */
      double sjj;
      /* diagonal element s[j,j] of scaling matrix S for this column;
         if the scaling is not used, s[j,j] is 1 */
      int stat;
      /* status of the structural variable:
         GLP_BS - basic variable
         GLP_NL - non-basic variable on lower bound
         GLP_NU - non-basic variable on upper bound
         GLP_NF - non-basic free variable
         GLP_NS - non-basic fixed variable */
      int bind;
      /* if the structural variable is basic, head[bind] refers to
         this column; otherwise, bind is 0; this attribute is valid only
         if the basis factorization is valid */
      double prim; /* non-scaled */
      /* primal value of the structural variable in basic solution */
      double dual; /* non-scaled */
      /* dual value of the structural variable in basic solution */
      double pval; /* non-scaled */
      /* primal value of the structural variable in interior solution */
      double dval; /* non-scaled */
      /* dual value of the structural variable in interior solution */
      double mipx; /* non-scaled */
      /* primal value of the structural variable in integer solution */
};

struct GLPAIJ
{     /* constraint coefficient a[i,j] */
      GLPROW *row;
      /* pointer to row, where this coefficient is placed */
      GLPCOL *col;
      /* pointer to column, where this coefficient is placed */
      double val;
      /* numeric (non-zero) value of this coefficient */
      GLPAIJ *r_prev;
      /* pointer to previous coefficient in the same row */
      GLPAIJ *r_next;
      /* pointer to next coefficient in the same row */
      GLPAIJ *c_prev;
      /* pointer to previous coefficient in the same column */
      GLPAIJ *c_next;
      /* pointer to next coefficient in the same column */
};

#endif

/* eof */
