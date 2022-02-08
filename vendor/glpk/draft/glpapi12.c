/* glpapi12.c (basis factorization and simplex tableau routines) */

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

#include "draft.h"
#include "env.h"
#include "prob.h"

/***********************************************************************
*  NAME
*
*  glp_bf_exists - check if the basis factorization exists
*
*  SYNOPSIS
*
*  int glp_bf_exists(glp_prob *lp);
*
*  RETURNS
*
*  If the basis factorization for the current basis associated with
*  the specified problem object exists and therefore is available for
*  computations, the routine glp_bf_exists returns non-zero. Otherwise
*  the routine returns zero. */

int glp_bf_exists(glp_prob *lp)
{     int ret;
      ret = (lp->m == 0 || lp->valid);
      return ret;
}

/***********************************************************************
*  NAME
*
*  glp_factorize - compute the basis factorization
*
*  SYNOPSIS
*
*  int glp_factorize(glp_prob *lp);
*
*  DESCRIPTION
*
*  The routine glp_factorize computes the basis factorization for the
*  current basis associated with the specified problem object.
*
*  RETURNS
*
*  0  The basis factorization has been successfully computed.
*
*  GLP_EBADB
*     The basis matrix is invalid, i.e. the number of basic (auxiliary
*     and structural) variables differs from the number of rows in the
*     problem object.
*
*  GLP_ESING
*     The basis matrix is singular within the working precision.
*
*  GLP_ECOND
*     The basis matrix is ill-conditioned. */

static int b_col(void *info, int j, int ind[], double val[])
{     glp_prob *lp = info;
      int m = lp->m;
      GLPAIJ *aij;
      int k, len;
      xassert(1 <= j && j <= m);
      /* determine the ordinal number of basic auxiliary or structural
         variable x[k] corresponding to basic variable xB[j] */
      k = lp->head[j];
      /* build j-th column of the basic matrix, which is k-th column of
         the scaled augmented matrix (I | -R*A*S) */
      if (k <= m)
      {  /* x[k] is auxiliary variable */
         len = 1;
         ind[1] = k;
         val[1] = 1.0;
      }
      else
      {  /* x[k] is structural variable */
         len = 0;
         for (aij = lp->col[k-m]->ptr; aij != NULL; aij = aij->c_next)
         {  len++;
            ind[len] = aij->row->i;
            val[len] = - aij->row->rii * aij->val * aij->col->sjj;
         }
      }
      return len;
}

int glp_factorize(glp_prob *lp)
{     int m = lp->m;
      int n = lp->n;
      GLPROW **row = lp->row;
      GLPCOL **col = lp->col;
      int *head = lp->head;
      int j, k, stat, ret;
      /* invalidate the basis factorization */
      lp->valid = 0;
      /* build the basis header */
      j = 0;
      for (k = 1; k <= m+n; k++)
      {  if (k <= m)
         {  stat = row[k]->stat;
            row[k]->bind = 0;
         }
         else
         {  stat = col[k-m]->stat;
            col[k-m]->bind = 0;
         }
         if (stat == GLP_BS)
         {  j++;
            if (j > m)
            {  /* too many basic variables */
               ret = GLP_EBADB;
               goto fini;
            }
            head[j] = k;
            if (k <= m)
               row[k]->bind = j;
            else
               col[k-m]->bind = j;
         }
      }
      if (j < m)
      {  /* too few basic variables */
         ret = GLP_EBADB;
         goto fini;
      }
      /* try to factorize the basis matrix */
      if (m > 0)
      {  if (lp->bfd == NULL)
         {  lp->bfd = bfd_create_it();
#if 0 /* 08/III-2014 */
            copy_bfcp(lp);
#endif
         }
         switch (bfd_factorize(lp->bfd, m, /*lp->head,*/ b_col, lp))
         {  case 0:
               /* ok */
               break;
            case BFD_ESING:
               /* singular matrix */
               ret = GLP_ESING;
               goto fini;
            case BFD_ECOND:
               /* ill-conditioned matrix */
               ret = GLP_ECOND;
               goto fini;
            default:
               xassert(lp != lp);
         }
         lp->valid = 1;
      }
      /* factorization successful */
      ret = 0;
fini: /* bring the return code to the calling program */
      return ret;
}

/***********************************************************************
*  NAME
*
*  glp_bf_updated - check if the basis factorization has been updated
*
*  SYNOPSIS
*
*  int glp_bf_updated(glp_prob *lp);
*
*  RETURNS
*
*  If the basis factorization has been just computed from scratch, the
*  routine glp_bf_updated returns zero. Otherwise, if the factorization
*  has been updated one or more times, the routine returns non-zero. */

int glp_bf_updated(glp_prob *lp)
{     int cnt;
      if (!(lp->m == 0 || lp->valid))
         xerror("glp_bf_update: basis factorization does not exist\n");
#if 0 /* 15/XI-2009 */
      cnt = (lp->m == 0 ? 0 : lp->bfd->upd_cnt);
#else
      cnt = (lp->m == 0 ? 0 : bfd_get_count(lp->bfd));
#endif
      return cnt;
}

/***********************************************************************
*  NAME
*
*  glp_get_bfcp - retrieve basis factorization control parameters
*
*  SYNOPSIS
*
*  void glp_get_bfcp(glp_prob *lp, glp_bfcp *parm);
*
*  DESCRIPTION
*
*  The routine glp_get_bfcp retrieves control parameters, which are
*  used on computing and updating the basis factorization associated
*  with the specified problem object.
*
*  Current values of control parameters are stored by the routine in
*  a glp_bfcp structure, which the parameter parm points to. */

#if 1 /* 08/III-2014 */
void glp_get_bfcp(glp_prob *P, glp_bfcp *parm)
{     if (P->bfd == NULL)
         P->bfd = bfd_create_it();
      bfd_get_bfcp(P->bfd, parm);
      return;
}
#endif

/***********************************************************************
*  NAME
*
*  glp_set_bfcp - change basis factorization control parameters
*
*  SYNOPSIS
*
*  void glp_set_bfcp(glp_prob *lp, const glp_bfcp *parm);
*
*  DESCRIPTION
*
*  The routine glp_set_bfcp changes control parameters, which are used
*  by internal GLPK routines in computing and updating the basis
*  factorization associated with the specified problem object.
*
*  New values of the control parameters should be passed in a structure
*  glp_bfcp, which the parameter parm points to.
*
*  The parameter parm can be specified as NULL, in which case all
*  control parameters are reset to their default values. */

#if 1 /* 08/III-2014 */
void glp_set_bfcp(glp_prob *P, const glp_bfcp *parm)
{     if (P->bfd == NULL)
         P->bfd = bfd_create_it();
      if (parm != NULL)
      {  if (!(parm->type == GLP_BF_LUF + GLP_BF_FT ||
               parm->type == GLP_BF_LUF + GLP_BF_BG ||
               parm->type == GLP_BF_LUF + GLP_BF_GR ||
               parm->type == GLP_BF_BTF + GLP_BF_BG ||
               parm->type == GLP_BF_BTF + GLP_BF_GR))
            xerror("glp_set_bfcp: type = 0x%02X; invalid parameter\n",
               parm->type);
         if (!(0.0 < parm->piv_tol && parm->piv_tol < 1.0))
            xerror("glp_set_bfcp: piv_tol = %g; invalid parameter\n",
               parm->piv_tol);
         if (parm->piv_lim < 1)
            xerror("glp_set_bfcp: piv_lim = %d; invalid parameter\n",
               parm->piv_lim);
         if (!(parm->suhl == GLP_ON || parm->suhl == GLP_OFF))
            xerror("glp_set_bfcp: suhl = %d; invalid parameter\n",
               parm->suhl);
         if (!(0.0 <= parm->eps_tol && parm->eps_tol <= 1e-6))
            xerror("glp_set_bfcp: eps_tol = %g; invalid parameter\n",
               parm->eps_tol);
         if (!(1 <= parm->nfs_max && parm->nfs_max <= 32767))
            xerror("glp_set_bfcp: nfs_max = %d; invalid parameter\n",
               parm->nfs_max);
         if (!(1 <= parm->nrs_max && parm->nrs_max <= 32767))
            xerror("glp_set_bfcp: nrs_max = %d; invalid parameter\n",
               parm->nrs_max);
      }
      bfd_set_bfcp(P->bfd, parm);
      return;
}
#endif

/***********************************************************************
*  NAME
*
*  glp_get_bhead - retrieve the basis header information
*
*  SYNOPSIS
*
*  int glp_get_bhead(glp_prob *lp, int k);
*
*  DESCRIPTION
*
*  The routine glp_get_bhead returns the basis header information for
*  the current basis associated with the specified problem object.
*
*  RETURNS
*
*  If xB[k], 1 <= k <= m, is i-th auxiliary variable (1 <= i <= m), the
*  routine returns i. Otherwise, if xB[k] is j-th structural variable
*  (1 <= j <= n), the routine returns m+j. Here m is the number of rows
*  and n is the number of columns in the problem object. */

int glp_get_bhead(glp_prob *lp, int k)
{     if (!(lp->m == 0 || lp->valid))
         xerror("glp_get_bhead: basis factorization does not exist\n");
      if (!(1 <= k && k <= lp->m))
         xerror("glp_get_bhead: k = %d; index out of range\n", k);
      return lp->head[k];
}

/***********************************************************************
*  NAME
*
*  glp_get_row_bind - retrieve row index in the basis header
*
*  SYNOPSIS
*
*  int glp_get_row_bind(glp_prob *lp, int i);
*
*  RETURNS
*
*  The routine glp_get_row_bind returns the index k of basic variable
*  xB[k], 1 <= k <= m, which is i-th auxiliary variable, 1 <= i <= m,
*  in the current basis associated with the specified problem object,
*  where m is the number of rows. However, if i-th auxiliary variable
*  is non-basic, the routine returns zero. */

int glp_get_row_bind(glp_prob *lp, int i)
{     if (!(lp->m == 0 || lp->valid))
         xerror("glp_get_row_bind: basis factorization does not exist\n"
            );
      if (!(1 <= i && i <= lp->m))
         xerror("glp_get_row_bind: i = %d; row number out of range\n",
            i);
      return lp->row[i]->bind;
}

/***********************************************************************
*  NAME
*
*  glp_get_col_bind - retrieve column index in the basis header
*
*  SYNOPSIS
*
*  int glp_get_col_bind(glp_prob *lp, int j);
*
*  RETURNS
*
*  The routine glp_get_col_bind returns the index k of basic variable
*  xB[k], 1 <= k <= m, which is j-th structural variable, 1 <= j <= n,
*  in the current basis associated with the specified problem object,
*  where m is the number of rows, n is the number of columns. However,
*  if j-th structural variable is non-basic, the routine returns zero.*/

int glp_get_col_bind(glp_prob *lp, int j)
{     if (!(lp->m == 0 || lp->valid))
         xerror("glp_get_col_bind: basis factorization does not exist\n"
            );
      if (!(1 <= j && j <= lp->n))
         xerror("glp_get_col_bind: j = %d; column number out of range\n"
            , j);
      return lp->col[j]->bind;
}

/***********************************************************************
*  NAME
*
*  glp_ftran - perform forward transformation (solve system B*x = b)
*
*  SYNOPSIS
*
*  void glp_ftran(glp_prob *lp, double x[]);
*
*  DESCRIPTION
*
*  The routine glp_ftran performs forward transformation, i.e. solves
*  the system B*x = b, where B is the basis matrix corresponding to the
*  current basis for the specified problem object, x is the vector of
*  unknowns to be computed, b is the vector of right-hand sides.
*
*  On entry elements of the vector b should be stored in dense format
*  in locations x[1], ..., x[m], where m is the number of rows. On exit
*  the routine stores elements of the vector x in the same locations.
*
*  SCALING/UNSCALING
*
*  Let A~ = (I | -A) is the augmented constraint matrix of the original
*  (unscaled) problem. In the scaled LP problem instead the matrix A the
*  scaled matrix A" = R*A*S is actually used, so
*
*     A~" = (I | A") = (I | R*A*S) = (R*I*inv(R) | R*A*S) =
*                                                                    (1)
*         = R*(I | A)*S~ = R*A~*S~,
*
*  is the scaled augmented constraint matrix, where R and S are diagonal
*  scaling matrices used to scale rows and columns of the matrix A, and
*
*     S~ = diag(inv(R) | S)                                          (2)
*
*  is an augmented diagonal scaling matrix.
*
*  By definition:
*
*     A~ = (B | N),                                                  (3)
*
*  where B is the basic matrix, which consists of basic columns of the
*  augmented constraint matrix A~, and N is a matrix, which consists of
*  non-basic columns of A~. From (1) it follows that:
*
*     A~" = (B" | N") = (R*B*SB | R*N*SN),                           (4)
*
*  where SB and SN are parts of the augmented scaling matrix S~, which
*  correspond to basic and non-basic variables, respectively. Therefore
*
*     B" = R*B*SB,                                                   (5)
*
*  which is the scaled basis matrix. */

void glp_ftran(glp_prob *lp, double x[])
{     int m = lp->m;
      GLPROW **row = lp->row;
      GLPCOL **col = lp->col;
      int i, k;
      /* B*x = b ===> (R*B*SB)*(inv(SB)*x) = R*b ===>
         B"*x" = b", where b" = R*b, x = SB*x" */
      if (!(m == 0 || lp->valid))
         xerror("glp_ftran: basis factorization does not exist\n");
      /* b" := R*b */
      for (i = 1; i <= m; i++)
         x[i] *= row[i]->rii;
      /* x" := inv(B")*b" */
      if (m > 0) bfd_ftran(lp->bfd, x);
      /* x := SB*x" */
      for (i = 1; i <= m; i++)
      {  k = lp->head[i];
         if (k <= m)
            x[i] /= row[k]->rii;
         else
            x[i] *= col[k-m]->sjj;
      }
      return;
}

/***********************************************************************
*  NAME
*
*  glp_btran - perform backward transformation (solve system B'*x = b)
*
*  SYNOPSIS
*
*  void glp_btran(glp_prob *lp, double x[]);
*
*  DESCRIPTION
*
*  The routine glp_btran performs backward transformation, i.e. solves
*  the system B'*x = b, where B' is a matrix transposed to the basis
*  matrix corresponding to the current basis for the specified problem
*  problem object, x is the vector of unknowns to be computed, b is the
*  vector of right-hand sides.
*
*  On entry elements of the vector b should be stored in dense format
*  in locations x[1], ..., x[m], where m is the number of rows. On exit
*  the routine stores elements of the vector x in the same locations.
*
*  SCALING/UNSCALING
*
*  See comments to the routine glp_ftran. */

void glp_btran(glp_prob *lp, double x[])
{     int m = lp->m;
      GLPROW **row = lp->row;
      GLPCOL **col = lp->col;
      int i, k;
      /* B'*x = b ===> (SB*B'*R)*(inv(R)*x) = SB*b ===>
         (B")'*x" = b", where b" = SB*b, x = R*x" */
      if (!(m == 0 || lp->valid))
         xerror("glp_btran: basis factorization does not exist\n");
      /* b" := SB*b */
      for (i = 1; i <= m; i++)
      {  k = lp->head[i];
         if (k <= m)
            x[i] /= row[k]->rii;
         else
            x[i] *= col[k-m]->sjj;
      }
      /* x" := inv[(B")']*b" */
      if (m > 0) bfd_btran(lp->bfd, x);
      /* x := R*x" */
      for (i = 1; i <= m; i++)
         x[i] *= row[i]->rii;
      return;
}

/***********************************************************************
*  NAME
*
*  glp_warm_up - "warm up" LP basis
*
*  SYNOPSIS
*
*  int glp_warm_up(glp_prob *P);
*
*  DESCRIPTION
*
*  The routine glp_warm_up "warms up" the LP basis for the specified
*  problem object using current statuses assigned to rows and columns
*  (that is, to auxiliary and structural variables).
*
*  This operation includes computing factorization of the basis matrix
*  (if it does not exist), computing primal and dual components of basic
*  solution, and determining the solution status.
*
*  RETURNS
*
*  0  The operation has been successfully performed.
*
*  GLP_EBADB
*     The basis matrix is invalid, i.e. the number of basic (auxiliary
*     and structural) variables differs from the number of rows in the
*     problem object.
*
*  GLP_ESING
*     The basis matrix is singular within the working precision.
*
*  GLP_ECOND
*     The basis matrix is ill-conditioned. */

int glp_warm_up(glp_prob *P)
{     GLPROW *row;
      GLPCOL *col;
      GLPAIJ *aij;
      int i, j, type, stat, ret;
      double eps, temp, *work;
      /* invalidate basic solution */
      P->pbs_stat = P->dbs_stat = GLP_UNDEF;
      P->obj_val = 0.0;
      P->some = 0;
      for (i = 1; i <= P->m; i++)
      {  row = P->row[i];
         row->prim = row->dual = 0.0;
      }
      for (j = 1; j <= P->n; j++)
      {  col = P->col[j];
         col->prim = col->dual = 0.0;
      }
      /* compute the basis factorization, if necessary */
      if (!glp_bf_exists(P))
      {  ret = glp_factorize(P);
         if (ret != 0) goto done;
      }
      /* allocate working array */
      work = xcalloc(1+P->m, sizeof(double));
      /* determine and store values of non-basic variables, compute
         vector (- N * xN) */
      for (i = 1; i <= P->m; i++)
         work[i] = 0.0;
      for (i = 1; i <= P->m; i++)
      {  row = P->row[i];
         if (row->stat == GLP_BS)
            continue;
         else if (row->stat == GLP_NL)
            row->prim = row->lb;
         else if (row->stat == GLP_NU)
            row->prim = row->ub;
         else if (row->stat == GLP_NF)
            row->prim = 0.0;
         else if (row->stat == GLP_NS)
            row->prim = row->lb;
         else
            xassert(row != row);
         /* N[j] is i-th column of matrix (I|-A) */
         work[i] -= row->prim;
      }
      for (j = 1; j <= P->n; j++)
      {  col = P->col[j];
         if (col->stat == GLP_BS)
            continue;
         else if (col->stat == GLP_NL)
            col->prim = col->lb;
         else if (col->stat == GLP_NU)
            col->prim = col->ub;
         else if (col->stat == GLP_NF)
            col->prim = 0.0;
         else if (col->stat == GLP_NS)
            col->prim = col->lb;
         else
            xassert(col != col);
         /* N[j] is (m+j)-th column of matrix (I|-A) */
         if (col->prim != 0.0)
         {  for (aij = col->ptr; aij != NULL; aij = aij->c_next)
               work[aij->row->i] += aij->val * col->prim;
         }
      }
      /* compute vector of basic variables xB = - inv(B) * N * xN */
      glp_ftran(P, work);
      /* store values of basic variables, check primal feasibility */
      P->pbs_stat = GLP_FEAS;
      for (i = 1; i <= P->m; i++)
      {  row = P->row[i];
         if (row->stat != GLP_BS)
            continue;
         row->prim = work[row->bind];
         type = row->type;
         if (type == GLP_LO || type == GLP_DB || type == GLP_FX)
         {  eps = 1e-6 + 1e-9 * fabs(row->lb);
            if (row->prim < row->lb - eps)
               P->pbs_stat = GLP_INFEAS;
         }
         if (type == GLP_UP || type == GLP_DB || type == GLP_FX)
         {  eps = 1e-6 + 1e-9 * fabs(row->ub);
            if (row->prim > row->ub + eps)
               P->pbs_stat = GLP_INFEAS;
         }
      }
      for (j = 1; j <= P->n; j++)
      {  col = P->col[j];
         if (col->stat != GLP_BS)
            continue;
         col->prim = work[col->bind];
         type = col->type;
         if (type == GLP_LO || type == GLP_DB || type == GLP_FX)
         {  eps = 1e-6 + 1e-9 * fabs(col->lb);
            if (col->prim < col->lb - eps)
               P->pbs_stat = GLP_INFEAS;
         }
         if (type == GLP_UP || type == GLP_DB || type == GLP_FX)
         {  eps = 1e-6 + 1e-9 * fabs(col->ub);
            if (col->prim > col->ub + eps)
               P->pbs_stat = GLP_INFEAS;
         }
      }
      /* compute value of the objective function */
      P->obj_val = P->c0;
      for (j = 1; j <= P->n; j++)
      {  col = P->col[j];
         P->obj_val += col->coef * col->prim;
      }
      /* build vector cB of objective coefficients at basic variables */
      for (i = 1; i <= P->m; i++)
         work[i] = 0.0;
      for (j = 1; j <= P->n; j++)
      {  col = P->col[j];
         if (col->stat == GLP_BS)
            work[col->bind] = col->coef;
      }
      /* compute vector of simplex multipliers pi = inv(B') * cB */
      glp_btran(P, work);
      /* compute and store reduced costs of non-basic variables d[j] =
         c[j] - N'[j] * pi, check dual feasibility */
      P->dbs_stat = GLP_FEAS;
      for (i = 1; i <= P->m; i++)
      {  row = P->row[i];
         if (row->stat == GLP_BS)
         {  row->dual = 0.0;
            continue;
         }
         /* N[j] is i-th column of matrix (I|-A) */
         row->dual = - work[i];
#if 0 /* 07/III-2013 */
         type = row->type;
         temp = (P->dir == GLP_MIN ? + row->dual : - row->dual);
         if ((type == GLP_FR || type == GLP_LO) && temp < -1e-5 ||
             (type == GLP_FR || type == GLP_UP) && temp > +1e-5)
            P->dbs_stat = GLP_INFEAS;
#else
         stat = row->stat;
         temp = (P->dir == GLP_MIN ? + row->dual : - row->dual);
         if ((stat == GLP_NF || stat == GLP_NL) && temp < -1e-5 ||
             (stat == GLP_NF || stat == GLP_NU) && temp > +1e-5)
            P->dbs_stat = GLP_INFEAS;
#endif
      }
      for (j = 1; j <= P->n; j++)
      {  col = P->col[j];
         if (col->stat == GLP_BS)
         {  col->dual = 0.0;
            continue;
         }
         /* N[j] is (m+j)-th column of matrix (I|-A) */
         col->dual = col->coef;
         for (aij = col->ptr; aij != NULL; aij = aij->c_next)
            col->dual += aij->val * work[aij->row->i];
#if 0 /* 07/III-2013 */
         type = col->type;
         temp = (P->dir == GLP_MIN ? + col->dual : - col->dual);
         if ((type == GLP_FR || type == GLP_LO) && temp < -1e-5 ||
             (type == GLP_FR || type == GLP_UP) && temp > +1e-5)
            P->dbs_stat = GLP_INFEAS;
#else
         stat = col->stat;
         temp = (P->dir == GLP_MIN ? + col->dual : - col->dual);
         if ((stat == GLP_NF || stat == GLP_NL) && temp < -1e-5 ||
             (stat == GLP_NF || stat == GLP_NU) && temp > +1e-5)
            P->dbs_stat = GLP_INFEAS;
#endif
      }
      /* free working array */
      xfree(work);
      ret = 0;
done: return ret;
}

/***********************************************************************
*  NAME
*
*  glp_eval_tab_row - compute row of the simplex tableau
*
*  SYNOPSIS
*
*  int glp_eval_tab_row(glp_prob *lp, int k, int ind[], double val[]);
*
*  DESCRIPTION
*
*  The routine glp_eval_tab_row computes a row of the current simplex
*  tableau for the basic variable, which is specified by the number k:
*  if 1 <= k <= m, x[k] is k-th auxiliary variable; if m+1 <= k <= m+n,
*  x[k] is (k-m)-th structural variable, where m is number of rows, and
*  n is number of columns. The current basis must be available.
*
*  The routine stores column indices and numerical values of non-zero
*  elements of the computed row using sparse format to the locations
*  ind[1], ..., ind[len] and val[1], ..., val[len], respectively, where
*  0 <= len <= n is number of non-zeros returned on exit.
*
*  Element indices stored in the array ind have the same sense as the
*  index k, i.e. indices 1 to m denote auxiliary variables and indices
*  m+1 to m+n denote structural ones (all these variables are obviously
*  non-basic by definition).
*
*  The computed row shows how the specified basic variable x[k] = xB[i]
*  depends on non-basic variables:
*
*     xB[i] = alfa[i,1]*xN[1] + alfa[i,2]*xN[2] + ... + alfa[i,n]*xN[n],
*
*  where alfa[i,j] are elements of the simplex table row, xN[j] are
*  non-basic (auxiliary and structural) variables.
*
*  RETURNS
*
*  The routine returns number of non-zero elements in the simplex table
*  row stored in the arrays ind and val.
*
*  BACKGROUND
*
*  The system of equality constraints of the LP problem is:
*
*     xR = A * xS,                                                   (1)
*
*  where xR is the vector of auxliary variables, xS is the vector of
*  structural variables, A is the matrix of constraint coefficients.
*
*  The system (1) can be written in homogenous form as follows:
*
*     A~ * x = 0,                                                    (2)
*
*  where A~ = (I | -A) is the augmented constraint matrix (has m rows
*  and m+n columns), x = (xR | xS) is the vector of all (auxiliary and
*  structural) variables.
*
*  By definition for the current basis we have:
*
*     A~ = (B | N),                                                  (3)
*
*  where B is the basis matrix. Thus, the system (2) can be written as:
*
*     B * xB + N * xN = 0.                                           (4)
*
*  From (4) it follows that:
*
*     xB = A^ * xN,                                                  (5)
*
*  where the matrix
*
*     A^ = - inv(B) * N                                              (6)
*
*  is called the simplex table.
*
*  It is understood that i-th row of the simplex table is:
*
*     e * A^ = - e * inv(B) * N,                                     (7)
*
*  where e is a unity vector with e[i] = 1.
*
*  To compute i-th row of the simplex table the routine first computes
*  i-th row of the inverse:
*
*     rho = inv(B') * e,                                             (8)
*
*  where B' is a matrix transposed to B, and then computes elements of
*  i-th row of the simplex table as scalar products:
*
*     alfa[i,j] = - rho * N[j]   for all j,                          (9)
*
*  where N[j] is a column of the augmented constraint matrix A~, which
*  corresponds to some non-basic auxiliary or structural variable. */

int glp_eval_tab_row(glp_prob *lp, int k, int ind[], double val[])
{     int m = lp->m;
      int n = lp->n;
      int i, t, len, lll, *iii;
      double alfa, *rho, *vvv;
      if (!(m == 0 || lp->valid))
         xerror("glp_eval_tab_row: basis factorization does not exist\n"
            );
      if (!(1 <= k && k <= m+n))
         xerror("glp_eval_tab_row: k = %d; variable number out of range"
            , k);
      /* determine xB[i] which corresponds to x[k] */
      if (k <= m)
         i = glp_get_row_bind(lp, k);
      else
         i = glp_get_col_bind(lp, k-m);
      if (i == 0)
         xerror("glp_eval_tab_row: k = %d; variable must be basic", k);
      xassert(1 <= i && i <= m);
      /* allocate working arrays */
      rho = xcalloc(1+m, sizeof(double));
      iii = xcalloc(1+m, sizeof(int));
      vvv = xcalloc(1+m, sizeof(double));
      /* compute i-th row of the inverse; see (8) */
      for (t = 1; t <= m; t++) rho[t] = 0.0;
      rho[i] = 1.0;
      glp_btran(lp, rho);
      /* compute i-th row of the simplex table */
      len = 0;
      for (k = 1; k <= m+n; k++)
      {  if (k <= m)
         {  /* x[k] is auxiliary variable, so N[k] is a unity column */
            if (glp_get_row_stat(lp, k) == GLP_BS) continue;
            /* compute alfa[i,j]; see (9) */
            alfa = - rho[k];
         }
         else
         {  /* x[k] is structural variable, so N[k] is a column of the
               original constraint matrix A with negative sign */
            if (glp_get_col_stat(lp, k-m) == GLP_BS) continue;
            /* compute alfa[i,j]; see (9) */
            lll = glp_get_mat_col(lp, k-m, iii, vvv);
            alfa = 0.0;
            for (t = 1; t <= lll; t++) alfa += rho[iii[t]] * vvv[t];
         }
         /* store alfa[i,j] */
         if (alfa != 0.0) len++, ind[len] = k, val[len] = alfa;
      }
      xassert(len <= n);
      /* free working arrays */
      xfree(rho);
      xfree(iii);
      xfree(vvv);
      /* return to the calling program */
      return len;
}

/***********************************************************************
*  NAME
*
*  glp_eval_tab_col - compute column of the simplex tableau
*
*  SYNOPSIS
*
*  int glp_eval_tab_col(glp_prob *lp, int k, int ind[], double val[]);
*
*  DESCRIPTION
*
*  The routine glp_eval_tab_col computes a column of the current simplex
*  table for the non-basic variable, which is specified by the number k:
*  if 1 <= k <= m, x[k] is k-th auxiliary variable; if m+1 <= k <= m+n,
*  x[k] is (k-m)-th structural variable, where m is number of rows, and
*  n is number of columns. The current basis must be available.
*
*  The routine stores row indices and numerical values of non-zero
*  elements of the computed column using sparse format to the locations
*  ind[1], ..., ind[len] and val[1], ..., val[len] respectively, where
*  0 <= len <= m is number of non-zeros returned on exit.
*
*  Element indices stored in the array ind have the same sense as the
*  index k, i.e. indices 1 to m denote auxiliary variables and indices
*  m+1 to m+n denote structural ones (all these variables are obviously
*  basic by the definition).
*
*  The computed column shows how basic variables depend on the specified
*  non-basic variable x[k] = xN[j]:
*
*     xB[1] = ... + alfa[1,j]*xN[j] + ...
*     xB[2] = ... + alfa[2,j]*xN[j] + ...
*              . . . . . .
*     xB[m] = ... + alfa[m,j]*xN[j] + ...
*
*  where alfa[i,j] are elements of the simplex table column, xB[i] are
*  basic (auxiliary and structural) variables.
*
*  RETURNS
*
*  The routine returns number of non-zero elements in the simplex table
*  column stored in the arrays ind and val.
*
*  BACKGROUND
*
*  As it was explained in comments to the routine glp_eval_tab_row (see
*  above) the simplex table is the following matrix:
*
*     A^ = - inv(B) * N.                                             (1)
*
*  Therefore j-th column of the simplex table is:
*
*     A^ * e = - inv(B) * N * e = - inv(B) * N[j],                   (2)
*
*  where e is a unity vector with e[j] = 1, B is the basis matrix, N[j]
*  is a column of the augmented constraint matrix A~, which corresponds
*  to the given non-basic auxiliary or structural variable. */

int glp_eval_tab_col(glp_prob *lp, int k, int ind[], double val[])
{     int m = lp->m;
      int n = lp->n;
      int t, len, stat;
      double *col;
      if (!(m == 0 || lp->valid))
         xerror("glp_eval_tab_col: basis factorization does not exist\n"
            );
      if (!(1 <= k && k <= m+n))
         xerror("glp_eval_tab_col: k = %d; variable number out of range"
            , k);
      if (k <= m)
         stat = glp_get_row_stat(lp, k);
      else
         stat = glp_get_col_stat(lp, k-m);
      if (stat == GLP_BS)
         xerror("glp_eval_tab_col: k = %d; variable must be non-basic",
            k);
      /* obtain column N[k] with negative sign */
      col = xcalloc(1+m, sizeof(double));
      for (t = 1; t <= m; t++) col[t] = 0.0;
      if (k <= m)
      {  /* x[k] is auxiliary variable, so N[k] is a unity column */
         col[k] = -1.0;
      }
      else
      {  /* x[k] is structural variable, so N[k] is a column of the
            original constraint matrix A with negative sign */
         len = glp_get_mat_col(lp, k-m, ind, val);
         for (t = 1; t <= len; t++) col[ind[t]] = val[t];
      }
      /* compute column of the simplex table, which corresponds to the
         specified non-basic variable x[k] */
      glp_ftran(lp, col);
      len = 0;
      for (t = 1; t <= m; t++)
      {  if (col[t] != 0.0)
         {  len++;
            ind[len] = glp_get_bhead(lp, t);
            val[len] = col[t];
         }
      }
      xfree(col);
      /* return to the calling program */
      return len;
}

/***********************************************************************
*  NAME
*
*  glp_transform_row - transform explicitly specified row
*
*  SYNOPSIS
*
*  int glp_transform_row(glp_prob *P, int len, int ind[], double val[]);
*
*  DESCRIPTION
*
*  The routine glp_transform_row performs the same operation as the
*  routine glp_eval_tab_row with exception that the row to be
*  transformed is specified explicitly as a sparse vector.
*
*  The explicitly specified row may be thought as a linear form:
*
*     x = a[1]*x[m+1] + a[2]*x[m+2] + ... + a[n]*x[m+n],             (1)
*
*  where x is an auxiliary variable for this row, a[j] are coefficients
*  of the linear form, x[m+j] are structural variables.
*
*  On entry column indices and numerical values of non-zero elements of
*  the row should be stored in locations ind[1], ..., ind[len] and
*  val[1], ..., val[len], where len is the number of non-zero elements.
*
*  This routine uses the system of equality constraints and the current
*  basis in order to express the auxiliary variable x in (1) through the
*  current non-basic variables (as if the transformed row were added to
*  the problem object and its auxiliary variable were basic), i.e. the
*  resultant row has the form:
*
*     x = alfa[1]*xN[1] + alfa[2]*xN[2] + ... + alfa[n]*xN[n],       (2)
*
*  where xN[j] are non-basic (auxiliary or structural) variables, n is
*  the number of columns in the LP problem object.
*
*  On exit the routine stores indices and numerical values of non-zero
*  elements of the resultant row (2) in locations ind[1], ..., ind[len']
*  and val[1], ..., val[len'], where 0 <= len' <= n is the number of
*  non-zero elements in the resultant row returned by the routine. Note
*  that indices (numbers) of non-basic variables stored in the array ind
*  correspond to original ordinal numbers of variables: indices 1 to m
*  mean auxiliary variables and indices m+1 to m+n mean structural ones.
*
*  RETURNS
*
*  The routine returns len', which is the number of non-zero elements in
*  the resultant row stored in the arrays ind and val.
*
*  BACKGROUND
*
*  The explicitly specified row (1) is transformed in the same way as it
*  were the objective function row.
*
*  From (1) it follows that:
*
*     x = aB * xB + aN * xN,                                         (3)
*
*  where xB is the vector of basic variables, xN is the vector of
*  non-basic variables.
*
*  The simplex table, which corresponds to the current basis, is:
*
*     xB = [-inv(B) * N] * xN.                                       (4)
*
*  Therefore substituting xB from (4) to (3) we have:
*
*     x = aB * [-inv(B) * N] * xN + aN * xN =
*                                                                    (5)
*       = rho * (-N) * xN + aN * xN = alfa * xN,
*
*  where:
*
*     rho = inv(B') * aB,                                            (6)
*
*  and
*
*     alfa = aN + rho * (-N)                                         (7)
*
*  is the resultant row computed by the routine. */

int glp_transform_row(glp_prob *P, int len, int ind[], double val[])
{     int i, j, k, m, n, t, lll, *iii;
      double alfa, *a, *aB, *rho, *vvv;
      if (!glp_bf_exists(P))
         xerror("glp_transform_row: basis factorization does not exist "
            "\n");
      m = glp_get_num_rows(P);
      n = glp_get_num_cols(P);
      /* unpack the row to be transformed to the array a */
      a = xcalloc(1+n, sizeof(double));
      for (j = 1; j <= n; j++) a[j] = 0.0;
      if (!(0 <= len && len <= n))
         xerror("glp_transform_row: len = %d; invalid row length\n",
            len);
      for (t = 1; t <= len; t++)
      {  j = ind[t];
         if (!(1 <= j && j <= n))
            xerror("glp_transform_row: ind[%d] = %d; column index out o"
               "f range\n", t, j);
         if (val[t] == 0.0)
            xerror("glp_transform_row: val[%d] = 0; zero coefficient no"
               "t allowed\n", t);
         if (a[j] != 0.0)
            xerror("glp_transform_row: ind[%d] = %d; duplicate column i"
               "ndices not allowed\n", t, j);
         a[j] = val[t];
      }
      /* construct the vector aB */
      aB = xcalloc(1+m, sizeof(double));
      for (i = 1; i <= m; i++)
      {  k = glp_get_bhead(P, i);
         /* xB[i] is k-th original variable */
         xassert(1 <= k && k <= m+n);
         aB[i] = (k <= m ? 0.0 : a[k-m]);
      }
      /* solve the system B'*rho = aB to compute the vector rho */
      rho = aB, glp_btran(P, rho);
      /* compute coefficients at non-basic auxiliary variables */
      len = 0;
      for (i = 1; i <= m; i++)
      {  if (glp_get_row_stat(P, i) != GLP_BS)
         {  alfa = - rho[i];
            if (alfa != 0.0)
            {  len++;
               ind[len] = i;
               val[len] = alfa;
            }
         }
      }
      /* compute coefficients at non-basic structural variables */
      iii = xcalloc(1+m, sizeof(int));
      vvv = xcalloc(1+m, sizeof(double));
      for (j = 1; j <= n; j++)
      {  if (glp_get_col_stat(P, j) != GLP_BS)
         {  alfa = a[j];
            lll = glp_get_mat_col(P, j, iii, vvv);
            for (t = 1; t <= lll; t++) alfa += vvv[t] * rho[iii[t]];
            if (alfa != 0.0)
            {  len++;
               ind[len] = m+j;
               val[len] = alfa;
            }
         }
      }
      xassert(len <= n);
      xfree(iii);
      xfree(vvv);
      xfree(aB);
      xfree(a);
      return len;
}

/***********************************************************************
*  NAME
*
*  glp_transform_col - transform explicitly specified column
*
*  SYNOPSIS
*
*  int glp_transform_col(glp_prob *P, int len, int ind[], double val[]);
*
*  DESCRIPTION
*
*  The routine glp_transform_col performs the same operation as the
*  routine glp_eval_tab_col with exception that the column to be
*  transformed is specified explicitly as a sparse vector.
*
*  The explicitly specified column may be thought as if it were added
*  to the original system of equality constraints:
*
*     x[1] = a[1,1]*x[m+1] + ... + a[1,n]*x[m+n] + a[1]*x
*     x[2] = a[2,1]*x[m+1] + ... + a[2,n]*x[m+n] + a[2]*x            (1)
*        .  .  .  .  .  .  .  .  .  .  .  .  .  .  .
*     x[m] = a[m,1]*x[m+1] + ... + a[m,n]*x[m+n] + a[m]*x
*
*  where x[i] are auxiliary variables, x[m+j] are structural variables,
*  x is a structural variable for the explicitly specified column, a[i]
*  are constraint coefficients for x.
*
*  On entry row indices and numerical values of non-zero elements of
*  the column should be stored in locations ind[1], ..., ind[len] and
*  val[1], ..., val[len], where len is the number of non-zero elements.
*
*  This routine uses the system of equality constraints and the current
*  basis in order to express the current basic variables through the
*  structural variable x in (1) (as if the transformed column were added
*  to the problem object and the variable x were non-basic), i.e. the
*  resultant column has the form:
*
*     xB[1] = ... + alfa[1]*x
*     xB[2] = ... + alfa[2]*x                                        (2)
*        .  .  .  .  .  .
*     xB[m] = ... + alfa[m]*x
*
*  where xB are basic (auxiliary and structural) variables, m is the
*  number of rows in the problem object.
*
*  On exit the routine stores indices and numerical values of non-zero
*  elements of the resultant column (2) in locations ind[1], ...,
*  ind[len'] and val[1], ..., val[len'], where 0 <= len' <= m is the
*  number of non-zero element in the resultant column returned by the
*  routine. Note that indices (numbers) of basic variables stored in
*  the array ind correspond to original ordinal numbers of variables:
*  indices 1 to m mean auxiliary variables and indices m+1 to m+n mean
*  structural ones.
*
*  RETURNS
*
*  The routine returns len', which is the number of non-zero elements
*  in the resultant column stored in the arrays ind and val.
*
*  BACKGROUND
*
*  The explicitly specified column (1) is transformed in the same way
*  as any other column of the constraint matrix using the formula:
*
*     alfa = inv(B) * a,                                             (3)
*
*  where alfa is the resultant column computed by the routine. */

int glp_transform_col(glp_prob *P, int len, int ind[], double val[])
{     int i, m, t;
      double *a, *alfa;
      if (!glp_bf_exists(P))
         xerror("glp_transform_col: basis factorization does not exist "
            "\n");
      m = glp_get_num_rows(P);
      /* unpack the column to be transformed to the array a */
      a = xcalloc(1+m, sizeof(double));
      for (i = 1; i <= m; i++) a[i] = 0.0;
      if (!(0 <= len && len <= m))
         xerror("glp_transform_col: len = %d; invalid column length\n",
            len);
      for (t = 1; t <= len; t++)
      {  i = ind[t];
         if (!(1 <= i && i <= m))
            xerror("glp_transform_col: ind[%d] = %d; row index out of r"
               "ange\n", t, i);
         if (val[t] == 0.0)
            xerror("glp_transform_col: val[%d] = 0; zero coefficient no"
               "t allowed\n", t);
         if (a[i] != 0.0)
            xerror("glp_transform_col: ind[%d] = %d; duplicate row indi"
               "ces not allowed\n", t, i);
         a[i] = val[t];
      }
      /* solve the system B*a = alfa to compute the vector alfa */
      alfa = a, glp_ftran(P, alfa);
      /* store resultant coefficients */
      len = 0;
      for (i = 1; i <= m; i++)
      {  if (alfa[i] != 0.0)
         {  len++;
            ind[len] = glp_get_bhead(P, i);
            val[len] = alfa[i];
         }
      }
      xfree(a);
      return len;
}

/***********************************************************************
*  NAME
*
*  glp_prim_rtest - perform primal ratio test
*
*  SYNOPSIS
*
*  int glp_prim_rtest(glp_prob *P, int len, const int ind[],
*     const double val[], int dir, double eps);
*
*  DESCRIPTION
*
*  The routine glp_prim_rtest performs the primal ratio test using an
*  explicitly specified column of the simplex table.
*
*  The current basic solution associated with the LP problem object
*  must be primal feasible.
*
*  The explicitly specified column of the simplex table shows how the
*  basic variables xB depend on some non-basic variable x (which is not
*  necessarily presented in the problem object):
*
*     xB[1] = ... + alfa[1] * x + ...
*     xB[2] = ... + alfa[2] * x + ...                                (*)
*         .  .  .  .  .  .  .  .
*     xB[m] = ... + alfa[m] * x + ...
*
*  The column (*) is specifed on entry to the routine using the sparse
*  format. Ordinal numbers of basic variables xB[i] should be placed in
*  locations ind[1], ..., ind[len], where ordinal number 1 to m denote
*  auxiliary variables, and ordinal numbers m+1 to m+n denote structural
*  variables. The corresponding non-zero coefficients alfa[i] should be
*  placed in locations val[1], ..., val[len]. The arrays ind and val are
*  not changed on exit.
*
*  The parameter dir specifies direction in which the variable x changes
*  on entering the basis: +1 means increasing, -1 means decreasing.
*
*  The parameter eps is an absolute tolerance (small positive number)
*  used by the routine to skip small alfa[j] of the row (*).
*
*  The routine determines which basic variable (among specified in
*  ind[1], ..., ind[len]) should leave the basis in order to keep primal
*  feasibility.
*
*  RETURNS
*
*  The routine glp_prim_rtest returns the index piv in the arrays ind
*  and val corresponding to the pivot element chosen, 1 <= piv <= len.
*  If the adjacent basic solution is primal unbounded and therefore the
*  choice cannot be made, the routine returns zero.
*
*  COMMENTS
*
*  If the non-basic variable x is presented in the LP problem object,
*  the column (*) can be computed with the routine glp_eval_tab_col;
*  otherwise it can be computed with the routine glp_transform_col. */

int glp_prim_rtest(glp_prob *P, int len, const int ind[],
      const double val[], int dir, double eps)
{     int k, m, n, piv, t, type, stat;
      double alfa, big, beta, lb, ub, temp, teta;
      if (glp_get_prim_stat(P) != GLP_FEAS)
         xerror("glp_prim_rtest: basic solution is not primal feasible "
            "\n");
      if (!(dir == +1 || dir == -1))
         xerror("glp_prim_rtest: dir = %d; invalid parameter\n", dir);
      if (!(0.0 < eps && eps < 1.0))
         xerror("glp_prim_rtest: eps = %g; invalid parameter\n", eps);
      m = glp_get_num_rows(P);
      n = glp_get_num_cols(P);
      /* initial settings */
      piv = 0, teta = DBL_MAX, big = 0.0;
      /* walk through the entries of the specified column */
      for (t = 1; t <= len; t++)
      {  /* get the ordinal number of basic variable */
         k = ind[t];
         if (!(1 <= k && k <= m+n))
            xerror("glp_prim_rtest: ind[%d] = %d; variable number out o"
               "f range\n", t, k);
         /* determine type, bounds, status and primal value of basic
            variable xB[i] = x[k] in the current basic solution */
         if (k <= m)
         {  type = glp_get_row_type(P, k);
            lb = glp_get_row_lb(P, k);
            ub = glp_get_row_ub(P, k);
            stat = glp_get_row_stat(P, k);
            beta = glp_get_row_prim(P, k);
         }
         else
         {  type = glp_get_col_type(P, k-m);
            lb = glp_get_col_lb(P, k-m);
            ub = glp_get_col_ub(P, k-m);
            stat = glp_get_col_stat(P, k-m);
            beta = glp_get_col_prim(P, k-m);
         }
         if (stat != GLP_BS)
            xerror("glp_prim_rtest: ind[%d] = %d; non-basic variable no"
               "t allowed\n", t, k);
         /* determine influence coefficient at basic variable xB[i]
            in the explicitly specified column and turn to the case of
            increasing the variable x in order to simplify the program
            logic */
         alfa = (dir > 0 ? + val[t] : - val[t]);
         /* analyze main cases */
         if (type == GLP_FR)
         {  /* xB[i] is free variable */
            continue;
         }
         else if (type == GLP_LO)
lo:      {  /* xB[i] has an lower bound */
            if (alfa > - eps) continue;
            temp = (lb - beta) / alfa;
         }
         else if (type == GLP_UP)
up:      {  /* xB[i] has an upper bound */
            if (alfa < + eps) continue;
            temp = (ub - beta) / alfa;
         }
         else if (type == GLP_DB)
         {  /* xB[i] has both lower and upper bounds */
            if (alfa < 0.0) goto lo; else goto up;
         }
         else if (type == GLP_FX)
         {  /* xB[i] is fixed variable */
            if (- eps < alfa && alfa < + eps) continue;
            temp = 0.0;
         }
         else
            xassert(type != type);
         /* if the value of the variable xB[i] violates its lower or
            upper bound (slightly, because the current basis is assumed
            to be primal feasible), temp is negative; we can think this
            happens due to round-off errors and the value is exactly on
            the bound; this allows replacing temp by zero */
         if (temp < 0.0) temp = 0.0;
         /* apply the minimal ratio test */
         if (teta > temp || teta == temp && big < fabs(alfa))
            piv = t, teta = temp, big = fabs(alfa);
      }
      /* return index of the pivot element chosen */
      return piv;
}

/***********************************************************************
*  NAME
*
*  glp_dual_rtest - perform dual ratio test
*
*  SYNOPSIS
*
*  int glp_dual_rtest(glp_prob *P, int len, const int ind[],
*     const double val[], int dir, double eps);
*
*  DESCRIPTION
*
*  The routine glp_dual_rtest performs the dual ratio test using an
*  explicitly specified row of the simplex table.
*
*  The current basic solution associated with the LP problem object
*  must be dual feasible.
*
*  The explicitly specified row of the simplex table is a linear form
*  that shows how some basic variable x (which is not necessarily
*  presented in the problem object) depends on non-basic variables xN:
*
*     x = alfa[1] * xN[1] + alfa[2] * xN[2] + ... + alfa[n] * xN[n]. (*)
*
*  The row (*) is specified on entry to the routine using the sparse
*  format. Ordinal numbers of non-basic variables xN[j] should be placed
*  in locations ind[1], ..., ind[len], where ordinal numbers 1 to m
*  denote auxiliary variables, and ordinal numbers m+1 to m+n denote
*  structural variables. The corresponding non-zero coefficients alfa[j]
*  should be placed in locations val[1], ..., val[len]. The arrays ind
*  and val are not changed on exit.
*
*  The parameter dir specifies direction in which the variable x changes
*  on leaving the basis: +1 means that x goes to its lower bound, and -1
*  means that x goes to its upper bound.
*
*  The parameter eps is an absolute tolerance (small positive number)
*  used by the routine to skip small alfa[j] of the row (*).
*
*  The routine determines which non-basic variable (among specified in
*  ind[1], ..., ind[len]) should enter the basis in order to keep dual
*  feasibility.
*
*  RETURNS
*
*  The routine glp_dual_rtest returns the index piv in the arrays ind
*  and val corresponding to the pivot element chosen, 1 <= piv <= len.
*  If the adjacent basic solution is dual unbounded and therefore the
*  choice cannot be made, the routine returns zero.
*
*  COMMENTS
*
*  If the basic variable x is presented in the LP problem object, the
*  row (*) can be computed with the routine glp_eval_tab_row; otherwise
*  it can be computed with the routine glp_transform_row. */

int glp_dual_rtest(glp_prob *P, int len, const int ind[],
      const double val[], int dir, double eps)
{     int k, m, n, piv, t, stat;
      double alfa, big, cost, obj, temp, teta;
      if (glp_get_dual_stat(P) != GLP_FEAS)
         xerror("glp_dual_rtest: basic solution is not dual feasible\n")
            ;
      if (!(dir == +1 || dir == -1))
         xerror("glp_dual_rtest: dir = %d; invalid parameter\n", dir);
      if (!(0.0 < eps && eps < 1.0))
         xerror("glp_dual_rtest: eps = %g; invalid parameter\n", eps);
      m = glp_get_num_rows(P);
      n = glp_get_num_cols(P);
      /* take into account optimization direction */
      obj = (glp_get_obj_dir(P) == GLP_MIN ? +1.0 : -1.0);
      /* initial settings */
      piv = 0, teta = DBL_MAX, big = 0.0;
      /* walk through the entries of the specified row */
      for (t = 1; t <= len; t++)
      {  /* get ordinal number of non-basic variable */
         k = ind[t];
         if (!(1 <= k && k <= m+n))
            xerror("glp_dual_rtest: ind[%d] = %d; variable number out o"
               "f range\n", t, k);
         /* determine status and reduced cost of non-basic variable
            x[k] = xN[j] in the current basic solution */
         if (k <= m)
         {  stat = glp_get_row_stat(P, k);
            cost = glp_get_row_dual(P, k);
         }
         else
         {  stat = glp_get_col_stat(P, k-m);
            cost = glp_get_col_dual(P, k-m);
         }
         if (stat == GLP_BS)
            xerror("glp_dual_rtest: ind[%d] = %d; basic variable not al"
               "lowed\n", t, k);
         /* determine influence coefficient at non-basic variable xN[j]
            in the explicitly specified row and turn to the case of
            increasing the variable x in order to simplify the program
            logic */
         alfa = (dir > 0 ? + val[t] : - val[t]);
         /* analyze main cases */
         if (stat == GLP_NL)
         {  /* xN[j] is on its lower bound */
            if (alfa < + eps) continue;
            temp = (obj * cost) / alfa;
         }
         else if (stat == GLP_NU)
         {  /* xN[j] is on its upper bound */
            if (alfa > - eps) continue;
            temp = (obj * cost) / alfa;
         }
         else if (stat == GLP_NF)
         {  /* xN[j] is non-basic free variable */
            if (- eps < alfa && alfa < + eps) continue;
            temp = 0.0;
         }
         else if (stat == GLP_NS)
         {  /* xN[j] is non-basic fixed variable */
            continue;
         }
         else
            xassert(stat != stat);
         /* if the reduced cost of the variable xN[j] violates its zero
            bound (slightly, because the current basis is assumed to be
            dual feasible), temp is negative; we can think this happens
            due to round-off errors and the reduced cost is exact zero;
            this allows replacing temp by zero */
         if (temp < 0.0) temp = 0.0;
         /* apply the minimal ratio test */
         if (teta > temp || teta == temp && big < fabs(alfa))
            piv = t, teta = temp, big = fabs(alfa);
      }
      /* return index of the pivot element chosen */
      return piv;
}

/***********************************************************************
*  NAME
*
*  glp_analyze_row - simulate one iteration of dual simplex method
*
*  SYNOPSIS
*
*  int glp_analyze_row(glp_prob *P, int len, const int ind[],
*     const double val[], int type, double rhs, double eps, int *piv,
*     double *x, double *dx, double *y, double *dy, double *dz);
*
*  DESCRIPTION
*
*  Let the current basis be optimal or dual feasible, and there be
*  specified a row (constraint), which is violated by the current basic
*  solution. The routine glp_analyze_row simulates one iteration of the
*  dual simplex method to determine some information on the adjacent
*  basis (see below), where the specified row becomes active constraint
*  (i.e. its auxiliary variable becomes non-basic).
*
*  The current basic solution associated with the problem object passed
*  to the routine must be dual feasible, and its primal components must
*  be defined.
*
*  The row to be analyzed must be previously transformed either with
*  the routine glp_eval_tab_row (if the row is in the problem object)
*  or with the routine glp_transform_row (if the row is external, i.e.
*  not in the problem object). This is needed to express the row only
*  through (auxiliary and structural) variables, which are non-basic in
*  the current basis:
*
*     y = alfa[1] * xN[1] + alfa[2] * xN[2] + ... + alfa[n] * xN[n],
*
*  where y is an auxiliary variable of the row, alfa[j] is an influence
*  coefficient, xN[j] is a non-basic variable.
*
*  The row is passed to the routine in sparse format. Ordinal numbers
*  of non-basic variables are stored in locations ind[1], ..., ind[len],
*  where numbers 1 to m denote auxiliary variables while numbers m+1 to
*  m+n denote structural variables. Corresponding non-zero coefficients
*  alfa[j] are stored in locations val[1], ..., val[len]. The arrays
*  ind and val are ot changed on exit.
*
*  The parameters type and rhs specify the row type and its right-hand
*  side as follows:
*
*     type = GLP_LO: y = sum alfa[j] * xN[j] >= rhs
*
*     type = GLP_UP: y = sum alfa[j] * xN[j] <= rhs
*
*  The parameter eps is an absolute tolerance (small positive number)
*  used by the routine to skip small coefficients alfa[j] on performing
*  the dual ratio test.
*
*  If the operation was successful, the routine stores the following
*  information to corresponding location (if some parameter is NULL,
*  its value is not stored):
*
*  piv   index in the array ind and val, 1 <= piv <= len, determining
*        the non-basic variable, which would enter the adjacent basis;
*
*  x     value of the non-basic variable in the current basis;
*
*  dx    difference between values of the non-basic variable in the
*        adjacent and current bases, dx = x.new - x.old;
*
*  y     value of the row (i.e. of its auxiliary variable) in the
*        current basis;
*
*  dy    difference between values of the row in the adjacent and
*        current bases, dy = y.new - y.old;
*
*  dz    difference between values of the objective function in the
*        adjacent and current bases, dz = z.new - z.old. Note that in
*        case of minimization dz >= 0, and in case of maximization
*        dz <= 0, i.e. in the adjacent basis the objective function
*        always gets worse (degrades). */

int _glp_analyze_row(glp_prob *P, int len, const int ind[],
      const double val[], int type, double rhs, double eps, int *_piv,
      double *_x, double *_dx, double *_y, double *_dy, double *_dz)
{     int t, k, dir, piv, ret = 0;
      double x, dx, y, dy, dz;
      if (P->pbs_stat == GLP_UNDEF)
         xerror("glp_analyze_row: primal basic solution components are "
            "undefined\n");
      if (P->dbs_stat != GLP_FEAS)
         xerror("glp_analyze_row: basic solution is not dual feasible\n"
            );
      /* compute the row value y = sum alfa[j] * xN[j] in the current
         basis */
      if (!(0 <= len && len <= P->n))
         xerror("glp_analyze_row: len = %d; invalid row length\n", len);
      y = 0.0;
      for (t = 1; t <= len; t++)
      {  /* determine value of x[k] = xN[j] in the current basis */
         k = ind[t];
         if (!(1 <= k && k <= P->m+P->n))
            xerror("glp_analyze_row: ind[%d] = %d; row/column index out"
               " of range\n", t, k);
         if (k <= P->m)
         {  /* x[k] is auxiliary variable */
            if (P->row[k]->stat == GLP_BS)
               xerror("glp_analyze_row: ind[%d] = %d; basic auxiliary v"
                  "ariable is not allowed\n", t, k);
            x = P->row[k]->prim;
         }
         else
         {  /* x[k] is structural variable */
            if (P->col[k-P->m]->stat == GLP_BS)
               xerror("glp_analyze_row: ind[%d] = %d; basic structural "
                  "variable is not allowed\n", t, k);
            x = P->col[k-P->m]->prim;
         }
         y += val[t] * x;
      }
      /* check if the row is primal infeasible in the current basis,
         i.e. the constraint is violated at the current point */
      if (type == GLP_LO)
      {  if (y >= rhs)
         {  /* the constraint is not violated */
            ret = 1;
            goto done;
         }
         /* in the adjacent basis y goes to its lower bound */
         dir = +1;
      }
      else if (type == GLP_UP)
      {  if (y <= rhs)
         {  /* the constraint is not violated */
            ret = 1;
            goto done;
         }
         /* in the adjacent basis y goes to its upper bound */
         dir = -1;
      }
      else
         xerror("glp_analyze_row: type = %d; invalid parameter\n",
            type);
      /* compute dy = y.new - y.old */
      dy = rhs - y;
      /* perform dual ratio test to determine which non-basic variable
         should enter the adjacent basis to keep it dual feasible */
      piv = glp_dual_rtest(P, len, ind, val, dir, eps);
      if (piv == 0)
      {  /* no dual feasible adjacent basis exists */
         ret = 2;
         goto done;
      }
      /* non-basic variable x[k] = xN[j] should enter the basis */
      k = ind[piv];
      xassert(1 <= k && k <= P->m+P->n);
      /* determine its value in the current basis */
      if (k <= P->m)
         x = P->row[k]->prim;
      else
         x = P->col[k-P->m]->prim;
      /* compute dx = x.new - x.old = dy / alfa[j] */
      xassert(val[piv] != 0.0);
      dx = dy / val[piv];
      /* compute dz = z.new - z.old = d[j] * dx, where d[j] is reduced
         cost of xN[j] in the current basis */
      if (k <= P->m)
         dz = P->row[k]->dual * dx;
      else
         dz = P->col[k-P->m]->dual * dx;
      /* store the analysis results */
      if (_piv != NULL) *_piv = piv;
      if (_x   != NULL) *_x   = x;
      if (_dx  != NULL) *_dx  = dx;
      if (_y   != NULL) *_y   = y;
      if (_dy  != NULL) *_dy  = dy;
      if (_dz  != NULL) *_dz  = dz;
done: return ret;
}

#if 0
int main(void)
{     /* example program for the routine glp_analyze_row */
      glp_prob *P;
      glp_smcp parm;
      int i, k, len, piv, ret, ind[1+100];
      double rhs, x, dx, y, dy, dz, val[1+100];
      P = glp_create_prob();
      /* read plan.mps (see glpk/examples) */
      ret = glp_read_mps(P, GLP_MPS_DECK, NULL, "plan.mps");
      glp_assert(ret == 0);
      /* and solve it to optimality */
      ret = glp_simplex(P, NULL);
      glp_assert(ret == 0);
      glp_assert(glp_get_status(P) == GLP_OPT);
      /* the optimal objective value is 296.217 */
      /* we would like to know what happens if we would add a new row
         (constraint) to plan.mps:
         .01 * bin1 + .01 * bin2 + .02 * bin4 + .02 * bin5 <= 12 */
      /* first, we specify this new row */
      glp_create_index(P);
      len = 0;
      ind[++len] = glp_find_col(P, "BIN1"), val[len] = .01;
      ind[++len] = glp_find_col(P, "BIN2"), val[len] = .01;
      ind[++len] = glp_find_col(P, "BIN4"), val[len] = .02;
      ind[++len] = glp_find_col(P, "BIN5"), val[len] = .02;
      rhs = 12;
      /* then we can compute value of the row (i.e. of its auxiliary
         variable) in the current basis to see if the constraint is
         violated */
      y = 0.0;
      for (k = 1; k <= len; k++)
         y += val[k] * glp_get_col_prim(P, ind[k]);
      glp_printf("y = %g\n", y);
      /* this prints y = 15.1372, so the constraint is violated, since
         we require that y <= rhs = 12 */
      /* now we transform the row to express it only through non-basic
         (auxiliary and artificial) variables */
      len = glp_transform_row(P, len, ind, val);
      /* finally, we simulate one step of the dual simplex method to
         obtain necessary information for the adjacent basis */
      ret = _glp_analyze_row(P, len, ind, val, GLP_UP, rhs, 1e-9, &piv,
         &x, &dx, &y, &dy, &dz);
      glp_assert(ret == 0);
      glp_printf("k = %d, x = %g; dx = %g; y = %g; dy = %g; dz = %g\n",
         ind[piv], x, dx, y, dy, dz);
      /* this prints dz = 5.64418 and means that in the adjacent basis
         the objective function would be 296.217 + 5.64418 = 301.861 */
      /* now we actually include the row into the problem object; note
         that the arrays ind and val are clobbered, so we need to build
         them once again */
      len = 0;
      ind[++len] = glp_find_col(P, "BIN1"), val[len] = .01;
      ind[++len] = glp_find_col(P, "BIN2"), val[len] = .01;
      ind[++len] = glp_find_col(P, "BIN4"), val[len] = .02;
      ind[++len] = glp_find_col(P, "BIN5"), val[len] = .02;
      rhs = 12;
      i = glp_add_rows(P, 1);
      glp_set_row_bnds(P, i, GLP_UP, 0, rhs);
      glp_set_mat_row(P, i, len, ind, val);
      /* and perform one dual simplex iteration */
      glp_init_smcp(&parm);
      parm.meth = GLP_DUAL;
      parm.it_lim = 1;
      glp_simplex(P, &parm);
      /* the current objective value is 301.861 */
      return 0;
}
#endif

/***********************************************************************
*  NAME
*
*  glp_analyze_bound - analyze active bound of non-basic variable
*
*  SYNOPSIS
*
*  void glp_analyze_bound(glp_prob *P, int k, double *limit1, int *var1,
*     double *limit2, int *var2);
*
*  DESCRIPTION
*
*  The routine glp_analyze_bound analyzes the effect of varying the
*  active bound of specified non-basic variable.
*
*  The non-basic variable is specified by the parameter k, where
*  1 <= k <= m means auxiliary variable of corresponding row while
*  m+1 <= k <= m+n means structural variable (column).
*
*  Note that the current basic solution must be optimal, and the basis
*  factorization must exist.
*
*  Results of the analysis have the following meaning.
*
*  value1 is the minimal value of the active bound, at which the basis
*  still remains primal feasible and thus optimal. -DBL_MAX means that
*  the active bound has no lower limit.
*
*  var1 is the ordinal number of an auxiliary (1 to m) or structural
*  (m+1 to n) basic variable, which reaches its bound first and thereby
*  limits further decreasing the active bound being analyzed.
*  if value1 = -DBL_MAX, var1 is set to 0.
*
*  value2 is the maximal value of the active bound, at which the basis
*  still remains primal feasible and thus optimal. +DBL_MAX means that
*  the active bound has no upper limit.
*
*  var2 is the ordinal number of an auxiliary (1 to m) or structural
*  (m+1 to n) basic variable, which reaches its bound first and thereby
*  limits further increasing the active bound being analyzed.
*  if value2 = +DBL_MAX, var2 is set to 0. */

void glp_analyze_bound(glp_prob *P, int k, double *value1, int *var1,
      double *value2, int *var2)
{     GLPROW *row;
      GLPCOL *col;
      int m, n, stat, kase, p, len, piv, *ind;
      double x, new_x, ll, uu, xx, delta, *val;
#if 0 /* 04/IV-2016 */
      /* sanity checks */
      if (P == NULL || P->magic != GLP_PROB_MAGIC)
         xerror("glp_analyze_bound: P = %p; invalid problem object\n",
            P);
#endif
      m = P->m, n = P->n;
      if (!(P->pbs_stat == GLP_FEAS && P->dbs_stat == GLP_FEAS))
         xerror("glp_analyze_bound: optimal basic solution required\n");
      if (!(m == 0 || P->valid))
         xerror("glp_analyze_bound: basis factorization required\n");
      if (!(1 <= k && k <= m+n))
         xerror("glp_analyze_bound: k = %d; variable number out of rang"
            "e\n", k);
      /* retrieve information about the specified non-basic variable
         x[k] whose active bound is to be analyzed */
      if (k <= m)
      {  row = P->row[k];
         stat = row->stat;
         x = row->prim;
      }
      else
      {  col = P->col[k-m];
         stat = col->stat;
         x = col->prim;
      }
      if (stat == GLP_BS)
         xerror("glp_analyze_bound: k = %d; basic variable not allowed "
            "\n", k);
      /* allocate working arrays */
      ind = xcalloc(1+m, sizeof(int));
      val = xcalloc(1+m, sizeof(double));
      /* compute column of the simplex table corresponding to the
         non-basic variable x[k] */
      len = glp_eval_tab_col(P, k, ind, val);
      xassert(0 <= len && len <= m);
      /* perform analysis */
      for (kase = -1; kase <= +1; kase += 2)
      {  /* kase < 0 means active bound of x[k] is decreasing;
            kase > 0 means active bound of x[k] is increasing */
         /* use the primal ratio test to determine some basic variable
            x[p] which reaches its bound first */
         piv = glp_prim_rtest(P, len, ind, val, kase, 1e-9);
         if (piv == 0)
         {  /* nothing limits changing the active bound of x[k] */
            p = 0;
            new_x = (kase < 0 ? -DBL_MAX : +DBL_MAX);
            goto store;
         }
         /* basic variable x[p] limits changing the active bound of
            x[k]; determine its value in the current basis */
         xassert(1 <= piv && piv <= len);
         p = ind[piv];
         if (p <= m)
         {  row = P->row[p];
            ll = glp_get_row_lb(P, row->i);
            uu = glp_get_row_ub(P, row->i);
            stat = row->stat;
            xx = row->prim;
         }
         else
         {  col = P->col[p-m];
            ll = glp_get_col_lb(P, col->j);
            uu = glp_get_col_ub(P, col->j);
            stat = col->stat;
            xx = col->prim;
         }
         xassert(stat == GLP_BS);
         /* determine delta x[p] = bound of x[p] - value of x[p] */
         if (kase < 0 && val[piv] > 0.0 ||
             kase > 0 && val[piv] < 0.0)
         {  /* delta x[p] < 0, so x[p] goes toward its lower bound */
            xassert(ll != -DBL_MAX);
            delta = ll - xx;
         }
         else
         {  /* delta x[p] > 0, so x[p] goes toward its upper bound */
            xassert(uu != +DBL_MAX);
            delta = uu - xx;
         }
         /* delta x[p] = alfa[p,k] * delta x[k], so new x[k] = x[k] +
            delta x[k] = x[k] + delta x[p] / alfa[p,k] is the value of
            x[k] in the adjacent basis */
         xassert(val[piv] != 0.0);
         new_x = x + delta / val[piv];
store:   /* store analysis results */
         if (kase < 0)
         {  if (value1 != NULL) *value1 = new_x;
            if (var1 != NULL) *var1 = p;
         }
         else
         {  if (value2 != NULL) *value2 = new_x;
            if (var2 != NULL) *var2 = p;
         }
      }
      /* free working arrays */
      xfree(ind);
      xfree(val);
      return;
}

/***********************************************************************
*  NAME
*
*  glp_analyze_coef - analyze objective coefficient at basic variable
*
*  SYNOPSIS
*
*  void glp_analyze_coef(glp_prob *P, int k, double *coef1, int *var1,
*     double *value1, double *coef2, int *var2, double *value2);
*
*  DESCRIPTION
*
*  The routine glp_analyze_coef analyzes the effect of varying the
*  objective coefficient at specified basic variable.
*
*  The basic variable is specified by the parameter k, where
*  1 <= k <= m means auxiliary variable of corresponding row while
*  m+1 <= k <= m+n means structural variable (column).
*
*  Note that the current basic solution must be optimal, and the basis
*  factorization must exist.
*
*  Results of the analysis have the following meaning.
*
*  coef1 is the minimal value of the objective coefficient, at which
*  the basis still remains dual feasible and thus optimal. -DBL_MAX
*  means that the objective coefficient has no lower limit.
*
*  var1 is the ordinal number of an auxiliary (1 to m) or structural
*  (m+1 to n) non-basic variable, whose reduced cost reaches its zero
*  bound first and thereby limits further decreasing the objective
*  coefficient being analyzed. If coef1 = -DBL_MAX, var1 is set to 0.
*
*  value1 is value of the basic variable being analyzed in an adjacent
*  basis, which is defined as follows. Let the objective coefficient
*  reaches its minimal value (coef1) and continues decreasing. Then the
*  reduced cost of the limiting non-basic variable (var1) becomes dual
*  infeasible and the current basis becomes non-optimal that forces the
*  limiting non-basic variable to enter the basis replacing there some
*  basic variable that leaves the basis to keep primal feasibility.
*  Should note that on determining the adjacent basis current bounds
*  of the basic variable being analyzed are ignored as if it were free
*  (unbounded) variable, so it cannot leave the basis. It may happen
*  that no dual feasible adjacent basis exists, in which case value1 is
*  set to -DBL_MAX or +DBL_MAX.
*
*  coef2 is the maximal value of the objective coefficient, at which
*  the basis still remains dual feasible and thus optimal. +DBL_MAX
*  means that the objective coefficient has no upper limit.
*
*  var2 is the ordinal number of an auxiliary (1 to m) or structural
*  (m+1 to n) non-basic variable, whose reduced cost reaches its zero
*  bound first and thereby limits further increasing the objective
*  coefficient being analyzed. If coef2 = +DBL_MAX, var2 is set to 0.
*
*  value2 is value of the basic variable being analyzed in an adjacent
*  basis, which is defined exactly in the same way as value1 above with
*  exception that now the objective coefficient is increasing. */

void glp_analyze_coef(glp_prob *P, int k, double *coef1, int *var1,
      double *value1, double *coef2, int *var2, double *value2)
{     GLPROW *row; GLPCOL *col;
      int m, n, type, stat, kase, p, q, dir, clen, cpiv, rlen, rpiv,
         *cind, *rind;
      double lb, ub, coef, x, lim_coef, new_x, d, delta, ll, uu, xx,
         *rval, *cval;
#if 0 /* 04/IV-2016 */
      /* sanity checks */
      if (P == NULL || P->magic != GLP_PROB_MAGIC)
         xerror("glp_analyze_coef: P = %p; invalid problem object\n",
            P);
#endif
      m = P->m, n = P->n;
      if (!(P->pbs_stat == GLP_FEAS && P->dbs_stat == GLP_FEAS))
         xerror("glp_analyze_coef: optimal basic solution required\n");
      if (!(m == 0 || P->valid))
         xerror("glp_analyze_coef: basis factorization required\n");
      if (!(1 <= k && k <= m+n))
         xerror("glp_analyze_coef: k = %d; variable number out of range"
            "\n", k);
      /* retrieve information about the specified basic variable x[k]
         whose objective coefficient c[k] is to be analyzed */
      if (k <= m)
      {  row = P->row[k];
         type = row->type;
         lb = row->lb;
         ub = row->ub;
         coef = 0.0;
         stat = row->stat;
         x = row->prim;
      }
      else
      {  col = P->col[k-m];
         type = col->type;
         lb = col->lb;
         ub = col->ub;
         coef = col->coef;
         stat = col->stat;
         x = col->prim;
      }
      if (stat != GLP_BS)
         xerror("glp_analyze_coef: k = %d; non-basic variable not allow"
            "ed\n", k);
      /* allocate working arrays */
      cind = xcalloc(1+m, sizeof(int));
      cval = xcalloc(1+m, sizeof(double));
      rind = xcalloc(1+n, sizeof(int));
      rval = xcalloc(1+n, sizeof(double));
      /* compute row of the simplex table corresponding to the basic
         variable x[k] */
      rlen = glp_eval_tab_row(P, k, rind, rval);
      xassert(0 <= rlen && rlen <= n);
      /* perform analysis */
      for (kase = -1; kase <= +1; kase += 2)
      {  /* kase < 0 means objective coefficient c[k] is decreasing;
            kase > 0 means objective coefficient c[k] is increasing */
         /* note that decreasing c[k] is equivalent to increasing dual
            variable lambda[k] and vice versa; we need to correctly set
            the dir flag as required by the routine glp_dual_rtest */
         if (P->dir == GLP_MIN)
            dir = - kase;
         else if (P->dir == GLP_MAX)
            dir = + kase;
         else
            xassert(P != P);
         /* use the dual ratio test to determine non-basic variable
            x[q] whose reduced cost d[q] reaches zero bound first */
         rpiv = glp_dual_rtest(P, rlen, rind, rval, dir, 1e-9);
         if (rpiv == 0)
         {  /* nothing limits changing c[k] */
            lim_coef = (kase < 0 ? -DBL_MAX : +DBL_MAX);
            q = 0;
            /* x[k] keeps its current value */
            new_x = x;
            goto store;
         }
         /* non-basic variable x[q] limits changing coefficient c[k];
            determine its status and reduced cost d[k] in the current
            basis */
         xassert(1 <= rpiv && rpiv <= rlen);
         q = rind[rpiv];
         xassert(1 <= q && q <= m+n);
         if (q <= m)
         {  row = P->row[q];
            stat = row->stat;
            d = row->dual;
         }
         else
         {  col = P->col[q-m];
            stat = col->stat;
            d = col->dual;
         }
         /* note that delta d[q] = new d[q] - d[q] = - d[q], because
            new d[q] = 0; delta d[q] = alfa[k,q] * delta c[k], so
            delta c[k] = delta d[q] / alfa[k,q] = - d[q] / alfa[k,q] */
         xassert(rval[rpiv] != 0.0);
         delta = - d / rval[rpiv];
         /* compute new c[k] = c[k] + delta c[k], which is the limiting
            value of the objective coefficient c[k] */
         lim_coef = coef + delta;
         /* let c[k] continue decreasing/increasing that makes d[q]
            dual infeasible and forces x[q] to enter the basis;
            to perform the primal ratio test we need to know in which
            direction x[q] changes on entering the basis; we determine
            that analyzing the sign of delta d[q] (see above), since
            d[q] may be close to zero having wrong sign */
         /* let, for simplicity, the problem is minimization */
         if (kase < 0 && rval[rpiv] > 0.0 ||
             kase > 0 && rval[rpiv] < 0.0)
         {  /* delta d[q] < 0, so d[q] being non-negative will become
               negative, so x[q] will increase */
            dir = +1;
         }
         else
         {  /* delta d[q] > 0, so d[q] being non-positive will become
               positive, so x[q] will decrease */
            dir = -1;
         }
         /* if the problem is maximization, correct the direction */
         if (P->dir == GLP_MAX) dir = - dir;
         /* check that we didn't make a silly mistake */
         if (dir > 0)
            xassert(stat == GLP_NL || stat == GLP_NF);
         else
            xassert(stat == GLP_NU || stat == GLP_NF);
         /* compute column of the simplex table corresponding to the
            non-basic variable x[q] */
         clen = glp_eval_tab_col(P, q, cind, cval);
         /* make x[k] temporarily free (unbounded) */
         if (k <= m)
         {  row = P->row[k];
            row->type = GLP_FR;
            row->lb = row->ub = 0.0;
         }
         else
         {  col = P->col[k-m];
            col->type = GLP_FR;
            col->lb = col->ub = 0.0;
         }
         /* use the primal ratio test to determine some basic variable
            which leaves the basis */
         cpiv = glp_prim_rtest(P, clen, cind, cval, dir, 1e-9);
         /* restore original bounds of the basic variable x[k] */
         if (k <= m)
         {  row = P->row[k];
            row->type = type;
            row->lb = lb, row->ub = ub;
         }
         else
         {  col = P->col[k-m];
            col->type = type;
            col->lb = lb, col->ub = ub;
         }
         if (cpiv == 0)
         {  /* non-basic variable x[q] can change unlimitedly */
            if (dir < 0 && rval[rpiv] > 0.0 ||
                dir > 0 && rval[rpiv] < 0.0)
            {  /* delta x[k] = alfa[k,q] * delta x[q] < 0 */
               new_x = -DBL_MAX;
            }
            else
            {  /* delta x[k] = alfa[k,q] * delta x[q] > 0 */
               new_x = +DBL_MAX;
            }
            goto store;
         }
         /* some basic variable x[p] limits changing non-basic variable
            x[q] in the adjacent basis */
         xassert(1 <= cpiv && cpiv <= clen);
         p = cind[cpiv];
         xassert(1 <= p && p <= m+n);
         xassert(p != k);
         if (p <= m)
         {  row = P->row[p];
            xassert(row->stat == GLP_BS);
            ll = glp_get_row_lb(P, row->i);
            uu = glp_get_row_ub(P, row->i);
            xx = row->prim;
         }
         else
         {  col = P->col[p-m];
            xassert(col->stat == GLP_BS);
            ll = glp_get_col_lb(P, col->j);
            uu = glp_get_col_ub(P, col->j);
            xx = col->prim;
         }
         /* determine delta x[p] = new x[p] - x[p] */
         if (dir < 0 && cval[cpiv] > 0.0 ||
             dir > 0 && cval[cpiv] < 0.0)
         {  /* delta x[p] < 0, so x[p] goes toward its lower bound */
            xassert(ll != -DBL_MAX);
            delta = ll - xx;
         }
         else
         {  /* delta x[p] > 0, so x[p] goes toward its upper bound */
            xassert(uu != +DBL_MAX);
            delta = uu - xx;
         }
         /* compute new x[k] = x[k] + alfa[k,q] * delta x[q], where
            delta x[q] = delta x[p] / alfa[p,q] */
         xassert(cval[cpiv] != 0.0);
         new_x = x + (rval[rpiv] / cval[cpiv]) * delta;
store:   /* store analysis results */
         if (kase < 0)
         {  if (coef1 != NULL) *coef1 = lim_coef;
            if (var1 != NULL) *var1 = q;
            if (value1 != NULL) *value1 = new_x;
         }
         else
         {  if (coef2 != NULL) *coef2 = lim_coef;
            if (var2 != NULL) *var2 = q;
            if (value2 != NULL) *value2 = new_x;
         }
      }
      /* free working arrays */
      xfree(cind);
      xfree(cval);
      xfree(rind);
      xfree(rval);
      return;
}

/* eof */
