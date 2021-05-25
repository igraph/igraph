/* glpspx02.c (dual simplex method) */

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
#pragma clang diagnostic ignored "-Wcomment"
#pragma clang diagnostic ignored "-Wsign-conversion"
#pragma clang diagnostic ignored "-Wsometimes-uninitialized"
#pragma clang diagnostic ignored "-Wlogical-op-parentheses"
#endif

#include "glpspx.h"

#define GLP_DEBUG 1

#if 0
#define GLP_LONG_STEP 1
#endif

struct csa
{     /* common storage area */
      /*--------------------------------------------------------------*/
      /* LP data */
      int m;
      /* number of rows (auxiliary variables), m > 0 */
      int n;
      /* number of columns (structural variables), n > 0 */
      char *type; /* char type[1+m+n]; */
      /* type[0] is not used;
         type[k], 1 <= k <= m+n, is the type of variable x[k]:
         GLP_FR - free variable
         GLP_LO - variable with lower bound
         GLP_UP - variable with upper bound
         GLP_DB - double-bounded variable
         GLP_FX - fixed variable */
      double *lb; /* double lb[1+m+n]; */
      /* lb[0] is not used;
         lb[k], 1 <= k <= m+n, is an lower bound of variable x[k];
         if x[k] has no lower bound, lb[k] is zero */
      double *ub; /* double ub[1+m+n]; */
      /* ub[0] is not used;
         ub[k], 1 <= k <= m+n, is an upper bound of variable x[k];
         if x[k] has no upper bound, ub[k] is zero;
         if x[k] is of fixed type, ub[k] is the same as lb[k] */
      double *coef; /* double coef[1+m+n]; */
      /* coef[0] is not used;
         coef[k], 1 <= k <= m+n, is an objective coefficient at
         variable x[k] */
      /*--------------------------------------------------------------*/
      /* original bounds of variables */
      char *orig_type; /* char orig_type[1+m+n]; */
      double *orig_lb; /* double orig_lb[1+m+n]; */
      double *orig_ub; /* double orig_ub[1+m+n]; */
      /*--------------------------------------------------------------*/
      /* original objective function */
      double *obj; /* double obj[1+n]; */
      /* obj[0] is a constant term of the original objective function;
         obj[j], 1 <= j <= n, is an original objective coefficient at
         structural variable x[m+j] */
      double zeta;
      /* factor used to scale original objective coefficients; its
         sign defines original optimization direction: zeta > 0 means
         minimization, zeta < 0 means maximization */
      /*--------------------------------------------------------------*/
      /* constraint matrix A; it has m rows and n columns and is stored
         by columns */
      int *A_ptr; /* int A_ptr[1+n+1]; */
      /* A_ptr[0] is not used;
         A_ptr[j], 1 <= j <= n, is starting position of j-th column in
         arrays A_ind and A_val; note that A_ptr[1] is always 1;
         A_ptr[n+1] indicates the position after the last element in
         arrays A_ind and A_val */
      int *A_ind; /* int A_ind[A_ptr[n+1]]; */
      /* row indices */
      double *A_val; /* double A_val[A_ptr[n+1]]; */
      /* non-zero element values */
#if 1 /* 06/IV-2009 */
      /* constraint matrix A stored by rows */
      int *AT_ptr; /* int AT_ptr[1+m+1];
      /* AT_ptr[0] is not used;
         AT_ptr[i], 1 <= i <= m, is starting position of i-th row in
         arrays AT_ind and AT_val; note that AT_ptr[1] is always 1;
         AT_ptr[m+1] indicates the position after the last element in
         arrays AT_ind and AT_val */
      int *AT_ind; /* int AT_ind[AT_ptr[m+1]]; */
      /* column indices */
      double *AT_val; /* double AT_val[AT_ptr[m+1]]; */
      /* non-zero element values */
#endif
      /*--------------------------------------------------------------*/
      /* basis header */
      int *head; /* int head[1+m+n]; */
      /* head[0] is not used;
         head[i], 1 <= i <= m, is the ordinal number of basic variable
         xB[i]; head[i] = k means that xB[i] = x[k] and i-th column of
         matrix B is k-th column of matrix (I|-A);
         head[m+j], 1 <= j <= n, is the ordinal number of non-basic
         variable xN[j]; head[m+j] = k means that xN[j] = x[k] and j-th
         column of matrix N is k-th column of matrix (I|-A) */
#if 1 /* 06/IV-2009 */
      int *bind; /* int bind[1+m+n]; */
      /* bind[0] is not used;
         bind[k], 1 <= k <= m+n, is the position of k-th column of the
         matrix (I|-A) in the matrix (B|N); that is, bind[k] = k' means
         that head[k'] = k */
#endif
      char *stat; /* char stat[1+n]; */
      /* stat[0] is not used;
         stat[j], 1 <= j <= n, is the status of non-basic variable
         xN[j], which defines its active bound:
         GLP_NL - lower bound is active
         GLP_NU - upper bound is active
         GLP_NF - free variable
         GLP_NS - fixed variable */
      /*--------------------------------------------------------------*/
      /* matrix B is the basis matrix; it is composed from columns of
         the augmented constraint matrix (I|-A) corresponding to basic
         variables and stored in a factorized (invertable) form */
      int valid;
      /* factorization is valid only if this flag is set */
      BFD *bfd; /* BFD bfd[1:m,1:m]; */
      /* factorized (invertable) form of the basis matrix */
#if 0 /* 06/IV-2009 */
      /*--------------------------------------------------------------*/
      /* matrix N is a matrix composed from columns of the augmented
         constraint matrix (I|-A) corresponding to non-basic variables
         except fixed ones; it is stored by rows and changes every time
         the basis changes */
      int *N_ptr; /* int N_ptr[1+m+1]; */
      /* N_ptr[0] is not used;
         N_ptr[i], 1 <= i <= m, is starting position of i-th row in
         arrays N_ind and N_val; note that N_ptr[1] is always 1;
         N_ptr[m+1] indicates the position after the last element in
         arrays N_ind and N_val */
      int *N_len; /* int N_len[1+m]; */
      /* N_len[0] is not used;
         N_len[i], 1 <= i <= m, is length of i-th row (0 to n) */
      int *N_ind; /* int N_ind[N_ptr[m+1]]; */
      /* column indices */
      double *N_val; /* double N_val[N_ptr[m+1]]; */
      /* non-zero element values */
#endif
      /*--------------------------------------------------------------*/
      /* working parameters */
      int phase;
      /* search phase:
         0 - not determined yet
         1 - search for dual feasible solution
         2 - search for optimal solution */
      glp_long tm_beg;
      /* time value at the beginning of the search */
      int it_beg;
      /* simplex iteration count at the beginning of the search */
      int it_cnt;
      /* simplex iteration count; it increases by one every time the
         basis changes */
      int it_dpy;
      /* simplex iteration count at the most recent display output */
      /*--------------------------------------------------------------*/
      /* basic solution components */
      double *bbar; /* double bbar[1+m]; */
      /* bbar[0] is not used on phase I; on phase II it is the current
         value of the original objective function;
         bbar[i], 1 <= i <= m, is primal value of basic variable xB[i]
         (if xB[i] is free, its primal value is not updated) */
      double *cbar; /* double cbar[1+n]; */
      /* cbar[0] is not used;
         cbar[j], 1 <= j <= n, is reduced cost of non-basic variable
         xN[j] (if xN[j] is fixed, its reduced cost is not updated) */
      /*--------------------------------------------------------------*/
      /* the following pricing technique options may be used:
         GLP_PT_STD - standard ("textbook") pricing;
         GLP_PT_PSE - projected steepest edge;
         GLP_PT_DVX - Devex pricing (not implemented yet);
         in case of GLP_PT_STD the reference space is not used, and all
         steepest edge coefficients are set to 1 */
      int refct;
      /* this count is set to an initial value when the reference space
         is defined and decreases by one every time the basis changes;
         once this count reaches zero, the reference space is redefined
         again */
      char *refsp; /* char refsp[1+m+n]; */
      /* refsp[0] is not used;
         refsp[k], 1 <= k <= m+n, is the flag which means that variable
         x[k] belongs to the current reference space */
      double *gamma; /* double gamma[1+m]; */
      /* gamma[0] is not used;
         gamma[i], 1 <= i <= n, is the steepest edge coefficient for
         basic variable xB[i]; if xB[i] is free, gamma[i] is not used
         and just set to 1 */
      /*--------------------------------------------------------------*/
      /* basic variable xB[p] chosen to leave the basis */
      int p;
      /* index of the basic variable xB[p] chosen, 1 <= p <= m;
         if the set of eligible basic variables is empty (i.e. if the
         current basic solution is primal feasible within a tolerance)
         and thus no variable has been chosen, p is set to 0 */
      double delta;
      /* change of xB[p] in the adjacent basis;
         delta > 0 means that xB[p] violates its lower bound and will
         increase to achieve it in the adjacent basis;
         delta < 0 means that xB[p] violates its upper bound and will
         decrease to achieve it in the adjacent basis */
      /*--------------------------------------------------------------*/
      /* pivot row of the simplex table corresponding to basic variable
         xB[p] chosen is the following vector:
            T' * e[p] = - N' * inv(B') * e[p] = - N' * rho,
         where B' is a matrix transposed to the current basis matrix,
         N' is a matrix, whose rows are columns of the matrix (I|-A)
         corresponding to non-basic non-fixed variables */
      int trow_nnz;
      /* number of non-zero components, 0 <= nnz <= n */
      int *trow_ind; /* int trow_ind[1+n]; */
      /* trow_ind[0] is not used;
         trow_ind[t], 1 <= t <= nnz, is an index of non-zero component,
         i.e. trow_ind[t] = j means that trow_vec[j] != 0 */
      double *trow_vec; /* int trow_vec[1+n]; */
      /* trow_vec[0] is not used;
         trow_vec[j], 1 <= j <= n, is a numeric value of j-th component
         of the row */
      double trow_max;
      /* infinity (maximum) norm of the row (max |trow_vec[j]|) */
      int trow_num;
      /* number of significant non-zero components, which means that:
         |trow_vec[j]| >= eps for j in trow_ind[1,...,num],
         |tcol_vec[j]| <  eps for j in trow_ind[num+1,...,nnz],
         where eps is a pivot tolerance */
      /*--------------------------------------------------------------*/
#ifdef GLP_LONG_STEP /* 07/IV-2009 */
      int nbps;
      /* number of breakpoints, 0 <= nbps <= n */
      struct bkpt
      {     int j;
            /* index of non-basic variable xN[j], 1 <= j <= n */
            double t;
            /* value of dual ray parameter at breakpoint, t >= 0 */
            double dz;
            /* dz = zeta(t = t[k]) - zeta(t = 0) */
      } *bkpt; /* struct bkpt bkpt[1+n]; */
      /* bkpt[0] is not used;
         bkpt[k], 1 <= k <= nbps, is k-th breakpoint of the dual
         objective */
#endif
      /*--------------------------------------------------------------*/
      /* non-basic variable xN[q] chosen to enter the basis */
      int q;
      /* index of the non-basic variable xN[q] chosen, 1 <= q <= n;
         if no variable has been chosen, q is set to 0 */
      double new_dq;
      /* reduced cost of xN[q] in the adjacent basis (it is the change
         of lambdaB[p]) */
      /*--------------------------------------------------------------*/
      /* pivot column of the simplex table corresponding to non-basic
         variable xN[q] chosen is the following vector:
            T * e[q] = - inv(B) * N * e[q] = - inv(B) * N[q],
         where B is the current basis matrix, N[q] is a column of the
         matrix (I|-A) corresponding to xN[q] */
      int tcol_nnz;
      /* number of non-zero components, 0 <= nnz <= m */
      int *tcol_ind; /* int tcol_ind[1+m]; */
      /* tcol_ind[0] is not used;
         tcol_ind[t], 1 <= t <= nnz, is an index of non-zero component,
         i.e. tcol_ind[t] = i means that tcol_vec[i] != 0 */
      double *tcol_vec; /* double tcol_vec[1+m]; */
      /* tcol_vec[0] is not used;
         tcol_vec[i], 1 <= i <= m, is a numeric value of i-th component
         of the column */
      /*--------------------------------------------------------------*/
      /* working arrays */
      double *work1; /* double work1[1+m]; */
      double *work2; /* double work2[1+m]; */
      double *work3; /* double work3[1+m]; */
      double *work4; /* double work4[1+m]; */
};

static const double kappa = 0.10;

/***********************************************************************
*  alloc_csa - allocate common storage area
*
*  This routine allocates all arrays in the common storage area (CSA)
*  and returns a pointer to the CSA. */

static struct csa *alloc_csa(glp_prob *lp)
{     struct csa *csa;
      int m = lp->m;
      int n = lp->n;
      int nnz = lp->nnz;
      csa = xmalloc(sizeof(struct csa));
      xassert(m > 0 && n > 0);
      csa->m = m;
      csa->n = n;
      csa->type = xcalloc(1+m+n, sizeof(char));
      csa->lb = xcalloc(1+m+n, sizeof(double));
      csa->ub = xcalloc(1+m+n, sizeof(double));
      csa->coef = xcalloc(1+m+n, sizeof(double));
      csa->orig_type = xcalloc(1+m+n, sizeof(char));
      csa->orig_lb = xcalloc(1+m+n, sizeof(double));
      csa->orig_ub = xcalloc(1+m+n, sizeof(double));
      csa->obj = xcalloc(1+n, sizeof(double));
      csa->A_ptr = xcalloc(1+n+1, sizeof(int));
      csa->A_ind = xcalloc(1+nnz, sizeof(int));
      csa->A_val = xcalloc(1+nnz, sizeof(double));
#if 1 /* 06/IV-2009 */
      csa->AT_ptr = xcalloc(1+m+1, sizeof(int));
      csa->AT_ind = xcalloc(1+nnz, sizeof(int));
      csa->AT_val = xcalloc(1+nnz, sizeof(double));
#endif
      csa->head = xcalloc(1+m+n, sizeof(int));
#if 1 /* 06/IV-2009 */
      csa->bind = xcalloc(1+m+n, sizeof(int));
#endif
      csa->stat = xcalloc(1+n, sizeof(char));
#if 0 /* 06/IV-2009 */
      csa->N_ptr = xcalloc(1+m+1, sizeof(int));
      csa->N_len = xcalloc(1+m, sizeof(int));
      csa->N_ind = NULL; /* will be allocated later */
      csa->N_val = NULL; /* will be allocated later */
#endif
      csa->bbar = xcalloc(1+m, sizeof(double));
      csa->cbar = xcalloc(1+n, sizeof(double));
      csa->refsp = xcalloc(1+m+n, sizeof(char));
      csa->gamma = xcalloc(1+m, sizeof(double));
      csa->trow_ind = xcalloc(1+n, sizeof(int));
      csa->trow_vec = xcalloc(1+n, sizeof(double));
#ifdef GLP_LONG_STEP /* 07/IV-2009 */
      csa->bkpt = xcalloc(1+n, sizeof(struct bkpt));
#endif
      csa->tcol_ind = xcalloc(1+m, sizeof(int));
      csa->tcol_vec = xcalloc(1+m, sizeof(double));
      csa->work1 = xcalloc(1+m, sizeof(double));
      csa->work2 = xcalloc(1+m, sizeof(double));
      csa->work3 = xcalloc(1+m, sizeof(double));
      csa->work4 = xcalloc(1+m, sizeof(double));
      return csa;
}

/***********************************************************************
*  init_csa - initialize common storage area
*
*  This routine initializes all data structures in the common storage
*  area (CSA). */

static void init_csa(struct csa *csa, glp_prob *lp)
{     int m = csa->m;
      int n = csa->n;
      char *type = csa->type;
      double *lb = csa->lb;
      double *ub = csa->ub;
      double *coef = csa->coef;
      char *orig_type = csa->orig_type;
      double *orig_lb = csa->orig_lb;
      double *orig_ub = csa->orig_ub;
      double *obj = csa->obj;
      int *A_ptr = csa->A_ptr;
      int *A_ind = csa->A_ind;
      double *A_val = csa->A_val;
#if 1 /* 06/IV-2009 */
      int *AT_ptr = csa->AT_ptr;
      int *AT_ind = csa->AT_ind;
      double *AT_val = csa->AT_val;
#endif
      int *head = csa->head;
#if 1 /* 06/IV-2009 */
      int *bind = csa->bind;
#endif
      char *stat = csa->stat;
      char *refsp = csa->refsp;
      double *gamma = csa->gamma;
      int i, j, k, loc;
      double cmax;
      /* auxiliary variables */
      for (i = 1; i <= m; i++)
      {  GLPROW *row = lp->row[i];
         type[i] = (char)row->type;
         lb[i] = row->lb * row->rii;
         ub[i] = row->ub * row->rii;
         coef[i] = 0.0;
      }
      /* structural variables */
      for (j = 1; j <= n; j++)
      {  GLPCOL *col = lp->col[j];
         type[m+j] = (char)col->type;
         lb[m+j] = col->lb / col->sjj;
         ub[m+j] = col->ub / col->sjj;
         coef[m+j] = col->coef * col->sjj;
      }
      /* original bounds of variables */
      memcpy(&orig_type[1], &type[1], (m+n) * sizeof(char));
      memcpy(&orig_lb[1], &lb[1], (m+n) * sizeof(double));
      memcpy(&orig_ub[1], &ub[1], (m+n) * sizeof(double));
      /* original objective function */
      obj[0] = lp->c0;
      memcpy(&obj[1], &coef[m+1], n * sizeof(double));
      /* factor used to scale original objective coefficients */
      cmax = 0.0;
      for (j = 1; j <= n; j++)
         if (cmax < fabs(obj[j])) cmax = fabs(obj[j]);
      if (cmax == 0.0) cmax = 1.0;
      switch (lp->dir)
      {  case GLP_MIN:
            csa->zeta = + 1.0 / cmax;
            break;
         case GLP_MAX:
            csa->zeta = - 1.0 / cmax;
            break;
         default:
            xassert(lp != lp);
      }
#if 1
      if (fabs(csa->zeta) < 1.0) csa->zeta *= 1000.0;
#endif
      /* scale working objective coefficients */
      for (j = 1; j <= n; j++) coef[m+j] *= csa->zeta;
      /* matrix A (by columns) */
      loc = 1;
      for (j = 1; j <= n; j++)
      {  GLPAIJ *aij;
         A_ptr[j] = loc;
         for (aij = lp->col[j]->ptr; aij != NULL; aij = aij->c_next)
         {  A_ind[loc] = aij->row->i;
            A_val[loc] = aij->row->rii * aij->val * aij->col->sjj;
            loc++;
         }
      }
      A_ptr[n+1] = loc;
      xassert(loc-1 == lp->nnz);
#if 1 /* 06/IV-2009 */
      /* matrix A (by rows) */
      loc = 1;
      for (i = 1; i <= m; i++)
      {  GLPAIJ *aij;
         AT_ptr[i] = loc;
         for (aij = lp->row[i]->ptr; aij != NULL; aij = aij->r_next)
         {  AT_ind[loc] = aij->col->j;
            AT_val[loc] = aij->row->rii * aij->val * aij->col->sjj;
            loc++;
         }
      }
      AT_ptr[m+1] = loc;
      xassert(loc-1 == lp->nnz);
#endif
      /* basis header */
      xassert(lp->valid);
      memcpy(&head[1], &lp->head[1], m * sizeof(int));
      k = 0;
      for (i = 1; i <= m; i++)
      {  GLPROW *row = lp->row[i];
         if (row->stat != GLP_BS)
         {  k++;
            xassert(k <= n);
            head[m+k] = i;
            stat[k] = (char)row->stat;
         }
      }
      for (j = 1; j <= n; j++)
      {  GLPCOL *col = lp->col[j];
         if (col->stat != GLP_BS)
         {  k++;
            xassert(k <= n);
            head[m+k] = m + j;
            stat[k] = (char)col->stat;
         }
      }
      xassert(k == n);
#if 1 /* 06/IV-2009 */
      for (k = 1; k <= m+n; k++)
         bind[head[k]] = k;
#endif
      /* factorization of matrix B */
      csa->valid = 1, lp->valid = 0;
      csa->bfd = lp->bfd, lp->bfd = NULL;
#if 0 /* 06/IV-2009 */
      /* matrix N (by rows) */
      alloc_N(csa);
      build_N(csa);
#endif
      /* working parameters */
      csa->phase = 0;
      csa->tm_beg = xtime();
      csa->it_beg = csa->it_cnt = lp->it_cnt;
      csa->it_dpy = -1;
      /* reference space and steepest edge coefficients */
      csa->refct = 0;
      memset(&refsp[1], 0, (m+n) * sizeof(char));
      for (i = 1; i <= m; i++) gamma[i] = 1.0;
      return;
}

#if 1 /* copied from primal */
/***********************************************************************
*  invert_B - compute factorization of the basis matrix
*
*  This routine computes factorization of the current basis matrix B.
*
*  If the operation is successful, the routine returns zero, otherwise
*  non-zero. */

static int inv_col(void *info, int i, int ind[], double val[])
{     /* this auxiliary routine returns row indices and numeric values
         of non-zero elements of i-th column of the basis matrix */
      struct csa *csa = info;
      int m = csa->m;
#ifdef GLP_DEBUG
      int n = csa->n;
#endif
      int *A_ptr = csa->A_ptr;
      int *A_ind = csa->A_ind;
      double *A_val = csa->A_val;
      int *head = csa->head;
      int k, len, ptr, t;
#ifdef GLP_DEBUG
      xassert(1 <= i && i <= m);
#endif
      k = head[i]; /* B[i] is k-th column of (I|-A) */
#ifdef GLP_DEBUG
      xassert(1 <= k && k <= m+n);
#endif
      if (k <= m)
      {  /* B[i] is k-th column of submatrix I */
         len = 1;
         ind[1] = k;
         val[1] = 1.0;
      }
      else
      {  /* B[i] is (k-m)-th column of submatrix (-A) */
         ptr = A_ptr[k-m];
         len = A_ptr[k-m+1] - ptr;
         memcpy(&ind[1], &A_ind[ptr], len * sizeof(int));
         memcpy(&val[1], &A_val[ptr], len * sizeof(double));
         for (t = 1; t <= len; t++) val[t] = - val[t];
      }
      return len;
}

static int invert_B(struct csa *csa)
{     int ret;
      ret = bfd_factorize(csa->bfd, csa->m, NULL, inv_col, csa);
      csa->valid = (ret == 0);
      return ret;
}
#endif

#if 1 /* copied from primal */
/***********************************************************************
*  update_B - update factorization of the basis matrix
*
*  This routine replaces i-th column of the basis matrix B by k-th
*  column of the augmented constraint matrix (I|-A) and then updates
*  the factorization of B.
*
*  If the factorization has been successfully updated, the routine
*  returns zero, otherwise non-zero. */

static int update_B(struct csa *csa, int i, int k)
{     int m = csa->m;
#ifdef GLP_DEBUG
      int n = csa->n;
#endif
      int ret;
#ifdef GLP_DEBUG
      xassert(1 <= i && i <= m);
      xassert(1 <= k && k <= m+n);
#endif
      if (k <= m)
      {  /* new i-th column of B is k-th column of I */
         int ind[1+1];
         double val[1+1];
         ind[1] = k;
         val[1] = 1.0;
         xassert(csa->valid);
         ret = bfd_update_it(csa->bfd, i, 0, 1, ind, val);
      }
      else
      {  /* new i-th column of B is (k-m)-th column of (-A) */
         int *A_ptr = csa->A_ptr;
         int *A_ind = csa->A_ind;
         double *A_val = csa->A_val;
         double *val = csa->work1;
         int beg, end, ptr, len;
         beg = A_ptr[k-m];
         end = A_ptr[k-m+1];
         len = 0;
         for (ptr = beg; ptr < end; ptr++)
            val[++len] = - A_val[ptr];
         xassert(csa->valid);
         ret = bfd_update_it(csa->bfd, i, 0, len, &A_ind[beg-1], val);
      }
      csa->valid = (ret == 0);
      return ret;
}
#endif

#if 1 /* copied from primal */
/***********************************************************************
*  error_ftran - compute residual vector r = h - B * x
*
*  This routine computes the residual vector r = h - B * x, where B is
*  the current basis matrix, h is the vector of right-hand sides, x is
*  the solution vector. */

static void error_ftran(struct csa *csa, double h[], double x[],
      double r[])
{     int m = csa->m;
#ifdef GLP_DEBUG
      int n = csa->n;
#endif
      int *A_ptr = csa->A_ptr;
      int *A_ind = csa->A_ind;
      double *A_val = csa->A_val;
      int *head = csa->head;
      int i, k, beg, end, ptr;
      double temp;
      /* compute the residual vector:
         r = h - B * x = h - B[1] * x[1] - ... - B[m] * x[m],
         where B[1], ..., B[m] are columns of matrix B */
      memcpy(&r[1], &h[1], m * sizeof(double));
      for (i = 1; i <= m; i++)
      {  temp = x[i];
         if (temp == 0.0) continue;
         k = head[i]; /* B[i] is k-th column of (I|-A) */
#ifdef GLP_DEBUG
         xassert(1 <= k && k <= m+n);
#endif
         if (k <= m)
         {  /* B[i] is k-th column of submatrix I */
            r[k] -= temp;
         }
         else
         {  /* B[i] is (k-m)-th column of submatrix (-A) */
            beg = A_ptr[k-m];
            end = A_ptr[k-m+1];
            for (ptr = beg; ptr < end; ptr++)
               r[A_ind[ptr]] += A_val[ptr] * temp;
         }
      }
      return;
}
#endif

#if 1 /* copied from primal */
/***********************************************************************
*  refine_ftran - refine solution of B * x = h
*
*  This routine performs one iteration to refine the solution of
*  the system B * x = h, where B is the current basis matrix, h is the
*  vector of right-hand sides, x is the solution vector. */

static void refine_ftran(struct csa *csa, double h[], double x[])
{     int m = csa->m;
      double *r = csa->work1;
      double *d = csa->work1;
      int i;
      /* compute the residual vector r = h - B * x */
      error_ftran(csa, h, x, r);
      /* compute the correction vector d = inv(B) * r */
      xassert(csa->valid);
      bfd_ftran(csa->bfd, d);
      /* refine the solution vector (new x) = (old x) + d */
      for (i = 1; i <= m; i++) x[i] += d[i];
      return;
}
#endif

#if 1 /* copied from primal */
/***********************************************************************
*  error_btran - compute residual vector r = h - B'* x
*
*  This routine computes the residual vector r = h - B'* x, where B'
*  is a matrix transposed to the current basis matrix, h is the vector
*  of right-hand sides, x is the solution vector. */

static void error_btran(struct csa *csa, double h[], double x[],
      double r[])
{     int m = csa->m;
#ifdef GLP_DEBUG
      int n = csa->n;
#endif
      int *A_ptr = csa->A_ptr;
      int *A_ind = csa->A_ind;
      double *A_val = csa->A_val;
      int *head = csa->head;
      int i, k, beg, end, ptr;
      double temp;
      /* compute the residual vector r = b - B'* x */
      for (i = 1; i <= m; i++)
      {  /* r[i] := b[i] - (i-th column of B)'* x */
         k = head[i]; /* B[i] is k-th column of (I|-A) */
#ifdef GLP_DEBUG
         xassert(1 <= k && k <= m+n);
#endif
         temp = h[i];
         if (k <= m)
         {  /* B[i] is k-th column of submatrix I */
            temp -= x[k];
         }
         else
         {  /* B[i] is (k-m)-th column of submatrix (-A) */
            beg = A_ptr[k-m];
            end = A_ptr[k-m+1];
            for (ptr = beg; ptr < end; ptr++)
               temp += A_val[ptr] * x[A_ind[ptr]];
         }
         r[i] = temp;
      }
      return;
}
#endif

#if 1 /* copied from primal */
/***********************************************************************
*  refine_btran - refine solution of B'* x = h
*
*  This routine performs one iteration to refine the solution of the
*  system B'* x = h, where B' is a matrix transposed to the current
*  basis matrix, h is the vector of right-hand sides, x is the solution
*  vector. */

static void refine_btran(struct csa *csa, double h[], double x[])
{     int m = csa->m;
      double *r = csa->work1;
      double *d = csa->work1;
      int i;
      /* compute the residual vector r = h - B'* x */
      error_btran(csa, h, x, r);
      /* compute the correction vector d = inv(B') * r */
      xassert(csa->valid);
      bfd_btran(csa->bfd, d);
      /* refine the solution vector (new x) = (old x) + d */
      for (i = 1; i <= m; i++) x[i] += d[i];
      return;
}
#endif

#if 1 /* copied from primal */
/***********************************************************************
*  get_xN - determine current value of non-basic variable xN[j]
*
*  This routine returns the current value of non-basic variable xN[j],
*  which is a value of its active bound. */

static double get_xN(struct csa *csa, int j)
{     int m = csa->m;
#ifdef GLP_DEBUG
      int n = csa->n;
#endif
      double *lb = csa->lb;
      double *ub = csa->ub;
      int *head = csa->head;
      char *stat = csa->stat;
      int k;
      double xN;
#ifdef GLP_DEBUG
      xassert(1 <= j && j <= n);
#endif
      k = head[m+j]; /* x[k] = xN[j] */
#ifdef GLP_DEBUG
      xassert(1 <= k && k <= m+n);
#endif
      switch (stat[j])
      {  case GLP_NL:
            /* x[k] is on its lower bound */
            xN = lb[k]; break;
         case GLP_NU:
            /* x[k] is on its upper bound */
            xN = ub[k]; break;
         case GLP_NF:
            /* x[k] is free non-basic variable */
            xN = 0.0; break;
         case GLP_NS:
            /* x[k] is fixed non-basic variable */
            xN = lb[k]; break;
         default:
            xassert(stat != stat);
      }
      return xN;
}
#endif

#if 1 /* copied from primal */
/***********************************************************************
*  eval_beta - compute primal values of basic variables
*
*  This routine computes current primal values of all basic variables:
*
*     beta = - inv(B) * N * xN,
*
*  where B is the current basis matrix, N is a matrix built of columns
*  of matrix (I|-A) corresponding to non-basic variables, and xN is the
*  vector of current values of non-basic variables. */

static void eval_beta(struct csa *csa, double beta[])
{     int m = csa->m;
      int n = csa->n;
      int *A_ptr = csa->A_ptr;
      int *A_ind = csa->A_ind;
      double *A_val = csa->A_val;
      int *head = csa->head;
      double *h = csa->work2;
      int i, j, k, beg, end, ptr;
      double xN;
      /* compute the right-hand side vector:
         h := - N * xN = - N[1] * xN[1] - ... - N[n] * xN[n],
         where N[1], ..., N[n] are columns of matrix N */
      for (i = 1; i <= m; i++)
         h[i] = 0.0;
      for (j = 1; j <= n; j++)
      {  k = head[m+j]; /* x[k] = xN[j] */
#ifdef GLP_DEBUG
         xassert(1 <= k && k <= m+n);
#endif
         /* determine current value of xN[j] */
         xN = get_xN(csa, j);
         if (xN == 0.0) continue;
         if (k <= m)
         {  /* N[j] is k-th column of submatrix I */
            h[k] -= xN;
         }
         else
         {  /* N[j] is (k-m)-th column of submatrix (-A) */
            beg = A_ptr[k-m];
            end = A_ptr[k-m+1];
            for (ptr = beg; ptr < end; ptr++)
               h[A_ind[ptr]] += xN * A_val[ptr];
         }
      }
      /* solve system B * beta = h */
      memcpy(&beta[1], &h[1], m * sizeof(double));
      xassert(csa->valid);
      bfd_ftran(csa->bfd, beta);
      /* and refine the solution */
      refine_ftran(csa, h, beta);
      return;
}
#endif

#if 1 /* copied from primal */
/***********************************************************************
*  eval_pi - compute vector of simplex multipliers
*
*  This routine computes the vector of current simplex multipliers:
*
*     pi = inv(B') * cB,
*
*  where B' is a matrix transposed to the current basis matrix, cB is
*  a subvector of objective coefficients at basic variables. */

static void eval_pi(struct csa *csa, double pi[])
{     int m = csa->m;
      double *c = csa->coef;
      int *head = csa->head;
      double *cB = csa->work2;
      int i;
      /* construct the right-hand side vector cB */
      for (i = 1; i <= m; i++)
         cB[i] = c[head[i]];
      /* solve system B'* pi = cB */
      memcpy(&pi[1], &cB[1], m * sizeof(double));
      xassert(csa->valid);
      bfd_btran(csa->bfd, pi);
      /* and refine the solution */
      refine_btran(csa, cB, pi);
      return;
}
#endif

#if 1 /* copied from primal */
/***********************************************************************
*  eval_cost - compute reduced cost of non-basic variable xN[j]
*
*  This routine computes the current reduced cost of non-basic variable
*  xN[j]:
*
*     d[j] = cN[j] - N'[j] * pi,
*
*  where cN[j] is the objective coefficient at variable xN[j], N[j] is
*  a column of the augmented constraint matrix (I|-A) corresponding to
*  xN[j], pi is the vector of simplex multipliers. */

static double eval_cost(struct csa *csa, double pi[], int j)
{     int m = csa->m;
#ifdef GLP_DEBUG
      int n = csa->n;
#endif
      double *coef = csa->coef;
      int *head = csa->head;
      int k;
      double dj;
#ifdef GLP_DEBUG
      xassert(1 <= j && j <= n);
#endif
      k = head[m+j]; /* x[k] = xN[j] */
#ifdef GLP_DEBUG
      xassert(1 <= k && k <= m+n);
#endif
      dj = coef[k];
      if (k <= m)
      {  /* N[j] is k-th column of submatrix I */
         dj -= pi[k];
      }
      else
      {  /* N[j] is (k-m)-th column of submatrix (-A) */
         int *A_ptr = csa->A_ptr;
         int *A_ind = csa->A_ind;
         double *A_val = csa->A_val;
         int beg, end, ptr;
         beg = A_ptr[k-m];
         end = A_ptr[k-m+1];
         for (ptr = beg; ptr < end; ptr++)
            dj += A_val[ptr] * pi[A_ind[ptr]];
      }
      return dj;
}
#endif

#if 1 /* copied from primal */
/***********************************************************************
*  eval_bbar - compute and store primal values of basic variables
*
*  This routine computes primal values of all basic variables and then
*  stores them in the solution array. */

static void eval_bbar(struct csa *csa)
{     eval_beta(csa, csa->bbar);
      return;
}
#endif

#if 1 /* copied from primal */
/***********************************************************************
*  eval_cbar - compute and store reduced costs of non-basic variables
*
*  This routine computes reduced costs of all non-basic variables and
*  then stores them in the solution array. */

static void eval_cbar(struct csa *csa)
{
#ifdef GLP_DEBUG
      int m = csa->m;
#endif
      int n = csa->n;
#ifdef GLP_DEBUG
      int *head = csa->head;
#endif
      double *cbar = csa->cbar;
      double *pi = csa->work3;
      int j;
#ifdef GLP_DEBUG
      int k;
#endif
      /* compute simplex multipliers */
      eval_pi(csa, pi);
      /* compute and store reduced costs */
      for (j = 1; j <= n; j++)
      {
#ifdef GLP_DEBUG
         k = head[m+j]; /* x[k] = xN[j] */
         xassert(1 <= k && k <= m+n);
#endif
         cbar[j] = eval_cost(csa, pi, j);
      }
      return;
}
#endif

/***********************************************************************
*  reset_refsp - reset the reference space
*
*  This routine resets (redefines) the reference space used in the
*  projected steepest edge pricing algorithm. */

static void reset_refsp(struct csa *csa)
{     int m = csa->m;
      int n = csa->n;
      int *head = csa->head;
      char *refsp = csa->refsp;
      double *gamma = csa->gamma;
      int i, k;
      xassert(csa->refct == 0);
      csa->refct = 1000;
      memset(&refsp[1], 0, (m+n) * sizeof(char));
      for (i = 1; i <= m; i++)
      {  k = head[i]; /* x[k] = xB[i] */
         refsp[k] = 1;
         gamma[i] = 1.0;
      }
      return;
}

/***********************************************************************
*  eval_gamma - compute steepest edge coefficients
*
*  This routine computes the vector of steepest edge coefficients for
*  all basic variables (except free ones) using its direct definition:
*
*     gamma[i] = eta[i] +  sum   alfa[i,j]^2,  i = 1,...,m,
*                         j in C
*
*  where eta[i] = 1 means that xB[i] is in the current reference space,
*  and 0 otherwise; C is a set of non-basic non-fixed variables xN[j],
*  which are in the current reference space; alfa[i,j] are elements of
*  the current simplex table.
*
*  NOTE: The routine is intended only for debugginig purposes. */

static void eval_gamma(struct csa *csa, double gamma[])
{     int m = csa->m;
      int n = csa->n;
      char *type = csa->type;
      int *head = csa->head;
      char *refsp = csa->refsp;
      double *alfa = csa->work3;
      double *h = csa->work3;
      int i, j, k;
      /* gamma[i] := eta[i] (or 1, if xB[i] is free) */
      for (i = 1; i <= m; i++)
      {  k = head[i]; /* x[k] = xB[i] */
#ifdef GLP_DEBUG
         xassert(1 <= k && k <= m+n);
#endif
         if (type[k] == GLP_FR)
            gamma[i] = 1.0;
         else
            gamma[i] = (refsp[k] ? 1.0 : 0.0);
      }
      /* compute columns of the current simplex table */
      for (j = 1; j <= n; j++)
      {  k = head[m+j]; /* x[k] = xN[j] */
#ifdef GLP_DEBUG
         xassert(1 <= k && k <= m+n);
#endif
         /* skip column, if xN[j] is not in C */
         if (!refsp[k]) continue;
#ifdef GLP_DEBUG
         /* set C must not contain fixed variables */
         xassert(type[k] != GLP_FX);
#endif
         /* construct the right-hand side vector h = - N[j] */
         for (i = 1; i <= m; i++)
            h[i] = 0.0;
         if (k <= m)
         {  /* N[j] is k-th column of submatrix I */
            h[k] = -1.0;
         }
         else
         {  /* N[j] is (k-m)-th column of submatrix (-A) */
            int *A_ptr = csa->A_ptr;
            int *A_ind = csa->A_ind;
            double *A_val = csa->A_val;
            int beg, end, ptr;
            beg = A_ptr[k-m];
            end = A_ptr[k-m+1];
            for (ptr = beg; ptr < end; ptr++)
               h[A_ind[ptr]] = A_val[ptr];
         }
         /* solve system B * alfa = h */
         xassert(csa->valid);
         bfd_ftran(csa->bfd, alfa);
         /* gamma[i] := gamma[i] + alfa[i,j]^2 */
         for (i = 1; i <= m; i++)
         {  k = head[i]; /* x[k] = xB[i] */
            if (type[k] != GLP_FR)
               gamma[i] += alfa[i] * alfa[i];
         }
      }
      return;
}

/***********************************************************************
*  chuzr - choose basic variable (row of the simplex table)
*
*  This routine chooses basic variable xB[p] having largest weighted
*  bound violation:
*
*     |r[p]| / sqrt(gamma[p]) = max  |r[i]| / sqrt(gamma[i]),
*                              i in I
*
*            / lB[i] - beta[i], if beta[i] < lB[i]
*            |
*     r[i] = < 0,               if lB[i] <= beta[i] <= uB[i]
*            |
*            \ uB[i] - beta[i], if beta[i] > uB[i]
*
*  where beta[i] is primal value of xB[i] in the current basis, lB[i]
*  and uB[i] are lower and upper bounds of xB[i], I is a subset of
*  eligible basic variables, which significantly violates their bounds,
*  gamma[i] is the steepest edge coefficient.
*
*  If |r[i]| is less than a specified tolerance, xB[i] is not included
*  in I and therefore ignored.
*
*  If I is empty and no variable has been chosen, p is set to 0. */

static void chuzr(struct csa *csa, double tol_bnd)
{     int m = csa->m;
#ifdef GLP_DEBUG
      int n = csa->n;
#endif
      char *type = csa->type;
      double *lb = csa->lb;
      double *ub = csa->ub;
      int *head = csa->head;
      double *bbar = csa->bbar;
      double *gamma = csa->gamma;
      int i, k, p;
      double delta, best, eps, ri, temp;
      /* nothing is chosen so far */
      p = 0, delta = 0.0, best = 0.0;
      /* look through the list of basic variables */
      for (i = 1; i <= m; i++)
      {  k = head[i]; /* x[k] = xB[i] */
#ifdef GLP_DEBUG
         xassert(1 <= k && k <= m+n);
#endif
         /* determine bound violation ri[i] */
         ri = 0.0;
         if (type[k] == GLP_LO || type[k] == GLP_DB ||
             type[k] == GLP_FX)
         {  /* xB[i] has lower bound */
            eps = tol_bnd * (1.0 + kappa * fabs(lb[k]));
            if (bbar[i] < lb[k] - eps)
            {  /* and significantly violates it */
               ri = lb[k] - bbar[i];
            }
         }
         if (type[k] == GLP_UP || type[k] == GLP_DB ||
             type[k] == GLP_FX)
         {  /* xB[i] has upper bound */
            eps = tol_bnd * (1.0 + kappa * fabs(ub[k]));
            if (bbar[i] > ub[k] + eps)
            {  /* and significantly violates it */
               ri = ub[k] - bbar[i];
            }
         }
         /* if xB[i] is not eligible, skip it */
         if (ri == 0.0) continue;
         /* xB[i] is eligible basic variable; choose one with largest
            weighted bound violation */
#ifdef GLP_DEBUG
         xassert(gamma[i] >= 0.0);
#endif
         temp = gamma[i];
         if (temp < DBL_EPSILON) temp = DBL_EPSILON;
         temp = (ri * ri) / temp;
         if (best < temp)
            p = i, delta = ri, best = temp;
      }
      /* store the index of basic variable xB[p] chosen and its change
         in the adjacent basis */
      csa->p = p;
      csa->delta = delta;
      return;
}

#if 1 /* copied from primal */
/***********************************************************************
*  eval_rho - compute pivot row of the inverse
*
*  This routine computes the pivot (p-th) row of the inverse inv(B),
*  which corresponds to basic variable xB[p] chosen:
*
*     rho = inv(B') * e[p],
*
*  where B' is a matrix transposed to the current basis matrix, e[p]
*  is unity vector. */

static void eval_rho(struct csa *csa, double rho[])
{     int m = csa->m;
      int p = csa->p;
      double *e = rho;
      int i;
#ifdef GLP_DEBUG
      xassert(1 <= p && p <= m);
#endif
      /* construct the right-hand side vector e[p] */
      for (i = 1; i <= m; i++)
         e[i] = 0.0;
      e[p] = 1.0;
      /* solve system B'* rho = e[p] */
      xassert(csa->valid);
      bfd_btran(csa->bfd, rho);
      return;
}
#endif

#if 1 /* copied from primal */
/***********************************************************************
*  refine_rho - refine pivot row of the inverse
*
*  This routine refines the pivot row of the inverse inv(B) assuming
*  that it was previously computed by the routine eval_rho. */

static void refine_rho(struct csa *csa, double rho[])
{     int m = csa->m;
      int p = csa->p;
      double *e = csa->work3;
      int i;
#ifdef GLP_DEBUG
      xassert(1 <= p && p <= m);
#endif
      /* construct the right-hand side vector e[p] */
      for (i = 1; i <= m; i++)
         e[i] = 0.0;
      e[p] = 1.0;
      /* refine solution of B'* rho = e[p] */
      refine_btran(csa, e, rho);
      return;
}
#endif

#if 1 /* 06/IV-2009 */
/***********************************************************************
*  eval_trow - compute pivot row of the simplex table
*
*  This routine computes the pivot row of the simplex table, which
*  corresponds to basic variable xB[p] chosen.
*
*  The pivot row is the following vector:
*
*     trow = T'* e[p] = - N'* inv(B') * e[p] = - N' * rho,
*
*  where rho is the pivot row of the inverse inv(B) previously computed
*  by the routine eval_rho.
*
*  Note that elements of the pivot row corresponding to fixed non-basic
*  variables are not computed.
*
*  NOTES
*
*  Computing pivot row of the simplex table is one of the most time
*  consuming operations, and for some instances it may take more than
*  50% of the total solution time.
*
*  In the current implementation there are two routines to compute the
*  pivot row. The routine eval_trow1 computes elements of the pivot row
*  as inner products of columns of the matrix N and the vector rho; it
*  is used when the vector rho is relatively dense. The routine
*  eval_trow2 computes the pivot row as a linear combination of rows of
*  the matrix N; it is used when the vector rho is relatively sparse. */

static void eval_trow1(struct csa *csa, double rho[])
{     int m = csa->m;
      int n = csa->n;
      int *A_ptr = csa->A_ptr;
      int *A_ind = csa->A_ind;
      double *A_val = csa->A_val;
      int *head = csa->head;
      char *stat = csa->stat;
      int *trow_ind = csa->trow_ind;
      double *trow_vec = csa->trow_vec;
      int j, k, beg, end, ptr, nnz;
      double temp;
      /* compute the pivot row as inner products of columns of the
         matrix N and vector rho: trow[j] = - rho * N[j] */
      nnz = 0;
      for (j = 1; j <= n; j++)
      {  if (stat[j] == GLP_NS)
         {  /* xN[j] is fixed */
            trow_vec[j] = 0.0;
            continue;
         }
         k = head[m+j]; /* x[k] = xN[j] */
         if (k <= m)
         {  /* N[j] is k-th column of submatrix I */
            temp = - rho[k];
         }
         else
         {  /* N[j] is (k-m)-th column of submatrix (-A) */
            beg = A_ptr[k-m], end = A_ptr[k-m+1];
            temp = 0.0;
            for (ptr = beg; ptr < end; ptr++)
               temp += rho[A_ind[ptr]] * A_val[ptr];
         }
         if (temp != 0.0)
            trow_ind[++nnz] = j;
         trow_vec[j] = temp;
      }
      csa->trow_nnz = nnz;
      return;
}

static void eval_trow2(struct csa *csa, double rho[])
{     int m = csa->m;
      int n = csa->n;
      int *AT_ptr = csa->AT_ptr;
      int *AT_ind = csa->AT_ind;
      double *AT_val = csa->AT_val;
      int *bind = csa->bind;
      char *stat = csa->stat;
      int *trow_ind = csa->trow_ind;
      double *trow_vec = csa->trow_vec;
      int i, j, beg, end, ptr, nnz;
      double temp;
      /* clear the pivot row */
      for (j = 1; j <= n; j++)
         trow_vec[j] = 0.0;
      /* compute the pivot row as a linear combination of rows of the
         matrix N: trow = - rho[1] * N'[1] - ... - rho[m] * N'[m] */
      for (i = 1; i <= m; i++)
      {  temp = rho[i];
         if (temp == 0.0) continue;
         /* trow := trow - rho[i] * N'[i] */
         j = bind[i] - m; /* x[i] = xN[j] */
         if (j >= 1 && stat[j] != GLP_NS)
            trow_vec[j] -= temp;
         beg = AT_ptr[i], end = AT_ptr[i+1];
         for (ptr = beg; ptr < end; ptr++)
         {  j = bind[m + AT_ind[ptr]] - m; /* x[k] = xN[j] */
            if (j >= 1 && stat[j] != GLP_NS)
               trow_vec[j] += temp * AT_val[ptr];
         }
      }
      /* construct sparse pattern of the pivot row */
      nnz = 0;
      for (j = 1; j <= n; j++)
      {  if (trow_vec[j] != 0.0)
            trow_ind[++nnz] = j;
      }
      csa->trow_nnz = nnz;
      return;
}

static void eval_trow(struct csa *csa, double rho[])
{     int m = csa->m;
      int i, nnz;
      double dens;
      /* determine the density of the vector rho */
      nnz = 0;
      for (i = 1; i <= m; i++)
         if (rho[i] != 0.0) nnz++;
      dens = (double)nnz / (double)m;
      if (dens >= 0.20)
      {  /* rho is relatively dense */
         eval_trow1(csa, rho);
      }
      else
      {  /* rho is relatively sparse */
         eval_trow2(csa, rho);
      }
      return;
}
#endif

/***********************************************************************
*  sort_trow - sort pivot row of the simplex table
*
*  This routine reorders the list of non-zero elements of the pivot
*  row to put significant elements, whose magnitude is not less than
*  a specified tolerance, in front of the list, and stores the number
*  of significant elements in trow_num. */

static void sort_trow(struct csa *csa, double tol_piv)
{
#ifdef GLP_DEBUG
      int n = csa->n;
      char *stat = csa->stat;
#endif
      int nnz = csa->trow_nnz;
      int *trow_ind = csa->trow_ind;
      double *trow_vec = csa->trow_vec;
      int j, num, pos;
      double big, eps, temp;
      /* compute infinity (maximum) norm of the row */
      big = 0.0;
      for (pos = 1; pos <= nnz; pos++)
      {
#ifdef GLP_DEBUG
         j = trow_ind[pos];
         xassert(1 <= j && j <= n);
         xassert(stat[j] != GLP_NS);
#endif
         temp = fabs(trow_vec[trow_ind[pos]]);
         if (big < temp) big = temp;
      }
      csa->trow_max = big;
      /* determine absolute pivot tolerance */
      eps = tol_piv * (1.0 + 0.01 * big);
      /* move significant row components to the front of the list */
      for (num = 0; num < nnz; )
      {  j = trow_ind[nnz];
         if (fabs(trow_vec[j]) < eps)
            nnz--;
         else
         {  num++;
            trow_ind[nnz] = trow_ind[num];
            trow_ind[num] = j;
         }
      }
      csa->trow_num = num;
      return;
}

#ifdef GLP_LONG_STEP /* 07/IV-2009 */
static int ls_func(const void *p1_, const void *p2_)
{     const struct bkpt *p1 = p1_, *p2 = p2_;
      if (p1->t < p2->t) return -1;
      if (p1->t > p2->t) return +1;
      return 0;
}

static int ls_func1(const void *p1_, const void *p2_)
{     const struct bkpt *p1 = p1_, *p2 = p2_;
      if (p1->dz < p2->dz) return -1;
      if (p1->dz > p2->dz) return +1;
      return 0;
}

static void long_step(struct csa *csa)
{     int m = csa->m;
#ifdef GLP_DEBUG
      int n = csa->n;
#endif
      char *type = csa->type;
      double *lb = csa->lb;
      double *ub = csa->ub;
      int *head = csa->head;
      char *stat = csa->stat;
      double *cbar = csa->cbar;
      double delta = csa->delta;
      int *trow_ind = csa->trow_ind;
      double *trow_vec = csa->trow_vec;
      int trow_num = csa->trow_num;
      struct bkpt *bkpt = csa->bkpt;
      int j, k, kk, nbps, pos;
      double alfa, s, slope, dzmax;
      /* delta > 0 means that xB[p] violates its lower bound, so to
         increase the dual objective lambdaB[p] must increase;
         delta < 0 means that xB[p] violates its upper bound, so to
         increase the dual objective lambdaB[p] must decrease */
      /* s := sign(delta) */
      s = (delta > 0.0 ? +1.0 : -1.0);
      /* determine breakpoints of the dual objective */
      nbps = 0;
      for (pos = 1; pos <= trow_num; pos++)
      {  j = trow_ind[pos];
#ifdef GLP_DEBUG
         xassert(1 <= j && j <= n);
         xassert(stat[j] != GLP_NS);
#endif
         /* if there is free non-basic variable, switch to the standard
            ratio test */
         if (stat[j] == GLP_NF)
         {  nbps = 0;
            goto done;
         }
         /* lambdaN[j] = ... - alfa * t - ..., where t = s * lambdaB[i]
            is the dual ray parameter, t >= 0 */
         alfa = s * trow_vec[j];
#ifdef GLP_DEBUG
         xassert(alfa != 0.0);
         xassert(stat[j] == GLP_NL || stat[j] == GLP_NU);
#endif
         if (alfa > 0.0 && stat[j] == GLP_NL ||
             alfa < 0.0 && stat[j] == GLP_NU)
         {  /* either lambdaN[j] >= 0 (if stat = GLP_NL) and decreases
               or lambdaN[j] <= 0 (if stat = GLP_NU) and increases; in
               both cases we have a breakpoint */
            nbps++;
#ifdef GLP_DEBUG
            xassert(nbps <= n);
#endif
            bkpt[nbps].j = j;
            bkpt[nbps].t = cbar[j] / alfa;
/*
if (stat[j] == GLP_NL && cbar[j] < 0.0 ||
    stat[j] == GLP_NU && cbar[j] > 0.0)
xprintf("%d %g\n", stat[j], cbar[j]);
*/
            /* if t is negative, replace it by exact zero (see comments
               in the routine chuzc) */
            if (bkpt[nbps].t < 0.0) bkpt[nbps].t = 0.0;
         }
      }
      /* if there are less than two breakpoints, switch to the standard
         ratio test */
      if (nbps < 2)
      {  nbps = 0;
         goto done;
      }
      /* sort breakpoints by ascending the dual ray parameter, t */
      qsort(&bkpt[1], nbps, sizeof(struct bkpt), ls_func);
      /* determine last breakpoint, at which the dual objective still
         greater than at t = 0 */
      dzmax = 0.0;
      slope = fabs(delta); /* initial slope */
      for (kk = 1; kk <= nbps; kk++)
      {  if (kk == 1)
            bkpt[kk].dz =
               0.0 + slope * (bkpt[kk].t - 0.0);
         else
            bkpt[kk].dz =
               bkpt[kk-1].dz + slope * (bkpt[kk].t - bkpt[kk-1].t);
         if (dzmax < bkpt[kk].dz)
            dzmax = bkpt[kk].dz;
         else if (bkpt[kk].dz < 0.05 * (1.0 + dzmax))
         {  nbps = kk - 1;
            break;
         }
         j = bkpt[kk].j;
         k = head[m+j]; /* x[k] = xN[j] */
         if (type[k] == GLP_DB)
            slope -= fabs(trow_vec[j]) * (ub[k] - lb[k]);
         else
         {  nbps = kk;
            break;
         }
      }
      /* if there are less than two breakpoints, switch to the standard
         ratio test */
      if (nbps < 2)
      {  nbps = 0;
         goto done;
      }
      /* sort breakpoints by ascending the dual change, dz */
      qsort(&bkpt[1], nbps, sizeof(struct bkpt), ls_func1);
/*
for (kk = 1; kk <= nbps; kk++)
xprintf("%d; t = %g; dz = %g\n", kk, bkpt[kk].t, bkpt[kk].dz);
*/
done: csa->nbps = nbps;
      return;
}
#endif

/***********************************************************************
*  chuzc - choose non-basic variable (column of the simplex table)
*
*  This routine chooses non-basic variable xN[q], which being entered
*  in the basis keeps dual feasibility of the basic solution.
*
*  The parameter rtol is a relative tolerance used to relax zero bounds
*  of reduced costs of non-basic variables. If rtol = 0, the routine
*  implements the standard ratio test. Otherwise, if rtol > 0, the
*  routine implements Harris' two-pass ratio test. In the latter case
*  rtol should be about three times less than a tolerance used to check
*  dual feasibility. */

static void chuzc(struct csa *csa, double rtol)
{
#ifdef GLP_DEBUG
      int m = csa->m;
      int n = csa->n;
#endif
      char *stat = csa->stat;
      double *cbar = csa->cbar;
#ifdef GLP_DEBUG
      int p = csa->p;
#endif
      double delta = csa->delta;
      int *trow_ind = csa->trow_ind;
      double *trow_vec = csa->trow_vec;
      int trow_num = csa->trow_num;
      int j, pos, q;
      double alfa, big, s, t, teta, tmax;
#ifdef GLP_DEBUG
      xassert(1 <= p && p <= m);
#endif
      /* delta > 0 means that xB[p] violates its lower bound and goes
         to it in the adjacent basis, so lambdaB[p] is increasing from
         its lower zero bound;
         delta < 0 means that xB[p] violates its upper bound and goes
         to it in the adjacent basis, so lambdaB[p] is decreasing from
         its upper zero bound */
#ifdef GLP_DEBUG
      xassert(delta != 0.0);
#endif
      /* s := sign(delta) */
      s = (delta > 0.0 ? +1.0 : -1.0);
      /*** FIRST PASS ***/
      /* nothing is chosen so far */
      q = 0, teta = DBL_MAX, big = 0.0;
      /* walk through significant elements of the pivot row */
      for (pos = 1; pos <= trow_num; pos++)
      {  j = trow_ind[pos];
#ifdef GLP_DEBUG
         xassert(1 <= j && j <= n);
#endif
         alfa = s * trow_vec[j];
#ifdef GLP_DEBUG
         xassert(alfa != 0.0);
#endif
         /* lambdaN[j] = ... - alfa * lambdaB[p] - ..., and due to s we
            need to consider only increasing lambdaB[p] */
         if (alfa > 0.0)
         {  /* lambdaN[j] is decreasing */
            if (stat[j] == GLP_NL || stat[j] == GLP_NF)
            {  /* lambdaN[j] has zero lower bound */
               t = (cbar[j] + rtol) / alfa;
            }
            else
            {  /* lambdaN[j] has no lower bound */
               continue;
            }
         }
         else
         {  /* lambdaN[j] is increasing */
            if (stat[j] == GLP_NU || stat[j] == GLP_NF)
            {  /* lambdaN[j] has zero upper bound */
               t = (cbar[j] - rtol) / alfa;
            }
            else
            {  /* lambdaN[j] has no upper bound */
               continue;
            }
         }
         /* t is a change of lambdaB[p], on which lambdaN[j] reaches
            its zero bound (possibly relaxed); since the basic solution
            is assumed to be dual feasible, t has to be non-negative by
            definition; however, it may happen that lambdaN[j] slightly
            (i.e. within a tolerance) violates its zero bound, that
            leads to negative t; in the latter case, if xN[j] is chosen,
            negative t means that lambdaB[p] changes in wrong direction
            that may cause wrong results on updating reduced costs;
            thus, if t is negative, we should replace it by exact zero
            assuming that lambdaN[j] is exactly on its zero bound, and
            violation appears due to round-off errors */
         if (t < 0.0) t = 0.0;
         /* apply minimal ratio test */
         if (teta > t || teta == t && big < fabs(alfa))
            q = j, teta = t, big = fabs(alfa);
      }
      /* the second pass is skipped in the following cases: */
      /* if the standard ratio test is used */
      if (rtol == 0.0) goto done;
      /* if no non-basic variable has been chosen on the first pass */
      if (q == 0) goto done;
      /* if lambdaN[q] prevents lambdaB[p] from any change */
      if (teta == 0.0) goto done;
      /*** SECOND PASS ***/
      /* here tmax is a maximal change of lambdaB[p], on which the
         solution remains dual feasible within a tolerance */
#if 0
      tmax = (1.0 + 10.0 * DBL_EPSILON) * teta;
#else
      tmax = teta;
#endif
      /* nothing is chosen so far */
      q = 0, teta = DBL_MAX, big = 0.0;
      /* walk through significant elements of the pivot row */
      for (pos = 1; pos <= trow_num; pos++)
      {  j = trow_ind[pos];
#ifdef GLP_DEBUG
         xassert(1 <= j && j <= n);
#endif
         alfa = s * trow_vec[j];
#ifdef GLP_DEBUG
         xassert(alfa != 0.0);
#endif
         /* lambdaN[j] = ... - alfa * lambdaB[p] - ..., and due to s we
            need to consider only increasing lambdaB[p] */
         if (alfa > 0.0)
         {  /* lambdaN[j] is decreasing */
            if (stat[j] == GLP_NL || stat[j] == GLP_NF)
            {  /* lambdaN[j] has zero lower bound */
               t = cbar[j] / alfa;
            }
            else
            {  /* lambdaN[j] has no lower bound */
               continue;
            }
         }
         else
         {  /* lambdaN[j] is increasing */
            if (stat[j] == GLP_NU || stat[j] == GLP_NF)
            {  /* lambdaN[j] has zero upper bound */
               t = cbar[j] / alfa;
            }
            else
            {  /* lambdaN[j] has no upper bound */
               continue;
            }
         }
         /* (see comments for the first pass) */
         if (t < 0.0) t = 0.0;
         /* t is a change of lambdaB[p], on which lambdaN[j] reaches
            its zero (lower or upper) bound; if t <= tmax, all reduced
            costs can violate their zero bounds only within relaxation
            tolerance rtol, so we can choose non-basic variable having
            largest influence coefficient to avoid possible numerical
            instability */
         if (t <= tmax && big < fabs(alfa))
            q = j, teta = t, big = fabs(alfa);
      }
      /* something must be chosen on the second pass */
      xassert(q != 0);
done: /* store the index of non-basic variable xN[q] chosen */
      csa->q = q;
      /* store reduced cost of xN[q] in the adjacent basis */
      csa->new_dq = s * teta;
      return;
}

#if 1 /* copied from primal */
/***********************************************************************
*  eval_tcol - compute pivot column of the simplex table
*
*  This routine computes the pivot column of the simplex table, which
*  corresponds to non-basic variable xN[q] chosen.
*
*  The pivot column is the following vector:
*
*     tcol = T * e[q] = - inv(B) * N * e[q] = - inv(B) * N[q],
*
*  where B is the current basis matrix, N[q] is a column of the matrix
*  (I|-A) corresponding to variable xN[q]. */

static void eval_tcol(struct csa *csa)
{     int m = csa->m;
#ifdef GLP_DEBUG
      int n = csa->n;
#endif
      int *head = csa->head;
      int q = csa->q;
      int *tcol_ind = csa->tcol_ind;
      double *tcol_vec = csa->tcol_vec;
      double *h = csa->tcol_vec;
      int i, k, nnz;
#ifdef GLP_DEBUG
      xassert(1 <= q && q <= n);
#endif
      k = head[m+q]; /* x[k] = xN[q] */
#ifdef GLP_DEBUG
      xassert(1 <= k && k <= m+n);
#endif
      /* construct the right-hand side vector h = - N[q] */
      for (i = 1; i <= m; i++)
         h[i] = 0.0;
      if (k <= m)
      {  /* N[q] is k-th column of submatrix I */
         h[k] = -1.0;
      }
      else
      {  /* N[q] is (k-m)-th column of submatrix (-A) */
         int *A_ptr = csa->A_ptr;
         int *A_ind = csa->A_ind;
         double *A_val = csa->A_val;
         int beg, end, ptr;
         beg = A_ptr[k-m];
         end = A_ptr[k-m+1];
         for (ptr = beg; ptr < end; ptr++)
            h[A_ind[ptr]] = A_val[ptr];
      }
      /* solve system B * tcol = h */
      xassert(csa->valid);
      bfd_ftran(csa->bfd, tcol_vec);
      /* construct sparse pattern of the pivot column */
      nnz = 0;
      for (i = 1; i <= m; i++)
      {  if (tcol_vec[i] != 0.0)
            tcol_ind[++nnz] = i;
      }
      csa->tcol_nnz = nnz;
      return;
}
#endif

#if 1 /* copied from primal */
/***********************************************************************
*  refine_tcol - refine pivot column of the simplex table
*
*  This routine refines the pivot column of the simplex table assuming
*  that it was previously computed by the routine eval_tcol. */

static void refine_tcol(struct csa *csa)
{     int m = csa->m;
#ifdef GLP_DEBUG
      int n = csa->n;
#endif
      int *head = csa->head;
      int q = csa->q;
      int *tcol_ind = csa->tcol_ind;
      double *tcol_vec = csa->tcol_vec;
      double *h = csa->work3;
      int i, k, nnz;
#ifdef GLP_DEBUG
      xassert(1 <= q && q <= n);
#endif
      k = head[m+q]; /* x[k] = xN[q] */
#ifdef GLP_DEBUG
      xassert(1 <= k && k <= m+n);
#endif
      /* construct the right-hand side vector h = - N[q] */
      for (i = 1; i <= m; i++)
         h[i] = 0.0;
      if (k <= m)
      {  /* N[q] is k-th column of submatrix I */
         h[k] = -1.0;
      }
      else
      {  /* N[q] is (k-m)-th column of submatrix (-A) */
         int *A_ptr = csa->A_ptr;
         int *A_ind = csa->A_ind;
         double *A_val = csa->A_val;
         int beg, end, ptr;
         beg = A_ptr[k-m];
         end = A_ptr[k-m+1];
         for (ptr = beg; ptr < end; ptr++)
            h[A_ind[ptr]] = A_val[ptr];
      }
      /* refine solution of B * tcol = h */
      refine_ftran(csa, h, tcol_vec);
      /* construct sparse pattern of the pivot column */
      nnz = 0;
      for (i = 1; i <= m; i++)
      {  if (tcol_vec[i] != 0.0)
            tcol_ind[++nnz] = i;
      }
      csa->tcol_nnz = nnz;
      return;
}
#endif

/***********************************************************************
*  update_cbar - update reduced costs of non-basic variables
*
*  This routine updates reduced costs of all (except fixed) non-basic
*  variables for the adjacent basis. */

static void update_cbar(struct csa *csa)
{
#ifdef GLP_DEBUG
      int n = csa->n;
#endif
      double *cbar = csa->cbar;
      int trow_nnz = csa->trow_nnz;
      int *trow_ind = csa->trow_ind;
      double *trow_vec = csa->trow_vec;
      int q = csa->q;
      double new_dq = csa->new_dq;
      int j, pos;
#ifdef GLP_DEBUG
      xassert(1 <= q && q <= n);
#endif
      /* set new reduced cost of xN[q] */
      cbar[q] = new_dq;
      /* update reduced costs of other non-basic variables */
      if (new_dq == 0.0) goto done;
      for (pos = 1; pos <= trow_nnz; pos++)
      {  j = trow_ind[pos];
#ifdef GLP_DEBUG
         xassert(1 <= j && j <= n);
#endif
         if (j != q)
            cbar[j] -= trow_vec[j] * new_dq;
      }
done: return;
}

/***********************************************************************
*  update_bbar - update values of basic variables
*
*  This routine updates values of all basic variables for the adjacent
*  basis. */

static void update_bbar(struct csa *csa)
{
#ifdef GLP_DEBUG
      int m = csa->m;
      int n = csa->n;
#endif
      double *bbar = csa->bbar;
      int p = csa->p;
      double delta = csa->delta;
      int q = csa->q;
      int tcol_nnz = csa->tcol_nnz;
      int *tcol_ind = csa->tcol_ind;
      double *tcol_vec = csa->tcol_vec;
      int i, pos;
      double teta;
#ifdef GLP_DEBUG
      xassert(1 <= p && p <= m);
      xassert(1 <= q && q <= n);
#endif
      /* determine the change of xN[q] in the adjacent basis */
#ifdef GLP_DEBUG
      xassert(tcol_vec[p] != 0.0);
#endif
      teta = delta / tcol_vec[p];
      /* set new primal value of xN[q] */
      bbar[p] = get_xN(csa, q) + teta;
      /* update primal values of other basic variables */
      if (teta == 0.0) goto done;
      for (pos = 1; pos <= tcol_nnz; pos++)
      {  i = tcol_ind[pos];
#ifdef GLP_DEBUG
         xassert(1 <= i && i <= m);
#endif
         if (i != p)
            bbar[i] += tcol_vec[i] * teta;
      }
done: return;
}

/***********************************************************************
*  update_gamma - update steepest edge coefficients
*
*  This routine updates steepest-edge coefficients for the adjacent
*  basis. */

static void update_gamma(struct csa *csa)
{     int m = csa->m;
#ifdef GLP_DEBUG
      int n = csa->n;
#endif
      char *type = csa->type;
      int *head = csa->head;
      char *refsp = csa->refsp;
      double *gamma = csa->gamma;
      int p = csa->p;
      int trow_nnz = csa->trow_nnz;
      int *trow_ind = csa->trow_ind;
      double *trow_vec = csa->trow_vec;
      int q = csa->q;
      int tcol_nnz = csa->tcol_nnz;
      int *tcol_ind = csa->tcol_ind;
      double *tcol_vec = csa->tcol_vec;
      double *u = csa->work3;
      int i, j, k,pos;
      double gamma_p, eta_p, pivot, t, t1, t2;
#ifdef GLP_DEBUG
      xassert(1 <= p && p <= m);
      xassert(1 <= q && q <= n);
#endif
      /* the basis changes, so decrease the count */
      xassert(csa->refct > 0);
      csa->refct--;
      /* recompute gamma[p] for the current basis more accurately and
         compute auxiliary vector u */
#ifdef GLP_DEBUG
      xassert(type[head[p]] != GLP_FR);
#endif
      gamma_p = eta_p = (refsp[head[p]] ? 1.0 : 0.0);
      for (i = 1; i <= m; i++) u[i] = 0.0;
      for (pos = 1; pos <= trow_nnz; pos++)
      {  j = trow_ind[pos];
#ifdef GLP_DEBUG
         xassert(1 <= j && j <= n);
#endif
         k = head[m+j]; /* x[k] = xN[j] */
#ifdef GLP_DEBUG
         xassert(1 <= k && k <= m+n);
         xassert(type[k] != GLP_FX);
#endif
         if (!refsp[k]) continue;
         t = trow_vec[j];
         gamma_p += t * t;
         /* u := u + N[j] * delta[j] * trow[j] */
         if (k <= m)
         {  /* N[k] = k-j stolbec submatrix I */
            u[k] += t;
         }
         else
         {  /* N[k] = k-m-k stolbec (-A) */
            int *A_ptr = csa->A_ptr;
            int *A_ind = csa->A_ind;
            double *A_val = csa->A_val;
            int beg, end, ptr;
            beg = A_ptr[k-m];
            end = A_ptr[k-m+1];
            for (ptr = beg; ptr < end; ptr++)
               u[A_ind[ptr]] -= t * A_val[ptr];
         }
      }
      xassert(csa->valid);
      bfd_ftran(csa->bfd, u);
      /* update gamma[i] for other basic variables (except xB[p] and
         free variables) */
      pivot = tcol_vec[p];
#ifdef GLP_DEBUG
      xassert(pivot != 0.0);
#endif
      for (pos = 1; pos <= tcol_nnz; pos++)
      {  i = tcol_ind[pos];
#ifdef GLP_DEBUG
         xassert(1 <= i && i <= m);
#endif
         k = head[i];
#ifdef GLP_DEBUG
         xassert(1 <= k && k <= m+n);
#endif
         /* skip xB[p] */
         if (i == p) continue;
         /* skip free basic variable */
         if (type[head[i]] == GLP_FR)
         {
#ifdef GLP_DEBUG
            xassert(gamma[i] == 1.0);
#endif
            continue;
         }
         /* compute gamma[i] for the adjacent basis */
         t = tcol_vec[i] / pivot;
         t1 = gamma[i] + t * t * gamma_p + 2.0 * t * u[i];
         t2 = (refsp[k] ? 1.0 : 0.0) + eta_p * t * t;
         gamma[i] = (t1 >= t2 ? t1 : t2);
         /* (though gamma[i] can be exact zero, because the reference
            space does not include non-basic fixed variables) */
         if (gamma[i] < DBL_EPSILON) gamma[i] = DBL_EPSILON;
      }
      /* compute gamma[p] for the adjacent basis */
      if (type[head[m+q]] == GLP_FR)
         gamma[p] = 1.0;
      else
      {  gamma[p] = gamma_p / (pivot * pivot);
         if (gamma[p] < DBL_EPSILON) gamma[p] = DBL_EPSILON;
      }
      /* if xB[p], which becomes xN[q] in the adjacent basis, is fixed
         and belongs to the reference space, remove it from there, and
         change all gamma's appropriately */
      k = head[p];
      if (type[k] == GLP_FX && refsp[k])
      {  refsp[k] = 0;
         for (pos = 1; pos <= tcol_nnz; pos++)
         {  i = tcol_ind[pos];
            if (i == p)
            {  if (type[head[m+q]] == GLP_FR) continue;
               t = 1.0 / tcol_vec[p];
            }
            else
            {  if (type[head[i]] == GLP_FR) continue;
               t = tcol_vec[i] / tcol_vec[p];
            }
            gamma[i] -= t * t;
            if (gamma[i] < DBL_EPSILON) gamma[i] = DBL_EPSILON;
         }
      }
      return;
}

#if 1 /* copied from primal */
/***********************************************************************
*  err_in_bbar - compute maximal relative error in primal solution
*
*  This routine returns maximal relative error:
*
*     max |beta[i] - bbar[i]| / (1 + |beta[i]|),
*
*  where beta and bbar are, respectively, directly computed and the
*  current (updated) values of basic variables.
*
*  NOTE: The routine is intended only for debugginig purposes. */

static double err_in_bbar(struct csa *csa)
{     int m = csa->m;
      double *bbar = csa->bbar;
      int i;
      double e, emax, *beta;
      beta = xcalloc(1+m, sizeof(double));
      eval_beta(csa, beta);
      emax = 0.0;
      for (i = 1; i <= m; i++)
      {  e = fabs(beta[i] - bbar[i]) / (1.0 + fabs(beta[i]));
         if (emax < e) emax = e;
      }
      xfree(beta);
      return emax;
}
#endif

#if 1 /* copied from primal */
/***********************************************************************
*  err_in_cbar - compute maximal relative error in dual solution
*
*  This routine returns maximal relative error:
*
*     max |cost[j] - cbar[j]| / (1 + |cost[j]|),
*
*  where cost and cbar are, respectively, directly computed and the
*  current (updated) reduced costs of non-basic non-fixed variables.
*
*  NOTE: The routine is intended only for debugginig purposes. */

static double err_in_cbar(struct csa *csa)
{     int m = csa->m;
      int n = csa->n;
      char *stat = csa->stat;
      double *cbar = csa->cbar;
      int j;
      double e, emax, cost, *pi;
      pi = xcalloc(1+m, sizeof(double));
      eval_pi(csa, pi);
      emax = 0.0;
      for (j = 1; j <= n; j++)
      {  if (stat[j] == GLP_NS) continue;
         cost = eval_cost(csa, pi, j);
         e = fabs(cost - cbar[j]) / (1.0 + fabs(cost));
         if (emax < e) emax = e;
      }
      xfree(pi);
      return emax;
}
#endif

/***********************************************************************
*  err_in_gamma - compute maximal relative error in steepest edge cff.
*
*  This routine returns maximal relative error:
*
*     max |gamma'[j] - gamma[j]| / (1 + |gamma'[j]),
*
*  where gamma'[j] and gamma[j] are, respectively, directly computed
*  and the current (updated) steepest edge coefficients for non-basic
*  non-fixed variable x[j].
*
*  NOTE: The routine is intended only for debugginig purposes. */

static double err_in_gamma(struct csa *csa)
{     int m = csa->m;
      char *type = csa->type;
      int *head = csa->head;
      double *gamma = csa->gamma;
      double *exact = csa->work4;
      int i;
      double e, emax, temp;
      eval_gamma(csa, exact);
      emax = 0.0;
      for (i = 1; i <= m; i++)
      {  if (type[head[i]] == GLP_FR)
         {  xassert(gamma[i] == 1.0);
            xassert(exact[i] == 1.0);
            continue;
         }
         temp = exact[i];
         e = fabs(temp - gamma[i]) / (1.0 + fabs(temp));
         if (emax < e) emax = e;
      }
      return emax;
}

/***********************************************************************
*  change_basis - change basis header
*
*  This routine changes the basis header to make it corresponding to
*  the adjacent basis. */

static void change_basis(struct csa *csa)
{     int m = csa->m;
#ifdef GLP_DEBUG
      int n = csa->n;
#endif
      char *type = csa->type;
      int *head = csa->head;
#if 1 /* 06/IV-2009 */
      int *bind = csa->bind;
#endif
      char *stat = csa->stat;
      int p = csa->p;
      double delta = csa->delta;
      int q = csa->q;
      int k;
      /* xB[p] leaves the basis, xN[q] enters the basis */
#ifdef GLP_DEBUG
      xassert(1 <= p && p <= m);
      xassert(1 <= q && q <= n);
#endif
      /* xB[p] <-> xN[q] */
      k = head[p], head[p] = head[m+q], head[m+q] = k;
#if 1 /* 06/IV-2009 */
      bind[head[p]] = p, bind[head[m+q]] = m + q;
#endif
      if (type[k] == GLP_FX)
         stat[q] = GLP_NS;
      else if (delta > 0.0)
      {
#ifdef GLP_DEBUG
         xassert(type[k] == GLP_LO || type[k] == GLP_DB);
#endif
         stat[q] = GLP_NL;
      }
      else /* delta < 0.0 */
      {
#ifdef GLP_DEBUG
         xassert(type[k] == GLP_UP || type[k] == GLP_DB);
#endif
         stat[q] = GLP_NU;
      }
      return;
}

/***********************************************************************
*  check_feas - check dual feasibility of basic solution
*
*  If the current basic solution is dual feasible within a tolerance,
*  this routine returns zero, otherwise it returns non-zero. */

static int check_feas(struct csa *csa, double tol_dj)
{     int m = csa->m;
      int n = csa->n;
      char *orig_type = csa->orig_type;
      int *head = csa->head;
      double *cbar = csa->cbar;
      int j, k;
      for (j = 1; j <= n; j++)
      {  k = head[m+j]; /* x[k] = xN[j] */
#ifdef GLP_DEBUG
         xassert(1 <= k && k <= m+n);
#endif
         if (cbar[j] < - tol_dj)
            if (orig_type[k] == GLP_LO || orig_type[k] == GLP_FR)
               return 1;
         if (cbar[j] > + tol_dj)
            if (orig_type[k] == GLP_UP || orig_type[k] == GLP_FR)
               return 1;
      }
      return 0;
}

/***********************************************************************
*  set_aux_bnds - assign auxiliary bounds to variables
*
*  This routine assigns auxiliary bounds to variables to construct an
*  LP problem solved on phase I. */

static void set_aux_bnds(struct csa *csa)
{     int m = csa->m;
      int n = csa->n;
      char *type = csa->type;
      double *lb = csa->lb;
      double *ub = csa->ub;
      char *orig_type = csa->orig_type;
      int *head = csa->head;
      char *stat = csa->stat;
      double *cbar = csa->cbar;
      int j, k;
      for (k = 1; k <= m+n; k++)
      {  switch (orig_type[k])
         {  case GLP_FR:
#if 0
               type[k] = GLP_DB, lb[k] = -1.0, ub[k] = +1.0;
#else
               /* to force free variables to enter the basis */
               type[k] = GLP_DB, lb[k] = -1e3, ub[k] = +1e3;
#endif
               break;
            case GLP_LO:
               type[k] = GLP_DB, lb[k] = 0.0, ub[k] = +1.0;
               break;
            case GLP_UP:
               type[k] = GLP_DB, lb[k] = -1.0, ub[k] = 0.0;
               break;
            case GLP_DB:
            case GLP_FX:
               type[k] = GLP_FX, lb[k] = ub[k] = 0.0;
               break;
            default:
               xassert(orig_type != orig_type);
         }
      }
      for (j = 1; j <= n; j++)
      {  k = head[m+j]; /* x[k] = xN[j] */
#ifdef GLP_DEBUG
         xassert(1 <= k && k <= m+n);
#endif
         if (type[k] == GLP_FX)
            stat[j] = GLP_NS;
         else if (cbar[j] >= 0.0)
            stat[j] = GLP_NL;
         else
            stat[j] = GLP_NU;
      }
      return;
}

/***********************************************************************
*  set_orig_bnds - restore original bounds of variables
*
*  This routine restores original types and bounds of variables and
*  determines statuses of non-basic variables assuming that the current
*  basis is dual feasible. */

static void set_orig_bnds(struct csa *csa)
{     int m = csa->m;
      int n = csa->n;
      char *type = csa->type;
      double *lb = csa->lb;
      double *ub = csa->ub;
      char *orig_type = csa->orig_type;
      double *orig_lb = csa->orig_lb;
      double *orig_ub = csa->orig_ub;
      int *head = csa->head;
      char *stat = csa->stat;
      double *cbar = csa->cbar;
      int j, k;
      memcpy(&type[1], &orig_type[1], (m+n) * sizeof(char));
      memcpy(&lb[1], &orig_lb[1], (m+n) * sizeof(double));
      memcpy(&ub[1], &orig_ub[1], (m+n) * sizeof(double));
      for (j = 1; j <= n; j++)
      {  k = head[m+j]; /* x[k] = xN[j] */
#ifdef GLP_DEBUG
         xassert(1 <= k && k <= m+n);
#endif
         switch (type[k])
         {  case GLP_FR:
               stat[j] = GLP_NF;
               break;
            case GLP_LO:
               stat[j] = GLP_NL;
               break;
            case GLP_UP:
               stat[j] = GLP_NU;
               break;
            case GLP_DB:
               if (cbar[j] >= +DBL_EPSILON)
                  stat[j] = GLP_NL;
               else if (cbar[j] <= -DBL_EPSILON)
                  stat[j] = GLP_NU;
               else if (fabs(lb[k]) <= fabs(ub[k]))
                  stat[j] = GLP_NL;
               else
                  stat[j] = GLP_NU;
               break;
            case GLP_FX:
               stat[j] = GLP_NS;
               break;
            default:
               xassert(type != type);
         }
      }
      return;
}

/***********************************************************************
*  check_stab - check numerical stability of basic solution
*
*  If the current basic solution is dual feasible within a tolerance,
*  this routine returns zero, otherwise it returns non-zero. */

static int check_stab(struct csa *csa, double tol_dj)
{     int n = csa->n;
      char *stat = csa->stat;
      double *cbar = csa->cbar;
      int j;
      for (j = 1; j <= n; j++)
      {  if (cbar[j] < - tol_dj)
            if (stat[j] == GLP_NL || stat[j] == GLP_NF) return 1;
         if (cbar[j] > + tol_dj)
            if (stat[j] == GLP_NU || stat[j] == GLP_NF) return 1;
      }
      return 0;
}

#if 1 /* copied from primal */
/***********************************************************************
*  eval_obj - compute original objective function
*
*  This routine computes the current value of the original objective
*  function. */

static double eval_obj(struct csa *csa)
{     int m = csa->m;
      int n = csa->n;
      double *obj = csa->obj;
      int *head = csa->head;
      double *bbar = csa->bbar;
      int i, j, k;
      double sum;
      sum = obj[0];
      /* walk through the list of basic variables */
      for (i = 1; i <= m; i++)
      {  k = head[i]; /* x[k] = xB[i] */
#ifdef GLP_DEBUG
         xassert(1 <= k && k <= m+n);
#endif
         if (k > m)
            sum += obj[k-m] * bbar[i];
      }
      /* walk through the list of non-basic variables */
      for (j = 1; j <= n; j++)
      {  k = head[m+j]; /* x[k] = xN[j] */
#ifdef GLP_DEBUG
         xassert(1 <= k && k <= m+n);
#endif
         if (k > m)
            sum += obj[k-m] * get_xN(csa, j);
      }
      return sum;
}
#endif

/***********************************************************************
*  display - display the search progress
*
*  This routine displays some information about the search progress. */

static void display(struct csa *csa, const glp_smcp *parm, int spec)
{     int m = csa->m;
      int n = csa->n;
      double *coef = csa->coef;
      char *orig_type = csa->orig_type;
      int *head = csa->head;
      char *stat = csa->stat;
      int phase = csa->phase;
      double *bbar = csa->bbar;
      double *cbar = csa->cbar;
      int i, j, cnt;
      double sum;
      if (parm->msg_lev < GLP_MSG_ON) goto skip;
      if (parm->out_dly > 0 &&
         1000.0 * xdifftime(xtime(), csa->tm_beg) < parm->out_dly)
         goto skip;
      if (csa->it_cnt == csa->it_dpy) goto skip;
      if (!spec && csa->it_cnt % parm->out_frq != 0) goto skip;
      /* compute the sum of dual infeasibilities */
      sum = 0.0;
      if (phase == 1)
      {  for (i = 1; i <= m; i++)
            sum -= coef[head[i]] * bbar[i];
         for (j = 1; j <= n; j++)
            sum -= coef[head[m+j]] * get_xN(csa, j);
      }
      else
      {  for (j = 1; j <= n; j++)
         {  if (cbar[j] < 0.0)
               if (stat[j] == GLP_NL || stat[j] == GLP_NF)
                  sum -= cbar[j];
            if (cbar[j] > 0.0)
               if (stat[j] == GLP_NU || stat[j] == GLP_NF)
                  sum += cbar[j];
         }
      }
      /* determine the number of basic fixed variables */
      cnt = 0;
      for (i = 1; i <= m; i++)
         if (orig_type[head[i]] == GLP_FX) cnt++;
      if (csa->phase == 1)
         xprintf(" %6d: %24s infeas = %10.3e (%d)\n",
            csa->it_cnt, "", sum, cnt);
      else
         xprintf("|%6d: obj = %17.9e  infeas = %10.3e (%d)\n",
            csa->it_cnt, eval_obj(csa), sum, cnt);
      csa->it_dpy = csa->it_cnt;
skip: return;
}

#if 1 /* copied from primal */
/***********************************************************************
*  store_sol - store basic solution back to the problem object
*
*  This routine stores basic solution components back to the problem
*  object. */

static void store_sol(struct csa *csa, glp_prob *lp, int p_stat,
      int d_stat, int ray)
{     int m = csa->m;
      int n = csa->n;
      double zeta = csa->zeta;
      int *head = csa->head;
      char *stat = csa->stat;
      double *bbar = csa->bbar;
      double *cbar = csa->cbar;
      int i, j, k;
#ifdef GLP_DEBUG
      xassert(lp->m == m);
      xassert(lp->n == n);
#endif
      /* basis factorization */
#ifdef GLP_DEBUG
      xassert(!lp->valid && lp->bfd == NULL);
      xassert(csa->valid && csa->bfd != NULL);
#endif
      lp->valid = 1, csa->valid = 0;
      lp->bfd = csa->bfd, csa->bfd = NULL;
      memcpy(&lp->head[1], &head[1], m * sizeof(int));
      /* basic solution status */
      lp->pbs_stat = p_stat;
      lp->dbs_stat = d_stat;
      /* objective function value */
      lp->obj_val = eval_obj(csa);
      /* simplex iteration count */
      lp->it_cnt = csa->it_cnt;
      /* unbounded ray */
      lp->some = ray;
      /* basic variables */
      for (i = 1; i <= m; i++)
      {  k = head[i]; /* x[k] = xB[i] */
#ifdef GLP_DEBUG
         xassert(1 <= k && k <= m+n);
#endif
         if (k <= m)
         {  GLPROW *row = lp->row[k];
            row->stat = GLP_BS;
            row->bind = i;
            row->prim = bbar[i] / row->rii;
            row->dual = 0.0;
         }
         else
         {  GLPCOL *col = lp->col[k-m];
            col->stat = GLP_BS;
            col->bind = i;
            col->prim = bbar[i] * col->sjj;
            col->dual = 0.0;
         }
      }
      /* non-basic variables */
      for (j = 1; j <= n; j++)
      {  k = head[m+j]; /* x[k] = xN[j] */
#ifdef GLP_DEBUG
         xassert(1 <= k && k <= m+n);
#endif
         if (k <= m)
         {  GLPROW *row = lp->row[k];
            row->stat = stat[j];
            row->bind = 0;
#if 0
            row->prim = get_xN(csa, j) / row->rii;
#else
            switch (stat[j])
            {  case GLP_NL:
                  row->prim = row->lb; break;
               case GLP_NU:
                  row->prim = row->ub; break;
               case GLP_NF:
                  row->prim = 0.0; break;
               case GLP_NS:
                  row->prim = row->lb; break;
               default:
                  xassert(stat != stat);
            }
#endif
            row->dual = (cbar[j] * row->rii) / zeta;
         }
         else
         {  GLPCOL *col = lp->col[k-m];
            col->stat = stat[j];
            col->bind = 0;
#if 0
            col->prim = get_xN(csa, j) * col->sjj;
#else
            switch (stat[j])
            {  case GLP_NL:
                  col->prim = col->lb; break;
               case GLP_NU:
                  col->prim = col->ub; break;
               case GLP_NF:
                  col->prim = 0.0; break;
               case GLP_NS:
                  col->prim = col->lb; break;
               default:
                  xassert(stat != stat);
            }
#endif
            col->dual = (cbar[j] / col->sjj) / zeta;
         }
      }
      return;
}
#endif

/***********************************************************************
*  free_csa - deallocate common storage area
*
*  This routine frees all the memory allocated to arrays in the common
*  storage area (CSA). */

static void free_csa(struct csa *csa)
{     xfree(csa->type);
      xfree(csa->lb);
      xfree(csa->ub);
      xfree(csa->coef);
      xfree(csa->orig_type);
      xfree(csa->orig_lb);
      xfree(csa->orig_ub);
      xfree(csa->obj);
      xfree(csa->A_ptr);
      xfree(csa->A_ind);
      xfree(csa->A_val);
#if 1 /* 06/IV-2009 */
      xfree(csa->AT_ptr);
      xfree(csa->AT_ind);
      xfree(csa->AT_val);
#endif
      xfree(csa->head);
#if 1 /* 06/IV-2009 */
      xfree(csa->bind);
#endif
      xfree(csa->stat);
#if 0 /* 06/IV-2009 */
      xfree(csa->N_ptr);
      xfree(csa->N_len);
      xfree(csa->N_ind);
      xfree(csa->N_val);
#endif
      xfree(csa->bbar);
      xfree(csa->cbar);
      xfree(csa->refsp);
      xfree(csa->gamma);
      xfree(csa->trow_ind);
      xfree(csa->trow_vec);
#ifdef GLP_LONG_STEP /* 07/IV-2009 */
      xfree(csa->bkpt);
#endif
      xfree(csa->tcol_ind);
      xfree(csa->tcol_vec);
      xfree(csa->work1);
      xfree(csa->work2);
      xfree(csa->work3);
      xfree(csa->work4);
      xfree(csa);
      return;
}

/***********************************************************************
*  spx_dual - core LP solver based on the dual simplex method
*
*  SYNOPSIS
*
*  #include "glpspx.h"
*  int spx_dual(glp_prob *lp, const glp_smcp *parm);
*
*  DESCRIPTION
*
*  The routine spx_dual is a core LP solver based on the two-phase dual
*  simplex method.
*
*  RETURNS
*
*  0  LP instance has been successfully solved.
*
*  GLP_EOBJLL
*     Objective lower limit has been reached (maximization).
*
*  GLP_EOBJUL
*     Objective upper limit has been reached (minimization).
*
*  GLP_EITLIM
*     Iteration limit has been exhausted.
*
*  GLP_ETMLIM
*     Time limit has been exhausted.
*
*  GLP_EFAIL
*     The solver failed to solve LP instance. */

int spx_dual(glp_prob *lp, const glp_smcp *parm)
{     struct csa *csa;
      int binv_st = 2;
      /* status of basis matrix factorization:
         0 - invalid; 1 - just computed; 2 - updated */
      int bbar_st = 0;
      /* status of primal values of basic variables:
         0 - invalid; 1 - just computed; 2 - updated */
      int cbar_st = 0;
      /* status of reduced costs of non-basic variables:
         0 - invalid; 1 - just computed; 2 - updated */
      int rigorous = 0;
      /* rigorous mode flag; this flag is used to enable iterative
         refinement on computing pivot rows and columns of the simplex
         table */
      int check = 0;
      int p_stat, d_stat, ret;
      /* allocate and initialize the common storage area */
      csa = alloc_csa(lp);
      init_csa(csa, lp);
      if (parm->msg_lev >= GLP_MSG_DBG)
         xprintf("Objective scale factor = %g\n", csa->zeta);
loop: /* main loop starts here */
      /* compute factorization of the basis matrix */
      if (binv_st == 0)
      {  ret = invert_B(csa);
         if (ret != 0)
         {  if (parm->msg_lev >= GLP_MSG_ERR)
            {  xprintf("Error: unable to factorize the basis matrix (%d"
                  ")\n", ret);
               xprintf("Sorry, basis recovery procedure not implemented"
                  " yet\n");
            }
            xassert(!lp->valid && lp->bfd == NULL);
            lp->bfd = csa->bfd, csa->bfd = NULL;
            lp->pbs_stat = lp->dbs_stat = GLP_UNDEF;
            lp->obj_val = 0.0;
            lp->it_cnt = csa->it_cnt;
            lp->some = 0;
            ret = GLP_EFAIL;
            goto done;
         }
         csa->valid = 1;
         binv_st = 1; /* just computed */
         /* invalidate basic solution components */
         bbar_st = cbar_st = 0;
      }
      /* compute reduced costs of non-basic variables */
      if (cbar_st == 0)
      {  eval_cbar(csa);
         cbar_st = 1; /* just computed */
         /* determine the search phase, if not determined yet */
         if (csa->phase == 0)
         {  if (check_feas(csa, 0.90 * parm->tol_dj) != 0)
            {  /* current basic solution is dual infeasible */
               /* start searching for dual feasible solution */
               csa->phase = 1;
               set_aux_bnds(csa);
            }
            else
            {  /* current basic solution is dual feasible */
               /* start searching for optimal solution */
               csa->phase = 2;
               set_orig_bnds(csa);
            }
            xassert(check_stab(csa, parm->tol_dj) == 0);
            /* some non-basic double-bounded variables might become
               fixed (on phase I) or vice versa (on phase II) */
#if 0 /* 06/IV-2009 */
            build_N(csa);
#endif
            csa->refct = 0;
            /* bounds of non-basic variables have been changed, so
               invalidate primal values */
            bbar_st = 0;
         }
         /* make sure that the current basic solution remains dual
            feasible */
         if (check_stab(csa, parm->tol_dj) != 0)
         {  if (parm->msg_lev >= GLP_MSG_ERR)
               xprintf("Warning: numerical instability (dual simplex, p"
                  "hase %s)\n", csa->phase == 1 ? "I" : "II");
#if 1
            if (parm->meth == GLP_DUALP)
            {  store_sol(csa, lp, GLP_UNDEF, GLP_UNDEF, 0);
               ret = GLP_EFAIL;
               goto done;
            }
#endif
            /* restart the search */
            csa->phase = 0;
            binv_st = 0;
            rigorous = 5;
            goto loop;
         }
      }
      xassert(csa->phase == 1 || csa->phase == 2);
      /* on phase I we do not need to wait until the current basic
         solution becomes primal feasible; it is sufficient to make
         sure that all reduced costs have correct signs */
      if (csa->phase == 1 && check_feas(csa, parm->tol_dj) == 0)
      {  /* the current basis is dual feasible; switch to phase II */
         display(csa, parm, 1);
         csa->phase = 2;
         if (cbar_st != 1)
         {  eval_cbar(csa);
            cbar_st = 1;
         }
         set_orig_bnds(csa);
#if 0 /* 06/IV-2009 */
         build_N(csa);
#endif
         csa->refct = 0;
         bbar_st = 0;
      }
      /* compute primal values of basic variables */
      if (bbar_st == 0)
      {  eval_bbar(csa);
         if (csa->phase == 2)
            csa->bbar[0] = eval_obj(csa);
         bbar_st = 1; /* just computed */
      }
      /* redefine the reference space, if required */
      switch (parm->pricing)
      {  case GLP_PT_STD:
            break;
         case GLP_PT_PSE:
            if (csa->refct == 0) reset_refsp(csa);
            break;
         default:
            xassert(parm != parm);
      }
      /* at this point the basis factorization and all basic solution
         components are valid */
      xassert(binv_st && bbar_st && cbar_st);
      /* check accuracy of current basic solution components (only for
         debugging) */
      if (check)
      {  double e_bbar = err_in_bbar(csa);
         double e_cbar = err_in_cbar(csa);
         double e_gamma =
            (parm->pricing == GLP_PT_PSE ? err_in_gamma(csa) : 0.0);
         xprintf("e_bbar = %10.3e; e_cbar = %10.3e; e_gamma = %10.3e\n",
            e_bbar, e_cbar, e_gamma);
         xassert(e_bbar <= 1e-5 && e_cbar <= 1e-5 && e_gamma <= 1e-3);
      }
      /* if the objective has to be maximized, check if it has reached
         its lower limit */
      if (csa->phase == 2 && csa->zeta < 0.0 &&
          parm->obj_ll > -DBL_MAX && csa->bbar[0] <= parm->obj_ll)
      {  if (bbar_st != 1 || cbar_st != 1)
         {  if (bbar_st != 1) bbar_st = 0;
            if (cbar_st != 1) cbar_st = 0;
            goto loop;
         }
         display(csa, parm, 1);
         if (parm->msg_lev >= GLP_MSG_ALL)
            xprintf("OBJECTIVE LOWER LIMIT REACHED; SEARCH TERMINATED\n"
               );
         store_sol(csa, lp, GLP_INFEAS, GLP_FEAS, 0);
         ret = GLP_EOBJLL;
         goto done;
      }
      /* if the objective has to be minimized, check if it has reached
         its upper limit */
      if (csa->phase == 2 && csa->zeta > 0.0 &&
          parm->obj_ul < +DBL_MAX && csa->bbar[0] >= parm->obj_ul)
      {  if (bbar_st != 1 || cbar_st != 1)
         {  if (bbar_st != 1) bbar_st = 0;
            if (cbar_st != 1) cbar_st = 0;
            goto loop;
         }
         display(csa, parm, 1);
         if (parm->msg_lev >= GLP_MSG_ALL)
            xprintf("OBJECTIVE UPPER LIMIT REACHED; SEARCH TERMINATED\n"
               );
         store_sol(csa, lp, GLP_INFEAS, GLP_FEAS, 0);
         ret = GLP_EOBJUL;
         goto done;
      }
      /* check if the iteration limit has been exhausted */
      if (parm->it_lim < INT_MAX &&
          csa->it_cnt - csa->it_beg >= parm->it_lim)
      {  if (csa->phase == 2 && bbar_st != 1 || cbar_st != 1)
         {  if (csa->phase == 2 && bbar_st != 1) bbar_st = 0;
            if (cbar_st != 1) cbar_st = 0;
            goto loop;
         }
         display(csa, parm, 1);
         if (parm->msg_lev >= GLP_MSG_ALL)
            xprintf("ITERATION LIMIT EXCEEDED; SEARCH TERMINATED\n");
         switch (csa->phase)
         {  case 1:
               d_stat = GLP_INFEAS;
               set_orig_bnds(csa);
               eval_bbar(csa);
               break;
            case 2:
               d_stat = GLP_FEAS;
               break;
            default:
               xassert(csa != csa);
         }
         store_sol(csa, lp, GLP_INFEAS, d_stat, 0);
         ret = GLP_EITLIM;
         goto done;
      }
      /* check if the time limit has been exhausted */
      if (parm->tm_lim < INT_MAX &&
          1000.0 * xdifftime(xtime(), csa->tm_beg) >= parm->tm_lim)
      {  if (csa->phase == 2 && bbar_st != 1 || cbar_st != 1)
         {  if (csa->phase == 2 && bbar_st != 1) bbar_st = 0;
            if (cbar_st != 1) cbar_st = 0;
            goto loop;
         }
         display(csa, parm, 1);
         if (parm->msg_lev >= GLP_MSG_ALL)
            xprintf("TIME LIMIT EXCEEDED; SEARCH TERMINATED\n");
         switch (csa->phase)
         {  case 1:
               d_stat = GLP_INFEAS;
               set_orig_bnds(csa);
               eval_bbar(csa);
               break;
            case 2:
               d_stat = GLP_FEAS;
               break;
            default:
               xassert(csa != csa);
         }
         store_sol(csa, lp, GLP_INFEAS, d_stat, 0);
         ret = GLP_ETMLIM;
         goto done;
      }
      /* display the search progress */
      display(csa, parm, 0);
      /* choose basic variable xB[p] */
      chuzr(csa, parm->tol_bnd);
      if (csa->p == 0)
      {  if (bbar_st != 1 || cbar_st != 1)
         {  if (bbar_st != 1) bbar_st = 0;
            if (cbar_st != 1) cbar_st = 0;
            goto loop;
         }
         display(csa, parm, 1);
         switch (csa->phase)
         {  case 1:
               if (parm->msg_lev >= GLP_MSG_ALL)
                  xprintf("PROBLEM HAS NO DUAL FEASIBLE SOLUTION\n");
               set_orig_bnds(csa);
               eval_bbar(csa);
               p_stat = GLP_INFEAS, d_stat = GLP_NOFEAS;
               break;
            case 2:
               if (parm->msg_lev >= GLP_MSG_ALL)
                  xprintf("OPTIMAL SOLUTION FOUND\n");
               p_stat = d_stat = GLP_FEAS;
               break;
            default:
               xassert(csa != csa);
         }
         store_sol(csa, lp, p_stat, d_stat, 0);
         ret = 0;
         goto done;
      }
      /* compute pivot row of the simplex table */
      {  double *rho = csa->work4;
         eval_rho(csa, rho);
         if (rigorous) refine_rho(csa, rho);
         eval_trow(csa, rho);
         sort_trow(csa, parm->tol_bnd);
      }
      /* unlike primal simplex there is no need to check accuracy of
         the primal value of xB[p] (which might be computed using the
         pivot row), since bbar is a result of FTRAN */
#ifdef GLP_LONG_STEP /* 07/IV-2009 */
      long_step(csa);
      if (csa->nbps > 0)
      {  csa->q = csa->bkpt[csa->nbps].j;
         if (csa->delta > 0.0)
            csa->new_dq = + csa->bkpt[csa->nbps].t;
         else
            csa->new_dq = - csa->bkpt[csa->nbps].t;
      }
      else
#endif
      /* choose non-basic variable xN[q] */
      switch (parm->r_test)
      {  case GLP_RT_STD:
            chuzc(csa, 0.0);
            break;
         case GLP_RT_HAR:
            chuzc(csa, 0.30 * parm->tol_dj);
            break;
         default:
            xassert(parm != parm);
      }
      if (csa->q == 0)
      {  if (bbar_st != 1 || cbar_st != 1 || !rigorous)
         {  if (bbar_st != 1) bbar_st = 0;
            if (cbar_st != 1) cbar_st = 0;
            rigorous = 1;
            goto loop;
         }
         display(csa, parm, 1);
         switch (csa->phase)
         {  case 1:
               if (parm->msg_lev >= GLP_MSG_ERR)
                  xprintf("Error: unable to choose basic variable on ph"
                     "ase I\n");
               xassert(!lp->valid && lp->bfd == NULL);
               lp->bfd = csa->bfd, csa->bfd = NULL;
               lp->pbs_stat = lp->dbs_stat = GLP_UNDEF;
               lp->obj_val = 0.0;
               lp->it_cnt = csa->it_cnt;
               lp->some = 0;
               ret = GLP_EFAIL;
               break;
            case 2:
               if (parm->msg_lev >= GLP_MSG_ALL)
                  xprintf("PROBLEM HAS NO FEASIBLE SOLUTION\n");
               store_sol(csa, lp, GLP_NOFEAS, GLP_FEAS,
                  csa->head[csa->p]);
               ret = 0;
               break;
            default:
               xassert(csa != csa);
         }
         goto done;
      }
      /* check if the pivot element is acceptable */
      {  double piv = csa->trow_vec[csa->q];
         double eps = 1e-5 * (1.0 + 0.01 * csa->trow_max);
         if (fabs(piv) < eps)
         {  if (parm->msg_lev >= GLP_MSG_DBG)
               xprintf("piv = %.12g; eps = %g\n", piv, eps);
            if (!rigorous)
            {  rigorous = 5;
               goto loop;
            }
         }
      }
      /* now xN[q] and xB[p] have been chosen anyhow */
      /* compute pivot column of the simplex table */
      eval_tcol(csa);
      if (rigorous) refine_tcol(csa);
      /* accuracy check based on the pivot element */
      {  double piv1 = csa->tcol_vec[csa->p]; /* more accurate */
         double piv2 = csa->trow_vec[csa->q]; /* less accurate */
         xassert(piv1 != 0.0);
         if (fabs(piv1 - piv2) > 1e-8 * (1.0 + fabs(piv1)) ||
             !(piv1 > 0.0 && piv2 > 0.0 || piv1 < 0.0 && piv2 < 0.0))
         {  if (parm->msg_lev >= GLP_MSG_DBG)
               xprintf("piv1 = %.12g; piv2 = %.12g\n", piv1, piv2);
            if (binv_st != 1 || !rigorous)
            {  if (binv_st != 1) binv_st = 0;
               rigorous = 5;
               goto loop;
            }
            /* (not a good idea; should be revised later) */
            if (csa->tcol_vec[csa->p] == 0.0)
            {  csa->tcol_nnz++;
               xassert(csa->tcol_nnz <= csa->m);
               csa->tcol_ind[csa->tcol_nnz] = csa->p;
            }
            csa->tcol_vec[csa->p] = piv2;
         }
      }
      /* update primal values of basic variables */
#ifdef GLP_LONG_STEP /* 07/IV-2009 */
      if (csa->nbps > 0)
      {  int kk, j, k;
         for (kk = 1; kk < csa->nbps; kk++)
         {  if (csa->bkpt[kk].t >= csa->bkpt[csa->nbps].t) continue;
            j = csa->bkpt[kk].j;
            k = csa->head[csa->m + j];
            xassert(csa->type[k] == GLP_DB);
            if (csa->stat[j] == GLP_NL)
               csa->stat[j] = GLP_NU;
            else
               csa->stat[j] = GLP_NL;
         }
      }
      bbar_st = 0;
#else
      update_bbar(csa);
      if (csa->phase == 2)
         csa->bbar[0] += (csa->cbar[csa->q] / csa->zeta) *
            (csa->delta / csa->tcol_vec[csa->p]);
      bbar_st = 2; /* updated */
#endif
      /* update reduced costs of non-basic variables */
      update_cbar(csa);
      cbar_st = 2; /* updated */
      /* update steepest edge coefficients */
      switch (parm->pricing)
      {  case GLP_PT_STD:
            break;
         case GLP_PT_PSE:
            if (csa->refct > 0) update_gamma(csa);
            break;
         default:
            xassert(parm != parm);
      }
      /* update factorization of the basis matrix */
      ret = update_B(csa, csa->p, csa->head[csa->m+csa->q]);
      if (ret == 0)
         binv_st = 2; /* updated */
      else
      {  csa->valid = 0;
         binv_st = 0; /* invalid */
      }
#if 0 /* 06/IV-2009 */
      /* update matrix N */
      del_N_col(csa, csa->q, csa->head[csa->m+csa->q]);
      if (csa->type[csa->head[csa->p]] != GLP_FX)
         add_N_col(csa, csa->q, csa->head[csa->p]);
#endif
      /* change the basis header */
      change_basis(csa);
      /* iteration complete */
      csa->it_cnt++;
      if (rigorous > 0) rigorous--;
      goto loop;
done: /* deallocate the common storage area */
      free_csa(csa);
      /* return to the calling program */
      return ret;
}

/* eof */
