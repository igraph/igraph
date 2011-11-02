/* glpspx01.c (primal simplex method) */

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

#include "glpspx.h"

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
         variable x[k] (note that on phase I auxiliary variables also
         may have non-zero objective coefficients) */
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
      /*--------------------------------------------------------------*/
      /* working parameters */
      int phase;
      /* search phase:
         0 - not determined yet
         1 - search for primal feasible solution
         2 - search for optimal solution */
      glp_long tm_beg;
      /* time value at the beginning of the search */
      int it_beg;
      /* simplex iteration count at the beginning of the search */
      int it_cnt;
      /* simplex iteration count; it increases by one every time the
         basis changes (including the case when a non-basic variable
         jumps to its opposite bound) */
      int it_dpy;
      /* simplex iteration count at the most recent display output */
      /*--------------------------------------------------------------*/
      /* basic solution components */
      double *bbar; /* double bbar[1+m]; */
      /* bbar[0] is not used;
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
      double *gamma; /* double gamma[1+n]; */
      /* gamma[0] is not used;
         gamma[j], 1 <= j <= n, is the steepest edge coefficient for
         non-basic variable xN[j]; if xN[j] is fixed, gamma[j] is not
         used and just set to 1 */
      /*--------------------------------------------------------------*/
      /* non-basic variable xN[q] chosen to enter the basis */
      int q;
      /* index of the non-basic variable xN[q] chosen, 1 <= q <= n;
         if the set of eligible non-basic variables is empty and thus
         no variable has been chosen, q is set to 0 */
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
      double tcol_max;
      /* infinity (maximum) norm of the column (max |tcol_vec[i]|) */
      int tcol_num;
      /* number of significant non-zero components, which means that:
         |tcol_vec[i]| >= eps for i in tcol_ind[1,...,num],
         |tcol_vec[i]| <  eps for i in tcol_ind[num+1,...,nnz],
         where eps is a pivot tolerance */
      /*--------------------------------------------------------------*/
      /* basic variable xB[p] chosen to leave the basis */
      int p;
      /* index of the basic variable xB[p] chosen, 1 <= p <= m;
         p = 0 means that no basic variable reaches its bound;
         p < 0 means that non-basic variable xN[q] reaches its opposite
         bound before any basic variable */
      int p_stat;
      /* new status (GLP_NL, GLP_NU, or GLP_NS) to be assigned to xB[p]
         once it has left the basis */
      double teta;
      /* change of non-basic variable xN[q] (see above), on which xB[p]
         (or, if p < 0, xN[q] itself) reaches its bound */
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
      csa->obj = xcalloc(1+n, sizeof(double));
      csa->A_ptr = xcalloc(1+n+1, sizeof(int));
      csa->A_ind = xcalloc(1+nnz, sizeof(int));
      csa->A_val = xcalloc(1+nnz, sizeof(double));
      csa->head = xcalloc(1+m+n, sizeof(int));
      csa->stat = xcalloc(1+n, sizeof(char));
      csa->N_ptr = xcalloc(1+m+1, sizeof(int));
      csa->N_len = xcalloc(1+m, sizeof(int));
      csa->N_ind = NULL; /* will be allocated later */
      csa->N_val = NULL; /* will be allocated later */
      csa->bbar = xcalloc(1+m, sizeof(double));
      csa->cbar = xcalloc(1+n, sizeof(double));
      csa->refsp = xcalloc(1+m+n, sizeof(char));
      csa->gamma = xcalloc(1+n, sizeof(double));
      csa->tcol_ind = xcalloc(1+m, sizeof(int));
      csa->tcol_vec = xcalloc(1+m, sizeof(double));
      csa->trow_ind = xcalloc(1+n, sizeof(int));
      csa->trow_vec = xcalloc(1+n, sizeof(double));
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

static void alloc_N(struct csa *csa);
static void build_N(struct csa *csa);

static void init_csa(struct csa *csa, glp_prob *lp)
{     int m = csa->m;
      int n = csa->n;
      char *type = csa->type;
      double *lb = csa->lb;
      double *ub = csa->ub;
      double *coef = csa->coef;
      double *obj = csa->obj;
      int *A_ptr = csa->A_ptr;
      int *A_ind = csa->A_ind;
      double *A_val = csa->A_val;
      int *head = csa->head;
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
      xassert(loc == lp->nnz+1);
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
      /* factorization of matrix B */
      csa->valid = 1, lp->valid = 0;
      csa->bfd = lp->bfd, lp->bfd = NULL;
      /* matrix N (by rows) */
      alloc_N(csa);
      build_N(csa);
      /* working parameters */
      csa->phase = 0;
      csa->tm_beg = xtime();
      csa->it_beg = csa->it_cnt = lp->it_cnt;
      csa->it_dpy = -1;
      /* reference space and steepest edge coefficients */
      csa->refct = 0;
      memset(&refsp[1], 0, (m+n) * sizeof(char));
      for (j = 1; j <= n; j++) gamma[j] = 1.0;
      return;
}

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

/***********************************************************************
*  alloc_N - allocate matrix N
*
*  This routine determines maximal row lengths of matrix N, sets its
*  row pointers, and then allocates arrays N_ind and N_val.
*
*  Note that some fixed structural variables may temporarily become
*  double-bounded, so corresponding columns of matrix A should not be
*  ignored on calculating maximal row lengths of matrix N. */

static void alloc_N(struct csa *csa)
{     int m = csa->m;
      int n = csa->n;
      int *A_ptr = csa->A_ptr;
      int *A_ind = csa->A_ind;
      int *N_ptr = csa->N_ptr;
      int *N_len = csa->N_len;
      int i, j, beg, end, ptr;
      /* determine number of non-zeros in each row of the augmented
         constraint matrix (I|-A) */
      for (i = 1; i <= m; i++)
         N_len[i] = 1;
      for (j = 1; j <= n; j++)
      {  beg = A_ptr[j];
         end = A_ptr[j+1];
         for (ptr = beg; ptr < end; ptr++)
            N_len[A_ind[ptr]]++;
      }
      /* determine maximal row lengths of matrix N and set its row
         pointers */
      N_ptr[1] = 1;
      for (i = 1; i <= m; i++)
      {  /* row of matrix N cannot have more than n non-zeros */
         if (N_len[i] > n) N_len[i] = n;
         N_ptr[i+1] = N_ptr[i] + N_len[i];
      }
      /* now maximal number of non-zeros in matrix N is known */
      csa->N_ind = xcalloc(N_ptr[m+1], sizeof(int));
      csa->N_val = xcalloc(N_ptr[m+1], sizeof(double));
      return;
}

/***********************************************************************
*  add_N_col - add column of matrix (I|-A) to matrix N
*
*  This routine adds j-th column to matrix N which is k-th column of
*  the augmented constraint matrix (I|-A). (It is assumed that old j-th
*  column was previously removed from matrix N.) */

static void add_N_col(struct csa *csa, int j, int k)
{     int m = csa->m;
#ifdef GLP_DEBUG
      int n = csa->n;
#endif
      int *N_ptr = csa->N_ptr;
      int *N_len = csa->N_len;
      int *N_ind = csa->N_ind;
      double *N_val = csa->N_val;
      int pos;
#ifdef GLP_DEBUG
      xassert(1 <= j && j <= n);
      xassert(1 <= k && k <= m+n);
#endif
      if (k <= m)
      {  /* N[j] is k-th column of submatrix I */
         pos = N_ptr[k] + (N_len[k]++);
#ifdef GLP_DEBUG
         xassert(pos < N_ptr[k+1]);
#endif
         N_ind[pos] = j;
         N_val[pos] = 1.0;
      }
      else
      {  /* N[j] is (k-m)-th column of submatrix (-A) */
         int *A_ptr = csa->A_ptr;
         int *A_ind = csa->A_ind;
         double *A_val = csa->A_val;
         int i, beg, end, ptr;
         beg = A_ptr[k-m];
         end = A_ptr[k-m+1];
         for (ptr = beg; ptr < end; ptr++)
         {  i = A_ind[ptr]; /* row number */
            pos = N_ptr[i] + (N_len[i]++);
#ifdef GLP_DEBUG
            xassert(pos < N_ptr[i+1]);
#endif
            N_ind[pos] = j;
            N_val[pos] = - A_val[ptr];
         }
      }
      return;
}

/***********************************************************************
*  del_N_col - remove column of matrix (I|-A) from matrix N
*
*  This routine removes j-th column from matrix N which is k-th column
*  of the augmented constraint matrix (I|-A). */

static void del_N_col(struct csa *csa, int j, int k)
{     int m = csa->m;
#ifdef GLP_DEBUG
      int n = csa->n;
#endif
      int *N_ptr = csa->N_ptr;
      int *N_len = csa->N_len;
      int *N_ind = csa->N_ind;
      double *N_val = csa->N_val;
      int pos, head, tail;
#ifdef GLP_DEBUG
      xassert(1 <= j && j <= n);
      xassert(1 <= k && k <= m+n);
#endif
      if (k <= m)
      {  /* N[j] is k-th column of submatrix I */
         /* find element in k-th row of N */
         head = N_ptr[k];
         for (pos = head; N_ind[pos] != j; pos++) /* nop */;
         /* and remove it from the row list */
         tail = head + (--N_len[k]);
#ifdef GLP_DEBUG
         xassert(pos <= tail);
#endif
         N_ind[pos] = N_ind[tail];
         N_val[pos] = N_val[tail];
      }
      else
      {  /* N[j] is (k-m)-th column of submatrix (-A) */
         int *A_ptr = csa->A_ptr;
         int *A_ind = csa->A_ind;
         int i, beg, end, ptr;
         beg = A_ptr[k-m];
         end = A_ptr[k-m+1];
         for (ptr = beg; ptr < end; ptr++)
         {  i = A_ind[ptr]; /* row number */
            /* find element in i-th row of N */
            head = N_ptr[i];
            for (pos = head; N_ind[pos] != j; pos++) /* nop */;
            /* and remove it from the row list */
            tail = head + (--N_len[i]);
#ifdef GLP_DEBUG
            xassert(pos <= tail);
#endif
            N_ind[pos] = N_ind[tail];
            N_val[pos] = N_val[tail];
         }
      }
      return;
}

/***********************************************************************
*  build_N - build matrix N for current basis
*
*  This routine builds matrix N for the current basis from columns
*  of the augmented constraint matrix (I|-A) corresponding to non-basic
*  non-fixed variables. */

static void build_N(struct csa *csa)
{     int m = csa->m;
      int n = csa->n;
      int *head = csa->head;
      char *stat = csa->stat;
      int *N_len = csa->N_len;
      int j, k;
      /* N := empty matrix */
      memset(&N_len[1], 0, m * sizeof(int));
      /* go through non-basic columns of matrix (I|-A) */
      for (j = 1; j <= n; j++)
      {  if (stat[j] != GLP_NS)
         {  /* xN[j] is non-fixed; add j-th column to matrix N which is
               k-th column of matrix (I|-A) */
            k = head[m+j]; /* x[k] = xN[j] */
#ifdef GLP_DEBUG
            xassert(1 <= k && k <= m+n);
#endif
            add_N_col(csa, j, k);
         }
      }
      return;
}

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

/***********************************************************************
*  eval_bbar - compute and store primal values of basic variables
*
*  This routine computes primal values of all basic variables and then
*  stores them in the solution array. */

static void eval_bbar(struct csa *csa)
{     eval_beta(csa, csa->bbar);
      return;
}

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
      int j, k;
      xassert(csa->refct == 0);
      csa->refct = 1000;
      memset(&refsp[1], 0, (m+n) * sizeof(char));
      for (j = 1; j <= n; j++)
      {  k = head[m+j]; /* x[k] = xN[j] */
         refsp[k] = 1;
         gamma[j] = 1.0;
      }
      return;
}

/***********************************************************************
*  eval_gamma - compute steepest edge coefficient
*
*  This routine computes the steepest edge coefficient for non-basic
*  variable xN[j] using its direct definition:
*
*     gamma[j] = delta[j] +  sum   alfa[i,j]^2,
*                           i in R
*
*  where delta[j] = 1, if xN[j] is in the current reference space,
*  and 0 otherwise; R is a set of basic variables xB[i], which are in
*  the current reference space; alfa[i,j] are elements of the current
*  simplex table.
*
*  NOTE: The routine is intended only for debugginig purposes. */

static double eval_gamma(struct csa *csa, int j)
{     int m = csa->m;
#ifdef GLP_DEBUG
      int n = csa->n;
#endif
      int *head = csa->head;
      char *refsp = csa->refsp;
      double *alfa = csa->work3;
      double *h = csa->work3;
      int i, k;
      double gamma;
#ifdef GLP_DEBUG
      xassert(1 <= j && j <= n);
#endif
      k = head[m+j]; /* x[k] = xN[j] */
#ifdef GLP_DEBUG
      xassert(1 <= k && k <= m+n);
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
      /* compute gamma */
      gamma = (refsp[k] ? 1.0 : 0.0);
      for (i = 1; i <= m; i++)
      {  k = head[i];
#ifdef GLP_DEBUG
         xassert(1 <= k && k <= m+n);
#endif
         if (refsp[k]) gamma += alfa[i] * alfa[i];
      }
      return gamma;
}

/***********************************************************************
*  chuzc - choose non-basic variable (column of the simplex table)
*
*  This routine chooses non-basic variable xN[q], which has largest
*  weighted reduced cost:
*
*     |d[q]| / sqrt(gamma[q]) = max  |d[j]| / sqrt(gamma[j]),
*                              j in J
*
*  where J is a subset of eligible non-basic variables xN[j], d[j] is
*  reduced cost of xN[j], gamma[j] is the steepest edge coefficient.
*
*  The working objective function is always minimized, so the sign of
*  d[q] determines direction, in which xN[q] has to change:
*
*     if d[q] < 0, xN[q] has to increase;
*
*     if d[q] > 0, xN[q] has to decrease.
*
*  If |d[j]| <= tol_dj, where tol_dj is a specified tolerance, xN[j]
*  is not included in J and therefore ignored. (It is assumed that the
*  working objective row is appropriately scaled, i.e. max|c[k]| = 1.)
*
*  If J is empty and no variable has been chosen, q is set to 0. */

static void chuzc(struct csa *csa, double tol_dj)
{     int n = csa->n;
      char *stat = csa->stat;
      double *cbar = csa->cbar;
      double *gamma = csa->gamma;
      int j, q;
      double dj, best, temp;
      /* nothing is chosen so far */
      q = 0, best = 0.0;
      /* look through the list of non-basic variables */
      for (j = 1; j <= n; j++)
      {  dj = cbar[j];
         switch (stat[j])
         {  case GLP_NL:
               /* xN[j] can increase */
               if (dj >= - tol_dj) continue;
               break;
            case GLP_NU:
               /* xN[j] can decrease */
               if (dj <= + tol_dj) continue;
               break;
            case GLP_NF:
               /* xN[j] can change in any direction */
               if (- tol_dj <= dj && dj <= + tol_dj) continue;
               break;
            case GLP_NS:
               /* xN[j] cannot change at all */
               continue;
            default:
               xassert(stat != stat);
         }
         /* xN[j] is eligible non-basic variable; choose one which has
            largest weighted reduced cost */
#ifdef GLP_DEBUG
         xassert(gamma[j] > 0.0);
#endif
         temp = (dj * dj) / gamma[j];
         if (best < temp)
            q = j, best = temp;
      }
      /* store the index of non-basic variable xN[q] chosen */
      csa->q = q;
      return;
}

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

/***********************************************************************
*  sort_tcol - sort pivot column of the simplex table
*
*  This routine reorders the list of non-zero elements of the pivot
*  column to put significant elements, whose magnitude is not less than
*  a specified tolerance, in front of the list, and stores the number
*  of significant elements in tcol_num. */

static void sort_tcol(struct csa *csa, double tol_piv)
{
#ifdef GLP_DEBUG
      int m = csa->m;
#endif
      int nnz = csa->tcol_nnz;
      int *tcol_ind = csa->tcol_ind;
      double *tcol_vec = csa->tcol_vec;
      int i, num, pos;
      double big, eps, temp;
      /* compute infinity (maximum) norm of the column */
      big = 0.0;
      for (pos = 1; pos <= nnz; pos++)
      {
#ifdef GLP_DEBUG
         i = tcol_ind[pos];
         xassert(1 <= i && i <= m);
#endif
         temp = fabs(tcol_vec[tcol_ind[pos]]);
         if (big < temp) big = temp;
      }
      csa->tcol_max = big;
      /* determine absolute pivot tolerance */
      eps = tol_piv * (1.0 + 0.01 * big);
      /* move significant column components to front of the list */
      for (num = 0; num < nnz; )
      {  i = tcol_ind[nnz];
         if (fabs(tcol_vec[i]) < eps)
            nnz--;
         else
         {  num++;
            tcol_ind[nnz] = tcol_ind[num];
            tcol_ind[num] = i;
         }
      }
      csa->tcol_num = num;
      return;
}

/***********************************************************************
*  chuzr - choose basic variable (row of the simplex table)
*
*  This routine chooses basic variable xB[p], which reaches its bound
*  first on changing non-basic variable xN[q] in valid direction.
*
*  The parameter rtol is a relative tolerance used to relax bounds of
*  basic variables. If rtol = 0, the routine implements the standard
*  ratio test. Otherwise, if rtol > 0, the routine implements Harris'
*  two-pass ratio test. In the latter case rtol should be about three
*  times less than a tolerance used to check primal feasibility. */

static void chuzr(struct csa *csa, double rtol)
{     int m = csa->m;
#ifdef GLP_DEBUG
      int n = csa->n;
#endif
      char *type = csa->type;
      double *lb = csa->lb;
      double *ub = csa->ub;
      double *coef = csa->coef;
      int *head = csa->head;
      int phase = csa->phase;
      double *bbar = csa->bbar;
      double *cbar = csa->cbar;
      int q = csa->q;
      int *tcol_ind = csa->tcol_ind;
      double *tcol_vec = csa->tcol_vec;
      int tcol_num = csa->tcol_num;
      int i, i_stat, k, p, p_stat, pos;
      double alfa, big, delta, s, t, teta, tmax;
#ifdef GLP_DEBUG
      xassert(1 <= q && q <= n);
#endif
      /* s := - sign(d[q]), where d[q] is reduced cost of xN[q] */
#ifdef GLP_DEBUG
      xassert(cbar[q] != 0.0);
#endif
      s = (cbar[q] > 0.0 ? -1.0 : +1.0);
      /*** FIRST PASS ***/
      k = head[m+q]; /* x[k] = xN[q] */
#ifdef GLP_DEBUG
      xassert(1 <= k && k <= m+n);
#endif
      if (type[k] == GLP_DB)
      {  /* xN[q] has both lower and upper bounds */
         p = -1, p_stat = 0, teta = ub[k] - lb[k], big = 1.0;
      }
      else
      {  /* xN[q] has no opposite bound */
         p = 0, p_stat = 0, teta = DBL_MAX, big = 0.0;
      }
      /* walk through significant elements of the pivot column */
      for (pos = 1; pos <= tcol_num; pos++)
      {  i = tcol_ind[pos];
#ifdef GLP_DEBUG
         xassert(1 <= i && i <= m);
#endif
         k = head[i]; /* x[k] = xB[i] */
#ifdef GLP_DEBUG
         xassert(1 <= k && k <= m+n);
#endif
         alfa = s * tcol_vec[i];
#ifdef GLP_DEBUG
         xassert(alfa != 0.0);
#endif
         /* xB[i] = ... + alfa * xN[q] + ..., and due to s we need to
            consider the only case when xN[q] is increasing */
         if (alfa > 0.0)
         {  /* xB[i] is increasing */
            if (phase == 1 && coef[k] < 0.0)
            {  /* xB[i] violates its lower bound, which plays the role
                  of an upper bound on phase I */
               delta = rtol * (1.0 + kappa * fabs(lb[k]));
               t = ((lb[k] + delta) - bbar[i]) / alfa;
               i_stat = GLP_NL;
            }
            else if (phase == 1 && coef[k] > 0.0)
            {  /* xB[i] violates its upper bound, which plays the role
                  of an lower bound on phase I */
               continue;
            }
            else if (type[k] == GLP_UP || type[k] == GLP_DB ||
                     type[k] == GLP_FX)
            {  /* xB[i] is within its bounds and has an upper bound */
               delta = rtol * (1.0 + kappa * fabs(ub[k]));
               t = ((ub[k] + delta) - bbar[i]) / alfa;
               i_stat = GLP_NU;
            }
            else
            {  /* xB[i] is within its bounds and has no upper bound */
               continue;
            }
         }
         else
         {  /* xB[i] is decreasing */
            if (phase == 1 && coef[k] > 0.0)
            {  /* xB[i] violates its upper bound, which plays the role
                  of an lower bound on phase I */
               delta = rtol * (1.0 + kappa * fabs(ub[k]));
               t = ((ub[k] - delta) - bbar[i]) / alfa;
               i_stat = GLP_NU;
            }
            else if (phase == 1 && coef[k] < 0.0)
            {  /* xB[i] violates its lower bound, which plays the role
                  of an upper bound on phase I */
               continue;
            }
            else if (type[k] == GLP_LO || type[k] == GLP_DB ||
                     type[k] == GLP_FX)
            {  /* xB[i] is within its bounds and has an lower bound */
               delta = rtol * (1.0 + kappa * fabs(lb[k]));
               t = ((lb[k] - delta) - bbar[i]) / alfa;
               i_stat = GLP_NL;
            }
            else
            {  /* xB[i] is within its bounds and has no lower bound */
               continue;
            }
         }
         /* t is a change of xN[q], on which xB[i] reaches its bound
            (possibly relaxed); since the basic solution is assumed to
            be primal feasible (or pseudo feasible on phase I), t has
            to be non-negative by definition; however, it may happen
            that xB[i] slightly (i.e. within a tolerance) violates its
            bound, that leads to negative t; in the latter case, if
            xB[i] is chosen, negative t means that xN[q] changes in
            wrong direction; if pivot alfa[i,q] is close to zero, even
            small bound violation of xB[i] may lead to a large change
            of xN[q] in wrong direction; let, for example, xB[i] >= 0
            and in the current basis its value be -5e-9; let also xN[q]
            be on its zero bound and should increase; from the ratio
            test rule it follows that the pivot alfa[i,q] < 0; however,
            if alfa[i,q] is, say, -1e-9, the change of xN[q] in wrong
            direction is 5e-9 / (-1e-9) = -5, and using it for updating
            values of other basic variables will give absolutely wrong
            results; therefore, if t is negative, we should replace it
            by exact zero assuming that xB[i] is exactly on its bound,
            and the violation appears due to round-off errors */
         if (t < 0.0) t = 0.0;
         /* apply minimal ratio test */
         if (teta > t || teta == t && big < fabs(alfa))
            p = i, p_stat = i_stat, teta = t, big = fabs(alfa);
      }
      /* the second pass is skipped in the following cases: */
      /* if the standard ratio test is used */
      if (rtol == 0.0) goto done;
      /* if xN[q] reaches its opposite bound or if no basic variable
         has been chosen on the first pass */
      if (p <= 0) goto done;
      /* if xB[p] is a blocking variable, i.e. if it prevents xN[q]
         from any change */
      if (teta == 0.0) goto done;
      /*** SECOND PASS ***/
      /* here tmax is a maximal change of xN[q], on which the solution
         remains primal feasible (or pseudo feasible on phase I) within
         a tolerance */
#if 0
      tmax = (1.0 + 10.0 * DBL_EPSILON) * teta;
#else
      tmax = teta;
#endif
      /* nothing is chosen so far */
      p = 0, p_stat = 0, teta = DBL_MAX, big = 0.0;
      /* walk through significant elements of the pivot column */
      for (pos = 1; pos <= tcol_num; pos++)
      {  i = tcol_ind[pos];
#ifdef GLP_DEBUG
         xassert(1 <= i && i <= m);
#endif
         k = head[i]; /* x[k] = xB[i] */
#ifdef GLP_DEBUG
         xassert(1 <= k && k <= m+n);
#endif
         alfa = s * tcol_vec[i];
#ifdef GLP_DEBUG
         xassert(alfa != 0.0);
#endif
         /* xB[i] = ... + alfa * xN[q] + ..., and due to s we need to
            consider the only case when xN[q] is increasing */
         if (alfa > 0.0)
         {  /* xB[i] is increasing */
            if (phase == 1 && coef[k] < 0.0)
            {  /* xB[i] violates its lower bound, which plays the role
                  of an upper bound on phase I */
               t = (lb[k] - bbar[i]) / alfa;
               i_stat = GLP_NL;
            }
            else if (phase == 1 && coef[k] > 0.0)
            {  /* xB[i] violates its upper bound, which plays the role
                  of an lower bound on phase I */
               continue;
            }
            else if (type[k] == GLP_UP || type[k] == GLP_DB ||
                     type[k] == GLP_FX)
            {  /* xB[i] is within its bounds and has an upper bound */
               t = (ub[k] - bbar[i]) / alfa;
               i_stat = GLP_NU;
            }
            else
            {  /* xB[i] is within its bounds and has no upper bound */
               continue;
            }
         }
         else
         {  /* xB[i] is decreasing */
            if (phase == 1 && coef[k] > 0.0)
            {  /* xB[i] violates its upper bound, which plays the role
                  of an lower bound on phase I */
               t = (ub[k] - bbar[i]) / alfa;
               i_stat = GLP_NU;
            }
            else if (phase == 1 && coef[k] < 0.0)
            {  /* xB[i] violates its lower bound, which plays the role
                  of an upper bound on phase I */
               continue;
            }
            else if (type[k] == GLP_LO || type[k] == GLP_DB ||
                     type[k] == GLP_FX)
            {  /* xB[i] is within its bounds and has an lower bound */
               t = (lb[k] - bbar[i]) / alfa;
               i_stat = GLP_NL;
            }
            else
            {  /* xB[i] is within its bounds and has no lower bound */
               continue;
            }
         }
         /* (see comments for the first pass) */
         if (t < 0.0) t = 0.0;
         /* t is a change of xN[q], on which xB[i] reaches its bound;
            if t <= tmax, all basic variables can violate their bounds
            only within relaxation tolerance delta; we can use this
            freedom and choose basic variable having largest influence
            coefficient to avoid possible numeric instability */
         if (t <= tmax && big < fabs(alfa))
            p = i, p_stat = i_stat, teta = t, big = fabs(alfa);
      }
      /* something must be chosen on the second pass */
      xassert(p != 0);
done: /* store the index and status of basic variable xB[p] chosen */
      csa->p = p;
      if (p > 0 && type[head[p]] == GLP_FX)
         csa->p_stat = GLP_NS;
      else
         csa->p_stat = p_stat;
      /* store corresponding change of non-basic variable xN[q] */
#ifdef GLP_DEBUG
      xassert(teta >= 0.0);
#endif
      csa->teta = s * teta;
      return;
}

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
*  variables are not computed. */

static void eval_trow(struct csa *csa, double rho[])
{     int m = csa->m;
      int n = csa->n;
#ifdef GLP_DEBUG
      char *stat = csa->stat;
#endif
      int *N_ptr = csa->N_ptr;
      int *N_len = csa->N_len;
      int *N_ind = csa->N_ind;
      double *N_val = csa->N_val;
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
         beg = N_ptr[i];
         end = beg + N_len[i];
         for (ptr = beg; ptr < end; ptr++)
         {
#ifdef GLP_DEBUG
            j = N_ind[ptr];
            xassert(1 <= j && j <= n);
            xassert(stat[j] != GLP_NS);
#endif
            trow_vec[N_ind[ptr]] -= temp * N_val[ptr];
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
      int q = csa->q;
      int tcol_nnz = csa->tcol_nnz;
      int *tcol_ind = csa->tcol_ind;
      double *tcol_vec = csa->tcol_vec;
      int p = csa->p;
      double teta = csa->teta;
      int i, pos;
#ifdef GLP_DEBUG
      xassert(1 <= q && q <= n);
      xassert(p < 0 || 1 <= p && p <= m);
#endif
      /* if xN[q] leaves the basis, compute its value in the adjacent
         basis, where it will replace xB[p] */
      if (p > 0)
         bbar[p] = get_xN(csa, q) + teta;
      /* update values of other basic variables (except xB[p], because
         it will be replaced by xN[q]) */
      if (teta == 0.0) goto done;
      for (pos = 1; pos <= tcol_nnz; pos++)
      {  i = tcol_ind[pos];
         /* skip xB[p] */
         if (i == p) continue;
         /* (change of xB[i]) = alfa[i,q] * (change of xN[q]) */
         bbar[i] += tcol_vec[i] * teta;
      }
done: return;
}

/***********************************************************************
*  reeval_cost - recompute reduced cost of non-basic variable xN[q]
*
*  This routine recomputes reduced cost of non-basic variable xN[q] for
*  the current basis more accurately using its direct definition:
*
*     d[q] = cN[q] - N'[q] * pi =
*
*          = cN[q] - N'[q] * (inv(B') * cB) =
*
*          = cN[q] - (cB' * inv(B) * N[q]) =
*
*          = cN[q] + cB' * (pivot column).
*
*  It is assumed that the pivot column of the simplex table is already
*  computed. */

static double reeval_cost(struct csa *csa)
{     int m = csa->m;
#ifdef GLP_DEBUG
      int n = csa->n;
#endif
      double *coef = csa->coef;
      int *head = csa->head;
      int q = csa->q;
      int tcol_nnz = csa->tcol_nnz;
      int *tcol_ind = csa->tcol_ind;
      double *tcol_vec = csa->tcol_vec;
      int i, pos;
      double dq;
#ifdef GLP_DEBUG
      xassert(1 <= q && q <= n);
#endif
      dq = coef[head[m+q]];
      for (pos = 1; pos <= tcol_nnz; pos++)
      {  i = tcol_ind[pos];
#ifdef GLP_DEBUG
         xassert(1 <= i && i <= m);
#endif
         dq += coef[head[i]] * tcol_vec[i];
      }
      return dq;
}

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
      int q = csa->q;
      int trow_nnz = csa->trow_nnz;
      int *trow_ind = csa->trow_ind;
      double *trow_vec = csa->trow_vec;
      int j, pos;
      double new_dq;
#ifdef GLP_DEBUG
      xassert(1 <= q && q <= n);
#endif
      /* compute reduced cost of xB[p] in the adjacent basis, where it
         will replace xN[q] */
#ifdef GLP_DEBUG
      xassert(trow_vec[q] != 0.0);
#endif
      new_dq = (cbar[q] /= trow_vec[q]);
      /* update reduced costs of other non-basic variables (except
         xN[q], because it will be replaced by xB[p]) */
      for (pos = 1; pos <= trow_nnz; pos++)
      {  j = trow_ind[pos];
         /* skip xN[q] */
         if (j == q) continue;
         cbar[j] -= trow_vec[j] * new_dq;
      }
      return;
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
      int *A_ptr = csa->A_ptr;
      int *A_ind = csa->A_ind;
      double *A_val = csa->A_val;
      int *head = csa->head;
      char *refsp = csa->refsp;
      double *gamma = csa->gamma;
      int q = csa->q;
      int tcol_nnz = csa->tcol_nnz;
      int *tcol_ind = csa->tcol_ind;
      double *tcol_vec = csa->tcol_vec;
      int p = csa->p;
      int trow_nnz = csa->trow_nnz;
      int *trow_ind = csa->trow_ind;
      double *trow_vec = csa->trow_vec;
      double *u = csa->work3;
      int i, j, k, pos, beg, end, ptr;
      double gamma_q, delta_q, pivot, s, t, t1, t2;
#ifdef GLP_DEBUG
      xassert(1 <= p && p <= m);
      xassert(1 <= q && q <= n);
#endif
      /* the basis changes, so decrease the count */
      xassert(csa->refct > 0);
      csa->refct--;
      /* recompute gamma[q] for the current basis more accurately and
         compute auxiliary vector u */
      gamma_q = delta_q = (refsp[head[m+q]] ? 1.0 : 0.0);
      for (i = 1; i <= m; i++) u[i] = 0.0;
      for (pos = 1; pos <= tcol_nnz; pos++)
      {  i = tcol_ind[pos];
         if (refsp[head[i]])
         {  u[i] = t = tcol_vec[i];
            gamma_q += t * t;
         }
         else
            u[i] = 0.0;
      }
      xassert(csa->valid);
      bfd_btran(csa->bfd, u);
      /* update gamma[k] for other non-basic variables (except fixed
         variables and xN[q], because it will be replaced by xB[p]) */
      pivot = trow_vec[q];
#ifdef GLP_DEBUG
      xassert(pivot != 0.0);
#endif
      for (pos = 1; pos <= trow_nnz; pos++)
      {  j = trow_ind[pos];
         /* skip xN[q] */
         if (j == q) continue;
         /* compute t */
         t = trow_vec[j] / pivot;
         /* compute inner product s = N'[j] * u */
         k = head[m+j]; /* x[k] = xN[j] */
         if (k <= m)
            s = u[k];
         else
         {  s = 0.0;
            beg = A_ptr[k-m];
            end = A_ptr[k-m+1];
            for (ptr = beg; ptr < end; ptr++)
               s -= A_val[ptr] * u[A_ind[ptr]];
         }
         /* compute gamma[k] for the adjacent basis */
         t1 = gamma[j] + t * t * gamma_q + 2.0 * t * s;
         t2 = (refsp[k] ? 1.0 : 0.0) + delta_q * t * t;
         gamma[j] = (t1 >= t2 ? t1 : t2);
         if (gamma[j] < DBL_EPSILON) gamma[j] = DBL_EPSILON;
      }
      /* compute gamma[q] for the adjacent basis */
      if (type[head[p]] == GLP_FX)
         gamma[q] = 1.0;
      else
      {  gamma[q] = gamma_q / (pivot * pivot);
         if (gamma[q] < DBL_EPSILON) gamma[q] = DBL_EPSILON;
      }
      return;
}

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
{     int n = csa->n;
      char *stat = csa->stat;
      double *gamma = csa->gamma;
      int j;
      double e, emax, temp;
      emax = 0.0;
      for (j = 1; j <= n; j++)
      {  if (stat[j] == GLP_NS)
         {  xassert(gamma[j] == 1.0);
            continue;
         }
         temp = eval_gamma(csa, j);
         e = fabs(temp - gamma[j]) / (1.0 + fabs(temp));
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
      char *type = csa->type;
#endif
      int *head = csa->head;
      char *stat = csa->stat;
      int q = csa->q;
      int p = csa->p;
      int p_stat = csa->p_stat;
      int k;
#ifdef GLP_DEBUG
      xassert(1 <= q && q <= n);
#endif
      if (p < 0)
      {  /* xN[q] goes to its opposite bound */
#ifdef GLP_DEBUG
         k = head[m+q]; /* x[k] = xN[q] */
         xassert(1 <= k && k <= m+n);
         xassert(type[k] == GLP_DB);
#endif
         switch (stat[q])
         {  case GLP_NL:
               /* xN[q] increases */
               stat[q] = GLP_NU;
               break;
            case GLP_NU:
               /* xN[q] decreases */
               stat[q] = GLP_NL;
               break;
            default:
               xassert(stat != stat);
         }
      }
      else
      {  /* xB[p] leaves the basis, xN[q] enters the basis */
#ifdef GLP_DEBUG
         xassert(1 <= p && p <= m);
         k = head[p]; /* x[k] = xB[p] */
         switch (p_stat)
         {  case GLP_NL:
               /* xB[p] goes to its lower bound */
               xassert(type[k] == GLP_LO || type[k] == GLP_DB);
               break;
            case GLP_NU:
               /* xB[p] goes to its upper bound */
               xassert(type[k] == GLP_UP || type[k] == GLP_DB);
               break;
            case GLP_NS:
               /* xB[p] goes to its fixed value */
               xassert(type[k] == GLP_NS);
               break;
            default:
               xassert(p_stat != p_stat);
         }
#endif
         /* xB[p] <-> xN[q] */
         k = head[p], head[p] = head[m+q], head[m+q] = k;
         stat[q] = (char)p_stat;
      }
      return;
}

/***********************************************************************
*  set_aux_obj - construct auxiliary objective function
*
*  The auxiliary objective function is a separable piecewise linear
*  convex function, which is the sum of primal infeasibilities:
*
*     z = t[1] + ... + t[m+n] -> minimize,
*
*  where:
*
*            / lb[k] - x[k], if x[k] < lb[k]
*            |
*     t[k] = <  0, if lb[k] <= x[k] <= ub[k]
*            |
*            \ x[k] - ub[k], if x[k] > ub[k]
*
*  This routine computes objective coefficients for the current basis
*  and returns the number of non-zero terms t[k]. */

static int set_aux_obj(struct csa *csa, double tol_bnd)
{     int m = csa->m;
      int n = csa->n;
      char *type = csa->type;
      double *lb = csa->lb;
      double *ub = csa->ub;
      double *coef = csa->coef;
      int *head = csa->head;
      double *bbar = csa->bbar;
      int i, k, cnt = 0;
      double eps;
      /* use a bit more restrictive tolerance */
      tol_bnd *= 0.90;
      /* clear all objective coefficients */
      for (k = 1; k <= m+n; k++)
         coef[k] = 0.0;
      /* walk through the list of basic variables */
      for (i = 1; i <= m; i++)
      {  k = head[i]; /* x[k] = xB[i] */
         if (type[k] == GLP_LO || type[k] == GLP_DB ||
             type[k] == GLP_FX)
         {  /* x[k] has lower bound */
            eps = tol_bnd * (1.0 + kappa * fabs(lb[k]));
            if (bbar[i] < lb[k] - eps)
            {  /* and violates it */
               coef[k] = -1.0;
               cnt++;
            }
         }
         if (type[k] == GLP_UP || type[k] == GLP_DB ||
             type[k] == GLP_FX)
         {  /* x[k] has upper bound */
            eps = tol_bnd * (1.0 + kappa * fabs(ub[k]));
            if (bbar[i] > ub[k] + eps)
            {  /* and violates it */
               coef[k] = +1.0;
               cnt++;
            }
         }
      }
      return cnt;
}

/***********************************************************************
*  set_orig_obj - restore original objective function
*
*  This routine assigns scaled original objective coefficients to the
*  working objective function. */

static void set_orig_obj(struct csa *csa)
{     int m = csa->m;
      int n = csa->n;
      double *coef = csa->coef;
      double *obj = csa->obj;
      double zeta = csa->zeta;
      int i, j;
      for (i = 1; i <= m; i++)
         coef[i] = 0.0;
      for (j = 1; j <= n; j++)
         coef[m+j] = zeta * obj[j];
      return;
}

/***********************************************************************
*  check_stab - check numerical stability of basic solution
*
*  If the current basic solution is primal feasible (or pseudo feasible
*  on phase I) within a tolerance, this routine returns zero, otherwise
*  it returns non-zero. */

static int check_stab(struct csa *csa, double tol_bnd)
{     int m = csa->m;
#ifdef GLP_DEBUG
      int n = csa->n;
#endif
      char *type = csa->type;
      double *lb = csa->lb;
      double *ub = csa->ub;
      double *coef = csa->coef;
      int *head = csa->head;
      int phase = csa->phase;
      double *bbar = csa->bbar;
      int i, k;
      double eps;
      /* walk through the list of basic variables */
      for (i = 1; i <= m; i++)
      {  k = head[i]; /* x[k] = xB[i] */
#ifdef GLP_DEBUG
         xassert(1 <= k && k <= m+n);
#endif
         if (phase == 1 && coef[k] < 0.0)
         {  /* x[k] must not be greater than its lower bound */
#ifdef GLP_DEBUG
            xassert(type[k] == GLP_LO || type[k] == GLP_DB ||
                    type[k] == GLP_FX);
#endif
            eps = tol_bnd * (1.0 + kappa * fabs(lb[k]));
            if (bbar[i] > lb[k] + eps) return 1;
         }
         else if (phase == 1 && coef[k] > 0.0)
         {  /* x[k] must not be less than its upper bound */
#ifdef GLP_DEBUG
            xassert(type[k] == GLP_UP || type[k] == GLP_DB ||
                    type[k] == GLP_FX);
#endif
            eps = tol_bnd * (1.0 + kappa * fabs(ub[k]));
            if (bbar[i] < ub[k] - eps) return 1;
         }
         else
         {  /* either phase = 1 and coef[k] = 0, or phase = 2 */
            if (type[k] == GLP_LO || type[k] == GLP_DB ||
                type[k] == GLP_FX)
            {  /* x[k] must not be less than its lower bound */
               eps = tol_bnd * (1.0 + kappa * fabs(lb[k]));
               if (bbar[i] < lb[k] - eps) return 1;
            }
            if (type[k] == GLP_UP || type[k] == GLP_DB ||
                type[k] == GLP_FX)
            {  /* x[k] must not be greater then its upper bound */
               eps = tol_bnd * (1.0 + kappa * fabs(ub[k]));
               if (bbar[i] > ub[k] + eps) return 1;
            }
         }
      }
      /* basic solution is primal feasible within a tolerance */
      return 0;
}

/***********************************************************************
*  check_feas - check primal feasibility of basic solution
*
*  If the current basic solution is primal feasible within a tolerance,
*  this routine returns zero, otherwise it returns non-zero. */

static int check_feas(struct csa *csa, double tol_bnd)
{     int m = csa->m;
#ifdef GLP_DEBUG
      int n = csa->n;
      char *type = csa->type;
#endif
      double *lb = csa->lb;
      double *ub = csa->ub;
      double *coef = csa->coef;
      int *head = csa->head;
      double *bbar = csa->bbar;
      int i, k;
      double eps;
      xassert(csa->phase == 1);
      /* walk through the list of basic variables */
      for (i = 1; i <= m; i++)
      {  k = head[i]; /* x[k] = xB[i] */
#ifdef GLP_DEBUG
         xassert(1 <= k && k <= m+n);
#endif
         if (coef[k] < 0.0)
         {  /* check if x[k] still violates its lower bound */
#ifdef GLP_DEBUG
            xassert(type[k] == GLP_LO || type[k] == GLP_DB ||
                    type[k] == GLP_FX);
#endif
            eps = tol_bnd * (1.0 + kappa * fabs(lb[k]));
            if (bbar[i] < lb[k] - eps) return 1;
         }
         else if (coef[k] > 0.0)
         {  /* check if x[k] still violates its upper bound */
#ifdef GLP_DEBUG
            xassert(type[k] == GLP_UP || type[k] == GLP_DB ||
                    type[k] == GLP_FX);
#endif
            eps = tol_bnd * (1.0 + kappa * fabs(ub[k]));
            if (bbar[i] > ub[k] + eps) return 1;
         }
      }
      /* basic solution is primal feasible within a tolerance */
      return 0;
}

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

/***********************************************************************
*  display - display the search progress
*
*  This routine displays some information about the search progress
*  that includes:
*
*  the search phase;
*
*  the number of simplex iterations performed by the solver;
*
*  the original objective value;
*
*  the sum of (scaled) primal infeasibilities;
*
*  the number of basic fixed variables. */

static void display(struct csa *csa, const glp_smcp *parm, int spec)
{     int m = csa->m;
#ifdef GLP_DEBUG
      int n = csa->n;
#endif
      char *type = csa->type;
      double *lb = csa->lb;
      double *ub = csa->ub;
      int phase = csa->phase;
      int *head = csa->head;
      double *bbar = csa->bbar;
      int i, k, cnt;
      double sum;
      if (parm->msg_lev < GLP_MSG_ON) goto skip;
      if (parm->out_dly > 0 &&
         1000.0 * xdifftime(xtime(), csa->tm_beg) < parm->out_dly)
         goto skip;
      if (csa->it_cnt == csa->it_dpy) goto skip;
      if (!spec && csa->it_cnt % parm->out_frq != 0) goto skip;
      /* compute the sum of primal infeasibilities and determine the
         number of basic fixed variables */
      sum = 0.0, cnt = 0;
      for (i = 1; i <= m; i++)
      {  k = head[i]; /* x[k] = xB[i] */
#ifdef GLP_DEBUG
         xassert(1 <= k && k <= m+n);
#endif
         if (type[k] == GLP_LO || type[k] == GLP_DB ||
             type[k] == GLP_FX)
         {  /* x[k] has lower bound */
            if (bbar[i] < lb[k])
               sum += (lb[k] - bbar[i]);
         }
         if (type[k] == GLP_UP || type[k] == GLP_DB ||
             type[k] == GLP_FX)
         {  /* x[k] has upper bound */
            if (bbar[i] > ub[k])
               sum += (bbar[i] - ub[k]);
         }
         if (type[k] == GLP_FX) cnt++;
      }
      xprintf("%c%6d: obj = %17.9e  infeas = %10.3e (%d)\n",
         phase == 1 ? ' ' : '*', csa->it_cnt, eval_obj(csa), sum, cnt);
      csa->it_dpy = csa->it_cnt;
skip: return;
}

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
      xfree(csa->obj);
      xfree(csa->A_ptr);
      xfree(csa->A_ind);
      xfree(csa->A_val);
      xfree(csa->head);
      xfree(csa->stat);
      xfree(csa->N_ptr);
      xfree(csa->N_len);
      xfree(csa->N_ind);
      xfree(csa->N_val);
      xfree(csa->bbar);
      xfree(csa->cbar);
      xfree(csa->refsp);
      xfree(csa->gamma);
      xfree(csa->tcol_ind);
      xfree(csa->tcol_vec);
      xfree(csa->trow_ind);
      xfree(csa->trow_vec);
      xfree(csa->work1);
      xfree(csa->work2);
      xfree(csa->work3);
      xfree(csa->work4);
      xfree(csa);
      return;
}

/***********************************************************************
*  spx_primal - core LP solver based on the primal simplex method
*
*  SYNOPSIS
*
*  #include "glpspx.h"
*  int spx_primal(glp_prob *lp, const glp_smcp *parm);
*
*  DESCRIPTION
*
*  The routine spx_primal is a core LP solver based on the two-phase
*  primal simplex method.
*
*  RETURNS
*
*  0  LP instance has been successfully solved.
*
*  GLP_EITLIM
*     Iteration limit has been exhausted.
*
*  GLP_ETMLIM
*     Time limit has been exhausted.
*
*  GLP_EFAIL
*     The solver failed to solve LP instance. */

int spx_primal(glp_prob *lp, const glp_smcp *parm)
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
      /* compute primal values of basic variables */
      if (bbar_st == 0)
      {  eval_bbar(csa);
         bbar_st = 1; /* just computed */
         /* determine the search phase, if not determined yet */
         if (csa->phase == 0)
         {  if (set_aux_obj(csa, parm->tol_bnd) > 0)
            {  /* current basic solution is primal infeasible */
               /* start to minimize the sum of infeasibilities */
               csa->phase = 1;
            }
            else
            {  /* current basic solution is primal feasible */
               /* start to minimize the original objective function */
               set_orig_obj(csa);
               csa->phase = 2;
            }
            xassert(check_stab(csa, parm->tol_bnd) == 0);
            /* working objective coefficients have been changed, so
               invalidate reduced costs */
            cbar_st = 0;
            display(csa, parm, 1);
         }
         /* make sure that the current basic solution remains primal
            feasible (or pseudo feasible on phase I) */
         if (check_stab(csa, parm->tol_bnd))
         {  /* there are excessive bound violations due to round-off
               errors */
            if (parm->msg_lev >= GLP_MSG_ERR)
               xprintf("Warning: numerical instability (primal simplex,"
                  " phase %s)\n", csa->phase == 1 ? "I" : "II");
            /* restart the search */
            csa->phase = 0;
            binv_st = 0;
            rigorous = 5;
            goto loop;
         }
      }
      xassert(csa->phase == 1 || csa->phase == 2);
      /* on phase I we do not need to wait until the current basic
         solution becomes dual feasible; it is sufficient to make sure
         that no basic variable violates its bounds */
      if (csa->phase == 1 && !check_feas(csa, parm->tol_bnd))
      {  /* the current basis is primal feasible; switch to phase II */
         csa->phase = 2;
         set_orig_obj(csa);
         cbar_st = 0;
         display(csa, parm, 1);
      }
      /* compute reduced costs of non-basic variables */
      if (cbar_st == 0)
      {  eval_cbar(csa);
         cbar_st = 1; /* just computed */
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
      /* check if the iteration limit has been exhausted */
      if (parm->it_lim < INT_MAX &&
          csa->it_cnt - csa->it_beg >= parm->it_lim)
      {  if (bbar_st != 1 || csa->phase == 2 && cbar_st != 1)
         {  if (bbar_st != 1) bbar_st = 0;
            if (csa->phase == 2 && cbar_st != 1) cbar_st = 0;
            goto loop;
         }
         display(csa, parm, 1);
         if (parm->msg_lev >= GLP_MSG_ALL)
            xprintf("ITERATION LIMIT EXCEEDED; SEARCH TERMINATED\n");
         switch (csa->phase)
         {  case 1:
               p_stat = GLP_INFEAS;
               set_orig_obj(csa);
               eval_cbar(csa);
               break;
            case 2:
               p_stat = GLP_FEAS;
               break;
            default:
               xassert(csa != csa);
         }
         chuzc(csa, parm->tol_dj);
         d_stat = (csa->q == 0 ? GLP_FEAS : GLP_INFEAS);
         store_sol(csa, lp, p_stat, d_stat, 0);
         ret = GLP_EITLIM;
         goto done;
      }
      /* check if the time limit has been exhausted */
      if (parm->tm_lim < INT_MAX &&
          1000.0 * xdifftime(xtime(), csa->tm_beg) >= parm->tm_lim)
      {  if (bbar_st != 1 || csa->phase == 2 && cbar_st != 1)
         {  if (bbar_st != 1) bbar_st = 0;
            if (csa->phase == 2 && cbar_st != 1) cbar_st = 0;
            goto loop;
         }
         display(csa, parm, 1);
         if (parm->msg_lev >= GLP_MSG_ALL)
            xprintf("TIME LIMIT EXCEEDED; SEARCH TERMINATED\n");
         switch (csa->phase)
         {  case 1:
               p_stat = GLP_INFEAS;
               set_orig_obj(csa);
               eval_cbar(csa);
               break;
            case 2:
               p_stat = GLP_FEAS;
               break;
            default:
               xassert(csa != csa);
         }
         chuzc(csa, parm->tol_dj);
         d_stat = (csa->q == 0 ? GLP_FEAS : GLP_INFEAS);
         store_sol(csa, lp, p_stat, d_stat, 0);
         ret = GLP_ETMLIM;
         goto done;
      }
      /* display the search progress */
      display(csa, parm, 0);
      /* choose non-basic variable xN[q] */
      chuzc(csa, parm->tol_dj);
      if (csa->q == 0)
      {  if (bbar_st != 1 || cbar_st != 1)
         {  if (bbar_st != 1) bbar_st = 0;
            if (cbar_st != 1) cbar_st = 0;
            goto loop;
         }
         display(csa, parm, 1);
         switch (csa->phase)
         {  case 1:
               if (parm->msg_lev >= GLP_MSG_ALL)
                  xprintf("PROBLEM HAS NO FEASIBLE SOLUTION\n");
               p_stat = GLP_NOFEAS;
               set_orig_obj(csa);
               eval_cbar(csa);
               chuzc(csa, parm->tol_dj);
               d_stat = (csa->q == 0 ? GLP_FEAS : GLP_INFEAS);
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
      /* compute pivot column of the simplex table */
      eval_tcol(csa);
      if (rigorous) refine_tcol(csa);
      sort_tcol(csa, parm->tol_piv);
      /* check accuracy of the reduced cost of xN[q] */
      {  double d1 = csa->cbar[csa->q]; /* less accurate */
         double d2 = reeval_cost(csa);  /* more accurate */
         xassert(d1 != 0.0);
         if (fabs(d1 - d2) > 1e-5 * (1.0 + fabs(d2)) ||
             !(d1 < 0.0 && d2 < 0.0 || d1 > 0.0 && d2 > 0.0))
         {  if (parm->msg_lev >= GLP_MSG_DBG)
               xprintf("d1 = %.12g; d2 = %.12g\n", d1, d2);
            if (cbar_st != 1 || !rigorous)
            {  if (cbar_st != 1) cbar_st = 0;
               rigorous = 5;
               goto loop;
            }
         }
         /* replace cbar[q] by more accurate value keeping its sign */
         if (d1 > 0.0)
            csa->cbar[csa->q] = (d2 > 0.0 ? d2 : +DBL_EPSILON);
         else
            csa->cbar[csa->q] = (d2 < 0.0 ? d2 : -DBL_EPSILON);
      }
      /* choose basic variable xB[p] */
      switch (parm->r_test)
      {  case GLP_RT_STD:
            chuzr(csa, 0.0);
            break;
         case GLP_RT_HAR:
            chuzr(csa, 0.30 * parm->tol_bnd);
            break;
         default:
            xassert(parm != parm);
      }
      if (csa->p == 0)
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
                  xprintf("PROBLEM HAS UNBOUNDED SOLUTION\n");
               store_sol(csa, lp, GLP_FEAS, GLP_NOFEAS,
                  csa->head[csa->m+csa->q]);
               ret = 0;
               break;
            default:
               xassert(csa != csa);
         }
         goto done;
      }
      /* check if the pivot element is acceptable */
      if (csa->p > 0)
      {  double piv = csa->tcol_vec[csa->p];
         double eps = 1e-5 * (1.0 + 0.01 * csa->tcol_max);
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
      /* compute pivot row of the simplex table */
      if (csa->p > 0)
      {  double *rho = csa->work4;
         eval_rho(csa, rho);
         if (rigorous) refine_rho(csa, rho);
         eval_trow(csa, rho);
      }
      /* accuracy check based on the pivot element */
      if (csa->p > 0)
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
            /* use more accurate version in the pivot row */
            if (csa->trow_vec[csa->q] == 0.0)
            {  csa->trow_nnz++;
               xassert(csa->trow_nnz <= csa->n);
               csa->trow_ind[csa->trow_nnz] = csa->q;
            }
            csa->trow_vec[csa->q] = piv1;
         }
      }
      /* update primal values of basic variables */
      update_bbar(csa);
      bbar_st = 2; /* updated */
      /* update reduced costs of non-basic variables */
      if (csa->p > 0)
      {  update_cbar(csa);
         cbar_st = 2; /* updated */
         /* on phase I objective coefficient of xB[p] in the adjacent
            basis becomes zero */
         if (csa->phase == 1)
         {  int k = csa->head[csa->p]; /* x[k] = xB[p] -> xN[q] */
            csa->cbar[csa->q] -= csa->coef[k];
            csa->coef[k] = 0.0;
         }
      }
      /* update steepest edge coefficients */
      if (csa->p > 0)
      {  switch (parm->pricing)
         {  case GLP_PT_STD:
               break;
            case GLP_PT_PSE:
               if (csa->refct > 0) update_gamma(csa);
               break;
            default:
               xassert(parm != parm);
         }
      }
      /* update factorization of the basis matrix */
      if (csa->p > 0)
      {  ret = update_B(csa, csa->p, csa->head[csa->m+csa->q]);
         if (ret == 0)
            binv_st = 2; /* updated */
         else
         {  csa->valid = 0;
            binv_st = 0; /* invalid */
         }
      }
      /* update matrix N */
      if (csa->p > 0)
      {  del_N_col(csa, csa->q, csa->head[csa->m+csa->q]);
         if (csa->type[csa->head[csa->p]] != GLP_FX)
            add_N_col(csa, csa->q, csa->head[csa->p]);
      }
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
