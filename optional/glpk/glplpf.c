/* glplpf.c (LP basis factorization, Schur complement version) */

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

#include "glplpf.h"
#include "glpenv.h"
#define xfault xerror

#define _GLPLPF_DEBUG 0

/* CAUTION: DO NOT CHANGE THE LIMIT BELOW */

#define M_MAX 100000000 /* = 100*10^6 */
/* maximal order of the basis matrix */

/***********************************************************************
*  NAME
*
*  lpf_create_it - create LP basis factorization
*
*  SYNOPSIS
*
*  #include "glplpf.h"
*  LPF *lpf_create_it(void);
*
*  DESCRIPTION
*
*  The routine lpf_create_it creates a program object, which represents
*  a factorization of LP basis.
*
*  RETURNS
*
*  The routine lpf_create_it returns a pointer to the object created. */

LPF *lpf_create_it(void)
{     LPF *lpf;
#if _GLPLPF_DEBUG
      xprintf("lpf_create_it: warning: debug mode enabled\n");
#endif
      lpf = xmalloc(sizeof(LPF));
      lpf->valid = 0;
      lpf->m0_max = lpf->m0 = 0;
      lpf->luf = luf_create_it();
      lpf->m = 0;
      lpf->B = NULL;
      lpf->n_max = 50;
      lpf->n = 0;
      lpf->R_ptr = lpf->R_len = NULL;
      lpf->S_ptr = lpf->S_len = NULL;
      lpf->scf = NULL;
      lpf->P_row = lpf->P_col = NULL;
      lpf->Q_row = lpf->Q_col = NULL;
      lpf->v_size = 1000;
      lpf->v_ptr = 0;
      lpf->v_ind = NULL;
      lpf->v_val = NULL;
      lpf->work1 = lpf->work2 = NULL;
      return lpf;
}

/***********************************************************************
*  NAME
*
*  lpf_factorize - compute LP basis factorization
*
*  SYNOPSIS
*
*  #include "glplpf.h"
*  int lpf_factorize(LPF *lpf, int m, const int bh[], int (*col)
*     (void *info, int j, int ind[], double val[]), void *info);
*
*  DESCRIPTION
*
*  The routine lpf_factorize computes the factorization of the basis
*  matrix B specified by the routine col.
*
*  The parameter lpf specified the basis factorization data structure
*  created with the routine lpf_create_it.
*
*  The parameter m specifies the order of B, m > 0.
*
*  The array bh specifies the basis header: bh[j], 1 <= j <= m, is the
*  number of j-th column of B in some original matrix. The array bh is
*  optional and can be specified as NULL.
*
*  The formal routine col specifies the matrix B to be factorized. To
*  obtain j-th column of A the routine lpf_factorize calls the routine
*  col with the parameter j (1 <= j <= n). In response the routine col
*  should store row indices and numerical values of non-zero elements
*  of j-th column of B to locations ind[1,...,len] and val[1,...,len],
*  respectively, where len is the number of non-zeros in j-th column
*  returned on exit. Neither zero nor duplicate elements are allowed.
*
*  The parameter info is a transit pointer passed to the routine col.
*
*  RETURNS
*
*  0  The factorization has been successfully computed.
*
*  LPF_ESING
*     The specified matrix is singular within the working precision.
*
*  LPF_ECOND
*     The specified matrix is ill-conditioned.
*
*  For more details see comments to the routine luf_factorize. */

int lpf_factorize(LPF *lpf, int m, const int bh[], int (*col)
      (void *info, int j, int ind[], double val[]), void *info)
{     int k, ret;
#if _GLPLPF_DEBUG
      int i, j, len, *ind;
      double *B, *val;
#endif
      xassert(bh == bh);
      if (m < 1)
         xfault("lpf_factorize: m = %d; invalid parameter\n", m);
      if (m > M_MAX)
         xfault("lpf_factorize: m = %d; matrix too big\n", m);
      lpf->m0 = lpf->m = m;
      /* invalidate the factorization */
      lpf->valid = 0;
      /* allocate/reallocate arrays, if necessary */
      if (lpf->R_ptr == NULL)
         lpf->R_ptr = xcalloc(1+lpf->n_max, sizeof(int));
      if (lpf->R_len == NULL)
         lpf->R_len = xcalloc(1+lpf->n_max, sizeof(int));
      if (lpf->S_ptr == NULL)
         lpf->S_ptr = xcalloc(1+lpf->n_max, sizeof(int));
      if (lpf->S_len == NULL)
         lpf->S_len = xcalloc(1+lpf->n_max, sizeof(int));
      if (lpf->scf == NULL)
         lpf->scf = scf_create_it(lpf->n_max);
      if (lpf->v_ind == NULL)
         lpf->v_ind = xcalloc(1+lpf->v_size, sizeof(int));
      if (lpf->v_val == NULL)
         lpf->v_val = xcalloc(1+lpf->v_size, sizeof(double));
      if (lpf->m0_max < m)
      {  if (lpf->P_row != NULL) xfree(lpf->P_row);
         if (lpf->P_col != NULL) xfree(lpf->P_col);
         if (lpf->Q_row != NULL) xfree(lpf->Q_row);
         if (lpf->Q_col != NULL) xfree(lpf->Q_col);
         if (lpf->work1 != NULL) xfree(lpf->work1);
         if (lpf->work2 != NULL) xfree(lpf->work2);
         lpf->m0_max = m + 100;
         lpf->P_row = xcalloc(1+lpf->m0_max+lpf->n_max, sizeof(int));
         lpf->P_col = xcalloc(1+lpf->m0_max+lpf->n_max, sizeof(int));
         lpf->Q_row = xcalloc(1+lpf->m0_max+lpf->n_max, sizeof(int));
         lpf->Q_col = xcalloc(1+lpf->m0_max+lpf->n_max, sizeof(int));
         lpf->work1 = xcalloc(1+lpf->m0_max+lpf->n_max, sizeof(double));
         lpf->work2 = xcalloc(1+lpf->m0_max+lpf->n_max, sizeof(double));
      }
      /* try to factorize the basis matrix */
      switch (luf_factorize(lpf->luf, m, col, info))
      {  case 0:
            break;
         case LUF_ESING:
            ret = LPF_ESING;
            goto done;
         case LUF_ECOND:
            ret = LPF_ECOND;
            goto done;
         default:
            xassert(lpf != lpf);
      }
      /* the basis matrix has been successfully factorized */
      lpf->valid = 1;
#if _GLPLPF_DEBUG
      /* store the basis matrix for debugging */
      if (lpf->B != NULL) xfree(lpf->B);
      xassert(m <= 32767);
      lpf->B = B = xcalloc(1+m*m, sizeof(double));
      ind = xcalloc(1+m, sizeof(int));
      val = xcalloc(1+m, sizeof(double));
      for (k = 1; k <= m * m; k++)
         B[k] = 0.0;
      for (j = 1; j <= m; j++)
      {  len = col(info, j, ind, val);
         xassert(0 <= len && len <= m);
         for (k = 1; k <= len; k++)
         {  i = ind[k];
            xassert(1 <= i && i <= m);
            xassert(B[(i - 1) * m + j] == 0.0);
            xassert(val[k] != 0.0);
            B[(i - 1) * m + j] = val[k];
         }
      }
      xfree(ind);
      xfree(val);
#endif
      /* B = B0, so there are no additional rows/columns */
      lpf->n = 0;
      /* reset the Schur complement factorization */
      scf_reset_it(lpf->scf);
      /* P := Q := I */
      for (k = 1; k <= m; k++)
      {  lpf->P_row[k] = lpf->P_col[k] = k;
         lpf->Q_row[k] = lpf->Q_col[k] = k;
      }
      /* make all SVA locations free */
      lpf->v_ptr = 1;
      ret = 0;
done: /* return to the calling program */
      return ret;
}

/***********************************************************************
*  The routine r_prod computes the product y := y + alpha * R * x,
*  where x is a n-vector, alpha is a scalar, y is a m0-vector.
*
*  Since matrix R is available by columns, the product is computed as
*  a linear combination:
*
*     y := y + alpha * (R[1] * x[1] + ... + R[n] * x[n]),
*
*  where R[j] is j-th column of R. */

static void r_prod(LPF *lpf, double y[], double a, const double x[])
{     int n = lpf->n;
      int *R_ptr = lpf->R_ptr;
      int *R_len = lpf->R_len;
      int *v_ind = lpf->v_ind;
      double *v_val = lpf->v_val;
      int j, beg, end, ptr;
      double t;
      for (j = 1; j <= n; j++)
      {  if (x[j] == 0.0) continue;
         /* y := y + alpha * R[j] * x[j] */
         t = a * x[j];
         beg = R_ptr[j];
         end = beg + R_len[j];
         for (ptr = beg; ptr < end; ptr++)
            y[v_ind[ptr]] += t * v_val[ptr];
      }
      return;
}

/***********************************************************************
*  The routine rt_prod computes the product y := y + alpha * R' * x,
*  where R' is a matrix transposed to R, x is a m0-vector, alpha is a
*  scalar, y is a n-vector.
*
*  Since matrix R is available by columns, the product components are
*  computed as inner products:
*
*     y[j] := y[j] + alpha * (j-th column of R) * x
*
*  for j = 1, 2, ..., n. */

static void rt_prod(LPF *lpf, double y[], double a, const double x[])
{     int n = lpf->n;
      int *R_ptr = lpf->R_ptr;
      int *R_len = lpf->R_len;
      int *v_ind = lpf->v_ind;
      double *v_val = lpf->v_val;
      int j, beg, end, ptr;
      double t;
      for (j = 1; j <= n; j++)
      {  /* t := (j-th column of R) * x */
         t = 0.0;
         beg = R_ptr[j];
         end = beg + R_len[j];
         for (ptr = beg; ptr < end; ptr++)
            t += v_val[ptr] * x[v_ind[ptr]];
         /* y[j] := y[j] + alpha * t */
         y[j] += a * t;
      }
      return;
}

/***********************************************************************
*  The routine s_prod computes the product y := y + alpha * S * x,
*  where x is a m0-vector, alpha is a scalar, y is a n-vector.
*
*  Since matrix S is available by rows, the product components are
*  computed as inner products:
*
*     y[i] = y[i] + alpha * (i-th row of S) * x
*
*  for i = 1, 2, ..., n. */

static void s_prod(LPF *lpf, double y[], double a, const double x[])
{     int n = lpf->n;
      int *S_ptr = lpf->S_ptr;
      int *S_len = lpf->S_len;
      int *v_ind = lpf->v_ind;
      double *v_val = lpf->v_val;
      int i, beg, end, ptr;
      double t;
      for (i = 1; i <= n; i++)
      {  /* t := (i-th row of S) * x */
         t = 0.0;
         beg = S_ptr[i];
         end = beg + S_len[i];
         for (ptr = beg; ptr < end; ptr++)
            t += v_val[ptr] * x[v_ind[ptr]];
         /* y[i] := y[i] + alpha * t */
         y[i] += a * t;
      }
      return;
}

/***********************************************************************
*  The routine st_prod computes the product y := y + alpha * S' * x,
*  where S' is a matrix transposed to S, x is a n-vector, alpha is a
*  scalar, y is m0-vector.
*
*  Since matrix R is available by rows, the product is computed as a
*  linear combination:
*
*     y := y + alpha * (S'[1] * x[1] + ... + S'[n] * x[n]),
*
*  where S'[i] is i-th row of S. */

static void st_prod(LPF *lpf, double y[], double a, const double x[])
{     int n = lpf->n;
      int *S_ptr = lpf->S_ptr;
      int *S_len = lpf->S_len;
      int *v_ind = lpf->v_ind;
      double *v_val = lpf->v_val;
      int i, beg, end, ptr;
      double t;
      for (i = 1; i <= n; i++)
      {  if (x[i] == 0.0) continue;
         /* y := y + alpha * S'[i] * x[i] */
         t = a * x[i];
         beg = S_ptr[i];
         end = beg + S_len[i];
         for (ptr = beg; ptr < end; ptr++)
            y[v_ind[ptr]] += t * v_val[ptr];
      }
      return;
}

#if _GLPLPF_DEBUG
/***********************************************************************
*  The routine check_error computes the maximal relative error between
*  left- and right-hand sides for the system B * x = b (if tr is zero)
*  or B' * x = b (if tr is non-zero), where B' is a matrix transposed
*  to B. (This routine is intended for debugging only.) */

static void check_error(LPF *lpf, int tr, const double x[],
      const double b[])
{     int m = lpf->m;
      double *B = lpf->B;
      int i, j;
      double  d, dmax = 0.0, s, t, tmax;
      for (i = 1; i <= m; i++)
      {  s = 0.0;
         tmax = 1.0;
         for (j = 1; j <= m; j++)
         {  if (!tr)
               t = B[m * (i - 1) + j] * x[j];
            else
               t = B[m * (j - 1) + i] * x[j];
            if (tmax < fabs(t)) tmax = fabs(t);
            s += t;
         }
         d = fabs(s - b[i]) / tmax;
         if (dmax < d) dmax = d;
      }
      if (dmax > 1e-8)
         xprintf("%s: dmax = %g; relative error too large\n",
            !tr ? "lpf_ftran" : "lpf_btran", dmax);
      return;
}
#endif

/***********************************************************************
*  NAME
*
*  lpf_ftran - perform forward transformation (solve system B*x = b)
*
*  SYNOPSIS
*
*  #include "glplpf.h"
*  void lpf_ftran(LPF *lpf, double x[]);
*
*  DESCRIPTION
*
*  The routine lpf_ftran performs forward transformation, i.e. solves
*  the system B*x = b, where B is the basis matrix, x is the vector of
*  unknowns to be computed, b is the vector of right-hand sides.
*
*  On entry elements of the vector b should be stored in dense format
*  in locations x[1], ..., x[m], where m is the number of rows. On exit
*  the routine stores elements of the vector x in the same locations.
*
*  BACKGROUND
*
*  Solution of the system B * x = b can be obtained by solving the
*  following augmented system:
*
*     ( B  F^) ( x )   ( b )
*     (      ) (   ) = (   )
*     ( G^ H^) ( y )   ( 0 )
*
*  which, using the main equality, can be written as follows:
*
*       ( L0 0 ) ( U0 R )   ( x )   ( b )
*     P (      ) (      ) Q (   ) = (   )
*       ( S  I ) ( 0  C )   ( y )   ( 0 )
*
*  therefore,
*
*     ( x )      ( U0 R )-1 ( L0 0 )-1    ( b )
*     (   ) = Q' (      )   (      )   P' (   )
*     ( y )      ( 0  C )   ( S  I )      ( 0 )
*
*  Thus, computing the solution includes the following steps:
*
*  1. Compute
*
*     ( f )      ( b )
*     (   ) = P' (   )
*     ( g )      ( 0 )
*
*  2. Solve the system
*
*     ( f1 )   ( L0 0 )-1 ( f )      ( L0 0 ) ( f1 )   ( f )
*     (    ) = (      )   (   )  =>  (      ) (    ) = (   )
*     ( g1 )   ( S  I )   ( g )      ( S  I ) ( g1 )   ( g )
*
*     from which it follows that:
*
*     { L0 * f1      = f      f1 = inv(L0) * f
*     {                   =>
*     {  S * f1 + g1 = g      g1 = g - S * f1
*
*  3. Solve the system
*
*     ( f2 )   ( U0 R )-1 ( f1 )      ( U0 R ) ( f2 )   ( f1 )
*     (    ) = (      )   (    )  =>  (      ) (    ) = (    )
*     ( g2 )   ( 0  C )   ( g1 )      ( 0  C ) ( g2 )   ( g1 )
*
*     from which it follows that:
*
*     { U0 * f2 + R * g2 = f1      f2 = inv(U0) * (f1 - R * g2)
*     {                        =>
*     {           C * g2 = g1      g2 = inv(C) * g1
*
*  4. Compute
*
*     ( x )      ( f2 )
*     (   ) = Q' (    )
*     ( y )      ( g2 )                                               */

void lpf_ftran(LPF *lpf, double x[])
{     int m0 = lpf->m0;
      int m = lpf->m;
      int n  = lpf->n;
      int *P_col = lpf->P_col;
      int *Q_col = lpf->Q_col;
      double *fg = lpf->work1;
      double *f = fg;
      double *g = fg + m0;
      int i, ii;
#if _GLPLPF_DEBUG
      double *b;
#endif
      if (!lpf->valid)
         xfault("lpf_ftran: the factorization is not valid\n");
      xassert(0 <= m && m <= m0 + n);
#if _GLPLPF_DEBUG
      /* save the right-hand side vector */
      b = xcalloc(1+m, sizeof(double));
      for (i = 1; i <= m; i++) b[i] = x[i];
#endif
      /* (f g) := inv(P) * (b 0) */
      for (i = 1; i <= m0 + n; i++)
         fg[i] = ((ii = P_col[i]) <= m ? x[ii] : 0.0);
      /* f1 := inv(L0) * f */
      luf_f_solve(lpf->luf, 0, f);
      /* g1 := g - S * f1 */
      s_prod(lpf, g, -1.0, f);
      /* g2 := inv(C) * g1 */
      scf_solve_it(lpf->scf, 0, g);
      /* f2 := inv(U0) * (f1 - R * g2) */
      r_prod(lpf, f, -1.0, g);
      luf_v_solve(lpf->luf, 0, f);
      /* (x y) := inv(Q) * (f2 g2) */
      for (i = 1; i <= m; i++)
         x[i] = fg[Q_col[i]];
#if _GLPLPF_DEBUG
      /* check relative error in solution */
      check_error(lpf, 0, x, b);
      xfree(b);
#endif
      return;
}

/***********************************************************************
*  NAME
*
*  lpf_btran - perform backward transformation (solve system B'*x = b)
*
*  SYNOPSIS
*
*  #include "glplpf.h"
*  void lpf_btran(LPF *lpf, double x[]);
*
*  DESCRIPTION
*
*  The routine lpf_btran performs backward transformation, i.e. solves
*  the system B'*x = b, where B' is a matrix transposed to the basis
*  matrix B, x is the vector of unknowns to be computed, b is the vector
*  of right-hand sides.
*
*  On entry elements of the vector b should be stored in dense format
*  in locations x[1], ..., x[m], where m is the number of rows. On exit
*  the routine stores elements of the vector x in the same locations.
*
*  BACKGROUND
*
*  Solution of the system B' * x = b, where B' is a matrix transposed
*  to B, can be obtained by solving the following augmented system:
*
*     ( B  F^)T ( x )   ( b )
*     (      )  (   ) = (   )
*     ( G^ H^)  ( y )   ( 0 )
*
*  which, using the main equality, can be written as follows:
*
*      T ( U0 R )T ( L0 0 )T  T ( x )   ( b )
*     Q  (      )  (      )  P  (   ) = (   )
*        ( 0  C )  ( S  I )     ( y )   ( 0 )
*
*  or, equivalently, as follows:
*
*        ( U'0 0 ) ( L'0 S')    ( x )   ( b )
*     Q' (       ) (       ) P' (   ) = (   )
*        ( R'  C') ( 0   I )    ( y )   ( 0 )
*
*  therefore,
*
*     ( x )     ( L'0 S')-1 ( U'0 0 )-1   ( b )
*     (   ) = P (       )   (       )   Q (   )
*     ( y )     ( 0   I )   ( R'  C')     ( 0 )
*
*  Thus, computing the solution includes the following steps:
*
*  1. Compute
*
*     ( f )     ( b )
*     (   ) = Q (   )
*     ( g )     ( 0 )
*
*  2. Solve the system
*
*     ( f1 )   ( U'0 0 )-1 ( f )      ( U'0 0 ) ( f1 )   ( f )
*     (    ) = (       )   (   )  =>  (       ) (    ) = (   )
*     ( g1 )   ( R'  C')   ( g )      ( R'  C') ( g1 )   ( g )
*
*     from which it follows that:
*
*     { U'0 * f1           = f      f1 = inv(U'0) * f
*     {                         =>
*     { R'  * f1 + C' * g1 = g      g1 = inv(C') * (g - R' * f1)
*
*  3. Solve the system
*
*     ( f2 )   ( L'0 S')-1 ( f1 )      ( L'0 S') ( f2 )   ( f1 )
*     (    ) = (       )   (    )  =>  (       ) (    ) = (    )
*     ( g2 )   ( 0   I )   ( g1 )      ( 0   I ) ( g2 )   ( g1 )
*
*     from which it follows that:
*
*     { L'0 * f2 + S' * g2 = f1
*     {                          =>  f2 = inv(L'0) * ( f1 - S' * g2)
*     {                 g2 = g1
*
*  4. Compute
*
*     ( x )     ( f2 )
*     (   ) = P (    )
*     ( y )     ( g2 )                                                */

void lpf_btran(LPF *lpf, double x[])
{     int m0 = lpf->m0;
      int m = lpf->m;
      int n = lpf->n;
      int *P_row = lpf->P_row;
      int *Q_row = lpf->Q_row;
      double *fg = lpf->work1;
      double *f = fg;
      double *g = fg + m0;
      int i, ii;
#if _GLPLPF_DEBUG
      double *b;
#endif
      if (!lpf->valid)
         xfault("lpf_btran: the factorization is not valid\n");
      xassert(0 <= m && m <= m0 + n);
#if _GLPLPF_DEBUG
      /* save the right-hand side vector */
      b = xcalloc(1+m, sizeof(double));
      for (i = 1; i <= m; i++) b[i] = x[i];
#endif
      /* (f g) := Q * (b 0) */
      for (i = 1; i <= m0 + n; i++)
         fg[i] = ((ii = Q_row[i]) <= m ? x[ii] : 0.0);
      /* f1 := inv(U'0) * f */
      luf_v_solve(lpf->luf, 1, f);
      /* g1 := inv(C') * (g - R' * f1) */
      rt_prod(lpf, g, -1.0, f);
      scf_solve_it(lpf->scf, 1, g);
      /* g2 := g1 */
      g = g;
      /* f2 := inv(L'0) * (f1 - S' * g2) */
      st_prod(lpf, f, -1.0, g);
      luf_f_solve(lpf->luf, 1, f);
      /* (x y) := P * (f2 g2) */
      for (i = 1; i <= m; i++)
         x[i] = fg[P_row[i]];
#if _GLPLPF_DEBUG
      /* check relative error in solution */
      check_error(lpf, 1, x, b);
      xfree(b);
#endif
      return;
}

/***********************************************************************
*  The routine enlarge_sva enlarges the Sparse Vector Area to new_size
*  locations by reallocating the arrays v_ind and v_val. */

static void enlarge_sva(LPF *lpf, int new_size)
{     int v_size = lpf->v_size;
      int used = lpf->v_ptr - 1;
      int *v_ind = lpf->v_ind;
      double *v_val = lpf->v_val;
      xassert(v_size < new_size);
      while (v_size < new_size) v_size += v_size;
      lpf->v_size = v_size;
      lpf->v_ind = xcalloc(1+v_size, sizeof(int));
      lpf->v_val = xcalloc(1+v_size, sizeof(double));
      xassert(used >= 0);
      memcpy(&lpf->v_ind[1], &v_ind[1], used * sizeof(int));
      memcpy(&lpf->v_val[1], &v_val[1], used * sizeof(double));
      xfree(v_ind);
      xfree(v_val);
      return;
}

/***********************************************************************
*  NAME
*
*  lpf_update_it - update LP basis factorization
*
*  SYNOPSIS
*
*  #include "glplpf.h"
*  int lpf_update_it(LPF *lpf, int j, int bh, int len, const int ind[],
*     const double val[]);
*
*  DESCRIPTION
*
*  The routine lpf_update_it updates the factorization of the basis
*  matrix B after replacing its j-th column by a new vector.
*
*  The parameter j specifies the number of column of B, which has been
*  replaced, 1 <= j <= m, where m is the order of B.
*
*  The parameter bh specifies the basis header entry for the new column
*  of B, which is the number of the new column in some original matrix.
*  This parameter is optional and can be specified as 0.
*
*  Row indices and numerical values of non-zero elements of the new
*  column of B should be placed in locations ind[1], ..., ind[len] and
*  val[1], ..., val[len], resp., where len is the number of non-zeros
*  in the column. Neither zero nor duplicate elements are allowed.
*
*  RETURNS
*
*  0  The factorization has been successfully updated.
*
*  LPF_ESING
*     New basis B is singular within the working precision.
*
*  LPF_ELIMIT
*     Maximal number of additional rows and columns has been reached.
*
*  BACKGROUND
*
*  Let j-th column of the current basis matrix B have to be replaced by
*  a new column a. This replacement is equivalent to removing the old
*  j-th column by fixing it at zero and introducing the new column as
*  follows:
*
*                   ( B   F^| a )
*     ( B  F^)      (       |   )
*     (      ) ---> ( G^  H^| 0 )
*     ( G^ H^)      (-------+---)
*                   ( e'j 0 | 0 )
*
*  where ej is a unit vector with 1 in j-th position which used to fix
*  the old j-th column of B (at zero). Then using the main equality we
*  have:
*
*     ( B   F^| a )            ( B0  F | f )
*     (       |   )   ( P  0 ) (       |   ) ( Q  0 )
*     ( G^  H^| 0 ) = (      ) ( G   H | g ) (      ) =
*     (-------+---)   ( 0  1 ) (-------+---) ( 0  1 )
*     ( e'j 0 | 0 )            ( v'  w'| 0 )
*
*       [   ( B0  F )|   ( f ) ]            [   ( B0 F )  |   ( f ) ]
*       [ P (       )| P (   ) ] ( Q  0 )   [ P (      ) Q| P (   ) ]
*     = [   ( G   H )|   ( g ) ] (      ) = [   ( G  H )  |   ( g ) ]
*       [------------+-------- ] ( 0  1 )   [-------------+---------]
*       [   ( v'  w')|     0   ]            [   ( v' w') Q|    0    ]
*
*  where:
*
*     ( a )     ( f )      ( f )        ( a )
*     (   ) = P (   )  =>  (   ) = P' * (   )
*     ( 0 )     ( g )      ( g )        ( 0 )
*
*                                 ( ej )      ( v )    ( v )     ( ej )
*     ( e'j  0 ) = ( v' w' ) Q => (    ) = Q' (   ) => (   ) = Q (    )
*                                 ( 0  )      ( w )    ( w )     ( 0  )
*
*  On the other hand:
*
*              ( B0| F  f )
*     ( P  0 ) (---+------) ( Q  0 )         ( B0    new F )
*     (      ) ( G | H  g ) (      ) = new P (             ) new Q
*     ( 0  1 ) (   |      ) ( 0  1 )         ( new G new H )
*              ( v'| w' 0 )
*
*  where:
*                               ( G )           ( H  g )
*     new F = ( F  f ), new G = (   ),  new H = (      ),
*                               ( v')           ( w' 0 )
*
*             ( P  0 )           ( Q  0 )
*     new P = (      ) , new Q = (      ) .
*             ( 0  1 )           ( 0  1 )
*
*  The factorization structure for the new augmented matrix remains the
*  same, therefore:
*
*           ( B0    new F )         ( L0     0 ) ( U0 new R )
*     new P (             ) new Q = (          ) (          )
*           ( new G new H )         ( new S  I ) ( 0  new C )
*
*  where:
*
*     new F = L0 * new R  =>
*
*     new R = inv(L0) * new F = inv(L0) * (F  f) = ( R  inv(L0)*f )
*
*     new G = new S * U0  =>
*
*                               ( G )             (     S      )
*     new S = new G * inv(U0) = (   ) * inv(U0) = (            )
*                               ( v')             ( v'*inv(U0) )
*
*     new H = new S * new R + new C  =>
*
*     new C = new H - new S * new R =
*
*             ( H  g )   (     S      )
*           = (      ) - (            ) * ( R  inv(L0)*f ) =
*             ( w' 0 )   ( v'*inv(U0) )
*
*             ( H - S*R           g - S*inv(L0)*f      )   ( C  x )
*           = (                                        ) = (      )
*             ( w'- v'*inv(U0)*R  -v'*inv(U0)*inv(L0)*f)   ( y' z )
*
*  Note that new C is resulted by expanding old C with new column x,
*  row y', and diagonal element z, where:
*
*     x = g - S * inv(L0) * f = g - S * (new column of R)
*
*     y = w - R'* inv(U'0)* v = w - R'* (new row of S)
*
*     z = - (new row of S) * (new column of R)
*
*  Finally, to replace old B by new B we have to permute j-th and last
*  (just added) columns of the matrix
*
*     ( B   F^| a )
*     (       |   )
*     ( G^  H^| 0 )
*     (-------+---)
*     ( e'j 0 | 0 )
*
*  and to keep the main equality do the same for matrix Q. */

int lpf_update_it(LPF *lpf, int j, int bh, int len, const int ind[],
      const double val[])
{     int m0 = lpf->m0;
      int m = lpf->m;
#if _GLPLPF_DEBUG
      double *B = lpf->B;
#endif
      int n = lpf->n;
      int *R_ptr = lpf->R_ptr;
      int *R_len = lpf->R_len;
      int *S_ptr = lpf->S_ptr;
      int *S_len = lpf->S_len;
      int *P_row = lpf->P_row;
      int *P_col = lpf->P_col;
      int *Q_row = lpf->Q_row;
      int *Q_col = lpf->Q_col;
      int v_ptr = lpf->v_ptr;
      int *v_ind = lpf->v_ind;
      double *v_val = lpf->v_val;
      double *a = lpf->work2; /* new column */
      double *fg = lpf->work1, *f = fg, *g = fg + m0;
      double *vw = lpf->work2, *v = vw, *w = vw + m0;
      double *x = g, *y = w, z;
      int i, ii, k, ret;
      xassert(bh == bh);
      if (!lpf->valid)
         xfault("lpf_update_it: the factorization is not valid\n");
      if (!(1 <= j && j <= m))
         xfault("lpf_update_it: j = %d; column number out of range\n",
            j);
      xassert(0 <= m && m <= m0 + n);
      /* check if the basis factorization can be expanded */
      if (n == lpf->n_max)
      {  lpf->valid = 0;
         ret = LPF_ELIMIT;
         goto done;
      }
      /* convert new j-th column of B to dense format */
      for (i = 1; i <= m; i++)
         a[i] = 0.0;
      for (k = 1; k <= len; k++)
      {  i = ind[k];
         if (!(1 <= i && i <= m))
            xfault("lpf_update_it: ind[%d] = %d; row number out of rang"
               "e\n", k, i);
         if (a[i] != 0.0)
            xfault("lpf_update_it: ind[%d] = %d; duplicate row index no"
               "t allowed\n", k, i);
         if (val[k] == 0.0)
            xfault("lpf_update_it: val[%d] = %g; zero element not allow"
               "ed\n", k, val[k]);
         a[i] = val[k];
      }
#if _GLPLPF_DEBUG
      /* change column in the basis matrix for debugging */
      for (i = 1; i <= m; i++)
         B[(i - 1) * m + j] = a[i];
#endif
      /* (f g) := inv(P) * (a 0) */
      for (i = 1; i <= m0+n; i++)
         fg[i] = ((ii = P_col[i]) <= m ? a[ii] : 0.0);
      /* (v w) := Q * (ej 0) */
      for (i = 1; i <= m0+n; i++) vw[i] = 0.0;
      vw[Q_col[j]] = 1.0;
      /* f1 := inv(L0) * f (new column of R) */
      luf_f_solve(lpf->luf, 0, f);
      /* v1 := inv(U'0) * v (new row of S) */
      luf_v_solve(lpf->luf, 1, v);
      /* we need at most 2 * m0 available locations in the SVA to store
         new column of matrix R and new row of matrix S */
      if (lpf->v_size < v_ptr + m0 + m0)
      {  enlarge_sva(lpf, v_ptr + m0 + m0);
         v_ind = lpf->v_ind;
         v_val = lpf->v_val;
      }
      /* store new column of R */
      R_ptr[n+1] = v_ptr;
      for (i = 1; i <= m0; i++)
      {  if (f[i] != 0.0)
            v_ind[v_ptr] = i, v_val[v_ptr] = f[i], v_ptr++;
      }
      R_len[n+1] = v_ptr - lpf->v_ptr;
      lpf->v_ptr = v_ptr;
      /* store new row of S */
      S_ptr[n+1] = v_ptr;
      for (i = 1; i <= m0; i++)
      {  if (v[i] != 0.0)
            v_ind[v_ptr] = i, v_val[v_ptr] = v[i], v_ptr++;
      }
      S_len[n+1] = v_ptr - lpf->v_ptr;
      lpf->v_ptr = v_ptr;
      /* x := g - S * f1 (new column of C) */
      s_prod(lpf, x, -1.0, f);
      /* y := w - R' * v1 (new row of C) */
      rt_prod(lpf, y, -1.0, v);
      /* z := - v1 * f1 (new diagonal element of C) */
      z = 0.0;
      for (i = 1; i <= m0; i++) z -= v[i] * f[i];
      /* update factorization of new matrix C */
      switch (scf_update_exp(lpf->scf, x, y, z))
      {  case 0:
            break;
         case SCF_ESING:
            lpf->valid = 0;
            ret = LPF_ESING;
            goto done;
         case SCF_ELIMIT:
            xassert(lpf != lpf);
         default:
            xassert(lpf != lpf);
      }
      /* expand matrix P */
      P_row[m0+n+1] = P_col[m0+n+1] = m0+n+1;
      /* expand matrix Q */
      Q_row[m0+n+1] = Q_col[m0+n+1] = m0+n+1;
      /* permute j-th and last (just added) column of matrix Q */
      i = Q_col[j], ii = Q_col[m0+n+1];
      Q_row[i] = m0+n+1, Q_col[m0+n+1] = i;
      Q_row[ii] = j, Q_col[j] = ii;
      /* increase the number of additional rows and columns */
      lpf->n++;
      xassert(lpf->n <= lpf->n_max);
      /* the factorization has been successfully updated */
      ret = 0;
done: /* return to the calling program */
      return ret;
}

/***********************************************************************
*  NAME
*
*  lpf_delete_it - delete LP basis factorization
*
*  SYNOPSIS
*
*  #include "glplpf.h"
*  void lpf_delete_it(LPF *lpf)
*
*  DESCRIPTION
*
*  The routine lpf_delete_it deletes LP basis factorization specified
*  by the parameter lpf and frees all memory allocated to this program
*  object. */

void lpf_delete_it(LPF *lpf)
{     luf_delete_it(lpf->luf);
#if _GLPLPF_DEBUG
      if (lpf->B != NULL) xfree(lpf->B);
#else
      xassert(lpf->B == NULL);
#endif
      if (lpf->R_ptr != NULL) xfree(lpf->R_ptr);
      if (lpf->R_len != NULL) xfree(lpf->R_len);
      if (lpf->S_ptr != NULL) xfree(lpf->S_ptr);
      if (lpf->S_len != NULL) xfree(lpf->S_len);
      if (lpf->scf != NULL) scf_delete_it(lpf->scf);
      if (lpf->P_row != NULL) xfree(lpf->P_row);
      if (lpf->P_col != NULL) xfree(lpf->P_col);
      if (lpf->Q_row != NULL) xfree(lpf->Q_row);
      if (lpf->Q_col != NULL) xfree(lpf->Q_col);
      if (lpf->v_ind != NULL) xfree(lpf->v_ind);
      if (lpf->v_val != NULL) xfree(lpf->v_val);
      if (lpf->work1 != NULL) xfree(lpf->work1);
      if (lpf->work2 != NULL) xfree(lpf->work2);
      xfree(lpf);
      return;
}

/* eof */
