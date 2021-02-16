/* glpscf.c (Schur complement factorization) */

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
#pragma clang diagnostic ignored "-Wself-assign"
#endif

#include "glpenv.h"
#include "glpscf.h"
#define xfault xerror

#define _GLPSCF_DEBUG 0

#define eps 1e-10

/***********************************************************************
*  NAME
*
*  scf_create_it - create Schur complement factorization
*
*  SYNOPSIS
*
*  #include "glpscf.h"
*  SCF *scf_create_it(int n_max);
*
*  DESCRIPTION
*
*  The routine scf_create_it creates the factorization of matrix C,
*  which initially has no rows and columns.
*
*  The parameter n_max specifies the maximal order of matrix C to be
*  factorized, 1 <= n_max <= 32767.
*
*  RETURNS
*
*  The routine scf_create_it returns a pointer to the structure SCF,
*  which defines the factorization. */

SCF *scf_create_it(int n_max)
{     SCF *scf;
#if _GLPSCF_DEBUG
      xprintf("scf_create_it: warning: debug mode enabled\n");
#endif
      if (!(1 <= n_max && n_max <= 32767))
         xfault("scf_create_it: n_max = %d; invalid parameter\n",
            n_max);
      scf = xmalloc(sizeof(SCF));
      scf->n_max = n_max;
      scf->n = 0;
      scf->f = xcalloc(1 + n_max * n_max, sizeof(double));
      scf->u = xcalloc(1 + n_max * (n_max + 1) / 2, sizeof(double));
      scf->p = xcalloc(1 + n_max, sizeof(int));
      scf->t_opt = SCF_TBG;
      scf->rank = 0;
#if _GLPSCF_DEBUG
      scf->c = xcalloc(1 + n_max * n_max, sizeof(double));
#else
      scf->c = NULL;
#endif
      scf->w = xcalloc(1 + n_max, sizeof(double));
      return scf;
}

/***********************************************************************
*  The routine f_loc determines location of matrix element F[i,j] in
*  the one-dimensional array f. */

static int f_loc(SCF *scf, int i, int j)
{     int n_max = scf->n_max;
      int n = scf->n;
      xassert(1 <= i && i <= n);
      xassert(1 <= j && j <= n);
      return (i - 1) * n_max + j;
}

/***********************************************************************
*  The routine u_loc determines location of matrix element U[i,j] in
*  the one-dimensional array u. */

static int u_loc(SCF *scf, int i, int j)
{     int n_max = scf->n_max;
      int n = scf->n;
      xassert(1 <= i && i <= n);
      xassert(i <= j && j <= n);
      return (i - 1) * n_max + j - i * (i - 1) / 2;
}

/***********************************************************************
*  The routine bg_transform applies Bartels-Golub version of gaussian
*  elimination to restore triangular structure of matrix U.
*
*  On entry matrix U has the following structure:
*
*        1       k         n
*     1  * * * * * * * * * *
*        . * * * * * * * * *
*        . . * * * * * * * *
*        . . . * * * * * * *
*     k  . . . . * * * * * *
*        . . . . . * * * * *
*        . . . . . . * * * *
*        . . . . . . . * * *
*        . . . . . . . . * *
*     n  . . . . # # # # # #
*
*  where '#' is a row spike to be eliminated.
*
*  Elements of n-th row are passed separately in locations un[k], ...,
*  un[n]. On exit the content of the array un is destroyed.
*
*  REFERENCES
*
*  R.H.Bartels, G.H.Golub, "The Simplex Method of Linear Programming
*  Using LU-decomposition", Comm. ACM, 12, pp. 266-68, 1969. */

static void bg_transform(SCF *scf, int k, double un[])
{     int n = scf->n;
      double *f = scf->f;
      double *u = scf->u;
      int j, k1, kj, kk, n1, nj;
      double t;
      xassert(1 <= k && k <= n);
      /* main elimination loop */
      for (k = k; k < n; k++)
      {  /* determine location of U[k,k] */
         kk = u_loc(scf, k, k);
         /* determine location of F[k,1] */
         k1 = f_loc(scf, k, 1);
         /* determine location of F[n,1] */
         n1 = f_loc(scf, n, 1);
         /* if |U[k,k]| < |U[n,k]|, interchange k-th and n-th rows to
            provide |U[k,k]| >= |U[n,k]| */
         if (fabs(u[kk]) < fabs(un[k]))
         {  /* interchange k-th and n-th rows of matrix U */
            for (j = k, kj = kk; j <= n; j++, kj++)
               t = u[kj], u[kj] = un[j], un[j] = t;
            /* interchange k-th and n-th rows of matrix F to keep the
               main equality F * C = U * P */
            for (j = 1, kj = k1, nj = n1; j <= n; j++, kj++, nj++)
               t = f[kj], f[kj] = f[nj], f[nj] = t;
         }
         /* now |U[k,k]| >= |U[n,k]| */
         /* if U[k,k] is too small in the magnitude, replace U[k,k] and
            U[n,k] by exact zero */
         if (fabs(u[kk]) < eps) u[kk] = un[k] = 0.0;
         /* if U[n,k] is already zero, elimination is not needed */
         if (un[k] == 0.0) continue;
         /* compute gaussian multiplier t = U[n,k] / U[k,k] */
         t = un[k] / u[kk];
         /* apply gaussian elimination to nullify U[n,k] */
         /* (n-th row of U) := (n-th row of U) - t * (k-th row of U) */
         for (j = k+1, kj = kk+1; j <= n; j++, kj++)
            un[j] -= t * u[kj];
         /* (n-th row of F) := (n-th row of F) - t * (k-th row of F)
            to keep the main equality F * C = U * P */
         for (j = 1, kj = k1, nj = n1; j <= n; j++, kj++, nj++)
            f[nj] -= t * f[kj];
      }
      /* if U[n,n] is too small in the magnitude, replace it by exact
         zero */
      if (fabs(un[n]) < eps) un[n] = 0.0;
      /* store U[n,n] in a proper location */
      u[u_loc(scf, n, n)] = un[n];
      return;
}

/***********************************************************************
*  The routine givens computes the parameters of Givens plane rotation
*  c = cos(teta) and s = sin(teta) such that:
*
*     ( c -s ) ( a )   ( r )
*     (      ) (   ) = (   ) ,
*     ( s  c ) ( b )   ( 0 )
*
*  where a and b are given scalars.
*
*  REFERENCES
*
*  G.H.Golub, C.F.Van Loan, "Matrix Computations", 2nd ed. */

static void givens(double a, double b, double *c, double *s)
{     double t;
      if (b == 0.0)
         (*c) = 1.0, (*s) = 0.0;
      else if (fabs(a) <= fabs(b))
         t = - a / b, (*s) = 1.0 / sqrt(1.0 + t * t), (*c) = (*s) * t;
      else
         t = - b / a, (*c) = 1.0 / sqrt(1.0 + t * t), (*s) = (*c) * t;
      return;
}

/*----------------------------------------------------------------------
*  The routine gr_transform applies Givens plane rotations to restore
*  triangular structure of matrix U.
*
*  On entry matrix U has the following structure:
*
*        1       k         n
*     1  * * * * * * * * * *
*        . * * * * * * * * *
*        . . * * * * * * * *
*        . . . * * * * * * *
*     k  . . . . * * * * * *
*        . . . . . * * * * *
*        . . . . . . * * * *
*        . . . . . . . * * *
*        . . . . . . . . * *
*     n  . . . . # # # # # #
*
*  where '#' is a row spike to be eliminated.
*
*  Elements of n-th row are passed separately in locations un[k], ...,
*  un[n]. On exit the content of the array un is destroyed.
*
*  REFERENCES
*
*  R.H.Bartels, G.H.Golub, "The Simplex Method of Linear Programming
*  Using LU-decomposition", Comm. ACM, 12, pp. 266-68, 1969. */

static void gr_transform(SCF *scf, int k, double un[])
{     int n = scf->n;
      double *f = scf->f;
      double *u = scf->u;
      int j, k1, kj, kk, n1, nj;
      double c, s;
      xassert(1 <= k && k <= n);
      /* main elimination loop */
      for (k = k; k < n; k++)
      {  /* determine location of U[k,k] */
         kk = u_loc(scf, k, k);
         /* determine location of F[k,1] */
         k1 = f_loc(scf, k, 1);
         /* determine location of F[n,1] */
         n1 = f_loc(scf, n, 1);
         /* if both U[k,k] and U[n,k] are too small in the magnitude,
            replace them by exact zero */
         if (fabs(u[kk]) < eps && fabs(un[k]) < eps)
            u[kk] = un[k] = 0.0;
         /* if U[n,k] is already zero, elimination is not needed */
         if (un[k] == 0.0) continue;
         /* compute the parameters of Givens plane rotation */
         givens(u[kk], un[k], &c, &s);
         /* apply Givens rotation to k-th and n-th rows of matrix U */
         for (j = k, kj = kk; j <= n; j++, kj++)
         {  double ukj = u[kj], unj = un[j];
            u[kj] = c * ukj - s * unj;
            un[j] = s * ukj + c * unj;
         }
         /* apply Givens rotation to k-th and n-th rows of matrix F
            to keep the main equality F * C = U * P */
         for (j = 1, kj = k1, nj = n1; j <= n; j++, kj++, nj++)
         {  double fkj = f[kj], fnj = f[nj];
            f[kj] = c * fkj - s * fnj;
            f[nj] = s * fkj + c * fnj;
         }
      }
      /* if U[n,n] is too small in the magnitude, replace it by exact
         zero */
      if (fabs(un[n]) < eps) un[n] = 0.0;
      /* store U[n,n] in a proper location */
      u[u_loc(scf, n, n)] = un[n];
      return;
}

/***********************************************************************
*  The routine transform restores triangular structure of matrix U.
*  It is a driver to the routines bg_transform and gr_transform (see
*  comments to these routines above). */

static void transform(SCF *scf, int k, double un[])
{     switch (scf->t_opt)
      {  case SCF_TBG:
            bg_transform(scf, k, un);
            break;
         case SCF_TGR:
            gr_transform(scf, k, un);
            break;
         default:
            xassert(scf != scf);
      }
      return;
}

/***********************************************************************
*  The routine estimate_rank estimates the rank of matrix C.
*
*  Since all transformations applied to matrix F are non-singular,
*  and F is assumed to be well conditioned, from the main equaility
*  F * C = U * P it follows that rank(C) = rank(U), where rank(U) is
*  estimated as the number of non-zero diagonal elements of U. */

static int estimate_rank(SCF *scf)
{     int n_max = scf->n_max;
      int n = scf->n;
      double *u = scf->u;
      int i, ii, inc, rank = 0;
      for (i = 1, ii = u_loc(scf, i, i), inc = n_max; i <= n;
         i++, ii += inc, inc--)
         if (u[ii] != 0.0) rank++;
      return rank;
}

#if _GLPSCF_DEBUG
/***********************************************************************
*  The routine check_error computes the maximal relative error between
*  left- and right-hand sides of the main equality F * C = U * P. (This
*  routine is intended only for debugging.) */

static void check_error(SCF *scf, const char *func)
{     int n = scf->n;
      double *f = scf->f;
      double *u = scf->u;
      int *p = scf->p;
      double *c = scf->c;
      int i, j, k;
      double d, dmax = 0.0, s, t;
      xassert(c != NULL);
      for (i = 1; i <= n; i++)
      {  for (j = 1; j <= n; j++)
         {  /* compute element (i,j) of product F * C */
            s = 0.0;
            for (k = 1; k <= n; k++)
               s += f[f_loc(scf, i, k)] * c[f_loc(scf, k, j)];
            /* compute element (i,j) of product U * P */
            k = p[j];
            t = (i <= k ? u[u_loc(scf, i, k)] : 0.0);
            /* compute the maximal relative error */
            d = fabs(s - t) / (1.0 + fabs(t));
            if (dmax < d) dmax = d;
         }
      }
      if (dmax > 1e-8)
         xprintf("%s: dmax = %g; relative error too large\n", func,
            dmax);
      return;
}
#endif

/***********************************************************************
*  NAME
*
*  scf_update_exp - update factorization on expanding C
*
*  SYNOPSIS
*
*  #include "glpscf.h"
*  int scf_update_exp(SCF *scf, const double x[], const double y[],
*     double z);
*
*  DESCRIPTION
*
*  The routine scf_update_exp updates the factorization of matrix C on
*  expanding it by adding a new row and column as follows:
*
*             ( C  x )
*     new C = (      )
*             ( y' z )
*
*  where x[1,...,n] is a new column, y[1,...,n] is a new row, and z is
*  a new diagonal element.
*
*  If on entry the factorization is empty, the parameters x and y can
*  be specified as NULL.
*
*  RETURNS
*
*  0  The factorization has been successfully updated.
*
*  SCF_ESING
*     The factorization has been successfully updated, however, new
*     matrix C is singular within working precision. Note that the new
*     factorization remains valid.
*
*  SCF_ELIMIT
*     There is not enough room to expand the factorization, because
*     n = n_max. The factorization remains unchanged.
*
*  ALGORITHM
*
*  We can see that:
*
*     ( F  0 ) ( C  x )   ( FC  Fx )   ( UP  Fx )
*     (      ) (      ) = (        ) = (        ) =
*     ( 0  1 ) ( y' z )   ( y'   z )   ( y'   z )
*
*        ( U   Fx ) ( P  0 )
*     =  (        ) (      ),
*        ( y'P' z ) ( 0  1 )
*
*  therefore to keep the main equality F * C = U * P we can take:
*
*             ( F  0 )           ( U   Fx )           ( P  0 )
*     new F = (      ),  new U = (        ),  new P = (      ),
*             ( 0  1 )           ( y'P' z )           ( 0  1 )
*
*  and eliminate the row spike y'P' in the last row of new U to restore
*  its upper triangular structure. */

int scf_update_exp(SCF *scf, const double x[], const double y[],
      double z)
{     int n_max = scf->n_max;
      int n = scf->n;
      double *f = scf->f;
      double *u = scf->u;
      int *p = scf->p;
#if _GLPSCF_DEBUG
      double *c = scf->c;
#endif
      double *un = scf->w;
      int i, ij, in, j, k, nj, ret = 0;
      double t;
      /* check if the factorization can be expanded */
      if (n == n_max)
      {  /* there is not enough room */
         ret = SCF_ELIMIT;
         goto done;
      }
      /* increase the order of the factorization */
      scf->n = ++n;
      /* fill new zero column of matrix F */
      for (i = 1, in = f_loc(scf, i, n); i < n; i++, in += n_max)
         f[in] = 0.0;
      /* fill new zero row of matrix F */
      for (j = 1, nj = f_loc(scf, n, j); j < n; j++, nj++)
         f[nj] = 0.0;
      /* fill new unity diagonal element of matrix F */
      f[f_loc(scf, n, n)] = 1.0;
      /* compute new column of matrix U, which is (old F) * x */
      for (i = 1; i < n; i++)
      {  /* u[i,n] := (i-th row of old F) * x */
         t = 0.0;
         for (j = 1, ij = f_loc(scf, i, 1); j < n; j++, ij++)
            t += f[ij] * x[j];
         u[u_loc(scf, i, n)] = t;
      }
      /* compute new (spiked) row of matrix U, which is (old P) * y */
      for (j = 1; j < n; j++) un[j] = y[p[j]];
      /* store new diagonal element of matrix U, which is z */
      un[n] = z;
      /* expand matrix P */
      p[n] = n;
#if _GLPSCF_DEBUG
      /* expand matrix C */
      /* fill its new column, which is x */
      for (i = 1, in = f_loc(scf, i, n); i < n; i++, in += n_max)
         c[in] = x[i];
      /* fill its new row, which is y */
      for (j = 1, nj = f_loc(scf, n, j); j < n; j++, nj++)
         c[nj] = y[j];
      /* fill its new diagonal element, which is z */
      c[f_loc(scf, n, n)] = z;
#endif
      /* restore upper triangular structure of matrix U */
      for (k = 1; k < n; k++)
         if (un[k] != 0.0) break;
      transform(scf, k, un);
      /* estimate the rank of matrices C and U */
      scf->rank = estimate_rank(scf);
      if (scf->rank != n) ret = SCF_ESING;
#if _GLPSCF_DEBUG
      /* check that the factorization is accurate enough */
      check_error(scf, "scf_update_exp");
#endif
done: return ret;
}

/***********************************************************************
*  The routine solve solves the system C * x = b.
*
*  From the main equation F * C = U * P it follows that:
*
*     C * x = b  =>  F * C * x = F * b  =>  U * P * x = F * b  =>
*
*     P * x = inv(U) * F * b  =>  x = P' * inv(U) * F * b.
*
*  On entry the array x contains right-hand side vector b. On exit this
*  array contains solution vector x. */

static void solve(SCF *scf, double x[])
{     int n = scf->n;
      double *f = scf->f;
      double *u = scf->u;
      int *p = scf->p;
      double *y = scf->w;
      int i, j, ij;
      double t;
      /* y := F * b */
      for (i = 1; i <= n; i++)
      {  /* y[i] = (i-th row of F) * b */
         t = 0.0;
         for (j = 1, ij = f_loc(scf, i, 1); j <= n; j++, ij++)
            t += f[ij] * x[j];
         y[i] = t;
      }
      /* y := inv(U) * y */
      for (i = n; i >= 1; i--)
      {  t = y[i];
         for (j = n, ij = u_loc(scf, i, n); j > i; j--, ij--)
            t -= u[ij] * y[j];
         y[i] = t / u[ij];
      }
      /* x := P' * y */
      for (i = 1; i <= n; i++) x[p[i]] = y[i];
      return;
}

/***********************************************************************
*  The routine tsolve solves the transposed system C' * x = b.
*
*  From the main equation F * C = U * P it follows that:
*
*     C' * F' = P' * U',
*
*  therefore:
*
*     C' * x = b  =>  C' * F' * inv(F') * x = b  =>
*
*     P' * U' * inv(F') * x = b  =>  U' * inv(F') * x = P * b  =>
*
*     inv(F') * x = inv(U') * P * b  =>  x = F' * inv(U') * P * b.
*
*  On entry the array x contains right-hand side vector b. On exit this
*  array contains solution vector x. */

static void tsolve(SCF *scf, double x[])
{     int n = scf->n;
      double *f = scf->f;
      double *u = scf->u;
      int *p = scf->p;
      double *y = scf->w;
      int i, j, ij;
      double t;
      /* y := P * b */
      for (i = 1; i <= n; i++) y[i] = x[p[i]];
      /* y := inv(U') * y */
      for (i = 1; i <= n; i++)
      {  /* compute y[i] */
         ij = u_loc(scf, i, i);
         t = (y[i] /= u[ij]);
         /* substitute y[i] in other equations */
         for (j = i+1, ij++; j <= n; j++, ij++)
            y[j] -= u[ij] * t;
      }
      /* x := F' * y (computed as linear combination of rows of F) */
      for (j = 1; j <= n; j++) x[j] = 0.0;
      for (i = 1; i <= n; i++)
      {  t = y[i]; /* coefficient of linear combination */
         for (j = 1, ij = f_loc(scf, i, 1); j <= n; j++, ij++)
            x[j] += f[ij] * t;
      }
      return;
}

/***********************************************************************
*  NAME
*
*  scf_solve_it - solve either system C * x = b or C' * x = b
*
*  SYNOPSIS
*
*  #include "glpscf.h"
*  void scf_solve_it(SCF *scf, int tr, double x[]);
*
*  DESCRIPTION
*
*  The routine scf_solve_it solves either the system C * x = b (if tr
*  is zero) or the system C' * x = b, where C' is a matrix transposed
*  to C (if tr is non-zero). C is assumed to be non-singular.
*
*  On entry the array x should contain the right-hand side vector b in
*  locations x[1], ..., x[n], where n is the order of matrix C. On exit
*  the array x contains the solution vector x in the same locations. */

void scf_solve_it(SCF *scf, int tr, double x[])
{     if (scf->rank < scf->n)
         xfault("scf_solve_it: singular matrix\n");
      if (!tr)
         solve(scf, x);
      else
         tsolve(scf, x);
      return;
}

void scf_reset_it(SCF *scf)
{     /* reset factorization for empty matrix C */
      scf->n = scf->rank = 0;
      return;
}

/***********************************************************************
*  NAME
*
*  scf_delete_it - delete Schur complement factorization
*
*  SYNOPSIS
*
*  #include "glpscf.h"
*  void scf_delete_it(SCF *scf);
*
*  DESCRIPTION
*
*  The routine scf_delete_it deletes the specified factorization and
*  frees all the memory allocated to this object. */

void scf_delete_it(SCF *scf)
{     xfree(scf->f);
      xfree(scf->u);
      xfree(scf->p);
#if _GLPSCF_DEBUG
      xfree(scf->c);
#endif
      xfree(scf->w);
      xfree(scf);
      return;
}

/* eof */
