/* ifu.c (dense updatable IFU-factorization) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2012-2013 Free Software Foundation, Inc.
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
#include "ifu.h"

/***********************************************************************
*  ifu_expand - expand IFU-factorization
*
*  This routine expands the IFU-factorization of the matrix A according
*  to the following expansion of A:
*
*             ( A  c )
*     new A = (      )
*             ( r' d )
*
*  where c[1,...,n] is a new column, r[1,...,n] is a new row, and d is
*  a new diagonal element.
*
*  From the main equality F * A = U it follows that:
*
*     ( F  0 ) ( A  c )   ( FA  Fc )   ( U  Fc )
*     (      ) (      ) = (        ) = (       ),
*     ( 0  1 ) ( r' d )   ( r'   d )   ( r'  d )
*
*  thus,
*
*             ( F  0 )           ( U  Fc )
*     new F = (      ),  new U = (       ).
*             ( 0  1 )           ( r'  d )
*
*  Note that the resulting matrix U loses its upper triangular form due
*  to row spike r', which should be eliminated. */

void ifu_expand(IFU *ifu, double c[/*1+n*/], double r[/*1+n*/],
      double d)
{     /* non-optimized version */
      int n_max = ifu->n_max;
      int n = ifu->n;
      double *f_ = ifu->f;
      double *u_ = ifu->u;
      int i, j;
      double t;
#     define f(i,j) f_[(i)*n_max+(j)]
#     define u(i,j) u_[(i)*n_max+(j)]
      xassert(0 <= n && n < n_max);
      /* adjust indexing */
      c++, r++;
      /* set new zero column of matrix F */
      for (i = 0; i < n; i++)
         f(i,n) = 0.0;
      /* set new zero row of matrix F */
      for (j = 0; j < n; j++)
         f(n,j) = 0.0;
      /* set new unity diagonal element of matrix F */
      f(n,n) = 1.0;
      /* set new column of matrix U to vector (old F) * c */
      for (i = 0; i < n; i++)
      {  /* u[i,n] := (i-th row of old F) * c */
         t = 0.0;
         for (j = 0; j < n; j++)
            t += f(i,j) * c[j];
         u(i,n) = t;
      }
      /* set new row of matrix U to vector r */
      for (j = 0; j < n; j++)
         u(n,j) = r[j];
      /* set new diagonal element of matrix U to scalar d */
      u(n,n) = d;
      /* increase factorization order */
      ifu->n++;
#     undef f
#     undef u
      return;
}

/***********************************************************************
*  ifu_bg_update - update IFU-factorization (Bartels-Golub)
*
*  This routine updates IFU-factorization of the matrix A according to
*  its expansion (see comments to the routine ifu_expand). The routine
*  is based on the method proposed by Bartels and Golub [1].
*
*  RETURNS
*
*  0  The factorization has been successfully updated.
*
*  1  On some elimination step diagional element u[k,k] to be used as
*     pivot is too small in magnitude.
*
*  2  Diagonal element u[n,n] is too small in magnitude (at the end of
*     update).
*
*  REFERENCES
*
*  1. R.H.Bartels, G.H.Golub, "The Simplex Method of Linear Programming
*     Using LU-decomposition", Comm. ACM, 12, pp. 266-68, 1969. */

int ifu_bg_update(IFU *ifu, double c[/*1+n*/], double r[/*1+n*/],
      double d)
{     /* non-optimized version */
      int n_max = ifu->n_max;
      int n = ifu->n;
      double *f_ = ifu->f;
      double *u_ = ifu->u;
#if 1 /* FIXME */
      double tol = 1e-5;
#endif
      int j, k;
      double t;
#     define f(i,j) f_[(i)*n_max+(j)]
#     define u(i,j) u_[(i)*n_max+(j)]
      /* expand factorization */
      ifu_expand(ifu, c, r, d);
      /* NOTE: n keeps its old value */
      /* eliminate spike (non-zero subdiagonal elements) in last row of
       * matrix U */
      for (k = 0; k < n; k++)
      {  /* if |u[k,k]| < |u[n,k]|, interchange k-th and n-th rows to
          * provide |u[k,k]| >= |u[n,k]| for numeric stability */
         if (fabs(u(k,k)) < fabs(u(n,k)))
         {  /* interchange k-th and n-th rows of matrix U */
            for (j = k; j <= n; j++)
               t = u(k,j), u(k,j) = u(n,j), u(n,j) = t;
            /* interchange k-th and n-th rows of matrix F to keep the
             * main equality F * A = U */
            for (j = 0; j <= n; j++)
               t = f(k,j), f(k,j) = f(n,j), f(n,j) = t;
         }
         /* now |u[k,k]| >= |u[n,k]| */
         /* check if diagonal element u[k,k] can be used as pivot */
         if (fabs(u(k,k)) < tol)
         {  /* u[k,k] is too small in magnitude */
            return 1;
         }
         /* if u[n,k] = 0, elimination is not needed */
         if (u(n,k) == 0.0)
            continue;
         /* compute gaussian multiplier t = u[n,k] / u[k,k] */
         t = u(n,k) / u(k,k);
         /* apply gaussian transformation to eliminate u[n,k] */
         /* (n-th row of U) := (n-th row of U) - t * (k-th row of U) */
         for (j = k+1; j <= n; j++)
            u(n,j) -= t * u(k,j);
         /* apply the same transformation to matrix F to keep the main
          * equality F * A = U */
         for (j = 0; j <= n; j++)
            f(n,j) -= t * f(k,j);
      }
      /* now matrix U is upper triangular */
      if (fabs(u(n,n)) < tol)
      {  /* u[n,n] is too small in magnitude */
         return 2;
      }
#     undef f
#     undef u
      return 0;
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
{     /* non-optimized version */
      double t;
      if (b == 0.0)
         (*c) = 1.0, (*s) = 0.0;
      else if (fabs(a) <= fabs(b))
         t = - a / b, (*s) = 1.0 / sqrt(1.0 + t * t), (*c) = (*s) * t;
      else
         t = - b / a, (*c) = 1.0 / sqrt(1.0 + t * t), (*s) = (*c) * t;
      return;
}

/***********************************************************************
*  ifu_gr_update - update IFU-factorization (Givens rotations)
*
*  This routine updates IFU-factorization of the matrix A according to
*  its expansion (see comments to the routine ifu_expand). The routine
*  is based on Givens plane rotations [1].
*
*  RETURNS
*
*  0  The factorization has been successfully updated.
*
*  1  On some elimination step both elements u[k,k] and u[n,k] are too
*     small in magnitude.
*
*  2  Diagonal element u[n,n] is too small in magnitude (at the end of
*     update).
*
*  REFERENCES
*
*  1. G.H.Golub, C.F.Van Loan, "Matrix Computations", 2nd ed. */

int ifu_gr_update(IFU *ifu, double c[/*1+n*/], double r[/*1+n*/],
      double d)
{     /* non-optimized version */
      int n_max = ifu->n_max;
      int n = ifu->n;
      double *f_ = ifu->f;
      double *u_ = ifu->u;
#if 1 /* FIXME */
      double tol = 1e-5;
#endif
      int j, k;
      double cs, sn;
#     define f(i,j) f_[(i)*n_max+(j)]
#     define u(i,j) u_[(i)*n_max+(j)]
      /* expand factorization */
      ifu_expand(ifu, c, r, d);
      /* NOTE: n keeps its old value */
      /* eliminate spike (non-zero subdiagonal elements) in last row of
       * matrix U */
      for (k = 0; k < n; k++)
      {  /* check if elements u[k,k] and u[n,k] are eligible */
         if (fabs(u(k,k)) < tol && fabs(u(n,k)) < tol)
         {  /* both u[k,k] and u[n,k] are too small in magnitude */
            return 1;
         }
         /* if u[n,k] = 0, elimination is not needed */
         if (u(n,k) == 0.0)
            continue;
         /* compute parameters of Givens plane rotation */
         givens(u(k,k), u(n,k), &cs, &sn);
         /* apply Givens rotation to k-th and n-th rows of matrix U to
          * eliminate u[n,k] */
         for (j = k; j <= n; j++)
         {  double ukj = u(k,j), unj = u(n,j);
            u(k,j) = cs * ukj - sn * unj;
            u(n,j) = sn * ukj + cs * unj;
         }
         /* apply the same transformation to matrix F to keep the main
          * equality F * A = U */
         for (j = 0; j <= n; j++)
         {  double fkj = f(k,j), fnj = f(n,j);
            f(k,j) = cs * fkj - sn * fnj;
            f(n,j) = sn * fkj + cs * fnj;
         }
      }
      /* now matrix U is upper triangular */
      if (fabs(u(n,n)) < tol)
      {  /* u[n,n] is too small in magnitude */
         return 2;
      }
#     undef f
#     undef u
      return 0;
}

/***********************************************************************
*  ifu_a_solve - solve system A * x = b
*
*  This routine solves the system A * x = b, where the matrix A is
*  specified by its IFU-factorization.
*
*  Using the main equality F * A = U we have:
*
*     A * x = b  =>  F * A * x = F * b  =>  U * x = F * b  =>
*
*     x = inv(U) * F * b.
*
*  On entry the array x should contain elements of the right-hand side
*  vector b in locations x[1], ..., x[n], where n is the order of the
*  matrix A. On exit this array will contain elements of the solution
*  vector x in the same locations.
*
*  The working array w should have at least 1+n elements (0-th element
*  is not used). */

void ifu_a_solve(IFU *ifu, double x[/*1+n*/], double w[/*1+n*/])
{     /* non-optimized version */
      int n_max = ifu->n_max;
      int n = ifu->n;
      double *f_ = ifu->f;
      double *u_ = ifu->u;
      int i, j;
      double t;
#     define f(i,j) f_[(i)*n_max+(j)]
#     define u(i,j) u_[(i)*n_max+(j)]
      xassert(0 <= n && n <= n_max);
      /* adjust indexing */
      x++, w++;
      /* y := F * b */
      memcpy(w, x, n * sizeof(double));
      for (i = 0; i < n; i++)
      {  /* y[i] := (i-th row of F) * b */
         t = 0.0;
         for (j = 0; j < n; j++)
            t += f(i,j) * w[j];
         x[i] = t;
      }
      /* x := inv(U) * y */
      for (i = n-1; i >= 0; i--)
      {  t = x[i];
         for (j = i+1; j < n; j++)
            t -= u(i,j) * x[j];
         x[i] = t / u(i,i);
      }
#     undef f
#     undef u
      return;
}

/***********************************************************************
*  ifu_at_solve - solve system A'* x = b
*
*  This routine solves the system A'* x = b, where A' is a matrix
*  transposed to the matrix A, specified by its IFU-factorization.
*
*  Using the main equality F * A = U, from which it follows that
*  A'* F' = U', we have:
*
*     A'* x = b  =>  A'* F'* inv(F') * x = b  =>
*
*     U'* inv(F') * x = b  =>  inv(F') * x = inv(U') * b  =>
*
*     x = F' * inv(U') * b.
*
*  On entry the array x should contain elements of the right-hand side
*  vector b in locations x[1], ..., x[n], where n is the order of the
*  matrix A. On exit this array will contain elements of the solution
*  vector x in the same locations.
*
*  The working array w should have at least 1+n elements (0-th element
*  is not used). */

void ifu_at_solve(IFU *ifu, double x[/*1+n*/], double w[/*1+n*/])
{     /* non-optimized version */
      int n_max = ifu->n_max;
      int n = ifu->n;
      double *f_ = ifu->f;
      double *u_ = ifu->u;
      int i, j;
      double t;
#     define f(i,j) f_[(i)*n_max+(j)]
#     define u(i,j) u_[(i)*n_max+(j)]
      xassert(0 <= n && n <= n_max);
      /* adjust indexing */
      x++, w++;
      /* y := inv(U') * b */
      for (i = 0; i < n; i++)
      {  t = (x[i] /= u(i,i));
         for (j = i+1; j < n; j++)
            x[j] -= u(i,j) * t;
      }
      /* x := F'* y */
      for (j = 0; j < n; j++)
      {  /* x[j] := (j-th column of F) * y */
         t = 0.0;
         for (i = 0; i < n; i++)
            t += f(i,j) * x[i];
         w[j] = t;
      }
      memcpy(x, w, n * sizeof(double));
#     undef f
#     undef u
      return;
}

/* eof */
