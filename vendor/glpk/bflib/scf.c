/* scf.c (sparse updatable Schur-complement-based factorization) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2013-2014 Free Software Foundation, Inc.
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
#include "scf.h"

/***********************************************************************
*  scf_r0_solve - solve system R0 * x = b or R0'* x = b
*
*  This routine solves the system R0 * x = b (if tr is zero) or the
*  system R0'* x = b (if tr is non-zero), where R0 is the left factor
*  of the initial matrix A0 = R0 * S0.
*
*  On entry the array x should contain elements of the right-hand side
*  vector b in locations x[1], ..., x[n0], where n0 is the order of the
*  matrix R0. On exit the array x will contain elements of the solution
*  vector in the same locations. */

void scf_r0_solve(SCF *scf, int tr, double x[/*1+n0*/])
{     switch (scf->type)
      {  case 1:
            /* A0 = F0 * V0, so R0 = F0 */
            if (!tr)
               luf_f_solve(scf->a0.luf, x);
            else
               luf_ft_solve(scf->a0.luf, x);
            break;
         case 2:
            /* A0 = I * A0, so R0 = I */
            break;
         default:
            xassert(scf != scf);
      }
      return;
}

/***********************************************************************
*  scf_s0_solve - solve system S0 * x = b or S0'* x = b
*
*  This routine solves the system S0 * x = b (if tr is zero) or the
*  system S0'* x = b (if tr is non-zero), where S0 is the right factor
*  of the initial matrix A0 = R0 * S0.
*
*  On entry the array x should contain elements of the right-hand side
*  vector b in locations x[1], ..., x[n0], where n0 is the order of the
*  matrix S0. On exit the array x will contain elements of the solution
*  vector in the same locations.
*
*  The routine uses locations [1], ..., [n0] of three working arrays
*  w1, w2, and w3. (In case of type = 1 arrays w2 and w3 are not used
*  and can be specified as NULL.) */

void scf_s0_solve(SCF *scf, int tr, double x[/*1+n0*/],
      double w1[/*1+n0*/], double w2[/*1+n0*/], double w3[/*1+n0*/])
{     int n0 = scf->n0;
      switch (scf->type)
      {  case 1:
            /* A0 = F0 * V0, so S0 = V0 */
            if (!tr)
               luf_v_solve(scf->a0.luf, x, w1);
            else
               luf_vt_solve(scf->a0.luf, x, w1);
            break;
         case 2:
            /* A0 = I * A0, so S0 = A0 */
            if (!tr)
               btf_a_solve(scf->a0.btf, x, w1, w2, w3);
            else
               btf_at_solve(scf->a0.btf, x, w1, w2, w3);
            break;
         default:
            xassert(scf != scf);
      }
      memcpy(&x[1], &w1[1], n0 * sizeof(double));
      return;
}

/***********************************************************************
*  scf_r_prod - compute product y := y + alpha * R * x
*
*  This routine computes the product y := y + alpha * R * x, where
*  x is a n0-vector, alpha is a scalar, y is a nn-vector.
*
*  Since matrix R is available by rows, the product components are
*  computed as inner products:
*
*     y[i] = y[i] + alpha * (i-th row of R) * x
*
*  for i = 1, 2, ..., nn. */

void scf_r_prod(SCF *scf, double y[/*1+nn*/], double a, const double
      x[/*1+n0*/])
{     int nn = scf->nn;
      SVA *sva = scf->sva;
      int *sv_ind = sva->ind;
      double *sv_val = sva->val;
      int rr_ref = scf->rr_ref;
      int *rr_ptr = &sva->ptr[rr_ref-1];
      int *rr_len = &sva->len[rr_ref-1];
      int i, ptr, end;
      double t;
      for (i = 1; i <= nn; i++)
      {  /* t := (i-th row of R) * x */
         t = 0.0;
         for (end = (ptr = rr_ptr[i]) + rr_len[i]; ptr < end; ptr++)
            t += sv_val[ptr] * x[sv_ind[ptr]];
         /* y[i] := y[i] + alpha * t */
         y[i] += a * t;
      }
      return;
}

/***********************************************************************
*  scf_rt_prod - compute product y := y + alpha * R'* x
*
*  This routine computes the product y := y + alpha * R'* x, where
*  R' is a matrix transposed to R, x is a nn-vector, alpha is a scalar,
*  y is a n0-vector.
*
*  Since matrix R is available by rows, the product is computed as a
*  linear combination:
*
*     y := y + alpha * (R'[1] * x[1] + ... + R'[nn] * x[nn]),
*
*  where R'[i] is i-th row of R. */

void scf_rt_prod(SCF *scf, double y[/*1+n0*/], double a, const double
      x[/*1+nn*/])
{     int nn = scf->nn;
      SVA *sva = scf->sva;
      int *sv_ind = sva->ind;
      double *sv_val = sva->val;
      int rr_ref = scf->rr_ref;
      int *rr_ptr = &sva->ptr[rr_ref-1];
      int *rr_len = &sva->len[rr_ref-1];
      int i, ptr, end;
      double t;
      for (i = 1; i <= nn; i++)
      {  if (x[i] == 0.0)
            continue;
         /* y := y + alpha * R'[i] * x[i] */
         t = a * x[i];
         for (end = (ptr = rr_ptr[i]) + rr_len[i]; ptr < end; ptr++)
            y[sv_ind[ptr]] += sv_val[ptr] * t;
      }
      return;
}

/***********************************************************************
*  scf_s_prod - compute product y := y + alpha * S * x
*
*  This routine computes the product y := y + alpha * S * x, where
*  x is a nn-vector, alpha is a scalar, y is a n0 vector.
*
*  Since matrix S is available by columns, the product is computed as
*  a linear combination:
*
*     y := y + alpha * (S[1] * x[1] + ... + S[nn] * x[nn]),
*
*  where S[j] is j-th column of S. */

void scf_s_prod(SCF *scf, double y[/*1+n0*/], double a, const double
      x[/*1+nn*/])
{     int nn = scf->nn;
      SVA *sva = scf->sva;
      int *sv_ind = sva->ind;
      double *sv_val = sva->val;
      int ss_ref = scf->ss_ref;
      int *ss_ptr = &sva->ptr[ss_ref-1];
      int *ss_len = &sva->len[ss_ref-1];
      int j, ptr, end;
      double t;
      for (j = 1; j <= nn; j++)
      {  if (x[j] == 0.0)
            continue;
         /* y := y + alpha * S[j] * x[j] */
         t = a * x[j];
         for (end = (ptr = ss_ptr[j]) + ss_len[j]; ptr < end; ptr++)
            y[sv_ind[ptr]] += sv_val[ptr] * t;
      }
      return;
}

/***********************************************************************
*  scf_st_prod - compute product y := y + alpha * S'* x
*
*  This routine computes the product y := y + alpha * S'* x, where
*  S' is a matrix transposed to S, x is a n0-vector, alpha is a scalar,
*  y is a nn-vector.
*
*  Since matrix S is available by columns, the product components are
*  computed as inner products:
*
*     y[j] := y[j] + alpha * (j-th column of S) * x
*
*  for j = 1, 2, ..., nn. */

void scf_st_prod(SCF *scf, double y[/*1+nn*/], double a, const double
      x[/*1+n0*/])
{     int nn = scf->nn;
      SVA *sva = scf->sva;
      int *sv_ind = sva->ind;
      double *sv_val = sva->val;
      int ss_ref = scf->ss_ref;
      int *ss_ptr = &sva->ptr[ss_ref-1];
      int *ss_len = &sva->len[ss_ref-1];
      int j, ptr, end;
      double t;
      for (j = 1; j <= nn; j++)
      {  /* t := (j-th column of S) * x */
         t = 0.0;
         for (end = (ptr = ss_ptr[j]) + ss_len[j]; ptr < end; ptr++)
            t += sv_val[ptr] * x[sv_ind[ptr]];
         /* y[j] := y[j] + alpha * t */
         y[j] += a * t;
      }
      return;
}

/***********************************************************************
*  scf_a_solve - solve system A * x = b
*
*  This routine solves the system A * x = b, where A is the current
*  matrix.
*
*  On entry the array x should contain elements of the right-hand side
*  vector b in locations x[1], ..., x[n], where n is the order of the
*  matrix A. On exit the array x will contain elements of the solution
*  vector in the same locations.
*
*  For details see the program documentation. */

void scf_a_solve(SCF *scf, double x[/*1+n*/],
      double w[/*1+n0+nn*/], double work1[/*1+max(n0,nn)*/],
      double work2[/*1+n*/], double work3[/*1+n*/])
{     int n = scf->n;
      int n0 = scf->n0;
      int nn = scf->nn;
      int *pp_ind = scf->pp_ind;
      int *qq_inv = scf->qq_inv;
      int i, ii;
      /* (u1, u2) := inv(P) * (b, 0) */
      for (ii = 1; ii <= n0+nn; ii++)
      {  i = pp_ind[ii];
#if 1 /* FIXME: currently P = I */
         xassert(i == ii);
#endif
         w[ii] = (i <= n ? x[i] : 0.0);
      }
      /* v1 := inv(R0) * u1 */
      scf_r0_solve(scf, 0, &w[0]);
      /* v2 := u2 - R * v1 */
      scf_r_prod(scf, &w[n0], -1.0, &w[0]);
      /* w2 := inv(C) * v2 */
      ifu_a_solve(&scf->ifu, &w[n0], work1);
      /* w1 := inv(S0) * (v1 - S * w2) */
      scf_s_prod(scf, &w[0], -1.0, &w[n0]);
      scf_s0_solve(scf, 0, &w[0], work1, work2, work3);
      /* (x, x~) := inv(Q) * (w1, w2); x~ is not needed */
      for (i = 1; i <= n; i++)
         x[i] = w[qq_inv[i]];
      return;
}

/***********************************************************************
*  scf_at_solve - solve system A'* x = b
*
*  This routine solves the system A'* x = b, where A' is a matrix
*  transposed to the current matrix A.
*
*  On entry the array x should contain elements of the right-hand side
*  vector b in locations x[1], ..., x[n], where n is the order of the
*  matrix A. On exit the array x will contain elements of the solution
*  vector in the same locations.
*
*  For details see the program documentation. */

void scf_at_solve(SCF *scf, double x[/*1+n*/],
      double w[/*1+n0+nn*/], double work1[/*1+max(n0,nn)*/],
      double work2[/*1+n*/], double work3[/*1+n*/])
{     int n = scf->n;
      int n0 = scf->n0;
      int nn = scf->nn;
      int *pp_inv = scf->pp_inv;
      int *qq_ind = scf->qq_ind;
      int i, ii;
      /* (u1, u2) := Q * (b, 0) */
      for (ii = 1; ii <= n0+nn; ii++)
      {  i = qq_ind[ii];
         w[ii] = (i <= n ? x[i] : 0.0);
      }
      /* v1 := inv(S0') * u1 */
      scf_s0_solve(scf, 1, &w[0], work1, work2, work3);
      /* v2 := inv(C') * (u2 - S'* v1) */
      scf_st_prod(scf, &w[n0], -1.0, &w[0]);
      ifu_at_solve(&scf->ifu, &w[n0], work1);
      /* w2 := v2 */
      /* nop */
      /* w1 := inv(R0') * (v1 - R'* w2) */
      scf_rt_prod(scf, &w[0], -1.0, &w[n0]);
      scf_r0_solve(scf, 1, &w[0]);
      /* compute (x, x~) := P * (w1, w2); x~ is not needed */
      for (i = 1; i <= n; i++)
      {
#if 1 /* FIXME: currently P = I */
         xassert(pp_inv[i] == i);
#endif
         x[i] = w[pp_inv[i]];
      }
      return;
}

/***********************************************************************
*  scf_add_r_row - add new row to matrix R
*
*  This routine adds new (nn+1)-th row to matrix R, whose elements are
*  specified in locations w[1,...,n0]. */

void scf_add_r_row(SCF *scf, const double w[/*1+n0*/])
{     int n0 = scf->n0;
      int nn = scf->nn;
      SVA *sva = scf->sva;
      int *sv_ind = sva->ind;
      double *sv_val = sva->val;
      int rr_ref = scf->rr_ref;
      int *rr_ptr = &sva->ptr[rr_ref-1];
      int *rr_len = &sva->len[rr_ref-1];
      int j, len, ptr;
      xassert(0 <= nn && nn < scf->nn_max);
      /* determine length of new row */
      len = 0;
      for (j = 1; j <= n0; j++)
      {  if (w[j] != 0.0)
            len++;
      }
      /* reserve locations for new row in static part of SVA */
      if (len > 0)
      {  if (sva->r_ptr - sva->m_ptr < len)
         {  sva_more_space(sva, len);
            sv_ind = sva->ind;
            sv_val = sva->val;
         }
         sva_reserve_cap(sva, rr_ref + nn, len);
      }
      /* store new row in sparse format */
      ptr = rr_ptr[nn+1];
      for (j = 1; j <= n0; j++)
      {  if (w[j] != 0.0)
         {  sv_ind[ptr] = j;
            sv_val[ptr] = w[j];
            ptr++;
         }
      }
      xassert(ptr - rr_ptr[nn+1] == len);
      rr_len[nn+1] = len;
#ifdef GLP_DEBUG
      sva_check_area(sva);
#endif
      return;
}

/***********************************************************************
*  scf_add_s_col - add new column to matrix S
*
*  This routine adds new (nn+1)-th column to matrix S, whose elements
*  are specified in locations v[1,...,n0]. */

void scf_add_s_col(SCF *scf, const double v[/*1+n0*/])
{     int n0 = scf->n0;
      int nn = scf->nn;
      SVA *sva = scf->sva;
      int *sv_ind = sva->ind;
      double *sv_val = sva->val;
      int ss_ref = scf->ss_ref;
      int *ss_ptr = &sva->ptr[ss_ref-1];
      int *ss_len = &sva->len[ss_ref-1];
      int i, len, ptr;
      xassert(0 <= nn && nn < scf->nn_max);
      /* determine length of new column */
      len = 0;
      for (i = 1; i <= n0; i++)
      {  if (v[i] != 0.0)
            len++;
      }
      /* reserve locations for new column in static part of SVA */
      if (len > 0)
      {  if (sva->r_ptr - sva->m_ptr < len)
         {  sva_more_space(sva, len);
            sv_ind = sva->ind;
            sv_val = sva->val;
         }
         sva_reserve_cap(sva, ss_ref + nn, len);
      }
      /* store new column in sparse format */
      ptr = ss_ptr[nn+1];
      for (i = 1; i <= n0; i++)
      {  if (v[i] != 0.0)
         {  sv_ind[ptr] = i;
            sv_val[ptr] = v[i];
            ptr++;
         }
      }
      xassert(ptr - ss_ptr[nn+1] == len);
      ss_len[nn+1] = len;
#ifdef GLP_DEBUG
      sva_check_area(sva);
#endif
      return;
}

/***********************************************************************
*  scf_update_aug - update factorization of augmented matrix
*
*  Given factorization of the current augmented matrix:
*
*     ( A0  A1 )   ( R0    ) ( S0  S )
*     (        ) = (       ) (       ),
*     ( A2  A3 )   ( R   I ) (     C )
*
*  this routine computes factorization of the new augmented matrix:
*
*     ( A0 | A1  b )
*     ( ---+------ )   ( A0  A1^ )   ( R0    ) ( S0  S^ )
*     ( A2 | A3  f ) = (         ) = (       ) (        ),
*     (    |       )   ( A2^ A3^ )   ( R^  I ) (     C^ )
*     ( d' | g'  h )
*
*  where b and d are specified n0-vectors, f and g are specified
*  nn-vectors, and h is a specified scalar. (Note that corresponding
*  arrays are clobbered on exit.)
*
*  The parameter upd specifies how to update factorization of the Schur
*  complement C:
*
*  1  Bartels-Golub updating.
*
*  2  Givens rotations updating.
*
*  The working arrays w1, w2, and w3 are used in the same way as in the
*  routine scf_s0_solve.
*
*  RETURNS
*
*  0  Factorization has been successfully updated.
*
*  1  Updating limit has been reached.
*
*  2  Updating IFU-factorization of matrix C failed.
*
*  For details see the program documentation. */

int scf_update_aug(SCF *scf, double b[/*1+n0*/], double d[/*1+n0*/],
      double f[/*1+nn*/], double g[/*1+nn*/], double h, int upd,
      double w1[/*1+n0*/], double w2[/*1+n0*/], double w3[/*1+n0*/])
{     int n0 = scf->n0;
      int k, ret;
      double *v, *w, *x, *y, z;
      if (scf->nn == scf->nn_max)
      {  /* updating limit has been reached */
         return 1;
      }
      /* v := inv(R0) * b */
      scf_r0_solve(scf, 0, (v = b));
      /* w := inv(S0') * d */
      scf_s0_solve(scf, 1, (w = d), w1, w2, w3);
      /* x := f - R * v */
      scf_r_prod(scf, (x = f), -1.0, v);
      /* y := g - S'* w */
      scf_st_prod(scf, (y = g), -1.0, w);
      /* z := h - v'* w */
      z = h;
      for (k = 1; k <= n0; k++)
         z -= v[k] * w[k];
      /* new R := R with row w added */
      scf_add_r_row(scf, w);
      /* new S := S with column v added */
      scf_add_s_col(scf, v);
      /* update IFU-factorization of C */
      switch (upd)
      {  case 1:
            ret = ifu_bg_update(&scf->ifu, x, y, z);
            break;
         case 2:
            ret = ifu_gr_update(&scf->ifu, x, y, z);
            break;
         default:
            xassert(upd != upd);
      }
      if (ret != 0)
      {  /* updating IFU-factorization failed */
         return 2;
      }
      /* increase number of additional rows and columns */
      scf->nn++;
      /* expand P and Q */
      k = n0 + scf->nn;
      scf->pp_ind[k] = scf->pp_inv[k] = k;
      scf->qq_ind[k] = scf->qq_inv[k] = k;
      /* factorization has been successfully updated */
      return 0;
}

/* eof */
