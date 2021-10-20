/* scfint.c (interface to Schur-complement-based factorization) */

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
#include "scfint.h"

SCFINT *scfint_create(int type)
{     /* create interface to SC-factorization */
      SCFINT *fi;
      fi = talloc(1, SCFINT);
      memset(fi, 0, sizeof(SCFINT));
      switch ((fi->scf.type = type))
      {  case 1:
            fi->u.lufi = lufint_create();
            break;
         case 2:
            fi->u.btfi = btfint_create();
            break;
         default:
            xassert(type != type);
      }
      return fi;
}

int scfint_factorize(SCFINT *fi, int n, int (*col)(void *info, int j,
      int ind[], double val[]), void *info)
{     /* compute SC-factorization of specified matrix A */
      int nn_max, old_n0_max, n0_max, k, ret;
      xassert(n > 0);
      fi->valid = 0;
      /* get required value of nn_max */
      nn_max = fi->nn_max;
      if (nn_max == 0)
         nn_max = 100;
      xassert(nn_max > 0);
      /* compute factorization of specified matrix A */
      switch (fi->scf.type)
      {  case 1:
            old_n0_max = fi->u.lufi->n_max;
            fi->u.lufi->sva_n_max = 4 * n + 2 * nn_max;
            ret = lufint_factorize(fi->u.lufi, n, col, info);
            n0_max = fi->u.lufi->n_max;
            fi->scf.sva = fi->u.lufi->sva;
            fi->scf.a0.luf = fi->u.lufi->luf;
            break;
         case 2:
            old_n0_max = fi->u.btfi->n_max;
            fi->u.btfi->sva_n_max = 6 * n + 2 * nn_max;
            ret = btfint_factorize(fi->u.btfi, n, col, info);
            n0_max = fi->u.btfi->n_max;
            fi->scf.sva = fi->u.btfi->sva;
            fi->scf.a0.btf = fi->u.btfi->btf;
            break;
         default:
            xassert(fi != fi);
      }
      /* allocate/reallocate arrays, if necessary */
      if (old_n0_max < n0_max)
      {  if (fi->w1 != NULL)
            tfree(fi->w1);
         if (fi->w2 != NULL)
            tfree(fi->w2);
         if (fi->w3 != NULL)
            tfree(fi->w3);
         fi->w1 = talloc(1+n0_max, double);
         fi->w2 = talloc(1+n0_max, double);
         fi->w3 = talloc(1+n0_max, double);
      }
      if (fi->scf.nn_max != nn_max)
      {  if (fi->scf.ifu.f != NULL)
            tfree(fi->scf.ifu.f);
         if (fi->scf.ifu.u != NULL)
            tfree(fi->scf.ifu.u);
         fi->scf.ifu.f = talloc(nn_max * nn_max, double);
         fi->scf.ifu.u = talloc(nn_max * nn_max, double);
      }
      if (old_n0_max < n0_max || fi->scf.nn_max != nn_max)
      {  if (fi->scf.pp_ind != NULL)
            tfree(fi->scf.pp_ind);
         if (fi->scf.pp_inv != NULL)
            tfree(fi->scf.pp_inv);
         if (fi->scf.qq_ind != NULL)
            tfree(fi->scf.qq_ind);
         if (fi->scf.qq_inv != NULL)
            tfree(fi->scf.qq_inv);
         if (fi->w4 != NULL)
            tfree(fi->w4);
         if (fi->w5 != NULL)
            tfree(fi->w5);
         fi->scf.pp_ind = talloc(1+n0_max+nn_max, int);
         fi->scf.pp_inv = talloc(1+n0_max+nn_max, int);
         fi->scf.qq_ind = talloc(1+n0_max+nn_max, int);
         fi->scf.qq_inv = talloc(1+n0_max+nn_max, int);
         fi->w4 = talloc(1+n0_max+nn_max, double);
         fi->w5 = talloc(1+n0_max+nn_max, double);
      }
      /* initialize SC-factorization */
      fi->scf.n = n;
      fi->scf.n0 = n;
      fi->scf.nn_max = nn_max;
      fi->scf.nn = 0;
      fi->scf.rr_ref = sva_alloc_vecs(fi->scf.sva, nn_max);
      fi->scf.ss_ref = sva_alloc_vecs(fi->scf.sva, nn_max);
      fi->scf.ifu.n_max = nn_max;
      fi->scf.ifu.n = 0;
      for (k = 1; k <= n; k++)
      {  fi->scf.pp_ind[k] = k;
         fi->scf.pp_inv[k] = k;
         fi->scf.qq_ind[k] = k;
         fi->scf.qq_inv[k] = k;
      }
      /* set validation flag */
      if (ret == 0)
         fi->valid = 1;
      return ret;
}

int scfint_update(SCFINT *fi, int upd, int j, int len, const int ind[],
      const double val[])
{     /* update SC-factorization after replacing j-th column of A */
      int n = fi->scf.n;
      int n0 = fi->scf.n0;
      int nn = fi->scf.nn;
      int *pp_ind = fi->scf.pp_ind;
      int *qq_ind = fi->scf.qq_ind;
      int *qq_inv = fi->scf.qq_inv;
      double *bf = fi->w4;
      double *dg = fi->w5;
      int k, t, ret;
      xassert(fi->valid);
      xassert(0 <= n && n <= n0+nn);
      /* (b, f) := inv(P) * (beta, 0) */
      for (k = 1; k <= n0+nn; k++)
         bf[k] = 0.0;
      for (t = 1; t <= len; t++)
      {  k = ind[t];
         xassert(1 <= k && k <= n);
#if 1 /* FIXME: currently P = I */
         xassert(pp_ind[k] == k);
#endif
         xassert(bf[k] == 0.0);
         xassert(val[t] != 0.0);
         bf[k] = val[t];
      }
      /* (d, g) := Q * (cj, 0) */
      for (k = 1; k <= n0+nn; k++)
         dg[k] = 0.0;
      xassert(1 <= j && j <= n);
      dg[fi->scf.qq_inv[j]] = 1;
      /* update factorization of augmented matrix */
      ret = scf_update_aug(&fi->scf, &bf[0], &dg[0], &bf[n0], &dg[n0],
         0.0, upd, fi->w1, fi->w2, fi->w3);
      if (ret == 0)
      {  /* swap j-th and last columns of new matrix Q */
         scf_swap_q_cols(j, n0+nn+1);
      }
      else
      {  /* updating failed */
         fi->valid = 0;
      }
      return ret;
}

void scfint_ftran(SCFINT *fi, double x[])
{     /* solve system A * x = b */
      xassert(fi->valid);
      scf_a_solve(&fi->scf, x, fi->w4, fi->w5, fi->w1, fi->w2);
      return;
}

void scfint_btran(SCFINT *fi, double x[])
{     /* solve system A'* x = b */
      xassert(fi->valid);
      scf_at_solve(&fi->scf, x, fi->w4, fi->w5, fi->w1, fi->w2);
      return;
}

double scfint_estimate(SCFINT *fi)
{     /* estimate 1-norm of inv(A) */
      double norm;
      xassert(fi->valid);
      xassert(fi->scf.n == fi->scf.n0);
      switch (fi->scf.type)
      {  case 1:
            norm = luf_estimate_norm(fi->scf.a0.luf, fi->w1, fi->w2);
            break;
         case 2:
            norm = btf_estimate_norm(fi->scf.a0.btf, fi->w1, fi->w2,
               fi->w3, fi->w4);
            break;
         default:
            xassert(fi != fi);
      }
      return norm;
}

void scfint_delete(SCFINT *fi)
{     /* delete interface to SC-factorization */
      switch (fi->scf.type)
      {  case 1:
            lufint_delete(fi->u.lufi);
            break;
         case 2:
            btfint_delete(fi->u.btfi);
            break;
         default:
            xassert(fi != fi);
      }
      if (fi->scf.ifu.f != NULL)
         tfree(fi->scf.ifu.f);
      if (fi->scf.ifu.u != NULL)
         tfree(fi->scf.ifu.u);
      if (fi->scf.pp_ind != NULL)
         tfree(fi->scf.pp_ind);
      if (fi->scf.pp_inv != NULL)
         tfree(fi->scf.pp_inv);
      if (fi->scf.qq_ind != NULL)
         tfree(fi->scf.qq_ind);
      if (fi->scf.qq_inv != NULL)
         tfree(fi->scf.qq_inv);
      if (fi->w1 != NULL)
         tfree(fi->w1);
      if (fi->w2 != NULL)
         tfree(fi->w2);
      if (fi->w3 != NULL)
         tfree(fi->w3);
      if (fi->w4 != NULL)
         tfree(fi->w4);
      if (fi->w5 != NULL)
         tfree(fi->w5);
      tfree(fi);
      return;
}

/* eof */
