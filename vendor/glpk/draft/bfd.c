/* bfd.c (LP basis factorization driver) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2007-2014 Free Software Foundation, Inc.
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

#include "glpk.h"
#include "env.h"
#include "bfd.h"
#include "fhvint.h"
#include "scfint.h"
#ifdef GLP_DEBUG
#include "glpspm.h"
#endif

struct BFD
{     /* LP basis factorization driver */
      int valid;
      /* factorization is valid only if this flag is set */
      int type;
      /* type of factorization used:
         0 - interface not established yet
         1 - FHV-factorization
         2 - Schur-complement-based factorization */
      union
      {  void *none;   /* type = 0 */
         FHVINT *fhvi; /* type = 1 */
         SCFINT *scfi; /* type = 2 */
      }  u;
      /* interface to factorization of LP basis */
      glp_bfcp parm;
      /* factorization control parameters */
#ifdef GLP_DEBUG
      SPM *B;
      /* current basis (for testing/debugging only) */
#endif
      int upd_cnt;
      /* factorization update count */
#if 1 /* 21/IV-2014 */
      double b_norm;
      /* 1-norm of matrix B */
      double i_norm;
      /* estimated 1-norm of matrix inv(B) */
#endif
};

BFD *bfd_create_it(void)
{     /* create LP basis factorization */
      BFD *bfd;
#ifdef GLP_DEBUG
      xprintf("bfd_create_it: warning: debugging version used\n");
#endif
      bfd = talloc(1, BFD);
      bfd->valid = 0;
      bfd->type = 0;
      bfd->u.none = NULL;
      bfd_set_bfcp(bfd, NULL);
#ifdef GLP_DEBUG
      bfd->B = NULL;
#endif
      bfd->upd_cnt = 0;
      return bfd;
}

#if 0 /* 08/III-2014 */
void bfd_set_parm(BFD *bfd, const void *parm)
{     /* change LP basis factorization control parameters */
      memcpy(&bfd->parm, parm, sizeof(glp_bfcp));
      return;
}
#endif

void bfd_get_bfcp(BFD *bfd, void /* glp_bfcp */ *parm)
{     /* retrieve LP basis factorization control parameters */
      memcpy(parm, &bfd->parm, sizeof(glp_bfcp));
      return;
}

void bfd_set_bfcp(BFD *bfd, const void /* glp_bfcp */ *parm)
{     /* change LP basis factorization control parameters */
      if (parm == NULL)
      {  /* reset to default */
         memset(&bfd->parm, 0, sizeof(glp_bfcp));
         bfd->parm.type = GLP_BF_LUF + GLP_BF_FT;
         bfd->parm.piv_tol = 0.10;
         bfd->parm.piv_lim = 4;
         bfd->parm.suhl = 1;
         bfd->parm.eps_tol = DBL_EPSILON;
         bfd->parm.nfs_max = 100;
         bfd->parm.nrs_max = 70;
      }
      else
         memcpy(&bfd->parm, parm, sizeof(glp_bfcp));
      return;
}

#if 1 /* 21/IV-2014 */
struct bfd_info
{     BFD *bfd;
      int (*col)(void *info, int j, int ind[], double val[]);
      void *info;
};

static int bfd_col(void *info_, int j, int ind[], double val[])
{     struct bfd_info *info = info_;
      int t, len;
      double sum;
      len = info->col(info->info, j, ind, val);
      sum = 0.0;
      for (t = 1; t <= len; t++)
      {  if (val[t] >= 0.0)
            sum += val[t];
         else
            sum -= val[t];
      }
      if (info->bfd->b_norm < sum)
         info->bfd->b_norm = sum;
      return len;
}
#endif

int bfd_factorize(BFD *bfd, int m, /*const int bh[],*/ int (*col1)
      (void *info, int j, int ind[], double val[]), void *info1)
{     /* compute LP basis factorization */
#if 1 /* 21/IV-2014 */
      struct bfd_info info;
#endif
      int type, ret;
      /*xassert(bh == bh);*/
      /* invalidate current factorization */
      bfd->valid = 0;
      /* determine required factorization type */
      switch (bfd->parm.type)
      {  case GLP_BF_LUF + GLP_BF_FT:
            type = 1;
            break;
         case GLP_BF_LUF + GLP_BF_BG:
         case GLP_BF_LUF + GLP_BF_GR:
         case GLP_BF_BTF + GLP_BF_BG:
         case GLP_BF_BTF + GLP_BF_GR:
            type = 2;
            break;
         default:
            xassert(bfd != bfd);
      }
      /* delete factorization interface, if necessary */
      switch (bfd->type)
      {  case 0:
            break;
         case 1:
            if (type != 1)
            {  bfd->type = 0;
               fhvint_delete(bfd->u.fhvi);
               bfd->u.fhvi = NULL;
            }
            break;
         case 2:
            if (type != 2)
            {  bfd->type = 0;
               scfint_delete(bfd->u.scfi);
               bfd->u.scfi = NULL;
            }
            break;
         default:
            xassert(bfd != bfd);
      }
      /* establish factorization interface, if necessary */
      if (bfd->type == 0)
      {  switch (type)
         {  case 1:
               bfd->type = 1;
               xassert(bfd->u.fhvi == NULL);
               bfd->u.fhvi = fhvint_create();
               break;
            case 2:
               bfd->type = 2;
               xassert(bfd->u.scfi == NULL);
               if (!(bfd->parm.type & GLP_BF_BTF))
                  bfd->u.scfi = scfint_create(1);
               else
                  bfd->u.scfi = scfint_create(2);
               break;
            default:
               xassert(type != type);
         }
      }
      /* try to compute factorization */
#if 1 /* 21/IV-2014 */
      bfd->b_norm = bfd->i_norm = 0.0;
      info.bfd = bfd;
      info.col = col1;
      info.info = info1;
#endif
      switch (bfd->type)
      {  case 1:
            bfd->u.fhvi->lufi->sgf_piv_tol = bfd->parm.piv_tol;
            bfd->u.fhvi->lufi->sgf_piv_lim = bfd->parm.piv_lim;
            bfd->u.fhvi->lufi->sgf_suhl = bfd->parm.suhl;
            bfd->u.fhvi->lufi->sgf_eps_tol = bfd->parm.eps_tol;
            bfd->u.fhvi->nfs_max = bfd->parm.nfs_max;
            ret = fhvint_factorize(bfd->u.fhvi, m, bfd_col, &info);
#if 1 /* FIXME */
            if (ret == 0)
               bfd->i_norm = fhvint_estimate(bfd->u.fhvi);
            else
               ret = BFD_ESING;
#endif
            break;
         case 2:
            if (bfd->u.scfi->scf.type == 1)
            {  bfd->u.scfi->u.lufi->sgf_piv_tol = bfd->parm.piv_tol;
               bfd->u.scfi->u.lufi->sgf_piv_lim = bfd->parm.piv_lim;
               bfd->u.scfi->u.lufi->sgf_suhl = bfd->parm.suhl;
               bfd->u.scfi->u.lufi->sgf_eps_tol = bfd->parm.eps_tol;
            }
            else if (bfd->u.scfi->scf.type == 2)
            {  bfd->u.scfi->u.btfi->sgf_piv_tol = bfd->parm.piv_tol;
               bfd->u.scfi->u.btfi->sgf_piv_lim = bfd->parm.piv_lim;
               bfd->u.scfi->u.btfi->sgf_suhl = bfd->parm.suhl;
               bfd->u.scfi->u.btfi->sgf_eps_tol = bfd->parm.eps_tol;
            }
            else
               xassert(bfd != bfd);
            bfd->u.scfi->nn_max = bfd->parm.nrs_max;
            ret = scfint_factorize(bfd->u.scfi, m, bfd_col, &info);
#if 1 /* FIXME */
            if (ret == 0)
               bfd->i_norm = scfint_estimate(bfd->u.scfi);
            else
               ret = BFD_ESING;
#endif
            break;
         default:
            xassert(bfd != bfd);
      }
#ifdef GLP_DEBUG
      /* save specified LP basis */
      if (bfd->B != NULL)
         spm_delete_mat(bfd->B);
      bfd->B = spm_create_mat(m, m);
      {  int *ind = talloc(1+m, int);
         double *val = talloc(1+m, double);
         int j, k, len;
         for (j = 1; j <= m; j++)
         {  len = col(info, j, ind, val);
            for (k = 1; k <= len; k++)
               spm_new_elem(bfd->B, ind[k], j, val[k]);
         }
         tfree(ind);
         tfree(val);
      }
#endif
      if (ret == 0)
      {  /* factorization has been successfully computed */
         double cond;
         bfd->valid = 1;
#ifdef GLP_DEBUG
         cond = bfd_condest(bfd);
         if (cond > 1e9)
            xprintf("bfd_factorize: warning: cond(B) = %g\n", cond);
#endif
      }
#ifdef GLP_DEBUG
      xprintf("bfd_factorize: m = %d; ret = %d\n", m, ret);
#endif
      bfd->upd_cnt = 0;
      return ret;
}

#if 0 /* 21/IV-2014 */
double bfd_estimate(BFD *bfd)
{     /* estimate 1-norm of inv(B) */
      double norm;
      xassert(bfd->valid);
      xassert(bfd->upd_cnt == 0);
      switch (bfd->type)
      {  case 1:
            norm = fhvint_estimate(bfd->u.fhvi);
            break;
         case 2:
            norm = scfint_estimate(bfd->u.scfi);
            break;
         default:
            xassert(bfd != bfd);
      }
      return norm;
}
#endif

#if 1 /* 21/IV-2014 */
double bfd_condest(BFD *bfd)
{     /* estimate condition of B */
      double cond;
      xassert(bfd->valid);
      /*xassert(bfd->upd_cnt == 0);*/
      cond = bfd->b_norm * bfd->i_norm;
      if (cond < 1.0)
         cond = 1.0;
      return cond;
}
#endif

void bfd_ftran(BFD *bfd, double x[])
{     /* perform forward transformation (solve system B * x = b) */
#ifdef GLP_DEBUG
      SPM *B = bfd->B;
      int m = B->m;
      double *b = talloc(1+m, double);
      SPME *e;
      int k;
      double s, relerr, maxerr;
      for (k = 1; k <= m; k++)
         b[k] = x[k];
#endif
      xassert(bfd->valid);
      switch (bfd->type)
      {  case 1:
            fhvint_ftran(bfd->u.fhvi, x);
            break;
         case 2:
            scfint_ftran(bfd->u.scfi, x);
            break;
         default:
            xassert(bfd != bfd);
      }
#ifdef GLP_DEBUG
      maxerr = 0.0;
      for (k = 1; k <= m; k++)
      {  s = 0.0;
         for (e = B->row[k]; e != NULL; e = e->r_next)
            s += e->val * x[e->j];
         relerr = (b[k] - s) / (1.0 + fabs(b[k]));
         if (maxerr < relerr)
            maxerr = relerr;
      }
      if (maxerr > 1e-8)
         xprintf("bfd_ftran: maxerr = %g; relative error too large\n",
            maxerr);
      tfree(b);
#endif
      return;
}

#if 1 /* 30/III-2016 */
void bfd_ftran_s(BFD *bfd, FVS *x)
{     /* sparse version of bfd_ftran */
      /* (sparse mode is not implemented yet) */
      int n = x->n;
      int *ind = x->ind;
      double *vec = x->vec;
      int j, nnz = 0;
      bfd_ftran(bfd, vec);
      for (j = n; j >= 1; j--)
      {  if (vec[j] != 0.0)
            ind[++nnz] = j;
      }
      x->nnz = nnz;
      return;
}
#endif

void bfd_btran(BFD *bfd, double x[])
{     /* perform backward transformation (solve system B'* x = b) */
#ifdef GLP_DEBUG
      SPM *B = bfd->B;
      int m = B->m;
      double *b = talloc(1+m, double);
      SPME *e;
      int k;
      double s, relerr, maxerr;
      for (k = 1; k <= m; k++)
         b[k] = x[k];
#endif
      xassert(bfd->valid);
      switch (bfd->type)
      {  case 1:
            fhvint_btran(bfd->u.fhvi, x);
            break;
         case 2:
            scfint_btran(bfd->u.scfi, x);
            break;
         default:
            xassert(bfd != bfd);
      }
#ifdef GLP_DEBUG
      maxerr = 0.0;
      for (k = 1; k <= m; k++)
      {  s = 0.0;
         for (e = B->col[k]; e != NULL; e = e->c_next)
            s += e->val * x[e->i];
         relerr = (b[k] - s) / (1.0 + fabs(b[k]));
         if (maxerr < relerr)
            maxerr = relerr;
      }
      if (maxerr > 1e-8)
         xprintf("bfd_btran: maxerr = %g; relative error too large\n",
            maxerr);
      tfree(b);
#endif
      return;
}

#if 1 /* 30/III-2016 */
void bfd_btran_s(BFD *bfd, FVS *x)
{     /* sparse version of bfd_btran */
      /* (sparse mode is not implemented yet) */
      int n = x->n;
      int *ind = x->ind;
      double *vec = x->vec;
      int j, nnz = 0;
      bfd_btran(bfd, vec);
      for (j = n; j >= 1; j--)
      {  if (vec[j] != 0.0)
            ind[++nnz] = j;
      }
      x->nnz = nnz;
      return;
}
#endif

int bfd_update(BFD *bfd, int j, int len, const int ind[], const double
      val[])
{     /* update LP basis factorization */
      int ret;
      xassert(bfd->valid);
      switch (bfd->type)
      {  case 1:
            ret = fhvint_update(bfd->u.fhvi, j, len, ind, val);
#if 1 /* FIXME */
            switch (ret)
            {  case 0:
                  break;
               case 1:
                  ret = BFD_ESING;
                  break;
               case 2:
               case 3:
                  ret = BFD_ECOND;
                  break;
               case 4:
                  ret = BFD_ELIMIT;
                  break;
               case 5:
                  ret = BFD_ECHECK;
                  break;
               default:
                  xassert(ret != ret);
            }
#endif
            break;
         case 2:
            switch (bfd->parm.type & 0x0F)
            {  case GLP_BF_BG:
                  ret = scfint_update(bfd->u.scfi, 1, j, len, ind, val);
                  break;
               case GLP_BF_GR:
                  ret = scfint_update(bfd->u.scfi, 2, j, len, ind, val);
                  break;
               default:
                  xassert(bfd != bfd);
            }
#if 1 /* FIXME */
            switch (ret)
            {  case 0:
                  break;
               case 1:
                  ret = BFD_ELIMIT;
                  break;
               case 2:
                  ret = BFD_ECOND;
                  break;
               default:
                  xassert(ret != ret);
            }
#endif
            break;
         default:
            xassert(bfd != bfd);
      }
      if (ret != 0)
      {  /* updating factorization failed */
         bfd->valid = 0;
      }
#ifdef GLP_DEBUG
      /* save updated LP basis */
      {  SPME *e;
         int k;
         for (e = bfd->B->col[j]; e != NULL; e = e->c_next)
            e->val = 0.0;
         spm_drop_zeros(bfd->B, 0.0);
         for (k = 1; k <= len; k++)
            spm_new_elem(bfd->B, ind[k], j, val[k]);
      }
#endif
      if (ret == 0)
         bfd->upd_cnt++;
      return ret;
}

int bfd_get_count(BFD *bfd)
{     /* determine factorization update count */
      return bfd->upd_cnt;
}

void bfd_delete_it(BFD *bfd)
{     /* delete LP basis factorization */
      switch (bfd->type)
      {  case 0:
            break;
         case 1:
            fhvint_delete(bfd->u.fhvi);
            break;
         case 2:
            scfint_delete(bfd->u.scfi);
            break;
         default:
            xassert(bfd != bfd);
      }
#ifdef GLP_DEBUG
      if (bfd->B != NULL)
         spm_delete_mat(bfd->B);
#endif
      tfree(bfd);
      return;
}

/* eof */
