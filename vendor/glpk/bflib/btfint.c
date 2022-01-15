/* btfint.c (interface to BT-factorization) */

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
#include "btfint.h"

BTFINT *btfint_create(void)
{     /* create interface to BT-factorization */
      BTFINT *fi;
      fi = talloc(1, BTFINT);
      fi->n_max = 0;
      fi->valid = 0;
      fi->sva = NULL;
      fi->btf = NULL;
      fi->sgf = NULL;
      fi->sva_n_max = fi->sva_size = 0;
      fi->delta_n0 = fi->delta_n = 0;
      fi->sgf_piv_tol = 0.10;
      fi->sgf_piv_lim = 4;
      fi->sgf_suhl = 1;
      fi->sgf_eps_tol = DBL_EPSILON;
      return fi;
}

static void factorize_triv(BTFINT *fi, int k, int (*col)(void *info,
      int j, int ind[], double val[]), void *info)
{     /* compute LU-factorization of diagonal block A~[k,k] and store
       * corresponding columns of matrix A except elements of A~[k,k]
       * (trivial case when the block has unity size) */
      SVA *sva = fi->sva;
      int *sv_ind = sva->ind;
      double *sv_val = sva->val;
      BTF *btf = fi->btf;
      int *pp_inv = btf->pp_inv;
      int *qq_ind = btf->qq_ind;
      int *beg = btf->beg;
      int ac_ref = btf->ac_ref;
      int *ac_ptr = &sva->ptr[ac_ref-1];
      int *ac_len = &sva->len[ac_ref-1];
      SGF *sgf = fi->sgf;
      int *ind = (int *)sgf->vr_max; /* working array */
      double *val = sgf->work; /* working array */
      int i, j, t, len, ptr, beg_k;
      /* diagonal block A~[k,k] has the only element in matrix A~,
       * which is a~[beg[k],beg[k]] = a[i,j] */
      beg_k = beg[k];
      i = pp_inv[beg_k];
      j = qq_ind[beg_k];
      /* get j-th column of A */
      len = col(info, j, ind, val);
      /* find element a[i,j] = a~[beg[k],beg[k]] in j-th column */
      for (t = 1; t <= len; t++)
      {  if (ind[t] == i)
            break;
      }
      xassert(t <= len);
      /* compute LU-factorization of diagonal block A~[k,k], where
       * F = (1), V = (a[i,j]), P = Q = (1) (see the module LUF) */
#if 1 /* FIXME */
      xassert(val[t] != 0.0);
#endif
      btf->vr_piv[beg_k] = val[t];
      btf->p1_ind[beg_k] = btf->p1_inv[beg_k] = 1;
      btf->q1_ind[beg_k] = btf->q1_inv[beg_k] = 1;
      /* remove element a[i,j] = a~[beg[k],beg[k]] from j-th column */
      memmove(&ind[t], &ind[t+1], (len-t) * sizeof(int));
      memmove(&val[t], &val[t+1], (len-t) * sizeof(double));
      len--;
      /* and store resulting j-th column of A into BTF */
      if (len > 0)
      {  /* reserve locations for j-th column of A */
         if (sva->r_ptr - sva->m_ptr < len)
         {  sva_more_space(sva, len);
            sv_ind = sva->ind;
            sv_val = sva->val;
         }
         sva_reserve_cap(sva, ac_ref+(j-1), len);
         /* store j-th column of A (except elements of A~[k,k]) */
         ptr = ac_ptr[j];
         memcpy(&sv_ind[ptr], &ind[1], len * sizeof(int));
         memcpy(&sv_val[ptr], &val[1], len * sizeof(double));
         ac_len[j] = len;
      }
      return;
}

static int factorize_block(BTFINT *fi, int k, int (*col)(void *info,
      int j, int ind[], double val[]), void *info)
{     /* compute LU-factorization of diagonal block A~[k,k] and store
       * corresponding columns of matrix A except elements of A~[k,k]
       * (general case) */
      SVA *sva = fi->sva;
      int *sv_ind = sva->ind;
      double *sv_val = sva->val;
      BTF *btf = fi->btf;
      int *pp_ind = btf->pp_ind;
      int *qq_ind = btf->qq_ind;
      int *beg = btf->beg;
      int ac_ref = btf->ac_ref;
      int *ac_ptr = &sva->ptr[ac_ref-1];
      int *ac_len = &sva->len[ac_ref-1];
      SGF *sgf = fi->sgf;
      int *ind = (int *)sgf->vr_max; /* working array */
      double *val = sgf->work; /* working array */
      LUF luf;
      int *vc_ptr, *vc_len, *vc_cap;
      int i, ii, j, jj, t, len, cnt, ptr, beg_k;
      /* construct fake LUF for LU-factorization of A~[k,k] */
      sgf->luf = &luf;
      luf.n = beg[k+1] - (beg_k = beg[k]);
      luf.sva = sva;
      luf.fr_ref = btf->fr_ref + (beg_k-1);
      luf.fc_ref = btf->fc_ref + (beg_k-1);
      luf.vr_ref = btf->vr_ref + (beg_k-1);
      luf.vr_piv = btf->vr_piv + (beg_k-1);
      luf.vc_ref = btf->vc_ref + (beg_k-1);
      luf.pp_ind = btf->p1_ind + (beg_k-1);
      luf.pp_inv = btf->p1_inv + (beg_k-1);
      luf.qq_ind = btf->q1_ind + (beg_k-1);
      luf.qq_inv = btf->q1_inv + (beg_k-1);
      /* process columns of k-th block of matrix A~ */
      vc_ptr = &sva->ptr[luf.vc_ref-1];
      vc_len = &sva->len[luf.vc_ref-1];
      vc_cap = &sva->cap[luf.vc_ref-1];
      for (jj = 1; jj <= luf.n; jj++)
      {  /* jj-th column of A~ = j-th column of A */
         j = qq_ind[jj + (beg_k-1)];
         /* get j-th column of A */
         len = col(info, j, ind, val);
         /* move elements of diagonal block A~[k,k] to the beginning of
          * the column list */
         cnt = 0;
         for (t = 1; t <= len; t++)
         {  /* i = row index of element a[i,j] */
            i = ind[t];
            /* i-th row of A = ii-th row of A~ */
            ii = pp_ind[i];
            if (ii >= beg_k)
            {  /* a~[ii,jj] = a[i,j] is in diagonal block A~[k,k] */
               double temp;
               cnt++;
               ind[t] = ind[cnt];
               ind[cnt] = ii - (beg_k-1); /* local index */
               temp = val[t], val[t] = val[cnt], val[cnt] = temp;
            }
         }
         /* first cnt elements in the column list give jj-th column of
          * diagonal block A~[k,k], which is initial matrix V in LUF */
         /* enlarge capacity of jj-th column of V = A~[k,k] */
         if (vc_cap[jj] < cnt)
         {  if (sva->r_ptr - sva->m_ptr < cnt)
            {  sva_more_space(sva, cnt);
               sv_ind = sva->ind;
               sv_val = sva->val;
            }
            sva_enlarge_cap(sva, luf.vc_ref+(jj-1), cnt, 0);
         }
         /* store jj-th column of V = A~[k,k] */
         ptr = vc_ptr[jj];
         memcpy(&sv_ind[ptr], &ind[1], cnt * sizeof(int));
         memcpy(&sv_val[ptr], &val[1], cnt * sizeof(double));
         vc_len[jj] = cnt;
         /* other (len-cnt) elements in the column list are stored in
          * j-th column of the original matrix A */
         len -= cnt;
         if (len > 0)
         {  /* reserve locations for j-th column of A */
            if (sva->r_ptr - sva->m_ptr < len)
            {  sva_more_space(sva, len);
               sv_ind = sva->ind;
               sv_val = sva->val;
            }
            sva_reserve_cap(sva, ac_ref-1+j, len);
            /* store j-th column of A (except elements of A~[k,k]) */
            ptr = ac_ptr[j];
            memcpy(&sv_ind[ptr], &ind[cnt+1], len * sizeof(int));
            memcpy(&sv_val[ptr], &val[cnt+1], len * sizeof(double));
            ac_len[j] = len;
         }
      }
      /* compute LU-factorization of diagonal block A~[k,k]; may note
       * that A~[k,k] is irreducible (strongly connected), so singleton
       * phase will have no effect */
      k = sgf_factorize(sgf, 0 /* disable singleton phase */);
      /* now left (dynamic) part of SVA should be empty (wichtig!) */
      xassert(sva->m_ptr == 1);
      return k;
}

int btfint_factorize(BTFINT *fi, int n, int (*col)(void *info, int j,
      int ind[], double val[]), void *info)
{     /* compute BT-factorization of specified matrix A */
      SVA *sva;
      BTF *btf;
      SGF *sgf;
      int k, rank;
      xassert(n > 0);
      fi->valid = 0;
      /* create sparse vector area (SVA), if necessary */
      sva = fi->sva;
      if (sva == NULL)
      {  int sva_n_max = fi->sva_n_max;
         int sva_size = fi->sva_size;
         if (sva_n_max == 0)
            sva_n_max = 6 * n;
         if (sva_size == 0)
            sva_size = 10 * n;
         sva = fi->sva = sva_create_area(sva_n_max, sva_size);
      }
      /* allocate/reallocate underlying objects, if necessary */
      if (fi->n_max < n)
      {  int n_max = fi->n_max;
         if (n_max == 0)
            n_max = fi->n_max = n + fi->delta_n0;
         else
            n_max = fi->n_max = n + fi->delta_n;
         xassert(n_max >= n);
         /* allocate/reallocate block triangular factorization (BTF) */
         btf = fi->btf;
         if (btf == NULL)
         {  btf = fi->btf = talloc(1, BTF);
            memset(btf, 0, sizeof(BTF));
            btf->sva = sva;
         }
         else
         {  tfree(btf->pp_ind);
            tfree(btf->pp_inv);
            tfree(btf->qq_ind);
            tfree(btf->qq_inv);
            tfree(btf->beg);
            tfree(btf->vr_piv);
            tfree(btf->p1_ind);
            tfree(btf->p1_inv);
            tfree(btf->q1_ind);
            tfree(btf->q1_inv);
         }
         btf->pp_ind = talloc(1+n_max, int);
         btf->pp_inv = talloc(1+n_max, int);
         btf->qq_ind = talloc(1+n_max, int);
         btf->qq_inv = talloc(1+n_max, int);
         btf->beg = talloc(1+n_max+1, int);
         btf->vr_piv = talloc(1+n_max, double);
         btf->p1_ind = talloc(1+n_max, int);
         btf->p1_inv = talloc(1+n_max, int);
         btf->q1_ind = talloc(1+n_max, int);
         btf->q1_inv = talloc(1+n_max, int);
         /* allocate/reallocate factorizer workspace (SGF) */
         /* (note that for SGF we could use the size of largest block
          * rather than n_max) */
         sgf = fi->sgf;
         sgf = fi->sgf;
         if (sgf == NULL)
         {  sgf = fi->sgf = talloc(1, SGF);
            memset(sgf, 0, sizeof(SGF));
         }
         else
         {  tfree(sgf->rs_head);
            tfree(sgf->rs_prev);
            tfree(sgf->rs_next);
            tfree(sgf->cs_head);
            tfree(sgf->cs_prev);
            tfree(sgf->cs_next);
            tfree(sgf->vr_max);
            tfree(sgf->flag);
            tfree(sgf->work);
         }
         sgf->rs_head = talloc(1+n_max, int);
         sgf->rs_prev = talloc(1+n_max, int);
         sgf->rs_next = talloc(1+n_max, int);
         sgf->cs_head = talloc(1+n_max, int);
         sgf->cs_prev = talloc(1+n_max, int);
         sgf->cs_next = talloc(1+n_max, int);
         sgf->vr_max = talloc(1+n_max, double);
         sgf->flag = talloc(1+n_max, char);
         sgf->work = talloc(1+n_max, double);
      }
      btf = fi->btf;
      btf->n = n;
      sgf = fi->sgf;
#if 1 /* FIXME */
      /* initialize SVA */
      sva->n = 0;
      sva->m_ptr = 1;
      sva->r_ptr = sva->size + 1;
      sva->head = sva->tail = 0;
#endif
      /* store pattern of original matrix A in column-wise format */
      btf->ac_ref = sva_alloc_vecs(btf->sva, btf->n);
      btf_store_a_cols(btf, col, info, btf->pp_ind, btf->vr_piv);
#ifdef GLP_DEBUG
      sva_check_area(sva);
#endif
      /* analyze pattern of original matrix A and determine permutation
       * matrices P and Q such that A = P * A~* Q, where A~ is an upper
       * block triangular matrix */
      rank = btf_make_blocks(btf);
      if (rank != n)
      {  /* original matrix A is structurally singular */
         return 1;
      }
#ifdef GLP_DEBUG
      btf_check_blocks(btf);
#endif
#if 1 /* FIXME */
      /* initialize SVA */
      sva->n = 0;
      sva->m_ptr = 1;
      sva->r_ptr = sva->size + 1;
      sva->head = sva->tail = 0;
#endif
      /* allocate sparse vectors in SVA */
      btf->ar_ref = sva_alloc_vecs(btf->sva, btf->n);
      btf->ac_ref = sva_alloc_vecs(btf->sva, btf->n);
      btf->fr_ref = sva_alloc_vecs(btf->sva, btf->n);
      btf->fc_ref = sva_alloc_vecs(btf->sva, btf->n);
      btf->vr_ref = sva_alloc_vecs(btf->sva, btf->n);
      btf->vc_ref = sva_alloc_vecs(btf->sva, btf->n);
      /* setup factorizer control parameters */
      sgf->updat = 0; /* wichtig! */
      sgf->piv_tol = fi->sgf_piv_tol;
      sgf->piv_lim = fi->sgf_piv_lim;
      sgf->suhl = fi->sgf_suhl;
      sgf->eps_tol = fi->sgf_eps_tol;
      /* compute LU-factorizations of diagonal blocks A~[k,k] and also
       * store corresponding columns of matrix A except elements of all
       * blocks A~[k,k] */
      for (k = 1; k <= btf->num; k++)
      {  if (btf->beg[k+1] - btf->beg[k] == 1)
         {  /* trivial case (A~[k,k] has unity order) */
            factorize_triv(fi, k, col, info);
         }
         else
         {  /* general case */
            if (factorize_block(fi, k, col, info) != 0)
               return 2; /* factorization of A~[k,k] failed */
         }
      }
#ifdef GLP_DEBUG
      sva_check_area(sva);
#endif
      /* build row-wise representation of matrix A */
      btf_build_a_rows(fi->btf, fi->sgf->rs_head);
#ifdef GLP_DEBUG
      sva_check_area(sva);
#endif
      /* BT-factorization has been successfully computed */
      fi->valid = 1;
      return 0;
}

void btfint_delete(BTFINT *fi)
{     /* delete interface to BT-factorization */
      SVA *sva = fi->sva;
      BTF *btf = fi->btf;
      SGF *sgf = fi->sgf;
      if (sva != NULL)
         sva_delete_area(sva);
      if (btf != NULL)
      {  tfree(btf->pp_ind);
         tfree(btf->pp_inv);
         tfree(btf->qq_ind);
         tfree(btf->qq_inv);
         tfree(btf->beg);
         tfree(btf->vr_piv);
         tfree(btf->p1_ind);
         tfree(btf->p1_inv);
         tfree(btf->q1_ind);
         tfree(btf->q1_inv);
         tfree(btf);
      }
      if (sgf != NULL)
      {  tfree(sgf->rs_head);
         tfree(sgf->rs_prev);
         tfree(sgf->rs_next);
         tfree(sgf->cs_head);
         tfree(sgf->cs_prev);
         tfree(sgf->cs_next);
         tfree(sgf->vr_max);
         tfree(sgf->flag);
         tfree(sgf->work);
         tfree(sgf);
      }
      tfree(fi);
      return;
}

/* eof */
