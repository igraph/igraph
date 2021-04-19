/* glpnpp01.c */

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
#pragma clang diagnostic ignored "-Wshorten-64-to-32"
#pragma clang diagnostic ignored "-Wsometimes-uninitialized"
#endif

#include "glpnpp.h"

NPP *npp_create_wksp(void)
{     /* create LP/MIP preprocessor workspace */
      NPP *npp;
      npp = xmalloc(sizeof(NPP));
      npp->orig_dir = 0;
      npp->orig_m = npp->orig_n = npp->orig_nnz = 0;
      npp->pool = dmp_create_pool();
      npp->name = npp->obj = NULL;
      npp->c0 = 0.0;
      npp->nrows = npp->ncols = 0;
      npp->r_head = npp->r_tail = NULL;
      npp->c_head = npp->c_tail = NULL;
      npp->stack = dmp_create_pool();
      npp->top = NULL;
#if 0 /* 16/XII-2009 */
      memset(&npp->count, 0, sizeof(npp->count));
#endif
      npp->m = npp->n = npp->nnz = 0;
      npp->row_ref = npp->col_ref = NULL;
      npp->sol = npp->scaling = 0;
      npp->p_stat = npp->d_stat = npp->t_stat = npp->i_stat = 0;
      npp->r_stat = NULL;
      /*npp->r_prim =*/ npp->r_pi = NULL;
      npp->c_stat = NULL;
      npp->c_value = /*npp->c_dual =*/ NULL;
      return npp;
}

void npp_insert_row(NPP *npp, NPPROW *row, int where)
{     /* insert row to the row list */
      if (where == 0)
      {  /* insert row to the beginning of the row list */
         row->prev = NULL;
         row->next = npp->r_head;
         if (row->next == NULL)
            npp->r_tail = row;
         else
            row->next->prev = row;
         npp->r_head = row;
      }
      else
      {  /* insert row to the end of the row list */
         row->prev = npp->r_tail;
         row->next = NULL;
         if (row->prev == NULL)
            npp->r_head = row;
         else
            row->prev->next = row;
         npp->r_tail = row;
      }
      return;
}

void npp_remove_row(NPP *npp, NPPROW *row)
{     /* remove row from the row list */
      if (row->prev == NULL)
         npp->r_head = row->next;
      else
         row->prev->next = row->next;
      if (row->next == NULL)
         npp->r_tail = row->prev;
      else
         row->next->prev = row->prev;
      return;
}

void npp_activate_row(NPP *npp, NPPROW *row)
{     /* make row active */
      if (!row->temp)
      {  row->temp = 1;
         /* move the row to the beginning of the row list */
         npp_remove_row(npp, row);
         npp_insert_row(npp, row, 0);
      }
      return;
}

void npp_deactivate_row(NPP *npp, NPPROW *row)
{     /* make row inactive */
      if (row->temp)
      {  row->temp = 0;
         /* move the row to the end of the row list */
         npp_remove_row(npp, row);
         npp_insert_row(npp, row, 1);
      }
      return;
}

void npp_insert_col(NPP *npp, NPPCOL *col, int where)
{     /* insert column to the column list */
      if (where == 0)
      {  /* insert column to the beginning of the column list */
         col->prev = NULL;
         col->next = npp->c_head;
         if (col->next == NULL)
            npp->c_tail = col;
         else
            col->next->prev = col;
         npp->c_head = col;
      }
      else
      {  /* insert column to the end of the column list */
         col->prev = npp->c_tail;
         col->next = NULL;
         if (col->prev == NULL)
            npp->c_head = col;
         else
            col->prev->next = col;
         npp->c_tail = col;
      }
      return;
}

void npp_remove_col(NPP *npp, NPPCOL *col)
{     /* remove column from the column list */
      if (col->prev == NULL)
         npp->c_head = col->next;
      else
         col->prev->next = col->next;
      if (col->next == NULL)
         npp->c_tail = col->prev;
      else
         col->next->prev = col->prev;
      return;
}

void npp_activate_col(NPP *npp, NPPCOL *col)
{     /* make column active */
      if (!col->temp)
      {  col->temp = 1;
         /* move the column to the beginning of the column list */
         npp_remove_col(npp, col);
         npp_insert_col(npp, col, 0);
      }
      return;
}

void npp_deactivate_col(NPP *npp, NPPCOL *col)
{     /* make column inactive */
      if (col->temp)
      {  col->temp = 0;
         /* move the column to the end of the column list */
         npp_remove_col(npp, col);
         npp_insert_col(npp, col, 1);
      }
      return;
}

NPPROW *npp_add_row(NPP *npp)
{     /* add new row to the current problem */
      NPPROW *row;
      row = dmp_get_atom(npp->pool, sizeof(NPPROW));
      row->i = ++(npp->nrows);
      row->name = NULL;
      row->lb = -DBL_MAX, row->ub = +DBL_MAX;
      row->ptr = NULL;
      row->temp = 0;
      npp_insert_row(npp, row, 1);
      return row;
}

NPPCOL *npp_add_col(NPP *npp)
{     /* add new column to the current problem */
      NPPCOL *col;
      col = dmp_get_atom(npp->pool, sizeof(NPPCOL));
      col->j = ++(npp->ncols);
      col->name = NULL;
#if 0
      col->kind = GLP_CV;
#else
      col->is_int = 0;
#endif
      col->lb = col->ub = col->coef = 0.0;
      col->ptr = NULL;
      col->temp = 0;
      npp_insert_col(npp, col, 1);
      return col;
}

NPPAIJ *npp_add_aij(NPP *npp, NPPROW *row, NPPCOL *col, double val)
{     /* add new element to the constraint matrix */
      NPPAIJ *aij;
      aij = dmp_get_atom(npp->pool, sizeof(NPPAIJ));
      aij->row = row;
      aij->col = col;
      aij->val = val;
      aij->r_prev = NULL;
      aij->r_next = row->ptr;
      aij->c_prev = NULL;
      aij->c_next = col->ptr;
      if (aij->r_next != NULL)
         aij->r_next->r_prev = aij;
      if (aij->c_next != NULL)
         aij->c_next->c_prev = aij;
      row->ptr = col->ptr = aij;
      return aij;
}

int npp_row_nnz(NPP *npp, NPPROW *row)
{     /* count number of non-zero coefficients in row */
      NPPAIJ *aij;
      int nnz;
      xassert(npp == npp);
      nnz = 0;
      for (aij = row->ptr; aij != NULL; aij = aij->r_next)
         nnz++;
      return nnz;
}

int npp_col_nnz(NPP *npp, NPPCOL *col)
{     /* count number of non-zero coefficients in column */
      NPPAIJ *aij;
      int nnz;
      xassert(npp == npp);
      nnz = 0;
      for (aij = col->ptr; aij != NULL; aij = aij->c_next)
         nnz++;
      return nnz;
}

void *npp_push_tse(NPP *npp, int (*func)(NPP *npp, void *info),
      int size)
{     /* push new entry to the transformation stack */
      NPPTSE *tse;
      tse = dmp_get_atom(npp->stack, sizeof(NPPTSE));
      tse->func = func;
      tse->info = dmp_get_atom(npp->stack, size);
      tse->link = npp->top;
      npp->top = tse;
      return tse->info;
}

#if 1 /* 23/XII-2009 */
void npp_erase_row(NPP *npp, NPPROW *row)
{     /* erase row content to make it empty */
      NPPAIJ *aij;
      while (row->ptr != NULL)
      {  aij = row->ptr;
         row->ptr = aij->r_next;
         if (aij->c_prev == NULL)
            aij->col->ptr = aij->c_next;
         else
            aij->c_prev->c_next = aij->c_next;
         if (aij->c_next == NULL)
            ;
         else
            aij->c_next->c_prev = aij->c_prev;
         dmp_free_atom(npp->pool, aij, sizeof(NPPAIJ));
      }
      return;
}
#endif

void npp_del_row(NPP *npp, NPPROW *row)
{     /* remove row from the current problem */
#if 0 /* 23/XII-2009 */
      NPPAIJ *aij;
#endif
      if (row->name != NULL)
         dmp_free_atom(npp->pool, row->name, strlen(row->name)+1);
#if 0 /* 23/XII-2009 */
      while (row->ptr != NULL)
      {  aij = row->ptr;
         row->ptr = aij->r_next;
         if (aij->c_prev == NULL)
            aij->col->ptr = aij->c_next;
         else
            aij->c_prev->c_next = aij->c_next;
         if (aij->c_next == NULL)
            ;
         else
            aij->c_next->c_prev = aij->c_prev;
         dmp_free_atom(npp->pool, aij, sizeof(NPPAIJ));
      }
#else
      npp_erase_row(npp, row);
#endif
      npp_remove_row(npp, row);
      dmp_free_atom(npp->pool, row, sizeof(NPPROW));
      return;
}

void npp_del_col(NPP *npp, NPPCOL *col)
{     /* remove column from the current problem */
      NPPAIJ *aij;
      if (col->name != NULL)
         dmp_free_atom(npp->pool, col->name, strlen(col->name)+1);
      while (col->ptr != NULL)
      {  aij = col->ptr;
         col->ptr = aij->c_next;
         if (aij->r_prev == NULL)
            aij->row->ptr = aij->r_next;
         else
            aij->r_prev->r_next = aij->r_next;
         if (aij->r_next == NULL)
            ;
         else
            aij->r_next->r_prev = aij->r_prev;
         dmp_free_atom(npp->pool, aij, sizeof(NPPAIJ));
      }
      npp_remove_col(npp, col);
      dmp_free_atom(npp->pool, col, sizeof(NPPCOL));
      return;
}

void npp_del_aij(NPP *npp, NPPAIJ *aij)
{     /* remove element from the constraint matrix */
      if (aij->r_prev == NULL)
         aij->row->ptr = aij->r_next;
      else
         aij->r_prev->r_next = aij->r_next;
      if (aij->r_next == NULL)
         ;
      else
         aij->r_next->r_prev = aij->r_prev;
      if (aij->c_prev == NULL)
         aij->col->ptr = aij->c_next;
      else
         aij->c_prev->c_next = aij->c_next;
      if (aij->c_next == NULL)
         ;
      else
         aij->c_next->c_prev = aij->c_prev;
      dmp_free_atom(npp->pool, aij, sizeof(NPPAIJ));
      return;
}

void npp_load_prob(NPP *npp, glp_prob *orig, int names, int sol,
      int scaling)
{     /* load original problem into the preprocessor workspace */
      int m = orig->m;
      int n = orig->n;
      NPPROW **link;
      int i, j;
      double dir;
      xassert(names == GLP_OFF || names == GLP_ON);
      xassert(sol == GLP_SOL || sol == GLP_IPT || sol == GLP_MIP);
      xassert(scaling == GLP_OFF || scaling == GLP_ON);
      if (sol == GLP_MIP) xassert(!scaling);
      npp->orig_dir = orig->dir;
      if (npp->orig_dir == GLP_MIN)
         dir = +1.0;
      else if (npp->orig_dir == GLP_MAX)
         dir = -1.0;
      else
         xassert(npp != npp);
      npp->orig_m = m;
      npp->orig_n = n;
      npp->orig_nnz = orig->nnz;
      if (names && orig->name != NULL)
      {  npp->name = dmp_get_atom(npp->pool, strlen(orig->name)+1);
         strcpy(npp->name, orig->name);
      }
      if (names && orig->obj != NULL)
      {  npp->obj = dmp_get_atom(npp->pool, strlen(orig->obj)+1);
         strcpy(npp->obj, orig->obj);
      }
      npp->c0 = dir * orig->c0;
      /* load rows */
      link = xcalloc(1+m, sizeof(NPPROW *));
      for (i = 1; i <= m; i++)
      {  GLPROW *rrr = orig->row[i];
         NPPROW *row;
         link[i] = row = npp_add_row(npp);
         xassert(row->i == i);
         if (names && rrr->name != NULL)
         {  row->name = dmp_get_atom(npp->pool, strlen(rrr->name)+1);
            strcpy(row->name, rrr->name);
         }
         if (!scaling)
         {  if (rrr->type == GLP_FR)
               row->lb = -DBL_MAX, row->ub = +DBL_MAX;
            else if (rrr->type == GLP_LO)
               row->lb = rrr->lb, row->ub = +DBL_MAX;
            else if (rrr->type == GLP_UP)
               row->lb = -DBL_MAX, row->ub = rrr->ub;
            else if (rrr->type == GLP_DB)
               row->lb = rrr->lb, row->ub = rrr->ub;
            else if (rrr->type == GLP_FX)
               row->lb = row->ub = rrr->lb;
            else
               xassert(rrr != rrr);
         }
         else
         {  double rii = rrr->rii;
            if (rrr->type == GLP_FR)
               row->lb = -DBL_MAX, row->ub = +DBL_MAX;
            else if (rrr->type == GLP_LO)
               row->lb = rrr->lb * rii, row->ub = +DBL_MAX;
            else if (rrr->type == GLP_UP)
               row->lb = -DBL_MAX, row->ub = rrr->ub * rii;
            else if (rrr->type == GLP_DB)
               row->lb = rrr->lb * rii, row->ub = rrr->ub * rii;
            else if (rrr->type == GLP_FX)
               row->lb = row->ub = rrr->lb * rii;
            else
               xassert(rrr != rrr);
         }
      }
      /* load columns and constraint coefficients */
      for (j = 1; j <= n; j++)
      {  GLPCOL *ccc = orig->col[j];
         GLPAIJ *aaa;
         NPPCOL *col;
         col = npp_add_col(npp);
         xassert(col->j == j);
         if (names && ccc->name != NULL)
         {  col->name = dmp_get_atom(npp->pool, strlen(ccc->name)+1);
            strcpy(col->name, ccc->name);
         }
         if (sol == GLP_MIP)
#if 0
            col->kind = ccc->kind;
#else
            col->is_int = (char)(ccc->kind == GLP_IV);
#endif
         if (!scaling)
         {  if (ccc->type == GLP_FR)
               col->lb = -DBL_MAX, col->ub = +DBL_MAX;
            else if (ccc->type == GLP_LO)
               col->lb = ccc->lb, col->ub = +DBL_MAX;
            else if (ccc->type == GLP_UP)
               col->lb = -DBL_MAX, col->ub = ccc->ub;
            else if (ccc->type == GLP_DB)
               col->lb = ccc->lb, col->ub = ccc->ub;
            else if (ccc->type == GLP_FX)
               col->lb = col->ub = ccc->lb;
            else
               xassert(ccc != ccc);
            col->coef = dir * ccc->coef;
            for (aaa = ccc->ptr; aaa != NULL; aaa = aaa->c_next)
               npp_add_aij(npp, link[aaa->row->i], col, aaa->val);
         }
         else
         {  double sjj = ccc->sjj;
            if (ccc->type == GLP_FR)
               col->lb = -DBL_MAX, col->ub = +DBL_MAX;
            else if (ccc->type == GLP_LO)
               col->lb = ccc->lb / sjj, col->ub = +DBL_MAX;
            else if (ccc->type == GLP_UP)
               col->lb = -DBL_MAX, col->ub = ccc->ub / sjj;
            else if (ccc->type == GLP_DB)
               col->lb = ccc->lb / sjj, col->ub = ccc->ub / sjj;
            else if (ccc->type == GLP_FX)
               col->lb = col->ub = ccc->lb / sjj;
            else
               xassert(ccc != ccc);
            col->coef = dir * ccc->coef * sjj;
            for (aaa = ccc->ptr; aaa != NULL; aaa = aaa->c_next)
               npp_add_aij(npp, link[aaa->row->i], col,
                  aaa->row->rii * aaa->val * sjj);
         }
      }
      xfree(link);
      /* keep solution indicator and scaling option */
      npp->sol = sol;
      npp->scaling = scaling;
      return;
}

void npp_build_prob(NPP *npp, glp_prob *prob)
{     /* build resultant (preprocessed) problem */
      NPPROW *row;
      NPPCOL *col;
      NPPAIJ *aij;
      int i, j, type, len, *ind;
      double dir, *val;
      glp_erase_prob(prob);
      glp_set_prob_name(prob, npp->name);
      glp_set_obj_name(prob, npp->obj);
      glp_set_obj_dir(prob, npp->orig_dir);
      if (npp->orig_dir == GLP_MIN)
         dir = +1.0;
      else if (npp->orig_dir == GLP_MAX)
         dir = -1.0;
      else
         xassert(npp != npp);
      glp_set_obj_coef(prob, 0, dir * npp->c0);
      /* build rows */
      for (row = npp->r_head; row != NULL; row = row->next)
      {  row->temp = i = glp_add_rows(prob, 1);
         glp_set_row_name(prob, i, row->name);
         if (row->lb == -DBL_MAX && row->ub == +DBL_MAX)
            type = GLP_FR;
         else if (row->ub == +DBL_MAX)
            type = GLP_LO;
         else if (row->lb == -DBL_MAX)
            type = GLP_UP;
         else if (row->lb != row->ub)
            type = GLP_DB;
         else
            type = GLP_FX;
         glp_set_row_bnds(prob, i, type, row->lb, row->ub);
      }
      /* build columns and the constraint matrix */
      ind = xcalloc(1+prob->m, sizeof(int));
      val = xcalloc(1+prob->m, sizeof(double));
      for (col = npp->c_head; col != NULL; col = col->next)
      {  j = glp_add_cols(prob, 1);
         glp_set_col_name(prob, j, col->name);
#if 0
         glp_set_col_kind(prob, j, col->kind);
#else
         glp_set_col_kind(prob, j, col->is_int ? GLP_IV : GLP_CV);
#endif
         if (col->lb == -DBL_MAX && col->ub == +DBL_MAX)
            type = GLP_FR;
         else if (col->ub == +DBL_MAX)
            type = GLP_LO;
         else if (col->lb == -DBL_MAX)
            type = GLP_UP;
         else if (col->lb != col->ub)
            type = GLP_DB;
         else
            type = GLP_FX;
         glp_set_col_bnds(prob, j, type, col->lb, col->ub);
         glp_set_obj_coef(prob, j, dir * col->coef);
         len = 0;
         for (aij = col->ptr; aij != NULL; aij = aij->c_next)
         {  len++;
            ind[len] = aij->row->temp;
            val[len] = aij->val;
         }
         glp_set_mat_col(prob, j, len, ind, val);
      }
      xfree(ind);
      xfree(val);
      /* resultant problem has been built */
      npp->m = prob->m;
      npp->n = prob->n;
      npp->nnz = prob->nnz;
      npp->row_ref = xcalloc(1+npp->m, sizeof(int));
      npp->col_ref = xcalloc(1+npp->n, sizeof(int));
      for (row = npp->r_head, i = 0; row != NULL; row = row->next)
         npp->row_ref[++i] = row->i;
      for (col = npp->c_head, j = 0; col != NULL; col = col->next)
         npp->col_ref[++j] = col->j;
      /* transformed problem segment is no longer needed */
      dmp_delete_pool(npp->pool), npp->pool = NULL;
      npp->name = npp->obj = NULL;
      npp->c0 = 0.0;
      npp->r_head = npp->r_tail = NULL;
      npp->c_head = npp->c_tail = NULL;
      return;
}

void npp_postprocess(NPP *npp, glp_prob *prob)
{     /* postprocess solution from the resultant problem */
      GLPROW *row;
      GLPCOL *col;
      NPPTSE *tse;
      int i, j, k;
      double dir;
      xassert(npp->orig_dir == prob->dir);
      if (npp->orig_dir == GLP_MIN)
         dir = +1.0;
      else if (npp->orig_dir == GLP_MAX)
         dir = -1.0;
      else
         xassert(npp != npp);
      xassert(npp->m == prob->m);
      xassert(npp->n == prob->n);
      xassert(npp->nnz == prob->nnz);
      /* copy solution status */
      if (npp->sol == GLP_SOL)
      {  npp->p_stat = prob->pbs_stat;
         npp->d_stat = prob->dbs_stat;
      }
      else if (npp->sol == GLP_IPT)
         npp->t_stat = prob->ipt_stat;
      else if (npp->sol == GLP_MIP)
         npp->i_stat = prob->mip_stat;
      else
         xassert(npp != npp);
      /* allocate solution arrays */
      if (npp->sol == GLP_SOL)
      {  if (npp->r_stat == NULL)
            npp->r_stat = xcalloc(1+npp->nrows, sizeof(char));
         for (i = 1; i <= npp->nrows; i++)
            npp->r_stat[i] = 0;
         if (npp->c_stat == NULL)
            npp->c_stat = xcalloc(1+npp->ncols, sizeof(char));
         for (j = 1; j <= npp->ncols; j++)
            npp->c_stat[j] = 0;
      }
#if 0
      if (npp->r_prim == NULL)
         npp->r_prim = xcalloc(1+npp->nrows, sizeof(double));
      for (i = 1; i <= npp->nrows; i++)
         npp->r_prim[i] = DBL_MAX;
#endif
      if (npp->c_value == NULL)
         npp->c_value = xcalloc(1+npp->ncols, sizeof(double));
      for (j = 1; j <= npp->ncols; j++)
         npp->c_value[j] = DBL_MAX;
      if (npp->sol != GLP_MIP)
      {  if (npp->r_pi == NULL)
            npp->r_pi = xcalloc(1+npp->nrows, sizeof(double));
         for (i = 1; i <= npp->nrows; i++)
            npp->r_pi[i] = DBL_MAX;
#if 0
         if (npp->c_dual == NULL)
            npp->c_dual = xcalloc(1+npp->ncols, sizeof(double));
         for (j = 1; j <= npp->ncols; j++)
            npp->c_dual[j] = DBL_MAX;
#endif
      }
      /* copy solution components from the resultant problem */
      if (npp->sol == GLP_SOL)
      {  for (i = 1; i <= npp->m; i++)
         {  row = prob->row[i];
            k = npp->row_ref[i];
            npp->r_stat[k] = (char)row->stat;
            /*npp->r_prim[k] = row->prim;*/
            npp->r_pi[k] = dir * row->dual;
         }
         for (j = 1; j <= npp->n; j++)
         {  col = prob->col[j];
            k = npp->col_ref[j];
            npp->c_stat[k] = (char)col->stat;
            npp->c_value[k] = col->prim;
            /*npp->c_dual[k] = dir * col->dual;*/
         }
      }
      else if (npp->sol == GLP_IPT)
      {  for (i = 1; i <= npp->m; i++)
         {  row = prob->row[i];
            k = npp->row_ref[i];
            /*npp->r_prim[k] = row->pval;*/
            npp->r_pi[k] = dir * row->dval;
         }
         for (j = 1; j <= npp->n; j++)
         {  col = prob->col[j];
            k = npp->col_ref[j];
            npp->c_value[k] = col->pval;
            /*npp->c_dual[k] = dir * col->dval;*/
         }
      }
      else if (npp->sol == GLP_MIP)
      {
#if 0
         for (i = 1; i <= npp->m; i++)
         {  row = prob->row[i];
            k = npp->row_ref[i];
            /*npp->r_prim[k] = row->mipx;*/
         }
#endif
         for (j = 1; j <= npp->n; j++)
         {  col = prob->col[j];
            k = npp->col_ref[j];
            npp->c_value[k] = col->mipx;
         }
      }
      else
         xassert(npp != npp);
      /* perform postprocessing to construct solution to the original
         problem */
      for (tse = npp->top; tse != NULL; tse = tse->link)
      {  xassert(tse->func != NULL);
         xassert(tse->func(npp, tse->info) == 0);
      }
      return;
}

void npp_unload_sol(NPP *npp, glp_prob *orig)
{     /* store solution to the original problem */
      GLPROW *row;
      GLPCOL *col;
      int i, j;
      double dir;
      xassert(npp->orig_dir == orig->dir);
      if (npp->orig_dir == GLP_MIN)
         dir = +1.0;
      else if (npp->orig_dir == GLP_MAX)
         dir = -1.0;
      else
         xassert(npp != npp);
      xassert(npp->orig_m == orig->m);
      xassert(npp->orig_n == orig->n);
      xassert(npp->orig_nnz == orig->nnz);
      if (npp->sol == GLP_SOL)
      {  /* store basic solution */
         orig->valid = 0;
         orig->pbs_stat = npp->p_stat;
         orig->dbs_stat = npp->d_stat;
         orig->obj_val = orig->c0;
         orig->some = 0;
         for (i = 1; i <= orig->m; i++)
         {  row = orig->row[i];
            row->stat = npp->r_stat[i];
            if (!npp->scaling)
            {  /*row->prim = npp->r_prim[i];*/
               row->dual = dir * npp->r_pi[i];
            }
            else
            {  /*row->prim = npp->r_prim[i] / row->rii;*/
               row->dual = dir * npp->r_pi[i] * row->rii;
            }
            if (row->stat == GLP_BS)
               row->dual = 0.0;
            else if (row->stat == GLP_NL)
            {  xassert(row->type == GLP_LO || row->type == GLP_DB);
               row->prim = row->lb;
            }
            else if (row->stat == GLP_NU)
            {  xassert(row->type == GLP_UP || row->type == GLP_DB);
               row->prim = row->ub;
            }
            else if (row->stat == GLP_NF)
            {  xassert(row->type == GLP_FR);
               row->prim = 0.0;
            }
            else if (row->stat == GLP_NS)
            {  xassert(row->type == GLP_FX);
               row->prim = row->lb;
            }
            else
               xassert(row != row);
         }
         for (j = 1; j <= orig->n; j++)
         {  col = orig->col[j];
            col->stat = npp->c_stat[j];
            if (!npp->scaling)
            {  col->prim = npp->c_value[j];
               /*col->dual = dir * npp->c_dual[j];*/
            }
            else
            {  col->prim = npp->c_value[j] * col->sjj;
               /*col->dual = dir * npp->c_dual[j] / col->sjj;*/
            }
            if (col->stat == GLP_BS)
               col->dual = 0.0;
#if 1
            else if (col->stat == GLP_NL)
            {  xassert(col->type == GLP_LO || col->type == GLP_DB);
               col->prim = col->lb;
            }
            else if (col->stat == GLP_NU)
            {  xassert(col->type == GLP_UP || col->type == GLP_DB);
               col->prim = col->ub;
            }
            else if (col->stat == GLP_NF)
            {  xassert(col->type == GLP_FR);
               col->prim = 0.0;
            }
            else if (col->stat == GLP_NS)
            {  xassert(col->type == GLP_FX);
               col->prim = col->lb;
            }
            else
               xassert(col != col);
#endif
            orig->obj_val += col->coef * col->prim;
         }
#if 1
         /* compute primal values of inactive rows */
         for (i = 1; i <= orig->m; i++)
         {  row = orig->row[i];
            if (row->stat == GLP_BS)
            {  GLPAIJ *aij;
               double temp;
               temp = 0.0;
               for (aij = row->ptr; aij != NULL; aij = aij->r_next)
                  temp += aij->val * aij->col->prim;
               row->prim = temp;
            }
         }
         /* compute reduced costs of active columns */
         for (j = 1; j <= orig->n; j++)
         {  col = orig->col[j];
            if (col->stat != GLP_BS)
            {  GLPAIJ *aij;
               double temp;
               temp = col->coef;
               for (aij = col->ptr; aij != NULL; aij = aij->c_next)
                  temp -= aij->val * aij->row->dual;
               col->dual = temp;
            }
         }
#endif
      }
      else if (npp->sol == GLP_IPT)
      {  /* store interior-point solution */
         orig->ipt_stat = npp->t_stat;
         orig->ipt_obj = orig->c0;
         for (i = 1; i <= orig->m; i++)
         {  row = orig->row[i];
            if (!npp->scaling)
            {  /*row->pval = npp->r_prim[i];*/
               row->dval = dir * npp->r_pi[i];
            }
            else
            {  /*row->pval = npp->r_prim[i] / row->rii;*/
               row->dval = dir * npp->r_pi[i] * row->rii;
            }
         }
         for (j = 1; j <= orig->n; j++)
         {  col = orig->col[j];
            if (!npp->scaling)
            {  col->pval = npp->c_value[j];
               /*col->dval = dir * npp->c_dual[j];*/
            }
            else
            {  col->pval = npp->c_value[j] * col->sjj;
               /*col->dval = dir * npp->c_dual[j] / col->sjj;*/
            }
            orig->ipt_obj += col->coef * col->pval;
         }
#if 1
         /* compute row primal values */
         for (i = 1; i <= orig->m; i++)
         {  row = orig->row[i];
            {  GLPAIJ *aij;
               double temp;
               temp = 0.0;
               for (aij = row->ptr; aij != NULL; aij = aij->r_next)
                  temp += aij->val * aij->col->pval;
               row->pval = temp;
            }
         }
         /* compute column dual values */
         for (j = 1; j <= orig->n; j++)
         {  col = orig->col[j];
            {  GLPAIJ *aij;
               double temp;
               temp = col->coef;
               for (aij = col->ptr; aij != NULL; aij = aij->c_next)
                  temp -= aij->val * aij->row->dval;
               col->dval = temp;
            }
         }
#endif
      }
      else if (npp->sol == GLP_MIP)
      {  /* store MIP solution */
         xassert(!npp->scaling);
         orig->mip_stat = npp->i_stat;
         orig->mip_obj = orig->c0;
#if 0
         for (i = 1; i <= orig->m; i++)
         {  row = orig->row[i];
            /*row->mipx = npp->r_prim[i];*/
         }
#endif
         for (j = 1; j <= orig->n; j++)
         {  col = orig->col[j];
            col->mipx = npp->c_value[j];
            if (col->kind == GLP_IV)
               xassert(col->mipx == floor(col->mipx));
            orig->mip_obj += col->coef * col->mipx;
         }
#if 1
         /* compute row primal values */
         for (i = 1; i <= orig->m; i++)
         {  row = orig->row[i];
            {  GLPAIJ *aij;
               double temp;
               temp = 0.0;
               for (aij = row->ptr; aij != NULL; aij = aij->r_next)
                  temp += aij->val * aij->col->mipx;
               row->mipx = temp;
            }
         }
#endif
      }
      else
         xassert(npp != npp);
      return;
}

void npp_delete_wksp(NPP *npp)
{     /* delete LP/MIP preprocessor workspace */
      if (npp->pool != NULL)
         dmp_delete_pool(npp->pool);
      if (npp->stack != NULL)
         dmp_delete_pool(npp->stack);
      if (npp->row_ref != NULL)
         xfree(npp->row_ref);
      if (npp->col_ref != NULL)
         xfree(npp->col_ref);
      if (npp->r_stat != NULL)
         xfree(npp->r_stat);
#if 0
      if (npp->r_prim != NULL)
         xfree(npp->r_prim);
#endif
      if (npp->r_pi != NULL)
         xfree(npp->r_pi);
      if (npp->c_stat != NULL)
         xfree(npp->c_stat);
      if (npp->c_value != NULL)
         xfree(npp->c_value);
#if 0
      if (npp->c_dual != NULL)
         xfree(npp->c_dual);
#endif
      xfree(npp);
      return;
}

/* eof */
