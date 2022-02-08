/* intfeas1.c (solve integer feasibility problem) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2011-2016 Free Software Foundation, Inc.
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
#include "npp.h"

int glp_intfeas1(glp_prob *P, int use_bound, int obj_bound)
{     /* solve integer feasibility problem */
      NPP *npp = NULL;
      glp_prob *mip = NULL;
      int *obj_ind = NULL;
      double *obj_val = NULL;
      int obj_row = 0;
      int i, j, k, obj_len, temp, ret;
#if 0 /* 04/IV-2016 */
      /* check the problem object */
      if (P == NULL || P->magic != GLP_PROB_MAGIC)
         xerror("glp_intfeas1: P = %p; invalid problem object\n",
            P);
#endif
      if (P->tree != NULL)
         xerror("glp_intfeas1: operation not allowed\n");
      /* integer solution is currently undefined */
      P->mip_stat = GLP_UNDEF;
      P->mip_obj = 0.0;
      /* check columns (variables) */
      for (j = 1; j <= P->n; j++)
      {  GLPCOL *col = P->col[j];
#if 0 /* binarization is not yet implemented */
         if (!(col->kind == GLP_IV || col->type == GLP_FX))
         {  xprintf("glp_intfeas1: column %d: non-integer non-fixed var"
               "iable not allowed\n", j);
#else
         if (!((col->kind == GLP_IV && col->lb == 0.0 && col->ub == 1.0)
            || col->type == GLP_FX))
         {  xprintf("glp_intfeas1: column %d: non-binary non-fixed vari"
               "able not allowed\n", j);
#endif
            ret = GLP_EDATA;
            goto done;
         }
         temp = (int)col->lb;
         if ((double)temp != col->lb)
         {  if (col->type == GLP_FX)
               xprintf("glp_intfeas1: column %d: fixed value %g is non-"
                  "integer or out of range\n", j, col->lb);
            else
               xprintf("glp_intfeas1: column %d: lower bound %g is non-"
                  "integer or out of range\n", j, col->lb);
            ret = GLP_EDATA;
            goto done;
         }
         temp = (int)col->ub;
         if ((double)temp != col->ub)
         {  xprintf("glp_intfeas1: column %d: upper bound %g is non-int"
               "eger or out of range\n", j, col->ub);
            ret = GLP_EDATA;
            goto done;
         }
         if (col->type == GLP_DB && col->lb > col->ub)
         {  xprintf("glp_intfeas1: column %d: lower bound %g is greater"
               " than upper bound %g\n", j, col->lb, col->ub);
            ret = GLP_EBOUND;
            goto done;
         }
      }
      /* check rows (constraints) */
      for (i = 1; i <= P->m; i++)
      {  GLPROW *row = P->row[i];
         GLPAIJ *aij;
         for (aij = row->ptr; aij != NULL; aij = aij->r_next)
         {  temp = (int)aij->val;
            if ((double)temp != aij->val)
            {  xprintf("glp_intfeas1: row = %d, column %d: constraint c"
                  "oefficient %g is non-integer or out of range\n",
                  i, aij->col->j, aij->val);
               ret = GLP_EDATA;
               goto done;
            }
         }
         temp = (int)row->lb;
         if ((double)temp != row->lb)
         {  if (row->type == GLP_FX)
               xprintf("glp_intfeas1: row = %d: fixed value %g is non-i"
                  "nteger or out of range\n", i, row->lb);
            else
               xprintf("glp_intfeas1: row = %d: lower bound %g is non-i"
                  "nteger or out of range\n", i, row->lb);
            ret = GLP_EDATA;
            goto done;
         }
         temp = (int)row->ub;
         if ((double)temp != row->ub)
         {  xprintf("glp_intfeas1: row = %d: upper bound %g is non-inte"
               "ger or out of range\n", i, row->ub);
            ret = GLP_EDATA;
            goto done;
         }
         if (row->type == GLP_DB && row->lb > row->ub)
         {  xprintf("glp_intfeas1: row %d: lower bound %g is greater th"
               "an upper bound %g\n", i, row->lb, row->ub);
            ret = GLP_EBOUND;
            goto done;
         }
      }
      /* check the objective function */
#if 1 /* 08/I-2017 by cmatraki & mao */
      if (!use_bound)
      {  /* skip check if no obj. bound is specified */
         goto skip;
      }
#endif
      temp = (int)P->c0;
      if ((double)temp != P->c0)
      {  xprintf("glp_intfeas1: objective constant term %g is non-integ"
            "er or out of range\n", P->c0);
         ret = GLP_EDATA;
         goto done;
      }
      for (j = 1; j <= P->n; j++)
      {  temp = (int)P->col[j]->coef;
         if ((double)temp != P->col[j]->coef)
         {  xprintf("glp_intfeas1: column %d: objective coefficient is "
               "non-integer or out of range\n", j, P->col[j]->coef);
            ret = GLP_EDATA;
            goto done;
         }
      }
#if 1 /* 08/I-2017 by cmatraki & mao */
skip: ;
#endif
      /* save the objective function and set it to zero */
      obj_ind = xcalloc(1+P->n, sizeof(int));
      obj_val = xcalloc(1+P->n, sizeof(double));
      obj_len = 0;
      obj_ind[0] = 0;
      obj_val[0] = P->c0;
      P->c0 = 0.0;
      for (j = 1; j <= P->n; j++)
      {  if (P->col[j]->coef != 0.0)
         {  obj_len++;
            obj_ind[obj_len] = j;
            obj_val[obj_len] = P->col[j]->coef;
            P->col[j]->coef = 0.0;
         }
      }
      /* add inequality to bound the objective function, if required */
      if (!use_bound)
         xprintf("Will search for ANY feasible solution\n");
      else
      {  xprintf("Will search only for solution not worse than %d\n",
            obj_bound);
         obj_row = glp_add_rows(P, 1);
         glp_set_mat_row(P, obj_row, obj_len, obj_ind, obj_val);
         if (P->dir == GLP_MIN)
            glp_set_row_bnds(P, obj_row,
               GLP_UP, 0.0, (double)obj_bound - obj_val[0]);
         else if (P->dir == GLP_MAX)
            glp_set_row_bnds(P, obj_row,
               GLP_LO, (double)obj_bound - obj_val[0], 0.0);
         else
            xassert(P != P);
      }
      /* create preprocessor workspace */
      xprintf("Translating to CNF-SAT...\n");
      xprintf("Original problem has %d row%s, %d column%s, and %d non-z"
         "ero%s\n", P->m, P->m == 1 ? "" : "s", P->n, P->n == 1 ? "" :
         "s", P->nnz, P->nnz == 1 ? "" : "s");
      npp = npp_create_wksp();
      /* load the original problem into the preprocessor workspace */
      npp_load_prob(npp, P, GLP_OFF, GLP_MIP, GLP_OFF);
      /* perform translation to SAT-CNF problem instance */
      ret = npp_sat_encode_prob(npp);
      if (ret == 0)
         ;
      else if (ret == GLP_ENOPFS)
         xprintf("PROBLEM HAS NO INTEGER FEASIBLE SOLUTION\n");
      else if (ret == GLP_ERANGE)
         xprintf("glp_intfeas1: translation to SAT-CNF failed because o"
            "f integer overflow\n");
      else
         xassert(ret != ret);
      if (ret != 0)
         goto done;
      /* build SAT-CNF problem instance and try to solve it */
      mip = glp_create_prob();
      npp_build_prob(npp, mip);
      ret = glp_minisat1(mip);
      /* only integer feasible solution can be postprocessed */
      if (!(mip->mip_stat == GLP_OPT || mip->mip_stat == GLP_FEAS))
      {  P->mip_stat = mip->mip_stat;
         goto done;
      }
      /* postprocess the solution found */
      npp_postprocess(npp, mip);
      /* the transformed problem is no longer needed */
      glp_delete_prob(mip), mip = NULL;
      /* store solution to the original problem object */
      npp_unload_sol(npp, P);
      /* change the solution status to 'integer feasible' */
      P->mip_stat = GLP_FEAS;
      /* check integer feasibility */
      for (i = 1; i <= P->m; i++)
      {  GLPROW *row;
         GLPAIJ *aij;
         double sum;
         row = P->row[i];
         sum = 0.0;
         for (aij = row->ptr; aij != NULL; aij = aij->r_next)
            sum += aij->val * aij->col->mipx;
         xassert(sum == row->mipx);
         if (row->type == GLP_LO || row->type == GLP_DB ||
             row->type == GLP_FX)
            xassert(sum >= row->lb);
         if (row->type == GLP_UP || row->type == GLP_DB ||
             row->type == GLP_FX)
            xassert(sum <= row->ub);
      }
      /* compute value of the original objective function */
      P->mip_obj = obj_val[0];
      for (k = 1; k <= obj_len; k++)
         P->mip_obj += obj_val[k] * P->col[obj_ind[k]]->mipx;
      xprintf("Objective value = %17.9e\n", P->mip_obj);
done: /* delete the transformed problem, if it exists */
      if (mip != NULL)
         glp_delete_prob(mip);
      /* delete the preprocessor workspace, if it exists */
      if (npp != NULL)
         npp_delete_wksp(npp);
      /* remove inequality used to bound the objective function */
      if (obj_row > 0)
      {  int ind[1+1];
         ind[1] = obj_row;
         glp_del_rows(P, 1, ind);
      }
      /* restore the original objective function */
      if (obj_ind != NULL)
      {  P->c0 = obj_val[0];
         for (k = 1; k <= obj_len; k++)
            P->col[obj_ind[k]]->coef = obj_val[k];
         xfree(obj_ind);
         xfree(obj_val);
      }
      return ret;
}

/* eof */
