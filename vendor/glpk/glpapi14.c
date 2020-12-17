/* glpapi14.c (processing models in GNU MathProg language) */

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

#define GLP_TRAN_DEFINED
typedef struct MPL glp_tran;

#include "glpmpl.h"
#include "glpapi.h"

glp_tran *glp_mpl_alloc_wksp(void)
{     /* allocate the MathProg translator workspace */
      glp_tran *tran;
      tran = mpl_initialize();
      return tran;
}

#if 1 /* 08/XII-2009 */
void _glp_mpl_init_rand(glp_tran *tran, int seed)
{     if (tran->phase != 0)
         xerror("glp_mpl_init_rand: invalid call sequence\n");
      rng_init_rand(tran->rand, seed);
      return;
}
#endif

int glp_mpl_read_model(glp_tran *tran, const char *fname, int skip)
{     /* read and translate model section */
      int ret;
      if (tran->phase != 0)
         xerror("glp_mpl_read_model: invalid call sequence\n");
      ret = mpl_read_model(tran, (char *)fname, skip);
      if (ret == 1 || ret == 2)
         ret = 0;
      else if (ret == 4)
         ret = 1;
      else
         xassert(ret != ret);
      return ret;
}

int glp_mpl_read_data(glp_tran *tran, const char *fname)
{     /* read and translate data section */
      int ret;
      if (!(tran->phase == 1 || tran->phase == 2))
         xerror("glp_mpl_read_data: invalid call sequence\n");
      ret = mpl_read_data(tran, (char *)fname);
      if (ret == 2)
         ret = 0;
      else if (ret == 4)
         ret = 1;
      else
         xassert(ret != ret);
      return ret;
}

int glp_mpl_generate(glp_tran *tran, const char *fname)
{     /* generate the model */
      int ret;
      if (!(tran->phase == 1 || tran->phase == 2))
         xerror("glp_mpl_generate: invalid call sequence\n");
      ret = mpl_generate(tran, (char *)fname);
      if (ret == 3)
         ret = 0;
      else if (ret == 4)
         ret = 1;
      return ret;
}

void glp_mpl_build_prob(glp_tran *tran, glp_prob *prob)
{     /* build LP/MIP problem instance from the model */
      int m, n, i, j, t, kind, type, len, *ind;
      double lb, ub, *val;
      if (tran->phase != 3)
         xerror("glp_mpl_build_prob: invalid call sequence\n");
      /* erase the problem object */
      glp_erase_prob(prob);
      /* set problem name */
      glp_set_prob_name(prob, mpl_get_prob_name(tran));
      /* build rows (constraints) */
      m = mpl_get_num_rows(tran);
      if (m > 0)
         glp_add_rows(prob, m);
      for (i = 1; i <= m; i++)
      {  /* set row name */
         glp_set_row_name(prob, i, mpl_get_row_name(tran, i));
         /* set row bounds */
         type = mpl_get_row_bnds(tran, i, &lb, &ub);
         switch (type)
         {  case MPL_FR: type = GLP_FR; break;
            case MPL_LO: type = GLP_LO; break;
            case MPL_UP: type = GLP_UP; break;
            case MPL_DB: type = GLP_DB; break;
            case MPL_FX: type = GLP_FX; break;
            default: xassert(type != type);
         }
         if (type == GLP_DB && fabs(lb - ub) < 1e-9 * (1.0 + fabs(lb)))
         {  type = GLP_FX;
            if (fabs(lb) <= fabs(ub)) ub = lb; else lb = ub;
         }
         glp_set_row_bnds(prob, i, type, lb, ub);
         /* warn about non-zero constant term */
         if (mpl_get_row_c0(tran, i) != 0.0)
            xprintf("glp_mpl_build_prob: row %s; constant term %.12g ig"
               "nored\n",
               mpl_get_row_name(tran, i), mpl_get_row_c0(tran, i));
      }
      /* build columns (variables) */
      n = mpl_get_num_cols(tran);
      if (n > 0)
         glp_add_cols(prob, n);
      for (j = 1; j <= n; j++)
      {  /* set column name */
         glp_set_col_name(prob, j, mpl_get_col_name(tran, j));
         /* set column kind */
         kind = mpl_get_col_kind(tran, j);
         switch (kind)
         {  case MPL_NUM:
               break;
            case MPL_INT:
            case MPL_BIN:
               glp_set_col_kind(prob, j, GLP_IV);
               break;
            default:
               xassert(kind != kind);
         }
         /* set column bounds */
         type = mpl_get_col_bnds(tran, j, &lb, &ub);
         switch (type)
         {  case MPL_FR: type = GLP_FR; break;
            case MPL_LO: type = GLP_LO; break;
            case MPL_UP: type = GLP_UP; break;
            case MPL_DB: type = GLP_DB; break;
            case MPL_FX: type = GLP_FX; break;
            default: xassert(type != type);
         }
         if (kind == MPL_BIN)
         {  if (type == GLP_FR || type == GLP_UP || lb < 0.0) lb = 0.0;
            if (type == GLP_FR || type == GLP_LO || ub > 1.0) ub = 1.0;
            type = GLP_DB;
         }
         if (type == GLP_DB && fabs(lb - ub) < 1e-9 * (1.0 + fabs(lb)))
         {  type = GLP_FX;
            if (fabs(lb) <= fabs(ub)) ub = lb; else lb = ub;
         }
         glp_set_col_bnds(prob, j, type, lb, ub);
      }
      /* load the constraint matrix */
      ind = xcalloc(1+n, sizeof(int));
      val = xcalloc(1+n, sizeof(double));
      for (i = 1; i <= m; i++)
      {  len = mpl_get_mat_row(tran, i, ind, val);
         glp_set_mat_row(prob, i, len, ind, val);
      }
      /* build objective function (the first objective is used) */
      for (i = 1; i <= m; i++)
      {  kind = mpl_get_row_kind(tran, i);
         if (kind == MPL_MIN || kind == MPL_MAX)
         {  /* set objective name */
            glp_set_obj_name(prob, mpl_get_row_name(tran, i));
            /* set optimization direction */
            glp_set_obj_dir(prob, kind == MPL_MIN ? GLP_MIN : GLP_MAX);
            /* set constant term */
            glp_set_obj_coef(prob, 0, mpl_get_row_c0(tran, i));
            /* set objective coefficients */
            len = mpl_get_mat_row(tran, i, ind, val);
            for (t = 1; t <= len; t++)
               glp_set_obj_coef(prob, ind[t], val[t]);
            break;
         }
      }
      /* free working arrays */
      xfree(ind);
      xfree(val);
      return;
}

int glp_mpl_postsolve(glp_tran *tran, glp_prob *prob, int sol)
{     /* postsolve the model */
      int i, j, m, n, stat, ret;
      double prim, dual;
      if (!(tran->phase == 3 && !tran->flag_p))
         xerror("glp_mpl_postsolve: invalid call sequence\n");
      if (!(sol == GLP_SOL || sol == GLP_IPT || sol == GLP_MIP))
         xerror("glp_mpl_postsolve: sol = %d; invalid parameter\n",
            sol);
      m = mpl_get_num_rows(tran);
      n = mpl_get_num_cols(tran);
      if (!(m == glp_get_num_rows(prob) &&
            n == glp_get_num_cols(prob)))
         xerror("glp_mpl_postsolve: wrong problem object\n");
      if (!mpl_has_solve_stmt(tran))
      {  ret = 0;
         goto done;
      }
      for (i = 1; i <= m; i++)
      {  if (sol == GLP_SOL)
         {  stat = glp_get_row_stat(prob, i);
            prim = glp_get_row_prim(prob, i);
            dual = glp_get_row_dual(prob, i);
         }
         else if (sol == GLP_IPT)
         {  stat = 0;
            prim = glp_ipt_row_prim(prob, i);
            dual = glp_ipt_row_dual(prob, i);
         }
         else if (sol == GLP_MIP)
         {  stat = 0;
            prim = glp_mip_row_val(prob, i);
            dual = 0.0;
         }
         else
            xassert(sol != sol);
         if (fabs(prim) < 1e-9) prim = 0.0;
         if (fabs(dual) < 1e-9) dual = 0.0;
         mpl_put_row_soln(tran, i, stat, prim, dual);
      }
      for (j = 1; j <= n; j++)
      {  if (sol == GLP_SOL)
         {  stat = glp_get_col_stat(prob, j);
            prim = glp_get_col_prim(prob, j);
            dual = glp_get_col_dual(prob, j);
         }
         else if (sol == GLP_IPT)
         {  stat = 0;
            prim = glp_ipt_col_prim(prob, j);
            dual = glp_ipt_col_dual(prob, j);
         }
         else if (sol == GLP_MIP)
         {  stat = 0;
            prim = glp_mip_col_val(prob, j);
            dual = 0.0;
         }
         else
            xassert(sol != sol);
         if (fabs(prim) < 1e-9) prim = 0.0;
         if (fabs(dual) < 1e-9) dual = 0.0;
         mpl_put_col_soln(tran, j, stat, prim, dual);
      }
      ret = mpl_postsolve(tran);
      if (ret == 3)
         ret = 0;
      else if (ret == 4)
         ret = 1;
done: return ret;
}

void glp_mpl_free_wksp(glp_tran *tran)
{     /* free the MathProg translator workspace */
      mpl_terminate(tran);
      return;
}

/* eof */
