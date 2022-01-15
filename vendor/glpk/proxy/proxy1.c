/* proxy1.c */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2013, 2018 Free Software Foundation, Inc.
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
#include "ios.h"
#include "proxy.h"

void ios_proxy_heur(glp_tree *T)
{     glp_prob *prob;
      int j, status;
      double *xstar, zstar;
      /* this heuristic is applied only once on the root level */
      if (!(T->curr->level == 0 && T->curr->solved == 1))
         goto done;
      prob = glp_create_prob();
      glp_copy_prob(prob, T->mip, 0);
      xstar = xcalloc(1+prob->n, sizeof(double));
      for (j = 1; j <= prob->n; j++)
         xstar[j] = 0.0;
      if (T->mip->mip_stat != GLP_FEAS)
         status = proxy(prob, &zstar, xstar, NULL, 0.0,
            T->parm->ps_tm_lim, 1);
      else
      {  double *xinit = xcalloc(1+prob->n, sizeof(double));
         for (j = 1; j <= prob->n; j++)
            xinit[j] = T->mip->col[j]->mipx;
         status = proxy(prob, &zstar, xstar, xinit, 0.0,
            T->parm->ps_tm_lim, 1);
         xfree(xinit);
      }
      if (status == 0)
#if 0 /* 17/III-2016 */
         glp_ios_heur_sol(T, xstar);
#else
      {  /* sometimes the proxy heuristic reports a wrong solution, so
          * make sure that the solution is really integer feasible */
         int i, feas1, feas2, ae_ind, re_ind;
         double ae_max, re_max;
         glp_copy_prob(prob, T->mip, 0);
         for (j = 1; j <= prob->n; j++)
            prob->col[j]->mipx = xstar[j];
         for (i = 1; i <= prob->m; i++)
         {  GLPROW *row;
            GLPAIJ *aij;
            row = prob->row[i];
            row->mipx = 0.0;
            for (aij = row->ptr; aij != NULL; aij = aij->r_next)
               row->mipx += aij->val * aij->col->mipx;
         }
         glp_check_kkt(prob, GLP_MIP, GLP_KKT_PE, &ae_max, &ae_ind,
            &re_max, &re_ind);
         feas1 = (re_max <= 1e-6);
         glp_check_kkt(prob, GLP_MIP, GLP_KKT_PB, &ae_max, &ae_ind,
            &re_max, &re_ind);
         feas2 = (re_max <= 1e-6);
         if (feas1 && feas2)
            glp_ios_heur_sol(T, xstar);
         else
            xprintf("WARNING: PROXY HEURISTIC REPORTED WRONG SOLUTION; "
               "SOLUTION REJECTED\n");
      }
#endif
      xfree(xstar);
      glp_delete_prob(prob);
done: return;
}

/* eof */
