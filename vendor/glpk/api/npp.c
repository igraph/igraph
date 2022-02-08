/* npp.c (LP/MIP preprocessing) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2017 Free Software Foundation, Inc.
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

glp_prep *glp_npp_alloc_wksp(void)
{     /* allocate the preprocessor workspace */
      glp_prep *prep;
      prep = npp_create_wksp();
      return prep;
}

void glp_npp_load_prob(glp_prep *prep, glp_prob *P, int sol, int names)
{     /* load original problem instance */
      if (prep->sol != 0)
         xerror("glp_npp_load_prob: invalid call sequence (original ins"
            "tance already loaded)\n");
      if (!(sol == GLP_SOL || sol == GLP_IPT || sol == GLP_MIP))
         xerror("glp_npp_load_prob: sol = %d; invalid parameter\n",
            sol);
      if (!(names == GLP_ON || names == GLP_OFF))
         xerror("glp_npp_load_prob: names = %d; invalid parameter\n",
            names);
      npp_load_prob(prep, P, names, sol, GLP_OFF);
      return;
}

int glp_npp_preprocess1(glp_prep *prep, int hard)
{     /* perform basic LP/MIP preprocessing */
      if (prep->sol == 0)
         xerror("glp_npp_preprocess1: invalid call sequence (original i"
            "nstance not loaded yet)\n");
      if (prep->pool == NULL)
         xerror("glp_npp_preprocess1: invalid call sequence (preprocess"
            "ing already finished)\n");
      if (!(hard == GLP_ON || hard == GLP_OFF))
         xerror("glp_npp_preprocess1: hard = %d; invalid parameter\n",
            hard);
      return npp_process_prob(prep, hard);
}

void glp_npp_build_prob(glp_prep *prep, glp_prob *Q)
{     /* build resultant problem instance */
      if (prep->sol == 0)
         xerror("glp_npp_build_prob: invalid call sequence (original in"
            "stance not loaded yet)\n");
      if (prep->pool == NULL)
         xerror("glp_npp_build_prob: invalid call sequence (resultant i"
            "nstance already built)\n");
      npp_build_prob(prep, Q);
      return;
}

void glp_npp_postprocess(glp_prep *prep, glp_prob *Q)
{     /* postprocess solution to resultant problem */
      if (prep->pool != NULL)
         xerror("glp_npp_postprocess: invalid call sequence (resultant "
            "instance not built yet)\n");
      if (!(prep->m == Q->m && prep->n == Q->n && prep->nnz == Q->nnz))
         xerror("glp_npp_postprocess: resultant instance mismatch\n");
      switch (prep->sol)
      {  case GLP_SOL:
            if (glp_get_status(Q) != GLP_OPT)
               xerror("glp_npp_postprocess: unable to recover non-optim"
                  "al basic solution\n");
            break;
         case GLP_IPT:
            if (glp_ipt_status(Q) != GLP_OPT)
               xerror("glp_npp_postprocess: unable to recover non-optim"
                  "al interior-point solution\n");
            break;
         case GLP_MIP:
            if (!(glp_mip_status(Q) == GLP_OPT || glp_mip_status(Q) ==
               GLP_FEAS))
               xerror("glp_npp_postprocess: unable to recover integer n"
                  "on-feasible solution\n");
            break;
         default:
            xassert(prep != prep);
      }
      npp_postprocess(prep, Q);
      return;
}

void glp_npp_obtain_sol(glp_prep *prep, glp_prob *P)
{     /* obtain solution to original problem */
      if (prep->pool != NULL)
         xerror("glp_npp_obtain_sol: invalid call sequence (resultant i"
            "nstance not built yet)\n");
      switch (prep->sol)
      {  case GLP_SOL:
            if (prep->p_stat == 0 || prep->d_stat == 0)
               xerror("glp_npp_obtain_sol: invalid call sequence (basic"
                  " solution not provided yet)\n");
            break;
         case GLP_IPT:
            if (prep->t_stat == 0)
               xerror("glp_npp_obtain_sol: invalid call sequence (inter"
                  "ior-point solution not provided yet)\n");
            break;
         case GLP_MIP:
            if (prep->i_stat == 0)
               xerror("glp_npp_obtain_sol: invalid call sequence (MIP s"
                  "olution not provided yet)\n");
            break;
         default:
            xassert(prep != prep);
      }
      if (!(prep->orig_dir == P->dir && prep->orig_m == P->m &&
            prep->orig_n == P->n && prep->orig_nnz == P->nnz))
         xerror("glp_npp_obtain_sol: original instance mismatch\n");
      npp_unload_sol(prep, P);
      return;
}

void glp_npp_free_wksp(glp_prep *prep)
{     /* free the preprocessor workspace */
      npp_delete_wksp(prep);
      return;
}

/* eof */
