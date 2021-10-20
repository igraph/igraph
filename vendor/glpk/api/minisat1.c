/* minisat1.c (driver to MiniSat solver) */

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
#include "minisat.h"
#include "prob.h"

int glp_minisat1(glp_prob *P)
{     /* solve CNF-SAT problem with MiniSat solver */
      solver *s;
      GLPAIJ *aij;
      int i, j, len, ret, *ind;
      double sum;
#if 0 /* 04/IV-2016 */
      /* check problem object */
      if (P == NULL || P->magic != GLP_PROB_MAGIC)
         xerror("glp_minisat1: P = %p; invalid problem object\n",
            P);
#endif
      if (P->tree != NULL)
         xerror("glp_minisat1: operation not allowed\n");
      /* integer solution is currently undefined */
      P->mip_stat = GLP_UNDEF;
      P->mip_obj = 0.0;
      /* check that problem object encodes CNF-SAT instance */
      if (glp_check_cnfsat(P) != 0)
      {  xprintf("glp_minisat1: problem object does not encode CNF-SAT "
            "instance\n");
         ret = GLP_EDATA;
         goto done;
      }
#if 0 /* 08/I-2017 by cmatraki */
#if 1 /* 07/XI-2015 */
      if (sizeof(void *) != sizeof(int))
      {  xprintf("glp_minisat1: sorry, MiniSat solver is not supported "
            "on 64-bit platforms\n");
         ret = GLP_EFAIL;
         goto done;
      }
#endif
#else
      if (sizeof(void *) != sizeof(size_t))
      {  xprintf("glp_minisat1: sorry, MiniSat solver is not supported "
            "on this platform\n");
         ret = GLP_EFAIL;
         goto done;
      }
#endif
      /* solve CNF-SAT problem */
      xprintf("Solving CNF-SAT problem...\n");
      xprintf("Instance has %d variable%s, %d clause%s, and %d literal%"
         "s\n", P->n, P->n == 1 ? "" : "s", P->m, P->m == 1 ? "" : "s",
         P->nnz, P->nnz == 1 ? "" : "s");
      /* if CNF-SAT has no clauses, it is satisfiable */
      if (P->m == 0)
      {  P->mip_stat = GLP_OPT;
         for (j = 1; j <= P->n; j++)
            P->col[j]->mipx = 0.0;
         goto fini;
      }
      /* if CNF-SAT has an empty clause, it is unsatisfiable */
      for (i = 1; i <= P->m; i++)
      {  if (P->row[i]->ptr == NULL)
         {  P->mip_stat = GLP_NOFEAS;
            goto fini;
         }
      }
      /* prepare input data for the solver */
      s = solver_new();
      solver_setnvars(s, P->n);
      ind = xcalloc(1+P->n, sizeof(int));
      for (i = 1; i <= P->m; i++)
      {  len = 0;
         for (aij = P->row[i]->ptr; aij != NULL; aij = aij->r_next)
         {  ind[++len] = toLit(aij->col->j-1);
            if (aij->val < 0.0)
               ind[len] = lit_neg(ind[len]);
         }
         xassert(len > 0);
#if 0 /* 08/I-2017 by cmatraki */
         xassert(solver_addclause(s, &ind[1], &ind[1+len]));
#else
         if (!solver_addclause(s, &ind[1], &ind[1+len]))
         {  /* found trivial conflict */
            xfree(ind);
            solver_delete(s);
            P->mip_stat = GLP_NOFEAS;
            goto fini;
         }
#endif
      }
      xfree(ind);
      /* call the solver */
      s->verbosity = 1;
      if (solver_solve(s, 0, 0))
      {  /* instance is reported as satisfiable */
         P->mip_stat = GLP_OPT;
         /* copy solution to the problem object */
         xassert(s->model.size == P->n);
         for (j = 1; j <= P->n; j++)
         {  P->col[j]->mipx =
               s->model.ptr[j-1] == l_True ? 1.0 : 0.0;
         }
         /* compute row values */
         for (i = 1; i <= P->m; i++)
         {  sum = 0;
            for (aij = P->row[i]->ptr; aij != NULL; aij = aij->r_next)
               sum += aij->val * aij->col->mipx;
            P->row[i]->mipx = sum;
         }
         /* check integer feasibility */
         for (i = 1; i <= P->m; i++)
         {  if (P->row[i]->mipx < P->row[i]->lb)
            {  /* solution is wrong */
               P->mip_stat = GLP_UNDEF;
               break;
            }
         }
      }
      else
      {  /* instance is reported as unsatisfiable */
         P->mip_stat = GLP_NOFEAS;
      }
      solver_delete(s);
fini: /* report the instance status */
      if (P->mip_stat == GLP_OPT)
      {  xprintf("SATISFIABLE\n");
         ret = 0;
      }
      else if (P->mip_stat == GLP_NOFEAS)
      {  xprintf("UNSATISFIABLE\n");
         ret = 0;
      }
      else
      {  xprintf("glp_minisat1: solver failed\n");
         ret = GLP_EFAIL;
      }
done: return ret;
}

/* eof */
