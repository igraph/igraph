/* wrcnf.c (write CNF-SAT problem data in DIMACS format) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2010-2016 Free Software Foundation, Inc.
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
#include "prob.h"

#define xfprintf        glp_format

int glp_write_cnfsat(glp_prob *P, const char *fname)
{     /* write CNF-SAT problem data in DIMACS format */
      glp_file *fp = NULL;
      GLPAIJ *aij;
      int i, j, len, count = 0, ret;
      char s[50];
#if 0 /* 04/IV-2016 */
      if (P == NULL || P->magic != GLP_PROB_MAGIC)
         xerror("glp_write_cnfsat: P = %p; invalid problem object\n",
            P);
#endif
      if (glp_check_cnfsat(P) != 0)
      {  xprintf("glp_write_cnfsat: problem object does not encode CNF-"
            "SAT instance\n");
         ret = 1;
         goto done;
      }
      xprintf("Writing CNF-SAT problem data to '%s'...\n", fname);
      fp = glp_open(fname, "w");
      if (fp == NULL)
      {  xprintf("Unable to create '%s' - %s\n", fname, get_err_msg());
         ret = 1;
         goto done;
      }
      xfprintf(fp, "c %s\n",
         P->name == NULL ? "unknown" : P->name), count++;
      xfprintf(fp, "p cnf %d %d\n", P->n, P->m), count++;
      for (i = 1; i <= P->m; i++)
      {  len = 0;
         for (aij = P->row[i]->ptr; aij != NULL; aij = aij->r_next)
         {  j = aij->col->j;
            if (aij->val < 0.0) j = -j;
            sprintf(s, "%d", j);
            if (len > 0 && len + 1 + strlen(s) > 72)
               xfprintf(fp, "\n"), count++, len = 0;
            xfprintf(fp, "%s%s", len == 0 ? "" : " ", s);
            if (len > 0) len++;
            len += strlen(s);
         }
         if (len > 0 && len + 1 + 1 > 72)
            xfprintf(fp, "\n"), count++, len = 0;
         xfprintf(fp, "%s0\n", len == 0 ? "" : " "), count++;
      }
      xfprintf(fp, "c eof\n"), count++;
#if 0 /* FIXME */
      xfflush(fp);
#endif
      if (glp_ioerr(fp))
      {  xprintf("Write error on '%s' - %s\n", fname, get_err_msg());
         ret = 1;
         goto done;
      }
      xprintf("%d lines were written\n", count);
      ret = 0;
done: if (fp != NULL) glp_close(fp);
      return ret;
}

/* eof */
