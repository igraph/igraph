/* rmfgen.c */

#include "env.h"
#include "glpk.h"

int glp_rmfgen(glp_graph *G_, int *s_, int *t_, int a_cap_,
      const int parm[1+5])
{     static const char func[] = "glp_rmfgen";
      xassert(G_ == G_);
      xassert(s_ == s_);
      xassert(t_ == t_);
      xassert(a_cap_ == a_cap_);
      xassert(parm == parm);
      xerror("%s: sorry, this routine is temporarily disabled due to li"
         "censing problems\n", func);
      abort();
      return -1;
}

/* eof */
