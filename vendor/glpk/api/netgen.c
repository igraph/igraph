/* netgen.c */

#include "env.h"
#include "glpk.h"

int glp_netgen(glp_graph *G_, int v_rhs_, int a_cap_, int a_cost_,
      const int parm[1+15])
{     static const char func[] = "glp_netgen";
      xassert(G_ == G_);
      xassert(v_rhs_ == v_rhs_);
      xassert(a_cap_ == a_cap_);
      xassert(a_cost_ == a_cost_);
      xassert(parm == parm);
      xerror("%s: sorry, this routine is temporarily disabled due to li"
         "censing problems\n", func);
      abort();
      return -1;
}

/* eof */
