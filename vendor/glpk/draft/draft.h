/* draft.h */

#ifndef DRAFT_H
#define DRAFT_H

#if 1 /* 28/III-2016 */
#define GLP_UNDOC 1
#endif
#include "glpk.h"

#if 1 /* 28/XI-2009 */
int _glp_analyze_row(glp_prob *P, int len, const int ind[],
      const double val[], int type, double rhs, double eps, int *_piv,
      double *_x, double *_dx, double *_y, double *_dy, double *_dz);
/* simulate one iteration of dual simplex method */
#endif

#endif

/* eof */
