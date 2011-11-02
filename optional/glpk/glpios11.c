/* glpios11.c (process cuts stored in the local cut pool) */

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

#include "glpios.h"

/***********************************************************************
*  NAME
*
*  ios_process_cuts - process cuts stored in the local cut pool
*
*  SYNOPSIS
*
*  #include "glpios.h"
*  void ios_process_cuts(glp_tree *T);
*
*  DESCRIPTION
*
*  The routine ios_process_cuts analyzes each cut currently stored in
*  the local cut pool, which must be non-empty, and either adds the cut
*  to the current subproblem or just discards it. All cuts are assumed
*  to be locally valid. On exit the local cut pool remains unchanged.
*
*  REFERENCES
*
*  1. E.Balas, S.Ceria, G.Cornuejols, "Mixed 0-1 Programming by
*     Lift-and-Project in a Branch-and-Cut Framework", Management Sc.,
*     42 (1996) 1229-1246.
*
*  2. G.Andreello, A.Caprara, and M.Fischetti, "Embedding Cuts in
*     a Branch&Cut Framework: a Computational Study with {0,1/2}-Cuts",
*     Preliminary Draft, October 28, 2003, pp.6-8. */

struct info
{     /* estimated cut efficiency */
      IOSCUT *cut;
      /* pointer to cut in the cut pool */
      char flag;
      /* if this flag is set, the cut is included into the current
         subproblem */
      double eff;
      /* cut efficacy (normalized residual) */
      double deg;
      /* lower bound to objective degradation */
};

static int fcmp(const void *arg1, const void *arg2)
{     const struct info *info1 = arg1, *info2 = arg2;
      if (info1->deg == 0.0 && info2->deg == 0.0)
      {  if (info1->eff > info2->eff) return -1;
         if (info1->eff < info2->eff) return +1;
      }
      else
      {  if (info1->deg > info2->deg) return -1;
         if (info1->deg < info2->deg) return +1;
      }
      return 0;
}

static double parallel(IOSCUT *a, IOSCUT *b, double work[]);

void ios_process_cuts(glp_tree *T)
{     IOSPOOL *pool;
      IOSCUT *cut;
      IOSAIJ *aij;
      struct info *info;
      int k, kk, max_cuts, len, ret, *ind;
      double *val, *work;
      /* the current subproblem must exist */
      xassert(T->curr != NULL);
      /* the pool must exist and be non-empty */
      pool = T->local;
      xassert(pool != NULL);
      xassert(pool->size > 0);
      /* allocate working arrays */
      info = xcalloc(1+pool->size, sizeof(struct info));
      ind = xcalloc(1+T->n, sizeof(int));
      val = xcalloc(1+T->n, sizeof(double));
      work = xcalloc(1+T->n, sizeof(double));
      for (k = 1; k <= T->n; k++) work[k] = 0.0;
      /* build the list of cuts stored in the cut pool */
      for (k = 0, cut = pool->head; cut != NULL; cut = cut->next)
         k++, info[k].cut = cut, info[k].flag = 0;
      xassert(k == pool->size);
      /* estimate efficiency of all cuts in the cut pool */
      for (k = 1; k <= pool->size; k++)
      {  double temp, dy, dz;
         cut = info[k].cut;
         /* build the vector of cut coefficients and compute its
            Euclidean norm */
         len = 0; temp = 0.0;
         for (aij = cut->ptr; aij != NULL; aij = aij->next)
         {  xassert(1 <= aij->j && aij->j <= T->n);
            len++, ind[len] = aij->j, val[len] = aij->val;
            temp += aij->val * aij->val;
         }
         if (temp < DBL_EPSILON * DBL_EPSILON) temp = DBL_EPSILON;
         /* transform the cut to express it only through non-basic
            (auxiliary and structural) variables */
         len = glp_transform_row(T->mip, len, ind, val);
         /* determine change in the cut value and in the objective
            value for the adjacent basis by simulating one step of the
            dual simplex */
         ret = _glp_analyze_row(T->mip, len, ind, val, cut->type,
            cut->rhs, 1e-9, NULL, NULL, NULL, NULL, &dy, &dz);
         /* determine normalized residual and lower bound to objective
            degradation */
         if (ret == 0)
         {  info[k].eff = fabs(dy) / sqrt(temp);
            /* if some reduced costs violates (slightly) their zero
               bounds (i.e. have wrong signs) due to round-off errors,
               dz also may have wrong sign being close to zero */
            if (T->mip->dir == GLP_MIN)
            {  if (dz < 0.0) dz = 0.0;
               info[k].deg = + dz;
            }
            else /* GLP_MAX */
            {  if (dz > 0.0) dz = 0.0;
               info[k].deg = - dz;
            }
         }
         else if (ret == 1)
         {  /* the constraint is not violated at the current point */
            info[k].eff = info[k].deg = 0.0;
         }
         else if (ret == 2)
         {  /* no dual feasible adjacent basis exists */
            info[k].eff = 1.0;
            info[k].deg = DBL_MAX;
         }
         else
            xassert(ret != ret);
         /* if the degradation is too small, just ignore it */
         if (info[k].deg < 0.01) info[k].deg = 0.0;
      }
      /* sort the list of cuts by decreasing objective degradation and
         then by decreasing efficacy */
      qsort(&info[1], pool->size, sizeof(struct info), fcmp);
      /* only first (most efficient) max_cuts in the list are qualified
         as candidates to be added to the current subproblem */
      max_cuts = (T->curr->level == 0 ? 90 : 10);
      if (max_cuts > pool->size) max_cuts = pool->size;
      /* add cuts to the current subproblem */
#if 0
      xprintf("*** adding cuts ***\n");
#endif
      for (k = 1; k <= max_cuts; k++)
      {  int i, len;
         /* if this cut seems to be inefficient, skip it */
         if (info[k].deg < 0.01 && info[k].eff < 0.01) continue;
         /* if the angle between this cut and every other cut included
            in the current subproblem is small, skip this cut */
         for (kk = 1; kk < k; kk++)
         {  if (info[kk].flag)
            {  if (parallel(info[k].cut, info[kk].cut, work) > 0.90)
                  break;
            }
         }
         if (kk < k) continue;
         /* add this cut to the current subproblem */
#if 0
         xprintf("eff = %g; deg = %g\n", info[k].eff, info[k].deg);
#endif
         cut = info[k].cut, info[k].flag = 1;
         i = glp_add_rows(T->mip, 1);
         if (cut->name != NULL)
            glp_set_row_name(T->mip, i, cut->name);
         xassert(T->mip->row[i]->origin == GLP_RF_CUT);
         T->mip->row[i]->klass = cut->klass;
         len = 0;
         for (aij = cut->ptr; aij != NULL; aij = aij->next)
            len++, ind[len] = aij->j, val[len] = aij->val;
         glp_set_mat_row(T->mip, i, len, ind, val);
         xassert(cut->type == GLP_LO || cut->type == GLP_UP);
         glp_set_row_bnds(T->mip, i, cut->type, cut->rhs, cut->rhs);
      }
      /* free working arrays */
      xfree(info);
      xfree(ind);
      xfree(val);
      xfree(work);
      return;
}

#if 0
/***********************************************************************
*  Given a cut a * x >= b (<= b) the routine efficacy computes the cut
*  efficacy as follows:
*
*     eff = d * (a * x~ - b) / ||a||,
*
*  where d is -1 (in case of '>= b') or +1 (in case of '<= b'), x~ is
*  the vector of values of structural variables in optimal solution to
*  LP relaxation of the current subproblem, ||a|| is the Euclidean norm
*  of the vector of cut coefficients.
*
*  If the cut is violated at point x~, the efficacy eff is positive,
*  and its value is the Euclidean distance between x~ and the cut plane
*  a * x = b in the space of structural variables.
*
*  Following geometrical intuition, it is quite natural to consider
*  this distance as a first-order measure of the expected efficacy of
*  the cut: the larger the distance the better the cut [1]. */

static double efficacy(glp_tree *T, IOSCUT *cut)
{     glp_prob *mip = T->mip;
      IOSAIJ *aij;
      double s = 0.0, t = 0.0, temp;
      for (aij = cut->ptr; aij != NULL; aij = aij->next)
      {  xassert(1 <= aij->j && aij->j <= mip->n);
         s += aij->val * mip->col[aij->j]->prim;
         t += aij->val * aij->val;
      }
      temp = sqrt(t);
      if (temp < DBL_EPSILON) temp = DBL_EPSILON;
      if (cut->type == GLP_LO)
         temp = (s >= cut->rhs ? 0.0 : (cut->rhs - s) / temp);
      else if (cut->type == GLP_UP)
         temp = (s <= cut->rhs ? 0.0 : (s - cut->rhs) / temp);
      else
         xassert(cut != cut);
      return temp;
}
#endif

/***********************************************************************
*  Given two cuts a1 * x >= b1 (<= b1) and a2 * x >= b2 (<= b2) the
*  routine parallel computes the cosine of angle between the cut planes
*  a1 * x = b1 and a2 * x = b2 (which is the acute angle between two
*  normals to these planes) in the space of structural variables as
*  follows:
*
*     cos phi = (a1' * a2) / (||a1|| * ||a2||),
*
*  where (a1' * a2) is a dot product of vectors of cut coefficients,
*  ||a1|| and ||a2|| are Euclidean norms of vectors a1 and a2.
*
*  Note that requirement cos phi = 0 forces the cuts to be orthogonal,
*  i.e. with disjoint support, while requirement cos phi <= 0.999 means
*  only avoiding duplicate (parallel) cuts [1]. */

static double parallel(IOSCUT *a, IOSCUT *b, double work[])
{     IOSAIJ *aij;
      double s = 0.0, sa = 0.0, sb = 0.0, temp;
      for (aij = a->ptr; aij != NULL; aij = aij->next)
      {  work[aij->j] = aij->val;
         sa += aij->val * aij->val;
      }
      for (aij = b->ptr; aij != NULL; aij = aij->next)
      {  s += work[aij->j] * aij->val;
         sb += aij->val * aij->val;
      }
      for (aij = a->ptr; aij != NULL; aij = aij->next)
         work[aij->j] = 0.0;
      temp = sqrt(sa) * sqrt(sb);
      if (temp < DBL_EPSILON * DBL_EPSILON) temp = DBL_EPSILON;
      return s / temp;
}

/* eof */
