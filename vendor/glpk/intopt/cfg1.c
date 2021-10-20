/* cfg1.c (conflict graph) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2012-2018 Free Software Foundation, Inc.
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

#include "cfg.h"
#include "env.h"
#include "prob.h"
#include "wclique.h"
#include "wclique1.h"

/***********************************************************************
*  cfg_build_graph - build conflict graph
*
*  This routine builds the conflict graph. It analyzes the specified
*  problem object to discover original and implied packing inequalities
*  and adds corresponding cliques to the conflict graph.
*
*  Packing inequality has the form:
*
*      sum z[j] <= 1,                                                (1)
*     j in J
*
*  where z[j] = x[j] or z[j] = 1 - x[j], x[j] is an original binary
*  variable. Every packing inequality (1) is equivalent to a set of
*  edge inequalities:
*
*     z[i] + z[j] <= 1   for all i, j in J, i != j,                  (2)
*
*  and since every edge inequality (2) defines an edge in the conflict
*  graph, corresponding packing inequality (1) defines a clique.
*
*  To discover packing inequalities the routine analyzes constraints
*  of the specified MIP. To simplify the analysis each constraint is
*  analyzed separately. The analysis is performed as follows.
*
*  Let some original constraint be the following:
*
*     L <= sum a[j] x[j] <= U.                                       (3)
*
*  To analyze it the routine analyzes two constraints of "not greater
*  than" type:
*
*     sum (-a[j]) x[j] <= -L,                                        (4)
*
*     sum (+a[j]) x[j] <= +U,                                        (5)
*
*  which are relaxations of the original constraint (3). (If, however,
*  L = -oo, or U = +oo, corresponding constraint being redundant is not
*  analyzed.)
*
*  Let a constraint of "not greater than" type be the following:
*
*      sum  a[j] x[j] + sum  a[j] x[j] <= b,                         (6)
*     j in J           j in J'
*
*  where J is a subset of binary variables, J' is a subset of other
*  (continues and non-binary integer) variables. The constraint (6) is
*  is relaxed as follows, to eliminate non-binary variables:
*
*      sum  a[j] x[j] <= b -  sum  a[j] x[j] <= b',                  (7)
*     j in J                 j in J'
*
*     b' = sup(b -  sum  a[j] x[j]) =
*                  j in J'
*
*        = b - inf(sum a[j] x[j]) =
*
*        = b - sum inf(a[j] x[j]) =                                  (8)
*
*        = b -  sum  a[j] inf(x[j]) -  sum  a[j] sup(x[j]) =
*              a[j]>0                 a[j]<0
*
*        = b -  sum  a[j] l[j] -  sum  a[j] u[j],
*              a[j]>0            a[j]<0
*
*  where l[j] and u[j] are, resp., lower and upper bounds of x[j].
*
*  Then the routine transforms the relaxed constraint containing only
*  binary variables:
*
*     sum a[j] x[j] <= b                                             (9)
*
*  to an equivalent 0-1 knapsack constraint as follows:
*
*     sum  a[j] x[j] + sum  a[j] x[j] <= b   ==>
*    a[j]>0           a[j]<0
*
*     sum  a[j] x[j] + sum  a[j] (1 - x[j]) <= b   ==>
*    a[j]>0           a[j]<0                                        (10)
*
*     sum  (+a[j]) x[j] + sum  (-a[j]) x[j] <= b + sum  (-a[j])   ==>
*    a[j]>0              a[j]<0                   a[j]<0
*
*     sum a'[j] z[j] <= b',
*
*  where a'[j] = |a[j]| > 0, and
*
*            ( x[j]      if a[j] > 0
*     z[j] = <
*            ( 1 - x[j]  if a[j] < 0
*
*  is a binary variable, which is either original binary variable x[j]
*  or its complement.
*
*  Finally, the routine analyzes the resultant 0-1 knapsack inequality:
*
*       sum a[j] z[j] <= b,                                         (11)
*     j in J
*
*  where all a[j] are positive, to discover clique inequalities (1),
*  which are valid for (11) and therefore valid for (3). (It is assumed
*  that the original MIP has been preprocessed, so it is not checked,
*  for example, that b > 0 or that a[j] <= b.)
*
*  In principle, to discover any edge inequalities valid for (11) it
*  is sufficient to check whether a[i] + a[j] > b for all i, j in J,
*  i < j. However, this way requires O(|J|^2) checks, so the routine
*  analyses (11) in the following way, which is much more efficient in
*  many practical cases.
*
*  1. Let a[p] and a[q] be two minimal coefficients:
*
*     a[p] = min a[j],                                              (12)
*
*     a[q] = min a[j], j != p,                                      (13)
*
*  such that
*
*     a[p] + a[q] > b.                                              (14)
*
*  This means that a[i] + a[j] > b for any i, j in J, i != j, so
*
*     z[i] + z[j] <= 1                                              (15)
*
*  are valid for (11) for any i, j in J, i != j. This case means that
*  J define a clique in the conflict graph.
*
*  2. Otherwise, let a[p] and [q] be two maximal coefficients:
*
*     a[p] = max a[j],                                              (16)
*
*     a[q] = max a[j], j != p,                                      (17)
*
*  such that
*
*     a[p] + a[q] <= b.                                             (18)
*
*  This means that a[i] + a[j] <= b for any i, j in J, i != j, so in
*  this case no valid edge inequalities for (11) exist.
*
*  3. Otherwise, let all a[j] be ordered by descending their values:
*
*     a[1] >= a[2] >= ... >= a[p-1] >= a[p] >= a[p+1] >= ...        (19)
*
*  where p is such that
*
*     a[p-1] + a[p] >  b,                                           (20)
*
*     a[p] + a[p+1] <= b.                                           (21)
*
*  (May note that due to the former two cases in this case we always
*  have 2 <= p <= |J|-1.)
*
*  Since a[p] and a[p-1] are two minimal coefficients in the set
*  J' = {1, ..., p}, J' define a clique in the conflict graph for the
*  same reason as in the first case. Similarly, since a[p] and a[p+1]
*  are two maximal coefficients in the set J" = {p, ..., |J|}, no edge
*  inequalities exist for all i, j in J" for the same reason as in the
*  second case. Thus, to discover other edge inequalities (15) valid
*  for (11), the routine checks if a[i] + a[j] > b for all i in J',
*  j in J", i != j. */

#define is_binary(j) \
      (P->col[j]->kind == GLP_IV && P->col[j]->type == GLP_DB && \
      P->col[j]->lb == 0.0 && P->col[j]->ub == 1.0)
/* check if x[j] is binary variable */

struct term { int ind; double val; };
/* term a[j] * z[j] used to sort a[j]'s */

static int CDECL fcmp(const void *e1, const void *e2)
{     /* auxiliary routine called from qsort */
      const struct term *t1 = e1, *t2 = e2;
      if (t1->val > t2->val)
         return -1;
      else if (t1->val < t2->val)
         return +1;
      else
         return 0;
}

static void analyze_ineq(glp_prob *P, CFG *G, int len, int ind[],
      double val[], double rhs, struct term t[])
{     /* analyze inequality constraint (6) */
      /* P is the original MIP
       * G is the conflict graph to be built
       * len is the number of terms in the constraint
       * ind[1], ..., ind[len] are indices of variables x[j]
       * val[1], ..., val[len] are constraint coefficients a[j]
       * rhs is the right-hand side b
       * t[1+len] is a working array */
      int j, k, kk, p, q, type, new_len;
      /* eliminate non-binary variables; see (7) and (8) */
      new_len = 0;
      for (k = 1; k <= len; k++)
      {  /* get index of variable x[j] */
         j = ind[k];
         if (is_binary(j))
         {  /* x[j] remains in relaxed constraint */
            new_len++;
            ind[new_len] = j;
            val[new_len] = val[k];
         }
         else if (val[k] > 0.0)
         {  /* eliminate non-binary x[j] in case a[j] > 0 */
            /* b := b - a[j] * l[j]; see (8) */
            type = P->col[j]->type;
            if (type == GLP_FR || type == GLP_UP)
            {  /* x[j] has no lower bound */
               goto done;
            }
            rhs -= val[k] * P->col[j]->lb;
         }
         else /* val[j] < 0.0 */
         {  /* eliminate non-binary x[j] in case a[j] < 0 */
            /* b := b - a[j] * u[j]; see (8) */
            type = P->col[j]->type;
            if (type == GLP_FR || type == GLP_LO)
            {  /* x[j] has no upper bound */
               goto done;
            }
            rhs -= val[k] * P->col[j]->ub;
         }
      }
      len = new_len;
      /* now we have the constraint (9) */
      if (len <= 1)
      {  /* at least two terms are needed */
         goto done;
      }
      /* make all constraint coefficients positive; see (10) */
      for (k = 1; k <= len; k++)
      {  if (val[k] < 0.0)
         {  /* a[j] < 0; substitute x[j] = 1 - x'[j], where x'[j] is
             * a complement binary variable */
            ind[k] = -ind[k];
            val[k] = -val[k];
            rhs += val[k];
         }
      }
      /* now we have 0-1 knapsack inequality (11) */
      /* increase the right-hand side a bit to avoid false checks due
       * to rounding errors */
      rhs += 0.001 * (1.0 + fabs(rhs));
      /*** first case ***/
      /* find two minimal coefficients a[p] and a[q] */
      p = 0;
      for (k = 1; k <= len; k++)
      {  if (p == 0 || val[p] > val[k])
            p = k;
      }
      q = 0;
      for (k = 1; k <= len; k++)
      {  if (k != p && (q == 0 || val[q] > val[k]))
            q = k;
      }
      xassert(p != 0 && q != 0 && p != q);
      /* check condition (14) */
      if (val[p] + val[q] > rhs)
      {  /* all z[j] define a clique in the conflict graph */
         cfg_add_clique(G, len, ind);
         goto done;
      }
      /*** second case ***/
      /* find two maximal coefficients a[p] and a[q] */
      p = 0;
      for (k = 1; k <= len; k++)
      {  if (p == 0 || val[p] < val[k])
            p = k;
      }
      q = 0;
      for (k = 1; k <= len; k++)
      {  if (k != p && (q == 0 || val[q] < val[k]))
            q = k;
      }
      xassert(p != 0 && q != 0 && p != q);
      /* check condition (18) */
      if (val[p] + val[q] <= rhs)
      {  /* no valid edge inequalities exist */
         goto done;
      }
      /*** third case ***/
      xassert(len >= 3);
      /* sort terms in descending order of coefficient values */
      for (k = 1; k <= len; k++)
      {  t[k].ind = ind[k];
         t[k].val = val[k];
      }
      qsort(&t[1], len, sizeof(struct term), fcmp);
      for (k = 1; k <= len; k++)
      {  ind[k] = t[k].ind;
         val[k] = t[k].val;
      }
      /* now a[1] >= a[2] >= ... >= a[len-1] >= a[len] */
      /* note that a[1] + a[2] > b and a[len-1] + a[len] <= b due two
       * the former two cases */
      xassert(val[1] + val[2] > rhs);
      xassert(val[len-1] + val[len] <= rhs);
      /* find p according to conditions (20) and (21) */
      for (p = 2; p < len; p++)
      {  if (val[p] + val[p+1] <= rhs)
            break;
      }
      xassert(p < len);
      /* z[1], ..., z[p] define a clique in the conflict graph */
      cfg_add_clique(G, p, ind);
      /* discover other edge inequalities */
      for (k = 1; k <= p; k++)
      {  for (kk = p; kk <= len; kk++)
         {  if (k != kk && val[k] + val[kk] > rhs)
            {  int iii[1+2];
               iii[1] = ind[k];
               iii[2] = ind[kk];
               cfg_add_clique(G, 2, iii);
            }
         }
      }
done: return;
}

CFG *cfg_build_graph(void *P_)
{     glp_prob *P = P_;
      int m = P->m;
      int n = P->n;
      CFG *G;
      int i, k, type, len, *ind;
      double *val;
      struct term *t;
      /* create the conflict graph (number of its vertices cannot be
       * greater than double number of binary variables) */
      G = cfg_create_graph(n, 2 * glp_get_num_bin(P));
      /* allocate working arrays */
      ind = talloc(1+n, int);
      val = talloc(1+n, double);
      t = talloc(1+n, struct term);
      /* analyze constraints to discover edge inequalities */
      for (i = 1; i <= m; i++)
      {  type = P->row[i]->type;
         if (type == GLP_LO || type == GLP_DB || type == GLP_FX)
         {  /* i-th row has lower bound */
            /* analyze inequality sum (-a[j]) * x[j] <= -lb */
            len = glp_get_mat_row(P, i, ind, val);
            for (k = 1; k <= len; k++)
               val[k] = -val[k];
            analyze_ineq(P, G, len, ind, val, -P->row[i]->lb, t);
         }
         if (type == GLP_UP || type == GLP_DB || type == GLP_FX)
         {  /* i-th row has upper bound */
            /* analyze inequality sum (+a[j]) * x[j] <= +ub */
            len = glp_get_mat_row(P, i, ind, val);
            analyze_ineq(P, G, len, ind, val, +P->row[i]->ub, t);
         }
      }
      /* free working arrays */
      tfree(ind);
      tfree(val);
      tfree(t);
      return G;
}

/***********************************************************************
*  cfg_find_clique - find maximum weight clique in conflict graph
*
*  This routine finds a maximum weight clique in the conflict graph
*  G = (V, E), where the weight of vertex v in V is the value of
*  corresponding binary variable z (which is either an original binary
*  variable or its complement) in the optimal solution to LP relaxation
*  provided in the problem object. The goal is to find a clique in G,
*  whose weight is greater than 1, in which case corresponding packing
*  inequality is violated at the optimal point.
*
*  On exit the routine stores vertex indices of the conflict graph
*  included in the clique found to locations ind[1], ..., ind[len], and
*  returns len, which is the clique size. The clique weight is stored
*  in location pointed to by the parameter sum. If no clique has been
*  found, the routine returns 0.
*
*  Since the conflict graph may have a big number of vertices and be
*  quite dense, the routine uses an induced subgraph G' = (V', E'),
*  which is constructed as follows:
*
*  1. If the weight of some vertex v in V is zero (close to zero), it
*     is not included in V'. Obviously, including in a clique
*     zero-weight vertices does not change its weight, so if in G there
*     exist a clique of a non-zero weight, in G' exists a clique of the
*     same weight. This point is extremely important, because dropping
*     out zero-weight vertices can be done without retrieving lists of
*     adjacent vertices whose size may be very large.
*
*  2. Cumulative weight of vertex v in V is the sum of the weight of v
*     and weights of all vertices in V adjacent to v. Obviously, if
*     a clique includes a vertex v, the clique weight cannot be greater
*     than the cumulative weight of v. Since we are interested only in
*     cliques whose weight is greater than 1, vertices of V, whose
*     cumulative weight is not greater than 1, are not included in V'.
*
*  May note that in many practical cases the size of the induced
*  subgraph G' is much less than the size of the original conflict
*  graph G due to many binary variables, whose optimal values are zero
*  or close to zero. For example, it may happen that |V| = 100,000 and
*  |E| = 1e9 while |V'| = 50 and |E'| = 1000. */

struct csa
{     /* common storage area */
      glp_prob *P;
      /* original MIP */
      CFG *G;
      /* original conflict graph G = (V, E), |V| = nv */
      int *ind; /* int ind[1+nv]; */
      /* working array */
      /*--------------------------------------------------------------*/
      /* induced subgraph G' = (V', E') of original conflict graph */
      int nn;
      /* number of vertices in V' */
      int *vtoi; /* int vtoi[1+nv]; */
      /* vtoi[v] = i, 1 <= v <= nv, means that vertex v in V is vertex
       * i in V'; vtoi[v] = 0 means that vertex v is not included in
       * the subgraph */
      int *itov; /* int itov[1+nv]; */
      /* itov[i] = v, 1 <= i <= nn, means that vertex i in V' is vertex
       * v in V */
      double *wgt; /* double wgt[1+nv]; */
      /* wgt[i], 1 <= i <= nn, is a weight of vertex i in V', which is
       * the value of corresponding binary variable in optimal solution
       * to LP relaxation */
};

static void build_subgraph(struct csa *csa)
{     /* build induced subgraph */
      glp_prob *P = csa->P;
      int n = P->n;
      CFG *G = csa->G;
      int *ind = csa->ind;
      int *pos = G->pos;
      int *neg = G->neg;
      int nv = G->nv;
      int *ref = G->ref;
      int *vtoi = csa->vtoi;
      int *itov = csa->itov;
      double *wgt = csa->wgt;
      int j, k, v, w, nn, len;
      double z, sum;
      /* initially induced subgraph is empty */
      nn = 0;
      /* walk thru vertices of original conflict graph */
      for (v = 1; v <= nv; v++)
      {  /* determine value of binary variable z[j] that corresponds to
          * vertex v */
         j = ref[v];
         xassert(1 <= j && j <= n);
         if (pos[j] == v)
         {  /* z[j] = x[j], where x[j] is original variable */
            z = P->col[j]->prim;
         }
         else if (neg[j] == v)
         {  /* z[j] = 1 - x[j], where x[j] is original variable */
            z = 1.0 - P->col[j]->prim;
         }
         else
            xassert(v != v);
         /* if z[j] is close to zero, do not include v in the induced
          * subgraph */
         if (z < 0.001)
         {  vtoi[v] = 0;
            continue;
         }
         /* calculate cumulative weight of vertex v */
         sum = z;
         /* walk thru all vertices adjacent to v */
         len = cfg_get_adjacent(G, v, ind);
         for (k = 1; k <= len; k++)
         {  /* there is an edge (v,w) in the conflict graph */
            w = ind[k];
            xassert(w != v);
            /* add value of z[j] that corresponds to vertex w */
            j = ref[w];
            xassert(1 <= j && j <= n);
            if (pos[j] == w)
               sum += P->col[j]->prim;
            else if (neg[j] == w)
               sum += 1.0 - P->col[j]->prim;
            else
               xassert(w != w);
         }
         /* cumulative weight of vertex v is an upper bound of weight
          * of any clique containing v; so if it not greater than 1, do
          * not include v in the induced subgraph */
         if (sum < 1.010)
         {  vtoi[v] = 0;
            continue;
         }
         /* include vertex v in the induced subgraph */
         nn++;
         vtoi[v] = nn;
         itov[nn] = v;
         wgt[nn] = z;
      }
      /* induced subgraph has been built */
      csa->nn = nn;
      return;
}

static int sub_adjacent(struct csa *csa, int i, int adj[])
{     /* retrieve vertices of induced subgraph adjacent to specified
       * vertex */
      CFG *G = csa->G;
      int nv = G->nv;
      int *ind = csa->ind;
      int nn = csa->nn;
      int *vtoi = csa->vtoi;
      int *itov = csa->itov;
      int j, k, v, w, len, len1;
      /* determine original vertex v corresponding to vertex i */
      xassert(1 <= i && i <= nn);
      v = itov[i];
      /* retrieve vertices adjacent to vertex v in original graph */
      len1 = cfg_get_adjacent(G, v, ind);
      /* keep only adjacent vertices which are in induced subgraph and
       * change their numbers appropriately */
      len = 0;
      for (k = 1; k <= len1; k++)
      {  /* there exists edge (v, w) in original graph */
         w = ind[k];
         xassert(1 <= w && w <= nv && w != v);
         j = vtoi[w];
         if (j != 0)
         {  /* vertex w is vertex j in induced subgraph */
            xassert(1 <= j && j <= nn && j != i);
            adj[++len] = j;
         }
      }
      return len;
}

static int find_clique(struct csa *csa, int c_ind[])
{     /* find maximum weight clique in induced subgraph with exact
       * Ostergard's algorithm */
      int nn = csa->nn;
      double *wgt = csa->wgt;
      int i, j, k, p, q, t, ne, nb, len, *iwt, *ind;
      unsigned char *a;
      xassert(nn >= 2);
      /* allocate working array */
      ind = talloc(1+nn, int);
      /* calculate the number of elements in lower triangle (without
       * diagonal) of adjacency matrix of induced subgraph */
      ne = (nn * (nn - 1)) / 2;
      /* calculate the number of bytes needed to store lower triangle
       * of adjacency matrix */
      nb = (ne + (CHAR_BIT - 1)) / CHAR_BIT;
      /* allocate lower triangle of adjacency matrix */
      a = talloc(nb, unsigned char);
      /* fill lower triangle of adjacency matrix */
      memset(a, 0, nb);
      for (p = 1; p <= nn; p++)
      {  /* retrieve vertices adjacent to vertex p */
         len = sub_adjacent(csa, p, ind);
         for (k = 1; k <= len; k++)
         {  /* there exists edge (p, q) in induced subgraph */
            q = ind[k];
            xassert(1 <= q && q <= nn && q != p);
            /* determine row and column indices of this edge in lower
             * triangle of adjacency matrix */
            if (p > q)
               i = p, j = q;
            else /* p < q */
               i = q, j = p;
            /* set bit a[i,j] to 1, i > j */
            t = ((i - 1) * (i - 2)) / 2 + (j - 1);
            a[t / CHAR_BIT] |=
               (unsigned char)(1 << ((CHAR_BIT - 1) - t % CHAR_BIT));
         }
      }
      /* scale vertex weights by 1000 and convert them to integers as
       * required by Ostergard's algorithm */
      iwt = ind;
      for (i = 1; i <= nn; i++)
      {  /* it is assumed that 0 <= wgt[i] <= 1 */
         t = (int)(1000.0 * wgt[i] + 0.5);
         if (t < 0)
            t = 0;
         else if (t > 1000)
            t = 1000;
         iwt[i] = t;
      }
      /* find maximum weight clique */
      len = wclique(nn, iwt, a, c_ind);
      /* free working arrays */
      tfree(ind);
      tfree(a);
      /* return clique size to calling routine */
      return len;
}

static int func(void *info, int i, int ind[])
{     /* auxiliary routine used by routine find_clique1 */
      struct csa *csa = info;
      xassert(1 <= i && i <= csa->nn);
      return sub_adjacent(csa, i, ind);
}

static int find_clique1(struct csa *csa, int c_ind[])
{     /* find maximum weight clique in induced subgraph with greedy
       * heuristic */
      int nn = csa->nn;
      double *wgt = csa->wgt;
      int len;
      xassert(nn >= 2);
      len = wclique1(nn, wgt, func, csa, c_ind);
      /* return clique size to calling routine */
      return len;
}

int cfg_find_clique(void *P, CFG *G, int ind[], double *sum_)
{     int nv = G->nv;
      struct csa csa;
      int i, k, len;
      double sum;
      /* initialize common storage area */
      csa.P = P;
      csa.G = G;
      csa.ind = talloc(1+nv, int);
      csa.nn = -1;
      csa.vtoi = talloc(1+nv, int);
      csa.itov = talloc(1+nv, int);
      csa.wgt = talloc(1+nv, double);
      /* build induced subgraph */
      build_subgraph(&csa);
#ifdef GLP_DEBUG
      xprintf("nn = %d\n", csa.nn);
#endif
      /* if subgraph has less than two vertices, do nothing */
      if (csa.nn < 2)
      {  len = 0;
         sum = 0.0;
         goto skip;
      }
      /* find maximum weight clique in induced subgraph */
#if 1 /* FIXME */
      if (csa.nn <= 50)
#endif
      {  /* induced subgraph is small; use exact algorithm */
         len = find_clique(&csa, ind);
      }
      else
      {  /* induced subgraph is large; use greedy heuristic */
         len = find_clique1(&csa, ind);
      }
      /* do not report clique, if it has less than two vertices */
      if (len < 2)
      {  len = 0;
         sum = 0.0;
         goto skip;
      }
      /* convert indices of clique vertices from induced subgraph to
       * original conflict graph and compute clique weight */
      sum = 0.0;
      for (k = 1; k <= len; k++)
      {  i = ind[k];
         xassert(1 <= i && i <= csa.nn);
         sum += csa.wgt[i];
         ind[k] = csa.itov[i];
      }
skip: /* free working arrays */
      tfree(csa.ind);
      tfree(csa.vtoi);
      tfree(csa.itov);
      tfree(csa.wgt);
      /* return to calling routine */
      *sum_ = sum;
      return len;
}

/* eof */
