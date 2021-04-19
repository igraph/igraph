/* glpnet04.c (grid-like network problem generator) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*
*  This code is a modified version of the program GRIDGEN, a grid-like
*  network problem generator developed by Yusin Lee and Jim Orlin.
*  The original code is publically available on the DIMACS ftp site at:
*  <ftp://dimacs.rutgers.edu/pub/netflow/generators/network/gridgen>.
*
*  All changes concern only the program interface, so this modified
*  version produces exactly the same instances as the original version.
*
*  Changes were made by Andrew Makhorin <mao@gnu.org>.
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

#ifdef __clang__
#pragma clang diagnostic ignored "-Wlogical-op-parentheses"
#endif

#include "glpapi.h"

/***********************************************************************
*  NAME
*
*  glp_gridgen - grid-like network problem generator
*
*  SYNOPSIS
*
*  int glp_gridgen(glp_graph *G, int v_rhs, int a_cap, int a_cost,
*     const int parm[1+14]);
*
*  DESCRIPTION
*
*  The routine glp_gridgen is a grid-like network problem generator
*  developed by Yusin Lee and Jim Orlin.
*
*  The parameter G specifies the graph object, to which the generated
*  problem data have to be stored. Note that on entry the graph object
*  is erased with the routine glp_erase_graph.
*
*  The parameter v_rhs specifies an offset of the field of type double
*  in the vertex data block, to which the routine stores the supply or
*  demand value. If v_rhs < 0, the value is not stored.
*
*  The parameter a_cap specifies an offset of the field of type double
*  in the arc data block, to which the routine stores the arc capacity.
*  If a_cap < 0, the capacity is not stored.
*
*  The parameter a_cost specifies an offset of the field of type double
*  in the arc data block, to which the routine stores the per-unit cost
*  if the arc flow. If a_cost < 0, the cost is not stored.
*
*  The array parm contains description of the network to be generated:
*
*  parm[0]  not used
*  parm[1]  two-ways arcs indicator:
*           1 - if links in both direction should be generated
*           0 - otherwise
*  parm[2]  random number seed (a positive integer)
*  parm[3]  number of nodes (the number of nodes generated might be
*           slightly different to make the network a grid)
*  parm[4]  grid width
*  parm[5]  number of sources
*  parm[6]  number of sinks
*  parm[7]  average degree
*  parm[8]  total flow
*  parm[9]  distribution of arc costs:
*           1 - uniform
*           2 - exponential
*  parm[10] lower bound for arc cost (uniform)
*           100 * lambda (exponential)
*  parm[11] upper bound for arc cost (uniform)
*           not used (exponential)
*  parm[12] distribution of arc capacities:
*           1 - uniform
*           2 - exponential
*  parm[13] lower bound for arc capacity (uniform)
*           100 * lambda (exponential)
*  parm[14] upper bound for arc capacity (uniform)
*           not used (exponential)
*
*  RETURNS
*
*  If the instance was successfully generated, the routine glp_gridgen
*  returns zero; otherwise, if specified parameters are inconsistent,
*  the routine returns a non-zero error code.
*
*  COMMENTS
*
*  This network generator generates a grid-like network plus a super
*  node. In additional to the arcs connecting the nodes in the grid,
*  there is an arc from each supply node to the super node and from the
*  super node to each demand node to guarantee feasiblity. These arcs
*  have very high costs and very big capacities.
*
*  The idea of this network generator is as follows: First, a grid of
*  n1 * n2 is generated. For example, 5 * 3. The nodes are numbered as
*  1 to 15, and the supernode is numbered as n1*n2+1. Then arcs between
*  adjacent nodes are generated. For these arcs, the user is allowed to
*  specify either to generate two-way arcs or one-way arcs. If two-way
*  arcs are to be generated, two arcs, one in each direction, will be
*  generated between each adjacent node pairs. Otherwise, only one arc
*  will be generated. If this is the case, the arcs will be generated
*  in alterntive directions as shown below.
*
*      1 ---> 2 ---> 3 ---> 4 ---> 5
*      |      ^      |      ^      |
*      |      |      |      |      |
*      V      |      V      |      V
*      6 <--- 7 <--- 8 <--- 9 <--- 10
*      |      ^      |      ^      |
*      |      |      |      |      |
*      V      |      V      |      V
*     11 --->12 --->13 --->14 ---> 15
*
*  Then the arcs between the super node and the source/sink nodes are
*  added as mentioned before. If the number of arcs still doesn't reach
*  the requirement, additional arcs will be added by uniformly picking
*  random node pairs. There is no checking to prevent multiple arcs
*  between any pair of nodes. However, there will be no self-arcs (arcs
*  that poins back to its tail node) in the network.
*
*  The source and sink nodes are selected uniformly in the network, and
*  the imbalances of each source/sink node are also assigned by uniform
*  distribution. */

struct stat_para
{     /* structure for statistical distributions */
      int distribution;
      /* the distribution: */
#define UNIFORM      1  /* uniform distribution */
#define EXPONENTIAL  2  /* exponential distribution */
      double parameter[5];
      /* the parameters of the distribution */
};

struct arcs
{     int from;
      /* the FROM node of that arc */
      int to;
      /* the TO node of that arc */
      int cost;
      /* original cost of that arc */
      int u;
      /* capacity of the arc */
};

struct imbalance
{     int node;
      /* Node ID */
      int supply;
      /* Supply of that node */
};

struct csa
{     /* common storage area */
      glp_graph *G;
      int v_rhs, a_cap, a_cost;
      int seed;
      /* random number seed */
      int seed_original;
      /* the original seed from input */
      int two_way;
      /* 0: generate arcs in both direction for the basic grid, except
         for the arcs to/from the super node.  1: o/w */
      int n_node;
      /* total number of nodes in the network, numbered 1 to n_node,
         including the super node, which is the last one */
      int n_arc;
      /* total number of arcs in the network, counting EVERY arc. */
      int n_grid_arc;
      /* number of arcs in the basic grid, including the arcs to/from
         the super node */
      int n_source, n_sink;
      /* number of source and sink nodes */
      int avg_degree;
      /* average degree, arcs to and from the super node are counted */
      int t_supply;
      /* total supply in the network */
      int n1, n2;
      /* the two edges of the network grid.  n1 >= n2 */
      struct imbalance *source_list, *sink_list;
      /* head of the array of source/sink nodes */
      struct stat_para arc_costs;
      /* the distribution of arc costs */
      struct stat_para capacities;
      /* distribution of the capacities of the arcs */
      struct arcs *arc_list;
      /* head of the arc list array.  Arcs in this array are in the
         order of grid_arcs, arcs to/from super node, and other arcs */
};

#define G (csa->G)
#define v_rhs (csa->v_rhs)
#define a_cap (csa->a_cap)
#define a_cost (csa->a_cost)
#define seed (csa->seed)
#define seed_original (csa->seed_original)
#define two_way (csa->two_way)
#define n_node (csa->n_node)
#define n_arc (csa->n_arc)
#define n_grid_arc (csa->n_grid_arc)
#define n_source (csa->n_source)
#define n_sink (csa->n_sink)
#define avg_degree (csa->avg_degree)
#define t_supply (csa->t_supply)
#define n1 (csa->n1)
#define n2 (csa->n2)
#define source_list (csa->source_list)
#define sink_list (csa->sink_list)
#define arc_costs (csa->arc_costs)
#define capacities (csa->capacities)
#define arc_list (csa->arc_list)

static void assign_capacities(struct csa *csa);
static void assign_costs(struct csa *csa);
static void assign_imbalance(struct csa *csa);
static int exponential(struct csa *csa, double lambda[1]);
static struct arcs *gen_additional_arcs(struct csa *csa, struct arcs
      *arc_ptr);
static struct arcs *gen_basic_grid(struct csa *csa, struct arcs
      *arc_ptr);
static void gen_more_arcs(struct csa *csa, struct arcs *arc_ptr);
static void generate(struct csa *csa);
static void output(struct csa *csa);
static double randy(struct csa *csa);
static void select_source_sinks(struct csa *csa);
static int uniform(struct csa *csa, double a[2]);

int glp_gridgen(glp_graph *G_, int _v_rhs, int _a_cap, int _a_cost,
      const int parm[1+14])
{     struct csa _csa, *csa = &_csa;
      int n, ret;
      G = G_;
      v_rhs = _v_rhs;
      a_cap = _a_cap;
      a_cost = _a_cost;
      if (G != NULL)
      {  if (v_rhs >= 0 && v_rhs > G->v_size - (int)sizeof(double))
            xerror("glp_gridgen: v_rhs = %d; invalid offset\n", v_rhs);
         if (a_cap >= 0 && a_cap > G->a_size - (int)sizeof(double))
            xerror("glp_gridgen: a_cap = %d; invalid offset\n", a_cap);
         if (a_cost >= 0 && a_cost > G->a_size - (int)sizeof(double))
            xerror("glp_gridgen: a_cost = %d; invalid offset\n", a_cost)
               ;
      }
      /* Check the parameters for consistency. */
      if (!(parm[1] == 0 || parm[1] == 1))
      {  ret = 1;
         goto done;
      }
      if (parm[2] < 1)
      {  ret = 2;
         goto done;
      }
      if (!(10 <= parm[3] && parm[3] <= 40000))
      {  ret = 3;
         goto done;
      }
      if (!(1 <= parm[4] && parm[4] <= 40000))
      {  ret = 4;
         goto done;
      }
      if (!(parm[5] >= 0 && parm[6] >= 0 && parm[5] + parm[6] <=
         parm[3]))
      {  ret = 5;
         goto done;
      }
      if (!(1 <= parm[7] && parm[7] <= parm[3]))
      {  ret = 6;
         goto done;
      }
      if (parm[8] < 0)
      {  ret = 7;
         goto done;
      }
      if (!(parm[9] == 1 || parm[9] == 2))
      {  ret = 8;
         goto done;
      }
      if (parm[9] == 1 && parm[10] > parm[11] ||
          parm[9] == 2 && parm[10] < 1)
      {  ret = 9;
         goto done;
      }
      if (!(parm[12] == 1 || parm[12] == 2))
      {  ret = 10;
         goto done;
      }
      if (parm[12] == 1 && !(0 <= parm[13] && parm[13] <= parm[14]) ||
          parm[12] == 2 && parm[13] < 1)
      {  ret = 11;
         goto done;
      }
      /* Initialize the graph object. */
      if (G != NULL)
      {  glp_erase_graph(G, G->v_size, G->a_size);
         glp_set_graph_name(G, "GRIDGEN");
      }
      /* Copy the generator parameters. */
      two_way = parm[1];
      seed_original = seed = parm[2];
      n_node = parm[3];
      n = parm[4];
      n_source = parm[5];
      n_sink = parm[6];
      avg_degree = parm[7];
      t_supply = parm[8];
      arc_costs.distribution = parm[9];
      if (parm[9] == 1)
      {  arc_costs.parameter[0] = parm[10];
         arc_costs.parameter[1] = parm[11];
      }
      else
      {  arc_costs.parameter[0] = (double)parm[10] / 100.0;
         arc_costs.parameter[1] = 0.0;
      }
      capacities.distribution = parm[12];
      if (parm[12] == 1)
      {  capacities.parameter[0] = parm[13];
         capacities.parameter[1] = parm[14];
      }
      else
      {  capacities.parameter[0] = (double)parm[13] / 100.0;
         capacities.parameter[1] = 0.0;
      }
      /* Calculate the edge lengths of the grid according to the
         input. */
      if (n * n >= n_node)
      {  n1 = n;
         n2 = (int)((double)n_node / (double)n + 0.5);
      }
      else
      {  n2 = n;
         n1 = (int)((double)n_node / (double)n + 0.5);
      }
      /* Recalculate the total number of nodes and plus 1 for the super
         node. */
      n_node = n1 * n2 + 1;
      n_arc = n_node * avg_degree;
      n_grid_arc = (two_way + 1) * ((n1 - 1) * n2 + (n2 - 1) * n1) +
         n_source + n_sink;
      if (n_grid_arc > n_arc) n_arc = n_grid_arc;
      arc_list = xcalloc(n_arc, sizeof(struct arcs));
      source_list = xcalloc(n_source, sizeof(struct imbalance));
      sink_list = xcalloc(n_sink, sizeof(struct imbalance));
      /* Generate a random network. */
      generate(csa);
      /* Output the network. */
      output(csa);
      /* Free all allocated memory. */
      xfree(arc_list);
      xfree(source_list);
      xfree(sink_list);
      /* The instance has been successfully generated. */
      ret = 0;
done: return ret;
}

#undef random

static void assign_capacities(struct csa *csa)
{     /* Assign a capacity to each arc. */
      struct arcs *arc_ptr = arc_list;
      int (*random)(struct csa *csa, double *);
      int i;
      /* Determine the random number generator to use. */
      switch (arc_costs.distribution)
      {  case UNIFORM:
            random = uniform;
            break;
         case EXPONENTIAL:
            random = exponential;
            break;
         default:
            xassert(csa != csa);
      }
      /* Assign capacities to grid arcs. */
      for (i = n_source + n_sink; i < n_grid_arc; i++, arc_ptr++)
         arc_ptr->u = random(csa, capacities.parameter);
      i = i - n_source - n_sink;
      /* Assign capacities to arcs to/from supernode. */
      for (; i < n_grid_arc; i++, arc_ptr++)
         arc_ptr->u = t_supply;
      /* Assign capacities to all other arcs. */
      for (; i < n_arc; i++, arc_ptr++)
         arc_ptr->u = random(csa, capacities.parameter);
      return;
}

static void assign_costs(struct csa *csa)
{     /* Assign a cost to each arc. */
      struct arcs *arc_ptr = arc_list;
      int (*random)(struct csa *csa, double *);
      int i;
      /* A high cost assigned to arcs to/from the supernode. */
      int high_cost;
      /* The maximum cost assigned to arcs in the base grid. */
      int max_cost = 0;
      /* Determine the random number generator to use. */
      switch (arc_costs.distribution)
      {  case UNIFORM:
            random = uniform;
            break;
         case EXPONENTIAL:
            random = exponential;
            break;
         default:
            xassert(csa != csa);
      }
      /* Assign costs to arcs in the base grid. */
      for (i = n_source + n_sink; i < n_grid_arc; i++, arc_ptr++)
      {  arc_ptr->cost = random(csa, arc_costs.parameter);
         if (max_cost < arc_ptr->cost) max_cost = arc_ptr->cost;
      }
      i = i - n_source - n_sink;
      /* Assign costs to arcs to/from the super node. */
      high_cost = max_cost * 2;
      for (; i < n_grid_arc; i++, arc_ptr++)
         arc_ptr->cost = high_cost;
      /* Assign costs to all other arcs. */
      for (; i < n_arc; i++, arc_ptr++)
         arc_ptr->cost = random(csa, arc_costs.parameter);
      return;
}

static void assign_imbalance(struct csa *csa)
{     /* Assign an imbalance to each node. */
      int total, i;
      double avg;
      struct imbalance *ptr;
      /* assign the supply nodes */
      avg = 2.0 * t_supply / n_source;
      do
      {  for (i = 1, total = t_supply, ptr = source_list + 1;
            i < n_source; i++, ptr++)
         {  ptr->supply = (int)(randy(csa) * avg + 0.5);
            total -= ptr->supply;
         }
         source_list->supply = total;
      }
      /* redo all if the assignment "overshooted" */
      while (total <= 0);
      /* assign the demand nodes */
      avg = -2.0 * t_supply / n_sink;
      do
      {  for (i = 1, total = t_supply, ptr = sink_list + 1;
            i < n_sink; i++, ptr++)
         {  ptr->supply = (int)(randy(csa) * avg - 0.5);
            total += ptr->supply;
         }
         sink_list->supply = - total;
      }
      while (total <= 0);
      return;
}

static int exponential(struct csa *csa, double lambda[1])
{     /* Returns an "exponentially distributed" integer with parameter
         lambda. */
      return ((int)(- lambda[0] * log((double)randy(csa)) + 0.5));
}

static struct arcs *gen_additional_arcs(struct csa *csa, struct arcs
      *arc_ptr)
{     /* Generate an arc from each source to the supernode and from
         supernode to each sink. */
      int i;
      for (i = 0; i < n_source; i++, arc_ptr++)
      {  arc_ptr->from = source_list[i].node;
         arc_ptr->to = n_node;
      }
      for (i = 0; i < n_sink; i++, arc_ptr++)
      {  arc_ptr->to = sink_list[i].node;
         arc_ptr->from = n_node;
      }
      return arc_ptr;
}

static struct arcs *gen_basic_grid(struct csa *csa, struct arcs
      *arc_ptr)
{     /* Generate the basic grid. */
      int direction = 1, i, j, k;
      if (two_way)
      {  /* Generate an arc in each direction. */
         for (i = 1; i < n_node; i += n1)
         {  for (j = i, k = j + n1 - 1; j < k; j++)
            {  arc_ptr->from = j;
               arc_ptr->to = j + 1;
               arc_ptr++;
               arc_ptr->from = j + 1;
               arc_ptr->to = j;
               arc_ptr++;
            }
         }
         for (i = 1; i <= n1; i++)
         {  for (j = i + n1; j < n_node; j += n1)
            {  arc_ptr->from = j;
               arc_ptr->to = j - n1;
               arc_ptr++;
               arc_ptr->from = j - n1;
               arc_ptr->to = j;
               arc_ptr++;
            }
         }
      }
      else
      {  /* Generate one arc in each direction. */
         for (i = 1; i < n_node; i += n1)
         {  if (direction == 1)
               j = i;
            else
               j = i + 1;
            for (k = j + n1 - 1; j < k; j++)
            {  arc_ptr->from = j;
               arc_ptr->to = j + direction;
               arc_ptr++;
            }
            direction = - direction;
         }
         for (i = 1; i <= n1; i++)
         {  j = i + n1;
            if (direction == 1)
            {  for (; j < n_node; j += n1)
               {  arc_ptr->from = j - n1;
                  arc_ptr->to = j;
                  arc_ptr++;
               }
            }
            else
            {  for (; j < n_node; j += n1)
               {  arc_ptr->from = j - n1;
                  arc_ptr->to = j;
                  arc_ptr++;
               }
            }
            direction = - direction;
         }
      }
      return arc_ptr;
}

static void gen_more_arcs(struct csa *csa, struct arcs *arc_ptr)
{     /* Generate random arcs to meet the specified density. */
      int i;
      double ab[2];
      ab[0] = 0.9;
      ab[1] = n_node - 0.99;  /* upper limit is n_node-1 because the
                                 supernode cannot be selected */
      for (i = n_grid_arc; i < n_arc; i++, arc_ptr++)
      {  arc_ptr->from = uniform(csa, ab);
         arc_ptr->to = uniform(csa, ab);
         if (arc_ptr->from == arc_ptr->to)
         {  arc_ptr--;
            i--;
         }
      }
      return;
}

static void generate(struct csa *csa)
{     /* Generate a random network. */
      struct arcs *arc_ptr = arc_list;
      arc_ptr = gen_basic_grid(csa, arc_ptr);
      select_source_sinks(csa);
      arc_ptr = gen_additional_arcs(csa, arc_ptr);
      gen_more_arcs(csa, arc_ptr);
      assign_costs(csa);
      assign_capacities(csa);
      assign_imbalance(csa);
      return;
}

static void output(struct csa *csa)
{     /* Output the network in DIMACS format. */
      struct arcs *arc_ptr;
      struct imbalance *imb_ptr;
      int i;
      if (G != NULL) goto skip;
      /* Output "c", "p" records. */
      xprintf("c generated by GRIDGEN\n");
      xprintf("c seed %d\n", seed_original);
      xprintf("c nodes %d\n", n_node);
      xprintf("c grid size %d X %d\n", n1, n2);
      xprintf("c sources %d sinks %d\n", n_source, n_sink);
      xprintf("c avg. degree %d\n", avg_degree);
      xprintf("c supply %d\n", t_supply);
      switch (arc_costs.distribution)
      {  case UNIFORM:
            xprintf("c arc costs: UNIFORM distr. min %d max %d\n",
               (int)arc_costs.parameter[0],
               (int)arc_costs.parameter[1]);
            break;
         case EXPONENTIAL:
            xprintf("c arc costs: EXPONENTIAL distr. lambda %d\n",
               (int)arc_costs.parameter[0]);
            break;
         default:
            xassert(csa != csa);
      }
      switch (capacities.distribution)
      {  case UNIFORM:
            xprintf("c arc caps :  UNIFORM distr. min %d max %d\n",
               (int)capacities.parameter[0],
               (int)capacities.parameter[1]);
            break;
         case EXPONENTIAL:
            xprintf("c arc caps :  EXPONENTIAL distr. %d lambda %d\n",
               (int)capacities.parameter[0]);
            break;
         default:
            xassert(csa != csa);
      }
skip: if (G == NULL)
         xprintf("p min %d %d\n", n_node, n_arc);
      else
      {  glp_add_vertices(G, n_node);
         if (v_rhs >= 0)
         {  double zero = 0.0;
            for (i = 1; i <= n_node; i++)
            {  glp_vertex *v = G->v[i];
               memcpy((char *)v->data + v_rhs, &zero, sizeof(double));
            }
         }
      }
      /* Output "n node supply". */
      for (i = 0, imb_ptr = source_list; i < n_source; i++, imb_ptr++)
      {  if (G == NULL)
            xprintf("n %d %d\n", imb_ptr->node, imb_ptr->supply);
         else
         {  if (v_rhs >= 0)
            {  double temp = (double)imb_ptr->supply;
               glp_vertex *v = G->v[imb_ptr->node];
               memcpy((char *)v->data + v_rhs, &temp, sizeof(double));
            }
         }
      }
      for (i = 0, imb_ptr = sink_list; i < n_sink; i++, imb_ptr++)
      {  if (G == NULL)
            xprintf("n %d %d\n", imb_ptr->node, imb_ptr->supply);
         else
         {  if (v_rhs >= 0)
            {  double temp = (double)imb_ptr->supply;
               glp_vertex *v = G->v[imb_ptr->node];
               memcpy((char *)v->data + v_rhs, &temp, sizeof(double));
            }
         }
      }
      /* Output "a from to lowcap=0 hicap cost". */
      for (i = 0, arc_ptr = arc_list; i < n_arc; i++, arc_ptr++)
      {  if (G == NULL)
            xprintf("a %d %d 0 %d %d\n", arc_ptr->from, arc_ptr->to,
               arc_ptr->u, arc_ptr->cost);
         else
         {  glp_arc *a = glp_add_arc(G, arc_ptr->from, arc_ptr->to);
            if (a_cap >= 0)
            {  double temp = (double)arc_ptr->u;
               memcpy((char *)a->data + a_cap, &temp, sizeof(double));
            }
            if (a_cost >= 0)
            {  double temp = (double)arc_ptr->cost;
               memcpy((char *)a->data + a_cost, &temp, sizeof(double));
            }
         }
      }
      return;
}

static double randy(struct csa *csa)
{     /* Returns a random number between 0.0 and 1.0.
         See Ward Cheney & David Kincaid, "Numerical Mathematics and
         Computing," 2Ed, pp. 335. */
      seed = 16807 * seed % 2147483647;
      if (seed < 0) seed = - seed;
      return seed * 4.6566128752459e-10;
}

static void select_source_sinks(struct csa *csa)
{     /* Randomly select the source nodes and sink nodes. */
      int i, *int_ptr;
      int *temp_list;   /* a temporary list of nodes */
      struct imbalance *ptr;
      double ab[2];     /* parameter for random number generator */
      ab[0] = 0.9;
      ab[1] = n_node - 0.99;  /* upper limit is n_node-1 because the
                                 supernode cannot be selected */
      temp_list = xcalloc(n_node, sizeof(int));
      for (i = 0, int_ptr = temp_list; i < n_node; i++, int_ptr++)
         *int_ptr = 0;
      /* Select the source nodes. */
      for (i = 0, ptr = source_list; i < n_source; i++, ptr++)
      {  ptr->node = uniform(csa, ab);
         if (temp_list[ptr->node] == 1) /* check for duplicates */
         {  ptr--;
            i--;
         }
         else
            temp_list[ptr->node] = 1;
      }
      /* Select the sink nodes. */
      for (i = 0, ptr = sink_list; i < n_sink; i++, ptr++)
      {  ptr->node = uniform(csa, ab);
         if (temp_list[ptr->node] == 1)
         {  ptr--;
            i--;
         }
         else
            temp_list[ptr->node] = 1;
      }
      xfree(temp_list);
      return;
}

int uniform(struct csa *csa, double a[2])
{     /* Generates an integer uniformly selected from [a[0],a[1]]. */
      return (int)((a[1] - a[0]) * randy(csa) + a[0] + 0.5);
}

/**********************************************************************/

#if 0
int main(void)
{     int parm[1+14];
      double temp;
      scanf("%d", &parm[1]);
      scanf("%d", &parm[2]);
      scanf("%d", &parm[3]);
      scanf("%d", &parm[4]);
      scanf("%d", &parm[5]);
      scanf("%d", &parm[6]);
      scanf("%d", &parm[7]);
      scanf("%d", &parm[8]);
      scanf("%d", &parm[9]);
      if (parm[9] == 1)
      {  scanf("%d", &parm[10]);
         scanf("%d", &parm[11]);
      }
      else
      {  scanf("%le", &temp);
         parm[10] = (int)(100.0 * temp + .5);
         parm[11] = 0;
      }
      scanf("%d", &parm[12]);
      if (parm[12] == 1)
      {  scanf("%d", &parm[13]);
         scanf("%d", &parm[14]);
      }
      else
      {  scanf("%le", &temp);
         parm[13] = (int)(100.0 * temp + .5);
         parm[14] = 0;
      }
      glp_gridgen(NULL, 0, 0, 0, parm);
      return 0;
}
#endif

/* eof */
