/* glpnet05.c (Goldfarb's maximum flow problem generator) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*
*  This code is a modified version of the program RMFGEN, a maxflow
*  problem generator developed by D.Goldfarb and M.Grigoriadis, and
*  originally implemented by Tamas Badics <badics@rutcor.rutgers.edu>.
*  The original code is publically available on the DIMACS ftp site at:
*  <ftp://dimacs.rutgers.edu/pub/netflow/generators/network/genrmf>.
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

#include "glpapi.h"
#include "glprng.h"

/***********************************************************************
*  NAME
*
*  glp_rmfgen - Goldfarb's maximum flow problem generator
*
*  SYNOPSIS
*
*  int glp_rmfgen(glp_graph *G, int *s, int *t, int a_cap,
*     const int parm[1+5]);
*
*  DESCRIPTION
*
*  The routine glp_rmfgen is a maximum flow problem generator developed
*  by D.Goldfarb and M.Grigoriadis.
*
*  The parameter G specifies the graph object, to which the generated
*  problem data have to be stored. Note that on entry the graph object
*  is erased with the routine glp_erase_graph.
*
*  The pointer s specifies a location, to which the routine stores the
*  source node number. If s is NULL, the node number is not stored.
*
*  The pointer t specifies a location, to which the routine stores the
*  sink node number. If t is NULL, the node number is not stored.
*
*  The parameter a_cap specifies an offset of the field of type double
*  in the arc data block, to which the routine stores the arc capacity.
*  If a_cap < 0, the capacity is not stored.
*
*  The array parm contains description of the network to be generated:
*
*  parm[0]  not used
*  parm[1]  (seed)   random number seed (a positive integer)
*  parm[2]  (a)      frame size
*  parm[3]  (b)      depth
*  parm[4]  (c1)     minimal arc capacity
*  parm[5]  (c2)     maximal arc capacity
*
*  RETURNS
*
*  If the instance was successfully generated, the routine glp_netgen
*  returns zero; otherwise, if specified parameters are inconsistent,
*  the routine returns a non-zero error code.
*
*  COMMENTS
*
*  The generated network is as follows. It has b pieces of frames of
*  size a * a. (So alltogether the number of vertices is a * a * b)
*
*  In each frame all the vertices are connected with their neighbours
*  (forth and back). In addition the vertices of a frame are connected
*  one to one with the vertices of next frame using a random permutation
*  of those vertices.
*
*  The source is the lower left vertex of the first frame, the sink is
*  the upper right vertex of the b'th frame.
*
*                    t
*           +-------+
*           |      .|
*           |     . |
*        /  |    /  |
*       +-------+/ -+ b
*       |    |  |/.
*     a |   -v- |/
*       |    |  |/
*       +-------+ 1
*      s    a
*
*  The capacities are randomly chosen integers from the range of [c1,c2]
*  in the case of interconnecting edges, and c2 * a * a for the in-frame
*  edges.
*
*  REFERENCES
*
*  D.Goldfarb and M.D.Grigoriadis, "A computational comparison of the
*  Dinic and network simplex methods for maximum flow." Annals of Op.
*  Res. 13 (1988), pp. 83-123.
*
*  U.Derigs and W.Meier, "Implementing Goldberg's max-flow algorithm:
*  A computational investigation." Zeitschrift fuer Operations Research
*  33 (1989), pp. 383-403. */

typedef struct VERTEX
{     struct EDGE **edgelist;
      /* Pointer to the list of pointers to the adjacent edges.
         (No matter that to or from edges) */
      struct EDGE **current;
      /* Pointer to the current edge */
      int degree;
      /* Number of adjacent edges (both direction) */
      int index;
} vertex;

typedef struct EDGE
{     int from;
      int to;
      int cap;
      /* Capacity */
} edge;

typedef struct NETWORK
{     struct NETWORK *next, *prev;
      int vertnum;
      int edgenum;
      vertex *verts;
      /* Vertex array[1..vertnum] */
      edge *edges;
      /* Edge array[1..edgenum] */
      int source;
      /* Pointer to the source */
      int sink;
      /* Pointer to the sink */
} network;

struct csa
{     /* common storage area */
      glp_graph *G;
      int *s, *t, a_cap;
      RNG *rand;
      network *N;
      int *Parr;
      int A, AA, C2AA, Ec;
};

#define G      (csa->G)
#define s      (csa->s)
#define t      (csa->t)
#define a_cap  (csa->a_cap)
#define N      (csa->N)
#define Parr   (csa->Parr)
#define A      (csa->A)
#define AA     (csa->AA)
#define C2AA   (csa->C2AA)
#define Ec     (csa->Ec)

#undef random
#define random(A) (int)(rng_unif_01(csa->rand) * (double)(A))
#define RANDOM(A, B) (int)(random((B) - (A) + 1) + (A))
#define sgn(A) (((A) > 0) ? 1 : ((A) == 0) ? 0 : -1)

static void make_edge(struct csa *csa, int from, int to, int c1, int c2)
{     Ec++;
      N->edges[Ec].from = from;
      N->edges[Ec].to = to;
      N->edges[Ec].cap = RANDOM(c1, c2);
      return;
}

static void permute(struct csa *csa)
{     int i, j, tmp;
      for (i = 1; i < AA; i++)
      {  j = RANDOM(i, AA);
         tmp = Parr[i];
         Parr[i] = Parr[j];
         Parr[j] = tmp;
      }
      return;
}

static void connect(struct csa *csa, int offset, int cv, int x1, int y1)
{     int cv1;
      cv1 = offset + (x1 - 1) * A + y1;
      Ec++;
      N->edges[Ec].from = cv;
      N->edges[Ec].to = cv1;
      N->edges[Ec].cap = C2AA;
      return;
}

static network *gen_rmf(struct csa *csa, int a, int b, int c1, int c2)
{     /* generates a network with a*a*b nodes and 6a*a*b-4ab-2a*a edges
         random_frame network:
         Derigs & Meier, Methods & Models of OR (1989), 33:383-403 */
      int x, y, z, offset, cv;
      A = a;
      AA = a * a;
      C2AA = c2 * AA;
      Ec = 0;
      N = (network *)xmalloc(sizeof(network));
      N->vertnum = AA * b;
      N->edgenum = 5 * AA * b - 4 * A * b - AA;
      N->edges = (edge *)xcalloc(N->edgenum + 1, sizeof(edge));
      N->source = 1;
      N->sink = N->vertnum;
      Parr = (int *)xcalloc(AA + 1, sizeof(int));
      for (x = 1; x <= AA; x++)
         Parr[x] = x;
      for (z = 1; z <= b; z++)
      {  offset = AA * (z - 1);
         if (z != b)
            permute(csa);
         for (x = 1; x <= A; x++)
         {  for (y = 1; y <= A; y++)
            {  cv = offset + (x - 1) * A + y;
               if (z != b)
                  make_edge(csa, cv, offset + AA + Parr[cv - offset],
                     c1, c2); /* the intermediate edges */
               if (y < A)
                  connect(csa, offset, cv, x, y + 1);
               if (y > 1)
                  connect(csa, offset, cv, x, y - 1);
               if (x < A)
                  connect(csa, offset, cv, x + 1, y);
               if (x > 1)
                  connect(csa, offset, cv, x - 1, y);
            }
         }
      }
      xfree(Parr);
      return N;
}

static void print_max_format(struct csa *csa, network *n, char *comm[],
      int dim)
{     /* prints a network heading with dim lines of comments (no \n
         needs at the ends) */
      int i, vnum, e_num;
      edge *e;
      vnum = n->vertnum;
      e_num = n->edgenum;
      if (G == NULL)
      {  for (i = 0; i < dim; i++)
            xprintf("c %s\n", comm[i]);
         xprintf("p max %7d %10d\n", vnum, e_num);
         xprintf("n %7d s\n", n->source);
         xprintf("n %7d t\n", n->sink);
      }
      else
      {  glp_add_vertices(G, vnum);
         if (s != NULL) *s = n->source;
         if (t != NULL) *t = n->sink;
      }
      for (i = 1; i <= e_num; i++)
      {  e = &n->edges[i];
         if (G == NULL)
            xprintf("a %7d %7d %10d\n", e->from, e->to, (int)e->cap);
         else
         {  glp_arc *a = glp_add_arc(G, e->from, e->to);
            if (a_cap >= 0)
            {  double temp = (double)e->cap;
               memcpy((char *)a->data + a_cap, &temp, sizeof(double));
            }
         }
      }
      return;
}

static void gen_free_net(network *n)
{     xfree(n->edges);
      xfree(n);
      return;
}

int glp_rmfgen(glp_graph *G_, int *_s, int *_t, int _a_cap,
      const int parm[1+5])
{     struct csa _csa, *csa = &_csa;
      network *n;
      char comm[10][80], *com1[10];
      int seed, a, b, c1, c2, ret;
      G = G_;
      s = _s;
      t = _t;
      a_cap = _a_cap;
      if (G != NULL)
      {  if (a_cap >= 0 && a_cap > G->a_size - (int)sizeof(double))
           xerror("glp_rmfgen: a_cap = %d; invalid offset\n", a_cap);
      }
      seed = parm[1];
      a = parm[2];
      b = parm[3];
      c1 = parm[4];
      c2 = parm[5];
      if (!(seed > 0 && 1 <= a && a <= 1000 && 1 <= b && b <= 1000 &&
            0 <= c1 && c1 <= c2 && c2 <= 1000))
      {  ret = 1;
         goto done;
      }
      if (G != NULL)
      {  glp_erase_graph(G, G->v_size, G->a_size);
         glp_set_graph_name(G, "RMFGEN");
      }
      csa->rand = rng_create_rand();
      rng_init_rand(csa->rand, seed);
      n = gen_rmf(csa, a, b, c1, c2);
      sprintf(comm[0], "This file was generated by genrmf.");
      sprintf(comm[1], "The parameters are: a: %d b: %d c1: %d c2: %d",
         a, b, c1, c2);
      com1[0] = comm[0];
      com1[1] = comm[1];
      print_max_format(csa, n, com1, 2);
      gen_free_net(n);
      rng_delete_rand(csa->rand);
      ret = 0;
done: return ret;
}

/**********************************************************************/

#if 0
int main(int argc, char *argv[])
{     int seed, a, b, c1, c2, i, parm[1+5];
      seed = 123;
      a = b = c1 = c2 = -1;
      for (i = 1; i < argc; i++)
      {  if (strcmp(argv[i], "-seed") == 0)
            seed = atoi(argv[++i]);
         else if (strcmp(argv[i], "-a") == 0)
            a = atoi(argv[++i]);
         else if (strcmp(argv[i], "-b") == 0)
            b = atoi(argv[++i]);
         else if (strcmp(argv[i], "-c1") == 0)
            c1 = atoi(argv[++i]);
         else if (strcmp(argv[i], "-c2") == 0)
            c2 = atoi(argv[++i]);
      }
      if (a < 0 || b < 0 || c1 < 0 || c2 < 0)
      {  xprintf("Usage:\n");
         xprintf("genrmf [-seed seed] -a frame_size -b depth\n");
         xprintf("        -c1 cap_range1 -c2 cap_range2\n");
      }
      else
      {  parm[1] = seed;
         parm[2] = a;
         parm[3] = b;
         parm[4] = c1;
         parm[5] = c2;
         glp_rmfgen(NULL, NULL, NULL, 0, parm);
      }
      return 0;
}
#endif

/* eof */
