/* -*- mode: C -*-  */
/* 
   IGraph R library.
   Copyright (C) 2003, 2004  Gabor Csardi <csardi@rmki.kfki.hu>
   MTA RMKI, Konkoly-Thege Miklos st. 29-33, Budapest 1121, Hungary
   
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
   02110-1301 USA

*/

#include "igraph.h"
#include "random.h"
#include "memory.h"

#include <math.h>

/**
 * \section about_games
 * 
 * <para>Games are randomized graph generators. Randomization means that
 * they generate a different graph every time you call them. </para>
 */

/**
 * \ingroup generators
 * \function igraph_barabasi_game
 * \brief Generates a graph based on the Barab&aacute;si-Albert model.
 *
 * \param graph An uninitialized graph object.
 * \param n The number of vertices in the graph.
 * \param m The number of outgoing edges generated for each 
 *        vertex. (Only if \p outseq is \c NULL.) 
 * \param outseq Gives the (out-)degrees of the vertices. If this is
 *        constant, this can be a NULL pointer or an empty (but
 *        initialized!) vector, in this case \p m contains
 *        the constant out-degree. The very first vertex has by definition 
 *        no outgoing edges, so the first number in this vector is 
 *        ignored.
 * \param outpref Boolean, if true not only the in- but also the out-degree
 *        of a vertex increases its citation probability. Ie. the
 *        citation probability is determined by the total degree of
 *        the vertices.
 * \param directed Boolean, whether to generate a directed graph.
 * \return Error code:
 *         \c IGRAPH_EINVAL: invalid \p n,
 *         \p m or \p outseq parameter.
 * 
 * Time complexity: O(|V|+|E|), the
 * number of vertices plus the number of edges.
 */

int igraph_barabasi_game(igraph_t *graph, igraph_integer_t n, igraph_integer_t m, 
			 const igraph_vector_t *outseq, igraph_bool_t outpref, 
			 igraph_bool_t directed) {

  long int no_of_nodes=n;
  long int no_of_neighbors=m;
  long int *bag;
  long int bagp;
  long int no_of_edges;
  igraph_vector_t edges=IGRAPH_VECTOR_NULL;
  
  long int resp=0;

  long int i,j;

  if (n<0) {
    IGRAPH_ERROR("Invalid number of vertices", IGRAPH_EINVAL);
  }
  if (outseq != 0 && igraph_vector_size(outseq) != 0 && igraph_vector_size(outseq) != n) {
    IGRAPH_ERROR("Invalid out degree sequence length", IGRAPH_EINVAL);
  }
  if ( (outseq == 0 || igraph_vector_size(outseq) == 0) && m<0) {
    IGRAPH_ERROR("Invalid out degree", IGRAPH_EINVAL);
  }

  if (outseq==0 || igraph_vector_size(outseq) == 0) {
    no_of_neighbors=m;
    bag=Calloc(no_of_nodes * no_of_neighbors + no_of_nodes +
	       outpref * no_of_nodes * no_of_neighbors,
	       long int);
    no_of_edges=(no_of_nodes-1)*no_of_neighbors;
  } else {
    no_of_edges=0;
    for (i=1; i<igraph_vector_size(outseq); i++) {
      no_of_edges+=VECTOR(*outseq)[i];
    }
    bag=Calloc(no_of_nodes + no_of_edges + outpref * no_of_edges,
	       long int);
  }
  
  if (bag==0) {
    IGRAPH_ERROR("barabasi_game failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(free, bag); 	/* TODO: hack */
  IGRAPH_VECTOR_INIT_FINALLY(&edges, no_of_edges*2);

  /* The first node */

  bagp=0;
  bag[bagp++]=0;
  
  RNG_BEGIN();

  /* and the others */
  
  for (i=1; i<no_of_nodes; i++) {
    /* draw edges */
    if (outseq != 0 && igraph_vector_size(outseq)!=0) { no_of_neighbors=VECTOR(*outseq)[i]; }
    for (j=0; j<no_of_neighbors; j++) {
      long int to=bag[RNG_INTEGER(0, bagp-1)];
      VECTOR(edges)[resp++] = i;
      VECTOR(edges)[resp++] = to;
    }
    /* update bag */
    bag[bagp++] = i;
    for (j=0; j<no_of_neighbors; j++) {
      bag[bagp++] = VECTOR(edges)[resp-2*j-1];
      if (outpref) {
	bag[bagp++] = i;
      }
    }
  }

  RNG_END();

  Free(bag);
  IGRAPH_CHECK(igraph_create(graph, &edges, 0, directed));
  igraph_vector_destroy(&edges);
  IGRAPH_FINALLY_CLEAN(2);

  return 0;
}

/**
 * \ingroup internal
 */

int igraph_erdos_renyi_game_gnp(igraph_t *graph, igraph_integer_t n, igraph_real_t p,
				igraph_bool_t directed, igraph_bool_t loops) {

  long int no_of_nodes=n;
  igraph_vector_t edges=IGRAPH_VECTOR_NULL;
  igraph_vector_t s=IGRAPH_VECTOR_NULL;
  int retval=0;

  if (n<0) {
    IGRAPH_ERROR("Invalid number of vertices", IGRAPH_EINVAL);
  }
  if (p<0.0 || p>1.0) {
    IGRAPH_ERROR("Invalid probability given", IGRAPH_EINVAL);
  }
  
  if (p==0.0 || no_of_nodes<=1) {
    IGRAPH_CHECK(retval=igraph_empty(graph, n, directed));
  } else if (p==1.0) { 
    IGRAPH_CHECK(retval=igraph_full(graph, n, directed, loops));
  } else {

    long int i;
    double maxedges;
    if (directed && loops) 
      { maxedges = n * n; }
    else if (directed && !loops)
      { maxedges = n * (n-1); }
    else if (!directed && loops) 
      { maxedges = n * (n+1)/2; }
    else 
      { maxedges = n * (n-1)/2; }
    
    IGRAPH_VECTOR_INIT_FINALLY(&s, 0);
    IGRAPH_CHECK(igraph_vector_reserve(&s, maxedges*p*1.1));

    RNG_BEGIN();

    IGRAPH_CHECK(igraph_vector_push_back(&s, RNG_GEOM(p)+1));
    while (igraph_vector_tail(&s) < maxedges) {
      IGRAPH_CHECK(igraph_vector_push_back(&s, igraph_vector_tail(&s)+RNG_GEOM(p)+1));
    }
    if (igraph_vector_tail(&s) > maxedges+1) {
      igraph_vector_pop_back(&s);
    }

    RNG_END();

    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
    IGRAPH_CHECK(igraph_vector_reserve(&edges, igraph_vector_size(&s)*2));

    if (directed && loops) {
      for (i=0; i<igraph_vector_size(&s); i++) {
	igraph_vector_push_back(&edges, ((long int)(VECTOR(s)[i]-1))/no_of_nodes);
	igraph_vector_push_back(&edges, ((long int)(VECTOR(s)[i]-1))%no_of_nodes);
      }
    } else if (directed && !loops) {
      for (i=0; i<igraph_vector_size(&s); i++) {
	long int from=((long int)(VECTOR(s)[i]-1))/(no_of_nodes-1);
	long int to=((long int)VECTOR(s)[i])%(no_of_nodes-1);
	if (from==to) {
	  to=no_of_nodes-1;
	}
	igraph_vector_push_back(&edges, from);
	igraph_vector_push_back(&edges, to);
      }
    } else if (!directed && loops) {
      for (i=0; i<igraph_vector_size(&s); i++) {
	igraph_real_t from=ceil((sqrt(8*(VECTOR(s)[i])+1)-1)/2);
	igraph_vector_push_back(&edges, from-1);
	igraph_vector_push_back(&edges, VECTOR(s)[i]-from*(from-1)/2-1);
      }
    } else {
      for (i=0; i<igraph_vector_size(&s); i++) {
	igraph_real_t from=ceil((sqrt(8*VECTOR(s)[i]+1)-1)/2)+1;
	igraph_vector_push_back(&edges, from-1);
	igraph_vector_push_back(&edges, VECTOR(s)[i]-(from-1)*(from-2)/2-1);
      }
    }      

    igraph_vector_destroy(&s);
    IGRAPH_FINALLY_CLEAN(1);
    IGRAPH_CHECK(retval=igraph_create(graph, &edges, n, directed));
    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);
  }

  return retval;
}

int igraph_erdos_renyi_game_gnm(igraph_t *graph, igraph_integer_t n, igraph_real_t m,
				igraph_bool_t directed, igraph_bool_t loops) {

  long int no_of_nodes=n;
  long int no_of_edges=m;
  igraph_vector_t edges=IGRAPH_VECTOR_NULL;
  igraph_vector_t s=IGRAPH_VECTOR_NULL;
  int retval=0;

  if (n<0) {
    IGRAPH_ERROR("Invalid number of vertices", IGRAPH_EINVAL);
  }
  if (m<0) {
    IGRAPH_ERROR("Invalid number of edges", IGRAPH_EINVAL);
  }
  
  if (m==0.0 || no_of_nodes<=1) {
    IGRAPH_CHECK(retval=igraph_empty(graph, n, directed));
  } else {

    long int i;    
    double maxedges;
    if (directed && loops) 
      { maxedges = n * n; }
    else if (directed && !loops)
      { maxedges = n * (n-1); }
    else if (!directed && loops) 
      { maxedges = n * (n+1)/2; }
    else 
      { maxedges = n * (n-1)/2; }
    
    if (no_of_edges > maxedges) {
      IGRAPH_ERROR("Invalid number (too large) of edges", IGRAPH_EINVAL);
    }
    
    if (maxedges == no_of_edges) {
      retval=igraph_full(graph, n, directed, loops);
    } else {
    
      IGRAPH_VECTOR_INIT_FINALLY(&s, 0);
      IGRAPH_CHECK(igraph_random_sample(&s, 1, maxedges, no_of_edges));
      
      IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
      IGRAPH_CHECK(igraph_vector_reserve(&edges, igraph_vector_size(&s)*2));
      
      if (directed && loops) {
	for (i=0; i<igraph_vector_size(&s); i++) {
	  igraph_vector_push_back(&edges, ((long int)(VECTOR(s)[i]-1))/no_of_nodes);
	  igraph_vector_push_back(&edges, ((long int)(VECTOR(s)[i]-1))%no_of_nodes);
	}
      } else if (directed && !loops) {
	for (i=0; i<igraph_vector_size(&s); i++) {
	  long int from=((long int)(VECTOR(s)[i]-1))/(no_of_nodes-1);
	  long int to=((long int)VECTOR(s)[i])%(no_of_nodes-1);
	  if (from==to) {
	    to=no_of_nodes-1;
	  }
	  igraph_vector_push_back(&edges, from);
	  igraph_vector_push_back(&edges, to);
	}
      } else if (!directed && loops) {
	for (i=0; i<igraph_vector_size(&s); i++) {
	  igraph_real_t from=ceil((sqrt(8*(VECTOR(s)[i])+1)-1)/2);
	  igraph_vector_push_back(&edges, from-1);
	  igraph_vector_push_back(&edges, VECTOR(s)[i]-from*(from-1)/2-1);
	}
      } else {
	for (i=0; i<igraph_vector_size(&s); i++) {
	  igraph_real_t from=ceil((sqrt(8*VECTOR(s)[i]+1)-1)/2)+1;
	  igraph_vector_push_back(&edges, from-1);
	  igraph_vector_push_back(&edges, VECTOR(s)[i]-(from-1)*(from-2)/2-1);
	}
      }  

      igraph_vector_destroy(&s);
      retval=igraph_create(graph, &edges, n, directed);
      igraph_vector_destroy(&edges);
      IGRAPH_FINALLY_CLEAN(2);
    }
  }
  
  return retval;
}

/**
 * \ingroup generators
 * \function igraph_erdos_renyi_game
 * \brief Generates a random (Erdos-Renyi) graph.
 * 
 * \param graph Pointer to an uninitialized graph object.
 * \param type The type of the random graph, possible values:
 *        \clist
 *        \cli IGRAPH_ERDOS_RENYI_GNM
 *          G(n,m) graph,  
 *          m edges are
 *          selected uniformly randomly in a graph with
 *          n vertices.
 *        \cli IGRAPH_ERDOS_RENYI_GNP
 *          G(n,p) graph,
 *          every possible edge is included in the graph with
 *          probability p.
 *        \endclist
 * \param n The number of vertices in the graph.
 * \param p_or_m This is the p parameter for
 *        G(n,p) graphs and the
 *        m 
 *        parameter for G(n,m) graphs.
 * \param directed Logical, whether to generate a directed graph.
 * \param loops Logical, whether to generate loops (self) edges.
 * \return Error code:
 *         \c IGRAPH_EINVAL: invalid
 *         \p type, \p n,
 *         \p p or \p m
 *          parameter.
 *         \c IGRAPH_ENOMEM: there is not enought
 *         memory for the operation.
 * 
 * Time complexity: O(|V|+|E|), the
 * number of vertices plus the number of edges in the graph.
 * 
 * \sa \ref igraph_barabasi_game(), \ref igraph_growing_random_game()
 */

int igraph_erdos_renyi_game(igraph_t *graph, igraph_erdos_renyi_t type,
			    igraph_integer_t n, igraph_real_t p_or_m,
			    igraph_bool_t directed, igraph_bool_t loops) {
  int retval=0;
  if (type == IGRAPH_ERDOS_RENYI_GNP) {
    retval=igraph_erdos_renyi_game_gnp(graph, n, p_or_m, directed, loops);
  } else if (type == IGRAPH_ERDOS_RENYI_GNM) {
    retval=igraph_erdos_renyi_game_gnm(graph, n, p_or_m, directed, loops);
  } else {
    IGRAPH_ERROR("Invalid type", IGRAPH_EINVAL);
  }
  
  return retval;
}

int igraph_degree_sequence_game_simple(igraph_t *graph, 
				       const igraph_vector_t *out_seq, 
				       const igraph_vector_t *in_seq) {

  long int outsum=0, insum=0;
  igraph_bool_t directed=(in_seq != 0 && igraph_vector_size(in_seq)!=0);
  long int no_of_nodes, no_of_edges;
  long int *bag1=0, *bag2=0;
  long int bagp1=0, bagp2=0;
  igraph_vector_t edges=IGRAPH_VECTOR_NULL;
  long int i,j;

  if (igraph_vector_any_smaller(out_seq, 0)) {
    IGRAPH_ERROR("Negative out degree", IGRAPH_EINVAL);
  }
  if (directed && igraph_vector_any_smaller(in_seq, 0)) {
    IGRAPH_ERROR("Negative in degree", IGRAPH_EINVAL);
  }
  if (directed && 
      igraph_vector_size(out_seq) != igraph_vector_size(in_seq)) { 
    IGRAPH_ERROR("Length of `out_seq' and `in_seq' differ for directed graph",
		 IGRAPH_EINVAL);
  }
  
  outsum=igraph_vector_sum(out_seq);
  insum=igraph_vector_sum(in_seq);
  
  if (!directed && outsum % 2 != 0) {
    IGRAPH_ERROR("Total degree not even for undirected graph", IGRAPH_EINVAL);
  }
  
  if (directed && outsum != insum) {
    IGRAPH_ERROR("Total in-degree and out-degree differ for directed graph",
		  IGRAPH_EINVAL);
  }

  no_of_nodes=igraph_vector_size(out_seq);
  if (directed) {
    no_of_edges=outsum;
  } else {
    no_of_edges=outsum/2;
  }

  bag1=Calloc(outsum, long int);
  if (bag1==0) {
    IGRAPH_ERROR("degree sequence game (simple)", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(free, bag1); 	/* TODO: hack */
    
  for (i=0; i<no_of_nodes; i++) {
    for (j=0; j<VECTOR(*out_seq)[i]; j++) {
      bag1[bagp1++]=i;
    }
  }
  if (directed) {
    bag2=Calloc(insum, long int);
    if (bag2==0) {
      IGRAPH_ERROR("degree sequence game (simple)", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(free, bag2);
    for (i=0; i<no_of_nodes; i++) {
      for (j=0; j<VECTOR(*in_seq)[i]; j++) {
	bag2[bagp2++]=i;
      }
    }
  }

  IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
  IGRAPH_CHECK(igraph_vector_reserve(&edges, no_of_edges));

  RNG_BEGIN();

  if (directed) {
    for (i=0; i<no_of_edges; i++) {
      long int from=RNG_INTEGER(0, bagp1-1);
      long int to=RNG_INTEGER(0, bagp2-1);
      igraph_vector_push_back(&edges, bag1[from]); /* safe, already reserved */
      igraph_vector_push_back(&edges, bag2[to]);   /* detto */
      bag1[from]=bag1[bagp1-1];
      bag2[to]=bag2[bagp2-1];
      bagp1--; bagp2--;
    }
  } else {
    for (i=0; i<no_of_edges; i++) {
      long int from=RNG_INTEGER(0, bagp1-1);
      long int to;
      igraph_vector_push_back(&edges, bag1[from]); /* safe, already reserved */
      bag1[from]=bag1[bagp1-1];
      bagp1--;
      to=RNG_INTEGER(0, bagp1-1);
      igraph_vector_push_back(&edges, bag1[to]);   /* detto */
      bag1[to]=bag1[bagp1-1];
      bagp1--;
    }
  }
  
  RNG_END();

  Free(bag1);
  IGRAPH_FINALLY_CLEAN(1);
  if (directed) {
    Free(bag2);
    IGRAPH_FINALLY_CLEAN(1);
  }

  IGRAPH_CHECK(igraph_create(graph, &edges, no_of_nodes, directed));
  igraph_vector_destroy(&edges);
  IGRAPH_FINALLY_CLEAN(1);
  
  return 0;
}

/**
 * \ingroup generators
 * \function igraph_degree_sequence_game
 * \brief Generates a random graph with a given degree sequence 
 * 
 * \param graph Pointer to an uninitialized graph object.
 * \param out_deg The degree sequence for an undirected graph (if
 *        \p in_seq is of length zero), or the out-degree
 *        sequence of a directed graph (if \p in_deq is not
 *        of length zero.
 * \param in_deg It is either a zero-length vector or
 *        \c NULL (if an undirected 
 *        graph is generated), or the in-degree sequence.
 * \param method The method to generate the graph. Possible values: 
 *        \c IGRAPH_DEGSEQ_SIMPLE, for undirected graphs this
 *        method puts all vertex ids in a bag, the multiplicity of a
 *        vertex in the bag is the same as its degree. Then it 
 *        draws pairs from the bag, until it is empty. This method can 
 *        generate both loop (self) edges and multiple edges.
 *        For directed graphs, the algorithm is basically the same,
 *        but two separate bags are used for the in- and out-degrees. 
 * \return Error code: 
 *          \c IGRAPH_ENOMEM: there is not enough
 *           memory to perform the operation.
 *          \c IGRAPH_EINVAL: invalid method parameter, or
 *           invalid in- and/or out-degree vectors. The degree vectors
 *           should be non-negative, \p out_deg should sum
 *           up to an even integer for undirected graphs; the length
 *           and sum of \p out_deg and
 *           \p in_deg 
 *           should match for directed graphs.
 * 
 * Time complexity: O(|V|+|E|), the
 * number of vertices plus the number of edges.
 * 
 * \sa \ref igraph_barabasi_game(), \ref igraph_erdos_renyi_game()
 */

int igraph_degree_sequence_game(igraph_t *graph, const igraph_vector_t *out_deg,
				const igraph_vector_t *in_deg, 
				igraph_degseq_t method) {

  int retval;

  if (method==IGRAPH_DEGSEQ_SIMPLE) {
    retval=igraph_degree_sequence_game_simple(graph, out_deg, in_deg);
  } else {
    IGRAPH_ERROR("Invalid degree sequence game method", IGRAPH_EINVAL);
  }

  return retval;
}

/**
 * \ingroup generators
 * \function igraph_growing_random_game
 * \brief Generates a growing random graph.
 *
 * </para><para>
 * This function simulates a growing random graph. In each discrete
 * time step a new vertex is added and a number of new edges are also
 * added. These graphs are known to be different from standard (not
 * growing) random graphs.
 * \param graph Uninitialized graph object.
 * \param n The number of vertices in the graph.
 * \param m The number of edges to add in a time step (ie. after
 *        adding a vertex).
 * \param directed Boolean, whether to generate a directed graph.
 * \param citation Boolean, if \c TRUE, the edges always
 *        originate from the most recently added vertex.
 * \return Error code:
 *          \c IGRAPH_EINVAL: invalid
 *          \p n or \p m
 *          parameter. 
 *
 * Time complexity: O(|V|+|E|), the
 * number of vertices plus the number of edges.
 */
int igraph_growing_random_game(igraph_t *graph, igraph_integer_t n, 
			       igraph_integer_t m, igraph_bool_t directed,
			       igraph_bool_t citation) {

  long int no_of_nodes=n;
  long int no_of_neighbors=m;
  long int no_of_edges;
  igraph_vector_t edges=IGRAPH_VECTOR_NULL;
  
  long int resp=0;

  long int i,j;

  if (n<0) {
    IGRAPH_ERROR("Invalid number of vertices", IGRAPH_EINVAL);
  }
  if (m<0) {
    IGRAPH_ERROR("Invalid number of edges per step (m)", IGRAPH_EINVAL);
  }

  no_of_edges=(no_of_nodes-1) * no_of_neighbors;

  IGRAPH_VECTOR_INIT_FINALLY(&edges, no_of_edges*2);  

  RNG_BEGIN();

  for (i=1; i<no_of_nodes; i++) {
    for (j=0; j<no_of_neighbors; j++) {
      if (citation) {
	long int to=RNG_INTEGER(0, i-1);
	VECTOR(edges)[resp++] = i;
	VECTOR(edges)[resp++] = to;
      } else {
	long int from=RNG_INTEGER(0, i);
	long int to=RNG_INTEGER(1,i);
	VECTOR(edges)[resp++] = from;
	VECTOR(edges)[resp++] = to;
      }
    }
  }

  RNG_END();

  IGRAPH_CHECK(igraph_create(graph, &edges, n, directed));
  igraph_vector_destroy(&edges);
  IGRAPH_FINALLY_CLEAN(1);

  return 0;
}

/**
 * \function igraph_callaway_traits_game
 * 
 * </para><para>
 * This function simulates a growing network with vertex types.
 * The different types of vertices prefer to connect other types of
 * vertices with a given probability.</para><para>
 * 
 * </para><para>
 * The simulation goes like this: in each discrete time step a new
 * vertex is added to the graph. The type of this vertex is generated
 * based on \p type_dist. Then two vertices are selected uniformly
 * randomly from the graph. The probability that they will be
 * connected depends on the types of these vertices and is taken from
 * \p pref_matrix. Then another two vertices are selected and this is
 * repeated \p edges_per_step times in each time step.
 * \param graph Pointer to an uninitialized graph.
 * \param nodes The number of nodes in the graph.
 * \param types Number of node types.
 * \param edges_per_step The number of edges to be add per time step.
 * \param type_dist Vector giving the distribution of the vertex
 * types.
 * \param pref_matrix Matrix giving the connection probabilities for
 * the vertex types.
 * \param directed Logical, whether to generate a directed graph.
 * \return Error code. 
 * 
 * Added in version 0.2.</para><para>
 * 
 * Time complexity: O(|V|e*log(|V|)), |V| is the number of vertices, e
 * is \p edges_per_step.
 */

int igraph_callaway_traits_game (igraph_t *graph, igraph_integer_t nodes, 
				igraph_integer_t types, igraph_integer_t edges_per_step, 
				igraph_vector_t *type_dist,
				igraph_matrix_t *pref_matrix,
				igraph_bool_t directed) {
  long int i, j;
  igraph_vector_t edges;
  igraph_vector_t cumdist;
  igraph_real_t maxcum;

  /* TODO: parameter checks */

  IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&cumdist, types+1);
  
  VECTOR(cumdist)[0]=0;
  for (i=0; i<types; i++) {
    VECTOR(cumdist)[i+1] = VECTOR(cumdist)[i]+VECTOR(*type_dist)[i];
  }
  maxcum=igraph_vector_tail(&cumdist);

  RNG_BEGIN();

  for (i=1; i<nodes; i++) {
    for (j=0; j<edges_per_step; j++) {
      long int node1=RNG_INTEGER(0, i);
      long int node2=RNG_INTEGER(0, i);
      igraph_real_t uni1=RNG_UNIF(0, maxcum), uni2=RNG_UNIF(0, maxcum);
      long int type1, type2;
      igraph_vector_binsearch(&cumdist, uni1, &type1);
      igraph_vector_binsearch(&cumdist, uni2, &type2);
/*    printf("unif: %f, %f, types: %li, %li\n", uni1, uni2, type1, type2); */
      if (RNG_UNIF01() < MATRIX(*pref_matrix, type1, type2)) {
	IGRAPH_CHECK(igraph_vector_push_back(&edges, node1));
	IGRAPH_CHECK(igraph_vector_push_back(&edges, node2));
      }
    }
  }

  RNG_END();

  igraph_vector_destroy(&cumdist);
  IGRAPH_FINALLY_CLEAN(1);
  IGRAPH_CHECK(igraph_create(graph, &edges, nodes, directed));
  igraph_vector_destroy(&edges);
  IGRAPH_FINALLY_CLEAN(1);
  return 0;
}

/**
 * \function igraph_establishment_game
 * \brief Generates a graph with a simple growing model with vertex
 * types 
 * 
 * </para><para>
 * The simulation goes like this: a single vertex is added at each
 * time step. This new vertex tries to connect to \p k vertices in the
 * graph. The probability that such a connection is realized depends
 * on the types of the vertices involved. 
 * 
 * \param graph Pointer to an uninitialized graph.
 * \param nodes The number of vertices in the graph.
 * \param types The number of vertex types.
 * \param k The number of connections tried in each time step.
 * \param type_dist Vector giving the distribution of vertex types.
 * \param pref_matrix Matrix giving the connection probabilities for
 * different vertex types.
 * \param directed Logical, whether to generate a directed graph.
 * \return Error code.
 *
 * Added in version 0.2.</para><para>
 *
 * Time complexity: O(|V|*k*log(|V|)), |V| is the number of vertices
 * and k is the \p k parameter.
 */

int igraph_establishment_game(igraph_t *graph, igraph_integer_t nodes,
			      igraph_integer_t types, igraph_integer_t k,
			      igraph_vector_t *type_dist,
			      igraph_matrix_t *pref_matrix,
			      igraph_bool_t directed) {
  
  long int i, j;
  igraph_vector_t edges;
  igraph_vector_t cumdist;
  igraph_vector_t potneis;
  igraph_real_t maxcum;
  
  IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&cumdist, types+1);
  IGRAPH_VECTOR_INIT_FINALLY(&potneis, k);
  
  VECTOR(cumdist)[0]=0;
  for (i=0; i<types; i++) {
    VECTOR(cumdist)[i+1] = VECTOR(cumdist)[i]+VECTOR(*type_dist)[i];
  }
  maxcum=igraph_vector_tail(&cumdist);

  RNG_BEGIN();
  
  for (i=k; i<nodes; i++) {
    igraph_random_sample(&potneis, 0, i-1, k);
    for (j=0; j<k; j++) {
      long int type1, type2;
      igraph_real_t uni1=RNG_UNIF(0, maxcum), uni2=RNG_UNIF(0, maxcum);
      igraph_vector_binsearch(&cumdist, uni1, &type1);
      igraph_vector_binsearch(&cumdist, uni2, &type2);
      if (RNG_UNIF01() < MATRIX(*pref_matrix, type1, type2)) {
	IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
	IGRAPH_CHECK(igraph_vector_push_back(&edges, VECTOR(potneis)[j]));
      }
    }
  }
  
  RNG_END();
  
  igraph_vector_destroy(&potneis);
  igraph_vector_destroy(&cumdist);
  IGRAPH_FINALLY_CLEAN(2);
  IGRAPH_CHECK(igraph_create(graph, &edges, nodes, directed));
  igraph_vector_destroy(&edges);
  IGRAPH_FINALLY_CLEAN(1);
  return 0;
}

/**
 * \function igraph_nonlinear_barabasi_game
 * \brief Generates graph with non-linear preferential attachment
 * 
 * </para><para>
 * This function is very similar to \ref igraph_barabasi_game(), only
 * in this game the probability that a new vertex attaches to a given
 * old vertex is not proportional to the degree of the old node, but
 * some power of the degree of the old node.
 * 
 * </para><para>
 * More precisely the attachment probability is the degree to the
 * power of \p power plus \p zeroappeal.
 * 
 * </para><para>
 * This function might generate graphs with multiple edges if the
 * value of \p m is at least two. You can call \ref igraph_simplify()
 * to get rid of the multiple edges. 
 * \param graph Pointer to an uninitialized graph object, the
 *        generated graph will be stored here.
 * \param n The number of vertices in the generated graph.
 * \param power The power of the preferential attachment.
 * \param m The number of edges to generate in each time step, if the
 *        \p outseq parameter is a null vector or a vector with length
 *        zero. It is ignored otherwise.
 * \param outseq The number of edges to generate in each time
 *        step. For directed graphs this is exactly the out-degree of
 *        the vertices. The first element of the vector is ignored. If
 *        this is a null vector or a vector of length zero then it is
 *        ignored and the value of the \p m argument is used.
 * \param outpref Logical constant, if TRUE then the preferential
 *        attachment is based on the total degree of the nodes instead
 *        of the in-degree.
 * \param zeroappeal Positive number, the attachment probability for
 *        vertices with degree zero. 
 * \param directed Logical constant, whether to generate a directed
 *        graph.
 * \return Error code.
 * 
 * Time complexity: O(|V|*m*log(|V|)+|E|), |V| is the number of
 * vertices, |E| is the number of edges and m is the average number of
 * edges added in a time step.
 * 
 * \sa \ref igraph_barabasi_game() for the slightly more efficient
 * implementation of the special case \p power=1.
 */

int igraph_nonlinear_barabasi_game(igraph_t *graph, igraph_integer_t n,
				   igraph_real_t power,
				   igraph_integer_t m,  
				   const igraph_vector_t *outseq,
				   igraph_bool_t outpref,
				   igraph_real_t zeroappeal,
				   igraph_bool_t directed) {
  long int no_of_nodes=n;
  long int no_of_neighbors=m;
  long int no_of_edges;
  igraph_vector_t edges;
  long int i, j;
  igraph_psumtree_t sumtree;
  long int edgeptr=0;
  igraph_vector_t degree;

  if (n<0) {
    IGRAPH_ERROR("Invalid number of vertices", IGRAPH_EINVAL);
  }
  if (outseq != 0 && igraph_vector_size(outseq) != 0 && igraph_vector_size(outseq) != n) {
    IGRAPH_ERROR("Invalid out degree sequence length", IGRAPH_EINVAL);
  }
  if ( (outseq == 0 || igraph_vector_size(outseq) == 0) && m<0) {
    IGRAPH_ERROR("Invalid out degree", IGRAPH_EINVAL);
  }

  if (outseq==0 || igraph_vector_size(outseq) == 0) {
    no_of_neighbors=m;
    no_of_edges=(no_of_nodes-1)*no_of_neighbors;
  } else {
    no_of_edges=0;
    for (i=1; i<igraph_vector_size(outseq); i++) {
      no_of_edges+=VECTOR(*outseq)[i];
    }
  }
  
  IGRAPH_VECTOR_INIT_FINALLY(&edges, no_of_edges*2);
  IGRAPH_CHECK(igraph_psumtree_init(&sumtree, no_of_nodes));
  IGRAPH_FINALLY(igraph_psumtree_destroy, &sumtree);
  IGRAPH_VECTOR_INIT_FINALLY(&degree, no_of_nodes);
  
  RNG_BEGIN();
  
  /* first node */
  igraph_psumtree_update(&sumtree, 0, zeroappeal);

  /* and the rest */
  for (i=1; i<no_of_nodes; i++) {
    igraph_real_t sum=igraph_psumtree_sum(&sumtree);
    long int to;
    if (outseq != 0 && igraph_vector_size(outseq)!=0) {
      no_of_neighbors=VECTOR(*outseq)[i];
    }
    for (j=0; j<no_of_neighbors; j++) {
      igraph_psumtree_search(&sumtree, &to, RNG_UNIF(0, sum));
      VECTOR(degree)[to]++;
      VECTOR(edges)[edgeptr++] = i;
      VECTOR(edges)[edgeptr++] = to;
    }
    /* update probabilities */
    for (j=0; j<no_of_neighbors; j++) {
      long int n=VECTOR(edges)[edgeptr-2*j-1];
      igraph_psumtree_update(&sumtree, n, 
			     pow(VECTOR(degree)[n], power)+zeroappeal);
    }
    if (outpref) {
      VECTOR(degree)[i] += no_of_neighbors;
      igraph_psumtree_update(&sumtree, i, 
			     pow(VECTOR(degree)[i], power)+zeroappeal);
    } else {
      igraph_psumtree_update(&sumtree, i, zeroappeal);
    }
  }
  
  RNG_END();

  igraph_psumtree_destroy(&sumtree);
  igraph_vector_destroy(&degree);
  IGRAPH_FINALLY_CLEAN(2);

  IGRAPH_CHECK(igraph_create(graph, &edges, 0, directed));
  igraph_vector_destroy(&edges);
  IGRAPH_FINALLY_CLEAN(1);

  return 0;
}

/**
 * \function igraph_recent_degree_game
 * \brief Stochastic graph generator based on the number of adjacent
 * edges a node has gained recently
 * 
 * \param graph Pointer to an uninitialized graph object.
 * \param n The number of vertices in the graph, this is the same as
 *        the number of time steps. 
 * \param power The exponent, the probability that a node gains a
 *        new edge is proportional to the number of edges it has
 *        gained recently (in the last \p window time steps) to \p
 *        power.
 * \param window Integer constant, the size of the time window to use
 *        to count the number of recent edges.
 * \param m Integer constant, the number of edges to add per time
 *        step if the \p outseq parameter is a null pointer or a
 *        zero-length vector.
 * \param outseq The number of edges to add in each time step. This
 *        argument is ignored if it is a null pointer or a zero length
 *        vector, is this case the constant \p m parameter is used. 
 * \param outpref Logical constant, if true the edges originated by a
 *        vertex also count as recent adjacent edges. It is false in
 *        most cases.
 * \param zero_appeal Constant giving the attractiveness of the
 *        vertices which haven't gained any edge recently. 
 * \param directed Logical constant, whether to generate a directed
 *        graph. 
 * \return Error code.
 * 
 * Time complexity: O(|V|*log(|V|)+|E|), |V| is the number of
 * vertices, |E| is the number of edges in the graph.
 *
 */ 

int igraph_recent_degree_game(igraph_t *graph, igraph_integer_t n,
			      igraph_real_t power,
			      igraph_integer_t window,
			      igraph_integer_t m,  
			      const igraph_vector_t *outseq,
			      igraph_bool_t outpref,
			      igraph_real_t zero_appeal,
			      igraph_bool_t directed) {
  long int no_of_nodes=n;
  long int no_of_neighbors=m;
  long int no_of_edges;
  igraph_vector_t edges;
  long int i, j;
  igraph_psumtree_t sumtree;
  long int edgeptr=0;
  igraph_vector_t degree;
  long int time_window=window;
  igraph_dqueue_t history;

  if (n<0) {
    IGRAPH_ERROR("Invalid number of vertices", IGRAPH_EINVAL);
  }
  if (outseq != 0 && igraph_vector_size(outseq) != 0 && igraph_vector_size(outseq) != n) {
    IGRAPH_ERROR("Invalid out degree sequence length", IGRAPH_EINVAL);
  }
  if ( (outseq == 0 || igraph_vector_size(outseq) == 0) && m<0) {
    IGRAPH_ERROR("Invalid out degree", IGRAPH_EINVAL);
  }

  if (outseq==0 || igraph_vector_size(outseq) == 0) {
    no_of_neighbors=m;
    no_of_edges=(no_of_nodes-1)*no_of_neighbors;
  } else {
    no_of_edges=0;
    for (i=1; i<igraph_vector_size(outseq); i++) {
      no_of_edges+=VECTOR(*outseq)[i];
    }
  }
  
  IGRAPH_VECTOR_INIT_FINALLY(&edges, no_of_edges*2);
  IGRAPH_CHECK(igraph_psumtree_init(&sumtree, no_of_nodes));
  IGRAPH_FINALLY(igraph_psumtree_destroy, &sumtree);
  IGRAPH_VECTOR_INIT_FINALLY(&degree, no_of_nodes);
  IGRAPH_CHECK(igraph_dqueue_init(&history, 
				  time_window*(no_of_neighbors+1)+10));
  IGRAPH_FINALLY(&history, igraph_dqueue_destroy);
  
  RNG_BEGIN();
  
  /* first node */
  igraph_psumtree_update(&sumtree, 0, zero_appeal);
  igraph_dqueue_push(&history, -1);

  /* and the rest */
  for (i=1; i<no_of_nodes; i++) {
    igraph_real_t sum;
    long int to;
    if (outseq != 0 && igraph_vector_size(outseq)!=0) {
      no_of_neighbors=VECTOR(*outseq)[i];
    }

    if (i>=time_window) {
      while ((j=igraph_dqueue_pop(&history)) != -1) {
	VECTOR(degree)[j] -= 1;
	igraph_psumtree_update(&sumtree, j, 
			       pow(VECTOR(degree)[j], power)+zero_appeal);
      }
    }
    
    sum=igraph_psumtree_sum(&sumtree);
    for (j=0; j<no_of_neighbors; j++) {
      igraph_psumtree_search(&sumtree, &to, RNG_UNIF(0, sum));
      VECTOR(degree)[to]++;
      VECTOR(edges)[edgeptr++] = i;
      VECTOR(edges)[edgeptr++] = to;
      igraph_dqueue_push(&history, to);
    }
    igraph_dqueue_push(&history, -1);

    /* update probabilities */
    for (j=0; j<no_of_neighbors; j++) {
      long int n=VECTOR(edges)[edgeptr-2*j-1];
      igraph_psumtree_update(&sumtree, n,
			     pow(VECTOR(degree)[n], power)+zero_appeal);
    }
    if (outpref) {
      VECTOR(degree)[i] += no_of_neighbors;
      igraph_psumtree_update(&sumtree, i, 
			     pow(VECTOR(degree)[i], power)+zero_appeal);
    } else {
      igraph_psumtree_update(&sumtree, i, zero_appeal);
    }
  }
  
  RNG_END();

  igraph_dqueue_destroy(&history);
  igraph_psumtree_destroy(&sumtree);
  igraph_vector_destroy(&degree);
  IGRAPH_FINALLY_CLEAN(3);

  IGRAPH_CHECK(igraph_create(graph, &edges, 0, directed));
  igraph_vector_destroy(&edges);
  IGRAPH_FINALLY_CLEAN(1);

  return 0;
}

/**
 * \function igraph_barabasi_aging_game
 * \brief Preferential attachment with aging of vertices
 * 
 * </para><para>
 * In this game, the probability that a node gains a new edge is
 * given by its (in-)degree (k) and age (l). This probability has a
 * degree dependent component multiplied by an age dependent
 * component. The degree dependent part is: \p deg_coef times k to the
 * power of \p pa_exp plus \p zero_deg_appeal; and the age dependent
 * part is \p age_coef times l to the power of \p aging_exp plus \p
 * zero_age_appeal. 
 * 
 * </para><para>
 * The age is based on the number of vertices in the
 * network and the \p aging_bin argument: vertices grew one unit older
 * after each \p aging_bin vertices added to the network.
 * \param graph Pointer to an uninitialized graph object.
 * \param nodes The number of vertices in the graph.
 * \param m The number of edges to add in each time step. If the \p
 *        outseq argument is not a null vector and not a zero-length
 *        vector. 
 * \param outseq The number of edges to add in each time step. If it
 *        is a null pointer or a zero-length vector then it is ignored
 *        and the \p m argument is used instead.
 * \param outpref Logical constant, whether the edges
 *        initiated by a vertex contribute to the probability to gain
 *        a new edge.
 * \param pa_exp The exponent of the preferential attachment, a small
 *        positive number usually, the value 1 yields the classic
 *        linear preferential attachment.
 * \param aging_exp The exponent of the aging, this is a negative
 *        number usually.
 * \param aging_bin Integer constant, the number of vertices to add
 *        before vertices in the network grew one unit older.
 * \param zero_deg_appeal The degree dependent part of the
 *        attractiveness of the zero degree vertices.
 * \param zero_age_appeal The age dependent part of the attractiveness
 *        of the vertices of age zero. This parameter is usually zero.
 * \param deg_coef The coefficient for the degree.
 * \param age_coef The coefficient for the age.
 * \param directed Logical constant, whether to generate a directed
 *        graph. 
 * \return Error code.
 * 
 * Time complexity: O((|V|+|V|/aging_bin)*log(|V|)+|E|). |V| is the number
 * of vertices, |E| the number of edges.
 */

int igraph_barabasi_aging_game(igraph_t *graph, 
			       igraph_integer_t nodes,
			       igraph_integer_t m,
			       const igraph_vector_t *outseq,
			       igraph_bool_t outpref,
			       igraph_real_t pa_exp,
			       igraph_real_t aging_exp,
			       igraph_integer_t aging_bin,
			       igraph_real_t zero_deg_appeal,
			       igraph_real_t zero_age_appeal,
			       igraph_real_t deg_coef,
			       igraph_real_t age_coef,
			       igraph_bool_t directed) {
  long int no_of_nodes=nodes;
  long int no_of_neighbors=m;
  long int binwidth=nodes/aging_bin+1;
  long int no_of_edges;
  igraph_vector_t edges;
  long int i, j, k;
  igraph_psumtree_t sumtree;
  long int edgeptr=0;
  igraph_vector_t degree;

  if (no_of_nodes<0) {
    IGRAPH_ERROR("Invalid number of vertices", IGRAPH_EINVAL);
  }
  if (outseq != 0 && igraph_vector_size(outseq) != 0 && igraph_vector_size(outseq) != no_of_nodes) {
    IGRAPH_ERROR("Invalid out degree sequence length", IGRAPH_EINVAL);
  }
  if ( (outseq == 0 || igraph_vector_size(outseq) == 0) && m<0) {
    IGRAPH_ERROR("Invalid out degree", IGRAPH_EINVAL);
  }
  if (aging_bin <= 0) { 
    IGRAPH_ERROR("Invalid aging bin", IGRAPH_EINVAL);
  }

  if (outseq==0 || igraph_vector_size(outseq) == 0) {
    no_of_neighbors=m;
    no_of_edges=(no_of_nodes-1)*no_of_neighbors;
  } else {
    no_of_edges=0;
    for (i=1; i<igraph_vector_size(outseq); i++) {
      no_of_edges+=VECTOR(*outseq)[i];
    }
  }

  IGRAPH_VECTOR_INIT_FINALLY(&edges, no_of_edges*2);
  IGRAPH_CHECK(igraph_psumtree_init(&sumtree, no_of_nodes));
  IGRAPH_FINALLY(igraph_psumtree_destroy, &sumtree);
  IGRAPH_VECTOR_INIT_FINALLY(&degree, no_of_nodes);
  
  RNG_BEGIN();
  
  /* first node */
  igraph_psumtree_update(&sumtree, 0, zero_deg_appeal*(1+zero_age_appeal));
  
  /* and the rest */
  for (i=1; i<no_of_nodes; i++) {
    igraph_real_t sum;
    long int to;
    if (outseq != 0 && igraph_vector_size(outseq)!=0) {
      no_of_neighbors=VECTOR(*outseq)[i];
    }
    sum=igraph_psumtree_sum(&sumtree);
    for (j=0; j<no_of_neighbors; j++) {
      igraph_psumtree_search(&sumtree, &to, RNG_UNIF(0, sum));
      VECTOR(degree)[to]++;
      VECTOR(edges)[edgeptr++] = i;
      VECTOR(edges)[edgeptr++] = to;
    }
    /* update probabilites */
    for (j=0; j<no_of_neighbors; j++) {
      long int n=VECTOR(edges)[edgeptr-2*j-1];
      long int age=(i-n)/binwidth;
      igraph_psumtree_update(&sumtree, n, 
			     (deg_coef*pow(VECTOR(degree)[n], pa_exp)
			      +zero_deg_appeal)*
			     (age_coef*pow(age+1,aging_exp)+zero_age_appeal));
    }
    if (outpref) {
      VECTOR(degree)[i] += no_of_neighbors;
      igraph_psumtree_update(&sumtree, i, (zero_age_appeal+1)*
			     (deg_coef*pow(VECTOR(degree)[i], pa_exp)
			      +zero_deg_appeal));
    } else { 
      igraph_psumtree_update(&sumtree, i, (1+zero_age_appeal)*zero_deg_appeal);
    }

    /* aging */
    for (k=1; i-binwidth*k+1 >= 1; k++) {
      long int shnode=i-binwidth*k;
      long int deg=VECTOR(degree)[shnode];
      long int age=(i-shnode)/binwidth;
      igraph_real_t old=igraph_psumtree_get(&sumtree, shnode);
      igraph_psumtree_update(&sumtree, shnode,
			     (deg_coef*pow(deg, pa_exp)+zero_deg_appeal) * 
			     (age_coef*pow(age+2, aging_exp)+zero_age_appeal));
    }
  }
  
  RNG_END();
  
  igraph_vector_destroy(&degree);
  igraph_psumtree_destroy(&sumtree);
  IGRAPH_FINALLY_CLEAN(2);

  IGRAPH_CHECK(igraph_create(graph, &edges, 0, directed));
  igraph_vector_destroy(&edges);
  IGRAPH_FINALLY_CLEAN(1);
  
  return 0;
}

/** 
 * \function igraph_recent_degree_aging_game
 * \brief Preferential attachment based on the number of edges gained
 * recently, with aging of vertices
 * 
 * </para><para>
 * This game is very similar to \ref igraph_barabasi_aging_game(),
 * except that instead of the total number of adjacent edges the
 * number of edges gained in the last \p time_window time steps are
 * counted. 
 * 
 * </para><para>The degree dependent part of the attractiveness is
 * given by k to the power of \p pa_exp plus \p zero_appeal; the age
 * dependent part is l to the power to \p aging_exp. 
 * k is the number of edges gained in the last \p time_window time
 * steps, l is the age of the vertex.
 * \param graph Pointer to an uninitialized graph object.
 * \param nodes The number of vertices in the graph.
 * \param m The number of edges to add in each time step. If the \p
 *        outseq argument is not a null vector or a zero-length vector
 *        then it is ignored.
 * \param outseq Vector giving the number of edges to add in each time
 *        step. If it is a null pointer or a zero-length vector then
 *        it is ignored and the \p m argument is used.
 * \param outpref Logical constant, if true the edges initiated by a
 *        vertex are also counted. Normally it is false.
 * \param pa_exp The exponent for the preferential attachment. 
 * \param aging_exp The exponent for the aging, normally it is
 *        negative: old vertices gain edges with less probability.
 * \param aging_bin Integer constant, gives the scale of the aging. 
 *        The age of the vertices is incremented by one after every \p
 *        aging_bin vertex added.
 * \param time_window The time window to use to count the number of 
 *        adjacent edges for the vertices.
 * \param zero_appeal The degree dependent part of the attractiveness
 *        for zero degree vertices.
 * \param directed Logical constant, whether to create a directed
 *        graph. 
 * \return Error code.
 * 
 * Time complexity: O((|V|+|V|/aging_bin)*log(|V|)+|E|). |V| is the number
 * of vertices, |E| the number of edges.
 */

int igraph_recent_degree_aging_game(igraph_t *graph,
				    igraph_integer_t nodes,
				    igraph_integer_t m, 
				    const igraph_vector_t *outseq,
				    igraph_bool_t outpref,
				    igraph_real_t pa_exp,
				    igraph_real_t aging_exp,
				    igraph_integer_t aging_bin,
				    igraph_integer_t time_window,
				    igraph_real_t zero_appeal,
				    igraph_bool_t directed) {
  
  long int no_of_nodes=nodes;
  long int no_of_neighbors=m;
  long int binwidth=nodes/aging_bin+1;
  long int no_of_edges;
  igraph_vector_t edges;
  long int i, j, k;
  igraph_psumtree_t sumtree;
  long int edgeptr=0;
  igraph_vector_t degree;
  igraph_dqueue_t history;
  
  if (no_of_nodes<0) {
    IGRAPH_ERROR("Invalid number of vertices", IGRAPH_EINVAL);
  }
  if (outseq != 0 && igraph_vector_size(outseq) != 0 && igraph_vector_size(outseq) != no_of_nodes) {
    IGRAPH_ERROR("Invalid out degree sequence length", IGRAPH_EINVAL);
  }
  if ( (outseq == 0 || igraph_vector_size(outseq) == 0) && m<0) {
    IGRAPH_ERROR("Invalid out degree", IGRAPH_EINVAL);
  }
  if (aging_bin <= 0) { 
    IGRAPH_ERROR("Invalid aging bin", IGRAPH_EINVAL);
  }

  if (outseq==0 || igraph_vector_size(outseq) == 0) {
    no_of_neighbors=m;
    no_of_edges=(no_of_nodes-1)*no_of_neighbors;
  } else {
    no_of_edges=0;
    for (i=1; i<igraph_vector_size(outseq); i++) {
      no_of_edges+=VECTOR(*outseq)[i];
    }
  }

  IGRAPH_VECTOR_INIT_FINALLY(&edges, no_of_edges*2);
  IGRAPH_CHECK(igraph_psumtree_init(&sumtree, no_of_nodes));
  IGRAPH_FINALLY(igraph_psumtree_destroy, &sumtree);
  IGRAPH_VECTOR_INIT_FINALLY(&degree, no_of_nodes);
  IGRAPH_CHECK(igraph_dqueue_init(&history, 
				  time_window*(no_of_neighbors+1)+10));
  IGRAPH_FINALLY(&history, igraph_dqueue_destroy);
  
  RNG_BEGIN();
  
  /* first node */
  igraph_psumtree_update(&sumtree, 0, zero_appeal);
  igraph_dqueue_push(&history, -1);
  
  /* and the rest */
  for (i=1; i<no_of_nodes; i++) {
    igraph_real_t sum;
    long int to;
    if (outseq != 0 && igraph_vector_size(outseq)!=0) {
      no_of_neighbors=VECTOR(*outseq)[i];
    }

    if (i>=time_window) {
      while ((j=igraph_dqueue_pop(&history)) != -1) {
	long int age=(i-j)/binwidth;
	VECTOR(degree)[j] -= 1;
	igraph_psumtree_update(&sumtree, j, 
			       (pow(VECTOR(degree)[j], pa_exp)+zero_appeal)*
			       pow(age+1, aging_exp));
      }
    }

    sum=igraph_psumtree_sum(&sumtree);
    for (j=0; j<no_of_neighbors; j++) {
      igraph_psumtree_search(&sumtree, &to, RNG_UNIF(0, sum));
      VECTOR(degree)[to]++;
      VECTOR(edges)[edgeptr++] = i;
      VECTOR(edges)[edgeptr++] = to;
      igraph_dqueue_push(&history, to);
    }
    igraph_dqueue_push(&history, -1);
    
    /* update probabilites */
    for (j=0; j<no_of_neighbors; j++) {
      long int n=VECTOR(edges)[edgeptr-2*j-1];
      long int age=(i-n)/binwidth;
      igraph_psumtree_update(&sumtree, n, 
			     (pow(VECTOR(degree)[n], pa_exp)+zero_appeal)*
			     pow(age+1,aging_exp));
    }
    if (outpref) {
      VECTOR(degree)[i] += no_of_neighbors;
      igraph_psumtree_update(&sumtree, i,
			     pow(VECTOR(degree)[i], pa_exp)+zero_appeal);
    } else { 
      igraph_psumtree_update(&sumtree, i, zero_appeal);
    }

    /* aging */
    for (k=1; i-binwidth*k+1 >= 1; k++) {
      long int shnode=i-binwidth*k;
      long int deg=VECTOR(degree)[shnode];
      long int age=(i-shnode)/binwidth;
      igraph_psumtree_update(&sumtree, shnode,
			     (pow(deg, pa_exp)+zero_appeal) *
			     pow(age+2, aging_exp));
    }
  }
  
  RNG_END();
  
  igraph_dqueue_destroy(&history);
  igraph_vector_destroy(&degree);
  igraph_psumtree_destroy(&sumtree);
  IGRAPH_FINALLY_CLEAN(3);

  IGRAPH_CHECK(igraph_create(graph, &edges, 0, directed));
  igraph_vector_destroy(&edges);
  IGRAPH_FINALLY_CLEAN(1);
  
  return 0;
}
				    
				   
