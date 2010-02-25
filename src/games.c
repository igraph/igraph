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
#include "config.h"

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
    bag=igraph_Calloc(no_of_nodes * no_of_neighbors + no_of_nodes +
	       outpref * no_of_nodes * no_of_neighbors,
	       long int);
    no_of_edges=(no_of_nodes-1)*no_of_neighbors;
  } else {
    no_of_edges=0;
    for (i=1; i<igraph_vector_size(outseq); i++) {
      no_of_edges+=VECTOR(*outseq)[i];
    }
    bag=igraph_Calloc(no_of_nodes + no_of_edges + outpref * no_of_edges,
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

  igraph_Free(bag);
  IGRAPH_CHECK(igraph_create(graph, &edges, n, directed));
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
    double maxedges, last;
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

    last=RNG_GEOM(p);
    while (last < maxedges) {
      IGRAPH_CHECK(igraph_vector_push_back(&s, last));
      last += RNG_GEOM(p);
      last += 1;
    }

    RNG_END();

    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
    IGRAPH_CHECK(igraph_vector_reserve(&edges, igraph_vector_size(&s)*2));

    if (directed && loops) {
      for (i=0; i<igraph_vector_size(&s); i++) {
	long int to=floor(VECTOR(s)[i]/no_of_nodes);
	long int from=VECTOR(s)[i]-((igraph_real_t)to)*no_of_nodes;
	igraph_vector_push_back(&edges, from);
	igraph_vector_push_back(&edges, to);
      }
    } else if (directed && !loops) {
      for (i=0; i<igraph_vector_size(&s); i++) {
	long int to=floor(VECTOR(s)[i]/no_of_nodes);
	long int from=VECTOR(s)[i]-((igraph_real_t)to)*no_of_nodes;
	if (from==to) {
	  to=no_of_nodes-1;
	}
	igraph_vector_push_back(&edges, from);
	igraph_vector_push_back(&edges, to);
      }
    } else if (!directed && loops) {
      for (i=0; i<igraph_vector_size(&s); i++) {
	long int to=floor((sqrt(8*VECTOR(s)[i]+1)-1)/2);
	long int from=VECTOR(s)[i]-(((igraph_real_t)to)*(to+1))/2;
	igraph_vector_push_back(&edges, from);
	igraph_vector_push_back(&edges, to);
      }
    } else /* !directed && !loops */ {
      for (i=0; i<igraph_vector_size(&s); i++) {
	long int to=floor((sqrt(8*VECTOR(s)[i]+1)+1)/2);
	long int from=VECTOR(s)[i]-(((igraph_real_t)to)*(to-1))/2;
	igraph_vector_push_back(&edges, from);
	igraph_vector_push_back(&edges, to);
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
    
      long int slen;
      
      IGRAPH_VECTOR_INIT_FINALLY(&s, 0);
      IGRAPH_CHECK(igraph_random_sample(&s, 0, maxedges-1, no_of_edges));
      
      IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
      IGRAPH_CHECK(igraph_vector_reserve(&edges, igraph_vector_size(&s)*2));

      slen=igraph_vector_size(&s);
      if (directed && loops) {
	for (i=0; i<slen; i++) {
	  long int to=floor(VECTOR(s)[i]/no_of_nodes);
	  long int from=VECTOR(s)[i]-((igraph_real_t)to)*no_of_nodes;
	  igraph_vector_push_back(&edges, from);
	  igraph_vector_push_back(&edges, to);
	}
      } else if (directed && !loops) {
	for (i=0; i<slen; i++) {
	  long int from=floor(VECTOR(s)[i]/(no_of_nodes-1));
	  long int to=VECTOR(s)[i]-((igraph_real_t)from)*(no_of_nodes-1);
	  if (from==to) {
	    to=no_of_nodes-1;
	  }
	  igraph_vector_push_back(&edges, from);
	  igraph_vector_push_back(&edges, to);
	}
      } else if (!directed && loops) {
	for (i=0; i<slen; i++) {
	  long int to=floor((sqrt(8*VECTOR(s)[i]+1)-1)/2);
	  long int from=VECTOR(s)[i]-(((igraph_real_t)to)*(to+1))/2;
	  igraph_vector_push_back(&edges, from);
	  igraph_vector_push_back(&edges, to);
	}
      } else /* !directed && !loops */ {
	for (i=0; i<slen; i++) {
	  long int to=floor((sqrt(8*VECTOR(s)[i]+1)+1)/2);
	  long int from=VECTOR(s)[i]-(((igraph_real_t)to)*(to-1))/2;
	  igraph_vector_push_back(&edges, from);
	  igraph_vector_push_back(&edges, to);
	}
      }  

      igraph_vector_destroy(&s);
      IGRAPH_FINALLY_CLEAN(1);
      retval=igraph_create(graph, &edges, n, directed);
      igraph_vector_destroy(&edges);
      IGRAPH_FINALLY_CLEAN(1);
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
  if (directed)
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

  bag1=igraph_Calloc(outsum, long int);
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
    bag2=igraph_Calloc(insum, long int);
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

  igraph_Free(bag1);
  IGRAPH_FINALLY_CLEAN(1);
  if (directed) {
    igraph_Free(bag2);
    IGRAPH_FINALLY_CLEAN(1);
  }

  IGRAPH_CHECK(igraph_create(graph, &edges, no_of_nodes, directed));
  igraph_vector_destroy(&edges);
  IGRAPH_FINALLY_CLEAN(1);
  
  return 0;
}

/* This is in gengraph_mr-connected.cpp */

int igraph_degree_sequence_game_vl(igraph_t *graph,
				   const igraph_vector_t *out_seq,
				   const igraph_vector_t *in_seq);
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
 *        \c IGRAPH_DEGSEQ_VL is a much more sophisticated generator, 
 *        that can sample undirected, connected simple graphs uniformly. 
 *        It uses Monte-Carlo methods to randomize the graphs. 
 *        This generator should be favoured if undirected and connected 
 *        graphs are to be generated. igraph uses the original implementation 
 *        Fabien Viger; see http://www-rp.lip6.fr/~latapy/FV/generation.html
 *        and the paper cited on it for the details of the algorithm.
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
  } else if (method==IGRAPH_DEGSEQ_VL) {
    retval=igraph_degree_sequence_game_vl(graph, out_deg, in_deg);
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
 * \brief Simulate a growing network with vertex types.
 * 
 * </para><para>
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
  igraph_vector_t nodetypes;

  /* TODO: parameter checks */

  IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&cumdist, types+1);
  IGRAPH_VECTOR_INIT_FINALLY(&nodetypes, nodes);
  
  VECTOR(cumdist)[0]=0;
  for (i=0; i<types; i++) {
    VECTOR(cumdist)[i+1] = VECTOR(cumdist)[i]+VECTOR(*type_dist)[i];
  }
  maxcum=igraph_vector_tail(&cumdist);

  RNG_BEGIN();

  for (i=0; i<nodes; i++) {
    igraph_real_t uni=RNG_UNIF(0, maxcum);
    long int type;
    igraph_vector_binsearch(&cumdist, uni, &type);
    VECTOR(nodetypes)[i]=type-1;
  }    

  for (i=1; i<nodes; i++) {
    for (j=0; j<edges_per_step; j++) {
      long int node1=RNG_INTEGER(0, i);
      long int node2=RNG_INTEGER(0, i);
      long int type1=VECTOR(nodetypes)[node1];
      long int type2=VECTOR(nodetypes)[node2];
/*    printf("unif: %f, %f, types: %li, %li\n", uni1, uni2, type1, type2); */
      if (RNG_UNIF01() < MATRIX(*pref_matrix, type1, type2)) {
	IGRAPH_CHECK(igraph_vector_push_back(&edges, node1));
	IGRAPH_CHECK(igraph_vector_push_back(&edges, node2));
      }
    }
  }

  RNG_END();

  igraph_vector_destroy(&nodetypes);
  igraph_vector_destroy(&cumdist);
  IGRAPH_FINALLY_CLEAN(2);
  IGRAPH_CHECK(igraph_create(graph, &edges, nodes, directed));
  igraph_vector_destroy(&edges);
  IGRAPH_FINALLY_CLEAN(1);
  return 0;
}

/**
 * \function igraph_establishment_game
 * \brief Generates a graph with a simple growing model with vertex types.
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
  igraph_vector_t nodetypes;
  
  IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&cumdist, types+1);
  IGRAPH_VECTOR_INIT_FINALLY(&potneis, k);
  IGRAPH_VECTOR_INIT_FINALLY(&nodetypes, nodes);
  
  VECTOR(cumdist)[0]=0;
  for (i=0; i<types; i++) {
    VECTOR(cumdist)[i+1] = VECTOR(cumdist)[i]+VECTOR(*type_dist)[i];
  }
  maxcum=igraph_vector_tail(&cumdist);

  RNG_BEGIN();

  for (i=0; i<nodes; i++) {
    igraph_real_t uni=RNG_UNIF(0, maxcum);
    long int type;
    igraph_vector_binsearch(&cumdist, uni, &type);
    VECTOR(nodetypes)[i]=type-1;
  }

  for (i=k; i<nodes; i++) {    
    long int type1=VECTOR(nodetypes)[i];
    igraph_random_sample(&potneis, 0, i-1, k);
    for (j=0; j<k; j++) {
      long int type2=VECTOR(nodetypes)[ (long int)VECTOR(potneis)[j] ];
      if (RNG_UNIF01() < MATRIX(*pref_matrix, type1, type2)) {
	IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
	IGRAPH_CHECK(igraph_vector_push_back(&edges, VECTOR(potneis)[j]));
      }
    }
  }
  
  RNG_END();
  
  igraph_vector_destroy(&nodetypes);
  igraph_vector_destroy(&potneis);
  igraph_vector_destroy(&cumdist);
  IGRAPH_FINALLY_CLEAN(3);
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

  IGRAPH_CHECK(igraph_create(graph, &edges, n, directed));
  igraph_vector_destroy(&edges);
  IGRAPH_FINALLY_CLEAN(1);

  return 0;
}

/**
 * \function igraph_recent_degree_game
 * \brief Stochastic graph generator based on the number of adjacent edges a node has gained recently
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
  IGRAPH_FINALLY(igraph_dqueue_destroy, &history);
  
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

  IGRAPH_CHECK(igraph_create(graph, &edges, n, directed));
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
      /* igraph_real_t old=igraph_psumtree_get(&sumtree, shnode); */
      igraph_psumtree_update(&sumtree, shnode,
			     (deg_coef*pow(deg, pa_exp)+zero_deg_appeal) * 
			     (age_coef*pow(age+2, aging_exp)+zero_age_appeal));
    }
  }
  
  RNG_END();
  
  igraph_vector_destroy(&degree);
  igraph_psumtree_destroy(&sumtree);
  IGRAPH_FINALLY_CLEAN(2);

  IGRAPH_CHECK(igraph_create(graph, &edges, nodes, directed));
  igraph_vector_destroy(&edges);
  IGRAPH_FINALLY_CLEAN(1);
  
  return 0;
}

/** 
 * \function igraph_recent_degree_aging_game
 * \brief Preferential attachment based on the number of edges gained recently, with aging of vertices
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
  IGRAPH_FINALLY(igraph_dqueue_destroy, &history);
  
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

  IGRAPH_CHECK(igraph_create(graph, &edges, nodes, directed));
  igraph_vector_destroy(&edges);
  IGRAPH_FINALLY_CLEAN(1);
  
  return 0;
}

/**
 * \function igraph_grg_game
 * \brief Generating geometric random graphs.
 *
 * A geometric random graph is created by dropping points (=vertices)
 * randomly to the unit square and then connecting all those pairs
 * which are less than \c radius apart in Euclidean norm.
 * 
 * </para><para>
 * Original code contributed by Keith Briggs, thanks Keith.
 * \param graph Pointer to an uninitialized graph object,
 * \param nodes The number of vertices in the graph.
 * \param radius The radius within which the vertices will be connected.
 * \param torus Logical constant, if true periodic boundary conditions
 *        will be used, ie. the vertices are assumed to be on a torus 
 *        instead of a square.
 * \return Error code.
 * 
 * Time complexity: TODO, less than O(|V|^2+|E|).
 */
				    
int igraph_grg_game(igraph_t *graph, igraph_integer_t nodes,
		    igraph_real_t radius, igraph_bool_t torus,
		    igraph_vector_t *x, igraph_vector_t *y) {
  
  long int i;
  igraph_vector_t myx, myy, *xx=&myx, *yy=&myy, edges;
  igraph_real_t r2=radius*radius;
  
  IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
  IGRAPH_CHECK(igraph_vector_reserve(&edges, nodes));

  if (x) { 
    xx=x;
    IGRAPH_CHECK(igraph_vector_resize(xx, nodes));
  } else {
    IGRAPH_VECTOR_INIT_FINALLY(xx, nodes);
  }
  if (y) {
    yy=y;
    IGRAPH_CHECK(igraph_vector_resize(yy, nodes));
  } else {
    IGRAPH_VECTOR_INIT_FINALLY(yy, nodes);
  }
  
  RNG_BEGIN();
  
  for (i=0; i<nodes; i++) {
    VECTOR(*xx)[i]=RNG_UNIF01();
    VECTOR(*yy)[i]=RNG_UNIF01();
  }
  
  RNG_END();

  igraph_vector_sort(xx);

  if (!torus) {
    for (i=0; i<nodes; i++) {
      igraph_real_t x1=VECTOR(*xx)[i];
      igraph_real_t y1=VECTOR(*yy)[i];
      long int j=i+1;
      igraph_real_t dx, dy;
      while ( j<nodes && (dx=VECTOR(*xx)[j] - x1) < radius) {
	dy=VECTOR(*yy)[j]-y1;
	if (dx*dx+dy*dy < r2) {
	  IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
	  IGRAPH_CHECK(igraph_vector_push_back(&edges, j));
	}
	j++;
      }
    }
  } else { 
    for (i=0; i<nodes; i++) {
      igraph_real_t x1=VECTOR(*xx)[i];
      igraph_real_t y1=VECTOR(*yy)[i];
      long int j=i+1;
      igraph_real_t dx, dy;
      while ( j<nodes && (dx=VECTOR(*xx)[j] - x1) < radius) {
	dy=fabs(VECTOR(*yy)[j]-y1);
	if (dx > 0.5) {
	  dx=1-dx;
	}
	if (dy > 0.5) {
	  dy=1-dy;
	}
	if (dx*dx+dy*dy < r2) {
	  IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
	  IGRAPH_CHECK(igraph_vector_push_back(&edges, j));
	}
	j++;
      }
      if (j==nodes) {
	j=0;      
	while (j<i && (dx=1-x1+VECTOR(*xx)[j]) < radius && 
	       x1-VECTOR(*xx)[j]>=radius) {
	  dy=fabs(VECTOR(*yy)[j]-y1);
	  if (dy > 0.5) {
	    dy=1-dy;
	  }
	  if (dx*dx+dy*dy < r2) {
	    IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
	    IGRAPH_CHECK(igraph_vector_push_back(&edges, j));
	  }
	  j++;
	}
      }
    }
  }
  
  if (!y) {
    igraph_vector_destroy(yy);
    IGRAPH_FINALLY_CLEAN(1);
  }
  if (!x) {
    igraph_vector_destroy(xx);
    IGRAPH_FINALLY_CLEAN(1);
  }

  IGRAPH_CHECK(igraph_create(graph, &edges, nodes, IGRAPH_UNDIRECTED));
  igraph_vector_destroy(&edges);
  IGRAPH_FINALLY_CLEAN(1);
  
  return 0;
}


void igraph_i_preference_game_free_vids_by_type(igraph_vector_ptr_t *vecs) {
  int i=0, n;
  igraph_vector_t *v;

  n = igraph_vector_ptr_size(vecs);
  for (i=0; i<n; i++) {
    v = (igraph_vector_t*)VECTOR(*vecs)[i];
	if (v) igraph_vector_destroy(v);
  }
  igraph_vector_ptr_destroy_all(vecs);
}

/**
 * \function igraph_preference_game
 * \brief Generates a graph with vertex types and connection preferences 
 * 
 * </para><para>
 * This is practically the nongrowing variant of \ref igraph_establishment_game .
 * A given number of vertices are generated. Every vertex is assigned to a
 * vertex type according to the given type probabilities. Finally, every
 * vertex pair is evaluated and an edge is created between them with a
 * probability depending on the types of the vertices involved.
 * 
 * \param graph Pointer to an uninitialized graph.
 * \param nodes The number of vertices in the graph.
 * \param types The number of vertex types.
 * \param type_dist Vector giving the distribution of vertex types. If \c NULL ,
 *   all vertex types will have equal probability.
 * \param pref_matrix Matrix giving the connection probabilities for
 *   different vertex types. This should be symmetric if the requested
 *   graph is undirected.
 * \param node_type_vec A vector where the individual generated vertex types
 *   will be stored. If \c NULL , the vertex types won't be saved.
 * \param directed Logical, whether to generate a directed graph. If undirected
 *   graphs are requested, only the lower left triangle of the preference
 *   matrix is considered.
 * \param loops Logical, whether loop edges are allowed.
 * \return Error code.
 *
 * Added in version 0.3.</para><para>
 *
 * Time complexity: O(|V|+|E|), the
 * number of vertices plus the number of edges in the graph.
 * 
 * \sa igraph_establishment_game()
 */

int igraph_preference_game(igraph_t *graph, igraph_integer_t nodes,
			   igraph_integer_t types,
			   igraph_vector_t *type_dist,
			   igraph_matrix_t *pref_matrix,
			   igraph_vector_t *node_type_vec,
			   igraph_bool_t directed,
			   igraph_bool_t loops) {
  
  long int i, j;
  igraph_vector_t edges, s;
  igraph_vector_t cumdist;
  igraph_vector_t* nodetypes;
  igraph_vector_ptr_t vids_by_type;
  igraph_real_t maxcum, maxedges;

  if (types < 1) IGRAPH_ERROR("types must be >= 1", IGRAPH_EINVAL);
  if (nodes < 0) IGRAPH_ERROR("nodes must be >= 0", IGRAPH_EINVAL);
  if (type_dist && igraph_vector_size(type_dist) != types) {
    if (igraph_vector_size(type_dist) > types)
      IGRAPH_WARNING("length of type_dist > types, type_dist will be trimmed");
    else
      IGRAPH_ERROR("type_dist vector too short", IGRAPH_EINVAL);
  }
  if (igraph_matrix_nrow(pref_matrix) < types ||
      igraph_matrix_ncol(pref_matrix) < types)
    IGRAPH_ERROR("pref_matrix too small", IGRAPH_EINVAL);

  IGRAPH_VECTOR_INIT_FINALLY(&cumdist, types+1);

  if (node_type_vec) {
    IGRAPH_CHECK(igraph_vector_resize(node_type_vec, nodes));
    nodetypes = node_type_vec;
  } else {
    nodetypes = igraph_Calloc(1, igraph_vector_t);
    if (nodetypes == 0) {
      IGRAPH_ERROR("preference_game failed", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, nodetypes);
    IGRAPH_VECTOR_INIT_FINALLY(nodetypes, nodes);
  }

  IGRAPH_CHECK(igraph_vector_ptr_init(&vids_by_type, types));
  IGRAPH_FINALLY(igraph_vector_ptr_destroy_all, &vids_by_type);
  for (i=0; i<types; i++) {
	  VECTOR(vids_by_type)[i] = igraph_Calloc(1, igraph_vector_t);
	  if (VECTOR(vids_by_type)[i] == 0) {
	    IGRAPH_ERROR("preference_game failed", IGRAPH_ENOMEM);
	  }
	  IGRAPH_CHECK(igraph_vector_init(VECTOR(vids_by_type)[i], 0));
  }
  IGRAPH_FINALLY_CLEAN(1);   /* removing igraph_vector_ptr_destroy_all */
  IGRAPH_FINALLY(igraph_i_preference_game_free_vids_by_type, &vids_by_type);

  VECTOR(cumdist)[0]=0;
  if (type_dist) {
    for (i=0; i<types; i++) 
      VECTOR(cumdist)[i+1] = VECTOR(cumdist)[i]+VECTOR(*type_dist)[i];
  } else {
    for (i=0; i<types; i++) VECTOR(cumdist)[i+1] = i+1;
  }
  maxcum=igraph_vector_tail(&cumdist);

  RNG_BEGIN();

  for (i=0; i<nodes; i++) {
    long int type1;
    igraph_real_t uni1=RNG_UNIF(0, maxcum);
    igraph_vector_binsearch(&cumdist, uni1, &type1);
    VECTOR(*nodetypes)[i] = type1-1;
    IGRAPH_CHECK(igraph_vector_push_back(
	    (igraph_vector_t*)VECTOR(vids_by_type)[type1-1], i));
  }

  igraph_vector_destroy(&cumdist);
  IGRAPH_FINALLY_CLEAN(1);

  IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&s, 0);

  for (i=0; i<types; i++) {
    for (j=0; j<types; j++) {
      /* Generating the random subgraph between vertices of type i and j */
      long int k, l;
      igraph_real_t p, last;
      igraph_vector_t *v1, *v2;
      long int v1_size, v2_size;
      
      IGRAPH_ALLOW_INTERRUPTION();

      v1 = (igraph_vector_t*)VECTOR(vids_by_type)[i];
      v2 = (igraph_vector_t*)VECTOR(vids_by_type)[j];
      v1_size = igraph_vector_size(v1);
      v2_size = igraph_vector_size(v2);

      p = MATRIX(*pref_matrix, i, j);
      igraph_vector_clear(&s);
      if (i != j) {
        /* The two vertex sets are disjoint, this is the easier case */
        if (i > j && !directed) continue;
        maxedges = v1_size * v2_size;
      } else {
        if (directed && loops) maxedges = v1_size * v1_size;
        else if (directed && !loops) maxedges = v1_size * (v1_size-1);
        else if (!directed && loops) maxedges = v1_size * (v1_size+1)/2;
        else maxedges = v1_size * (v1_size-1)/2;
      }

      IGRAPH_CHECK(igraph_vector_reserve(&s, maxedges*p*1.1));
        
      last=RNG_GEOM(p);
      while (last < maxedges) {
        IGRAPH_CHECK(igraph_vector_push_back(&s, last));
        last += RNG_GEOM(p);
        last += 1;
      }
      l = igraph_vector_size(&s);
        
      IGRAPH_CHECK(igraph_vector_reserve(&edges, igraph_vector_size(&edges)+l*2));
      
      if (i != j) {
        /* Generating the subgraph between vertices of type i and j */
        for (k=0; k<l; k++) {
          long int to=floor(VECTOR(s)[k]/v1_size);
          long int from=VECTOR(s)[k]-((igraph_real_t)to)*v1_size;
          igraph_vector_push_back(&edges, VECTOR(*v1)[from]);
          igraph_vector_push_back(&edges, VECTOR(*v2)[to]);
        }
      } else {
        /* Generating the subgraph among vertices of type i */
        if (directed && loops) {
          for (k=0; k<l; k++) {
            long int to=floor(VECTOR(s)[k]/v1_size);
            long int from=VECTOR(s)[k]-((igraph_real_t)to)*v1_size;
            igraph_vector_push_back(&edges, VECTOR(*v1)[from]);
            igraph_vector_push_back(&edges, VECTOR(*v1)[to]);
          }
        } else if (directed && !loops) {
          for (k=0; k<l; k++) {
            long int to=floor(VECTOR(s)[k]/v1_size);
            long int from=VECTOR(s)[k]-((igraph_real_t)to)*v1_size;
            if (from==to) to=v1_size-1;
            igraph_vector_push_back(&edges, VECTOR(*v1)[from]); 
            igraph_vector_push_back(&edges, VECTOR(*v1)[to]);
          }
        } else if (!directed && loops) {
          for (k=0; k<l; k++) {
            long int to=floor((sqrt(8*VECTOR(s)[k]+1)-1)/2);
            long int from=VECTOR(s)[k]-(((igraph_real_t)to)*(to+1))/2;
            igraph_vector_push_back(&edges, VECTOR(*v1)[from]);
            igraph_vector_push_back(&edges, VECTOR(*v1)[to]);
          }
        } else {
          for (k=0; k<l; k++) {
            long int to=floor((sqrt(8*VECTOR(s)[k]+1)+1)/2);
            long int from=VECTOR(s)[k]-(((igraph_real_t)to)*(to-1))/2;
            igraph_vector_push_back(&edges, VECTOR(*v1)[from]);
            igraph_vector_push_back(&edges, VECTOR(*v1)[to]);
          }
        }
      }
    }
  }

  RNG_END();

  igraph_vector_destroy(&s);
  igraph_i_preference_game_free_vids_by_type(&vids_by_type);
  IGRAPH_FINALLY_CLEAN(2);

  if (node_type_vec == 0) {
    igraph_vector_destroy(nodetypes);
    igraph_Free(nodetypes);
    IGRAPH_FINALLY_CLEAN(2);
  }

  IGRAPH_CHECK(igraph_create(graph, &edges, nodes, directed));
  igraph_vector_destroy(&edges);
  IGRAPH_FINALLY_CLEAN(1);
  return 0;
}

/**
 * \function igraph_asymmetric_preference_game
 * \brief Generates a graph with asymmetric vertex types and connection preferences 
 * 
 * </para><para>
 * This is the asymmetric variant of \ref igraph_preference_game() .
 * A given number of vertices are generated. Every vertex is assigned to an
 * "incoming" and an "outgoing" vertex type according to the given joint
 * type probabilities. Finally, every vertex pair is evaluated and a
 * directed edge is created between them with a probability depending on the
 * "outgoing" type of the source vertex and the "incoming" type of the target
 * vertex.
 * 
 * \param graph Pointer to an uninitialized graph.
 * \param nodes The number of vertices in the graph.
 * \param types The number of vertex types.
 * \param type_dist_matrix Matrix giving the joint distribution of vertex types.
 *   If null, incoming and outgoing vertex types are independent and uniformly
 *   distributed.
 * \param pref_matrix Matrix giving the connection probabilities for
 *   different vertex types.
 * \param node_type_in_vec A vector where the individual generated "incoming"
 *   vertex types will be stored. If NULL, the vertex types won't be saved.
 * \param node_type_out_vec A vector where the individual generated "outgoing"
 *   vertex types will be stored. If NULL, the vertex types won't be saved.
 * \param loops Logical, whether loop edges are allowed.
 * \return Error code.
 *
 * Added in version 0.3.</para><para>
 *
 * Time complexity: O(|V|+|E|), the
 * number of vertices plus the number of edges in the graph.
 * 
 * \sa \ref igraph_preference_game()
 */

int igraph_asymmetric_preference_game(igraph_t *graph, igraph_integer_t nodes,
				      igraph_integer_t types,
				      igraph_matrix_t *type_dist_matrix,
				      igraph_matrix_t *pref_matrix,
				      igraph_vector_t *node_type_in_vec,
				      igraph_vector_t *node_type_out_vec,
				      igraph_bool_t loops) {
  
  long int i, j, k;
  igraph_vector_t edges, cumdist, s, intersect;
  igraph_vector_t *nodetypes_in;
  igraph_vector_t *nodetypes_out;
  igraph_vector_ptr_t vids_by_intype, vids_by_outtype;
  igraph_real_t maxcum, maxedges;
  
  if (types < 1) IGRAPH_ERROR("types must be >= 1", IGRAPH_EINVAL);
  if (nodes < 0) IGRAPH_ERROR("nodes must be >= 0", IGRAPH_EINVAL);
  if (type_dist_matrix) {
    if (igraph_matrix_nrow(type_dist_matrix) < types ||
        igraph_matrix_ncol(type_dist_matrix) < types)
      IGRAPH_ERROR("type_dist_matrix too small", IGRAPH_EINVAL);
    else if (igraph_matrix_nrow(type_dist_matrix) > types ||
        igraph_matrix_ncol(type_dist_matrix) > types)
      IGRAPH_WARNING("type_dist_matrix will be trimmed");
  }
  if (igraph_matrix_nrow(pref_matrix) < types ||
      igraph_matrix_ncol(pref_matrix) < types)
    IGRAPH_ERROR("pref_matrix too small", IGRAPH_EINVAL);

  IGRAPH_VECTOR_INIT_FINALLY(&cumdist, types*types+1);

  if (node_type_in_vec) {
    nodetypes_in=node_type_in_vec;
    IGRAPH_CHECK(igraph_vector_resize(nodetypes_in, nodes));
  } else {
    nodetypes_in = igraph_Calloc(1, igraph_vector_t);
    if (nodetypes_in == 0) {
      IGRAPH_ERROR("asymmetric_preference_game failed", IGRAPH_ENOMEM);
    }
    IGRAPH_VECTOR_INIT_FINALLY(nodetypes_in, nodes);
  }
   
  if (node_type_out_vec) {
    nodetypes_out=node_type_out_vec;
    IGRAPH_CHECK(igraph_vector_resize(nodetypes_out, nodes));
  } else {
    nodetypes_out = igraph_Calloc(1, igraph_vector_t);
    if (nodetypes_out == 0) {
      IGRAPH_ERROR("asymmetric_preference_game failed", IGRAPH_ENOMEM);
    }
    IGRAPH_VECTOR_INIT_FINALLY(nodetypes_out, nodes);
  }
   
  IGRAPH_CHECK(igraph_vector_ptr_init(&vids_by_intype, types));
  IGRAPH_FINALLY(igraph_vector_ptr_destroy_all, &vids_by_intype);
  IGRAPH_CHECK(igraph_vector_ptr_init(&vids_by_outtype, types));
  IGRAPH_FINALLY(igraph_vector_ptr_destroy_all, &vids_by_outtype);
  for (i=0; i<types; i++) {
	  VECTOR(vids_by_intype)[i] = igraph_Calloc(1, igraph_vector_t);
	  VECTOR(vids_by_outtype)[i] = igraph_Calloc(1, igraph_vector_t);
	  if (VECTOR(vids_by_intype)[i] == 0 || VECTOR(vids_by_outtype)[i] == 0) {
	    IGRAPH_ERROR("asymmetric_preference_game failed", IGRAPH_ENOMEM);
	  }
	  IGRAPH_CHECK(igraph_vector_init(VECTOR(vids_by_intype)[i], 0));
	  IGRAPH_CHECK(igraph_vector_init(VECTOR(vids_by_outtype)[i], 0));
  }
  IGRAPH_FINALLY_CLEAN(2);   /* removing igraph_vector_ptr_destroy_all */
  IGRAPH_FINALLY(igraph_i_preference_game_free_vids_by_type, &vids_by_intype);
  IGRAPH_FINALLY(igraph_i_preference_game_free_vids_by_type, &vids_by_outtype);

  VECTOR(cumdist)[0]=0;
  if (type_dist_matrix) {
    for (i=0, k=0; i<types; i++) {
      for (j=0; j<types; j++, k++) {
        VECTOR(cumdist)[k+1] = VECTOR(cumdist)[k]+MATRIX(*type_dist_matrix, i, j);
      }
    }
  } else {
    for (i=0; i<types*types; i++) VECTOR(cumdist)[i+1] = i+1;
  }
  maxcum=igraph_vector_tail(&cumdist);

  RNG_BEGIN();

  for (i=0; i<nodes; i++) {
    long int type1, type2;
    igraph_real_t uni1=RNG_UNIF(0, maxcum);
    igraph_vector_binsearch(&cumdist, uni1, &type1);
    type2 = (type1-1) % (int)types;
    type1 = (type1-1) / (int)types;
    VECTOR(*nodetypes_in)[i] = type1; 
    VECTOR(*nodetypes_out)[i] = type2; 
    IGRAPH_CHECK(igraph_vector_push_back(
	    (igraph_vector_t*)VECTOR(vids_by_intype)[type1], i));
    IGRAPH_CHECK(igraph_vector_push_back(
	    (igraph_vector_t*)VECTOR(vids_by_outtype)[type2], i));
  }

  igraph_vector_destroy(&cumdist);
  IGRAPH_FINALLY_CLEAN(1);

  IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&s, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&intersect, 0);
  for (i=0; i<types; i++) {
    for (j=0; j<types; j++) {
      long int k, l, c;
      igraph_real_t p, last;
      igraph_vector_t *v1, *v2;
      long int v1_size, v2_size;

      IGRAPH_ALLOW_INTERRUPTION();
      
      v1 = (igraph_vector_t*)VECTOR(vids_by_outtype)[i];
      v2 = (igraph_vector_t*)VECTOR(vids_by_intype)[j];
      v1_size = igraph_vector_size(v1);
      v2_size = igraph_vector_size(v2);

      maxedges = v1_size * v2_size;
      if (!loops) {
        IGRAPH_CHECK(igraph_vector_intersect_sorted(v1, v2, &intersect, 0));
        c = igraph_vector_size(&intersect);
        maxedges -= c; 
      }

      p = MATRIX(*pref_matrix, i, j);
      igraph_vector_clear(&s);
      IGRAPH_CHECK(igraph_vector_reserve(&s, maxedges*p*1.1));
        
      last=RNG_GEOM(p);
      while (last < maxedges) {
        IGRAPH_CHECK(igraph_vector_push_back(&s, last));
        last += RNG_GEOM(p);
        last += 1;
      }
      l = igraph_vector_size(&s);
        
      IGRAPH_CHECK(igraph_vector_reserve(&edges, igraph_vector_size(&edges)+l*2));

      if (!loops && c>0) {
          for (k=0; k<l; k++) {
            long int to=floor(VECTOR(s)[k]/v1_size);
            long int from=VECTOR(s)[k]-((igraph_real_t)to)*v1_size;
            if (VECTOR(*v1)[from] == VECTOR(*v2)[to]) {
              /* remap loop edges */
              to = v2_size-1;
              igraph_vector_binsearch(&intersect, VECTOR(*v1)[from], &c);
              from = v1_size-1;
              if (VECTOR(*v1)[from] == VECTOR(*v2)[to]) from--;
              while (c>0) {
                c--; from--;
                if (VECTOR(*v1)[from] == VECTOR(*v2)[to]) from--;
              }
            }
            igraph_vector_push_back(&edges, VECTOR(*v1)[from]); 
            igraph_vector_push_back(&edges, VECTOR(*v2)[to]);
          }
      } else {
        for (k=0; k<l; k++) {
          long int to=floor(VECTOR(s)[k]/v1_size);
          long int from=VECTOR(s)[k]-((igraph_real_t)to)*v1_size;
          igraph_vector_push_back(&edges, VECTOR(*v1)[from]);
          igraph_vector_push_back(&edges, VECTOR(*v2)[to]);
        }
      }
    }
  }

  RNG_END();
  
  igraph_vector_destroy(&s);
  igraph_vector_destroy(&intersect);
  igraph_i_preference_game_free_vids_by_type(&vids_by_intype);
  igraph_i_preference_game_free_vids_by_type(&vids_by_outtype);
  IGRAPH_FINALLY_CLEAN(4);

  if (node_type_out_vec == 0) {
    igraph_vector_destroy(nodetypes_out);
    igraph_Free(nodetypes_out);
    IGRAPH_FINALLY_CLEAN(1);
  }
   
  if (node_type_in_vec == 0) {
    igraph_vector_destroy(nodetypes_in);
    igraph_Free(nodetypes_in);
    IGRAPH_FINALLY_CLEAN(1);
  }
  
  IGRAPH_CHECK(igraph_create(graph, &edges, nodes, 1));
  igraph_vector_destroy(&edges);
  IGRAPH_FINALLY_CLEAN(1);
  return 0;
}

/**
 * \function igraph_rewire_edges
 * \brief Rewire the edges of a graph with constant probability
 * 
 * This function rewires the edges of a graph with a constant
 * probability. More precisely each end point of each edge is rewired
 * to an uniformly randomly chosen vertex with constant probability \p
 * prob.
 * 
 * </para><para> Note that this function modifies the input \p graph,
 * call \ref igraph_copy() if you want to keep it.
 * 
 * \param graph The input graph, this will be rewired, it can be
 *    directed or undirected.
 * \param prob The rewiring probability a constant between zero and
 *    one (inclusive).
 * \return Error code.
 * 
 * \sa \ref igraph_watts_strogatz_game() uses this function for the
 * rewiring.
 * 
 * Time complexity: O(|V|+|E|).
 */

int igraph_rewire_edges(igraph_t *graph, igraph_real_t prob) {

  igraph_t newgraph;
  igraph_vector_t edges;
  long int no_of_edges=igraph_ecount(graph);
  long int no_of_nodes=igraph_vcount(graph);
  long int endpoints=no_of_edges*2;
  long int to_rewire;

  if (prob < 0 || prob > 1) {
    IGRAPH_ERROR("Rewiring probability should be between zero and one",
		 IGRAPH_EINVAL);
  }
  
  IGRAPH_VECTOR_INIT_FINALLY(&edges, endpoints);
  IGRAPH_CHECK(igraph_get_edgelist(graph, &edges, 0));
  
  /* Now do the rewiring, fast method. 
     Each endpoint of an edge is rewired with probability p.
     So the "skips" between the really rewired endpoints follow a 
     geometric distribution.
  */

  RNG_BEGIN();

  if (prob != 0) {
    to_rewire=RNG_GEOM(prob)+1;
    while (to_rewire <= endpoints) {
      VECTOR(edges)[ to_rewire-1 ] = RNG_INTEGER(0, no_of_nodes-1);
      to_rewire += RNG_GEOM(prob)+1;
    }
  }
  
  RNG_END();
  
  IGRAPH_CHECK(igraph_create(&newgraph, &edges, no_of_nodes, 
			     igraph_is_directed(graph)));
  igraph_vector_destroy(&edges);
  IGRAPH_FINALLY_CLEAN(1); 
  
  IGRAPH_FINALLY(igraph_destroy, &newgraph);
  IGRAPH_I_ATTRIBUTE_DESTROY(&newgraph);
  IGRAPH_I_ATTRIBUTE_COPY(&newgraph, graph, 1,1,1);
  IGRAPH_FINALLY_CLEAN(1);
  igraph_destroy(graph);
  *graph=newgraph;
  
  return 0;
}

/**
 * \function igraph_watts_strogatz_game
 * \brief The Watts-Strogatz small-world model
 * 
 * This function generates a graph according to the Watts-Strogatz
 * model of small-world networks. The graph is obtained by creating a
 * circular undirected lattice and then rewire the edges randomly with
 * a constant probability.
 * 
 * </para><para>See also: Duncan J Watts and Steven H Strogatz:
 * Collective dynamics of <quote>small world</quote> networks, Nature
 * 393, 440-442, 1998.
 * \param graph The graph to initialize.
 * \param dim The dimension of the lattice.
 * \param size The size of the lattice along each dimension.
 * \param nei The size of the neighborhood for each vertex. This is
 *    the same as the \p nei argument of \ref
 *    igraph_connect_neighborhood(). 
 * \param p The rewiring probability. A real number between zero and
 *   one (inclusive). 
 * \return Error code.
 * 
 * \sa \ref igraph_lattice(), \ref igraph_connect_neighborhood() and
 * \ref igraph_rewire_edges() can be used if more flexibility is
 * needed, eg. a different type of lattice.
 * 
 * Time complexity: O(|V|*d^o+|E|), |V| ans |E| are the number of
 * vertices and edges, d is the average degree, o is the \p nei
 * argument.
 */

int igraph_watts_strogatz_game(igraph_t *graph, igraph_integer_t dim,
			       igraph_integer_t size, igraph_integer_t nei,
			       igraph_real_t p) {
  
  igraph_vector_t dimvector;
  long int i;

  if (dim < 1) {
    IGRAPH_ERROR("WS game: dimension should be at least one", IGRAPH_EINVAL);
  }
  if (size < 1) { 
    IGRAPH_ERROR("WS game: lattice size should be at least one", 
		 IGRAPH_EINVAL);
  }
  if (p < 0 || p > 1) {
    IGRAPH_ERROR("WS game: rewiring probability should be between 0 and 1",
		 IGRAPH_EINVAL);
  }

  /* Create the lattice first */

  IGRAPH_VECTOR_INIT_FINALLY(&dimvector, dim);
  for (i=0; i<dim; i++) {
    VECTOR(dimvector)[i] = size;
  }
  
  IGRAPH_CHECK(igraph_lattice(graph, &dimvector, nei, IGRAPH_UNDIRECTED,
			      0 /* mutual */, 1 /* circular */));
  igraph_vector_destroy(&dimvector);
  IGRAPH_FINALLY_CLEAN(1);
  IGRAPH_FINALLY(igraph_destroy, graph);
  
  /* Rewire the edges then */

  IGRAPH_CHECK(igraph_rewire_edges(graph, p));

  IGRAPH_FINALLY_CLEAN(1);
  return 0;
}

/**
 * \function igraph_lastcit_game
 * \brief Simulate citation network, based on time passed since the last citation.
 * 
 * This is a quite special stochastic graph generator, it models an
 * evolving graph. In each time step a single vertex is added to the 
 * network and it cites a number of other vertices (as specified by 
 * the \p edges_per_step argument). The cited vertices are selected
 * based on the last time they were cited. Time is mesaured by the 
 * addition of vertices and it is binned into \p pagebins bins. 
 * So if the current time step is \c t and the last citation to a 
 * given \c i vertex was made in time step \c t0, then \c
 * (t-t0)/binwidth is calculated where binwidth is \c nodes/pagebins+1,
 * in the last expression '/' denotes integer division, so the
 * fraction part is omitted. 
 * 
 * </para><para>
 * The \p preference argument specifies the preferences for the
 * citation lags, ie. its first elements contains the attractivity 
 * of the very recently cited vertices, etc. The last element is
 * special, it contains the attractivity of the vertices which were
 * never cited. This element should be bigger than zero.
 * 
 * </para><para>
 * Note that this function generates networks with multiple edges if 
 * \p edges_per_step is bigger than one, call \ref igraph_simplify()
 * on the result to get rid of these edges.
 * \param graph Pointer to an uninitialized graph object, the result
 *     will be stored here.
 * \param node The number of vertices in the network.
 * \param edges_per_node The number of edges to add in each time
 *     step. 
 * \param pagebins The number of age bins to use.
 * \param preference Pointer to an initialized vector of length
 *     \c pagebins+1. This contains the `attractivity' of the various
 *     age bins, the last element is the attractivity of the vertices 
 *     which were never cited, and it should be greater than zero.
 *     It is a good idea to have all positive values in this vector.
 * \param directed Logical constant, whether to create directed
 *      networks. 
 * \return Error code.
 * 
 * \sa \ref igraph_barabasi_aging_game().
 * 
 * Time complexity: O(|V|*a+|E|*log|V|), |V| is the number of vertices,
 * |E| is the total number of edges, a is the \p pagebins parameter.
 */
  
int igraph_lastcit_game(igraph_t *graph, 
			igraph_integer_t nodes, igraph_integer_t edges_per_node, 
			igraph_integer_t pagebins,
			const igraph_vector_t *preference,
			igraph_bool_t directed) {

  long int no_of_nodes=nodes;
  igraph_psumtree_t sumtree;
  igraph_vector_t edges;
  long int i, j, k;
  long int *lastcit;
  long int *index;
  long int agebins=pagebins;
  long int binwidth=no_of_nodes/agebins+1;

  if (agebins != igraph_vector_size(preference)-1) {
    IGRAPH_ERROR("`preference' vector should be of length `agebins' plus one",
		 IGRAPH_EINVAL);
  }
  if (agebins <=1 ) {
    IGRAPH_ERROR("at least two age bins are need for lastcit game",
		 IGRAPH_EINVAL);
  }
  if (VECTOR(*preference)[agebins] <= 0) {
    IGRAPH_ERROR("the last element of the `preference' vector needs to be positive",
		 IGRAPH_EINVAL);
  }
  
  IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);

  lastcit=igraph_Calloc(no_of_nodes, long int);
  if (!lastcit) {
    IGRAPH_ERROR("lastcit game failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, lastcit);

  index=igraph_Calloc(no_of_nodes+1, long int);
  if (!index) {
    IGRAPH_ERROR("lastcit game failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, index);

  IGRAPH_CHECK(igraph_psumtree_init(&sumtree, nodes));
  IGRAPH_FINALLY(igraph_psumtree_destroy, &sumtree);
  IGRAPH_CHECK(igraph_vector_reserve(&edges, nodes*edges_per_node));  
  
  /* The first node */
  igraph_psumtree_update(&sumtree, 0, VECTOR(*preference)[agebins]);
  index[0]=0;
  index[1]=0;

  RNG_BEGIN();

  for (i=1; i<no_of_nodes; i++) {

    /* Add new edges */
    for (j=0; j<edges_per_node; j++) {
      long int to;
      igraph_real_t sum=igraph_psumtree_sum(&sumtree);
      igraph_psumtree_search(&sumtree, &to, RNG_UNIF(0, sum));
      igraph_vector_push_back(&edges, i);
      igraph_vector_push_back(&edges, to);
      lastcit[to]=i+1;
      igraph_psumtree_update(&sumtree, to, VECTOR(*preference)[0]);
    }

    /* Add the node itself */
    igraph_psumtree_update(&sumtree, i, VECTOR(*preference)[agebins]);
    index[i+1]=index[i]+edges_per_node;

    /* Update the preference of some vertices if they got to another bin.
       We need to know the citations of some older vertices, this is in the index. */
    for (k=1; i-binwidth*k >= 1; k++) {
      long int shnode=i-binwidth*k;
      long int m=index[shnode], n=index[shnode+1];
      for (j=2*m; j<2*n; j+=2) {
	long int cnode=VECTOR(edges)[j+1];
	if (lastcit[cnode]==shnode+1) {
	  igraph_psumtree_update(&sumtree, cnode, VECTOR(*preference)[k]);
	}
      }
    }
    
  }
  
  RNG_END();
  
  igraph_psumtree_destroy(&sumtree);
  igraph_free(index);
  igraph_free(lastcit);
  IGRAPH_FINALLY_CLEAN(3);
  
  IGRAPH_CHECK(igraph_create(graph, &edges, nodes, directed));
  igraph_vector_destroy(&edges);
  IGRAPH_FINALLY_CLEAN(1);

  return 0;
}

/**
 * \function igraph_cited_type_game
 * \brief Simulate a citation based on vertex types.
 * 
 * Function to create a network based on some vertex categories. This
 * function creates a citation network, in each step a single vertex
 * and \p edges_per_step citating edges are added, nodes with
 * different categories (may) have different probabilities to get
 * cited, as given by the \p pref vector.
 *
 * </para><para>
 * Note that this function might generate networks with multiple edges 
 * if \p edges_per_step is greater than one. You might want to call
 * \ref igraph_simplify() on the result to remove multiple edges.
 * \param graph Pointer to an uninitialized graph object.
 * \param nodes The number of vertices in the network.
 * \param types Numeric vector giving the categories of the vertices,
 *     so it should contain \p nodes non-negative integer
 *     numbers. Types are numbered from zero.
 * \param pref The attractivity of the different vertex categories in
 *     a vector. Its length should be the maximum element in \p types
 *     plus one (types are numbered from zero).
 * \param edges_per_step Integer constant, the number of edges to add
 *     in each time step.
 * \param directed Logical constant, whether to create a directed
 *     network. 
 * \return Error code.
 * 
 * \sa \ref igraph_citing_cited_type_game() for a bit more general
 * game. 
 * 
 * Time complexity: O((|V|+|E|)log|V|), |V| and |E| are number of
 * vertices and edges, respectively.
 */

int igraph_cited_type_game(igraph_t *graph, igraph_integer_t nodes,
			   const igraph_vector_t *types,
			   const igraph_vector_t *pref,
			   igraph_integer_t edges_per_step,
			   igraph_bool_t directed) {
  
  igraph_vector_t edges;
  igraph_vector_t cumsum;
  igraph_real_t sum;
  long int i,j;
  
  IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&cumsum, 2);
  IGRAPH_CHECK(igraph_vector_reserve(&cumsum, nodes+1));
  IGRAPH_CHECK(igraph_vector_reserve(&edges, nodes*edges_per_step));
  
  /* first node */
  VECTOR(cumsum)[0]=0;
  sum=VECTOR(cumsum)[1]=VECTOR(*pref)[ (long int) VECTOR(*types)[0] ];
  
  RNG_BEGIN();

  for (i=1; i<nodes; i++) {
    for (j=0; j<edges_per_step; j++) {
      long int to;
      igraph_real_t r=RNG_UNIF(0,sum);
      igraph_vector_binsearch(&cumsum, r, &to);
        igraph_vector_push_back(&edges, i);
      igraph_vector_push_back(&edges, to-1);
    }
    sum+=VECTOR(*pref)[(long int) VECTOR(*types)[i] ];
    igraph_vector_push_back(&cumsum, sum);
  }
  
  RNG_END();

  igraph_vector_destroy(&cumsum);
  IGRAPH_FINALLY_CLEAN(1);
  IGRAPH_CHECK(igraph_create(graph, &edges, nodes, directed));
  igraph_vector_destroy(&edges);
  IGRAPH_FINALLY_CLEAN(1);

  return 0;
}


typedef struct {
  long int no;
  igraph_psumtree_t *sumtrees;
} igraph_i_citing_cited_type_game_struct_t;

void igraph_i_citing_cited_type_game_free(igraph_i_citing_cited_type_game_struct_t *s) {
  long int i;
  if (!s->sumtrees) { return; }
  for (i=0; i<s->no; i++) {
    igraph_psumtree_destroy(&s->sumtrees[i]);
  }
}

/**
 * \function igraph_citing_cited_type_game
 * \brief Simulate a citation network based on vertex types.
 *
 * This game is similar to \ref igraph_cited_type_game() but here the
 * category of the citing vertex is also considered. 
 * 
 * </para><para>
 * An evolving citation network is modeled here, a single vertex and
 * its \p edges_per_step citation are added in each time step. The 
 * odds the a given vertex is cited by the new vertex depends on the 
 * category of both the citing and the cited vertex and is given in
 * the \p pref matrix. The categories of the citing vertex correspond
 * to the rows, the categories of the cited vertex to the columns of
 * this matrix. Ie. the element in row \c i and column \c j gives the
 * probability that a \c j vertex is cited, if the category of the
 * citing vertex is \c i.
 * 
 * </para><para>
 * Note that this function might generate networks with multiple edges 
 * if \p edges_per_step is greater than one. You might want to call
 * \ref igraph_simplify() on the result to remove multiple edges.
 * \param graph Pointer to an uninitialized graph object.
 * \param nodes The number of vertices in the network.
 * \param types A numeric matrix of length \p nodes, containing the
 *    categories of the vertices. The categories are numbered from
 *    zero.
 * \param pref The preference matrix, a square matrix is required, 
 *     both the number of rows and columns should be the maximum
 *     element in \p types plus one (types are numbered from zero).
 * \param directed Logical constant, whether to create a directed
 *     network.
 * \return Error code.
 * 
 * Time complexity: O((|V|+|E|)log|V|), |V| and |E| are number of
 * vertices and edges, respectively.
 */

int igraph_citing_cited_type_game(igraph_t *graph, igraph_integer_t nodes,
				  const igraph_vector_t *types,
				  const igraph_matrix_t *pref,
				  igraph_integer_t edges_per_step,
				  igraph_bool_t directed) {

  igraph_vector_t edges;
  igraph_i_citing_cited_type_game_struct_t str = { 0, 0 };
  igraph_psumtree_t *sumtrees;
  igraph_vector_t sums;
  long int nocats=igraph_matrix_ncol(pref);
  long int i, j;
  
  IGRAPH_VECTOR_INIT_FINALLY(&edges,0);
  str.sumtrees=sumtrees=igraph_Calloc(nocats, igraph_psumtree_t);  
  if (!sumtrees) {
    IGRAPH_ERROR("Citing-cited type game failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_i_citing_cited_type_game_free, &str);
  
  for (i=0; i<nocats; i++) {
    IGRAPH_CHECK(igraph_psumtree_init(&sumtrees[i], nodes));
    str.no++;    
  }
  IGRAPH_VECTOR_INIT_FINALLY(&sums, nocats);
					       
  IGRAPH_CHECK(igraph_vector_reserve(&edges, nodes*edges_per_step));

  /* First node */
  for (i=0; i<nocats; i++) {
    long int type=VECTOR(*types)[0];
    igraph_psumtree_update(&sumtrees[i], 0, MATRIX(*pref, i, type));
    VECTOR(sums)[i]=MATRIX(*pref, i, type);
  }

  RNG_BEGIN();
    
  for (i=1; i<nodes; i++) {
    long int type=VECTOR(*types)[i];
    igraph_real_t sum=VECTOR(sums)[type];
    for (j=0; j<edges_per_step; j++) {
      long int to;
      igraph_psumtree_search(&sumtrees[type], &to, RNG_UNIF(0, sum));
      igraph_vector_push_back(&edges, i);
      igraph_vector_push_back(&edges, to);
    }
    
    /* add i */
    for (j=0; j<nocats; j++) {
      igraph_psumtree_update(&sumtrees[j], i, MATRIX(*pref, j,  type));
      VECTOR(sums)[j] += MATRIX(*pref, j, type);
    }
  }

  RNG_END();
  
  igraph_i_citing_cited_type_game_free(&str);
  IGRAPH_FINALLY_CLEAN(1);
  
  igraph_create(graph, &edges, nodes, directed);
  igraph_vector_destroy(&edges);
  IGRAPH_FINALLY_CLEAN(1);
  return 0;
}
