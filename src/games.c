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
   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

*/

#include "igraph.h"
#include "random.h"
#include "memory.h"

#include <math.h>

/**
 * \ingroup generators
 * \brief Generates a graph based on the Barabási-Albert model.
 *
 * @param graph An uninitialized graph object.
 * @param n The number of vertices in the graph.
 * @param m The number of outgoing edges generated for each
 *        vertex. (Only if <code>outseq</code> is not given.)
 * @param outseq Gives the (out-)degrees of the vertices. If this is
 *        constant, this can be an empty (but initialized!) vector,
 *        in this case <code>m</code> contains the constant
 *        out-degree. 
 * @param outpref Boolean, if true not only the in- but also the out-degree
 *        of a vertex increases its citation probability. Ie. the
 *        citation probability is determined by the total degree of
 *        the vertices.
 * @param directed Boolean, whether to generate a directed graph.
 * @return Error code.
 * 
 * Time complexity: <code>O(|V|+|E|)</code>, the number of vertices
 * plus the number of edges.
 */

int igraph_barabasi_game(igraph_t *graph, integer_t n, integer_t m, 
			 vector_t *outseq, bool_t outpref, bool_t directed) {

  long int no_of_nodes=n;
  long int no_of_neighbors=m;
  long int *bag;
  long int bagp;
  long int no_of_edges;
  vector_t edges;
  
  long int resp=0;

  long int i,j;

  if (vector_size(outseq) == 0) {
    no_of_neighbors=m;
    bag=Calloc(no_of_nodes * no_of_neighbors + no_of_nodes +
	       outpref * no_of_nodes * no_of_neighbors,
	       long int);
    no_of_edges=(no_of_nodes-1)*no_of_neighbors;
  } else {
    no_of_edges=0;
    for (i=1; i<vector_size(outseq); i++) {
      no_of_edges+=VECTOR(*outseq)[i];
    }
    bag=Calloc(no_of_nodes + no_of_edges + outpref * no_of_edges,
	       long int);
  }
  
  vector_init(&edges, no_of_edges*2);
  
  /* The first node */

  bagp=0;
  bag[bagp++]=0;
  
  RNG_BEGIN();

  /* and the others */
  
  for (i=1; i<no_of_nodes; i++) {
    /* draw edges */
    if (vector_size(outseq)!=0) { no_of_neighbors=VECTOR(*outseq)[i]; }
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
  igraph_create(graph, &edges, 0, directed);
  vector_destroy(&edges);

  return 0;
}

/**
 * \ingroup generators
 * \brief Generates a random (Erdos-Renyi) graph.
 * 
 * @param graph Pointer to an uninitialized graph object.
 * @param n The number of vertices in the graph.
 * @param p The probability that an edge is present.
 * @param directed Logical, whether to generate a directed graph.
 * @param loops Logical, whether to generate loops (self) edges.
 * @return Error code.
 * 
 * Time complexity: <code>O(|V|+|E|)</code>, the number of vertices
 * plus the number of edges in the graph.
 * 
 * \sa barabasi_game(), growing_random_game()
 */
 

int igraph_erdos_renyi_game(igraph_t *graph, integer_t n, real_t p,
			    bool_t directed, bool_t loops) {

  long int no_of_nodes=n;
  vector_t edges;
  vector_t s;
  int retval=0;
  
  if (p==0.0 || no_of_nodes<=1) {
    retval=igraph_empty(graph, n, directed);
  } else if (p==1.0) { 
    retval=igraph_full(graph, n, directed, loops);
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
    
    vector_init(&s, 0);
    vector_reserve(&s, maxedges*p*1.1);

    RNG_BEGIN();

    vector_push_back(&s, RNG_GEOM(p)+1);
    while (vector_tail(&s) < maxedges) {
      vector_push_back(&s, vector_tail(&s)+RNG_GEOM(p)+1);
    }
    if (vector_tail(&s) > maxedges+1) {
      vector_pop_back(&s);
    }

    RNG_END();

    vector_init(&edges, 0);
    vector_reserve(&edges, vector_size(&s)*2);

    if (directed && loops) {
      for (i=0; i<vector_size(&s); i++) {
	vector_push_back(&edges, ((long int)(VECTOR(s)[i]-1))/no_of_nodes);
	vector_push_back(&edges, ((long int)(VECTOR(s)[i]-1))%no_of_nodes);
      }
    } else if (directed && !loops) {
      for (i=0; i<vector_size(&s); i++) {
	long int from=((long int)(VECTOR(s)[i]-1))/(no_of_nodes-1);
	long int to=((long int)VECTOR(s)[i])%(no_of_nodes-1);
	if (from==to) {
	  to=no_of_nodes-1;
	}
	vector_push_back(&edges, from);
	vector_push_back(&edges, to);
      }
    } else if (!directed && loops) {
      for (i=0; i<vector_size(&s); i++) {
	real_t from=ceil((sqrt(8*(VECTOR(s)[i])+1)-1)/2);
	vector_push_back(&edges, from-1);
	vector_push_back(&edges, VECTOR(s)[i]-from*(from-1)/2-1);
      }
    } else {
      for (i=0; i<vector_size(&s); i++) {
	real_t from=ceil((sqrt(8*VECTOR(s)[i]+1)-1)/2)+1;
	vector_push_back(&edges, from-1);
	vector_push_back(&edges, VECTOR(s)[i]-(from-1)*(from-2)/2-1);
      }
    }      

    vector_destroy(&s);
    retval=igraph_create(graph, &edges, n, directed);
    vector_destroy(&edges);
  }

  return retval;
}

int igraph_degree_sequence_game(igraph_t *graph, vector_t *out_deg,
				vector_t *in_deg, integer_t method) {
  /* TODO */
  return 0;
}

/**
 * \ingroup generators
 * \brief Generates a growing random graph.
 *
 * This function simulates a growing random graph. In each discrete
 * time step a new vertex is added and a number of new edges are also
 * added. These graphs are known to be different from standard (not
 * growing) random graphs.
 * @param graph Uninitialized graph object.
 * @param n The number of vertices in the graph.
 * @param m The number of edges to add in a time step (ie. after
 *        adding a vertex.
 * @param directed Boolean, whether to generate a directed graph.
 * @param citation Boolean, if <code>TRUE</code> the edges always
 *        originate from the most recently added vertex.
 * @return Error code.
 *
 * Time complexity: <code>O(|V|+|E|)</code>, the number of vertices
 * plus the number of edges.
 */
int igraph_growing_random_game(igraph_t *graph, integer_t n, 
			       integer_t m, bool_t directed,
			       bool_t citation) {

  long int no_of_nodes=n;
  long int no_of_neighbors=m;
  long int no_of_edges;
  vector_t edges;
  
  long int resp=0;

  long int i,j;

  no_of_edges=(no_of_nodes-1) * no_of_neighbors;
    
  vector_init(&edges, no_of_edges*2);  

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

  igraph_create(graph, &edges, n, directed);
  vector_destroy(&edges);

  return 0;
}

int igraph_aging_prefatt_game(igraph_t *graph, integer_t n, integer_t m,
			      integer_t aging_type, real_t aging_exp) {
  /* TODO */
  return 0;
}
