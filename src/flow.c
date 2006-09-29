/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2006  Gabor Csardi <csardi@rmki.kfki.hu>
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
#include "error.h"

#include <limits.h>
#include <stdio.h>

/**
 * \function igraph_maxflow
 * \brief Maximum flow in a network
 * 
 * </para><para>This function implements the Goldberg-Tarjan algorithm for
 * calculating the maximum flow in a network. The algorithm was given
 * in Andrew V. Goldberg, Robert E. Tarjan: A New Approach to the
 * Maximum-Flow Problem, Journal of the ACM, 35(4), 921-940, 1988.
 * \param graph The input graph.
 * \param value Pointer to a real number, the result will be placed here.
 * \param source The id of the source vertex.
 * \param target The id of the target vertex.
 * \param capacity Vector containing the capacity of the edges.
 * \return Error code.
 * 
 * Time complexity: O(n^3).
 */

int igraph_maxflow(const igraph_t *graph, igraph_real_t *value,
		   igraph_integer_t source, igraph_integer_t target,
		   const igraph_vector_t *capacity) {

  long int no_of_nodes=igraph_vcount(graph);
  long int no_of_orig_edges=igraph_ecount(graph);
  long int no_of_edges=2*no_of_orig_edges;

  igraph_vector_t from, to, rev, cap, rescap, excess, distance;
  igraph_vector_t edges, rank;
  igraph_vector_t current, first;
  igraph_dqueue_t q;

  long int i, j, k, l, idx;

  if (igraph_vector_size(capacity) != no_of_orig_edges) {
    IGRAPH_ERROR("Invalid capacity vector", IGRAPH_EINVAL);
  }
  if (source<0 || source>no_of_nodes || target<0 || target>no_of_nodes) {
    IGRAPH_ERROR("Invalid source or target vertex", IGRAPH_EINVAL);
  }
  if (!igraph_is_directed(graph)) {
    IGRAPH_ERROR("The maximum flow algorithm works only on directed graphs!",
		 IGRAPH_EINVAL);
  }
  
  IGRAPH_VECTOR_INIT_FINALLY(&to,       no_of_edges);
  IGRAPH_VECTOR_INIT_FINALLY(&rev,      no_of_edges);
  IGRAPH_VECTOR_INIT_FINALLY(&rescap,   no_of_edges);
  IGRAPH_VECTOR_INIT_FINALLY(&excess,   no_of_nodes);
  IGRAPH_VECTOR_INIT_FINALLY(&distance, no_of_nodes);
  IGRAPH_VECTOR_INIT_FINALLY(&first,    no_of_nodes+1);

  IGRAPH_VECTOR_INIT_FINALLY(&from,     no_of_edges);
  IGRAPH_VECTOR_INIT_FINALLY(&edges,    no_of_edges);
  IGRAPH_VECTOR_INIT_FINALLY(&rank,     no_of_edges);
  
  /* Create the basic data structure */
  IGRAPH_CHECK(igraph_get_edgelist(graph, &edges, 0));
  IGRAPH_CHECK(igraph_vector_rank(&edges, &rank, no_of_nodes));
  
  for (i=0; i<no_of_edges; i+=2) {
    long int pos=VECTOR(rank)[i];
    long int pos2=VECTOR(rank)[i+1];
    VECTOR(from)[pos] = VECTOR(edges)[i];
    VECTOR(to)[pos]   = VECTOR(edges)[i+1];
    VECTOR(from)[pos2] = VECTOR(edges)[i+1];
    VECTOR(to)[pos2]   = VECTOR(edges)[i];
    VECTOR(rev)[pos] = pos2;
    VECTOR(rev)[pos2] = pos;
    VECTOR(rescap)[pos] = VECTOR(*capacity)[i/2];
    VECTOR(rescap)[pos2] = 0.0;
  }  
 
  igraph_vector_destroy(&rank);
  igraph_vector_destroy(&edges);
  IGRAPH_FINALLY_CLEAN(2);

  /* The first pointers */
  
  idx=-1;
  for (i=0; i<=VECTOR(from)[0]; i++) {
    idx++; VECTOR(first)[idx]=0;
  }
  for (i=1; i<no_of_edges; i++) {
    long int n=VECTOR(from)[i]-VECTOR(from)[ (long int) VECTOR(first)[idx] ];
    for (j=0; j<n; j++) {
      idx++; VECTOR(first)[idx]=i;
    }
  }
  idx++;
  while (idx < no_of_nodes+1) {
    VECTOR(first)[idx++] = no_of_edges;
  }

  igraph_vector_destroy(&from);
  IGRAPH_FINALLY_CLEAN(1);

  /* And the current pointers, initially the same as the first */
  IGRAPH_VECTOR_INIT_FINALLY(&current, no_of_nodes);
  for (i=0; i<no_of_nodes; i++) {
    VECTOR(current)[i] = VECTOR(first)[i];
  }

  /* Some useful macros */
  
#define FIRST(i)       ((long int)VECTOR(first)[(long int)(i)])
#define LAST(i)        ((long int)VECTOR(first)[(long int)(i)+1])
#define CURRENT(i)     (VECTOR(current)[(i)])
#define RESCAP(i)      (VECTOR(rescap)[(i)])
#define REV(i)         ((long int)VECTOR(rev)[(i)])
#define HEAD(i)        ((long int)VECTOR(to)[(i)])
#define EXCESS(i)      (VECTOR(excess)[(long int)(i)])
#define DIST(i)        (VECTOR(distance)[(long int)(i)])
  
  /* OK, the graph is set up, initialization */
  IGRAPH_DQUEUE_INIT_FINALLY(&q, 100);
  for (i=FIRST(source), j=LAST(source); i<j; i++) {
    if (HEAD(i) != source) {
      igraph_real_t delta=RESCAP(i);
      RESCAP(i) -= delta;
      RESCAP(REV(i)) += delta;
      EXCESS(HEAD(i)) += delta;
    }
  }

  for (i=0; i<no_of_nodes; i++) {
    DIST(i) = 1;
  }
  DIST(source)=no_of_nodes;
  DIST(target)=0;
    
  for (i=0; i<no_of_nodes; i++) {
    if (EXCESS(i) > 0.0 && i != target) {
      IGRAPH_CHECK(igraph_dqueue_push(&q, i));
    }
  }

  /* The main part comes here */
  while (!igraph_dqueue_empty(&q)) {
    long int vertex=igraph_dqueue_pop(&q);
    igraph_bool_t endoflist=0;
    /* DISCHARGE(vertex) comes here */
    do {
      for (i=CURRENT(vertex), j=LAST(vertex); i<j; i++) {
	if (RESCAP(i) > 0) {
	  long int nei=HEAD(i);
	  
	  if (DIST(vertex) == DIST(nei)+1) {
	    igraph_real_t delta=
	      RESCAP(i) < EXCESS(vertex) ? RESCAP(i) : EXCESS(vertex);
	    RESCAP(i) -= delta;
	    RESCAP(REV(i)) += delta;
	    
	    if (nei != target && EXCESS(nei) == 0.0 &&
		DIST(nei) != no_of_nodes) {
	      IGRAPH_CHECK(igraph_dqueue_push(&q, nei));
	    }
	    
	    EXCESS(nei) += delta;
	    EXCESS(vertex) -= delta;
	    
	    if (EXCESS(vertex) == 0) break;
	    
	  }
	}
      }
      
      if (i==j) {
	
	/* RELABEL(vertex) comes here */
	igraph_real_t min;
	long int min_edge;
	DIST(vertex)=min=no_of_nodes;
	for (k=FIRST(vertex), l=LAST(vertex); k<l; k++) {
	  if (RESCAP(k) > 0) {
	    if (DIST(HEAD(k)) < min) {
	      min=DIST(HEAD(k));
	      min_edge=k;
	    }
	  }
	}
	
	min++;

	if (min < no_of_nodes) {
	  DIST(vertex)=min;
	  CURRENT(vertex)=min_edge;
	  /* Vertex is still active */
	  IGRAPH_CHECK(igraph_dqueue_push(&q, vertex));
	}
	
	if (DIST(vertex) == no_of_nodes) break;

	/* TODO: gap heuristics here */
	
      } else {
	CURRENT(vertex) = FIRST(vertex);
	break;
      }
    } while (1);
  }

  /* Store the result */
  if (value) {
    *value=EXCESS(target);
  }

  igraph_dqueue_destroy(&q);
  igraph_vector_destroy(&current);
  igraph_vector_destroy(&first);
  igraph_vector_destroy(&distance);
  igraph_vector_destroy(&excess);
  igraph_vector_destroy(&rescap);
  igraph_vector_destroy(&rev);
  igraph_vector_destroy(&to);
  IGRAPH_FINALLY_CLEAN(8);

  return 0;
}

int igraph_edge_connectivity_pair(const igraph_t *graph, igraph_integer_t *res,
				  igraph_integer_t source, 
				  igraph_integer_t target) {
  
  long int no_of_edges=igraph_ecount(graph);
  igraph_vector_t capacity;
  igraph_real_t flow;
  long int i;

  IGRAPH_VECTOR_INIT_FINALLY(&capacity, no_of_edges);
  for (i=0; i<no_of_edges; i++) {
    VECTOR(capacity)[i]=1.0;
  }
  
  IGRAPH_CHECK(igraph_maxflow(graph, &flow, source, target, &capacity));
  *res = flow;

  igraph_vector_destroy(&capacity);
  IGRAPH_FINALLY_CLEAN(1);
  return 0;
}

int igraph_edge_connectivity(const igraph_t *graph, igraph_integer_t *res) {
  
  long int no_of_nodes=igraph_vcount(graph);
  long int no_of_edges=igraph_ecount(graph);
  igraph_vector_t capacity;
  igraph_real_t minmaxflow, flow;
  long int i, j;

  IGRAPH_VECTOR_INIT_FINALLY(&capacity, no_of_edges);
  for (i=0; i<no_of_edges; i++) {
    VECTOR(capacity)[i]=1.0;
  }
  
  minmaxflow=1.0/0.0;

  for (i=1; i<no_of_nodes; i++) {
    IGRAPH_CHECK(igraph_maxflow(graph, &flow, 0, i, &capacity));
    if (flow < minmaxflow) {
      minmaxflow = flow;
      if (flow==0) break;
    }
  }

  if (res) {
    *res=minmaxflow;
  }
  
  igraph_vector_destroy(&capacity);
  IGRAPH_FINALLY_CLEAN(1);
  return 0;
}

int igraph_vertex_connectivity_pair(const igraph_t *graph,
				    igraph_integer_t *res,
				    igraph_integer_t source, 
				    igraph_integer_t target) {

  long int no_of_nodes=igraph_vcount(graph);
  long int no_of_edges=igraph_ecount(graph);
  igraph_vector_t edges;
  igraph_vector_t capacity;
  igraph_t newgraph;
  long int i;
  
  /* Create the new graph */

  IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
  IGRAPH_CHECK(igraph_vector_reserve(&edges, 2*(no_of_edges+no_of_nodes)));
  IGRAPH_CHECK(igraph_get_edgelist(graph, &edges, 0));
  IGRAPH_CHECK(igraph_vector_resize(&edges, 2*(no_of_edges+no_of_nodes)));
  
  for (i=0; i<2*no_of_edges; i+=2) {
    igraph_integer_t to=VECTOR(edges)[i+1];
    if (to != source && to != target) {
      VECTOR(edges)[i+1] = no_of_nodes + to;
    }
  }
  
  for (i=0; i<no_of_nodes; i++) {
    VECTOR(edges)[ 2*(no_of_edges+i)   ] = no_of_nodes+i;
    VECTOR(edges)[ 2*(no_of_edges+i)+1 ] = i;
  }

  IGRAPH_CHECK(igraph_create(&newgraph, &edges, 2*no_of_nodes, 
			     igraph_is_directed(graph)));

  igraph_vector_destroy(&edges);
  IGRAPH_FINALLY_CLEAN(1);
  IGRAPH_FINALLY(igraph_destroy, &newgraph);
  
  /* Do the maximum flow */

  no_of_nodes=igraph_vcount(&newgraph);
  no_of_edges=igraph_ecount(&newgraph);
  
  IGRAPH_VECTOR_INIT_FINALLY(&capacity, no_of_edges);
  for (i=0; i<no_of_edges; i++) {
    VECTOR(capacity)[i] = 1.0;
  }
  
  IGRAPH_CHECK(igraph_maxflow(&newgraph, res, source, target, &capacity));
  
  igraph_vector_destroy(&capacity);
  igraph_destroy(&newgraph);
  IGRAPH_FINALLY_CLEAN(2);
  return 0;
}

int igraph_vertex_connectivity(const igraph_t *graph, igraph_integer_t *res) {

  long int no_of_nodes=igraph_vcount(graph);
  long int i;
  igraph_integer_t minconn=LONG_MAX, conn;

  for (i=1; i<no_of_nodes; i++) {
    IGRAPH_CHECK(igraph_vertex_connectivity_pair(graph, &conn, 0, i));
    if (conn < minconn) {
      minconn = conn;
      if (conn == 0) { break; }
    }
  }

  if (res) {
    *res = minconn;
  }

  return 0;
}
