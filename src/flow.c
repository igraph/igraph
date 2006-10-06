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
 * \brief Maximum flow in a network with the push/relabel algorithm
 * 
 * </para><para>This function implements the Goldberg-Tarjan algorithm for
 * calculating the maximum flow in a directed graph. The algorithm was given
 * in Andrew V. Goldberg, Robert E. Tarjan: A New Approach to the
 * Maximum-Flow Problem, Journal of the ACM, 35(4), 921-940, 1988.
 * </para>
 * 
 * <para> The input of the function is a directed graph, a vector
 * of real numbers giving the capacity of the edges and two vertices
 * of the graph, the source and the target. A flow is a function 
 * assigning positive real numbers to the edges and satisfying two
 * requirements: (1) the flow value is less than the capacity of the
 * edge and (2) at each vertex except the source and the target, the
 * incoming flow (ie. the sum of the flow on the incoming edges) is
 * the same as the outgoing flow (ie. the sum of the flow on the
 * outgoing edges). The value of the flow is the incoming flow at the
 * target vertex. The maximum flow is the flow with the maximum
 * value. </para>
 * 
 * <para> This function can only calculate the value of the maximum
 * flow, but not the flow itself (may be added later). </para>
 * 
 * <para> Note that the value of the maximum flow is the same as the
 * minimum cut in the graph.
 * \param graph The input graph, it has to be directed.
 * \param value Pointer to a real number, the result will be placed here.
 * \param source The id of the source vertex.
 * \param target The id of the target vertex.
 * \param capacity Vector containing the capacity of the edges.
 * \return Error code.
 * 
 * Time complexity: O(|V|^3). In practice it is much faster, but i
 * cannot prove a better lower bound for the data structure i've
 * used. In fact, this implementation runs much faster than the
 * \c hi_pr implementation discussed in
 * B. V. Cherkassky and A. V. Goldberg: On implementing the
 * push-relabel method for the maximum flow problem, (TODO: citation)
 * on all the graph classes i've tried.
 * 
 * \sa \ref igraph_edge_connectivity_pair(), \ref
 * igraph_edge_connectivity, \ref igraph_vertex_connectivity_pair(),
 * \ref igraph_vertex_connectivity() for 
 * properties based on the maximum flow.
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
  igraph_buckets_t buckets;

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

  IGRAPH_CHECK(igraph_buckets_init(&buckets, no_of_nodes+1, no_of_nodes));
  IGRAPH_FINALLY(igraph_buckets_destroy, &buckets);

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
      igraph_buckets_add(&buckets, DIST(i), i);
    }
  }

  /* The main part comes here */
  while (!igraph_buckets_empty(&buckets)) {
    long int vertex=igraph_buckets_popmax(&buckets);
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
	      igraph_buckets_add(&buckets, DIST(nei), nei);
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
	  igraph_buckets_add(&buckets, DIST(vertex), vertex);
	}
	
	/* TODO: gap heuristics here ??? */

      } else {
	CURRENT(vertex) = FIRST(vertex);
      }

      break;

    } while (1);
  }

  /* Store the result */
  if (value) {
    *value=EXCESS(target);
  }

  igraph_buckets_destroy(&buckets);
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

/**
 * \function igraph_edge_connectivity_pair
 * \brief Edge connectivity of a pair of vertices
 * 
 * </para><para> The edge connectivity of two vertices (\c source and
 * \c target) in a directed graph is the minimum number of edges that
 * have to be deleted from the graph to eliminate all paths from \c
 * source to \c target.</para>
 * 
 * <para>This function uses the maximum flow algorithm to calculate
 * the edge connectivity.
 * \param graph The input graph, it has to be directed.
 * \param res Pointer to an integer, the result will be stored here.
 * \param source The id of the source vertex.
 * \param target The id of the target vertex.
 * \return Error code.
 *
 * Time complexity: O(|V|^3). 
 * 
 * \sa \ref igraph_maxflow(), \refigraph_edge_connectivity_pair(),
 * \ref igraph_vertex_connectivity_pair(), \ref
 * igraph_vertex_connectivity(). 
 */
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

/**
 * \function igraph_edge_connectivity
 * \brief The minimum edge connectivity in a graph.
 * 
 * </para><para>
 * The edge connectivity of a graph is the same as group adhesion as
 * defined in Douglas R. White and Frank Harary: The cohesiveness of
 * blocks in social networks: node connectivity and conditional
 * density, TODO: citation
 * \param graph The input graph, it must be directed.
 * \param res Pointer to an integer, the result will be stored here.
 * \return Error code.
 * 
 * Time complexity: O(|V|^4).
 * 
 * \sa \ref igraph_edge_connectivity_pair(), \ref igraph_maxflow(), 
 * \ref igraph_vertex_connectivity().
 */

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
    IGRAPH_CHECK(igraph_maxflow(graph, &flow, i, 0, &capacity));
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

/**
 * \function igraph_vertex_connectivity_pair
 * The vertex connectivity of a pair of vertices
 * 
 * </para><para>The vertex connectivity of two vertices (\c source and
 * \c target) is the minimum number of vertices that have to be
 * deleted to eliminate all paths from \c source to \c target.</para>
 * 
 * <para>The vertex connectivity of a pair is the same as the number
 * of different (ie. node-independent) paths from source to
 * target.</para> 
 *
 * <para>The current implementation uses a maximum flow calculation to
 * obtain the result.
 * \param graph The input graph, it must be directed.
 * \param res Pointer to an integer, the result will be stored here.
 * \param source The id of the source vertex.
 * \param target The id of the target vertex.
 * \return Error code.
 * 
 * Time complexity: O(|V|^3).
 * 
 * \sa \ref igraph_vertex_connectivity(), \ref
 * igraph_edge_connectivity_pair(), \ref igraph_edge_connectivity(), 
 * \ref igraph_maxflow().
 */

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

/**
 * \function igraph_vertex_connectivity
 * The vertex connectivity of a directed graph
 * 
 * </para><para> The vertex connectivity of a graph is the minimum
 * vertex connectivity along each pairs of vertices in the graph.
 * </para>
 * <para> The vertex connectivity of a graph is the same as group
 * cohesion as defined in Douglas R. White and Frank Harary: The
 * cohesiveness of blocks in social networks: node connectivity and
 * conditional density, TODO: citation
 * \param graph The input graph, it must be directed.
 * \param res Pointer to an integer, the result will be stored here. 
 * \return Error code.
 * 
 * Time complecity: O(|V|^4).
 * 
 * \sa \ref igraph_vertex_connectivity_pair(), \ref igraph_maxflow(),
 * \ref igraph_edge_connectivity_pair() and \ref
 * igraph_edge_connectivity(). 
 */

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
    IGRAPH_CHECK(igraph_vertex_connectivity_pair(graph, &conn, i, 0));
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

/* void print_vector(igraph_vector_t *v, FILE *f) { */
/*   long int i; */
/*   for (i=0; i<igraph_vector_size(v); i++) { */
/*     fprintf(f, " %li", (long int) VECTOR(*v)[i]); */
/*   } */
/*   fprintf(f, "\n"); */
/* } */

/* void print_adjlist(igraph_i_adjlist_t *a, long int n) { */
/*   long int i; */
/*   for (i=0; i<n; i++) { */
/*     igraph_vector_t *v=igraph_i_adjlist_get(a, i); */
/*     long int j, m=igraph_vector_size(v); */
/*     printf("%li: ", i); */
/*     for (j=0; j<m; j++) { */
/*       printf("%li ", (long int) VECTOR(*v)[j]); */
/*     } */
/*     printf("\n"); */
/*   } */
/*   printf("-------------\n"); */
/* } */

int igraph_minimum_cut(const igraph_t *graph, igraph_integer_t *res,
		       igraph_vector_t *capacity) {
  
  long int no_of_nodes=igraph_vcount(graph);
  long int no_of_edges=igraph_ecount(graph);

  igraph_i_cutheap_t heap;
  igraph_real_t mincut=1.0/0.0;	/* infinity */
  long int i;
  
  igraph_i_adjlist_t adjlist;
  igraph_i_adjedgelist_t adjedgelist;
  
  if (igraph_is_directed(graph)) {
    IGRAPH_ERROR("The minimum cut function expects undirected graphs",
		 IGRAPH_EINVAL);
  }
  
  if (igraph_vector_size(capacity) != no_of_edges) {
    IGRAPH_ERROR("Invalid capacity vector size", IGRAPH_EINVAL);
  }
  
  IGRAPH_CHECK(igraph_i_cutheap_init(&heap, no_of_nodes));
  IGRAPH_FINALLY(igraph_i_cutheap_destroy, &heap);

  IGRAPH_CHECK(igraph_i_adjedgelist_init(graph, &adjedgelist, IGRAPH_OUT));
  IGRAPH_FINALLY(igraph_i_adjedgelist_destroy, &adjedgelist);

  IGRAPH_CHECK(igraph_i_adjlist_init(graph, &adjlist, IGRAPH_OUT));
  IGRAPH_FINALLY(igraph_i_adjlist_destroy, &adjlist);

  while (igraph_i_cutheap_size(&heap) >= 2) {

    long int last;
    igraph_real_t acut;
    long int a, n;

    igraph_vector_t *edges, *neis, *edges2, *neis2;
   
/*     print_adjlist(&adjlist, no_of_nodes); */
 
    do {
      a=igraph_i_cutheap_popmax(&heap);

      /* update the weights of the active vertices connected to a */
      edges=igraph_i_adjedgelist_get(&adjedgelist, a);
      neis=igraph_i_adjlist_get(&adjlist, a);
      n=igraph_vector_size(edges);
      for (i=0; i<n; i++) {
	igraph_integer_t edge=VECTOR(*edges)[i];
	igraph_integer_t to=VECTOR(*neis)[i];
	igraph_real_t weight=VECTOR(*capacity)[(long int)edge];
	igraph_i_cutheap_update(&heap, to, weight);
/* 	print_vector(&heap.heap, stdout); */
/* 	print_vector(&heap.hptr, stdout); */
      }
            
    } while (igraph_i_cutheap_active_size(&heap) > 1);

    /* Now, there is only one active vertex left, 
       calculate the cut of the phase */
    acut=igraph_i_cutheap_maxvalue(&heap);
    last=igraph_i_cutheap_popmax(&heap);

/*     printf("acut: %g\n", acut); */
    if (acut < mincut) {
      mincut=acut;
    }    

    if (mincut == 0) {
      break;
    }

    /* And contract the last and the remaining vertex (a and last) */
    /* First remove the a--last edge if there is one, a is still the
       last deactivated vertex */
/*     printf("contracting %li into %li \n", (long int)last, (long int)a); */
    edges=igraph_i_adjedgelist_get(&adjedgelist, a);
    neis=igraph_i_adjlist_get(&adjlist, a);
    n=igraph_vector_size(edges);
    for (i=0; i<n; ) {
      if (VECTOR(*neis)[i]==last) {
	VECTOR(*neis)[i] = VECTOR(*neis)[n-1];
	VECTOR(*edges)[i] = VECTOR(*edges)[n-1];
	igraph_vector_pop_back(neis);
	igraph_vector_pop_back(edges);
	n--;
      } else {
	i++;
      }
    }
    
    edges=igraph_i_adjedgelist_get(&adjedgelist, last);
    neis=igraph_i_adjlist_get(&adjlist, last);
    n=igraph_vector_size(edges);
    for (i=0; i<n; ) {
      if (VECTOR(*neis)[i] == a) {
	VECTOR(*neis)[i] = VECTOR(*neis)[n-1];
	VECTOR(*edges)[i] = VECTOR(*edges)[n-1];
	igraph_vector_pop_back(neis);
	igraph_vector_pop_back(edges);
	n--;
      } else {
	i++;
      }
    }

    /* Now rewrite the edge lists of last's neighbors */
    neis=igraph_i_adjlist_get(&adjlist, last);
    n=igraph_vector_size(neis);    
    for (i=0; i<n; i++) {     
      igraph_integer_t nei=VECTOR(*neis)[i];
      neis2=igraph_i_adjlist_get(&adjlist, nei);
      long int n2=igraph_vector_size(neis2);
      long int j;
      for (j=0; j<n2; j++) {
	if (VECTOR(*neis2)[j] == last) {
	  VECTOR(*neis2)[j] = a;
	}
      }
    }
    
    /* And append the lists of last to the lists of a */
    edges=igraph_i_adjedgelist_get(&adjedgelist, a);
    neis=igraph_i_adjlist_get(&adjlist, a);
    edges2=igraph_i_adjedgelist_get(&adjedgelist, last);
    neis2=igraph_i_adjlist_get(&adjlist, last);
    IGRAPH_CHECK(igraph_vector_append(edges, edges2));
    IGRAPH_CHECK(igraph_vector_append(neis, neis2));
    igraph_vector_clear(edges2); /* TODO: free it */
    igraph_vector_clear(neis2);	 /* TODO: free it */

    /* Remove the deleted vertex from the heap entirely */
    igraph_i_cutheap_reset_undefine(&heap, last);    
  }

  *res=mincut;

  igraph_i_adjedgelist_destroy(&adjedgelist);
  igraph_i_cutheap_destroy(&heap);
  IGRAPH_FINALLY_CLEAN(2);
  return 0;
}

int igraph_minimum_vertex_cut_pair(const igraph_t *graph, 
				   igraph_integer_t *res,
				   igraph_integer_t source,
				   igraph_integer_t target) {
  
  long int no_of_nodes=igraph_vcount(graph);
  igraph_t newgraph;

  IGRAPH_CHECK(igraph_copy(&newgraph, graph));
  IGRAPH_FINALLY(igraph_destroy, &newgraph);
  IGRAPH_CHECK(igraph_to_directed(&newgraph, IGRAPH_TO_DIRECTED_MUTUAL));
  
  IGRAPH_CHECK(igraph_vertex_connectivity_pair(&newgraph, res, 
					       source, target));
  
  igraph_destroy(&newgraph);
  IGRAPH_FINALLY_CLEAN(1);

  return 0;
}

int igraph_minimum_vertex_cut(const igraph_t *graph,
			      igraph_integer_t *res) {
  
  long int no_of_nodes=igraph_vcount(graph);
  igraph_t newgraph;
  
  IGRAPH_CHECK(igraph_copy(&newgraph, graph));
  IGRAPH_FINALLY(igraph_destroy, &newgraph);
  IGRAPH_CHECK(igraph_to_directed(&newgraph, IGRAPH_TO_DIRECTED_MUTUAL));
  
  IGRAPH_CHECK(igraph_vertex_connectivity(&newgraph, res));
  
  igraph_destroy(&newgraph);
  IGRAPH_FINALLY_CLEAN(1);

  return 0;
}
