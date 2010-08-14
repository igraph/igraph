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

#include "igraph_flow.h"
#include "igraph_error.h"
#include "igraph_memory.h"
#include "igraph_constants.h"
#include "igraph_interface.h"
#include "igraph_adjlist.h"
#include "igraph_conversion.h"
#include "igraph_constructors.h"
#include "igraph_structural.h"
#include "igraph_components.h"
#include "igraph_types_internal.h"
#include "config.h"
#include "igraph_math.h"
#include "igraph_dqueue.h"
#include "igraph_visitor.h"
#include "igraph_marked_queue.h"
#include "igraph_stack.h"

#include <limits.h>
#include <stdio.h>

/*
 * Some general remarks about the functions in this file.
 *
 * The following measures can be calculated:
 * ( 1) s-t maximum flow value, directed graph
 * ( 2) s-t maximum flow value, undirected graph
 * ( 3) s-t maximum flow, directed graph
 * ( 4) s-t maximum flow, undirected graph
 * ( 5) s-t minimum cut value, directed graph
 * ( 6) s-t minimum cut value, undirected graph
 * ( 7) minimum cut value, directed graph
 * ( 8) minimum cut value, undirected graph
 * ( 9) s-t minimum cut, directed graph
 * (10) s-t minimum cut, undirected graph
 * (11) minimum cut, directed graph
 * (12) minimum cut, undirected graph
 * (13) s-t edge connectivity, directed graph
 * (14) s-t edge connectivity, undirected graph
 * (15) edge connectivity, directed graph
 * (16) edge connectivity, undirected graph
 * (17) s-t vertex connectivity, directed graph
 * (18) s-t vertex connectivity, undireced graph
 * (19) vertex connectivity, directed graph
 * (20) vertex connectivity, undireced graph
 * (21) s-t number of edge disjoint paths, directed graph
 * (22) s-t number of edge disjoint paths, undirected graph
 * (23) s-t number of vertex disjoint paths, directed graph
 * (24) s-t number of vertex disjoint paths, undirected graph
 * (25) graph adhesion, directed graph
 * (26) graph adhesion, undirected graph
 * (27) graph cohesion, directed graph
 * (28) graph cohesion, undirected graph
 * 
 * This is how they are calculated:
 * ( 1) igraph_maxflow_value, calls igraph_maxflow.
 * ( 2) igraph_maxflow_value, calls igraph_maxflow, this calls
 *      igraph_i_maxflow_undirected. This transforms the graph into a
 *      directed graph, including two mutual edges instead of every
 *      undirected edge, then igraph_maxflow is called again with the
 *      directed graph.
 * ( 3) igraph_maxflow, does the push-relabel algorithm, optionally
 *      calculates the cut, the partitions and the flow itself.
 * ( 4) igraph_maxflow calls igraph_i_maxflow_undirected, this converts 
 *      the undirected graph into a directed one, adding two mutual edges
 *      for each undirected edge, then igraph_maxflow is called again, 
 *      with the directed graph. After igraph_maxflow returns, we need 
 *      to edit the flow (and the cut) to make it sense for the
 *      original graph.
 * ( 5) igraph_st_mincut_value, we just call igraph_maxflow_value
 * ( 6) igraph_st_mincut_value, we just call igraph_maxflow_value
 * ( 7) igraph_mincut_value, we call igraph_maxflow_value (|V|-1)*2
 *      times, from vertex 0 to all other vertices and from all other
 *      vertices to vertex 0
 * ( 8) We call igraph_i_mincut_value_undirected, that calls 
 *      igraph_i_mincut_undirected with partition=partition2=cut=NULL
 *      The Stoer-Wagner algorithm is used.
 * ( 9) igraph_st_mincut, just calls igraph_maxflow.
 * (10) igraph_st_mincut, just calls igraph_maxflow.
 * (11) igraph_mincut, calls igraph_i_mincut_directed, which runs 
 *      the maximum flow algorithm 2(|V|-1) times, from vertex zero to
 *      and from all other vertices and stores the smallest cut.
 * (12) igraph_mincut, igraph_i_mincut_undirected is called, 
 *      this is the Stoer-Wagner algorithm
 * (13) We just call igraph_maxflow_value, back to (1)
 * (14) We just call igraph_maxflow_value, back to (2)
 * (15) We just call igraph_mincut_value (possibly after some basic
 *      checks). Back to (7)
 * (16) We just call igraph_mincut_value (possibly after some basic
 *      checks). Back to (8).
 * (17) We call igraph_i_st_vertex_connectivity_directed.
 *      That creates a new graph with 2*|V| vertices and smartly chosen
 *      edges, so that the s-t edge connectivity of this graph is the
 *      same as the s-t vertex connectivity of the original graph.
 *      So finally it calls igraph_maxflow_value, go to (1)
 * (18) We call igraph_i_st_vertex_connectivity_undirected.
 *      We convert the graph to a directed one,
 *      IGRAPH_TO_DIRECTED_MUTUAL method. Then we call 
 *      igraph_i_st_vertex_connectivity_directed, see (17).
 * (19) We call igraph_i_vertex_connectivity_directed.
 *      That calls igraph_st_vertex_connectivity for all pairs of
 *      vertices. Back to (17).
 * (20) We call igraph_i_vertex_connectivity_undirected.
 *      That converts the graph into a directed one
 *      (IGRAPH_TO_DIRECTED_MUTUAL) and calls the directed version,
 *      igraph_i_vertex_connectivity_directed, see (19).
 * (21) igraph_edge_disjoint_paths, we just call igraph_maxflow_value, (1).
 * (22) igraph_edge_disjoint_paths, we just call igraph_maxflow_value, (2).
 * (23) igraph_vertex_disjoint_paths, if there is a connection between
 *      the two vertices, then we remove that (or all of them if there
 *      are many), as this could mess up vertex connectivity
 *      calculation. The we call
 *      igraph_i_st_vertex_connectivity_directed, see (19).
 * (24) igraph_vertex_disjoint_paths, if there is a connection between
 *      the two vertices, then we remove that (or all of them if there
 *      are many), as this could mess up vertex connectivity
 *      calculation. The we call
 *      igraph_i_st_vertex_connectivity_undirected, see (20).
 * (25) We just call igraph_edge_connectivity, see (15). 
 * (26) We just call igraph_edge_connectivity, see (16).
 * (27) We just call igraph_vertex_connectivity, see (19).
 * (28) We just call igraph_vertex_connectivity, see (20).
 */

/*
 * This is an internal function that calculates the maximum flow value
 * on undirected graphs, either for an s-t vertex pair or for the
 * graph (i.e. all vertex pairs). 
 * 
 * It does it by converting the undirected graph to a corresponding
 * directed graph, including reciprocal directed edges instead of each
 * undirected edge.
 */

int igraph_i_maxflow_undirected(const igraph_t *graph, 
				igraph_real_t *value,
				igraph_vector_t *flow,
				igraph_vector_t *cut,
				igraph_vector_t *partition,
				igraph_vector_t *partition2,
				igraph_integer_t source, 
				igraph_integer_t target,
				const igraph_vector_t *capacity) {
  long int no_of_edges=igraph_ecount(graph);
  long int no_of_nodes=igraph_vcount(graph);
  igraph_vector_t edges;
  igraph_vector_t newcapacity;
  igraph_t newgraph;
  long int i;
  
  /* We need to convert this to directed by hand, since we need to be
     sure that the edge ids will be handled properly to build the new
     capacity vector. */

  IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&newcapacity, no_of_edges*2);
  IGRAPH_CHECK(igraph_vector_reserve(&edges, no_of_edges*4));
  IGRAPH_CHECK(igraph_get_edgelist(graph, &edges, 0));
  IGRAPH_CHECK(igraph_vector_resize(&edges, no_of_edges*4));
  for (i=0; i<no_of_edges; i++) {
    VECTOR(edges)[no_of_edges*2+i*2] = VECTOR(edges)[i*2+1];
    VECTOR(edges)[no_of_edges*2+i*2+1] = VECTOR(edges)[i*2];
    VECTOR(newcapacity)[i] = VECTOR(newcapacity)[no_of_edges+i] = 
      capacity ? VECTOR(*capacity)[i] : 1.0;
  }
  
  IGRAPH_CHECK(igraph_create(&newgraph, &edges, no_of_nodes, IGRAPH_DIRECTED));
  IGRAPH_FINALLY(igraph_destroy, &newgraph);
  
  IGRAPH_CHECK(igraph_maxflow(&newgraph, value, flow, cut, partition,
			      partition2, source, target, &newcapacity));

  /* Edges are ordered by ids in the cut, we just remove the
     additional ones */
  if (cut) {
    IGRAPH_CHECK(igraph_vector_resize(cut,
				      igraph_vector_size(cut)/2));
  }
  
  /* The flow has one non-zero value for each real-nonreal edge pair,
     by definition, we convert it to a positive-negative vector. If
     for an edge the flow is negative that means that it is going
     from the bigger vertex id to the smaller one. For positive
     values the direction is the opposite. */
  if (flow) {
    long int i;
    for (i=0; i<no_of_edges; i++) {
      if (VECTOR(*flow)[i] == 0) {
	VECTOR(*flow)[i] = - VECTOR(*flow)[i+no_of_edges];
      }
    }
    IGRAPH_CHECK(igraph_vector_resize(flow, no_of_edges));
  }
  
  igraph_destroy(&newgraph);
  igraph_vector_destroy(&edges);
  igraph_vector_destroy(&newcapacity);
  IGRAPH_FINALLY_CLEAN(3);
  
  return 0;
}

/**
 * \function igraph_maxflow
 * Maximum network flow between a pair of vertices
 * 
 * </para><para>This function implements the Goldberg-Tarjan algorithm for
 * calculating value of the maximum flow in a directed or undirected
 * graph. The algorithm was given in Andrew V. Goldberg, Robert
 * E. Tarjan: A New Approach to the Maximum-Flow Problem, Journal of
 * the ACM, 35(4), 921-940, 1988. </para>
 * 
 * <para> The input of the function is a graph, a vector
 * of real numbers giving the capacity of the edges and two vertices
 * of the graph, the source and the target. A flow is a function 
 * assigning positive real numbers to the edges and satisfying two
 * requirements: (1) the flow value is less than the capacity of the
 * edge and (2) at each vertex except the source and the target, the
 * incoming flow (ie. the sum of the flow on the incoming edges) is
 * the same as the outgoing flow (ie. the sum of the flow on the
 * outgoing edges). The value of the flow is the incoming flow at the
 * target vertex. The maximum flow is the flow with the maximum
 * value.
 * 
 * \param graph The input graph, either directed or undirected.
 * \param value Pointer to a real number, the value of the maximum
 *        will be placed here, unless it is a null pointer.
 * \param flow If not a null pointer, then it must be a pointer to an
 *        initialized vector. The vector will be resized, and the flow
 *        on each edge will be placed in it, in the order of the edge
 *        ids. For undirected graphs this argument is bit trickier,
 *        since for these the flow direction is not predetermined by
 *        the edge direction. For these graphs the elements of the
 *        \p flow vector can be negative, this means that the flow
 *        goes from the bigger vertex id to the smaller one. Positive
 *        values mean that the flow goes from the smaller vertex id to
 *        the bigger one.
 * \param cut A null pointer or a pointer to an initialized vector. 
 *        If not a null pointer, then the minimum cut corresponding to
 *        the maximum flow is stored here, i.e. all edge ids that are
 *        part of the minimum cut are stored in the vector.
 * \param partition A null pointer or a pointer to an initialized
 *        vector. If not a null pointer, then the first partition of
 *        the minimum cut that corresponds to the maximum flow will be
 *        placed here.
 * \param partition2 A null pointer or a pointer to an initialized
 *        vector. If not a null pointer, then the second partition of
 *        the minimum cut that corresponds to the maximum flow will be
 *        placed here.
 * \param source The id of the source vertex.
 * \param target The id of the target vertex.
 * \param capacity Vector containing the capacity of the edges. If NULL, then
 *        every edge is considered to have capacity 1.0.
 * \return Error code.
 * 
 * Time complexity: O(|V|^3). In practice it is much faster, but i
 * cannot prove a better lower bound for the data structure i've
 * used. In fact, this implementation runs much faster than the
 * \c hi_pr implementation discussed in
 * B. V. Cherkassky and A. V. Goldberg: On implementing the 
 * push-relabel method for the maximum flow problem, (Algorithmica, 
 * 19:390--410, 1997) on all the graph classes i've tried.
 * 
 * \sa \ref igraph_mincut_value(), \ref igraph_edge_connectivity(),
 * \ref igraph_vertex_connectivity() for 
 * properties based on the maximum flow.
 */

int igraph_maxflow(const igraph_t *graph, igraph_real_t *value,
		   igraph_vector_t *flow, igraph_vector_t *cut,
		   igraph_vector_t *partition, igraph_vector_t *partition2,
		   igraph_integer_t source, igraph_integer_t target,
		   const igraph_vector_t *capacity) {

  long int no_of_nodes=igraph_vcount(graph);
  long int no_of_orig_edges=igraph_ecount(graph);
  long int no_of_edges=2*no_of_orig_edges;

  igraph_vector_t from, to, rev, rescap, excess, distance;
  igraph_vector_t edges, rank;
  igraph_vector_t current, first;
  igraph_buckets_t buckets;

  long int i, j, k, l, idx;

  if (!igraph_is_directed(graph)) {
    IGRAPH_CHECK(igraph_i_maxflow_undirected(graph, value, flow, cut,
					partition, partition2, source, 
					target, capacity));
    return 0;
  }

  if (capacity && igraph_vector_size(capacity) != no_of_orig_edges) {
    IGRAPH_ERROR("Invalid capacity vector", IGRAPH_EINVAL);
  }
  if (source<0 || source>=no_of_nodes || target<0 || target>=no_of_nodes) {
    IGRAPH_ERROR("Invalid source or target vertex", IGRAPH_EINVAL);
  }

  /* 
   * The data structure:
   * - First of all, we consider every edge twice, first the edge
   *   itself, but also its opposite.
   * - (from, to) contain all edges (original + opposite), ordered by 
   *   the id of the source vertex. During the algorithm we just need
   *   'to', so from is destroyed soon. We only need it in the
   *   beginning, to create the 'first' pointers.
   * - 'first' is a pointer vector for 'to', first[i] points to the
   *   first neighbor of vertex i and first[i+1]-1 is the last
   *   neighbor of vertex i. (Unless vertex i is isolate, in which
   *   case first[i]==first[i+1]).
   * - 'rev' contains a mapping from an edge to its opposite pair
   * - 'rescap' contains the residual capacities of the edges, this is
   *   initially equal to the capacity of the edges for the original
   *   edges and it is zero for the opposite edges.
   * - 'excess' contains the excess flow for the vertices. I.e. the flow
   *   that is coming in, but it is not going out.
   * - 'current' stores the next neighboring vertex to check, for every
   *   vertex, when excess flow is being pushed to neighbors.
   * - 'distance' stores the distance of the vertices from the source.
   * - 'rank' and 'edges' are only needed temporarily, for ordering and
   *   storing the edges.
   * - we use an igraph_buckets_t data structure ('buckets') to find
   *   the vertices with the highest 'distance' values quickly.
   *   This always contains the vertices that have a positive excess
   *   flow.
   */
  
  IGRAPH_VECTOR_INIT_FINALLY(&to,       no_of_edges);
  IGRAPH_VECTOR_INIT_FINALLY(&rev,      no_of_edges);
  IGRAPH_VECTOR_INIT_FINALLY(&rescap,   no_of_edges);
  IGRAPH_VECTOR_INIT_FINALLY(&excess,   no_of_nodes);
  IGRAPH_VECTOR_INIT_FINALLY(&distance, no_of_nodes);
  IGRAPH_VECTOR_INIT_FINALLY(&first,    no_of_nodes+1);

  IGRAPH_VECTOR_INIT_FINALLY(&rank,     no_of_edges);
  IGRAPH_VECTOR_INIT_FINALLY(&from,     no_of_edges);
  IGRAPH_VECTOR_INIT_FINALLY(&edges,    no_of_edges);
  
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
    VECTOR(rescap)[pos] = capacity ? VECTOR(*capacity)[i/2] : 1.0;
    VECTOR(rescap)[pos2] = 0.0;
  }  
 
  /* The first pointers. This is a but trickier, than one would
     think, because of the possible isolate vertices. */
  
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
  igraph_vector_destroy(&edges);
  IGRAPH_FINALLY_CLEAN(2);

  if (!flow) {
    igraph_vector_destroy(&rank);
    IGRAPH_FINALLY_CLEAN(1);
  }

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

  /* Send as much flow as possible from the source to its neighbors */
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

  /* It would be enough to do this for the neighbors of the source,
     the rest have EXCESS(i) == 0.0 by definition. */
  for (i=0; i<no_of_nodes; i++) {
    if (EXCESS(i) > 0.0 && i != target) {
      igraph_buckets_add(&buckets, DIST(i), i);
    }
  }

  /* The main part comes here */
  while (!igraph_buckets_empty(&buckets)) {
    long int vertex=igraph_buckets_popmax(&buckets);
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
	long int min_edge=0;
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

  /* If we also need the minimum cut */
  if (cut || partition || partition2) {
    /* We need to find all vertices from which the target is reachable 
       in the residual graph. We do a breadth-first search, going
       backwards. */
    igraph_dqueue_t Q;
    igraph_vector_bool_t added;
    long int marked=0;

    IGRAPH_CHECK(igraph_vector_bool_init(&added, no_of_nodes));
    IGRAPH_FINALLY(igraph_vector_bool_destroy, &added);

    IGRAPH_CHECK(igraph_dqueue_init(&Q, 100));
    IGRAPH_FINALLY(igraph_dqueue_destroy, &Q);

    igraph_dqueue_push(&Q, target);
    VECTOR(added)[(long int)target]=1;
    marked++;
    while (!igraph_dqueue_empty(&Q)) {
      long int actnode=igraph_dqueue_pop(&Q);
      for (i=FIRST(actnode), j=LAST(actnode); i<j; i++) {
	long int nei=HEAD(i);
	if (!VECTOR(added)[nei] && RESCAP(REV(i)) > 0.0) {
	  VECTOR(added)[nei]=1;
	  marked++;
	  IGRAPH_CHECK(igraph_dqueue_push(&Q, nei));
	}
      }
    }    
    igraph_dqueue_destroy(&Q);
    IGRAPH_FINALLY_CLEAN(1);

    /* Now we marked each vertex that is on one side of the cut,
       check the crossing edges */

    if (cut) {
      igraph_vector_clear(cut);
      for (i=0; i<no_of_orig_edges; i++) {
	long int v1=IGRAPH_FROM(graph, i);
	long int v2=IGRAPH_TO(graph, i);
	char p1=VECTOR(added)[v1];
	char p2=VECTOR(added)[v2];
	if (p1 ^ p2) {
	  IGRAPH_CHECK(igraph_vector_push_back(cut, i));
	}
      }
    }

    if (partition) {
      long int x=0;
      IGRAPH_CHECK(igraph_vector_resize(partition, marked));
      for (i=0; i<no_of_nodes; i++) {
	if (VECTOR(added)[i]) {
	  VECTOR(*partition)[x++]=i;
	}
      }
    }

    if (partition2) {
      long int x=0;
      IGRAPH_CHECK(igraph_vector_resize(partition2,
					no_of_nodes-marked));
      for (i=0; i<no_of_nodes; i++) {
	if (!VECTOR(added)[i]) {
	  VECTOR(*partition2)[x++]=i;
	}
      }
    }
    
    igraph_vector_bool_destroy(&added);
    IGRAPH_FINALLY_CLEAN(1);
  }

  if (flow) {
    /* Initialize the backward distances, with a breadth-first search 
       from the source */ 
    igraph_dqueue_t Q;
    igraph_vector_bool_t added;
    long int j;

    IGRAPH_CHECK(igraph_vector_bool_init(&added, no_of_nodes));
    IGRAPH_FINALLY(igraph_vector_bool_destroy, &added);
    IGRAPH_CHECK(igraph_dqueue_init(&Q, 100));
    IGRAPH_FINALLY(igraph_dqueue_destroy, &added);
    
    igraph_dqueue_push(&Q, source);
    igraph_dqueue_push(&Q, 0);
    VECTOR(added)[(long int)source]=1;
    while (!igraph_dqueue_empty(&Q)) {
      long int actnode=igraph_dqueue_pop(&Q);
      long int actdist=igraph_dqueue_pop(&Q);
      DIST(actnode)=actdist;
	  
      for (i=FIRST(actnode), j=LAST(actnode); i<j; i++) {
	long int nei=HEAD(i);
	if (!VECTOR(added)[nei] && RESCAP(REV(i)) > 0.0) {
	  VECTOR(added)[nei]=1;
	  IGRAPH_CHECK(igraph_dqueue_push(&Q, nei));
	  IGRAPH_CHECK(igraph_dqueue_push(&Q, actdist+1));
	}
      }
    } /* !igraph_dqueue_empty(&Q) */
	  
    igraph_vector_bool_destroy(&added);
    igraph_dqueue_destroy(&Q);
    IGRAPH_FINALLY_CLEAN(2);

    /* Reinitialize the buckets */
    igraph_buckets_clear(&buckets);
    for (i=0; i<no_of_nodes; i++) {
      if (EXCESS(i) > 0.0 && i != source && i != target) {
	igraph_buckets_add(&buckets, DIST(i), i);
      }
    }

    /* Now we return the flow to the source */
    while (!igraph_buckets_empty(&buckets)) {
      long int vertex=igraph_buckets_popmax(&buckets);

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
	      
	      if (nei != source && EXCESS(nei) == 0.0 &&
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
	  long int min_edge=0;
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
    
    IGRAPH_CHECK(igraph_vector_resize(flow, no_of_orig_edges));
    for (i=0, j=0; i<no_of_edges; i+=2, j++) {
      long int pos=VECTOR(rank)[i];
      VECTOR(*flow)[j] = (capacity ? VECTOR(*capacity)[j] : 1.0) - 
	RESCAP(pos);
    }
    
    igraph_vector_destroy(&rank);
    IGRAPH_FINALLY_CLEAN(1);
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
 * \function igraph_maxflow_value
 * \brief Maximum flow in a network with the push/relabel algorithm
 * 
 * </para><para>This function implements the Goldberg-Tarjan algorithm for
 * calculating value of the maximum flow in a directed or undirected
 * graph. The algorithm was given in Andrew V. Goldberg, Robert
 * E. Tarjan: A New Approach to the Maximum-Flow Problem, Journal of
 * the ACM, 35(4), 921-940, 1988. </para>
 * 
 * <para> The input of the function is a graph, a vector
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
 * <para> According to a theorem by Ford and Furkelson 
 * (L. R. Ford Jr. and D. R. Fulkerson. Maximal flow through a
 * network. Canadian J. Math., 8:399-404, 1956.) the maximum flow
 * between two vertices is the same as the 
 * minimum cut between them (also called the minimum s-t cut). So \ref
 * igraph_st_mincut_value() gives the same result in all cases as \c
 * igraph_maxflow_value().</para>
 * 
 * <para> Note that the value of the maximum flow is the same as the
 * minimum cut in the graph.
 * \param graph The input graph, either directed or undirected.
 * \param value Pointer to a real number, the result will be placed here.
 * \param source The id of the source vertex.
 * \param target The id of the target vertex.
 * \param capacity Vector containing the capacity of the edges. If NULL, then
 *        every edge is considered to have capacity 1.0.
 * \return Error code.
 * 
 * Time complexity: O(|V|^3). In practice it is much faster, but i
 * cannot prove a better lower bound for the data structure i've
 * used. In fact, this implementation runs much faster than the
 * \c hi_pr implementation discussed in
 * B. V. Cherkassky and A. V. Goldberg: On implementing the 
 * push-relabel method for the maximum flow problem, (Algorithmica, 
 * 19:390--410, 1997) on all the graph classes i've tried.
 * 
 * \sa \ref igraph_maxflow() to calculate the actual flow. 
 * \ref igraph_mincut_value(), \ref igraph_edge_connectivity(),
 * \ref igraph_vertex_connectivity() for 
 * properties based on the maximum flow.
 */

int igraph_maxflow_value(const igraph_t *graph, igraph_real_t *value,
			 igraph_integer_t source, igraph_integer_t target,
			 const igraph_vector_t *capacity) {

  return igraph_maxflow(graph, value, /*flow=*/ 0, /*cut=*/ 0, 
			/*partition=*/ 0, /*partition1=*/ 0,
			source, target, capacity);
}

/**
 * \function igraph_st_mincut_value
 * \brief The minimum s-t cut in a graph
 * 
 * </para><para> The minimum s-t cut in a weighted (=valued) graph is the
 * total minimum edge weight needed to remove from the graph to
 * eliminate all paths from a given vertex (\c source) to
 * another vertex (\c target). Directed paths are considered in
 * directed graphs, and undirected paths in undirected graphs.  </para>
 * 
 * <para> The minimum s-t cut between two vertices is known to be same
 * as the maximum flow between these two vertices. So this function
 * calls \ref igraph_maxflow_value() to do the calculation.
 * \param graph The input graph.
 * \param value Pointer to a real variable, the result will be stored
 *        here. 
 * \param source The id of the source vertex.
 * \param target The id of the target vertex.
 * \param capacity Pointer to the capacity vector, it should contain
 *        non-negative numbers and its length should be the same the
 *        the number of edges in the graph. It can be a null pointer, then
 *        every edge has unit capacity.
 * \return Error code.
 * 
 * Time complexity: O(|V|^3), see also the discussion for \ref
 * igraph_maxflow_value(), |V| is the number of vertices. 
 */

int igraph_st_mincut_value(const igraph_t *graph, igraph_real_t *value,
			   igraph_integer_t source, igraph_integer_t target,
			   const igraph_vector_t *capacity) {
  
  if (source == target) {
    IGRAPH_ERROR("source and target vertices are the same", IGRAPH_EINVAL);
  }
  
  IGRAPH_CHECK(igraph_maxflow_value(graph, value, source, target, capacity));
  return 0;
}			    

/** 
 * \function igraph_st_mincut
 * Minimum cut between a source and a target vertex
 * 
 * Finds the edge set that has the smallest total capacity among all
 * edge sets that disconnect the source and target vertices.
 * 
 * </para><para>The calculation is performed using maximum flow
 * techniques, by calling \ref igraph_maxflow().
 * \param graph The input graph.
 * \param value Pointer to a real variable, the value of the cut is
 *        stored here.
 * \param cut Pointer to a real vector, the edge ids that are included
 *        in the cut are stored here. This argument is ignored if it
 *        is a null pointer.
 * \param partition Pointer to a real vector, the vertex ids of the
 *        vertices in the first partition of the cut are stored
 *        here. This argument is ignored if it is a null pointer.
 * \param partition2 Pointer to a real vector, the vertex ids of the
 *        vertices in the second partinio of the cut are stored here.
 *        This argument is ignored if it is a null pointer.
 * \param source Integer, the id of the source vertex.
 * \param target Integer, the id of the target vertex.
 * \param capacity Vector containing the capacity of the edges. If a
 *        null pointer, then every edge is considered to have capacity
 *        1.0.
 * \return Error code.
 * 
 * \sa \ref igraph_maxflow().
 * 
 * Time complecity: see \ref igraph_maxflow().
 */

int igraph_st_mincut(const igraph_t *graph, igraph_real_t *value,
		     igraph_vector_t *cut, igraph_vector_t *partition,
		     igraph_vector_t *partition2,
		     igraph_integer_t source, igraph_integer_t target,
		     const igraph_vector_t *capacity) {
  
  return igraph_maxflow(graph, value, /*flow=*/ 0, 
			cut, partition, partition2, 
			source, target, capacity);
}

/* This is a flow-based version, but there is a better one
   for undirected graphs */

/* int igraph_i_mincut_value_undirected(const igraph_t *graph, */
/* 				     igraph_real_t *res, */
/* 				     const igraph_vector_t *capacity) { */
  
/*   long int no_of_edges=igraph_ecount(graph); */
/*   long int no_of_nodes=igraph_vcount(graph); */
/*   igraph_vector_t edges; */
/*   igraph_vector_t newcapacity; */
/*   igraph_t newgraph; */
/*   long int i; */
  
/*   /\* We need to convert this to directed by hand, since we need to be */
/*      sure that the edge ids will be handled properly to build the new */
/*      capacity vector. *\/ */

/*   IGRAPH_VECTOR_INIT_FINALLY(&edges, 0); */
/*   IGRAPH_VECTOR_INIT_FINALLY(&newcapacity, no_of_edges*2); */
/*   IGRAPH_CHECK(igraph_vector_reserve(&edges, no_of_edges*4)); */
/*   IGRAPH_CHECK(igraph_get_edgelist(graph, &edges, 0)); */
/*   IGRAPH_CHECK(igraph_vector_resize(&edges, no_of_edges*4)); */
/*   for (i=0; i<no_of_edges; i++) { */
/*     VECTOR(edges)[no_of_edges*2+i*2] = VECTOR(edges)[i*2+1]; */
/*     VECTOR(edges)[no_of_edges*2+i*2+1] = VECTOR(edges)[i*2]; */
/*     VECTOR(newcapacity)[i] = VECTOR(newcapacity)[no_of_edges+i] =  */
/*       capacity ? VECTOR(*capacity)[i] : 1.0 ; */
/*   } */
  
/*   IGRAPH_CHECK(igraph_create(&newgraph, &edges, no_of_nodes, IGRAPH_DIRECTED)); */
/*   IGRAPH_FINALLY(igraph_destroy, &newgraph); */
  
/*   IGRAPH_CHECK(igraph_mincut_value(&newgraph, res, &newcapacity)); */
  
/*   igraph_destroy(&newgraph); */
/*   igraph_vector_destroy(&edges); */
/*   igraph_vector_destroy(&newcapacity); */
/*   IGRAPH_FINALLY_CLEAN(3); */
  
/*   return 0; */
/* } */

/*
 * This is the Stoer-Wagner algorithm, it works for calcuating the
 * minimum cut for undirected graphs, for the whole graph. 
 * I.e. this is basically the edge-connectivity of the graph. 
 * It can also calculate the cut itself, not just the cut value.
 */

int igraph_i_mincut_undirected(const igraph_t *graph, 
			       igraph_integer_t *res,
			       igraph_vector_t *partition,
			       igraph_vector_t *partition2,
			       igraph_vector_t *cut,
			       const igraph_vector_t *capacity) {
  
  long int no_of_nodes=igraph_vcount(graph);
  long int no_of_edges=igraph_ecount(graph);  

  igraph_i_cutheap_t heap;
  igraph_real_t mincut=IGRAPH_INFINITY;	/* infinity */
  long int i;
  
  igraph_adjlist_t adjlist;
  igraph_adjedgelist_t adjedgelist;

  igraph_vector_t mergehist;
  igraph_bool_t calc_cut=partition || partition2 || cut;
  long int act_step=0, mincut_step=0;
  
  if (capacity && igraph_vector_size(capacity) != no_of_edges) {
    IGRAPH_ERROR("Invalid capacity vector size", IGRAPH_EINVAL);
  }

  if (calc_cut) {
    IGRAPH_VECTOR_INIT_FINALLY(&mergehist, 0);
    IGRAPH_CHECK(igraph_vector_reserve(&mergehist, no_of_nodes*2));
  }

  IGRAPH_CHECK(igraph_i_cutheap_init(&heap, no_of_nodes));
  IGRAPH_FINALLY(igraph_i_cutheap_destroy, &heap);

  IGRAPH_CHECK(igraph_adjedgelist_init(graph, &adjedgelist, IGRAPH_OUT));
  IGRAPH_FINALLY(igraph_adjedgelist_destroy, &adjedgelist);

  IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist, IGRAPH_OUT));
  IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist);

  while (igraph_i_cutheap_size(&heap) >= 2) {

    long int last;
    igraph_real_t acut;
    long int a, n;

    igraph_vector_t *edges, *neis, *edges2, *neis2;
   
    do {
      a=igraph_i_cutheap_popmax(&heap);

      /* update the weights of the active vertices connected to a */
      edges=igraph_adjedgelist_get(&adjedgelist, a);
      neis=igraph_adjlist_get(&adjlist, a);
      n=igraph_vector_size(edges);
      for (i=0; i<n; i++) {
	igraph_integer_t edge=VECTOR(*edges)[i];
	igraph_integer_t to=VECTOR(*neis)[i];
	igraph_real_t weight=capacity ? VECTOR(*capacity)[(long int)edge] : 1.0;
	igraph_i_cutheap_update(&heap, to, weight);
      }
            
    } while (igraph_i_cutheap_active_size(&heap) > 1);

    /* Now, there is only one active vertex left, 
       calculate the cut of the phase */
    acut=igraph_i_cutheap_maxvalue(&heap);
    last=igraph_i_cutheap_popmax(&heap);

    if (acut < mincut) {
      mincut=acut;
      mincut_step=act_step;
    }    

    if (mincut == 0) {
      break;
    }

    /* And contract the last and the remaining vertex (a and last) */
    /* Before actually doing that, make some notes */
    act_step++;
    if (calc_cut) {
      IGRAPH_CHECK(igraph_vector_push_back(&mergehist, a));
      IGRAPH_CHECK(igraph_vector_push_back(&mergehist, last));
    }
    /* First remove the a--last edge if there is one, a is still the
       last deactivated vertex */
    edges=igraph_adjedgelist_get(&adjedgelist, a);
    neis=igraph_adjlist_get(&adjlist, a);
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
    
    edges=igraph_adjedgelist_get(&adjedgelist, last);
    neis=igraph_adjlist_get(&adjlist, last);
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
    neis=igraph_adjlist_get(&adjlist, last);
    n=igraph_vector_size(neis);    
    for (i=0; i<n; i++) {     
      igraph_integer_t nei=VECTOR(*neis)[i];
      long int n2, j;
      neis2=igraph_adjlist_get(&adjlist, nei);
      n2=igraph_vector_size(neis2);
      for (j=0; j<n2; j++) {
	if (VECTOR(*neis2)[j] == last) {
	  VECTOR(*neis2)[j] = a;
	}
      }
    }
    
    /* And append the lists of last to the lists of a */
    edges=igraph_adjedgelist_get(&adjedgelist, a);
    neis=igraph_adjlist_get(&adjlist, a);
    edges2=igraph_adjedgelist_get(&adjedgelist, last);
    neis2=igraph_adjlist_get(&adjlist, last);
    IGRAPH_CHECK(igraph_vector_append(edges, edges2));
    IGRAPH_CHECK(igraph_vector_append(neis, neis2));
    igraph_vector_clear(edges2); /* TODO: free it */
    igraph_vector_clear(neis2);	 /* TODO: free it */

    /* Remove the deleted vertex from the heap entirely */
    igraph_i_cutheap_reset_undefine(&heap, last);    
  }

  *res=mincut;

  igraph_adjedgelist_destroy(&adjedgelist);
  igraph_adjlist_destroy(&adjlist);
  igraph_i_cutheap_destroy(&heap);
  IGRAPH_FINALLY_CLEAN(3);

  if (calc_cut) {
    long int bignode=VECTOR(mergehist)[2*mincut_step+1];
    long int i, idx;
    long int size=1;
    char *mark;
    mark=igraph_Calloc(no_of_nodes, char);
    if (!mark) { 
      IGRAPH_ERROR("Not enough memory for minumum cut", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, mark);
    
    /* first count the vertices in the partition */
    mark[bignode]=1;
    for (i=mincut_step-1; i>=0; i--) {
      if ( mark[ (long int) VECTOR(mergehist)[2*i] ] ) {
	size++;
	mark [ (long int) VECTOR(mergehist)[2*i+1] ]=1;
      }
    }
	
    /* now store them, if requested */
    if (partition) {
      IGRAPH_CHECK(igraph_vector_resize(partition, size));
      idx=0;
      VECTOR(*partition)[idx++]=bignode;
      for (i=mincut_step-1; i>=0; i--) {
	if (mark[ (long int) VECTOR(mergehist)[2*i] ]) {
	  VECTOR(*partition)[idx++] = VECTOR(mergehist)[2*i+1];
	}
      }
    }

    /* The other partition too? */
    if (partition2) {
      IGRAPH_CHECK(igraph_vector_resize(partition2, no_of_nodes-size));
      idx=0;
      for (i=0; i<no_of_nodes; i++) {
	if (!mark[i]) {
	  VECTOR(*partition2)[idx++]=i;
	}
      }
    }
    
    /* The edges in the cut are also requested? */
    /* We want as few memory allocated for 'cut' as possible,
       so we first collect the edges in mergehist, we don't 
       need that anymore. Then we copy it to 'cut';  */
    if (cut) {
      igraph_integer_t from, to;
      igraph_vector_clear(&mergehist);
      for (i=0; i<no_of_edges; i++) {
	igraph_edge(graph, i, &from, &to);
	if ((mark[(long int)from] && !mark[(long int)to]) ||
	    (mark[(long int)to] && !mark[(long int)from])) {
	  IGRAPH_CHECK(igraph_vector_push_back(&mergehist, i));
	}
      }
      igraph_vector_clear(cut);
      IGRAPH_CHECK(igraph_vector_append(cut, &mergehist));
    }

    igraph_free(mark);
    igraph_vector_destroy(&mergehist);
    IGRAPH_FINALLY_CLEAN(2);
  }
  
  return 0;
}

int igraph_i_mincut_directed(const igraph_t *graph,
			     igraph_real_t *value,
			     igraph_vector_t *partition,
			     igraph_vector_t *partition2,
			     igraph_vector_t *cut,
			     const igraph_vector_t *capacity) {
  long int i;
  long int no_of_nodes=igraph_vcount(graph);
  igraph_real_t flow;
  igraph_real_t minmaxflow=IGRAPH_INFINITY;
  igraph_vector_t mypartition, mypartition2, mycut;
  igraph_vector_t *ppartition=0, *ppartition2=0, *pcut=0;
  igraph_vector_t bestpartition, bestpartition2, bestcut;
  igraph_vector_t *pbestpartition=0, *pbestpartition2=0, *pbestcut=0;

  if (partition) {
    IGRAPH_VECTOR_INIT_FINALLY(&bestpartition, 0);
    pbestpartition=&bestpartition;
  }
  if (partition2) {
    IGRAPH_VECTOR_INIT_FINALLY(&bestpartition2, 0);
    pbestpartition2=&bestpartition2;
  }
  if (cut) {
    IGRAPH_VECTOR_INIT_FINALLY(&bestcut, 0);
    pbestcut=&bestcut;
  }
  
  if (partition) {
    IGRAPH_VECTOR_INIT_FINALLY(&mypartition, 0);
    ppartition=&mypartition;
  }
  if (partition2) {
    IGRAPH_VECTOR_INIT_FINALLY(&mypartition2, 0);
    ppartition2=&mypartition2;
  }
  if (cut) {
    IGRAPH_VECTOR_INIT_FINALLY(&mycut, 0);
    pcut=&mycut;
  }

  for (i=1; i<no_of_nodes; i++) {
    IGRAPH_CHECK(igraph_maxflow(graph, /*value=*/ &flow, /*flow=*/ 0, 
				pcut, ppartition, ppartition2, /*source=*/ 0,
				/*target=*/ i, capacity));
    if (flow < minmaxflow) {
      minmaxflow = flow;
      if (cut) { 
	IGRAPH_CHECK(igraph_vector_update(&bestcut, &mycut)); 
      }
      if (partition) { 
	IGRAPH_CHECK(igraph_vector_update(&bestpartition, &mypartition)); 
      }
      if (partition2) { 
	IGRAPH_CHECK(igraph_vector_update(&bestpartition2, &mypartition2)); 
      }

      if (minmaxflow == 0) { break; }
    }
    IGRAPH_CHECK(igraph_maxflow(graph, /*value=*/ &flow, /*flow=*/ 0,
				pcut, ppartition, ppartition2, /*source=*/ i,
				/*target=*/ 0, capacity));
    if (flow < minmaxflow) {
      minmaxflow = flow;
      if (cut) { 
	IGRAPH_CHECK(igraph_vector_update(&bestcut, &mycut)); 
      }
      if (partition) { 
	IGRAPH_CHECK(igraph_vector_update(&bestpartition, &mypartition)); 
      }
      if (partition2) { 
	IGRAPH_CHECK(igraph_vector_update(&bestpartition2, &mypartition2)); 
      }

      if (minmaxflow == 0) { break; }
    }
  }
  
  if (value) {
    *value = minmaxflow;
  }

  if (cut) {
    igraph_vector_destroy(&mycut);
    IGRAPH_FINALLY_CLEAN(1);
  }
  if (partition) {
    igraph_vector_destroy(&mypartition);
    IGRAPH_FINALLY_CLEAN(1);
  }
  if (partition2) {
    igraph_vector_destroy(&mypartition2);
    IGRAPH_FINALLY_CLEAN(1);
  }
  if (cut) {
    IGRAPH_CHECK(igraph_vector_update(cut, &bestcut));
    igraph_vector_destroy(&bestcut);
    IGRAPH_FINALLY_CLEAN(1);
  }  
  if (partition2) {
    IGRAPH_CHECK(igraph_vector_update(partition2, &bestpartition2));
    igraph_vector_destroy(&bestpartition2);
    IGRAPH_FINALLY_CLEAN(1);
  }
  if (partition) {
    IGRAPH_CHECK(igraph_vector_update(partition, &bestpartition));
    igraph_vector_destroy(&bestpartition);
    IGRAPH_FINALLY_CLEAN(1);
  }
  
  return 0;
}

/** 
 * \function igraph_mincut
 * \brief Calculates the minimum cut in a graph.
 * 
 * This function calculates the minimum cut in a graph. 
 * The minimum cut is the minimum set of edges which needs to be
 * removed to disconnect the graph. The minimum is calculated using
 * the weigths (\p capacity) of the edges, so the cut with the minimum
 * total capacity is calculated. 
 * 
 * </para><para> For directed graphs an implementation based on
 * calculating 2|V|-2 maximum flows is used. 
 * For undirected graphs we use the Stoer-Wagner
 * algorithm, as described in M. Stoer and F. Wagner: A simple min-cut
 * algorithm, Journal of the ACM, 44 585-591, 1997.
 * 
 * </para><para>
 * The first implementation of the actual cut calculation for
 * undirected graphs was made by Gregory Benison, thanks Greg.
 * \param graph The input graph.
 * \param value Pointer to an integer, the value of the cut will be
 *    stored here.
 * \param partition Pointer to an initialized vector, the ids
 *    of the vertices in the first partition after separating the
 *    graph will be stored here. The vector will be resized as
 *    needed. This argument is ignored if it is a NULL pointer. 
 * \param partition2 Pointer to an initialized vector the ids
 *    of the vertices in the second partition will be stored here. 
 *    The vector will be resized as needed. This argument is ignored
 *    if it is a NULL pointer.
 * \param cut Pointer to an initialized vector, the ids of the edges
 *    in the cut will be stored here. This argument is ignored if it
 *    is a NULL pointer.
 * \param capacity A numeric vector giving the capacities of the
 *    edges. If a null pointer then all edges have unit capacity.
 * \return Error code.
 * 
 * \sa \ref igraph_mincut_value(), a simpler interface for calculating
 * the value of the cut only.
 *
 * Time complexity: for directed graphs it is O(|V|^4), but see the
 * remarks at \ref igraph_maxflow(). For undirected graphs it is
 * O(|V||E|+|V|^2 log|V|). |V| and |E| are the number of vertices and
 * edges respectively.
 */

int igraph_mincut(const igraph_t *graph,
		  igraph_real_t *value,
		  igraph_vector_t *partition,
		  igraph_vector_t *partition2,
		  igraph_vector_t *cut,
		  const igraph_vector_t *capacity) {
  
  if (igraph_is_directed(graph)) {
    if (partition || partition2 || cut) {
      igraph_i_mincut_directed(graph, value, partition, partition2, cut, 
			       capacity);
    } else {
      return igraph_mincut_value(graph, value, capacity);      
    }
  } else {
    return igraph_i_mincut_undirected(graph, value, partition, 
				      partition2, cut, capacity);
  }
  
  return 0;
}
    

int igraph_i_mincut_value_undirected(const igraph_t *graph, 
				     igraph_integer_t *res,
				     const igraph_vector_t *capacity) {
  return igraph_i_mincut_undirected(graph, res, 0, 0, 0, capacity);
}

/** 
 * \function igraph_mincut_value
 * \brief The minimum edge cut in a graph
 * 
 * </para><para> The minimum edge cut in a graph is the total minimum
 * weight of the edges needed to remove from the graph to make the
 * graph \em not strongly connected. (If the original graph is not
 * strongly connected then this is zero.) Note that in undirected
 * graphs strong connectedness is the same as weak connectedness. </para>
 * 
 * <para> The minimum cut can be calculated with maximum flow
 * techniques, although the current implementation does this only for
 * directed graphs and a separate non-flow based implementation is
 * used for undirected graphs. See Mechthild Stoer and Frank Wagner: A
 * simple min-cut algorithm, Journal of the ACM 44 585--591, 1997.
 * For directed graphs
 * the maximum flow is calculated between a fixed vertex and all the
 * other vertices in the graph and this is done in both
 * directions. Then the minimum is taken to get the minimum cut.
 * 
 * \param graph The input graph. 
 * \param res Pointer to a real variable, the result will be stored
 *    here.
 * \param capacity Pointer to the capacity vector, it should contain
 *    the same number of non-negative numbers as the number of edges in
 *    the graph. If a null pointer then all edges will have unit capacity.
 * \return Error code.
 *
 * \sa \ref igraph_mincut(), \ref igraph_maxflow_value(), \ref
 * igraph_st_mincut_value(). 
 * 
 * Time complexity: O(log(|V|)*|V|^2) for undirected graphs and 
 * O(|V|^4) for directed graphs, but see also the discussion at the
 * documentation of \ref igraph_maxflow_value().
 */

int igraph_mincut_value(const igraph_t *graph, igraph_real_t *res, 
			const igraph_vector_t *capacity) {

  long int no_of_nodes=igraph_vcount(graph);
  igraph_real_t minmaxflow, flow;
  long int i;

  minmaxflow=IGRAPH_INFINITY;

  if (!igraph_is_directed(graph)) {
    IGRAPH_CHECK(igraph_i_mincut_value_undirected(graph, res, capacity));
    return 0;
  }    

  for (i=1; i<no_of_nodes; i++) {
    IGRAPH_CHECK(igraph_maxflow_value(graph, &flow, 0, i, capacity));
    if (flow < minmaxflow) {
      minmaxflow = flow;
      if (flow==0) break;
    }
    IGRAPH_CHECK(igraph_maxflow_value(graph, &flow, i, 0, capacity));
    if (flow < minmaxflow) {
      minmaxflow = flow;
      if (flow==0) break;
    }
  }

  if (res) {
    *res=minmaxflow;
  }
  
  return 0;
}

int igraph_i_st_vertex_connectivity_directed(const igraph_t *graph,
					     igraph_integer_t *res,
					     igraph_integer_t source, 
					     igraph_integer_t target,
					     igraph_vconn_nei_t neighbors) {
  
  long int no_of_nodes=igraph_vcount(graph);
  long int no_of_edges=igraph_ecount(graph);
  igraph_vector_t edges;
  igraph_t newgraph;
  long int i;
  igraph_bool_t conn1;
  
  if (source<0 || source>=no_of_nodes || target<0 || target>=no_of_nodes) {
    IGRAPH_ERROR("Invalid source or target vertex", IGRAPH_EINVAL);
  }

  switch (neighbors) {
  case IGRAPH_VCONN_NEI_ERROR:
    IGRAPH_CHECK(igraph_are_connected(graph, source, target, &conn1));
    if (conn1) {
      IGRAPH_ERROR("vertices connected", IGRAPH_EINVAL);
      return 0;
    }
    break;
  case IGRAPH_VCONN_NEI_INFINITY:
    IGRAPH_CHECK(igraph_are_connected(graph, source, target, &conn1));
    if (conn1) {
/*       fprintf(stderr, "%li -> %li connected\n", (long int)source, (long int) target); */
      *res=IGRAPH_INFINITY;
      return 0;
/*     } else { */
/*       fprintf(stderr, "not connected\n"); */
    }
    break;
  case IGRAPH_VCONN_NEI_IGNORE:
    break;
  default:
    IGRAPH_ERROR("Unknown `igraph_vconn_nei_t'", IGRAPH_EINVAL);
    break;
  }
  
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
  
  IGRAPH_CHECK(igraph_maxflow_value(&newgraph, res, 
				    source, target, 0));
  
  igraph_destroy(&newgraph);
  IGRAPH_FINALLY_CLEAN(1);
  
  return 0;
}

int igraph_i_st_vertex_connectivity_undirected(const igraph_t *graph, 
					       igraph_integer_t *res,
					       igraph_integer_t source,
					       igraph_integer_t target,
					       igraph_vconn_nei_t neighbors){
  
  long int no_of_nodes=igraph_vcount(graph);
  igraph_t newgraph;
  igraph_bool_t conn;
  
  if (source<0 || source>=no_of_nodes || target<0 || target>=no_of_nodes) {
    IGRAPH_ERROR("Invalid source or target vertex", IGRAPH_EINVAL);
  }

  switch (neighbors) {
  case IGRAPH_VCONN_NEI_ERROR:
    IGRAPH_CHECK(igraph_are_connected(graph, source, target, &conn));
    if (conn) { 
      IGRAPH_ERROR("vertices connected", IGRAPH_EINVAL);      
      return 0;
    }
    break;
  case IGRAPH_VCONN_NEI_INFINITY:
    IGRAPH_CHECK(igraph_are_connected(graph, source, target, &conn));
    if (conn) {
      *res=IGRAPH_INFINITY;
      return 0;
    }
    break;
  case IGRAPH_VCONN_NEI_IGNORE:
    break;
  default:
    IGRAPH_ERROR("Unknown `igraph_vconn_nei_t'", IGRAPH_EINVAL);
    break;
  }

  IGRAPH_CHECK(igraph_copy(&newgraph, graph));
  IGRAPH_FINALLY(igraph_destroy, &newgraph);
  IGRAPH_CHECK(igraph_to_directed(&newgraph, IGRAPH_TO_DIRECTED_MUTUAL));
  
  IGRAPH_CHECK(igraph_i_st_vertex_connectivity_directed(&newgraph, res, 
							source, target, 
							IGRAPH_VCONN_NEI_IGNORE));
  
  igraph_destroy(&newgraph);
  IGRAPH_FINALLY_CLEAN(1);

  return 0;
}

/**
 * \function igraph_st_vertex_connectivity
 * \brief The vertex connectivity of a pair of vertices
 * 
 * </para><para>The vertex connectivity of two vertices (\c source and
 * \c target) is the minimum number of vertices that have to be
 * deleted to eliminate all paths from \c source to \c
 * target. Directed paths are considered in directed graphs.</para>
 * 
 * <para>The vertex connectivity of a pair is the same as the number
 * of different (ie. node-independent) paths from source to
 * target.</para> 
 *
 * <para>The current implementation uses maximum flow calculations to
 * obtain the result.
 * \param graph The input graph.
 * \param res Pointer to an integer, the result will be stored here.
 * \param source The id of the source vertex.
 * \param target The id of the target vertex.
 * \param neighbors A constant giving what to do if the two vertices
 *     are connected. Possible values: 
 *     \c IGRAPH_VCONN_NEI_ERROR, stop with an error message,
 *     \c IGRAPH_VCONN_INFINITY, return infinity (ie. 1.0/0.0).
 *     \c IGRAPH_VCONN_IGNORE, ignore the fact that the two vertices
 *        are connected and calculated the number of vertices needed
 *        to aliminate all paths except for the trivial (direct) paths
 *        between \c source and \c vertex. TOOD: what about neighbors?
 * \return Error code.
 * 
 * Time complexity: O(|V|^3), but see the discussion at \ref
 * igraph_maxflow_value(). 
 * 
 * \sa \ref igraph_vertex_connectivity(),
 * \ref igraph_edge_connectivity(),
 * \ref igraph_maxflow_value().
 */

int igraph_st_vertex_connectivity(const igraph_t *graph, 
				  igraph_integer_t *res,
				  igraph_integer_t source,
				  igraph_integer_t target,
				  igraph_vconn_nei_t neighbors) {

  if (source == target) { 
    IGRAPH_ERROR("source and target vertices are the same", IGRAPH_EINVAL);
  }
  
  if (igraph_is_directed(graph)) {
    IGRAPH_CHECK(igraph_i_st_vertex_connectivity_directed(graph, res,
							  source, target,
							  neighbors));
  } else {
    IGRAPH_CHECK(igraph_i_st_vertex_connectivity_undirected(graph, res,
							    source, target,
							    neighbors));
  }
  
  return 0;
}

int igraph_i_vertex_connectivity_directed(const igraph_t *graph, 
					igraph_integer_t *res) {

  long int no_of_nodes=igraph_vcount(graph);
  long int i, j;
  igraph_integer_t minconn=no_of_nodes-1, conn;

  for (i=0; i<no_of_nodes; i++) {
    for (j=0; j<no_of_nodes; j++) {
      if (i==j) { continue; }
      IGRAPH_CHECK(igraph_st_vertex_connectivity(graph, &conn, i, j, 
						 IGRAPH_VCONN_NEI_INFINITY));
      if (conn < minconn) {
	minconn = conn;
	if (conn == 0) { break; }
      }
    }
    if (conn == 0) { break; }
  }

  if (res) {
    *res = minconn;
  }

  return 0;
}

int igraph_i_vertex_connectivity_undirected(const igraph_t *graph, 
					    igraph_integer_t *res) {
  igraph_t newgraph;

  IGRAPH_CHECK(igraph_copy(&newgraph, graph));
  IGRAPH_FINALLY(igraph_destroy, &newgraph);
  IGRAPH_CHECK(igraph_to_directed(&newgraph, IGRAPH_TO_DIRECTED_MUTUAL));
  
  IGRAPH_CHECK(igraph_i_vertex_connectivity_directed(&newgraph, res));
  
  igraph_destroy(&newgraph);
  IGRAPH_FINALLY_CLEAN(1);

  return 0;  
}

/* Use that vertex.connectivity(G) <= edge.connectivity(G) <= min(degree(G)) */
int igraph_i_connectivity_checks(const igraph_t *graph,
				 igraph_integer_t *res,
				 igraph_bool_t *found) {
  igraph_bool_t conn;
  *found=0;
  IGRAPH_CHECK(igraph_is_connected(graph, &conn, IGRAPH_STRONG));
  if (!conn) {
    *res=0;
    *found=1;
  } else {
    igraph_vector_t degree;
    IGRAPH_VECTOR_INIT_FINALLY(&degree, 0);
    if (!igraph_is_directed(graph)) {
      IGRAPH_CHECK(igraph_degree(graph, &degree, igraph_vss_all(),
				 IGRAPH_OUT, IGRAPH_LOOPS));
      if (igraph_vector_min(&degree)==1) {
	*res=1;
	*found=1;
      }
    } else {
      /* directed, check both in- & out-degree */
      IGRAPH_CHECK(igraph_degree(graph, &degree, igraph_vss_all(),
				 IGRAPH_OUT, IGRAPH_LOOPS));
      if (igraph_vector_min(&degree)==1) {
	*res=1;
	*found=1;
      } else {
	IGRAPH_CHECK(igraph_degree(graph, &degree, igraph_vss_all(),
				   IGRAPH_IN, IGRAPH_LOOPS));
	if (igraph_vector_min(&degree)==1) {
	  *res=1;
	  *found=1;
	}
      }
    }
    igraph_vector_destroy(&degree);
    IGRAPH_FINALLY_CLEAN(1);
  }
  return 0;
}

/**
 * \function igraph_vertex_connectivity
 * The vertex connectivity of a graph
 * 
 * </para><para> The vertex connectivity of a graph is the minimum
 * vertex connectivity along each pairs of vertices in the graph.
 * </para>
 * <para> The vertex connectivity of a graph is the same as group
 * cohesion as defined in Douglas R. White and Frank Harary: The
 * cohesiveness of blocks in social networks: node connectivity and
 * conditional density, Sociological Methodology 31:305--359, 2001.
 * \param graph The input graph.
 * \param res Pointer to an integer, the result will be stored here. 
 * \param checks Logical constant. Whether to check that the graph is
 *    connected and also the degree of the vertices. If the graph is
 *    not (strongly) connected then the connectivity is obviously zero. Otherwise
 *    if the minimum degree is one then the vertex connectivity is also
 *    one. It is a good idea to perform these checks, as they can be
 *    done quickly compared to the connectivity calculation itself. 
 *    They were suggested by Peter McMahan, thanks Peter.
 * \return Error code.
 * 
 * Time complecity: O(|V|^5).
 * 
 * \sa \ref igraph_st_vertex_connectivity(), \ref igraph_maxflow_value(),
 * and \ref igraph_edge_connectivity(). 
 */

int igraph_vertex_connectivity(const igraph_t *graph, igraph_integer_t *res, 
			       igraph_bool_t checks) {

  igraph_bool_t ret=0;

  if (checks) {
    IGRAPH_CHECK(igraph_i_connectivity_checks(graph, res, &ret));
  }
  
  /* Are we done yet? */
  if (!ret) {
    if (igraph_is_directed(graph)) {
      IGRAPH_CHECK(igraph_i_vertex_connectivity_directed(graph, res));
    } else {
      IGRAPH_CHECK(igraph_i_vertex_connectivity_undirected(graph, res));
    }
  }

  return 0;
}

/**
 * \function igraph_st_edge_connectivity
 * \brief Edge connectivity of a pair of vertices
 * 
 * </para><para> The edge connectivity of two vertices (\c source and
 * \c target) in a graph is the minimum number of edges that
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
 * \sa \ref igraph_maxflow_value(), \ref igraph_edge_connectivity(),
 * \ref igraph_st_vertex_connectivity(), \ref
 * igraph_vertex_connectivity().
 */

int igraph_st_edge_connectivity(const igraph_t *graph, igraph_integer_t *res,
				igraph_integer_t source, 
				igraph_integer_t target) {
  igraph_real_t flow;

  if (source == target) {
    IGRAPH_ERROR("source and target vertices are the same", IGRAPH_EINVAL);
  }

  IGRAPH_CHECK(igraph_maxflow_value(graph, &flow, source, target, 0));
  *res = flow;

  return 0;
}


/**
 * \function igraph_edge_connectivity
 * \brief The minimum edge connectivity in a graph.
 * 
 * </para><para> This is the minimum of the edge connectivity over all
 * pairs of vertices in the graph. </para>
 * 
 * <para>
 * The edge connectivity of a graph is the same as group adhesion as
 * defined in Douglas R. White and Frank Harary: The cohesiveness of
 * blocks in social networks: node connectivity and conditional
 * density, Sociological Methodology 31:305--359, 2001.
 * \param graph The input graph.
 * \param res Pointer to an integer, the result will be stored here.
 * \param checks Logical constant. Whether to check that the graph is
 *    connected and also the degree of the vertices. If the graph is
 *    not (strongly) connected then the connectivity is obviously zero. Otherwise
 *    if the minimum degree is one then the edge connectivity is also
 *    one. It is a good idea to perform these checks, as they can be
 *    done quickly compared to the connectivity calculation itself. 
 *    They were suggested by Peter McMahan, thanks Peter.
 * \return Error code.
 * 
 * Time complexity: O(log(|V|)*|V|^2) for undirected graphs and 
 * O(|V|^4) for directed graphs, but see also the discussion at the
 * documentation of \ref igraph_maxflow_value().
 * 
 * \sa \ref igraph_st_edge_connectivity(), \ref igraph_maxflow_value(), 
 * \ref igraph_vertex_connectivity().
 */

int igraph_edge_connectivity(const igraph_t *graph, igraph_integer_t *res,
			     igraph_bool_t checks) {
  
  igraph_bool_t ret;
  
  /* Use that vertex.connectivity(G) <= edge.connectivity(G) <= min(degree(G)) */
  if (checks) {
    IGRAPH_CHECK(igraph_i_connectivity_checks(graph, res, &ret));
  }  

  if (!ret) {
    IGRAPH_CHECK(igraph_mincut_value(graph, res, 0));
  }

  return 0;
}

/**
 * \function igraph_edge_disjoint_paths
 * \brief The maximum number of edge-disjoint paths between two vertices. 
 * 
 * </para><para> A set of paths between two vertices is called
 * edge-disjoint if they do not share any edges. The maximum number of
 * edge-disjoint paths are calculated by this function using maximum
 * flow techniques. Directed paths are considered in directed
 * graphs. </para>
 * 
 * <para> Note that the number of disjoint paths is the same as the
 * edge connectivity of the two vertices using uniform edge weights.
 * \param graph The input graph, can be directed or undirected.
 * \param res Pointer to an integer variable, the result will be
 *        stored here. 
 * \param source The id of the source vertex.
 * \param target The id of the target vertex.
 * \return Error code.
 * 
 * Time complecity: O(|V|^3), but see the discussion at \ref
 * igraph_maxflow_value().
 * 
 * \sa \ref igraph_vertex_disjoint_paths(), \ref
 * igraph_st_edge_connectivity(), \ref igraph_maxflow_value().
 */

int igraph_edge_disjoint_paths(const igraph_t *graph, igraph_integer_t *res,
			       igraph_integer_t source, 
			       igraph_integer_t target) {

  igraph_real_t flow;

  if (source == target) {
    IGRAPH_ERROR("Not implemented for source=target", IGRAPH_UNIMPLEMENTED);
  }

  IGRAPH_CHECK(igraph_maxflow_value(graph, &flow, source, target, 0));

  *res = flow;
  
  return 0;
}

/**
 * \function igraph_vertex_disjoint_paths
 * \brief Maximum number of vertex-disjoint paths between two vertices.
 * 
 * </para><para> A set of paths between two vertices is called
 * vertex-disjoint if they share no vertices. The calculation is
 * performed by using maximum flow techniques. </para>
 * 
 * <para> Note that the number of vertex-disjoint paths is the same as
 * the vertex connectivity of the two vertices in most cases (if the
 * two vertices are not connected by an edge).
 * \param graph The input graph.
 * \param res Pointer to an integer variable, the result will be
 *        stored here. 
 * \param source The id of the source vertex.
 * \param target The id of the target vertex.
 * \return Error code.
 * 
 * Time complexity: O(|V|^3).
 * 
 * \sa \ref igraph_edge_disjoint_paths(), \ref
 * igraph_vertex_connectivity(), \ref igraph_maxflow_value().
 */

int igraph_vertex_disjoint_paths(const igraph_t *graph, igraph_integer_t *res,
				 igraph_integer_t source,
				 igraph_integer_t target) {

  igraph_bool_t conn;

  if (source==target) {
    IGRAPH_ERROR("The source==target case is not implemented",
		 IGRAPH_UNIMPLEMENTED);
  }

  igraph_are_connected(graph, source, target, &conn);
  if (conn) { 
    /* We need to remove every (possibly directed) edge between source
       and target and calculate the disjoint paths on the new
       graph. Finally we add 1 for the removed connection(s).  */
    igraph_es_t es;
    igraph_vector_t v;
    igraph_t newgraph;
    IGRAPH_VECTOR_INIT_FINALLY(&v, 2);
    VECTOR(v)[0]=source;
    VECTOR(v)[1]=target;
    IGRAPH_CHECK(igraph_es_multipairs(&es, &v, IGRAPH_DIRECTED));
    IGRAPH_FINALLY(igraph_es_destroy, &es);
    
    IGRAPH_CHECK(igraph_copy(&newgraph, graph));    
    IGRAPH_FINALLY(igraph_destroy, &newgraph);
    IGRAPH_CHECK(igraph_delete_edges(&newgraph, es));

    if (igraph_is_directed(graph)) {
      IGRAPH_CHECK(igraph_i_st_vertex_connectivity_directed(&newgraph, res,
							    source, target,
							    IGRAPH_VCONN_NEI_IGNORE));
    } else {
      IGRAPH_CHECK(igraph_i_st_vertex_connectivity_undirected(&newgraph, res,
							      source, target, 
							      IGRAPH_VCONN_NEI_IGNORE));
    }

    if (res) {
      *res += 1;
    }
    
    IGRAPH_FINALLY_CLEAN(3);
    igraph_destroy(&newgraph);
    igraph_es_destroy(&es);
    igraph_vector_destroy(&v);
  }

  /* These do nothing if the two vertices are connected, 
     so it is safe to call them. */

  if (igraph_is_directed(graph)) {
    IGRAPH_CHECK(igraph_i_st_vertex_connectivity_directed(graph, res,
							  source, target,
							  IGRAPH_VCONN_NEI_IGNORE));
  } else {
    IGRAPH_CHECK(igraph_i_st_vertex_connectivity_undirected(graph, res,
							    source, target,
							    IGRAPH_VCONN_NEI_IGNORE));
  }    
  
  return 0;
}

/**
 * \function igraph_adhesion
 * \brief Graph adhesion, this is (almost) the same as edge connectivity.
 * 
 * </para><para> This quantity is defined by White and Harary in
 * The cohesiveness of blocks in social networks: node connectivity and
 * conditional density, (Sociological Methodology 31:305--359, 2001)
 * and basically it is the edge connectivity of the graph
 * with uniform edge weights.
 * \param graph The input graph, either directed or undirected.
 * \param res Pointer to an integer, the result will be stored here.
 * \param checks Logical constant. Whether to check that the graph is
 *    connected and also the degree of the vertices. If the graph is
 *    not (strongly) connected then the adhesion is obviously zero. Otherwise
 *    if the minimum degree is one then the adhesion is also
 *    one. It is a good idea to perform these checks, as they can be
 *    done quickly compared to the edge connectivity calculation itself. 
 *    They were suggested by Peter McMahan, thanks Peter.
* \return Error code.
 * 
 * Time complexity: O(log(|V|)*|V|^2) for undirected graphs and 
 * O(|V|^4) for directed graphs, but see also the discussion at the
 * documentation of \ref igraph_maxflow_value().
 *
 * \sa \ref igraph_cohesion(), \ref igraph_maxflow_value(), \ref
 * igraph_edge_connectivity(), \ref igraph_mincut_value().
 */

int igraph_adhesion(const igraph_t *graph, igraph_integer_t *res,
		    igraph_bool_t checks) {
  return igraph_edge_connectivity(graph, res, checks);
}

/**
 * \function igraph_cohesion
 * \brief Graph cohesion, this is the same as vertex connectivity. 
 * 
 * </para><para> This quantity was defined by White and Harary in <quote>The
 * cohesiveness of blocks in social networks: node connectivity and
 * conditional density</quote>, (Sociological Methodology 31:305--359, 2001)
 * and it is the same as the vertex connectivity of a 
 * graph. 
 * \param graph The input graph.
 * \param res Pointer to an integer variable, the result will be
 *        stored here.
 * \param checks Logical constant. Whether to check that the graph is
 *    connected and also the degree of the vertices. If the graph is
 *    not (strongly) connected then the cohesion is obviously zero. Otherwise
 *    if the minimum degree is one then the cohesion is also
 *    one. It is a good idea to perform these checks, as they can be
 *    done quickly compared to the vertex connectivity calculation itself. 
 *    They were suggested by Peter McMahan, thanks Peter.
 * \return Error code.
 * 
 * Time complexity: O(|V|^4), |V| is the number of vertices. In
 * practice it is more like O(|V|^2), see \ref igraph_maxflow_value().
 * 
 * \sa \ref igraph_vertex_connectivity(), \ref igraph_adhesion(), 
 * \ref igraph_maxflow_value().
 */

int igraph_cohesion(const igraph_t *graph, igraph_integer_t *res,
		    igraph_bool_t checks) {
  
  IGRAPH_CHECK(igraph_vertex_connectivity(graph, res, checks));
  return 0;
}

/**
 * \function igraph_even_tarjan_reduction
 * Even-Tarjan reduction of a graph
 */

int igraph_even_tarjan_reduction(const igraph_t *graph, igraph_t *graphbar) {

  long int no_of_nodes=igraph_vcount(graph);
  long int no_of_edges=igraph_ecount(graph);
  
  long int new_no_of_nodes=no_of_nodes*2;
  long int new_no_of_edges=no_of_nodes + no_of_edges * 2;
  
  igraph_vector_t edges;
  long int edgeptr=0;
  long int i;

  IGRAPH_VECTOR_INIT_FINALLY(&edges, new_no_of_edges * 2);

  /* Every vertex 'i' is replaced by two vertices, i' and i'' */
  /* id[i'] := id[i] ; id[i''] := id[i] + no_of_nodes */

  /* One edge for each original vertex, for i, we add (i',i'') */
  for (i=0; i<no_of_nodes; i++) {
    VECTOR(edges)[edgeptr++] = i;
    VECTOR(edges)[edgeptr++] = i + no_of_nodes;
  }
  
  /* Two news edges for each original edge 
     (from,to) becomes (from'',to'), (to'',from') */
  for (i=0; i<no_of_edges; i++) {
    long int from=IGRAPH_FROM(graph, i);
    long int to=IGRAPH_TO(graph, i);
    VECTOR(edges)[edgeptr++] = from + no_of_nodes;
    VECTOR(edges)[edgeptr++] = to;
    VECTOR(edges)[edgeptr++] = to + no_of_nodes;
    VECTOR(edges)[edgeptr++] = from;
  }
  
  IGRAPH_CHECK(igraph_create(graphbar, &edges, new_no_of_nodes, 
			     IGRAPH_DIRECTED));

  igraph_vector_destroy(&edges);
  IGRAPH_FINALLY_CLEAN(1);
  
  return 0;
}

int igraph_i_residual_graph(const igraph_t *graph,
			    const igraph_vector_t *capacity,
			    igraph_t *residual,
			    igraph_vector_t *residual_capacity,
			    const igraph_vector_t *flow, 
			    igraph_vector_t *tmp) {
  
  long int no_of_nodes=igraph_vcount(graph);
  long int no_of_edges=igraph_ecount(graph);
  long int i, no_new_edges=0;
  long int edgeptr=0, capptr=0;
  
  for (i=0; i<no_of_edges; i++) {
    if (VECTOR(*flow)[i] < VECTOR(*capacity)[i]) {
      no_new_edges++;
    }
  }
  
  IGRAPH_CHECK(igraph_vector_resize(tmp, no_new_edges*2));
  if (residual_capacity) {
    IGRAPH_CHECK(igraph_vector_resize(residual_capacity, no_new_edges));
  }
  
  for (i=0; i<no_of_edges; i++) {
    if (VECTOR(*capacity)[i] - VECTOR(*flow)[i] > 0) {
      long int from=IGRAPH_FROM(graph, i);
      long int to=IGRAPH_TO(graph, i);
      igraph_real_t c=VECTOR(*capacity)[i];
      VECTOR(*tmp)[edgeptr++] = from;
      VECTOR(*tmp)[edgeptr++] = to;
      if (residual_capacity) {
	VECTOR(*residual_capacity)[capptr++] = c;
      }
    }
  }

  IGRAPH_CHECK(igraph_create(residual, tmp, no_of_nodes, IGRAPH_DIRECTED));  

  return 0;
}

int igraph_residual_graph(const igraph_t *graph,
			  const igraph_vector_t *capacity,
			  igraph_t *residual,
			  igraph_vector_t *residual_capacity,
			  const igraph_vector_t *flow) {
  
  igraph_vector_t tmp;
  long int no_of_edges=igraph_ecount(graph);
  
  if (igraph_vector_size(capacity) != no_of_edges) {
    IGRAPH_ERROR("Invalid `capacity' vector size", IGRAPH_EINVAL);
  }
  if (igraph_vector_size(flow) != no_of_edges) {
    IGRAPH_ERROR("Invalid `flow' vector size", IGRAPH_EINVAL);
  }

  IGRAPH_VECTOR_INIT_FINALLY(&tmp, 0);

  IGRAPH_CHECK(igraph_i_residual_graph(graph, capacity, residual, 
				       residual_capacity, flow, &tmp));
  
  igraph_vector_destroy(&tmp);
  IGRAPH_FINALLY_CLEAN(1);
  
  return 0;
}

int igraph_i_inverse_residual_graph(const igraph_t *graph,
				    const igraph_vector_t *capacity,
				    igraph_t *residual,
				    const igraph_vector_t *flow,
				    igraph_vector_t *tmp) {

  long int no_of_nodes=igraph_vcount(graph);
  long int no_of_edges=igraph_ecount(graph);
  long int i, no_new_edges=0;
  long int edgeptr=0;

  for (i=0; i<no_of_edges; i++) {
    if (VECTOR(*flow)[i] > 0) {
      no_new_edges++;
    }
    if (VECTOR(*flow)[i] < (capacity ? VECTOR(*capacity)[i] : 1.0)) {
      no_new_edges++;
    }
  }
  
  IGRAPH_CHECK(igraph_vector_resize(tmp, no_new_edges*2));

  for (i=0; i<no_of_edges; i++) {
    if (VECTOR(*flow)[i] > 0) {
      long int from=IGRAPH_FROM(graph, i);
      long int to=IGRAPH_TO(graph, i);
      VECTOR(*tmp)[edgeptr++] = from;
      VECTOR(*tmp)[edgeptr++] = to;
    }
    if (VECTOR(*flow)[i] < (capacity ? VECTOR(*capacity)[i] : 1.0)) {
      long int from=IGRAPH_TO(graph, i);
      long int to=IGRAPH_FROM(graph, i);
      VECTOR(*tmp)[edgeptr++] = from;
      VECTOR(*tmp)[edgeptr++] = to;
    }
  }
  
  IGRAPH_CHECK(igraph_create(residual, tmp, no_of_nodes, IGRAPH_DIRECTED));
  
  return 0;
}
  
int igraph_inverse_residual_graph(const igraph_t *graph,
				  const igraph_vector_t *capacity,
				  igraph_t *residual,
				  const igraph_vector_t *flow) {
  igraph_vector_t tmp;
  long int no_of_edges=igraph_ecount(graph);
  
  if (capacity && igraph_vector_size(capacity) != no_of_edges) {
    IGRAPH_ERROR("Invalid `capacity' vector size", IGRAPH_EINVAL);
  }
  if (igraph_vector_size(flow) != no_of_edges) {
    IGRAPH_ERROR("Invalid `flow' vector size", IGRAPH_EINVAL);
  }
  IGRAPH_VECTOR_INIT_FINALLY(&tmp, 0);
  
  IGRAPH_CHECK(igraph_i_inverse_residual_graph(graph, capacity, residual,
					       flow, &tmp));
  
  igraph_vector_destroy(&tmp);
  IGRAPH_FINALLY_CLEAN(1);
  
  return 0;
}

typedef struct igraph_i_dbucket_t {
  igraph_vector_long_t head;
  igraph_vector_long_t next;
} igraph_i_dbucket_t;

int igraph_i_dbucket_init(igraph_i_dbucket_t *buck, long int size) {
  IGRAPH_CHECK(igraph_vector_long_init(&buck->head, size));
  IGRAPH_FINALLY(igraph_vector_long_destroy, &buck->head);
  IGRAPH_CHECK(igraph_vector_long_init(&buck->next, size));
  IGRAPH_FINALLY_CLEAN(1);
  return 0;
}

void igraph_i_dbucket_destroy(igraph_i_dbucket_t *buck) {
  igraph_vector_long_destroy(&buck->head);
  igraph_vector_long_destroy(&buck->next);
}

int igraph_i_dbucket_insert(igraph_i_dbucket_t *buck, long int bid, 
			    long int elem) {
  /* Note: we can do this, since elem is not in any buckets */
  VECTOR(buck->next)[elem]=VECTOR(buck->head)[bid];
  VECTOR(buck->head)[bid]=elem+1;
  return 0;
}

long int igraph_i_dbucket_empty(const igraph_i_dbucket_t *buck, 
				long int bid) {
  return VECTOR(buck->head)[bid] == 0;
}

long int igraph_i_dbucket_delete(igraph_i_dbucket_t *buck, long int bid) {
  long int elem=VECTOR(buck->head)[bid]-1;
  VECTOR(buck->head)[bid]=VECTOR(buck->next)[elem];
  return elem;
}

int igraph_i_dominator_LINK(long int v, long int w,
			    igraph_vector_long_t *ancestor) {
  VECTOR(*ancestor)[w] = v+1;
  return 0;
}

/* TODO: don't always reallocate path */

int igraph_i_dominator_COMPRESS(long int v,
				igraph_vector_long_t *ancestor,
				igraph_vector_long_t *label,
				igraph_vector_long_t *semi) {
  igraph_stack_long_t path;
  long int w=v;
  long int top, pretop;
    
  IGRAPH_CHECK(igraph_stack_long_init(&path, 10));
  IGRAPH_FINALLY(igraph_stack_long_destroy, &path);
  
  while (VECTOR(*ancestor)[w] != 0) {
    IGRAPH_CHECK(igraph_stack_long_push(&path, w));
    w=VECTOR(*ancestor)[w]-1;
  }
  
  top=igraph_stack_long_pop(&path);
  while (!igraph_stack_long_empty(&path)) {
    pretop=igraph_stack_long_pop(&path);

    if (VECTOR(*semi)[VECTOR(*label)[top]] < 
	VECTOR(*semi)[VECTOR(*label)[pretop]]) {
      VECTOR(*label)[pretop] = VECTOR(*label)[top];
    }
    VECTOR(*ancestor)[pretop]=VECTOR(*ancestor)[top];

    top=pretop;
  }

  igraph_stack_long_destroy(&path);
  IGRAPH_FINALLY_CLEAN(1);

  return 0;
}

long int igraph_i_dominator_EVAL(long int v,
				 igraph_vector_long_t *ancestor,
				 igraph_vector_long_t *label,
				 igraph_vector_long_t *semi) {
  if (VECTOR(*ancestor)[v] == 0) { 
    return v;
  } else {
    igraph_i_dominator_COMPRESS(v, ancestor, label, semi);
    return VECTOR(*label)[v];
  }
}

/* TODO: implement the faster version. */

int igraph_dominator_tree(const igraph_t *graph,
			  igraph_integer_t root,
			  igraph_vector_t *dom,
			  igraph_t *domtree,
			  igraph_vector_t *leftout,
			  igraph_neimode_t mode) {

  long int no_of_nodes=igraph_vcount(graph);

  igraph_adjlist_t succ, pred;
  igraph_vector_t parent;
  igraph_vector_long_t semi;	/* +1 always */
  igraph_vector_t vertex;	/* +1 always */
  igraph_i_dbucket_t bucket;
  igraph_vector_long_t ancestor;
  igraph_vector_long_t label;

  igraph_neimode_t invmode= mode==IGRAPH_IN ? IGRAPH_OUT: IGRAPH_IN;

  long int i;

  igraph_vector_t vdom, *mydom=dom;

  long int component_size=0;

  if (root < 0 || root >= no_of_nodes) {
    IGRAPH_ERROR("Invalid root vertex id for dominator tree", 
		 IGRAPH_EINVAL);
  }

  if (!igraph_is_directed(graph)) {
    IGRAPH_ERROR("Dominator tree of an undirected graph requested",
		 IGRAPH_EINVAL);
  }
  
  if (mode == IGRAPH_ALL) {
    IGRAPH_ERROR("Invalid neighbor mode for dominator tree",
		 IGRAPH_EINVAL);
  }

  if (dom) {
    IGRAPH_CHECK(igraph_vector_resize(dom, no_of_nodes));
  } else {
    mydom=&vdom;
    IGRAPH_VECTOR_INIT_FINALLY(mydom, no_of_nodes);
  }
  igraph_vector_fill(mydom, IGRAPH_NAN);

  IGRAPH_CHECK(igraph_vector_init(&parent, no_of_nodes));
  IGRAPH_FINALLY(igraph_vector_destroy, &parent);
  IGRAPH_CHECK(igraph_vector_long_init(&semi, no_of_nodes));
  IGRAPH_FINALLY(igraph_vector_long_destroy, &semi);
  IGRAPH_CHECK(igraph_vector_init(&vertex, no_of_nodes));
  IGRAPH_FINALLY(igraph_vector_destroy, &vertex);
  IGRAPH_CHECK(igraph_vector_long_init(&ancestor, no_of_nodes));
  IGRAPH_FINALLY(igraph_vector_long_destroy, &ancestor);
  IGRAPH_CHECK(igraph_vector_long_init_seq(&label, 0, no_of_nodes-1));
  IGRAPH_FINALLY(igraph_vector_long_destroy, &label);
  IGRAPH_CHECK(igraph_adjlist_init(graph, &succ, mode));
  IGRAPH_FINALLY(igraph_adjlist_destroy, &succ);
  IGRAPH_CHECK(igraph_adjlist_init(graph, &pred, invmode));
  IGRAPH_FINALLY(igraph_adjlist_destroy, &pred);
  IGRAPH_CHECK(igraph_i_dbucket_init(&bucket, no_of_nodes));
  IGRAPH_FINALLY(igraph_i_dbucket_destroy, &bucket);

  /* DFS first, to set semi, vertex and parent, step 1 */
  
  IGRAPH_CHECK(igraph_dfs(graph, root, mode, /*unreachable=*/ 0,
			  /*order=*/ &vertex,
			  /*order_out=*/ 0, /*father=*/ &parent,
			  /*dist=*/ 0, /*in_callback=*/ 0, 
			  /*out_callback=*/ 0, /*extra=*/ 0));

  for (i=0; i<no_of_nodes; i++) {
    if (IGRAPH_FINITE(VECTOR(vertex)[i])) {
      long int t=VECTOR(vertex)[i];
      VECTOR(semi)[t] = component_size+1;
      VECTOR(vertex)[component_size] = t+1;
      component_size++;
    }
  }
  if (leftout) {
    long int n=no_of_nodes-component_size;
    long int p=0, j;
    IGRAPH_CHECK(igraph_vector_resize(leftout, n));
    for (j=0; j<no_of_nodes && p<n; j++) {
      if (!IGRAPH_FINITE(VECTOR(parent)[j])) {
	VECTOR(*leftout)[p++] = j;
      }
    }
  }

  /* We need to go over 'pred' because it should contain only the
     edges towards the target vertex. */
  for (i=0; i<no_of_nodes; i++) {
    igraph_vector_t *v=igraph_adjlist_get(&pred, i);
    long int j, n=igraph_vector_size(v);
    for (j=0; j<n; ) {
      long int v2=VECTOR(*v)[j];
      if (IGRAPH_FINITE(VECTOR(parent)[v2])) {
	j++; 
      } else {
	VECTOR(*v)[j]=VECTOR(*v)[n-1];
	igraph_vector_pop_back(v);
	n--;
      }
    }
  }

  /* Now comes the main algorithm, steps 2 & 3 */

  for (i=component_size-1; i>0; i--) {
    long int w=VECTOR(vertex)[i]-1;
    igraph_vector_t *predw=igraph_adjlist_get(&pred, w);
    long int j, n=igraph_vector_size(predw);
    for (j=0; j<n; j++) {
      long int v=VECTOR(*predw)[j];
      long int u=igraph_i_dominator_EVAL(v, &ancestor, &label, &semi);
      if (VECTOR(semi)[u] < VECTOR(semi)[w]) {
	VECTOR(semi)[w]=VECTOR(semi)[u];
      }
    }
    igraph_i_dbucket_insert(&bucket, 
			    VECTOR(vertex)[ VECTOR(semi)[w]-1 ]-1, w);
    igraph_i_dominator_LINK(VECTOR(parent)[w], w, &ancestor);
    while (!igraph_i_dbucket_empty(&bucket, VECTOR(parent)[w])) {
      long int v=igraph_i_dbucket_delete(&bucket, VECTOR(parent)[w]);
      long int u=igraph_i_dominator_EVAL(v, &ancestor, &label, &semi);
      VECTOR(*mydom)[v] = VECTOR(semi)[u] < VECTOR(semi)[v] ? u : 
	VECTOR(parent)[w];
    }
  }

  /* Finally, step 4 */

  for (i=1; i<component_size; i++) {
    long int w=VECTOR(vertex)[i]-1;
    if (VECTOR(*mydom)[w] != VECTOR(vertex)[VECTOR(semi)[w]-1]-1) {
      VECTOR(*mydom)[w] = VECTOR(*mydom)[(long int)VECTOR(*mydom)[w]];
    }
  }
  VECTOR(*mydom)[(long int)root]=-1;

  igraph_i_dbucket_destroy(&bucket);
  igraph_adjlist_destroy(&pred);
  igraph_adjlist_destroy(&succ);
  igraph_vector_long_destroy(&label);
  igraph_vector_long_destroy(&ancestor);
  igraph_vector_destroy(&vertex);
  igraph_vector_long_destroy(&semi);
  igraph_vector_destroy(&parent);
  IGRAPH_FINALLY_CLEAN(8);

  if (domtree) {
    igraph_vector_t edges;
    long int ptr=0;
    IGRAPH_VECTOR_INIT_FINALLY(&edges, component_size*2-2);
    for (i=0; i<no_of_nodes; i++) {
      if (i!=root && IGRAPH_FINITE(VECTOR(*mydom)[i])) {
	if (mode==IGRAPH_OUT) {
	  VECTOR(edges)[ptr++] = VECTOR(*mydom)[i];
	  VECTOR(edges)[ptr++] = i;
	} else {
	  VECTOR(edges)[ptr++] = i;
	  VECTOR(edges)[ptr++] = VECTOR(*mydom)[i];
	}
      }
    }
    IGRAPH_CHECK(igraph_create(domtree, &edges, no_of_nodes,
			       IGRAPH_DIRECTED));
    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);
  }

  if (!dom) {
    igraph_vector_destroy(&vdom);
    IGRAPH_FINALLY_CLEAN(1);
  }

  return 0;
}

typedef struct igraph_i_all_st_cuts_minimal_dfs_data_t {
  igraph_stack_t *stack;
  igraph_vector_bool_t *nomark;  
  const igraph_vector_bool_t *GammaX;
  long int root;
  const igraph_vector_t *map;
} igraph_i_all_st_cuts_minimal_dfs_data_t;

igraph_bool_t igraph_i_all_st_cuts_minimal_dfs_incb(const igraph_t *graph,
						    igraph_integer_t vid,
						    igraph_integer_t dist,
						    void *extra) {

  igraph_i_all_st_cuts_minimal_dfs_data_t *data=extra;
  igraph_stack_t *stack=data->stack;
  igraph_vector_bool_t *nomark=data->nomark;
  const igraph_vector_bool_t *GammaX=data->GammaX;
  const igraph_vector_t *map=data->map;
  long int realvid=VECTOR(*map)[(long int)vid];

  if (VECTOR(*GammaX)[(long int)realvid]) {
    if (!igraph_stack_empty(stack)) {
      long int top=igraph_stack_top(stack);
      VECTOR(*nomark)[top]=1;	/* we just found a smaller one */
    }
    igraph_stack_push(stack, realvid); /* TODO: error check */
  }

  return 0;
}

igraph_bool_t igraph_i_all_st_cuts_minimal_dfs_otcb(const igraph_t *graph,
						    igraph_integer_t vid,
						    igraph_integer_t dist,
						    void *extra) {
  igraph_i_all_st_cuts_minimal_dfs_data_t *data=extra;
  igraph_stack_t *stack=data->stack;
  const igraph_vector_t *map=data->map;
  long int realvid=VECTOR(*map)[(long int)vid];

  if (!igraph_stack_empty(stack) && 
      igraph_stack_top(stack) == realvid) {
    igraph_stack_pop(stack);
  }

  return 0;
}

int igraph_i_all_st_cuts_minimal(const igraph_t *graph,
				 const igraph_t *domtree,
				 long int root,
				 const igraph_marked_queue_t *X,
				 const igraph_vector_bool_t *GammaX,
				 const igraph_vector_t *invmap,
				 igraph_vector_t *minimal) {

  long int no_of_nodes=igraph_vcount(graph);
  igraph_stack_t stack;
  igraph_vector_bool_t nomark;
  igraph_i_all_st_cuts_minimal_dfs_data_t data;
  long int i;
  
  IGRAPH_CHECK(igraph_stack_init(&stack, 10));
  IGRAPH_FINALLY(igraph_stack_destroy, &stack);
  IGRAPH_CHECK(igraph_vector_bool_init(&nomark, no_of_nodes));
  IGRAPH_FINALLY(igraph_vector_bool_destroy, &nomark);

  data.stack=&stack;
  data.nomark=&nomark;
  data.GammaX=GammaX;
  data.root=root;
  data.map=invmap;

  /* We mark all GammaX elements as minimal first.
     TODO: actually, we could just use GammaX to return the minimal
     elements. */
  for (i=0; i<no_of_nodes; i++) {
    VECTOR(nomark)[i] = VECTOR(*GammaX)[i] == 0 ? 1 : 0;
  }

  /* We do a reverse DFS from root. If, along a path we find a GammaX
     vertex after (=below) another GammaX vertex, we mark the higher
     one as non-minimal. */

  IGRAPH_CHECK(igraph_dfs(domtree, root, IGRAPH_IN, /*unreachable=*/ 0,
			  /*order=*/ 0,
			  /*order_out=*/ 0, /*father=*/ 0, 
			  /*dist=*/ 0, /*in_callback=*/ 
			  igraph_i_all_st_cuts_minimal_dfs_incb,
			  /*out_callback=*/ 
			  igraph_i_all_st_cuts_minimal_dfs_otcb,
			  /*extra=*/ &data));
  
  igraph_vector_clear(minimal);
  for (i=0; i<no_of_nodes; i++) {
    if (!VECTOR(nomark)[i]) { 
      IGRAPH_CHECK(igraph_vector_push_back(minimal, i));
    }
  }

  igraph_vector_bool_destroy(&nomark);
  igraph_stack_destroy(&stack);
  IGRAPH_FINALLY_CLEAN(2);

  return 0;
}

int igraph_i_all_st_cuts_pivot(const igraph_t *graph,
			       const igraph_marked_queue_t *S,
			       const igraph_stack_t *T,
			       const igraph_vector_bool_t *TV,
			       long int source,
			       long int target,
			       long int *v,
			       igraph_vector_t *Isv) {

  long int no_of_nodes=igraph_vcount(graph);
  igraph_t Sbar;
  igraph_vector_t Sbar_map, Sbar_invmap;
  igraph_vector_t keep;
  igraph_t domtree;
  igraph_vector_t leftout;
  long int i, nomin;
  long int root;
  igraph_vector_t M;
  igraph_vector_bool_t GammaS;
  long int GammaS_len;
  igraph_vector_t Nuv;
  igraph_vector_t Isv_min;
  igraph_vector_t GammaS_vec;

  /* We need to create the graph induced by Sbar */
  IGRAPH_VECTOR_INIT_FINALLY(&Sbar_map, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&Sbar_invmap, 0);

  IGRAPH_VECTOR_INIT_FINALLY(&keep, 0);
  for (i=0; i<no_of_nodes; i++) {
    if (!igraph_marked_queue_iselement(S, i)) {
      IGRAPH_CHECK(igraph_vector_push_back(&keep, i));
    }
  }

  IGRAPH_CHECK(igraph_induced_subgraph_map(graph, &Sbar,
					   igraph_vss_vector(&keep),
					   IGRAPH_SUBGRAPH_AUTO,
					   /* map= */ &Sbar_map, 
					   /* invmap= */ &Sbar_invmap));  
  igraph_vector_destroy(&keep);
  IGRAPH_FINALLY_CLEAN(1);
  IGRAPH_FINALLY(igraph_destroy, &Sbar);

  root=VECTOR(Sbar_map)[target]-1;

  /* -------------------------------------------------------------*/
  /* Construct the dominator tree of Sbar */

  IGRAPH_VECTOR_INIT_FINALLY(&leftout, 0);
  IGRAPH_CHECK(igraph_dominator_tree(&Sbar, root, /*dom=*/ 0, &domtree,
				     &leftout, IGRAPH_IN));
  IGRAPH_FINALLY(igraph_destroy, &domtree);

  /* -------------------------------------------------------------*/
  /* Identify the set M of minimal elements of Gamma(S) with respect
     to the dominator relation. */
  
  /* First we create GammaS */
  /* TODO: use the adjacency list, instead of neighbors() */
  IGRAPH_CHECK(igraph_vector_bool_init(&GammaS, no_of_nodes));
  IGRAPH_FINALLY(igraph_vector_bool_destroy, &GammaS);
  if (igraph_marked_queue_size(S)==0) {
    VECTOR(GammaS)[(long int) VECTOR(Sbar_map)[source]-1]=1;
  } else {
    for (i=0; i<no_of_nodes; i++) {
      if (igraph_marked_queue_iselement(S, i)) {
	igraph_vector_t neis;
	long int j, n;
	IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
	IGRAPH_CHECK(igraph_neighbors(graph, &neis, i, IGRAPH_OUT));
	n=igraph_vector_size(&neis);
	for (j=0; j<n; j++) {
	  long int nei=VECTOR(neis)[j];
	  if (!igraph_marked_queue_iselement(S, nei)) {
	    VECTOR(GammaS)[nei]=1;
	  }
	}
	igraph_vector_destroy(&neis);
	IGRAPH_FINALLY_CLEAN(1);
      }
    }
  }

  IGRAPH_VECTOR_INIT_FINALLY(&M, 0);
  if (igraph_ecount(&domtree)>0) {
    IGRAPH_CHECK(igraph_i_all_st_cuts_minimal(graph, &domtree, root, S,
					      &GammaS, &Sbar_invmap, &M));
  }

  igraph_vector_clear(Isv);
  IGRAPH_VECTOR_INIT_FINALLY(&Nuv, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&Isv_min, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&GammaS_vec, 0);
  for (i=0; i<no_of_nodes; i++) {
    if (VECTOR(GammaS)[i]) {
      IGRAPH_CHECK(igraph_vector_push_back(&GammaS_vec, i));
    }
  }
  GammaS_len=igraph_vector_size(&GammaS_vec);

  nomin=igraph_vector_size(&M);
  for (i=0; i<nomin; i++) {
    /* -------------------------------------------------------------*/
    /* For each v in M find the set Nu(v)=dom(Sbar, v)-K
       Nu(v) contains all vertices that are dominated by v, for every
       v, this is a subtree of the dominator tree, rooted at v. The
       different subtrees are disjoint. */
    long int min=VECTOR(Sbar_map)[(long int) VECTOR(M)[i] ]-1;
    long int nuvsize, isvlen, j;
    IGRAPH_CHECK(igraph_dfs(&domtree, min, IGRAPH_IN, /*unreachable=*/ 0, 
			    /*order=*/ &Nuv, 
			    /*order_out=*/ 0, /*father=*/ 0, /*dist=*/ 0,
			    /*in_callback=*/ 0, /*out_callback=*/ 0, 
			    /*extra=*/ 0));
    /* Remove the NAN values from the end of the vector */
    for (nuvsize=0; nuvsize<no_of_nodes; nuvsize++) {
      igraph_real_t t=VECTOR(Nuv)[nuvsize];
      if (IGRAPH_FINITE(t)) {
	VECTOR(Nuv)[nuvsize]=VECTOR(Sbar_invmap)[(long int) t];
      } else {
	break;
      }
    }
    igraph_vector_resize(&Nuv, nuvsize);
    
    /* -------------------------------------------------------------*/
    /* By a BFS search of <Nu(v)> determine I(S,v)-K.
       I(S,v) contains all vertices that are in Nu(v) and that are
       reachable from Gamma(S) via a path in Nu(v). */
    IGRAPH_CHECK(igraph_bfs(graph, /*root=*/ -1, /*roots=*/ &GammaS_vec,
			    /*mode=*/ IGRAPH_OUT, /*unreachable=*/ 0,
			    /*restricted=*/ &Nuv, 
			    /*order=*/ &Isv_min, /*rank=*/ 0,
			    /*father=*/ 0, /*pred=*/ 0, /*succ=*/ 0,
			    /*dist=*/ 0, /*callback=*/ 0, /*extra=*/ 0));
    for (isvlen=0; isvlen<no_of_nodes; isvlen++) {
      if (!IGRAPH_FINITE(VECTOR(Isv_min)[isvlen])) { break; }
    }
    igraph_vector_resize(&Isv_min, isvlen);
    
    /* -------------------------------------------------------------*/
    /* For each c in M check whether Isv-K is included in Tbar. If 
       such a v is found, compute Isv={x|v[Nu(v) U K]x} and return v and 
       Isv; otherwise return Isv={}. */
    for (j=0; j<isvlen; j++) {
      long int v=VECTOR(Isv_min)[j];
      if (VECTOR(*TV)[v]) { break; }
    }
    /* We might have found one */
    if (j==isvlen) {
      *v=VECTOR(M)[i];
      /* Calculate real Isv */
      IGRAPH_CHECK(igraph_vector_append(&Nuv, &leftout));
      IGRAPH_CHECK(igraph_bfs(graph, /*root=*/ *v, /*roots=*/ 0,
			      /*mode=*/ IGRAPH_OUT, /*unreachable=*/ 0,
			      /*restricted=*/ &Nuv, 
			      /*order=*/ &Isv_min, /*rank=*/ 0,
			      /*father=*/ 0, /*pred=*/ 0, /*succ=*/ 0,
			      /*dist=*/ 0, /*callback=*/ 0, /*extra=*/ 0));
      for (isvlen=0; isvlen<no_of_nodes; isvlen++) {
	if (!IGRAPH_FINITE(VECTOR(Isv_min)[isvlen])) { break; }
      }
      igraph_vector_resize(&Isv_min, isvlen);
      igraph_vector_update(Isv, &Isv_min);

      break;
    }
  }

  igraph_vector_destroy(&GammaS_vec);
  igraph_vector_destroy(&Isv_min);
  igraph_vector_destroy(&Nuv);
  IGRAPH_FINALLY_CLEAN(3);

  igraph_vector_destroy(&M);
  igraph_vector_bool_destroy(&GammaS);
  igraph_destroy(&domtree);
  igraph_vector_destroy(&leftout);  
  igraph_destroy(&Sbar);
  igraph_vector_destroy(&Sbar_map);
  igraph_vector_destroy(&Sbar_invmap);
  IGRAPH_FINALLY_CLEAN(7);

  return 0;
}

/* TODO: This is a temporary resursive version, without proper error
   handling */

int igraph_i_all_st_cuts_list(const igraph_t *graph,
			      igraph_marked_queue_t *S,
			      igraph_stack_t *T,
			      igraph_vector_bool_t *TV, 
			      long int source,
			      long int target,
			      igraph_vector_ptr_t *closed_sets) {

  igraph_vector_t Isv;
  long int v=0;
  long int i, n;

  igraph_vector_init(&Isv, 0);
  
  igraph_i_all_st_cuts_pivot(graph, S, T, TV, source, target, &v, &Isv);
  if (igraph_vector_size(&Isv)==0) {
    if (igraph_marked_queue_size(S) != 0) {
      igraph_vector_t *vec=igraph_Calloc(1, igraph_vector_t);
      igraph_vector_init(vec, igraph_marked_queue_size(S));
      igraph_marked_queue_as_vector(S, vec);
      igraph_vector_ptr_push_back(closed_sets, vec);
    }
  } else {
    /* Put v into T */
    igraph_stack_push(T, v);
    VECTOR(*TV)[v] = 1;

    /* Go down left in the search tree */
    igraph_i_all_st_cuts_list(graph, S, T, TV, source, target, 
			      closed_sets);

    /* Take out v from T */
    igraph_stack_pop(T);
    VECTOR(*TV)[v] = 0;
    
    /* Add Isv to S */
    igraph_marked_queue_start_batch(S);
    n=igraph_vector_size(&Isv);
    for (i=0; i<n; i++) {
      if (!igraph_marked_queue_iselement(S, VECTOR(Isv)[i])) {
	igraph_marked_queue_push(S, VECTOR(Isv)[i]);
      }
    }

    /* Go down right in the search tree */
    
    igraph_i_all_st_cuts_list(graph, S, T, TV, source, target, 
			      closed_sets);

    /* Take out Isv from S */
    igraph_marked_queue_pop_back_batch(S);
  }
    
  igraph_vector_destroy(&Isv);

  return 0;
}

int igraph_all_st_cuts(const igraph_t *graph,
		       igraph_vector_ptr_t *cuts,
		       igraph_vector_ptr_t *partition1s,
		       igraph_integer_t source,
		       igraph_integer_t target) {

  /* S is a special stack, in which elements are pushed in batches. 
     It is then possible to remove the whole batch in one step.

     T & TV together form a stack, that has a fast "is-element" 
     operation. Every element has to be included at most once.
  */

  long int no_of_nodes=igraph_vcount(graph);
  long int no_of_edges=igraph_ecount(graph);
  igraph_marked_queue_t S;
  igraph_stack_t T;
  igraph_vector_bool_t TV;

  if (!igraph_is_directed(graph)) {
    IGRAPH_ERROR("Listing all s-t cuts only implemented for "
		 "directed graphs", IGRAPH_UNIMPLEMENTED);
  }
  
  if (!partition1s) { 
    IGRAPH_ERROR("`partition1s' must not be a null pointer", 
		 IGRAPH_UNIMPLEMENTED);
  }

  IGRAPH_CHECK(igraph_marked_queue_init(&S, no_of_nodes));
  IGRAPH_FINALLY(igraph_marked_queue_destroy, &S);
  IGRAPH_CHECK(igraph_stack_init(&T, 0));
  IGRAPH_FINALLY(igraph_stack_destroy, &T);
  IGRAPH_CHECK(igraph_vector_bool_init(&TV, no_of_nodes));
  IGRAPH_FINALLY(igraph_vector_bool_destroy, &TV);

  if (cuts)        { igraph_vector_ptr_clear(cuts);        }
  if (partition1s) { igraph_vector_ptr_clear(partition1s); }    
  
  /* We call it with S={}, T={} */
  igraph_marked_queue_push(&S, source);
  IGRAPH_CHECK(igraph_i_all_st_cuts_list(graph, &S, &T, &TV, 
					 source, target, partition1s));
  
  if (cuts) {
    igraph_vector_long_t inS;
    long int i, nocuts=igraph_vector_ptr_size(partition1s);
    IGRAPH_CHECK(igraph_vector_long_init(&inS, no_of_nodes));
    IGRAPH_FINALLY(igraph_vector_long_destroy, &inS);
    IGRAPH_CHECK(igraph_vector_ptr_resize(cuts, nocuts));
    for (i=0; i<nocuts; i++) {
      igraph_vector_t *cut;
      igraph_vector_t *part=VECTOR(*partition1s)[i];
      long int cutsize=0;
      long int j, partlen=igraph_vector_size(part);
      /* Mark elements */
      for (j=0; j<partlen; j++) {
	long int v=VECTOR(*part)[j];
	VECTOR(inS)[v] = i+1;
      }
      /* Check how many edges */
      for (j=0; j<no_of_edges; j++) {
	long int from=IGRAPH_FROM(graph, j);
	long int to=IGRAPH_TO(graph, j);
	long int pfrom=VECTOR(inS)[from];
	long int pto=VECTOR(inS)[to];
	if ((pfrom == i+1 && pto != i+1) ||
	    (pfrom != i+1 && pto == i+1)) {
	  cutsize++;
	}
      }
      /* Add the edges */
      cut=igraph_Calloc(1, igraph_vector_t);
      if (!cut) {
	IGRAPH_ERROR("Cannot calculate s-t cuts", IGRAPH_ENOMEM);
      }
      IGRAPH_VECTOR_INIT_FINALLY(cut, cutsize);
      cutsize=0;
      for (j=0; j<no_of_edges; j++) {
	long int from=IGRAPH_FROM(graph, j);
	long int to=IGRAPH_TO(graph, j);
	long int pfrom=VECTOR(inS)[from];
	long int pto=VECTOR(inS)[to];
	if ((pfrom == i+1 && pto != i+1) ||
	    (pfrom != i+1 && pto == i+1)) {
	  VECTOR(*cut)[cutsize++]=j;
	}
      }
      VECTOR(*cuts)[i]=cut;
      IGRAPH_FINALLY_CLEAN(1);
    }
    
    igraph_vector_long_destroy(&inS);
    IGRAPH_FINALLY_CLEAN(1);
  }

  igraph_vector_bool_destroy(&TV);
  igraph_stack_destroy(&T);
  igraph_marked_queue_destroy(&S);
  IGRAPH_FINALLY_CLEAN(3);

  return 0;
}

int igraph_all_st_mincuts(const igraph_t *graph, igraph_real_t *value,
			  igraph_vector_ptr_t *cuts,
			  igraph_vector_ptr_t *partition1s,
			  igraph_integer_t source,
			  igraph_integer_t target,
			  const igraph_vector_t *capacity) {
  /* TODO */
  return 0;
}
