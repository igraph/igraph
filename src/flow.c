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

  /* TODO: what if not directed? Gives correct result??? */
  /* TODO: what if not connected? Gives correctly zero??? */
  /* TODO: what if not acyclic? Gives correct results??? */
  /* TODO: return the vertices in the cut. */

  long int no_of_nodes=igraph_vcount(graph);
  long int no_of_edges=igraph_ecount(graph);

  igraph_vector_t flow;
  igraph_vector_t current;
  igraph_dqueue_t queue;
  igraph_vector_t distance;
  igraph_vector_t excess;
  
  igraph_i_adjlist_t out_adjlist;
  igraph_i_adjedgelist_t out_adjedgelist;

  igraph_vector_t *neis, *neis2;

  long int i;

  if (igraph_vector_size(capacity) != no_of_edges) {
    IGRAPH_ERROR("Invalid length of capacity vector", IGRAPH_EINVAL);
  }

  IGRAPH_VECTOR_INIT_FINALLY(&flow, no_of_edges);
  IGRAPH_VECTOR_INIT_FINALLY(&current, no_of_nodes);
  IGRAPH_DQUEUE_INIT_FINALLY(&queue, 100);
  IGRAPH_VECTOR_INIT_FINALLY(&distance, no_of_nodes);
  IGRAPH_VECTOR_INIT_FINALLY(&excess, no_of_nodes);
  
  IGRAPH_CHECK(igraph_i_adjlist_init(graph, &out_adjlist, IGRAPH_OUT));
  IGRAPH_FINALLY(igraph_i_adjlist_destroy, &out_adjlist);

  IGRAPH_CHECK(igraph_i_adjedgelist_init(graph, &out_adjedgelist, IGRAPH_OUT));
  IGRAPH_FINALLY(igraph_i_adjedgelist_destroy, &out_adjedgelist);  

  /* Initialization */
  /* flowvalue[] == 0 initialized already */
  neis=igraph_i_adjedgelist_get(&out_adjedgelist, source);
  neis2=igraph_i_adjlist_get(&out_adjlist, source);
  for (i=0; i<igraph_vector_size(neis); i++) {
    long int edge=VECTOR(*neis)[i];
    long int nei=VECTOR(*neis2)[i];
    VECTOR(flow)[edge] = VECTOR(*capacity)[edge];
    VECTOR(excess)[nei] = VECTOR(*capacity)[edge];
    if (nei != target) {
      IGRAPH_CHECK(igraph_dqueue_push(&queue, nei));
    }
  }
  VECTOR(distance)[ (long int) source ] = no_of_nodes;  

  /* The main part comes here */
  while (!igraph_dqueue_empty(&queue)) {
    long int vertex=igraph_dqueue_pop(&queue);
    igraph_real_t origdist=VECTOR(distance)[vertex];
    do {
      /* PUSH/RELABEL */
      long int cedge, nei;
      igraph_real_t res;
      neis=igraph_i_adjedgelist_get(&out_adjedgelist, vertex);
      neis2=igraph_i_adjlist_get(&out_adjlist, vertex);
      cedge=VECTOR(*neis)[ (long int) VECTOR(current)[vertex] ];
      nei=VECTOR(*neis2)[ (long int) VECTOR(current)[vertex] ];
      res=VECTOR(*capacity)[cedge]-VECTOR(flow)[cedge];
      /* Is the push applicable? */
      if (res > 0.0 && VECTOR(distance)[vertex] == VECTOR(distance)[nei]) {
	/* PUSH */
	igraph_real_t delta=
	  VECTOR(excess)[vertex] < res ? VECTOR(excess)[vertex] : res;
	VECTOR(flow)[cedge] += delta;
	VECTOR(excess)[vertex] -= delta;
	VECTOR(excess)[nei] += delta;
	if (VECTOR(excess)[nei] == delta && nei != target) {
	  /* becomes active */
	  IGRAPH_CHECK(igraph_dqueue_push(&queue, nei));
	}
      } else {
	if (VECTOR(current)[vertex] < igraph_vector_size(neis)-1) {
	  VECTOR(current)[vertex] += 1;
	} else {
	  igraph_bool_t l=0;
	  igraph_real_t min;
	  VECTOR(current)[vertex] = 0;
	  /* RELABEL */	  
	  for (i=0; i<igraph_vector_size(neis2); i++) {
	    long int ee=VECTOR(*neis)[i];
	    long int nn=VECTOR(*neis2)[i];
	    if (VECTOR(flow)[ee] < VECTOR(*capacity)[ee]) {
	      if (!l) {
		min=VECTOR(distance)[nn];
		l=1;
	      } else if (min > VECTOR(distance)[nn]) {
		min=VECTOR(distance)[nn];
	      }
	    }
	  } /* for (neis2) */
	  if (!l) { 
	    VECTOR(distance)[vertex] = 1.0/0.0;	/* infinity */
	  } else {
	    VECTOR(distance)[vertex] = min+1;
	  }
	}
      }
    } while (VECTOR(excess)[vertex] != 0.0 && 
	     VECTOR(distance)[vertex] == origdist);
  }

  /* Store the result */
  if (value) {
    *value=VECTOR(excess)[(long int)target];
  }

  /* Clean up everything */
  igraph_i_adjedgelist_destroy(&out_adjedgelist);
  igraph_i_adjlist_destroy(&out_adjlist);
  
  igraph_vector_destroy(&excess);
  igraph_vector_destroy(&distance);
  igraph_dqueue_destroy(&queue);
  igraph_vector_destroy(&current);
  igraph_vector_destroy(&flow);
  IGRAPH_FINALLY_CLEAN(7);
  return 0;
}
