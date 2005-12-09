/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2003, 2004, 2005  Gabor Csardi <csardi@rmki.kfki.hu>
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
#include "memory.h"

int igraph_clusters_weak(igraph_t *graph, vector_t *membership,
			 vector_t *csize);

int igraph_clusters_strong(igraph_t *graph, vector_t *membership,
			   vector_t *csize);

/**
 * \ingroup structural
 * \brief Calculates the (weakly or strongly) connected components in
 * a graph. 
 *
 * @param graph The graph object to analyze.
 * @param membership First half of the result will be stored here. For
 *        every vertex the id of its component is given. The vector
 *        has to be preinitialized and will be resized.
 * @param csize The second half of the result. For every component it
 *        gives its size, the order is defined by the component ids.
 *        The vector has to be preinitialized and will be resized.
 * @param mode For directed graph this specifies whether to calculate
 *        weakly or strongly connected components. Possible values: 
 *        <b>IGRAPH_WEAK</B>, <b>IGRAPH_STRONG</b>. This argument is
 *        igrored for undirected graphs.
 * @return Error code:
 *         - <b>IGRAPH_EINVAL</b>: invalid mode argument.
 * 
 * Time complexity: <code>O(|V|+|E|)</code>, <code>|V|</code> and
 * <code>|E|</code> are the number of vertices and edges in the graph.
 */

#include <string.h>

int igraph_clusters(igraph_t *graph, vector_t *membership, vector_t *csize, 
		    igraph_connectedness_t mode) {
  if (mode==IGRAPH_WEAK || !igraph_is_directed(graph)) {
    return igraph_clusters_weak(graph, membership, csize);
  } else if (mode==IGRAPH_STRONG) {
    return igraph_clusters_strong(graph, membership, csize);
  } else {
    IGRAPH_FERROR("", IGRAPH_EINVAL);
  }
  
  return 1;
}

int igraph_clusters_weak(igraph_t *graph, vector_t *membership,
			 vector_t *csize) {

  long int no_of_nodes=igraph_vcount(graph);
  char *already_added;
  long int first_node, act_cluster_size=0, no_of_clusters=1;
  
  dqueue_t q=DQUEUE_NULL;
  
  long int i;
  vector_t neis=VECTOR_NULL;

  int ret1, ret2;
  
  already_added=Calloc(no_of_nodes,char);
  if (already_added==0) {
    IGRAPH_FERROR("", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(free, already_added); /* TODO: hack */

  DQUEUE_INIT_FINALLY(&q, no_of_nodes > 100000 ? 10000 : no_of_nodes/10);
  VECTOR_INIT_FINALLY(&neis, 0);

  /* Memory for result, csize is dynamically allocated */
  IGRAPH_CHECK(vector_resize(membership, no_of_nodes));
  vector_clear(csize);

  /* The algorithm */

  for (first_node=0; first_node < no_of_nodes; ++first_node) {
    if (already_added[first_node]==1) continue;

    already_added[first_node]=1;
    act_cluster_size=1;
    VECTOR(*membership)[first_node]=no_of_clusters-1;
    IGRAPH_CHECK(dqueue_push(&q, first_node));
    
    while ( !dqueue_empty(&q) ) {
      long int act_node=dqueue_pop(&q);
      IGRAPH_CHECK(igraph_neighbors(graph, &neis, act_node, IGRAPH_ALL));
      for (i=0; i<vector_size(&neis); i++) {
	long int neighbor=VECTOR(neis)[i];
	if (already_added[neighbor]==1) { continue; }
	IGRAPH_CHECK(dqueue_push(&q, neighbor));
	already_added[neighbor]=1;
	act_cluster_size++;
	VECTOR(*membership)[neighbor]=no_of_clusters-1;
      }
    }
    no_of_clusters++;
    IGRAPH_CHECK(vector_push_back(csize, act_cluster_size));
  }
  
  /* Cleaning up */
  
  Free(already_added);
  dqueue_destroy(&q);
  vector_destroy(&neis);
  IGRAPH_FINALLY_CLEAN(3);
  
  return 0;
}

int igraph_clusters_strong(igraph_t *graph, vector_t *membership,
			   vector_t *csize) {

  long int no_of_nodes=igraph_vcount(graph);
  vector_t next_nei=VECTOR_NULL;
  
  long int i;
  dqueue_t q=DQUEUE_NULL;
  
  long int no_of_clusters=1;
  long int act_cluster_size;

  vector_t out=VECTOR_NULL;
  vector_t tmp=VECTOR_NULL;

  int ret1, ret2, ret3;

  /* The result */

  VECTOR_INIT_FINALLY(&next_nei, no_of_nodes);
  VECTOR_INIT_FINALLY(&out, 0);
  DQUEUE_INIT_FINALLY(&q, 100);
  VECTOR_INIT_FINALLY(&tmp, 0);

  IGRAPH_CHECK(vector_resize(membership, no_of_nodes));
  IGRAPH_CHECK(vector_reserve(&out, no_of_nodes));

  vector_null(&out);
  vector_clear(csize);
  
  for (i=0; i<no_of_nodes; i++) {
    IGRAPH_CHECK(igraph_neighbors(graph, &tmp, i, IGRAPH_OUT));
    if (VECTOR(next_nei)[i] > vector_size(&tmp)) { continue; }
    
    IGRAPH_CHECK(dqueue_push(&q, i));
    while (!dqueue_empty(&q)) {
      long int act_node=dqueue_back(&q);
      IGRAPH_CHECK(igraph_neighbors(graph, &tmp, act_node, IGRAPH_OUT));
      if (VECTOR(next_nei)[act_node]==0) {
	/* this is the first time we've met this vertex */
	VECTOR(next_nei)[act_node]++;
      } else if (VECTOR(next_nei)[act_node] <= vector_size(&tmp)) {
	/* we've already met this vertex but it has more children */
	long int neighbor=VECTOR(tmp)[(long int)VECTOR(next_nei)[act_node]-1];
	if (VECTOR(next_nei)[neighbor] == 0) {
	  IGRAPH_CHECK(dqueue_push(&q, neighbor));
	}
	VECTOR(next_nei)[act_node]++;
      } else {
	/* we've met this vertex and it has no more children */
	IGRAPH_CHECK(vector_push_back(&out, act_node));
	dqueue_pop_back(&q);
      }
    } /* while q */
  }  /* for */

  /* OK, we've the 'out' values for the nodes, let's use them in
     descreasing order with the help of a heap */

  vector_null(&next_nei);                            /* mark already
							added vertices */
  while (!vector_empty(&out)) {
    long int grandfather=vector_pop_back(&out);
    if (VECTOR(next_nei)[grandfather] != 0) { continue; }
    VECTOR(next_nei)[grandfather]=1;
    act_cluster_size=1;
    VECTOR(*membership)[grandfather]=no_of_clusters-1;
    IGRAPH_CHECK(dqueue_push(&q, grandfather));
    
    while (!dqueue_empty(&q)) {
      long int act_node=dqueue_pop_back(&q);
      IGRAPH_CHECK(igraph_neighbors(graph, &tmp, act_node, IGRAPH_IN));
      for (i=0; i<vector_size(&tmp); i++) {
	long int neighbor=VECTOR(tmp)[i];
	if (VECTOR(next_nei)[neighbor] != 0) { continue; }
	IGRAPH_CHECK(dqueue_push(&q, neighbor));
	VECTOR(next_nei)[neighbor]=1;
	act_cluster_size++;
	VECTOR(*membership)[neighbor]=no_of_clusters-1;
      }
    }
    no_of_clusters++;
    IGRAPH_CHECK(vector_push_back(csize, act_cluster_size));
  }
  
  /* Clean up, return */

  vector_destroy(&out);
  vector_destroy(&tmp);
  dqueue_destroy(&q);
  vector_destroy(&next_nei);
  IGRAPH_FINALLY_CLEAN(4);

  return 0;
}

int igraph_is_connected_weak(igraph_t *graph, bool_t *res);

/**
 * \ingroup structural
 * \brief Decides whether the graph is (weakly or strongly) connected.
 * 
 * @param graph The graph object to analyze.
 * @param res Pointer to a logical variable, the result will be stored
 *        here. 
 * @param mode For directed graph this specifies whether to calculate
 *        weak or strong connectedness. Possible values: 
 *        <b>IGRAPH_WEAK</B>, <b>IGRAPH_STRONG</b>. This argument is
 *        igrored for undirected graphs.
 * @return Error code:
 *         - <b>IGRAPH_EINVAL</b>: invalid mode argument.
 *
 * Time complexity: <code>O(|V|+|E|)</code>, the number of vertices
 * plus the number of edges in the graph.
 */


int igraph_is_connected(igraph_t *graph, bool_t *res, 
			igraph_connectedness_t mode) {
  if (mode==IGRAPH_WEAK || !igraph_is_directed(graph)) {
    return igraph_is_connected_weak(graph, res);
  } else if (mode==IGRAPH_STRONG) {
    vector_t membership=VECTOR_NULL;
    vector_t csize=VECTOR_NULL;
    int retval;
    VECTOR_INIT_FINALLY(&membership, 0);
    VECTOR_INIT_FINALLY(&csize, 0);
    retval = igraph_clusters_strong(graph, &membership, &csize);
    *res = (vector_size(&csize)==1);
    vector_destroy(&membership);
    vector_destroy(&csize);
    IGRAPH_FINALLY_CLEAN(2);
    return retval;
  } else {
    IGRAPH_FERROR("mode argument", IGRAPH_EINVAL);
  }
  return 0;
}

int igraph_is_connected_weak(igraph_t *graph, bool_t *res) {

  long int no_of_nodes=igraph_vcount(graph);
  char *already_added;
  vector_t neis=VECTOR_NULL;
  dqueue_t q=DQUEUE_NULL;
  
  long int i, j;

  int ret1, ret2;

  already_added=Calloc(no_of_nodes, char);
  if (already_added==0) {
    IGRAPH_FERROR("is connected (weak) failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(free, already_added); /* TODO: hack */

  DQUEUE_INIT_FINALLY(&q, 10);
  VECTOR_INIT_FINALLY(&neis, 0);
  
  /* Try to find at least two clusters */
  already_added[0]=1;
  IGRAPH_CHECK(dqueue_push(&q, 0));
  
  j=1;
  while ( !dqueue_empty(&q)) {
    long int actnode=dqueue_pop(&q);
    IGRAPH_CHECK(igraph_neighbors(graph, &neis, actnode, IGRAPH_ALL));
    for (i=0; i <vector_size(&neis); i++) {
      long int neighbor=VECTOR(neis)[i];
      if (already_added[neighbor] != 0) { continue; }
      IGRAPH_CHECK(dqueue_push(&q, neighbor));
      j++;
      already_added[neighbor]++;
    }
  }
  
  /* Connected? */
  *res = (j == no_of_nodes);

  Free(already_added);
  dqueue_destroy(&q);
  vector_destroy(&neis);
  IGRAPH_FINALLY_CLEAN(3);

  return 0;
}
