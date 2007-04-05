/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2007  Gabor Csardi <csardi@rmki.kfki.hu>
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
#include "memory.h"

#include <string.h>

/**
 * \function igraph_community_eb_get_merges
 * \brief Calculating the merges, ie. the dendrogram for an edge betweenness
 * community structure
 * 
 * </para><para> 
 * This function is handy if you have a sequence of edge which are
 * gradually removed from the network and you would like to know how
 * the network falls apart into separate components. The edge sequence
 * may come from the \ref igraph_community_edge_betweenness()
 * function, but this is not neccessary. Note that \ref
 * igraph_community_edge_betweenness can also calculate the
 * dendrogram, via its \p merges argument.
 *
 * \param graph The input graph.
 * \param edges Vector containing the edges to be removed from the
 *    network, all edges are expected to appear exactly once in the
 *    vector.
 * \param res Pointer to an initialized matrix, if not NULL then the 
 *    dendrogram will be stored here, in the same form as for the \ref
 *    igraph_community_walktrap() function: the matrix has two columns
 *    and each line is a merge given by the ids of the merged
 *    components. The component ids are number from zero and
 *    component ids smaller than the number of vertices in the graph
 *    belong to individual vertices. The non-trivial components
 *    containing at least two vertices are numbered from \c n, \c n is
 *    the number of vertices in the graph. So if the first line
 *    contains \c a and \c b that means that components \c a and \c b
 *    are merged into component \c n, the second line creates
 *    component \c n+1, etc. The matrix will be resized as needed.
 * \param bridges Pointer to an initialized vector or NULL. If not
 *    null then the index of the edge removals which split the network
 *    will be stored here. The vector will be resized as needed.
 * \return Error code.
 * 
 * \sa \ref igraph_community_edge_betweenness().
 * 
 * Time complexity: O(|E|+|V|log|V|), |V| is the number of vertices,
 * |E| is the number of edges.
 */

int igraph_community_eb_get_merges(const igraph_t *graph, 
				   const igraph_vector_t *edges,
				   igraph_matrix_t *res,
				   igraph_vector_t *bridges) {

  long int no_of_nodes=igraph_vcount(graph);
  igraph_vector_t ptr;
  long int i, midx=0;
  
  IGRAPH_VECTOR_INIT_FINALLY(&ptr, no_of_nodes*2-1);
  if (res) { 
    IGRAPH_CHECK(igraph_matrix_resize(res, no_of_nodes-1, 2));
  }
  if (bridges) {
    IGRAPH_CHECK(igraph_vector_resize(bridges, no_of_nodes-1));
  }
  
  for (i=igraph_vector_size(edges)-1; i>=0; i--) {
    long int edge=VECTOR(*edges)[i];
    igraph_integer_t from, to;
    long int c1, c2, idx;
    igraph_edge(graph, edge, &from, &to);
    idx=from+1;
    while (VECTOR(ptr)[idx-1] != 0) {
      idx=VECTOR(ptr)[idx-1];
    }
    c1=idx-1;
    idx=to+1;
    while (VECTOR(ptr)[idx-1] != 0) {
      idx=VECTOR(ptr)[idx-1];
    }
    c2=idx-1;
    if (c1 != c2) {		/* this is a merge */
      if (res) {
	MATRIX(*res, midx, 0)=c1;
	MATRIX(*res, midx, 1)=c2;
      }
      if (bridges) {
	VECTOR(*bridges)[midx]=i+1;
      }
      
      VECTOR(ptr)[c1]=no_of_nodes+midx+1;
      VECTOR(ptr)[c2]=no_of_nodes+midx+1;
      VECTOR(ptr)[(long int)from]=no_of_nodes+midx+1;
      VECTOR(ptr)[(long int)to]=no_of_nodes+midx+1;
      
      midx++;
    }
  }

  igraph_vector_destroy(&ptr);
  IGRAPH_FINALLY_CLEAN(1);
  
  return 0;
}

/* Find the smallest active element in the vector */
long int igraph_i_vector_which_max_not_null(const igraph_vector_t *v, 
					    const char *passive) {
  long int which, i=0, size=igraph_vector_size(v);
  igraph_real_t max;
  while (passive[i]) {
    i++;
  }
  which=i;
  max=VECTOR(*v)[which];
  for (i++; i<size; i++) {
    igraph_real_t elem=VECTOR(*v)[i];
    if (!passive[i] && elem > max) {
      max=elem;
      which=i;
    }
  }
  
  return which;
}

/**
 * \function igraph_community_edge_betweenness
 * 
 * Community structure detection based on the betweenness of the edges
 * in the network. The algorithm was invented by M. Girvan ans
 * M. Newman, see TODO: citation.
 * 
 * </para><para>
 * The idea is that the betweenness of the edges connecting two
 * communities is typically high, as many of the shortest paths
 * between nodes in separate communities go through them. So we
 * gradually remove the edge with highest betweenness from the
 * network, and recalculate edge betweenness after every removal. 
 * This way sooner or later the network falls off to two components,
 * then after a while one of these components falls off to two smaller 
 * components, etc. until all edges are removed. This is a divisive
 * hieararchical approach, the result is a dendrogram.
 * \param graph The input graph.
 * \param result Pointer to an initialized vector, the result will be
 *     stored here, the ids of the removed edges in the order of their 
 *     removal. It will be resized as needed.
 * \param edge_betweenness Pointer to an initialized vector or
 *     NULL. In the former case the edge betweenness of the removed
 *     edge is stored here. The vector will be resized as needed.
 * \param merges Pointer to an initialized matrix or NULL. If not NULL
 *     then merges performed by the algorithm are stored here. Even if
 *     this is a divisive algorithm, we can replay it backwards and
 *     note which two clusters were merged. Clusters are numbered from
 *     zero, see the \p merges argument of \ref
 *     igraph_community_walktrap() for details. The matrix will be
 *     resized as needed.
 * \param bridges Pointer to an initialized vector of NULL. If not
 *     NULL then all edge removals which separated the network into
 *     more components are marked here.
 * \param directed Logical constant, whether to calculate directed
 *    betweenness (ie. directed paths) for directed graphs. It is
 *    ignored for undirected graphs.
 * \return Error code.
 * 
 * \sa \ref igraph_community_eb_get_merges(), \ref
 * igraph_community_spinglass(), \ref igraph_community_walktrap().
 * 
 * Time complexity: O(|V|^3), as the betweenness calculation requires
 * O(|V|^2) and we do it |V|-1 times.
 */
  
int igraph_community_edge_betweenness(const igraph_t *graph, 
				      igraph_vector_t *result,
				      igraph_vector_t *edge_betweenness,
				      igraph_matrix_t *merges,
				      igraph_vector_t *bridges,
				      igraph_bool_t directed) {
  
  long int no_of_nodes=igraph_vcount(graph);
  long int no_of_edges=igraph_ecount(graph);
  igraph_dqueue_t q=IGRAPH_DQUEUE_NULL;
  long int *distance, *nrgeo;
  double *tmpscore;
  igraph_stack_t stack=IGRAPH_STACK_NULL;
  long int source, i, e;
  
  igraph_i_adjedgelist_t elist_out, elist_in;
  igraph_i_adjedgelist_t *elist_out_p, *elist_in_p;
  igraph_vector_t *neip;
  long int neino;
  igraph_integer_t modein, modeout;
  igraph_vector_t eb;
  long int maxedge, pos;
  igraph_integer_t from, to;

  char *passive;

  directed=directed && igraph_is_directed(graph);
  if (directed) {
    modeout=IGRAPH_OUT;
    modeout=IGRAPH_IN;
    IGRAPH_CHECK(igraph_i_adjedgelist_init(graph, &elist_out, IGRAPH_OUT));
    IGRAPH_FINALLY(igraph_i_adjedgelist_destroy, &elist_out);
    IGRAPH_CHECK(igraph_i_adjedgelist_init(graph, &elist_in, IGRAPH_IN));
    IGRAPH_FINALLY(igraph_i_adjedgelist_destroy, &elist_in);
    elist_out_p=&elist_out;
    elist_in_p=&elist_in;
  } else {
    modeout=modein=IGRAPH_ALL;
    IGRAPH_CHECK(igraph_i_adjedgelist_init(graph, &elist_out, IGRAPH_ALL));
    IGRAPH_FINALLY(igraph_i_adjedgelist_destroy, &elist_out);
    elist_out_p=elist_in_p=&elist_out;
  }
  
  distance=Calloc(no_of_nodes, long int);
  if (distance==0) {
    IGRAPH_ERROR("edge betweenness community struture failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, distance);
  nrgeo=Calloc(no_of_nodes, long int);
  if (nrgeo==0) {
    IGRAPH_ERROR("edge betweenness community struture failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, nrgeo);
  tmpscore=Calloc(no_of_nodes, double);
  if (tmpscore==0) {
    IGRAPH_ERROR("edge betweenness community struture failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, tmpscore);

  IGRAPH_DQUEUE_INIT_FINALLY(&q, 100);
  IGRAPH_CHECK(igraph_stack_init(&stack, no_of_nodes));
  IGRAPH_FINALLY(igraph_stack_destroy, &stack);
  
  IGRAPH_CHECK(igraph_vector_resize(result, no_of_edges));
  if (edge_betweenness) {
    IGRAPH_CHECK(igraph_vector_resize(edge_betweenness, no_of_edges));
    VECTOR(*edge_betweenness)[no_of_edges-1]=0;
  }

  IGRAPH_VECTOR_INIT_FINALLY(&eb, no_of_edges);
  
  passive=Calloc(no_of_edges, char);
  if (!passive) {
    IGRAPH_ERROR("edge betweenness community struture failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, passive);

  for (e=0; e<no_of_edges; e++) {
    
    igraph_vector_null(&eb);

    for (source=0; source<no_of_nodes; source++) {

      /* This will contain the edge betweenness in the current step */
      IGRAPH_ALLOW_INTERRUPTION();

      memset(distance, 0, no_of_nodes*sizeof(long int));
      memset(nrgeo, 0, no_of_nodes*sizeof(long int));
      memset(tmpscore, 0, no_of_nodes*sizeof(double));
      igraph_stack_clear(&stack); /* it should be empty anyway... */
      
      IGRAPH_CHECK(igraph_dqueue_push(&q, source));
      
      nrgeo[source]=1;
      distance[source]=0;
      
      while (!igraph_dqueue_empty(&q)) {
	long int actnode=igraph_dqueue_pop(&q);
	
	neip=igraph_i_adjedgelist_get(elist_out_p, actnode);
	neino=igraph_vector_size(neip);
	for (i=0; i<neino; i++) {
	  igraph_integer_t edge=VECTOR(*neip)[i], from, to;
	  long int neighbor;
	  igraph_edge(graph, edge, &from, &to);
	  neighbor = actnode!=from ? from : to;
	  if (nrgeo[neighbor] != 0) {
	    /* we've already seen this node, another shortest path? */
	    if (distance[neighbor]==distance[actnode]+1) {
	      nrgeo[neighbor]+=nrgeo[actnode];
	    }
	  } else {
	    /* we haven't seen this node yet */
	    nrgeo[neighbor]+=nrgeo[actnode];
	    distance[neighbor]=distance[actnode]+1;
	    IGRAPH_CHECK(igraph_dqueue_push(&q, neighbor));
	    IGRAPH_CHECK(igraph_stack_push(&stack, neighbor));
	  }
	}
      } /* while !igraph_dqueue_empty */
      
      /* Ok, we've the distance of each node and also the number of
	 shortest paths to them. Now we do an inverse search, starting
	 with the farthest nodes. */
      while (!igraph_stack_empty(&stack)) {
	long int actnode=igraph_stack_pop(&stack);
	if (distance[actnode]<1) { continue; } /* skip source node */
	
	/* set the temporary score of the friends */
	neip=igraph_i_adjedgelist_get(elist_in_p, actnode);
	neino=igraph_vector_size(neip);
	for (i=0; i<neino; i++) {
	  igraph_integer_t from, to;
	  long int neighbor;
	  long int edgeno=VECTOR(*neip)[i];
	  igraph_edge(graph, edgeno, &from, &to);
	  neighbor= actnode != from ? from : to;
	  if (distance[neighbor]==distance[actnode]-1 &&
	      nrgeo[neighbor] != 0) {
	    tmpscore[neighbor] +=
	      (tmpscore[actnode]+1)*nrgeo[neighbor]/nrgeo[actnode];
	    VECTOR(eb)[edgeno] +=
	      (tmpscore[actnode]+1)*nrgeo[neighbor]/nrgeo[actnode];
	  }
	}
      }
      /* Ok, we've the scores for this source */
    } /* for source <= no_of_nodes */
    
    /* Now look for the smallest edge betweenness */
    /* and eliminate that edge from the network */
    maxedge=igraph_i_vector_which_max_not_null(&eb, passive);
    VECTOR(*result)[e]=maxedge;
    if (edge_betweenness) {
      VECTOR(*edge_betweenness)[e]=VECTOR(eb)[maxedge];
      if (!directed) { 
	VECTOR(*edge_betweenness)[e] /= 2.0;
      }
    }
    passive[maxedge]=1;
    igraph_edge(graph, maxedge, &from, &to);

    neip=igraph_i_adjedgelist_get(elist_in_p, to);
    neino=igraph_vector_size(neip);
    igraph_vector_search(neip, 0, maxedge, &pos);
    VECTOR(*neip)[pos]=VECTOR(*neip)[neino-1];
    igraph_vector_pop_back(neip);
    
    neip=igraph_i_adjedgelist_get(elist_out_p, from);
    neino=igraph_vector_size(neip);
    igraph_vector_search(neip, 0, maxedge, &pos);
    VECTOR(*neip)[pos]=VECTOR(*neip)[neino-1];
    igraph_vector_pop_back(neip);
  }

  igraph_free(passive);
  igraph_vector_destroy(&eb);
  igraph_stack_destroy(&stack);
  igraph_dqueue_destroy(&q);
  igraph_free(tmpscore);
  igraph_free(nrgeo);
  igraph_free(distance);
  IGRAPH_FINALLY_CLEAN(6);
  
  if (directed) {
    igraph_i_adjedgelist_destroy(&elist_out);
    igraph_i_adjedgelist_destroy(&elist_in);
    IGRAPH_FINALLY_CLEAN(2);
  } else {
    igraph_i_adjedgelist_destroy(&elist_out);
    IGRAPH_FINALLY_CLEAN(1);
  }

  if (merges || bridges) {
    IGRAPH_CHECK(igraph_community_eb_get_merges(graph, result, merges, bridges));
  }
  
  return 0;
}
