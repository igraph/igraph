/* -*- mode: C -*-  */
/* vim:set ts=2 sw=2 sts=2 et: */
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
#include "random.h"
#include "arpack.h"
#include "arpack_internal.h"
#include "config.h"

#include <string.h>
#include <math.h>

/**
 * \function igraph_community_eb_get_merges
 * \brief Calculating the merges, ie. the dendrogram for an edge betweenness community structure
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
  igraph_integer_t no_comps;
  
  IGRAPH_CHECK(igraph_clusters(graph, 0, 0, &no_comps, IGRAPH_WEAK));
  
  IGRAPH_VECTOR_INIT_FINALLY(&ptr, no_of_nodes*2-1);
  if (res) { 
    IGRAPH_CHECK(igraph_matrix_resize(res, no_of_nodes-no_comps, 2));
  }
  if (bridges) {
    IGRAPH_CHECK(igraph_vector_resize(bridges, no_of_nodes-no_comps));
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
 * \brief Community findinf based on edge betweenness
 * 
 * Community structure detection based on the betweenness of the edges
 * in the network. The algorithm was invented by M. Girvan and
 * M. Newman, see: M. Girvan and M. E. J. Newman: Community structure in
 * social and biological networks, Proc. Nat. Acad. Sci. USA 99, 7821-7826
 * (2002).
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
  
  igraph_adjedgelist_t elist_out, elist_in;
  igraph_adjedgelist_t *elist_out_p, *elist_in_p;
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
    IGRAPH_CHECK(igraph_adjedgelist_init(graph, &elist_out, IGRAPH_OUT));
    IGRAPH_FINALLY(igraph_adjedgelist_destroy, &elist_out);
    IGRAPH_CHECK(igraph_adjedgelist_init(graph, &elist_in, IGRAPH_IN));
    IGRAPH_FINALLY(igraph_adjedgelist_destroy, &elist_in);
    elist_out_p=&elist_out;
    elist_in_p=&elist_in;
  } else {
    modeout=modein=IGRAPH_ALL;
    IGRAPH_CHECK(igraph_adjedgelist_init(graph, &elist_out, IGRAPH_ALL));
    IGRAPH_FINALLY(igraph_adjedgelist_destroy, &elist_out);
    elist_out_p=elist_in_p=&elist_out;
  }
  
  distance=igraph_Calloc(no_of_nodes, long int);
  if (distance==0) {
    IGRAPH_ERROR("edge betweenness community structure failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, distance);
  nrgeo=igraph_Calloc(no_of_nodes, long int);
  if (nrgeo==0) {
    IGRAPH_ERROR("edge betweenness community structure failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, nrgeo);
  tmpscore=igraph_Calloc(no_of_nodes, double);
  if (tmpscore==0) {
    IGRAPH_ERROR("edge betweenness community structure failed", IGRAPH_ENOMEM);
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
  
  passive=igraph_Calloc(no_of_edges, char);
  if (!passive) {
    IGRAPH_ERROR("edge betweenness community structure failed", IGRAPH_ENOMEM);
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
	
	neip=igraph_adjedgelist_get(elist_out_p, actnode);
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
	neip=igraph_adjedgelist_get(elist_in_p, actnode);
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

    neip=igraph_adjedgelist_get(elist_in_p, to);
    neino=igraph_vector_size(neip);
    igraph_vector_search(neip, 0, maxedge, &pos);
    VECTOR(*neip)[pos]=VECTOR(*neip)[neino-1];
    igraph_vector_pop_back(neip);
    
    neip=igraph_adjedgelist_get(elist_out_p, from);
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
  IGRAPH_FINALLY_CLEAN(7);
  
  if (directed) {
    igraph_adjedgelist_destroy(&elist_out);
    igraph_adjedgelist_destroy(&elist_in);
    IGRAPH_FINALLY_CLEAN(2);
  } else {
    igraph_adjedgelist_destroy(&elist_out);
    IGRAPH_FINALLY_CLEAN(1);
  }

  if (merges || bridges) {
    IGRAPH_CHECK(igraph_community_eb_get_merges(graph, result, merges, bridges));
  }
  
  return 0;
}


/**
 * \function igraph_community_to_membership 
 * \brief Create membership vector from community structure dendrogram
 * 
 * This function creates a membership vector from a community
 * structure dendrogram. A membership vector contains for each vertex
 * the id of its graph component, the graph components are numbered
 * from zero, see the same argument of \ref igraph_clusters() for an
 * example of a membership vector.
 * 
 * </para><para>
 * Many community detection algorithms return with a \em merges
 * matrix, \ref igraph_community_walktrap() an \ref
 * igraph_community_edge_betweenness() are two examples. The matrix
 * contains the merge operations performed while mapping the
 * hierarchical structure of a network. If the matrix has \c n-1 rows,
 * where \c n is the number of vertices in the graph, then it contains
 * the hierarchical structure of the whole network and it is called a
 * dendrogram. 
 * 
 * </para><para>
 * This function performs \p steps merge operations as prescribed by
 * the \p merges matrix and returns the current state of the network.
 * 
 * </para><para>
 * If if \p merges is not a complete dendrogram, it is possible to
 * take \p steps steps if \p steps is not bigger than the number 
 * lines in \p merges.
 * \param merges The two-column matrix containing the merge
 *    operations. See \ref igraph_community_walktrap() for the
 *    detailed syntax.
 * \param nodes The number of leaf nodes in the dendrogram
 * \param steps Integer constant, the number of steps to take.
 * \param membership Pointer to an initialied vector, the membership
 *    results will be stored here, if not NULL. The vector will be
 *    resized as needed.
 * \param csize Pointer to an initialized vector, or NULL. If not NULL
 *    then the sizes of the components will be stored here, the vector
 *    will be resized as needed.
 * 
 * \sa \ref igraph_community_walktrap(), \ref
 * igraph_community_edge_betweenness(), \ref
 * igraph_community_fastgreedy() for community structure detection
 * algorithms.
 * 
 * Time complexity: O(|V|), the number of vertices in the graph.
 */

int igraph_community_to_membership(const igraph_matrix_t *merges,
				   igraph_integer_t nodes,
				   igraph_integer_t steps,
				   igraph_vector_t *membership,
				   igraph_vector_t *csize) {
  
  long int no_of_nodes=nodes;
  long int components=no_of_nodes-steps;
  long int i, found=0;
  igraph_vector_t tmp;
  
  if (steps > igraph_matrix_nrow(merges)) {
    IGRAPH_ERROR("`steps' to big or `merges' matrix too short", IGRAPH_EINVAL);
  }

  if (membership) {
    IGRAPH_CHECK(igraph_vector_resize(membership, no_of_nodes));
    igraph_vector_null(membership);
  }
  if (csize) {
    IGRAPH_CHECK(igraph_vector_resize(csize, components));
    igraph_vector_null(csize);
  }
  
  IGRAPH_VECTOR_INIT_FINALLY(&tmp, steps);
  
  for (i=steps-1; i>=0; i--) {
    long int c1=MATRIX(*merges, i, 0);
    long int c2=MATRIX(*merges, i, 1);

    /* new component? */
    if (VECTOR(tmp)[i]==0) {
      found++;
      VECTOR(tmp)[i]=found;
    }

    if (c1<no_of_nodes) {
      long int cid=VECTOR(tmp)[i]-1;
      if (membership) { VECTOR(*membership)[c1]=cid+1; }
      if (csize) { VECTOR(*csize)[cid] += 1; }
    } else {
      VECTOR(tmp)[c1-no_of_nodes]=VECTOR(tmp)[i];
    }
    
    if (c2<no_of_nodes) { 
      long int cid=VECTOR(tmp)[i]-1;
      if (membership) {	VECTOR(*membership)[c2]=cid+1; }
      if (csize) { VECTOR(*csize)[cid] += 1; }
    } else {
      VECTOR(tmp)[c2-no_of_nodes]=VECTOR(tmp)[i];
    }
    
  }
  
  if (membership || csize) {
    for (i=0; i<no_of_nodes; i++) {
      long int tmp=VECTOR(*membership)[i];
      if (tmp!=0) {
	if (membership) {
	  VECTOR(*membership)[i]=tmp-1;
	}
      } else {
	if (csize) {
	  VECTOR(*csize)[found]+=1;
	}
	if (membership) {
	  VECTOR(*membership)[i]=found;
	}
	found++;
      }
    }
  }
  
  igraph_vector_destroy(&tmp);
  IGRAPH_FINALLY_CLEAN(1);
  
  return 0;
}

/**
 * \function igraph_modularity
 * \brief Calculate the modularity of a graph with respect to some vertex types
 * 
 * The modularity of a graph with respect to some division (or vertex
 * types) measures how good the division is, or how separated are the 
 * different vertex types from each other. It defined as 
 * Q=1/(2m) * sum(Aij-ki*kj/(2m)delta(ci,cj),i,j), here `m' is the
 * number of edges, `Aij' is the element of the `A' adjacency matrix
 * in row `i' and column `j', `ki' is the degree of `i', `kj' is the
 * degree of `j', `ci' is the type (or component) of `i', `cj' that of
 * `j', the sum goes over all `i' and `j' pairs of vertices, and
 * `delta(x,y)' is one if x=y and zero otherwise.
 *
 * </para><para>
 * Modularity on weighted graphs is also meaningful. When taking edge
 * weights into account, `Aij' becomes the weight of the corresponding
 * edge (or 0 if there is no edge), `ki' is the total weight of edges
 * adjacent to vertex `i', `kj' is the total weight of edges adjacent
 * to vertex `j' and `m' is the total weight of all edges.
 * 
 * </para><para>
 * See also MEJ Newman and M Girvan: Finding and evaluating community
 * structure in networks. Physical Review E 69 026113, 2004.
 * \param graph The input graph.
 * \param membership Numeric vector which gives the type of each
 *     vertex, ie. the component to which it belongs.
 * \param modularity Pointer to a real number, the result will be
 *     stored here.
 * \param weights Weight vector or NULL if no weights are specified.
 * \return Error code.
 * 
 * Time complexity: O(|V|+|E|), the number of vertices plus the number
 * of edges.
 */

int igraph_modularity(const igraph_t *graph, 
		      const igraph_vector_t *membership,
		      igraph_real_t *modularity,
			  const igraph_vector_t *weights) {
  
  igraph_vector_t e, a;
  long int types=igraph_vector_max(membership)+1;
  long int no_of_edges=igraph_ecount(graph);
  long int i;
  igraph_integer_t from, to, m;
  long int c1, c2;
  
  IGRAPH_VECTOR_INIT_FINALLY(&e, types);
  IGRAPH_VECTOR_INIT_FINALLY(&a, types);
  
  if (weights) {
    if (igraph_vector_size(weights) < no_of_edges)
      IGRAPH_ERROR("cannot calculate modularity, weight vector too short",
        IGRAPH_EINVAL);
    m=igraph_vector_sum(weights);
    for (i=0; i<no_of_edges; i++) {
      igraph_real_t w=VECTOR(*weights)[i];
      igraph_edge(graph, i, &from, &to);
      c1=VECTOR(*membership)[(long int)from];
      c2=VECTOR(*membership)[(long int)to];
      if (c1==c2) VECTOR(e)[c1] += 2*w;
      VECTOR(a)[c1] += w;
      VECTOR(a)[c2] += w;
    }
  } else {
    m=no_of_edges;
    for (i=0; i<no_of_edges; i++) {
      igraph_edge(graph, i, &from, &to);
      c1=VECTOR(*membership)[(long int)from];
      c2=VECTOR(*membership)[(long int)to];
      if (c1==c2) VECTOR(e)[c1] += 2;
      VECTOR(a)[c1] += 1;
      VECTOR(a)[c2] += 1;
    }
  }

  *modularity=0.0;
  for (i=0; i<types; i++) {
    igraph_real_t tmp=VECTOR(a)[i]/2/m;
    *modularity += VECTOR(e)[i]/2/m;
    *modularity -= tmp*tmp;
  }
  
  igraph_vector_destroy(&e);
  igraph_vector_destroy(&a);
  IGRAPH_FINALLY_CLEAN(2);
  
  return 0;
}

/**
 * \section about_leading_eigenvector_methods
 * 
 * <para>
 * The functions documented in these section implement the 
 * <quote>leading eigenvector</quote> method developed by Mark Newman and 
 * published in MEJ Newman: Finding community structure using the
 * eigenvectors of matrices, arXiv:physics/0605087. TODO: proper
 * citation.</para>
 * 
 * <para>
 * The heart of the method is the definition of the modularity matrix,
 * B, which is B=A-P, A being the adjacency matrix of the (undirected)
 * network, and P contains the probability that certain edges are
 * present according to the <quote>configuration model</quote> In
 * other words, a Pij element of P is the probability that there is an
 * edge between vertices i and j in a random network in which the
 * degrees of all vertices are the same as in the input graph.</para>
 * 
 * <para>
 * The leading eigenvector method works by calculating the eigenvector
 * of the modularity matrix for the largest positive eigenvalue and
 * then separating vertices into two community based on the sign of
 * the corresponding element in the eigenvector. If all elements in
 * the eigenvector are of the same sign that means that the network
 * has no underlying comuunity structure.
 * Check Newman's paper to understand why this is a good method for
 * detecting community structure. </para>
 * 
 * <para>
 * Three function are implemented, they all work accoding to the same
 * principles. The simplest is perhaps \ref
 * igraph_community_leading_eigenvector_naive(). This function splits
 * the network as described above and then recursively splits the 
 * two components after the split as individual networks, if possible.
 * This however is not a good way for maximizing moduilarity, again
 * see the paper for explanation and the proper definition of
 * modularity.</para>
 * 
 * <para>
 * The correct recursive community structure detection method is 
 * implemented in \ref igraph_community_leading_eigenvector(). 
 * Here, after the initial split, the following splits are done in a
 * way to optimize modularity regarding the original network. 
 * I can't say it enough, see the paper, particularly section VI.
 * </para>
 * 
 * <para>
 * The third function is \ref
 * igraph_community_leading_eigenvector_step(), this starts from a
 * division of the network and tries to split a given community into 
 * two subcommunities via the same (correct) method as \ref
 * igraph_community_leading_eigenvector().
 * </para>
 */

typedef struct igraph_i_community_leading_eigenvector_naive_data_t {
  igraph_vector_t *idx;
  igraph_adjlist_t *adjlist;
} igraph_i_community_leading_eigenvector_naive_data_t;

int igraph_i_community_leading_eigenvector_naive(igraph_real_t *to,
						 const igraph_real_t *from,
						 long int n, void *extra) {

  igraph_i_community_leading_eigenvector_naive_data_t *data=extra;
  long int j, k, nlen, size=n;
  igraph_vector_t *idx=data->idx;
  igraph_adjlist_t *adjlist=data->adjlist;
  igraph_real_t ktx, ktx2;

 /* Calculate Ax first */
  for (j=0; j<size; j++) {
    long int oldid=VECTOR(*idx)[j];
    igraph_vector_t *neis=igraph_adjlist_get(adjlist, oldid);
    nlen=igraph_vector_size(neis);
    to[j]=0.0;
    for (k=0; k<nlen; k++) {
      long int nei=VECTOR(*neis)[k];
      to[j] += from[nei];
    }
  }
  
  /* Now calculate k^Tx/2m */
  ktx=0.0; ktx2=0;
  for (j=0; j<size; j++) {
    long int oldid=VECTOR(*idx)[j];
    igraph_vector_t *neis=igraph_adjlist_get(adjlist, oldid);
    long int degree=igraph_vector_size(neis);
    ktx2 += degree;
    ktx += from[j] * degree;
  }
  ktx = ktx/ ktx2; 
  
  /* Now calculate Bx */
  for (j=0; j<size; j++) {
    long int oldid=VECTOR(*idx)[j];
    igraph_vector_t *neis=igraph_adjlist_get(adjlist, oldid);
    long int degree=igraph_vector_size(neis);
    to[j] = to[j] - ktx*degree + degree*degree*from[j]/ktx2;
  }

  return 0;
}

/**
 * \ingroup communities
 * \function igraph_community_leading_eigenvector_naive
 * \brief Leading eigenvector community finding (naive version).
 * 
 * A naive implementation of Newman's eigenvector community structure
 * detection. This function splits the network into two components 
 * according to the leading eigenvector of the modularity matrix and
 * then recursively takes \p steps steps by splitting the components
 * as individual network. This is not the correct way however, see
 * MEJ Newman: Finding community structure in networks using the
 * eigenvectors of matrices, arXiv:physics/0605087. Consider using the 
 * correct \ref igraph_community_leading_eigenvector() function instead.
 * \param graph The input graph, should be undirected to make sense.
 * \param merges The merge matrix. The splits done by the algorithm
 *    are stored here, its structure is the same ad for \ref
 *    igraph_community_leading_eigenvector(). This argument is ignored
 *    if it is \c NULL.
 * \param membership The membership vector, for each vertex it gives
 *    the id of its community after all the splits are performed.
 *    This argument is ignored if it is \c NULL.
 * \param steps The number of splits to do, if possible. Supply the
 *    number of vertices in the network here to perform as many steps 
 *    as possible.
 * \param options The options for ARPACK. \c n is always
 *    overwritten. \c ncv is set to at least 3.
 * \return Error code.
 * 
 * \sa \ref igraph_community_leading_eigenvector() for the proper way, 
 * \ref igraph_community_leading_eigenvector_step() to do just one split.
 * 
 * Time complexity: O(E|+|V|^2*steps), |V| is the number of vertices,
 * |E| is the number of edges.
 */ 

int igraph_community_leading_eigenvector_naive(const igraph_t *graph,
					       igraph_matrix_t *merges,
					       igraph_vector_t *membership,
					       igraph_integer_t steps,
					       igraph_arpack_options_t *options) {
  
  long int no_of_nodes=igraph_vcount(graph);
  igraph_dqueue_t tosplit;
  igraph_vector_t mymerges;
  igraph_vector_t idx;
  long int staken=0;
  igraph_adjlist_t adjlist;
  long int i, j, k, l, m;
  long int communities=1;
  igraph_vector_t vmembership, *mymembership=membership;  
  igraph_i_community_leading_eigenvector_naive_data_t extra;
  igraph_arpack_storage_t storage;

  if (igraph_is_directed(graph)) {
    IGRAPH_WARNING("This method was developed for undirected graphs");
  }
  
  if (steps < 0 || steps > no_of_nodes-1) {
    steps=no_of_nodes-1;
  }
  
  if (!membership) {
    mymembership=&vmembership;
    IGRAPH_VECTOR_INIT_FINALLY(mymembership, 0);
  }

  IGRAPH_VECTOR_INIT_FINALLY(&mymerges, 0);
  IGRAPH_CHECK(igraph_vector_reserve(&mymerges, steps*2));
  IGRAPH_VECTOR_INIT_FINALLY(&idx, no_of_nodes);
  IGRAPH_DQUEUE_INIT_FINALLY(&tosplit, 100);
  IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist, IGRAPH_ALL));
  IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist);
  
  IGRAPH_CHECK(igraph_vector_resize(mymembership, no_of_nodes));
  igraph_vector_null(mymembership);  

  igraph_dqueue_push(&tosplit, 0);

  if (options->ncv<3) { options->ncv=3; }

  /* Memory for ARPACK */
  IGRAPH_CHECK(igraph_arpack_storage_init(&storage, no_of_nodes, options->ncv, 
					  no_of_nodes, 1));
  IGRAPH_FINALLY(igraph_arpack_storage_destroy, &storage);
  extra.idx=&idx;
  extra.adjlist=&adjlist;
  
  while (!igraph_dqueue_empty(&tosplit) && staken < steps) {
    long int comm=igraph_dqueue_pop_back(&tosplit); /* depth first search */
    long int size=0;

    IGRAPH_ALLOW_INTERRUPTION();

    for (i=0; i<no_of_nodes; i++) {
      if (VECTOR(*mymembership)[i]==comm) {
	VECTOR(idx)[size++]=i;
      }
    }

    /* now 'size' is the size of the current community and 
       idx[0:(size-1)] contains the original ids of the vertices in 
       the current community. We need this to index the neighbor list.  */

    staken++;
    if (size==1) {
      continue;			/* nothing to do */
    }

    options->start=0;
    options->n=size;
    options->which[0]='L'; options->which[1]='A';
    if (options->ncv<3) { options->ncv=3; }
    if (options->ncv > options->n) { options->ncv=options->n; }
    
    /* Call ARPACK solver */
    IGRAPH_CHECK(igraph_arpack_rssolve(igraph_i_community_leading_eigenvector_naive,
				       &extra, options, &storage, 0, 0));

    if (options->noiter > options->mxiter) {
      IGRAPH_WARNING("Maximum number of ARPACK iterations reached");
    }

    /* just to have the always the same result, we multiply by -1
       if the first (nonzero) element is not positive  */
    for (i=0; i<size; i++) {
      if (storage.v[i] != 0) { break; }
    }
    if (storage.v[i]<0) {
      for (; i<size; i++) {
	storage.v[i] = - storage.v[i];
      }
    }

    /* Ok, we have the eigenvector */

    /* Non-positive eigenvalue */
/*     printf("%f\n", storage.d[0]); */
    if (storage.d[0] <= 0) { continue; }
    
    /* We create an index vector in workd to renumber the vertices */
    l=0; m=0;
    for (j=0; j<size; j++) {
      if (storage.v[j] <= 0) {
	storage.workd[j]=l++;
      } else {
	storage.workd[j]=m++;
      }
    }    
    /* if l==0 or m==0 then there was no split */
    if (l==0 || m==0) {
      continue;
    }
    communities++;
    
    /* Rewrite the adjacency lists */
    for (j=0; j<size; j++) {
      long int oldid=VECTOR(idx)[j];
      igraph_vector_t *neis=igraph_adjlist_get(&adjlist, oldid);
      long int n=igraph_vector_size(neis);
      for (k=0; k<n; ) {
	long int nei=VECTOR(*neis)[k];
	if ( (storage.v[j] <= 0 && storage.v[nei] <= 0) ||
	     (storage.v[j] > 0 && storage.v[nei] > 0)) {
	  /* they remain in the same community */
	  VECTOR(*neis)[k] = storage.workd[nei];
	  k++;
	} else {
	  /* nei in the other community, remove from neighbor list */
	  VECTOR(*neis)[k] = VECTOR(*neis)[n-1];
	  igraph_vector_pop_back(neis);
	  n--;
	}
      }
    }    

    /* Also rewrite the mymembership vector */
    for (j=0; j<size; j++) {
      if (storage.v[j] <= 0) {
	long int oldid=VECTOR(idx)[j];
	VECTOR(*mymembership)[oldid]=communities-1;
      }
    }

    /* Record merge */
    IGRAPH_CHECK(igraph_vector_push_back(&mymerges, comm));
    IGRAPH_CHECK(igraph_vector_push_back(&mymerges, communities-1));

    /* Store the resulting communities in the queue if needed */
    if (l > 1) {
      IGRAPH_CHECK(igraph_dqueue_push(&tosplit, communities-1));
    }
    if (m > 1) {
      IGRAPH_CHECK(igraph_dqueue_push(&tosplit, comm));
    }

  }
  
  igraph_arpack_storage_destroy(&storage);
  igraph_adjlist_destroy(&adjlist);
  igraph_dqueue_destroy(&tosplit);
  IGRAPH_FINALLY_CLEAN(3);

  /* reform the mymerges vector into merges matrix */
  if (merges) {
    igraph_vector_null(&idx);
    l=igraph_vector_size(&mymerges);
    k=communities;
    j=0;
    IGRAPH_CHECK(igraph_matrix_resize(merges, l/2, 2));
    for (i=l; i>0; i-=2) {
      long int from=VECTOR(mymerges)[i-1];
      long int to=VECTOR(mymerges)[i-2];
      MATRIX(*merges, j, 0)=VECTOR(mymerges)[i-2];
      MATRIX(*merges, j, 1)=VECTOR(mymerges)[i-1];    
      if (VECTOR(idx)[from]!=0) {
	MATRIX(*merges, j, 1)=VECTOR(idx)[from]-1;
      }
      if (VECTOR(idx)[to]!=0) {
	MATRIX(*merges, j, 0)=VECTOR(idx)[to]-1;
      }
      VECTOR(idx)[to]=++k;
      j++;
    }  
  }

  igraph_vector_destroy(&idx);
  igraph_vector_destroy(&mymerges);
  IGRAPH_FINALLY_CLEAN(2);

  if (!membership) {
    igraph_vector_destroy(mymembership);
    IGRAPH_FINALLY_CLEAN(1);
  }
 
  return 0;
}

typedef struct igraph_i_community_leading_eigenvector_data_t {
  igraph_vector_t *idx;
  igraph_vector_t *idx2;
  igraph_adjlist_t *adjlist;
  igraph_vector_t *tmp;
  long int no_of_edges;
  igraph_vector_t *mymembership;
  long int comm;
} igraph_i_community_leading_eigenvector_data_t;

int igraph_i_community_leading_eigenvector(igraph_real_t *to,
					   const igraph_real_t *from,
					   long int n, void *extra) {
  
  igraph_i_community_leading_eigenvector_data_t *data=extra;
  long int j, k, nlen, size=n;
  igraph_vector_t *idx=data->idx;
  igraph_vector_t *idx2=data->idx2;
  igraph_vector_t *tmp=data->tmp;
  igraph_adjlist_t *adjlist=data->adjlist;
  igraph_real_t ktx, ktx2;
  long int no_of_edges=data->no_of_edges;
  igraph_vector_t *mymembership=data->mymembership;
  long int comm=data->comm;
  
  for (j=0; j<size; j++) {
    long int oldid=VECTOR(*idx)[j];
    igraph_vector_t *neis=igraph_adjlist_get(adjlist, oldid);
    nlen=igraph_vector_size(neis);
    to[j]=0.0;
    VECTOR(*tmp)[j]=0.0;
    for (k=0; k<nlen; k++) {
      long int nei=VECTOR(*neis)[k];
      if (VECTOR(*mymembership)[nei]==comm) {
	to[j] += from[ (long int) VECTOR(*idx2)[nei] ];
	VECTOR(*tmp)[j] += 1;
      }
    }
  }
  
  /* Now calculate k^Tx/2m */
  ktx=0.0; ktx2=0.0;
  for (j=0; j<size; j++) {
    long int oldid=VECTOR(*idx)[j];
    igraph_vector_t *neis=igraph_adjlist_get(adjlist, oldid);
    long int degree=igraph_vector_size(neis);
    ktx += from[j] * degree;
    ktx2 += degree;
  }
  ktx = ktx / no_of_edges/2.0;
  ktx2 = ktx2 / no_of_edges/2.0;
  
  /* Now calculate Bx */
  for (j=0; j<size; j++) {
    long int oldid=VECTOR(*idx)[j];
    igraph_vector_t *neis=igraph_adjlist_get(adjlist, oldid);
    igraph_real_t degree=igraph_vector_size(neis);
    to[j] = to[j] - ktx*degree + degree*degree*from[j]/no_of_edges/2.0;
    VECTOR(*tmp)[j] = VECTOR(*tmp)[j] - ktx2*degree + 
      degree*degree*1.0/no_of_edges/2.0;
  }
  
  /* -d_ij summa l in G B_il */
  for (j=0; j<size; j++) {
    to[j] -= VECTOR(*tmp)[j] * from[j];
  }

  return 0;
}

/**
 * \ingroup communities
 * \function igraph_community_leading_eigenvector
 * \brief Leading eigenvector community finding (proper version).
 * 
 * Newman's leading eigenvector method for detecting community
 * structure. This is the proper implementation of the recursive,
 * divisive algorithm: each split is done by maximizing the modularity 
 * regarding the original network, see MEJ Newman: Finding community
 * structure in networks using the eigenvectors of matrices,
 * arXiv:physics/0605087.
 * 
 * \param graph The undirected input graph.
 * \param merges The result of the algorithm, a matrix containing the
 *    information about the splits performed. The matrix is built in
 *    the opposite way however, it is like the result of an
 *    agglomerative algorithm. If at the end of the algorithm (after
 *    \p steps steps was done) there are <quote>p</quote> communities,
 *    then these are numbered from zero to <quote>p-1</quote>. The
 *    first line of the matrix contains the first <quote>merge</quote>
 *    (which is in reality the last split) of two communities into
 *    community <quote>p</quote>, the merge in the second line forms 
 *    community <quote>p+1</quote>, etc. The matrix should be
 *    initialized before calling and will be resized as needed.
 *    This argument is ignored of it is \c NULL.
 * \param membership The membership of the vertices after all the
 *    splits were performed will be stored here. The vector must be
 *    initialized  before calling and will be resized as needed.
 *    This argument is ignored if it is \c NULL.
 * \param steps The maximum number of steps to perform. It might
 *    happen that some component (or the whole network) has no
 *    underlying community structure and no further steps can be
 *    done. If you wany as many steps as possible then supply the 
 *    number of vertices in the network here.
 * \param options The options for ARPACK. \c n is always
 *    overwritten. \c ncv is set to at least 3.
 * \return Error code.
 * 
 * \sa \ref igraph_community_walktrap() and \ref
 * igraph_community_spinglass() for other community structure
 * detection methods.
 * 
 * Time complexity: O(|E|+|V|^2*steps), |V| is the number of vertices,
 * |E| the number of edges, <quote>steps</quote> the number of splits
 * performed.
 */

int igraph_community_leading_eigenvector(const igraph_t *graph,
					 igraph_matrix_t *merges,
					 igraph_vector_t *membership,
					 igraph_integer_t steps,
					 igraph_arpack_options_t *options) {

  long int no_of_nodes=igraph_vcount(graph);
  long int no_of_edges=igraph_ecount(graph);
  igraph_dqueue_t tosplit;
  igraph_vector_t idx, idx2, mymerges;
  igraph_vector_t tmp;
  long int staken=0;
  igraph_adjlist_t adjlist;
  long int i, j, k, l;
  long int communities=1;
  igraph_vector_t vmembership, *mymembership=membership;
  igraph_i_community_leading_eigenvector_data_t extra;
  igraph_arpack_storage_t storage;
  
  if (igraph_is_directed(graph)) {
    IGRAPH_WARNING("This method was developed for undirected graphs");
  }
  
  if (steps < 0 || steps > no_of_nodes-1) {
    steps=no_of_nodes-1;
  }
  
  if (steps > no_of_nodes-1) {
    steps=no_of_nodes-1;
  }
  
  if (!membership) {
    mymembership=&vmembership;
    IGRAPH_VECTOR_INIT_FINALLY(mymembership, 0);
  }
  
  IGRAPH_VECTOR_INIT_FINALLY(&mymerges, 0);
  IGRAPH_CHECK(igraph_vector_reserve(&mymerges, steps*2));

  IGRAPH_VECTOR_INIT_FINALLY(&idx, no_of_nodes);
  IGRAPH_VECTOR_INIT_FINALLY(&idx2, no_of_nodes);
  IGRAPH_VECTOR_INIT_FINALLY(&tmp, no_of_nodes);
  IGRAPH_DQUEUE_INIT_FINALLY(&tosplit, 100);
  IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist, IGRAPH_ALL));
  IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist);
  
  IGRAPH_CHECK(igraph_vector_resize(mymembership, no_of_nodes));
  igraph_vector_null(mymembership);  

  igraph_dqueue_push(&tosplit, 0);

  if (options->ncv<3) { options->ncv=3; }

  /* Memory for ARPACK */
  IGRAPH_CHECK(igraph_arpack_storage_init(&storage, no_of_nodes, options->ncv, 
					  no_of_nodes, 1));
  IGRAPH_FINALLY(igraph_arpack_storage_destroy, &storage);
  extra.idx=&idx;
  extra.idx2=&idx2;
  extra.tmp=&tmp;
  extra.adjlist=&adjlist;
  extra.no_of_edges=no_of_edges;
  extra.mymembership=mymembership;

  while (!igraph_dqueue_empty(&tosplit) && staken < steps) {
    long int comm=igraph_dqueue_pop_back(&tosplit); /* depth first search */
    long int size=0;

    IGRAPH_ALLOW_INTERRUPTION();

    for (i=0; i<no_of_nodes; i++) {
      if (VECTOR(*mymembership)[i]==comm) {
	VECTOR(idx)[size]=i;
	VECTOR(idx2)[i]=size++;
      }
    }

    staken++;
    if (size==1) {
      continue;
    }
    
    options->start=0;
    options->n=size;
    options->which[0]='L'; options->which[1]='A';
    if (options->ncv < 3) { options->ncv=3; }
    if (options->ncv > options->n) { options->ncv=options->n; }
    extra.comm=comm;
    
    /* Call ARPACK solver */
    IGRAPH_CHECK(igraph_arpack_rssolve(igraph_i_community_leading_eigenvector,
				       &extra, options, &storage, 0, 0));
    
    if (options->noiter > options->mxiter) {
      IGRAPH_WARNING("Maximum number of ARPACK iterations reached");
    }

    /* just to have the always the same result, we multiply by -1
       if the first (nonzero) element is not positive  */
    for (i=0; i<size; i++) {
      if (storage.v[i] != 0) { break; }
    }
    if (storage.v[i]<0) {
      for (; i<size; i++) {
	storage.v[i] = - storage.v[i];
      }
    }

    /* Ok, we have the eigenvector */

    /* Non-positive eigenvalue */
/*     printf("%f\n", storage.d[0]); */
/*     for (j=0; j<size; j++) { printf("%g ", storage.v[j]); } */
/*     printf("\n"); */
    if (storage.d[0] <= 0) { continue; }

    /* Count the number of vertices in each community after the split */
    l=0;
    for (j=0; j<size; j++) {
      if (storage.v[j] <= 0) {
	l++;
      }
    }
    if (l==0 || l==size) {
      continue;
    }
    communities++;
    
    /* Rewrite the mymembership vector */
    for (j=0; j<size; j++) {
      if (storage.v[j] <= 0) {
	long int oldid=VECTOR(idx)[j];
	VECTOR(*mymembership)[oldid]=communities-1;
      }
    }

    /* Record merge */
    IGRAPH_CHECK(igraph_vector_push_back(&mymerges, comm));
    IGRAPH_CHECK(igraph_vector_push_back(&mymerges, communities-1));

    /* Store the resulting communities in the queue if needed */
    if (l > 1) {
      IGRAPH_CHECK(igraph_dqueue_push(&tosplit, communities-1));
    }
    if (size-l > 1) {
      IGRAPH_CHECK(igraph_dqueue_push(&tosplit, comm));
    }
    
  }
  
  igraph_arpack_storage_destroy(&storage);
  igraph_adjlist_destroy(&adjlist);
  igraph_dqueue_destroy(&tosplit);
  igraph_vector_destroy(&tmp);
  igraph_vector_destroy(&idx2);
  IGRAPH_FINALLY_CLEAN(5);

  /* reform the mymerges vector */
  if (merges) {
    igraph_vector_null(&idx);
    l=igraph_vector_size(&mymerges);
    k=communities;
    j=0;
    IGRAPH_CHECK(igraph_matrix_resize(merges, l/2, 2));
    for (i=l; i>0; i-=2) {
      long int from=VECTOR(mymerges)[i-1];
      long int to=VECTOR(mymerges)[i-2];
      MATRIX(*merges, j, 0)=VECTOR(mymerges)[i-2];
      MATRIX(*merges, j, 1)=VECTOR(mymerges)[i-1];    
      if (VECTOR(idx)[from]!=0) {
	MATRIX(*merges, j, 1)=VECTOR(idx)[from]-1;
      }
      if (VECTOR(idx)[to]!=0) {
	MATRIX(*merges, j, 0)=VECTOR(idx)[to]-1;
      }
      VECTOR(idx)[to]=++k;
      j++;
    }      
  }
  
  igraph_vector_destroy(&idx);
  igraph_vector_destroy(&mymerges);
  IGRAPH_FINALLY_CLEAN(2);
  
  if (!membership) {
    igraph_vector_destroy(mymembership);
    IGRAPH_FINALLY_CLEAN(1);
  }
  
  return 0;
}

typedef struct igraph_i_community_leading_eigenvector_step_data_t {
  igraph_vector_t *idx;
  igraph_vector_t *idx2;
  igraph_lazy_adjlist_t *adjlist;
  igraph_vector_t *tmp;
  long int no_of_edges;
  igraph_vector_t *mymembership;
  long int comm;
} igraph_i_community_leading_eigenvector_step_data_t;

int igraph_i_community_leading_eigenvector_step(igraph_real_t *to,
						const igraph_real_t *from,
						long int n, void *extra) {
  
  igraph_i_community_leading_eigenvector_step_data_t *data=extra;
  long int j, k, nlen, size=n;
  igraph_vector_t *idx=data->idx;
  igraph_vector_t *idx2=data->idx2;
  igraph_vector_t *tmp=data->tmp;
  igraph_lazy_adjlist_t *adjlist=data->adjlist;
  igraph_real_t ktx, ktx2;
  long int no_of_edges=data->no_of_edges;
  igraph_vector_t *mymembership=data->mymembership;
  long int comm=data->comm;
  
  for (j=0; j<size; j++) {
    long int oldid=VECTOR(*idx)[j];
    igraph_vector_t *neis=igraph_lazy_adjlist_get(adjlist, oldid);
    nlen=igraph_vector_size(neis);
    to[j]=0.0;
    VECTOR(*tmp)[j]=0.0;
    for (k=0; k<nlen; k++) {
      long int nei=VECTOR(*neis)[k];
      if (VECTOR(*mymembership)[nei]==comm) {
	to[j] += from[ (long int) VECTOR(*idx2)[nei] ];
	VECTOR(*tmp)[j] += 1;
      }
    }
  }
  
  /* Now calculate k^Tx/2m */
  ktx=0.0; ktx2=0.0;
  for (j=0; j<size; j++) {
    long int oldid=VECTOR(*idx)[j];
    igraph_vector_t *neis=igraph_lazy_adjlist_get(adjlist, oldid);
    long int degree=igraph_vector_size(neis);
    ktx += from[j] * degree;
    ktx2 += degree;
  }
  ktx = ktx / no_of_edges/2.0;
  ktx2 = ktx2 / no_of_edges/2.0;
  
  /* Now calculate Bx */
  for (j=0; j<size; j++) {
    long int oldid=VECTOR(*idx)[j];
    igraph_vector_t *neis=igraph_lazy_adjlist_get(adjlist, oldid);
    igraph_real_t degree=igraph_vector_size(neis);
    to[j] = to[j] - ktx*degree + degree*degree*from[j]/no_of_edges/2.0;
    VECTOR(*tmp)[j] = VECTOR(*tmp)[j] - ktx2*degree + 
      degree*degree*1.0/no_of_edges/2.0;
  }
  
  /* -d_ij summa l in G B_il */
  for (j=0; j<size; j++) {
    to[j] -= VECTOR(*tmp)[j] * from[j];
  }

  return 0;
}

/**
 * \ingroup communities
 * \function igraph_community_leading_eigenvector_step
 * \brief Leading eigenvector community finding (make one step).
 * 
 * Do one split according to Mark Newman's leading eigenvector
 * community detection method. See MEJ Newman: Finding community
 * structure in networks using the eigenvectors of matrices,
 * arXiv:phyisics/0605087 for the details.
 * 
 * </para><para>Use this function instead of \ref
 * igraph_community_leading_eigenvector() if you want to have full
 * control over and information about each split performed along
 * community structure detection. \ref
 * igraph_community_leading_eigenvector() can be simulated by
 * repeatedly calling this function.
 * 
 * \param graph The undirected input graph.
 * \param membership Numeric vector giving a division of \p graph.
 *    The result will be also stored here. The vector contains the
 *    community ids for each vertex, these are numbered from 0.
 * \param community The id of the community to split.
 * \param split Pointer to a logical variable, if it was possible to
 *    split community \p community then 1, otherwise 0 will be stored
 *    here. This argument is ignored if it is \c NULL.
 * \param eigenvector Pointer to an initialized vector, the
 *    eigenvector on which the split was done will be stored here. 
 *    It will be resised to have the same length as the number of
 *    vertices in community \p community. This argument is ignored 
 *    if it is \c NULL.
 * \param eigenvalue Pointer to a real variable, the eigenvalue
 *    associated with \p eigenvector will be stored here.
 *    This argument is ignored if it is \c NULL.
 * \return Error code.
 * 
 * \sa \ref igraph_community_leading_eigenvector().
 * 
 * Time complexity: O(|E|+|V|^2), |E| is the number of edges, |V| is
 * the number of vertices.
 */

int igraph_community_leading_eigenvector_step(const igraph_t *graph,
					      igraph_vector_t *membership,
					      igraph_integer_t community,
					      igraph_bool_t *split,
					      igraph_vector_t *eigenvector,
					      igraph_real_t *eigenvalue, 
					      igraph_arpack_options_t *options,
					      igraph_arpack_storage_t *storage) {

  long int no_of_nodes=igraph_vcount(graph);
  long int no_of_edges=igraph_ecount(graph);
  igraph_vector_t tmp;
  igraph_vector_t idx, idx2;
  long int i, j, k;
  long int communities=1;
  igraph_lazy_adjlist_t adjlist;
  long int size=0;
  igraph_i_community_leading_eigenvector_step_data_t extra;
  igraph_arpack_storage_t real_storage, *mystorage= 
    storage ? storage : &real_storage;
  long int comm=community;
  
  if (igraph_vector_size(membership) != no_of_nodes) {
    IGRAPH_ERROR("Invalid membership vector length", IGRAPH_EINVAL);
  }
  
  if (igraph_is_directed(graph)) {
    IGRAPH_WARNING("This method was developed for undirected graphs");
  }

  IGRAPH_VECTOR_INIT_FINALLY(&idx, no_of_nodes);
  IGRAPH_VECTOR_INIT_FINALLY(&idx2, no_of_nodes);

  for (i=0; i<no_of_nodes; i++) {
    if (VECTOR(*membership)[i]==comm) {
      VECTOR(idx)[size]=i;
      VECTOR(idx2)[i]=size;
      size++;
    }
    if (VECTOR(*membership)[i] > communities-1) {
      communities=VECTOR(*membership)[i]+1;
    }
  }
  
  if (split) { *split=0; }
  if (size != 1) {
    IGRAPH_CHECK(igraph_lazy_adjlist_init(graph, &adjlist, IGRAPH_ALL, 
					  IGRAPH_DONT_SIMPLIFY));  
    IGRAPH_FINALLY(igraph_lazy_adjlist_destroy, &adjlist);  
    if (!storage) {
      IGRAPH_CHECK(igraph_arpack_storage_init(mystorage, no_of_nodes, 3, no_of_nodes, 1));
      IGRAPH_FINALLY(igraph_arpack_storage_destroy, mystorage);
    }
    IGRAPH_VECTOR_INIT_FINALLY(&tmp, size);

    extra.idx=&idx;
    extra.idx2=&idx2;
    extra.tmp=&tmp;
    extra.adjlist=&adjlist;
    extra.no_of_edges=no_of_edges;
    extra.mymembership=membership;
    extra.comm=comm;
   
    options->start=0;
    options->n=size;
    if (options->ncv < 3) { options->ncv=3; }
    options->which[0]='L'; options->which[1]='A';
    if (options->ncv > options->n) { options->ncv=options->n; }

    IGRAPH_CHECK(igraph_arpack_rssolve(igraph_i_community_leading_eigenvector_step,
				       &extra, options, mystorage, 0, 0));
    
    if (options->noiter > options->mxiter) {
      IGRAPH_WARNING("Maximum number of ARPACK iterations reached");
    }
    
    /* just to have the always the same result, we multiply by -1
       if the first (nonzero) element is not positive  */
    for (i=0; i<size; i++) {
      if (mystorage->v[i] != 0) { break; }
    }
    if (mystorage->v[i]<0) {
      for (; i<size; i++) {
	mystorage->v[i] = - mystorage->v[i];
      }
    }

    /* Ok, we have the eigenvector */
    
    /* Save eigenvalue/vector if requested */
    if (eigenvalue) {
      *eigenvalue=mystorage->d[0];
    }
    if (eigenvector) {
      IGRAPH_CHECK(igraph_vector_resize(eigenvector, size));
      for (i=0; i<size; i++) {
	VECTOR(*eigenvector)[i] = mystorage->v[i];
      }
    }
    
    /* Positive eigenvalue? */
    if (mystorage->d[0] > 0) {
      
      /* Rewrite the membership vector, check if there was a split */
      for (j=0, k=0; j<size; j++) {
	if (VECTOR(*eigenvector)[j] <= 0) {
	  long int oldid=VECTOR(idx)[j];
	  VECTOR(*membership)[oldid]=communities;
	  k++;
	}
      }
      
      if (split && k>0) {
	*split=1;
      }
    }
    
    igraph_vector_destroy(&tmp);
    IGRAPH_FINALLY_CLEAN(1);
    if (!storage) { 
      igraph_arpack_storage_destroy(mystorage); 
      IGRAPH_FINALLY_CLEAN(1);
    }
    igraph_lazy_adjlist_destroy(&adjlist);
    IGRAPH_FINALLY_CLEAN(1);
  
  } /* size != 1 */

  igraph_vector_destroy(&idx2);
  igraph_vector_destroy(&idx);
  IGRAPH_FINALLY_CLEAN(2);
  
  return 0;
}

/**
 * \function igraph_le_community_to_membership
 * Vertex membership from the leading eigenvector community structure
 * 
 * This function creates a membership vector from the
 * result of \ref igraph_community_leading_eigenvector() or
 * \ref igraph_community_leading_eigenvector_naive(). It takes \c membership
 * and permformes \c steps merges, according to the supplied
 * \c merges matrix.
 * \param merges The matrix defining the merges to make. 
 *     This is usually from the output of the leading eigenvector community
 *     structure detection routines.
 * \param steps The number of steps to make according to \c merges.
 * \param membership Initially the starting membership vector, 
 *     on output the resulting membership vector, after performing \c steps merges.
 * \param csize Optionally the sizes of the commmunities is stored here, 
 *     if this is not a null pointer, but an initialized vector.
 * \return Error code.
 * 
 * Time complexity: O(|V|), the number of vertices.
 */

int igraph_le_community_to_membership(const igraph_matrix_t *merges,
				      igraph_integer_t steps,
				      igraph_vector_t *membership,
				      igraph_vector_t *csize) {

  long int no_of_nodes=igraph_vector_size(membership);
  igraph_vector_t fake_memb;
  long int components, i;

  if (igraph_matrix_nrow(merges) < steps) {
    IGRAPH_ERROR("`steps' to big or `merges' matrix too short", IGRAPH_EINVAL);
  }    
  
  components=igraph_vector_max(membership)+1;
  if (components > no_of_nodes) { 
    IGRAPH_ERROR("Invalid membership vector, too many components", IGRAPH_EINVAL);
  }
  if (steps >= components) {
    IGRAPH_ERROR("Cannot make `steps' steps from supplied membership vector",
		 IGRAPH_EINVAL);
  }
  
  IGRAPH_VECTOR_INIT_FINALLY(&fake_memb, components);
  
  /* Check membership vector */
  for (i=0; i<no_of_nodes; i++) {
    if (VECTOR(*membership)[i] < 0) {
      IGRAPH_ERROR("Invalid membership vector, negative id", IGRAPH_EINVAL);
    }
    VECTOR(fake_memb)[ (long int) VECTOR(*membership)[i] ] += 1;
  }
  for (i=0; i<components; i++) {
    if (VECTOR(fake_memb)[i] == 0) {
      IGRAPH_ERROR("Invalid membership vector, empty cluster", IGRAPH_EINVAL);
    }
  }
  
  IGRAPH_CHECK(igraph_community_to_membership(merges, components, steps, 
					      &fake_memb, 0));
  
  /* Ok, now we have the membership of the initial components, 
     rewrite the original membership vector. */

  if (csize) {
    IGRAPH_CHECK(igraph_vector_resize(csize, components-steps));
    igraph_vector_null(csize);
  }

  for (i=0; i<no_of_nodes; i++) {
    VECTOR(*membership)[i] = VECTOR(fake_memb)[ (long int) VECTOR(*membership)[i] ];
    if (csize) {
      VECTOR(*csize)[ (long int) VECTOR(*membership)[i] ] += 1;
    }
  }

  igraph_vector_destroy(&fake_memb);
  IGRAPH_FINALLY_CLEAN(1);
  return 0;  
}

/**
 * \function igraph_community_label_propagation
 * \ingroup communities
 * \brief Community detection based on label propagation
 *
 * This function implements the community detection method described in:
 * Raghavan, U.N. and Albert, R. and Kumara, S.: Near linear time algorithm
 * to detect community structures in large-scale networks. Phys Rev E
 * 76, 036106. (2007). This version extends the original method by
 * the ability to take edge weights into consideration and also
 * by allowing some labels to be fixed.
 * 
 * \param graph The input graph, should be undirected to make sense.
 * \param membership The membership vector, the result is returned here.
 *    For each vertex it gives the ID of its community (label).
 * \param weights The weight vector, it should contain a positive
 *    weight for all the edges.
 * \param initial The initial state. If NULL, every vertex will have
 *   a different label at the beginning. Otherwise it must be a vector
 *   with an entry for each vertex. Non-negative values denote different
 *   labels, negative entries denote vertices without labels.
 * \param fixed Boolean vector denoting which labels are fixed. Of course
 *   this makes sense only if you provided an initial state, otherwise
 *   this element will be ignored. Also note that vertices without labels
 *   cannot be fixed.
 * \return Error code.
 * 
 * Time complexity: O(m+n)
 */
int igraph_community_label_propagation(const igraph_t *graph,
                                       igraph_vector_t *membership,
                                       const igraph_vector_t *weights,
                                       const igraph_vector_t *initial,
                                       igraph_vector_bool_t *fixed) {
  long int no_of_nodes=igraph_vcount(graph);
  long int no_of_edges=igraph_ecount(graph);
  long int no_of_not_fixed_nodes=no_of_nodes;
  long int i, j, k;
  igraph_adjlist_t al;
  igraph_adjedgelist_t ael;
  igraph_bool_t running = 1;

  igraph_vector_t label_counters, dominant_labels, node_order;

  /* The implementation uses a trick to avoid negative array indexing:
   * elements of the membership vector are increased by 1 at the start
   * of the algorithm; this to allow us to denote unlabeled vertices
   * (if any) by zeroes. The membership vector is shifted back in the end
   */

  /* Do some initial checks */
  if (fixed && igraph_vector_bool_size(fixed) != no_of_nodes) {
    IGRAPH_ERROR("Invalid fixed labeling vector length", IGRAPH_EINVAL);
  }
  if (weights) {
    if (igraph_vector_size(weights) != no_of_edges) {
      IGRAPH_ERROR("Invalid weight vector length", IGRAPH_EINVAL);
    } else if (igraph_vector_min(weights) < 0) {
      IGRAPH_ERROR("Weights must be non-negative", IGRAPH_EINVAL);
    }
  }
  if (fixed && !initial) {
    IGRAPH_WARNING("Ignoring fixed vertices as no initial labeling given");
  }

  IGRAPH_CHECK(igraph_vector_resize(membership, no_of_nodes));

  if (initial) {
    if (igraph_vector_size(initial) != no_of_nodes) {
      IGRAPH_ERROR("Invalid initial labeling vector length", IGRAPH_EINVAL);
    }
    /* Check if the labels used are valid, initialize membership vector */
    for (i=0; i<no_of_nodes; i++) {
      if (VECTOR(*initial)[i] < -1) {
        VECTOR(*membership)[i] = 0;
      } else {
        VECTOR(*membership)[i] = VECTOR(*initial)[i] + 1;
      }
    }
    if (fixed) {
      for (i=0; i<no_of_nodes; i++) {
        if (VECTOR(*fixed)[i]) {
          if (VECTOR(*membership)[i] == 0) {
            IGRAPH_WARNING("Fixed nodes cannot be unlabeled, ignoring them");
            VECTOR(*fixed)[i] = 0;
          } else {
            no_of_not_fixed_nodes--;
          }
        }
      }
    }

    i = igraph_vector_max(membership);
    if (i > no_of_nodes) {
      IGRAPH_ERROR("elements of the initial labeling vector must be between 0 and |V|-1", IGRAPH_EINVAL);
    }
    if (i <= 0) {
      IGRAPH_ERROR("at least one vertex must be labeled in the initial labeling", IGRAPH_EINVAL);
    }
  } else {
    for (i=0; i<no_of_nodes; i++) {
      VECTOR(*membership)[i] = i;
    }
  }

  /* Create an adjacency (edge) list representation for efficiency.
   * For the unweighted case, the adjacency list is enough. For the
   * weighted case, we need the adjacency edge list */
  if (weights) {
    IGRAPH_CHECK(igraph_adjedgelist_init(graph, &ael, IGRAPH_IN));
    IGRAPH_FINALLY(igraph_adjedgelist_destroy, &ael);
  } else {
    IGRAPH_CHECK(igraph_adjlist_init(graph, &al, IGRAPH_IN));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &al);
  }

  /* Create storage space for counting distinct labels and dominant ones */
  IGRAPH_VECTOR_INIT_FINALLY(&label_counters, no_of_nodes+1);
  IGRAPH_VECTOR_INIT_FINALLY(&dominant_labels, 0);
  IGRAPH_CHECK(igraph_vector_reserve(&dominant_labels, 2));

  RNG_BEGIN();

  /* Initialize node ordering vector with only the not fixed nodes */
  if (fixed) {
    IGRAPH_VECTOR_INIT_FINALLY(&node_order, no_of_not_fixed_nodes);
    for (i=0, j=0; i<no_of_nodes; i++) {
      if (!VECTOR(*fixed)[i]) {
        VECTOR(node_order)[j] = i;
        j++;
      }
    }
  } else {
    IGRAPH_CHECK(igraph_vector_init_seq(&node_order, 0, no_of_nodes-1));
    IGRAPH_FINALLY(igraph_vector_destroy, &node_order);
  }

  running = 1;
  while (running) {
    long int v1, num_neis;
    igraph_real_t max_count;
    igraph_vector_t *neis;

    running = 0;

    /* Shuffle the node ordering vector */
    IGRAPH_CHECK(igraph_vector_shuffle(&node_order));
    /* In the prescribed order, loop over the vertices and reassign labels */
    for (i=0; i<no_of_not_fixed_nodes; i++) {
      v1 = VECTOR(node_order)[i];

      /* Count the weights corresponding to different labels */
      igraph_vector_null(&label_counters);
      igraph_vector_clear(&dominant_labels);
      max_count = 0.0;
      if (weights) {
        neis = igraph_adjedgelist_get(&ael, v1);
        num_neis = igraph_vector_size(neis);
        for (j=0; j<num_neis; j++) {
          k = VECTOR(*membership)[(long)IGRAPH_OTHER(graph, VECTOR(*neis)[j], v1)];
          if (k == 0) continue;   /* skip if it has no label yet */
          VECTOR(label_counters)[k] += VECTOR(*weights)[(long)VECTOR(*neis)[j]];
          if (max_count < VECTOR(label_counters)[k]) {
            max_count = VECTOR(label_counters)[k];
            igraph_vector_resize(&dominant_labels, 1);
            VECTOR(dominant_labels)[0] = k;
          } else if (max_count == VECTOR(label_counters)[k]) {
            igraph_vector_push_back(&dominant_labels, k);
          }
        }
      } else {
        neis = igraph_adjlist_get(&al, v1);
        num_neis = igraph_vector_size(neis);
        for (j=0; j<num_neis; j++) {
          k = VECTOR(*membership)[(long)VECTOR(*neis)[j]];
          if (k == 0) continue;   /* skip if it has no label yet */
          VECTOR(label_counters)[k]++;
          if (max_count < VECTOR(label_counters)[k]) {
            max_count = VECTOR(label_counters)[k];
            igraph_vector_resize(&dominant_labels, 1);
            VECTOR(dominant_labels)[0] = k;
          } else if (max_count == VECTOR(label_counters)[k]) {
            igraph_vector_push_back(&dominant_labels, k);
          }
        }
      }

      if (igraph_vector_size(&dominant_labels) > 0) {
        /* Select randomly from the dominant labels */
        k = RNG_INTEGER(0, igraph_vector_size(&dominant_labels)-1); 
        k = VECTOR(dominant_labels)[k];
        /* Check if the _current_ label of the node is also dominant */
        if (VECTOR(label_counters)[(long)VECTOR(*membership)[v1]]!=max_count) {
          /* Nope, we need at least one more iteration */
          running = 1;
        }
        VECTOR(*membership)[v1] = k;
      }
    }
  }

  RNG_END();

  /* Shift back the membership vector, permute labels in increasing order */
  /* We recycle label_counters here :) */
  igraph_vector_fill(&label_counters, -1);
  j = 0;
  for (i=0; i<no_of_nodes; i++) {
    k = (long)VECTOR(*membership)[i]-1;
    if (k >= 0) {
      if (VECTOR(label_counters)[k] == -1) {
        /* We have seen this label for the first time */
        VECTOR(label_counters)[k] = j;
        k = j;
        j++;
      } else {
        k = VECTOR(label_counters)[k];
      }
    } else {
      /* This is an unlabeled vertex */
    }
    VECTOR(*membership)[i] = k;
  }

  if (weights)
    igraph_adjedgelist_destroy(&ael);
  else
    igraph_adjlist_destroy(&al);

  igraph_vector_destroy(&node_order);
  igraph_vector_destroy(&label_counters);
  igraph_vector_destroy(&dominant_labels);
  IGRAPH_FINALLY_CLEAN(4);

  return 0;
}

