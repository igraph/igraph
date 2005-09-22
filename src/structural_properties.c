/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2005  Gabor Csardi <csardi@rmki.kfki.hu>
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

/**
 * \ingroup structural
 * \brief Calculates the diameter of a graph (longest geodesic).
 *
 * @param graph The graph object.
 * @param res Pointer to an integer, this will contain the result.
 * @param directed Boolean, whether to consider directed
 *        paths. Ignored for undirected graphs.
 * @param unconn What to do if the graph is not connected. If
 *        <code>TRUE</code> the longest geodesic within a component
 *        will be returned, otherwise the number of vertices is
 *        returned. (The ratio behind the latter is that this is
 *        always longer than the longest possible diameter in a
 *        graph.) 
 * @return Error code.
 *
 * Time complexity: <code>O(|V||E|)</code>, the number of vertices
 * times the number of edges.
 */

int igraph_diameter(igraph_t *graph, integer_t *res, 
		    bool_t directed, bool_t unconn) {

  long int no_of_nodes=igraph_vcount(graph);
  long int i, j;
  long int *already_added;
  long int nodes_reached;

  dqueue_t q;
  vector_t neis;
  integer_t dirmode;

  *res=0;  
  if (directed) { dirmode=IGRAPH_OUT; } else { dirmode=IGRAPH_ALL; }
  already_added=Calloc(no_of_nodes, long int);
  dqueue_init(&q, 100);
  vector_init(&neis, 0);
  
  for (i=0; i<no_of_nodes; i++) {
    nodes_reached=1;
    dqueue_push(&q, i);
    dqueue_push(&q, 0);
    already_added[i]=i+1;
    
    while (!dqueue_empty(&q)) {
      long int actnode=dqueue_pop(&q);
      long int actdist=dqueue_pop(&q);
      if (actdist>*res) { *res=actdist; }
      
      igraph_neighbors(graph, &neis, actnode, dirmode);
      for (j=0; j<vector_size(&neis); j++) {
	long int neighbor=VECTOR(neis)[j];
	if (already_added[neighbor] == i+1) { continue; }
	already_added[neighbor]=i+1;
	nodes_reached++;
	dqueue_push(&q, neighbor);
	dqueue_push(&q, actdist+1);
      }
    } /* while !dqueue_empty */
    
    /* not connected, return largest possible */
    if (nodes_reached != no_of_nodes && !unconn) {
      *res=no_of_nodes;
      break;
    }
  } /* for i<no_of_nodes */
  
  /* clean */
  Free(already_added);
  dqueue_destroy(&q);
  vector_destroy(&neis);

  return 0;
}

/**
 * \ingroup structural
 * \brief Calculates the average geodesic length in a graph.
 *
 * @param graph The graph object.
 * @param res Pointer to a real number, this will contain the result.
 * @param directed Boolean, whether to consider directed
 *        paths. Ignored for undirected graphs.
 * @param unconn What to do if the graph is not connected. If
 *        <code>TRUE</code> the average of thr geodesics within the components
 *        will be returned, otherwise the number of vertices is
 *        used for the length of non-existing geodesics. (The ratio
 *        behind this is that this is always longer than the longest
 *        possible geodesic in a graph.) 
 * @return Error code.
 *
 * Time complexity: <code>O(|V||E|)</code>, the number of vertices
 * times the number of edges.
 */

int igraph_average_path_length(igraph_t *graph, real_t *res,
			       bool_t directed, bool_t unconn) {
  long int no_of_nodes=igraph_vcount(graph);
  long int i, j;
  long int *already_added;
  long int nodes_reached=0;
  long int normfact=0;

  dqueue_t q;
  vector_t neis;
  integer_t dirmode;

  *res=0;  
  if (directed) { dirmode=IGRAPH_OUT; } else { dirmode=IGRAPH_ALL; }
  already_added=Calloc(no_of_nodes, long int);
  dqueue_init(&q, 100);
  vector_init(&neis, 0);
  
  for (i=0; i<no_of_nodes; i++) {
    nodes_reached=0;
    dqueue_push(&q, i);
    dqueue_push(&q, 0);
    already_added[i]=i+1;
    
    while (!dqueue_empty(&q)) {
      long int actnode=dqueue_pop(&q);
      long int actdist=dqueue_pop(&q);
    
      igraph_neighbors(graph, &neis, actnode, dirmode);
      for (j=0; j<vector_size(&neis); j++) {
	long int neighbor=VECTOR(neis)[j];
	if (already_added[neighbor] == i+1) { continue; }
	already_added[neighbor]=i+1;
	nodes_reached++;
	*res += actdist+1;
	normfact++;
	dqueue_push(&q, neighbor);
	dqueue_push(&q, actdist+1);
      }
    } /* while !dqueue_empty */
    
    /* not connected, return largest possible */
    if (!unconn) {
      *res += (no_of_nodes * (no_of_nodes-1-nodes_reached));
      normfact += no_of_nodes-1-nodes_reached;
    }    
  } /* for i<no_of_nodes */

  printf("%f / %li\n", *res, normfact);
  *res /= normfact;

  /* clean */
  Free(already_added);
  dqueue_destroy(&q);
  vector_destroy(&neis);

  return 0;
}

/**
 * \ingroup structural
 * \brief Calculates one minimum spanning tree of an unweighted graph.
 * 
 * If the graph has more minimum spanning trees (this is always the
 * case, except if it is a forest) this implementation returns only
 * the same one.
 * 
 * Directed graphs are considered as undirected for this computation.
 *
 * If the graph is not connected then its minimum spanning forest is
 * returned. This is the set of the minimum spanning trees of each
 * component.
 * @param graph The graph object.
 * @param mst The minimum spanning tree, another graph object. Do
 *        <em>not</em> initialize this object before passing it to
 *        this function, but be sure to call igraph_destroy() on it if
 *        you don't need it any more.
 * @return Error code.
 *
 * Time complexity: <code>O(|V|+|E|)</code>, <code>|V|</code> is the
 * number of vertices, <code>|E|</code> the number of edges in the
 * graph. 
 *
 * \sa igraph_minimum_spanning_tree_prim() for weighted graphs.
 */

int igraph_minimum_spanning_tree_unweighted(igraph_t *graph, igraph_t *mst) {

  long int no_of_nodes=igraph_vcount(graph);
  char *already_added;
  
  dqueue_t q;
  vector_t edges;
  vector_t tmp;
  long int i, j;

  already_added=Calloc(no_of_nodes, char);
  vector_init(&tmp, 0);

  dqueue_init(&q, 100);
  vector_init(&edges, 0);
  vector_reserve(&edges, (no_of_nodes-1)*2);
  
  for (i=0; i<no_of_nodes; i++) {
    if (already_added[i]>0) { continue; }

    already_added[i]=1;
    dqueue_push(&q, i);
    while (! dqueue_empty(&q)) {
      long int act_node=dqueue_pop(&q);
      igraph_neighbors(graph, &tmp, act_node, 3);
      for (j=0; j<vector_size(&tmp); j++) {
	long int neighbor=VECTOR(tmp)[j];
	if (already_added[neighbor]==0) {
	  already_added[neighbor]=1;
	  dqueue_push(&q, neighbor);
	  vector_push_back(&edges, act_node);
	  vector_push_back(&edges, neighbor);
	}
      }
    }
  }
  
  dqueue_destroy(&q);
  Free(already_added);
  vector_destroy(&tmp);

  igraph_create(mst, &edges, 0, igraph_is_directed(graph));
  vector_destroy(&edges);

  return 0;
}

/**
 * \ingroup structural
 * \brief Calculates one minimum spanning tree of a weighted graph.
 *
 * This function uses Prim's method for carrying out the computation,
 * see Prim, R.C.: Shortest connection networks and some
 * generalizations, <em>Bell System Technical Journal</em>, Vol. 36,
 * 1957, 1389--1401.
 * 
 * If the graph has more than one minimum spanning tree, the current
 * implementation returns always the same one.
 *
 * Directed graphs are considered as undirected for this computation. 
 * 
 * If the graph is not connected then its minimum spanning forest is
 * returned. This is the set of the minimum spanning trees of each
 * component.
 * 
 * @param graph The graph object.
 * @param mst The result of the computation, a graph object containing
 *        the minimum spanning tree of the graph.
 *        Do <em>not</em> initialize this object before passing it to
 *        this function, but be sure to call igraph_destroy() on it if
 *        you don't need it any more.
 * @param weights A vector containing the weights of the the edges.
 *        in the same order as the simple edge iterator visits them.
 * @return Error code.
 *
 * Time complexity: <code>O(|V|+|E|)</code>, <code>|V|</code> is the
 * number of vertices, <code>|E|</code> the number of edges in the
 * graph. 
 *
 * \sa igraph_minimum_spanning_tree_unweighted() for unweighted graphs.
 */

int igraph_minimum_spanning_tree_prim(igraph_t *graph, igraph_t *mst,
				      vector_t *weights) {

  long int no_of_nodes=igraph_vcount(graph);
  char *already_added;

  d_indheap_t heap;
  vector_t edges;
  integer_t mode=3;
  
  igraph_iterator_t it;

  long int i, j;

  already_added=Calloc(no_of_nodes, char);
  
  d_indheap_init(&heap, 0);
  vector_init(&edges, 0);
  vector_reserve(&edges, (no_of_nodes-1)*2);
  igraph_iterator_eneis(graph, &it, 0, mode);

  for (i=0; i<no_of_nodes; i++) {
    if (already_added[i]>0) { continue; }

    already_added[i]=1;
    /* add all edges of the first vertex */
    igraph_iterator_eneis_set(graph, &it, i, mode);
    while( !igraph_end(graph, &it)) {
      long int neighbor=igraph_get_vertex_nei(graph, &it);
      long int edgeno=igraph_get_edge(graph, &it);
      if (already_added[neighbor] == 0) {
	d_indheap_push(&heap, -VECTOR(*weights)[edgeno], i, neighbor);
      }
      igraph_next(graph, &it);
    }

    while(! d_indheap_empty(&heap)) {
      /* Get minimal edge */
      long int from, to;
      d_indheap_max_index(&heap, &from, &to);

      /* Erase it */
      d_indheap_delete_max(&heap);
      
      /* Does it point to a visited node? */
      if (already_added[to]==0) {
	already_added[to]=1;
	vector_push_back(&edges, from);
	vector_push_back(&edges, to);
	/* add all outgoing edges */
	igraph_iterator_eneis_set(graph, &it, to, mode);
	while ( !igraph_end(graph, &it)) {
	  long int neighbor=igraph_get_vertex_nei(graph, &it);
	  long int edgeno=igraph_get_edge(graph, &it);
	  if (already_added[neighbor] == 0) {
	    d_indheap_push(&heap, -VECTOR(*weights)[edgeno], to, neighbor);
	  }
	  igraph_next(graph, &it);
	} /* for */
      } /* if !already_added */
    } /* while in the same component */
  } /* for all nodes */

  d_indheap_destroy(&heap);
  Free(already_added);
  igraph_iterator_destroy(graph, &it);

  igraph_create(mst, &edges, 0, igraph_is_directed(graph));
  vector_destroy(&edges);

  return 0;
}

/**
 * \ingroup structural
 * \brief Cloness centrality calculations for some vertices.
 *
 * The closeness centerality of a vertex measures how easily other
 * vertices can be reached from it (or the other way: how easily it
 * can be reached from the other vertices). It is defined as the
 * number of the number of vertices minus one divided by the sum of the
 * lengths of all geodesics from/to the given vertex.
 *
 * If the graph is not connected, and there is no path between two
 * vertices, the number of vertices is used instead the length of the
 * geodesic. This is always longer than the longest possible geodesic.
 * 
 * @param graph The graph object.
 * @param res The result of the computatuion, a vector containing the
 *        closeness centrality scores for the given vertices.
 * @param vids Vector giving the vertices for which the closeness
 *        centrality scores will be computed.
 * @param mode The type of shortest paths to be use for the
 *        calculation in directed graphs. Possible values: 
 *        - <b>IGRAPH_OUT</b>, the lengths of the outgoing paths are
 *          calculated. 
 *        - <b>IGRAPH_IN</b>, the lengths of the incoming paths are
 *          calculated. 
 *        - <b>IGRAPH_ALL</b>, the directed graph is considered as an
 *          undirected one for the computation.
 * @return Error code.
 *
 * Time complexity: <code>O(n|E|)</code>, <code>n</code> is the number
 * of vertices for which the calculation is done and |E| is the number
 * of edges in the graph.
 *
 * \sa Other centrality types: igraph_degree(), igraph_betweenness().
 */

int igraph_closeness(igraph_t *graph, vector_t *res, vector_t *vids, 
		     igraph_neimode_t mode) {

  long int no_of_nodes=igraph_vcount(graph);
  long int *already_counted;
  long int i, j;
  long int nodes_reached;

  dqueue_t q;
  
  long int nodes_to_calc=vector_size(vids);
  vector_t tmp;

  vector_resize(res, nodes_to_calc);
  vector_null(res);
  
  already_counted=Calloc(no_of_nodes, long int);
  vector_init(&tmp, 0);
  dqueue_init(&q, 100);

  for (i=0; i<nodes_to_calc; i++) {
    dqueue_push(&q, VECTOR(*vids)[i]);
    dqueue_push(&q, 0);
    nodes_reached=1;
    already_counted[(long int)VECTOR(*vids)[i]]=i+1;
    
    while (!dqueue_empty(&q)) {
      long int act=dqueue_pop(&q);
      long int actdist=dqueue_pop(&q);
      VECTOR(*res)[i] += actdist;

      igraph_neighbors(graph, &tmp, act, mode);
      for (j=0; j<vector_size(&tmp); j++) {
	long int neighbor=VECTOR(tmp)[j];
	if (already_counted[neighbor] == i+1) { continue; }
	already_counted[neighbor] = i+1;
	nodes_reached++;
	dqueue_push(&q, neighbor);
	dqueue_push(&q, actdist+1);
      }
    }
    VECTOR(*res)[i] += (no_of_nodes * (no_of_nodes-nodes_reached));
    VECTOR(*res)[i] = (no_of_nodes-1) / VECTOR(*res)[i];
  }
  
  /* Clean */
  dqueue_destroy(&q);
  vector_destroy(&tmp);
  Free(already_counted);

  return 0;
}

/**
 * \ingroup structural
 * \brief The length of the shortest paths between vertices.
 *
 * @param graph The graph object.
 * @param res The result of the calculation, a matrix. It has the same
 *        number of rows as the length of the <code>from</code>
 *        argument, and its number of columns is the number of
 *        vertices in the graph. One row of the matrix shows the
 *        distances from/to a given vertex to all the others in the
 *        graph, the order is fixed by the vertex ids.
 * @param from Vector of the vertex ids for which the path length
 *        calculations are done.
 * @param mode The type of shortest paths to be use for the
 *        calculation in directed graphs. Possible values: 
 *        - <b>IGRAPH_OUT</b>, the lengths of the outgoing paths are
 *          calculated. 
 *        - <b>IGRAPH_IN</b>, the lengths of the incoming paths are
 *          calculated. 
 *        - <b>IGRAPH_ALL</b>, the directed graph is considered as an
 *          undirected one for the computation.
 * @return Error code.
 * 
 * Time complexity: <code>O(n(|V|+|E|))</code>, <code>n</code> is the
 * number of vertices to calculate, <code>|V|</code> and
 * <code>|E|</code> are the number of vertices and edges in the graph.
 *
 * \sa igraph_get_shortest_paths() to get the paths themselves.
 */

int igraph_shortest_paths(igraph_t *graph, matrix_t *res, 
			  vector_t *from, igraph_neimode_t mode) {

  long int no_of_nodes=igraph_vcount(graph);
  long int no_of_from=vector_size(from);
  long int *already_counted;
  
  dqueue_t q;

  long int i, j;
  vector_t tmp;

  already_counted=Calloc(no_of_nodes, long int);
  vector_init(&tmp, 0);
  matrix_resize(res, no_of_from, no_of_nodes);
  matrix_null(res);

  dqueue_init(&q, 100);

  for (i=0; i<no_of_from; i++) {
    long int reached=1;
    dqueue_push(&q, VECTOR(*from)[i]);
    dqueue_push(&q, 0);
    already_counted[ (long int) VECTOR(*from)[i] ] = i+1;
    
    while (!dqueue_empty(&q)) {
      long int act=dqueue_pop(&q);
      long int actdist=dqueue_pop(&q);
      MATRIX(*res, i, act)=actdist;
      
      igraph_neighbors(graph, &tmp, act, mode);
      for (j=0; j<vector_size(&tmp); j++) {
	long int neighbor=VECTOR(tmp)[j];
	if (already_counted[neighbor] == i+1) { continue; }
	already_counted[neighbor] = i+1;
	reached++;
	dqueue_push(&q, neighbor);
	dqueue_push(&q, actdist+1);
      }
    }
    /* Plus the unreachable nodes */
    j=0;
    while (reached < no_of_nodes) {
      if (MATRIX(*res, i, j) == 0 && j != VECTOR(*from)[i]) {
	MATRIX(*res, i, j)=no_of_nodes;
	reached++;
      }
      j++;
    }
  }

  /* Clean */
  vector_destroy(&tmp);
  Free(already_counted);
  dqueue_destroy(&q);

  return 0;
}

/**
 * \ingroup structural
 * \brief Calculates the shortest paths from/to one vertex.
 * 
 * If there is more than one geodesic between two vertices, this
 * function gives only one of them. 
 * @param graph The graph object.
 * @param res The result, this is a pointer array to vector
 *        objects. These should be initialized before passing them to
 *        the function, which will properly erase and/or resize them
 *        and fill the ids of the vertices along the geodesics from/to
 *        the vertices.
 * @param from The id of the vertex from/to which the geodesics are
 *        calculated. 
 * @param mode The type of shortest paths to be use for the
 *        calculation in directed graphs. Possible values: 
 *        - <b>IGRAPH_OUT</b>, the lengths of the outgoing paths are
 *          calculated. 
 *        - <b>IGRAPH_IN</b>, the lengths of the incoming paths are
 *          calculated. 
 *        - <b>IGRAPH_ALL</b>, the directed graph is considered as an
 *          undirected one for the computation.
 * @return Error code.
 * 
 * Time complexity: <code>O(|V|+|E|)</code>, <code>|V|</code> is the
 * number of vertices, <code>|E|</code> the number of edges in the
 * graph. 
 *
 * \sa igraph_shortest_paths() if you only need the path length but
 * not the paths themselves.
 */
 

int igraph_get_shortest_paths(igraph_t *graph, vector_t *res,
			      integer_t from, igraph_neimode_t mode) {

  long int no_of_nodes=igraph_vcount(graph);
  long int *father;
  
  dqueue_t q;

  long int i, j;
  vector_t tmp;
  
  father=Calloc(no_of_nodes, long int);
  vector_init(&tmp, 0);
  dqueue_init(&q, 100);

  dqueue_push(&q, from+1);
  father[ (long int)from ] = from+1;
  
  while (!dqueue_empty(&q)) {
    long int act=dqueue_pop(&q);
    
    igraph_neighbors(graph, &tmp, act-1, mode);
    for (j=0; j<vector_size(&tmp); j++) {
      long int neighbor=VECTOR(tmp)[j]+1;
      if (father[neighbor-1] != 0) { continue; }
      father[neighbor-1] = act;
      dqueue_push(&q, neighbor);
    }
  }
  
    
  for (j=0; j<no_of_nodes; j++) {
    vector_clear(&res[j]);
    if (father[j]!=0) {
      long int act=j+1;
      long int size=0;
      while (father[act-1] != act) {
	size++;
	act=father[act-1];
      }
      size++;
      vector_resize(&res[j], size);
      VECTOR(res[j])[--size]=j;
      act=j+1;
      while (father[act-1] != act) {
	VECTOR(res[j])[--size]=father[act-1]-1;
	act=father[act-1];
      }
    }
  }
  
  /* Clean */
  Free(father);
  dqueue_destroy(&q);
  vector_destroy(&tmp);
  
  return 0;
}

/** 
 * \ingroup structural
 * \brief The vertices in the same component as a given vertex.
 *
 * @param graph The graph object.
 * @param res The result, vector with the ids of the vertices in the
 *        same component. 
 * @param vertex The id of the vertex of which the component is
 *        searched. 
 * @param mode Type of the component for directed graphs, possible values:
 *        - <b>IGRAPH_OUT</b>, the set of vertices reachable
 *          <em>from</em> the <code>vertex</code>,
 *        - <b>IGRAPH_IN</b>, the set of vertices from which the
 *          <code>vertex</code> is reachable.
 *        - <b>IGRAPH_ALL</b>, the graph is considered as an
 *          undirected graph. Note that this is <em>not</em> the same
 *          as the union of the previous two.
 * @return Error code.
 * 
 * Time complexity: <code>O(|V|+|E|)</code>, <code><|V|</code> and
 * <code>|E|</code> are the number of vertices and edges in the graph.
 * 
 * \sa igraph_subgraph() if you want a graph object consisting only
 * a given set of vertices and the edges between them.
 */

int igraph_subcomponent(igraph_t *graph, vector_t *res, real_t vertex, 
			igraph_neimode_t mode) {

  long int no_of_nodes=igraph_vcount(graph);
  dqueue_t q;
  char *already_added;
  long int i,j;
  vector_t tmp;
  
  vector_init(&tmp, 0);
  already_added=Calloc(no_of_nodes, char);
  dqueue_init(&q, 100);
  
  dqueue_push(&q, vertex);
  vector_push_back(res, vertex);
  already_added[(long int)vertex]=1;
  
  while (!dqueue_empty(&q)) {
    long int actnode=dqueue_pop(&q);

    igraph_neighbors(graph, &tmp, actnode, mode);
    for (i=0; i<vector_size(&tmp); i++) {
      long int neighbor=VECTOR(tmp)[i];
      
      if (already_added[neighbor]) { continue; }
      already_added[neighbor]=1;
      vector_push_back(res, neighbor);
      dqueue_push(&q, neighbor);
    }
  }

  dqueue_destroy(&q);
  vector_destroy(&tmp);

  return 0;
}

/**
 * \ingroup structural
 * \brief Betweenness centrality of some vertices.
 * 
 * The betweenness centrality of a vertex is the number of geodesics
 * going through it. If there are more than one geodesics between two
 * vertices, the value of these geodesics are weighted by one over the 
 * number of geodesics.
 * @param graph The graph object.
 * @param res The result of the computation, vector containing the
 *        betweenness scores for the specified vertices.
 * @param vids The vertices of which the betweenness centrality scores
 *        will be calculated.
 * @param directed Logical, if true directed paths will be considered
 *        for directed graphs. It is ignored for undirected graphs.
 * @return Error code.
 *
 * Time complexity: <code>O(|V||E|)</code>, <code>|V|</code> and
 * <code>|E|</code> are the number of vertices and edges in the graph.
 * Note that the time complexy is independent of the number of
 * vertices for which the score is calculated.
 *
 * \sa Other centrality types: igraph_degree(), igraph_closeness().
 *     See igraph_edge_betweenness() for calculating the betweenness score
 *     of the edges in a graph.
 */

int igraph_betweenness (igraph_t *graph, vector_t *res, vector_t *vids, 
			bool_t directed) {

  long int no_of_nodes=igraph_vcount(graph);
  dqueue_t q;
  long int *distance;
  long int *nrgeo;
  double *tmpscore;
  stack_t stack;
  long int source;
  long int j;
  vector_t tmp;
  integer_t modein, modeout;
  real_t *tmpres;

  if (directed) 
    { modeout=IGRAPH_OUT; modein=IGRAPH_IN; } 
  else 
    { modeout=modein=IGRAPH_ALL; }

  distance=Calloc(no_of_nodes, long int);
  nrgeo=Calloc(no_of_nodes, long int);
  tmpscore=Calloc(no_of_nodes, double);
  tmpres=Calloc(no_of_nodes, real_t);

  vector_init(&tmp, 0);
  vector_resize(res, vector_size(vids));
  vector_null(res);

  dqueue_init(&q, 100);
  stack_init(&stack, no_of_nodes);

  /* here we go */
  
  for (source=0; source<no_of_nodes; source++) {
    
    memset(distance, 0, no_of_nodes*sizeof(long int));
    memset(nrgeo, 0, no_of_nodes*sizeof(long int));
    memset(tmpscore, 0, no_of_nodes*sizeof(double));
    stack_clear(&stack); /* it should be empty anyway... */
    
    dqueue_push(&q, source);
    nrgeo[source]=1;
    distance[source]=0;
    
    while (!dqueue_empty(&q)) {
      long int actnode=dqueue_pop(&q);

      igraph_neighbors(graph, &tmp, actnode, modeout);
      for (j=0; j<vector_size(&tmp); j++) {
	long int neighbor=VECTOR(tmp)[j];
	if (nrgeo[neighbor] != 0) {
	  /* we've already seen this node, another shortest path? */
	  if (distance[neighbor]==distance[actnode]+1) {
	    nrgeo[neighbor]++;
	  }
	} else {
	  /* we haven't seen this node yet */
	  nrgeo[neighbor]++;
	  distance[neighbor]=distance[actnode]+1;
	  dqueue_push(&q, neighbor);
	  stack_push(&stack, neighbor);
	}
      }
    } /* while !dqueue_empty */
    
    /* Ok, we've the distance of each node and also the number of
       shortest paths to them. Now we do an inverse search, starting
       with the farthest nodes. */
    while (!stack_empty(&stack)) {
      long int actnode=stack_pop(&stack);
      long int friends=0;
      if (distance[actnode]<=1) { continue; } /* skip source node */
      
      /* search for the neighbors on the geodesics */
      igraph_neighbors(graph, &tmp, actnode, modein);
      for (j=0; j<vector_size(&tmp); j++) {
	long int neighbor=VECTOR(tmp)[j];
	if (distance[neighbor]==distance[actnode]-1) { friends++; }
      }

      /* set the temporary score of the friends */
      for (j=0; j<vector_size(&tmp); j++) {
	long int neighbor=VECTOR(tmp)[j];
	if (distance[neighbor]==distance[actnode]-1) {
	  tmpscore[neighbor] += (tmpscore[actnode]+1)/friends;
	}
      }
    }
    
    /* Ok, we've the scores for this source */
    for (j=0; j<vector_size(vids); j++) {
      long int node=VECTOR(*vids)[j];
      VECTOR(*res)[node] += tmpscore[node];
      tmpscore[node] = 0.0; /* in case a node is in vids multiple times */
    }

  } /* for source < no_of_nodes */

  /* divide by 2 for undirected graph */
  if (!directed || !igraph_is_directed(graph)) {
    for (j=0; j<vector_size(res); j++) {
      VECTOR(*res)[j] /= 2.0;
    }
  }
  
  /* clean  */
  Free(distance);
  Free(nrgeo);
  Free(tmpscore);
  dqueue_destroy(&q);
  stack_destroy(&stack);
  vector_destroy(&tmp);

  return 0;
}

/**
 * \ingroup structural
 * \brief Betweenness centrality of the edges.
 * 
 * The betweenness centrality of an edge is the number of geodesics
 * going through it. If there are more than one geodesics between two
 * vertices, the value of these geodesics are weighted by one over the 
 * number of geodesics.
 * @param graph The graph object.
 * @param result The result of the computation, vector containing the
 *        betweenness scores for the edges.
 * @param directed Logical, if true directed paths will be considered
 *        for directed graphs. It is ignored for undirected graphs.
 * @return Error code.
 *
 * Time complexity: <code>O(|V||E|)</code>, <code>|V|</code> and
 * <code>|E|</code> are the number of vertices and edges in the graph.
 *
 * \sa Other centrality types: igraph_degree(), igraph_closeness().
 *     See igraph_edge_betweenness() for calculating the betweenness score
 *     of the edges in a graph.
 */

int igraph_edge_betweenness (igraph_t *graph, vector_t *result, 
			     bool_t directed) {
  
  long int no_of_nodes=igraph_vcount(graph);
  long int no_of_edges=igraph_ecount(graph);
  dqueue_t q;
  long int *distance;
  long int *nrgeo;
  double *tmpscore;
  stack_t stack;
  long int source;
  long int j;

  igraph_iterator_t it;
  integer_t modein, modeout;

  if (directed) {
    modeout=IGRAPH_OUT;
    modein=IGRAPH_IN;
  } else {
    modeout=modein=IGRAPH_ALL;
  }
  
  distance=Calloc(no_of_nodes, long int);
  nrgeo=Calloc(no_of_nodes, long int);
  tmpscore=Calloc(no_of_nodes, double);

  vector_resize(result, no_of_edges);
  vector_null(result);

  dqueue_init(&q, 100);
  stack_init(&stack, no_of_nodes);
  igraph_iterator_eneis(graph, &it, 0, modeout);

  /* here we go */
  
  for (source=0; source<no_of_nodes; source++) {

    memset(distance, 0, no_of_nodes*sizeof(long int));
    memset(nrgeo, 0, no_of_nodes*sizeof(long int));
    memset(tmpscore, 0, no_of_nodes*sizeof(double));
    stack_clear(&stack); /* it should be empty anyway... */
    
    dqueue_push(&q, source);
    nrgeo[source]=1;
    distance[source]=0;
    
    while (!dqueue_empty(&q)) {
      long int actnode=dqueue_pop(&q);
    
      igraph_iterator_eneis_set(graph, &it, actnode, modeout);
      while ( !igraph_end(graph, &it)) {
	long int neighbor=igraph_get_vertex_nei(graph, &it);
	if (nrgeo[neighbor] != 0) {
	  /* we've already seen this node, another shortest path? */
	  if (distance[neighbor]==distance[actnode]+1) {
	    nrgeo[neighbor]++;
	  }
	} else {
	  /* we haven't seen this node yet */
	  nrgeo[neighbor]++;
	  distance[neighbor]=distance[actnode]+1;
	  dqueue_push(&q, neighbor);
	  stack_push(&stack, neighbor);
	}
	igraph_next(graph, &it);
      }
    } /* while !dqueue_empty */
    
    /* Ok, we've the distance of each node and also the number of
       shortest paths to them. Now we do an inverse search, starting
       with the farthest nodes. */
    while (!stack_empty(&stack)) {
      long int actnode=stack_pop(&stack);
      long int friends=0;
      if (distance[actnode]<1) { continue; } /* skip source node */
      
      /* search for the neighbors on the geodesics */
      igraph_iterator_eneis_set(graph, &it, actnode, modein);
      while ( !igraph_end(graph,  &it)) {
	long int neighbor=igraph_get_vertex_nei(graph, &it);
	if (distance[neighbor]==distance[actnode]-1) { friends++; }
	igraph_next(graph, &it);
      }

      /* set the temporary score of the friends */
      igraph_iterator_eneis_set(graph, &it, actnode, modein);
      while ( !igraph_end(graph,  &it)) {
	long int neighbor=igraph_get_vertex_nei(graph, &it);
	long int edgeno=igraph_get_edge(graph, &it);
	if (distance[neighbor]==distance[actnode]-1) {
	  tmpscore[neighbor] += (tmpscore[actnode]+1)/friends;
	  VECTOR(*result)[edgeno] += (tmpscore[actnode]+1)/friends;
	}
	igraph_next(graph, &it);	
      }
    }
    /* Ok, we've the scores for this source */
  } /* for source <= no_of_nodes */
  
  /* clean and return */
  igraph_iterator_destroy(graph, &it);
  Free(distance);
  Free(nrgeo);
  Free(tmpscore);
  dqueue_destroy(&q);
  stack_destroy(&stack);

  /* divide by 2 for undirected graph */
  if (!directed || !igraph_is_directed(graph)) {
    for (j=0; j<vector_size(result); j++) {
      VECTOR(*result)[j] /= 2.0;
    }
  }
  
  return 0;
}

/**
 * \ingroup structural
 * \brief Creates a subgraph with the specified vertices.
 * 
 * This function collects the specified vertices and all edges between
 * them to a new graph.
 * As the vertex ids in a graph always start with one, this function
 * very likely needs to reassign ids to the vertices.
 * @param graph The graph object.
 * @param res The subgraph, another graph object will be stored here,
 *        do <em>not</em> initialize this object before calling this
 *        function, and call igraph_destroy() on it if you don't need
 *        it any more.
 * @param vids Vector with the vertex ids to put in the subgraph.
 * @return Error code.
 * 
 * Time complexity: <code>O(|V|+|E|)</code>, <code>|V|</code> and
 * <code>|E|</code> are the number of vertices and edges in the
 * original graph.
 *
 * \sa igraph_delete_vertices() to delete the specified set of
 * vertices from a graph, the opposite of this function.
 */

int igraph_subgraph(igraph_t *graph, igraph_t *res, vector_t *vids) {
  
  long int no_of_nodes=igraph_vcount(graph);
  vector_t delete;
  char *remain;
  long int i;

  vector_init(&delete, 0);
  vector_reserve(&delete, no_of_nodes-vector_size(vids));
  remain=Calloc(no_of_nodes, char);
  
  for (i=0; i<vector_size(vids); i++) {
    remain[ (long int) VECTOR(*vids)[i] ] = 1;
  }

  for (i=0; i<no_of_nodes; i++) {
    if (remain[i] == 0) {
      vector_push_back(&delete, i);
    }
  }

  Free(remain);
  igraph_copy(res, graph);
  igraph_delete_vertices(res, &delete);
  
  vector_destroy(&delete);
  return 0;
}
