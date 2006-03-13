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
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
   02110-1301 USA

*/

#include "igraph.h"
#include "memory.h"
#include "random.h"

#include <string.h>

/** 
 * \section about_structural
 *
 * <para>These functions usually calculate some structural property
 * of a graph, like its diameter, the degree of the nodes, etc.</para>
 */

/**
 * \ingroup structural
 * \function igraph_diameter
 * \brief Calculates the diameter of a graph (longest geodesic).
 *
 * \param graph The graph object.
 * \param res Pointer to an integer, this will contain the result.
 * \param directed Boolean, whether to consider directed
 *        paths. Ignored for undirected graphs.
 * \param unconn What to do if the graph is not connected. If
 *        \c TRUE the longest geodesic within a component
 *        will be returned, otherwise the number of vertices is
 *        returned. (The rationale behind the latter is that this is
 *        always longer than the longest possible diameter in a
 *        graph.) 
 * \return Error code:
 *         \c IGRAPH_ENOMEM, not enough memory for
 *         temporary data.  
 *
 * Time complexity: O(|V||E|), the
 * number of vertices times the number of edges.
 */

int igraph_diameter(const igraph_t *graph, integer_t *pres, 
		    integer_t *pfrom, integer_t *pto, 
		    igraph_vector_t *path,
		    bool_t directed, bool_t unconn) {

  long int no_of_nodes=igraph_vcount(graph);
  long int i, j, n;
  long int *already_added;
  long int nodes_reached;
  long int from=0, to=0;
  long int res=0;

  igraph_dqueue_t q=IGRAPH_DQUEUE_NULL;
  igraph_vector_t *neis;
  integer_t dirmode;
  igraph_i_adjlist_t allneis;
  
  if (directed) { dirmode=IGRAPH_OUT; } else { dirmode=IGRAPH_ALL; }
  already_added=Calloc(no_of_nodes, long int);
  if (already_added==0) {
    IGRAPH_ERROR("diameter failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, already_added);
  IGRAPH_DQUEUE_INIT_FINALLY(&q, 100);
  
  igraph_i_adjlist_init(graph, &allneis, dirmode);
  IGRAPH_FINALLY(igraph_i_adjlist_destroy, &allneis);
  
  for (i=0; i<no_of_nodes; i++) {
    nodes_reached=1;
    IGRAPH_CHECK(igraph_dqueue_push(&q, i));
    IGRAPH_CHECK(igraph_dqueue_push(&q, 0));
    already_added[i]=i+1;
    
    while (!igraph_dqueue_empty(&q)) {
      long int actnode=igraph_dqueue_pop(&q);
      long int actdist=igraph_dqueue_pop(&q);
      if (actdist>res) { 
	res=actdist; 
	from=i;
	to=actnode;
      }
      
      neis=igraph_i_adjlist_get(&allneis, actnode);
      n=igraph_vector_size(neis);
      for (j=0; j<n; j++) {
	long int neighbor=VECTOR(*neis)[j];
	if (already_added[neighbor] == i+1) { continue; }
	already_added[neighbor]=i+1;
	nodes_reached++;
	IGRAPH_CHECK(igraph_dqueue_push(&q, neighbor));
	IGRAPH_CHECK(igraph_dqueue_push(&q, actdist+1));
      }
    } /* while !igraph_dqueue_empty */
    
    /* not connected, return largest possible */
    if (nodes_reached != no_of_nodes && !unconn) {
      res=no_of_nodes;
      from=-1;
      to=-1;
      break;
    }
  } /* for i<no_of_nodes */

  /* return the requested info */
  if (pres != 0) {
    *pres=res;
  }
  if (pfrom != 0) {
    *pfrom=from;
  }
  if (pto != 0) {
    *pto=to;
  }
  if (path != 0) {
    if (res==no_of_nodes) {
      igraph_vector_clear(path);
    } else {
      igraph_vector_ptr_t tmpptr;
      igraph_vector_ptr_init(&tmpptr, 1);
      IGRAPH_FINALLY(igraph_vector_ptr_destroy, &tmpptr);
      VECTOR(tmpptr)[0]=path;
      IGRAPH_CHECK(igraph_get_shortest_paths(graph, &tmpptr, from, 
					     IGRAPH_VS_1(graph, to), dirmode));
      igraph_vector_ptr_destroy(&tmpptr);
      IGRAPH_FINALLY_CLEAN(1);
    }
  }
  
  /* clean */
  Free(already_added);
  igraph_dqueue_destroy(&q);
  igraph_i_adjlist_destroy(&allneis);
  IGRAPH_FINALLY_CLEAN(3);

  return 0;
}

/**
 * \ingroup structural
 * \function igraph_average_path_length
 * \brief Calculates the average geodesic length in a graph.
 *
 * \param graph The graph object.
 * \param res Pointer to a real number, this will contain the result.
 * \param directed Boolean, whether to consider directed
 *        paths. Ignored for undirected graphs.
 * \param unconn What to do if the graph is not connected. If
 *        \c TRUE the average of thr geodesics
 *        within the components 
 *        will be returned, otherwise the number of vertices is
 *        used for the length of non-existing geodesics. (The rationale
 *        behind this is that this is always longer than the longest
 *        possible geodesic in a graph.) 
 * \return Error code:
 *         \c IGRAPH_ENOMEM, not enough memory for
 *         data structures 
 *
 * Time complexity: O(|V||E|), the
 * number of vertices times the number of edges.
 */

int igraph_average_path_length(const igraph_t *graph, real_t *res,
			       bool_t directed, bool_t unconn) {
  long int no_of_nodes=igraph_vcount(graph);
  long int i, j, n;
  long int *already_added;
  long int nodes_reached=0;
  real_t normfact=0.0;

  igraph_dqueue_t q=IGRAPH_DQUEUE_NULL;
  igraph_vector_t *neis;
  integer_t dirmode;
  igraph_i_adjlist_t allneis;

  *res=0;  
  if (directed) { dirmode=IGRAPH_OUT; } else { dirmode=IGRAPH_ALL; }
  already_added=Calloc(no_of_nodes, long int);
  if (already_added==0) {
    IGRAPH_ERROR("average path length failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(free, already_added); /* TODO: hack */
  IGRAPH_DQUEUE_INIT_FINALLY(&q, 100);

  igraph_i_adjlist_init(graph, &allneis, dirmode);
  IGRAPH_FINALLY(igraph_i_adjlist_destroy, &allneis);

  for (i=0; i<no_of_nodes; i++) {
    nodes_reached=0;
    IGRAPH_CHECK(igraph_dqueue_push(&q, i));
    IGRAPH_CHECK(igraph_dqueue_push(&q, 0));
    already_added[i]=i+1;
    
    while (!igraph_dqueue_empty(&q)) {
      long int actnode=igraph_dqueue_pop(&q);
      long int actdist=igraph_dqueue_pop(&q);
    
      neis=igraph_i_adjlist_get(&allneis, actnode);
      n=igraph_vector_size(neis);
      for (j=0; j<n; j++) {
	long int neighbor=VECTOR(*neis)[j];
	if (already_added[neighbor] == i+1) { continue; }
	already_added[neighbor]=i+1;
	nodes_reached++;
	*res += actdist+1;
	normfact+=1;
	IGRAPH_CHECK(igraph_dqueue_push(&q, neighbor));
	IGRAPH_CHECK(igraph_dqueue_push(&q, actdist+1));
      }
    } /* while !igraph_dqueue_empty */
    
    /* not connected, return largest possible */
    if (!unconn) {
      *res += (no_of_nodes * (no_of_nodes-1-nodes_reached));
      normfact += no_of_nodes-1-nodes_reached;
    }    
  } /* for i<no_of_nodes */

  *res /= normfact;

  /* clean */
  Free(already_added);
  igraph_dqueue_destroy(&q);
  igraph_i_adjlist_destroy(&allneis);
  IGRAPH_FINALLY_CLEAN(3);

  return 0;
}

/**
 * \ingroup structural
 * \function igraph_minimum_spanning_tree_unweighted
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
 * \param graph The graph object.
 * \param mst The minimum spanning tree, another graph object. Do
 *        \em not initialize this object before passing it to
 *        this function, but be sure to call \ref igraph_destroy() on it if
 *        you don't need it any more.
 * \return Error code:
 *         \c IGRAPH_ENOMEM, not enough memory for
 *         temporary data. 
 *
 * Time complexity: O(|V|+|E|),
 * |V| is the 
 * number of vertices, |E| the number
 * of edges in the graph. 
 *
 * \sa \ref igraph_minimum_spanning_tree_prim() for weighted graphs.
 */

int igraph_minimum_spanning_tree_unweighted(const igraph_t *graph, 
					    igraph_t *mst) {

  long int no_of_nodes=igraph_vcount(graph);
  char *already_added;
  
  igraph_dqueue_t q=IGRAPH_DQUEUE_NULL;
  igraph_vector_t edges=IGRAPH_VECTOR_NULL;
  igraph_vector_t tmp=IGRAPH_VECTOR_NULL;
  long int i, j;

  IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
  already_added=Calloc(no_of_nodes, char);
  if (already_added==0) {
    IGRAPH_ERROR("unweighted spanning tree failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(free, already_added); /* TODO: hack */
  IGRAPH_VECTOR_INIT_FINALLY(&tmp, 0);
  IGRAPH_DQUEUE_INIT_FINALLY(&q, 100);
  IGRAPH_CHECK(igraph_vector_reserve(&edges, (no_of_nodes-1)*2));
  
  for (i=0; i<no_of_nodes; i++) {
    if (already_added[i]>0) { continue; }

    already_added[i]=1;
    IGRAPH_CHECK(igraph_dqueue_push(&q, i));
    while (! igraph_dqueue_empty(&q)) {
      long int act_node=igraph_dqueue_pop(&q);
      IGRAPH_CHECK(igraph_neighbors(graph, &tmp, act_node, 3));
      for (j=0; j<igraph_vector_size(&tmp); j++) {
	long int neighbor=VECTOR(tmp)[j];
	if (already_added[neighbor]==0) {
	  already_added[neighbor]=1;
	  IGRAPH_CHECK(igraph_dqueue_push(&q, neighbor));
	  IGRAPH_CHECK(igraph_vector_push_back(&edges, act_node));
	  IGRAPH_CHECK(igraph_vector_push_back(&edges, neighbor));
	}
      }
    }
  }
  
  igraph_dqueue_destroy(&q);
  Free(already_added);
  igraph_vector_destroy(&tmp);
  IGRAPH_FINALLY_CLEAN(3);

  IGRAPH_CHECK(igraph_create(mst, &edges, 0, igraph_is_directed(graph)));
  igraph_vector_destroy(&edges);
  IGRAPH_FINALLY_CLEAN(1);

  return 0;
}

/**
 * \ingroup structural
 * \function igraph_minimum_spanning_tree_prim
 * \brief Calculates one minimum spanning tree of a weighted graph.
 *
 * This function uses Prim's method for carrying out the computation,
 * see Prim, R.C.: Shortest connection networks and some
 * generalizations, Bell System Technical
 * Journal, Vol. 36, 
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
 * \param graph The graph object.
 * \param mst The result of the computation, a graph object containing
 *        the minimum spanning tree of the graph.
 *        Do \em not initialize this object before passing it to
 *        this function, but be sure to call \ref igraph_destroy() on it if
 *        you don't need it any more.
 * \param weights A vector containing the weights of the the edges.
 *        in the same order as the simple edge iterator visits them.
 * \return Error code:
 *         \c IGRAPH_ENOMEM, not enough memory.
 *         \c IGRAPH_EINVAL, length of weight vector does not
 *           match number of edges.
 *
 * Time complexity: O(|V|+|E|),
 * |V| is the number of vertices,
 * |E| the number of edges in the 
 * graph. 
 *
 * \sa \ref igraph_minimum_spanning_tree_unweighted() for unweighted graphs.
 */

int igraph_minimum_spanning_tree_prim(const igraph_t *graph, igraph_t *mst,
				      const igraph_vector_t *weights) {

  long int no_of_nodes=igraph_vcount(graph);
  char *already_added;

  igraph_d_indheap_t heap=IGRAPH_D_INDHEAP_NULL;
  igraph_vector_t edges=IGRAPH_VECTOR_NULL;
  integer_t mode=3;
  
  igraph_es_t it;

  long int i;

  if (igraph_vector_size(weights) != igraph_ecount(graph)) {
    IGRAPH_ERROR("Invalid weights length", IGRAPH_EINVAL);
  }

  IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
  already_added=Calloc(no_of_nodes, char);
  if (already_added == 0) {
    IGRAPH_ERROR("prim spanning tree failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(free, already_added); /* TODO: hack */
  IGRAPH_CHECK(igraph_d_indheap_init(&heap, 0));
  IGRAPH_FINALLY(igraph_d_indheap_destroy, &heap);
  IGRAPH_CHECK(igraph_vector_reserve(&edges, (no_of_nodes-1)*2));
  IGRAPH_CHECK(igraph_es_adj(graph, &it, 0, mode));
  IGRAPH_FINALLY(igraph_es_destroy, &it);

  for (i=0; i<no_of_nodes; i++) {
    if (already_added[i]>0) { continue; }

    already_added[i]=1;
    /* add all edges of the first vertex */
    igraph_es_adj_set(graph, &it, i, mode);
    while( !igraph_es_end(graph, &it)) {
      long int neighbor=igraph_es_adj_vertex(graph, &it);
      long int edgeno=igraph_es_get(graph, &it);
      if (already_added[neighbor] == 0) {
	IGRAPH_CHECK(igraph_d_indheap_push(&heap, -VECTOR(*weights)[edgeno], i, 
				    neighbor));
      }
      igraph_es_next(graph, &it);
    }

    while(! igraph_d_indheap_empty(&heap)) {
      /* Get minimal edge */
      long int from, to;
      igraph_d_indheap_max_index(&heap, &from, &to);

      /* Erase it */
      igraph_d_indheap_delete_max(&heap);
      
      /* Does it point to a visited node? */
      if (already_added[to]==0) {
	already_added[to]=1;
	IGRAPH_CHECK(igraph_vector_push_back(&edges, from));
	IGRAPH_CHECK(igraph_vector_push_back(&edges, to));
	/* add all outgoing edges */
	igraph_es_adj_set(graph, &it, to, mode);
	while ( !igraph_es_end(graph, &it)) {
	  long int neighbor=igraph_es_adj_vertex(graph, &it);
	  long int edgeno=igraph_es_get(graph, &it);
	  if (already_added[neighbor] == 0) {
	    IGRAPH_CHECK(igraph_d_indheap_push(&heap, -VECTOR(*weights)[edgeno], to, 
					neighbor));
	  }
	  igraph_es_next(graph, &it);
	} /* for */
      } /* if !already_added */
    } /* while in the same component */
  } /* for all nodes */

  igraph_d_indheap_destroy(&heap);
  Free(already_added);
  IGRAPH_FINALLY_CLEAN(2);

  IGRAPH_CHECK(igraph_create(mst, &edges, 0, igraph_is_directed(graph)));
  igraph_vector_destroy(&edges);
  igraph_es_destroy(&it);
  IGRAPH_FINALLY_CLEAN(2);
  
  return 0;
}

/**
 * \ingroup structural
 * \function igraph_closeness
 * \brief Closeness centrality calculations for some vertices.
 *
 * The closeness centrality of a vertex measures how easily other
 * vertices can be reached from it (or the other way: how easily it
 * can be reached from the other vertices). It is defined as the
 * number of the number of vertices minus one divided by the sum of the
 * lengths of all geodesics from/to the given vertex.
 *
 * If the graph is not connected, and there is no path between two
 * vertices, the number of vertices is used instead the length of the
 * geodesic. This is always longer than the longest possible geodesic.
 * 
 * \param graph The graph object.
 * \param res The result of the computation, a vector containing the
 *        closeness centrality scores for the given vertices.
 * \param vids Vector giving the vertices for which the closeness
 *        centrality scores will be computed.
 * \param mode The type of shortest paths to be used for the
 *        calculation in directed graphs. Possible values: 
 *        \clist
 *        \cli IGRAPH_OUT 
 *          the lengths of the outgoing paths are calculated. 
 *        \cli IGRAPH_IN 
 *          the lengths of the incoming paths are calculated. 
 *        \cli IGRAPH_ALL
 *          the directed graph is considered as an
 *          undirected one for the computation.
 *        \endclist
 * \return Error code:
 *        \clist
 *        \cli IGRAPH_ENOMEM
 *           not enough memory for temporary data.
 *        \cli IGRAPH_EINVVID
 *           invalid vertex id passed.
 *        \cli IGRAPH_EINVMODE
 *           invalid mode argument.
 *        \endclist
 *
 * Time complexity: O(n|E|),
 * n is the number 
 * of vertices for which the calculation is done and
 * |E| is the number 
 * of edges in the graph.
 *
 * \sa Other centrality types: \ref igraph_degree(), \ref igraph_betweenness().
 */

int igraph_closeness(const igraph_t *graph, igraph_vector_t *res, 
		     const igraph_vs_t *vids, 
		     igraph_neimode_t mode) {

  long int no_of_nodes=igraph_vcount(graph);
  igraph_vector_t already_counted;
  long int i, j;
  long int nodes_reached;

  igraph_dqueue_t q;
  
  long int nodes_to_calc;
  igraph_vector_t tmp;
  igraph_vs_t myvids;
  const igraph_vector_t *myvidsv;

  IGRAPH_CHECK(igraph_vs_vectorview_it(graph, vids, &myvids));
  IGRAPH_FINALLY(igraph_vs_destroy, &myvids);
  myvidsv=igraph_vs_vector_getvector(graph, &myvids);

  nodes_to_calc=igraph_vector_size(myvidsv);
  
  if (mode != IGRAPH_OUT && mode != IGRAPH_IN && 
      mode != IGRAPH_ALL) {
    IGRAPH_ERROR("calculating closeness", IGRAPH_EINVMODE);
  }
  if (!igraph_vector_isininterval(myvidsv, 0, no_of_nodes-1)) {
    IGRAPH_ERROR("calculating closeness", IGRAPH_EINVVID);
  }

  IGRAPH_VECTOR_INIT_FINALLY(&already_counted, no_of_nodes);
  IGRAPH_VECTOR_INIT_FINALLY(&tmp, 0);
  IGRAPH_DQUEUE_INIT_FINALLY(&q, 100);

  IGRAPH_CHECK(igraph_vector_resize(res, nodes_to_calc));
  igraph_vector_null(res);
  
  for (i=0; i<nodes_to_calc; i++) {
    IGRAPH_CHECK(igraph_dqueue_push(&q, VECTOR(*myvidsv)[i]));
    IGRAPH_CHECK(igraph_dqueue_push(&q, 0));
    nodes_reached=1;
    VECTOR(already_counted)[(long int)VECTOR(*myvidsv)[i]]=i+1;
    
    while (!igraph_dqueue_empty(&q)) {
      long int act=igraph_dqueue_pop(&q);
      long int actdist=igraph_dqueue_pop(&q);
      VECTOR(*res)[i] += actdist;

      IGRAPH_CHECK(igraph_neighbors(graph, &tmp, act, mode));
      for (j=0; j<igraph_vector_size(&tmp); j++) {
	long int neighbor=VECTOR(tmp)[j];
	if (VECTOR(already_counted)[neighbor] == i+1) { continue; }
	VECTOR(already_counted)[neighbor] = i+1;
	nodes_reached++;
	IGRAPH_CHECK(igraph_dqueue_push(&q, neighbor));
	IGRAPH_CHECK(igraph_dqueue_push(&q, actdist+1));
      }
    }
    VECTOR(*res)[i] += (no_of_nodes * (no_of_nodes-nodes_reached));
    VECTOR(*res)[i] = (no_of_nodes-1) / VECTOR(*res)[i];
  }
  
  /* Clean */
  igraph_dqueue_destroy(&q);
  igraph_vector_destroy(&tmp);
  igraph_vector_destroy(&already_counted);
  igraph_vs_destroy(&myvids);
  IGRAPH_FINALLY_CLEAN(4);
  
  if (vids->shorthand) { igraph_vs_destroy((igraph_vs_t*)vids); }
  return 0;
}

/**
 * \ingroup structural
 * \function igraph_shortest_paths
 * \brief The length of the shortest paths between vertices.
 *
 * \param graph The graph object.
 * \param res The result of the calculation, a matrix. It has the same
 *        number of rows as the length of the \c from
 *        argument, and its number of columns is the number of
 *        vertices in the graph. One row of the matrix shows the
 *        distances from/to a given vertex to all the others in the
 *        graph, the order is fixed by the vertex ids.
 * \param from Vector of the vertex ids for which the path length
 *        calculations are done.
 * \param mode The type of shortest paths to be use for the
 *        calculation in directed graphs. Possible values: 
 *        \clist
 *        \cli IGRAPH_OUT 
 *          the lengths of the outgoing paths are calculated. 
 *        \cli IGRAPH_IN 
 *          the lengths of the incoming paths are calculated. 
 *        \cli IGRAPH_ALL 
 *          the directed graph is considered as an undirected one for
 *          the computation. 
 *        \endclist
 * \return Error code:
 *        \clist
 *        \cli IGRAPH_ENOMEM 
 *           not enough memory for temporary
 *           data.
 *        \cli IGRAPH_EINVVID
 *           invalid vertex id passed.
 *        \cli IGRAPH_EINVMODE 
 *           invalid mode argument.
 *        \endclist
 * 
 * Time complexity: O(n(|V|+|E|)),
 * n is the 
 * number of vertices to calculate, |V| and
 * |E| are the number of vertices and
 * edges in the graph. 
 *
 * \sa \ref igraph_get_shortest_paths() to get the paths themselves.
 */

int igraph_shortest_paths(const igraph_t *graph, igraph_matrix_t *res, 
			  const igraph_vs_t *from, igraph_neimode_t mode) {

  long int no_of_nodes=igraph_vcount(graph);
  long int no_of_from;
  long int *already_counted;
  
  igraph_dqueue_t q=IGRAPH_DQUEUE_NULL;

  long int i, j;
  igraph_vector_t tmp=IGRAPH_VECTOR_NULL;
  igraph_vs_t myfrom;
  const igraph_vector_t *myfromv;

  IGRAPH_CHECK(igraph_vs_vectorview_it(graph, from, &myfrom));
  IGRAPH_FINALLY(igraph_vs_destroy, &myfrom);
  myfromv=igraph_vs_vector_getvector(graph, &myfrom);

  no_of_from=igraph_vector_size(myfromv);

  if (!igraph_vector_isininterval(myfromv, 0, no_of_nodes-1)) {
    IGRAPH_ERROR("shortest paths failed", IGRAPH_EINVVID);
  }
  if (mode != IGRAPH_OUT && mode != IGRAPH_IN && 
      mode != IGRAPH_ALL) {
    IGRAPH_ERROR("Invalid mode argument", IGRAPH_EINVMODE);
  }
  already_counted=Calloc(no_of_nodes, long int);
  if (already_counted==0) {
    IGRAPH_ERROR("shortest paths failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(free, already_counted);
  IGRAPH_VECTOR_INIT_FINALLY(&tmp, 0);
  IGRAPH_DQUEUE_INIT_FINALLY(&q, 100);

  IGRAPH_CHECK(igraph_matrix_resize(res, no_of_from, no_of_nodes));
  igraph_matrix_null(res);

  for (i=0; i<no_of_from; i++) {
    long int reached=1;
    IGRAPH_CHECK(igraph_dqueue_push(&q, VECTOR(*myfromv)[i]));
    IGRAPH_CHECK(igraph_dqueue_push(&q, 0));
    already_counted[ (long int) VECTOR(*myfromv)[i] ] = i+1;
    
    while (!igraph_dqueue_empty(&q)) {
      long int act=igraph_dqueue_pop(&q);
      long int actdist=igraph_dqueue_pop(&q);
      MATRIX(*res, i, act)=actdist;
      
      IGRAPH_CHECK(igraph_neighbors(graph, &tmp, act, mode));
      for (j=0; j<igraph_vector_size(&tmp); j++) {
	long int neighbor=VECTOR(tmp)[j];
	if (already_counted[neighbor] == i+1) { continue; }
	already_counted[neighbor] = i+1;
	reached++;
	IGRAPH_CHECK(igraph_dqueue_push(&q, neighbor));
	IGRAPH_CHECK(igraph_dqueue_push(&q, actdist+1));
      }
    }
    /* Plus the unreachable nodes */
    j=0;
    while (reached < no_of_nodes) {
      if (MATRIX(*res, i, j) == 0 && j != VECTOR(*myfromv)[i]) {
	MATRIX(*res, i, j)=no_of_nodes;
	reached++;
      }
      j++;
    }
  }

  /* Clean */
  igraph_vector_destroy(&tmp);
  Free(already_counted);
  igraph_dqueue_destroy(&q);
  igraph_vs_destroy(&myfrom);
  IGRAPH_FINALLY_CLEAN(4);

  if (from->shorthand) { igraph_vs_destroy((igraph_vs_t*)from); }
  return 0;
}

/**
 * \ingroup structural
 * \function igraph_get_shortest_paths
 * \brief Calculates the shortest paths from/to one vertex.
 * 
 * If there is more than one geodesic between two vertices, this
 * function gives only one of them. 
 * \param graph The graph object.
 * \param res The result, this is a pointer vector, each element points 
 *        to a vector
 *        object. These should be initialized before passing them to
 *        the function, which will properly clear and/or resize them
 *        and fill the ids of the vertices along the geodesics from/to
 *        the vertices.
 * \param from The id of the vertex from/to which the geodesics are
 *        calculated. 
 * \param to Vertex sequence with the ids of the vertices to/from which the 
 *        shortest paths will be calculated. A vertex might be given multiple
 *        times.
 * \param mode The type of shortest paths to be use for the
 *        calculation in directed graphs. Possible values: 
 *        \clist
 *        \cli IGRAPH_OUT 
 *          the outgoing paths are calculated. 
 *        \cli IGRAPH_IN 
 *          the incoming paths are calculated. 
 *        \cli IGRAPH_ALL 
 *          the directed graph is considered as an
 *          undirected one for the computation.
 *        \endclist
 * \return Error code:
 *        \clist
 *        \cli IGRAPH_ENOMEM 
 *           not enough memory for temporary data.
 *        \cli IGRAPH_EINVVID
 *           \p from is invalid vertex id, or the length of \p to is 
 *           not the same as the length of \p res.
 *        \cli IGRAPH_EINVMODE 
 *           invalid mode argument.
 *        \endclist
 * 
 * Time complexity: O(|V|+|E|),
 * |V| is the number of vertices,
 * |E| the number of edges in the
 * graph.  
 *
 * \sa \ref igraph_shortest_paths() if you only need the path length but
 * not the paths themselves.
 */
 

int igraph_get_shortest_paths(const igraph_t *graph, igraph_vector_ptr_t *res,
			      integer_t from, const igraph_vs_t *to, 
			      igraph_neimode_t mode) {

  /* TODO: use adjlist_t if to is long (longer than 1?) */

  long int no_of_nodes=igraph_vcount(graph);
  long int *father;
  
  igraph_dqueue_t q=IGRAPH_DQUEUE_NULL;

  long int j;
  igraph_vector_t tmp=IGRAPH_VECTOR_NULL;

  igraph_vs_t myto;
  const igraph_vector_t *mytov;
  
  long int to_reach;
  long int reached=0;

  if (from<0 || from>=no_of_nodes) {
    IGRAPH_ERROR("cannot get shortest paths", IGRAPH_EINVVID);
  }
  if (mode != IGRAPH_OUT && mode != IGRAPH_IN && 
      mode != IGRAPH_ALL) {
    IGRAPH_ERROR("Invalid mode argument", IGRAPH_EINVMODE);
  }

  IGRAPH_CHECK(igraph_vs_vectorview_it(graph, to, &myto));
  IGRAPH_FINALLY(igraph_vs_destroy, &myto);
  mytov=igraph_vs_vector_getvector(graph, &myto);

  if (igraph_vector_size(mytov) != igraph_vector_ptr_size(res)) {
    IGRAPH_ERROR("Size of the `res' and the `to' should match", IGRAPH_EINVAL);
  }

  father=Calloc(no_of_nodes, long int);
  if (father==0) {
    IGRAPH_ERROR("cannot get shortest paths", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(free, father);	/* TODO: hack */
  IGRAPH_VECTOR_INIT_FINALLY(&tmp, 0);
  IGRAPH_DQUEUE_INIT_FINALLY(&q, 100);

  /* Mark the vertices we need to reach */
  to_reach=igraph_vector_size(mytov);
  for (j=0; j<igraph_vector_size(mytov); j++) {
    if (father[ (long int) VECTOR(*mytov)[j] ] == 0) {
      father[ (long int) VECTOR(*mytov)[j] ] = -1;
    } else {
      to_reach--;		/* this node was given multiple times */
    }
  }

  IGRAPH_CHECK(igraph_dqueue_push(&q, from+1));
  if (father[ (long int) from ] < 0) { reached++; }
  father[ (long int)from ] = from+1;
  
  while (!igraph_dqueue_empty(&q) && reached < to_reach) {
    long int act=igraph_dqueue_pop(&q);
    
    IGRAPH_CHECK(igraph_neighbors(graph, &tmp, act-1, mode));
    for (j=0; j<igraph_vector_size(&tmp); j++) {
      long int neighbor=VECTOR(tmp)[j]+1;
      if (father[neighbor-1] > 0) { 
	continue; 
      } else if (father[neighbor-1] < 0) { 
	reached++; 
      }
      father[neighbor-1] = act;
      IGRAPH_CHECK(igraph_dqueue_push(&q, neighbor));
    }
  }

  if (reached < to_reach) {
    IGRAPH_WARNING("Couldn't reach some vertices");
  }
  
  for (j=0; j<igraph_vector_size(mytov); j++) {
    long int node=VECTOR(*mytov)[j];
    igraph_vector_t *vec=VECTOR(*res)[j];
    igraph_vector_clear(vec);
    if (father[node]>0) {
      long int act=node+1;
      long int size=0;
      while (father[act-1] != act) {
	size++;
	act=father[act-1];
      }
      size++;
      IGRAPH_CHECK(igraph_vector_resize(vec, size));
      VECTOR(*vec)[--size]=node;
      act=node+1;
      while (father[act-1] != act) {
	VECTOR(*vec)[--size]=father[act-1]-1;
	act=father[act-1];
      }
    }
  }
  
  /* Clean */
  Free(father);
  igraph_dqueue_destroy(&q);
  igraph_vector_destroy(&tmp);
  igraph_vs_destroy(&myto);
  IGRAPH_FINALLY_CLEAN(4);

  if (to->shorthand) { igraph_vs_destroy((igraph_vs_t*)to); }
  return 0;
}

void igraph_i_gasp_paths_destroy(igraph_vector_ptr_t *v) {
  long int i;
  for (i=0; i<igraph_vector_ptr_size(v); i++) {
    if (VECTOR(*v)[i] != 0) {
      igraph_vector_destroy(VECTOR(*v)[i]);
      Free(VECTOR(*v)[i]);
    }
  }
  igraph_vector_ptr_destroy(v);
}

/**
 * \function igraph_get_all_shortest_paths
 * \brief Finds all shortest paths (geodesics) from a vertex to all
 * other vertices 
 * 
 * \param graph The graph object.
 * \param res Pointer to an initialized pointer vector, the result
 *   will be stored here in igraph_vector_t objects. Each vector
 *   object contains the vertices along a shortest path from \p from
 *   to another vertex. The vectors are ordered according to their
 *   target vertex: first the shortest paths to vertex 0, then to
 *   vertex 1, etc. No data is included for unreachable vertices.
 * \param nrgeo Pointer to an initialized igraph_vector_t object or
 *   NULL. If not NULL the number of shortest paths from \p from are
 * is stored here for every vertex in the graph.
 * \param from The id of the vertex from/to which the geodesics are
 *        calculated. 
 * \param mode The type of shortest paths to be use for the
 *        calculation in directed graphs. Possible values: 
 *        \clist
 *        \cli IGRAPH_OUT 
 *          the lengths of the outgoing paths are calculated. 
 *        \cli IGRAPH_IN 
 *          the lengths of the incoming paths are calculated. 
 *        \cli IGRAPH_ALL 
 *          the directed graph is considered as an
 *          undirected one for the computation.
 *        \endclist
 * \return Error code:
 *        \clist
 *        \cli IGRAPH_ENOMEM 
 *           not enough memory for temporary data.
 *        \cli IGRAPH_EINVVID
 *           \p from is invalid vertex id.
 *        \cli IGRAPH_EINVMODE 
 *           invalid mode argument.
 *        \endclist
 *
 * Added in version 0.2.</para><para>
 *
 * Time complexity: O(|V|+|E|) for most graphs, O(|V|^2) in the worst
 * case. 
 */

int igraph_get_all_shortest_paths(const igraph_t *graph,
				  igraph_vector_ptr_t *res, 
				  igraph_vector_t *nrgeo,
				  integer_t from, const igraph_vs_t *to,
				  igraph_neimode_t mode) {
  
  long int no_of_nodes=igraph_vcount(graph);
  long int *geodist;
  igraph_vector_ptr_t paths;
  igraph_dqueue_t q;
  igraph_vector_t *vptr;
  igraph_vector_t neis;
  igraph_vector_t ptrlist;
  igraph_vector_t ptrhead;
  long int n, j, i;

  igraph_vs_t myto;
  const igraph_vector_t *mytov;

  if (from<0 || from>=no_of_nodes) {
    IGRAPH_ERROR("cannot get shortest paths", IGRAPH_EINVVID);
  }
  if (mode != IGRAPH_OUT && mode != IGRAPH_IN && 
      mode != IGRAPH_ALL) {
    IGRAPH_ERROR("Invalid mode argument", IGRAPH_EINVMODE);
  }

  IGRAPH_CHECK(igraph_vs_vectorview_it(graph, to, &myto));
  IGRAPH_FINALLY(igraph_vs_destroy, &myto);
  mytov=igraph_vs_vector_getvector(graph, &myto);

  IGRAPH_CHECK(igraph_vector_ptr_init(&paths, 0));
  IGRAPH_FINALLY(igraph_i_gasp_paths_destroy, &paths);
  IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&ptrlist, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&ptrhead, no_of_nodes);
  geodist=Calloc(no_of_nodes, long int);
  if (geodist==0) {
    IGRAPH_ERROR("Cannot calculate shortest paths", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, geodist);
  IGRAPH_CHECK(igraph_dqueue_init(&q, 100));
  IGRAPH_FINALLY(igraph_dqueue_destroy, &q);

  if (nrgeo) { 
    IGRAPH_CHECK(igraph_vector_resize(nrgeo, no_of_nodes));
    igraph_vector_null(nrgeo);
  }
  
  /* from -> from */
  vptr=Calloc(1, igraph_vector_t); /* TODO: dirty */
  IGRAPH_CHECK(igraph_vector_ptr_push_back(&paths, vptr));
  IGRAPH_CHECK(igraph_vector_init(vptr, 1));
  VECTOR(*vptr)[0]=from;
  geodist[(long int)from]=1;
  VECTOR(ptrhead)[(long int)from]=1;
  IGRAPH_CHECK(igraph_vector_push_back(&ptrlist, 0));
  if (nrgeo) { VECTOR(*nrgeo)[(long int)from]=1; }

  /* Init queue */
  IGRAPH_CHECK(igraph_dqueue_push(&q, from));
  IGRAPH_CHECK(igraph_dqueue_push(&q, 0.0));
  while (!igraph_dqueue_empty(&q)) {
    long int actnode=igraph_dqueue_pop(&q);
    long int actdist=igraph_dqueue_pop(&q);
    
    IGRAPH_CHECK(igraph_neighbors(graph, &neis, actnode, mode));
    n=igraph_vector_size(&neis);
    for (j=0; j<n; j++) {
      long int neighbor=VECTOR(neis)[j];
      long int fatherptr=VECTOR(ptrhead)[actnode];
      if (geodist[neighbor] != 0 && 
	  geodist[neighbor]-1 < actdist+1) { continue; }
      if (nrgeo) { VECTOR(*nrgeo)[neighbor] += VECTOR(*nrgeo)[actnode]; }
      if (geodist[neighbor] == 0) {
	IGRAPH_CHECK(igraph_dqueue_push(&q, neighbor));
	IGRAPH_CHECK(igraph_dqueue_push(&q, actdist+1));
      }
      geodist[neighbor]=actdist+2;

      /* copy all existing paths to the parent */
      while (fatherptr != 0) {
	long int tmp;
	vptr=Calloc(1, igraph_vector_t);
	IGRAPH_CHECK(igraph_vector_ptr_push_back(&paths, vptr));
	IGRAPH_CHECK(igraph_vector_copy(vptr, VECTOR(paths)[fatherptr-1]));
	IGRAPH_CHECK(igraph_vector_reserve(vptr, actdist+2));
	igraph_vector_push_back(vptr, neighbor);

	IGRAPH_CHECK(igraph_vector_push_back(&ptrlist, 
					     VECTOR(ptrhead)[neighbor]));
	VECTOR(ptrhead)[neighbor]=igraph_vector_size(&ptrlist);
	
	fatherptr=VECTOR(ptrlist)[fatherptr-1];
      }
    }
  }

  igraph_dqueue_destroy(&q);
  IGRAPH_FINALLY_CLEAN(1);

  /* Copy to the result */
  memset(geodist, 0, sizeof(long int)*no_of_nodes);
  for (i=0; i<igraph_vector_size(mytov); i++) {
    geodist[ (long int) VECTOR(*mytov)[i] ] = 1;
  }
  
  n=0;
  for (i=0; i<no_of_nodes; i++) {
    long int fatherptr=VECTOR(ptrhead)[i];
    if (geodist[i] > 0) {
      while (fatherptr != 0) {
	n++;
	fatherptr=VECTOR(ptrlist)[fatherptr-1];
      }
    }
  }

  IGRAPH_CHECK(igraph_vector_ptr_resize(res, n));
  j=0;
  for (i=0; i<no_of_nodes; i++) {
    long int fatherptr=VECTOR(ptrhead)[i];
    if (geodist[i] > 0) {
      while (fatherptr != 0) {
	VECTOR(*res)[j++]=VECTOR(paths)[fatherptr-1];
	fatherptr=VECTOR(ptrlist)[fatherptr-1];
      }
    } else {
      while (fatherptr != 0) {
	igraph_vector_destroy(VECTOR(paths)[fatherptr-1]);
	Free(VECTOR(paths)[fatherptr-1]);
	fatherptr=VECTOR(ptrlist)[fatherptr-1];
      }
    }
  }

  Free(geodist);
  igraph_vector_destroy(&ptrlist);
  igraph_vector_destroy(&ptrhead);
  igraph_vector_destroy(&neis);
  igraph_vector_ptr_destroy(&paths);
  igraph_vs_destroy(&myto);
  IGRAPH_FINALLY_CLEAN(6);

  if (to->shorthand) { igraph_vs_destroy((igraph_vs_t*)to); }
  return 0;
}
				  

/** 
 * \ingroup structural
 * \function igraph_subcomponent
 * \brief The vertices in the same component as a given vertex.
 *
 * \param graph The graph object.
 * \param res The result, vector with the ids of the vertices in the
 *        same component. 
 * \param vertex The id of the vertex of which the component is
 *        searched. 
 * \param mode Type of the component for directed graphs, possible
 *        values:
 *        \clist
 *        \cli IGRAPH_OUT 
 *          the set of vertices reachable \em from the
 *          \p vertex, 
 *        \cli IGRAPH_IN
 *          the set of vertices from which the
 *          \p vertex is reachable.
 *        \cli IGRAPH_ALL 
 *          the graph is considered as an
 *          undirected graph. Note that this is \em not the same
 *          as the union of the previous two.
 *        \endclist
 * \return Error code:
 *        \clist
 *        \cli IGRAPH_ENOMEM 
 *          not enough memory for temporary data.
 *        \cli IGRAPH_EINVVID 
 *           \p vertex is an invalid vertex id
 *        \cli IGRAPH_EINVMODE
 *           invalid mode argument passed.
 *        \endclist
 * 
 * Time complexity: O(|V|+|E|),
 * |V| and
 * |E| are the number of vertices and
 * edges in the graph. 
 * 
 * \sa \ref igraph_subgraph() if you want a graph object consisting only
 * a given set of vertices and the edges between them.
 */

int igraph_subcomponent(const igraph_t *graph, igraph_vector_t *res, real_t vertex, 
			igraph_neimode_t mode) {

  long int no_of_nodes=igraph_vcount(graph);
  igraph_dqueue_t q=IGRAPH_DQUEUE_NULL;
  char *already_added;
  long int i;
  igraph_vector_t tmp=IGRAPH_VECTOR_NULL;

  if (vertex<0 || vertex>=no_of_nodes) {
    IGRAPH_ERROR("subcomponent failed", IGRAPH_EINVVID);
  }
  if (mode != IGRAPH_OUT && mode != IGRAPH_IN && 
      mode != IGRAPH_ALL) {
    IGRAPH_ERROR("invalid mode argument", IGRAPH_EINVMODE);
  }

  already_added=Calloc(no_of_nodes, char);
  if (already_added==0) {
    IGRAPH_ERROR("subcomponent failed",IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(free, already_added); /* TODO: hack */

  IGRAPH_VECTOR_INIT_FINALLY(&tmp, 0);
  IGRAPH_DQUEUE_INIT_FINALLY(&q, 100);
  
  IGRAPH_CHECK(igraph_dqueue_push(&q, vertex));
  IGRAPH_CHECK(igraph_vector_push_back(res, vertex));
  already_added[(long int)vertex]=1;
  
  while (!igraph_dqueue_empty(&q)) {
    long int actnode=igraph_dqueue_pop(&q);

    IGRAPH_CHECK(igraph_neighbors(graph, &tmp, actnode, mode));
    for (i=0; i<igraph_vector_size(&tmp); i++) {
      long int neighbor=VECTOR(tmp)[i];
      
      if (already_added[neighbor]) { continue; }
      already_added[neighbor]=1;
      IGRAPH_CHECK(igraph_vector_push_back(res, neighbor));
      IGRAPH_CHECK(igraph_dqueue_push(&q, neighbor));
    }
  }

  igraph_dqueue_destroy(&q);
  igraph_vector_destroy(&tmp);
  Free(already_added);
  IGRAPH_FINALLY_CLEAN(3);
   
  return 0;
}

/**
 * \ingroup structural
 * \function igraph_betweenness
 * \brief Betweenness centrality of some vertices.
 * 
 * The betweenness centrality of a vertex is the number of geodesics
 * going through it. If there are more than one geodesic between two
 * vertices, the value of these geodesics are weighted by one over the 
 * number of geodesics.
 * \param graph The graph object.
 * \param res The result of the computation, a vector containing the
 *        betweenness scores for the specified vertices.
 * \param vids The vertices of which the betweenness centrality scores
 *        will be calculated.
 * \param directed Logical, if true directed paths will be considered
 *        for directed graphs. It is ignored for undirected graphs.
 * \return Error code:
 *        \c IGRAPH_ENOMEM, not enough memory for
 *        temporary data. 
 *        \c IGRAPH_EINVVID, invalid vertex id passed in
 *        \p vids. 
 *
 * Time complexity: O(|V||E|),
 * |V| and 
 * |E| are the number of vertices and
 * edges in the graph. 
 * Note that the time complexity is independent of the number of
 * vertices for which the score is calculated.
 *
 * \sa Other centrality types: \ref igraph_degree(), \ref igraph_closeness().
 *     See \ref igraph_edge_betweenness() for calculating the betweenness score
 *     of the edges in a graph.
 */

int igraph_betweenness (const igraph_t *graph, igraph_vector_t *res, 
			const igraph_vs_t *vids, 
			bool_t directed) {

  long int no_of_nodes=igraph_vcount(graph);
  igraph_dqueue_t q=IGRAPH_DQUEUE_NULL;
  long int *distance;
  long int *nrgeo;
  double *tmpscore;
  igraph_stack_t stack=IGRAPH_STACK_NULL;
  long int source;
  long int j;
  igraph_vector_t tmp=IGRAPH_VECTOR_NULL;
  integer_t modein, modeout;
  igraph_vs_t myvids;
  const igraph_vector_t *myvidsv;

  IGRAPH_CHECK(igraph_vs_vectorview_it(graph, vids, &myvids));
  IGRAPH_FINALLY(igraph_vs_destroy, &myvids);
  myvidsv=igraph_vs_vector_getvector(graph, &myvids);

  if (!igraph_vector_isininterval(myvidsv, 0, no_of_nodes-1)) {
    IGRAPH_ERROR("betweenness failed", IGRAPH_EINVVID);
  }

  if (directed) 
    { modeout=IGRAPH_OUT; modein=IGRAPH_IN; } 
  else 
    { modeout=modein=IGRAPH_ALL; }

  distance=Calloc(no_of_nodes, long int);
  if (distance==0) {
    IGRAPH_ERROR("betweenness failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(free, distance); /* TODO: hack */
  nrgeo=Calloc(no_of_nodes, long int);
  if (nrgeo==0) {
    IGRAPH_ERROR("betweenness failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(free, nrgeo);	/* TODO: hack */
  tmpscore=Calloc(no_of_nodes, double);
  if (tmpscore==0) {
    IGRAPH_ERROR("betweenness failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(free, tmpscore); /* TODO: hack */

  IGRAPH_VECTOR_INIT_FINALLY(&tmp, 0);
  IGRAPH_DQUEUE_INIT_FINALLY(&q, 100);
  igraph_stack_init(&stack, no_of_nodes);
  IGRAPH_FINALLY(igraph_stack_destroy, &stack);
    
  IGRAPH_CHECK(igraph_vector_resize(res, igraph_vector_size(myvidsv)));
  igraph_vector_null(res);

  /* here we go */
  
  for (source=0; source<no_of_nodes; source++) {
    
    memset(distance, 0, no_of_nodes*sizeof(long int));
    memset(nrgeo, 0, no_of_nodes*sizeof(long int));
    memset(tmpscore, 0, no_of_nodes*sizeof(double));
    igraph_stack_clear(&stack); /* it should be empty anyway... */
    
    IGRAPH_CHECK(igraph_dqueue_push(&q, source));
    nrgeo[source]=1;
    distance[source]=0;
    
    while (!igraph_dqueue_empty(&q)) {
      long int actnode=igraph_dqueue_pop(&q);

      IGRAPH_CHECK(igraph_neighbors(graph, &tmp, actnode, modeout));
      for (j=0; j<igraph_vector_size(&tmp); j++) {
	long int neighbor=VECTOR(tmp)[j];
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
      if (distance[actnode]<=1) { continue; } /* skip source node */
      
      /* set the temporary score of the friends */
      IGRAPH_CHECK(igraph_neighbors(graph, &tmp, actnode, modein));
      for (j=0; j<igraph_vector_size(&tmp); j++) {
	long int neighbor=VECTOR(tmp)[j];
	if (distance[neighbor]==distance[actnode]-1 &&
	    nrgeo[neighbor] != 0) {
	  tmpscore[neighbor] += 
	    (tmpscore[actnode]+1)*nrgeo[neighbor]/nrgeo[actnode];
	}
      }
    }
    
    /* Ok, we've the scores for this source */
    for (j=0; j<igraph_vector_size(myvidsv); j++) {
      long int node=VECTOR(*myvidsv)[j];
      VECTOR(*res)[node] += tmpscore[node];
      tmpscore[node] = 0.0; /* in case a node is in vids multiple times */
    }

  } /* for source < no_of_nodes */

  /* divide by 2 for undirected graph */
  if (!directed || !igraph_is_directed(graph)) {
    for (j=0; j<igraph_vector_size(res); j++) {
      VECTOR(*res)[j] /= 2.0;
    }
  }
  
  /* clean  */
  Free(distance);
  Free(nrgeo);
  Free(tmpscore);
  
  igraph_dqueue_destroy(&q);
  igraph_stack_destroy(&stack);
  igraph_vector_destroy(&tmp);
  igraph_vs_destroy(&myvids);
  IGRAPH_FINALLY_CLEAN(7);

  if (vids->shorthand) { igraph_vs_destroy((igraph_vs_t*)vids); }
  return 0;
}

/**
 * \ingroup structural
 * \function igraph_edge_betweenness
 * \brief Betweenness centrality of the edges.
 * 
 * The betweenness centrality of an edge is the number of geodesics
 * going through it. If there are more than one geodesics between two
 * vertices, the value of these geodesics are weighted by one over the 
 * number of geodesics.
 * \param graph The graph object.
 * \param result The result of the computation, vector containing the
 *        betweenness scores for the edges.
 * \param directed Logical, if true directed paths will be considered
 *        for directed graphs. It is ignored for undirected graphs.
 * \return Error code:
 *        \c IGRAPH_ENOMEM, not enough memory for
 *        temporary data. 
 *
 * Time complexity: O(|V||E|),
 * |V| and
 * |E| are the number of vertices and
 * edges in the graph. 
 *
 * \sa Other centrality types: \ref igraph_degree(), \ref igraph_closeness().
 *     See \ref igraph_edge_betweenness() for calculating the betweenness score
 *     of the edges in a graph.
 */

int igraph_edge_betweenness (const igraph_t *graph, igraph_vector_t *result, 
			     bool_t directed) {
  
  long int no_of_nodes=igraph_vcount(graph);
  long int no_of_edges=igraph_ecount(graph);
  igraph_dqueue_t q=IGRAPH_DQUEUE_NULL;
  long int *distance;
  long int *nrgeo;
  double *tmpscore;
  igraph_stack_t stack=IGRAPH_STACK_NULL;
  long int source;
  long int j;

  igraph_es_t it;
  integer_t modein, modeout;

  if (directed) {
    modeout=IGRAPH_OUT;
    modein=IGRAPH_IN;
  } else {
    modeout=modein=IGRAPH_ALL;
  }
  
  distance=Calloc(no_of_nodes, long int);
  if (distance==0) {
    IGRAPH_ERROR("edge betweenness failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(free, distance); /* TODO: hack */
  nrgeo=Calloc(no_of_nodes, long int);
  if (nrgeo==0) {
    IGRAPH_ERROR("edge betweenness failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(free, nrgeo);	/* TODO: hack */
  tmpscore=Calloc(no_of_nodes, double);
  if (tmpscore==0) {
    IGRAPH_ERROR("edge betweenness failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(free, tmpscore); /* TODO: hack */

  IGRAPH_DQUEUE_INIT_FINALLY(&q, 100);
  IGRAPH_CHECK(igraph_stack_init(&stack, no_of_nodes));
  IGRAPH_FINALLY(igraph_stack_destroy, &stack);

  IGRAPH_CHECK(igraph_es_adj(graph, &it, 0, modeout));
  IGRAPH_FINALLY(igraph_es_destroy, &it);
  IGRAPH_CHECK(igraph_vector_resize(result, no_of_edges));

  igraph_vector_null(result);

  /* here we go */
  
  for (source=0; source<no_of_nodes; source++) {

    memset(distance, 0, no_of_nodes*sizeof(long int));
    memset(nrgeo, 0, no_of_nodes*sizeof(long int));
    memset(tmpscore, 0, no_of_nodes*sizeof(double));
    igraph_stack_clear(&stack); /* it should be empty anyway... */
    
    IGRAPH_CHECK(igraph_dqueue_push(&q, source));
      
    nrgeo[source]=1;
    distance[source]=0;
    
    while (!igraph_dqueue_empty(&q)) {
      long int actnode=igraph_dqueue_pop(&q);
    
      igraph_es_adj_set(graph, &it, actnode, modeout);
      while ( !igraph_es_end(graph, &it)) {
	long int neighbor=igraph_es_adj_vertex(graph, &it);
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
	igraph_es_next(graph, &it);
      }
    } /* while !igraph_dqueue_empty */
    
    /* Ok, we've the distance of each node and also the number of
       shortest paths to them. Now we do an inverse search, starting
       with the farthest nodes. */
    while (!igraph_stack_empty(&stack)) {
      long int actnode=igraph_stack_pop(&stack);
      if (distance[actnode]<1) { continue; } /* skip source node */
      
      /* set the temporary score of the friends */
      igraph_es_adj_set(graph, &it, actnode, modein);
      while ( !igraph_es_end(graph,  &it)) {
	long int neighbor=igraph_es_adj_vertex(graph, &it);
	long int edgeno=igraph_es_get(graph, &it);
	if (distance[neighbor]==distance[actnode]-1 &&
	    nrgeo[neighbor] != 0) {
	  tmpscore[neighbor] += 
	    (tmpscore[actnode]+1)*nrgeo[neighbor]/nrgeo[actnode];
	  VECTOR(*result)[edgeno] += 
	    (tmpscore[actnode]+1)*nrgeo[neighbor]/nrgeo[actnode];
	}
	igraph_es_next(graph, &it);	
      }
    }
    /* Ok, we've the scores for this source */
  } /* for source <= no_of_nodes */
  
  /* clean and return */
  Free(distance);
  Free(nrgeo);
  Free(tmpscore);
  igraph_dqueue_destroy(&q);
  igraph_stack_destroy(&stack);
  igraph_es_destroy(&it);
  IGRAPH_FINALLY_CLEAN(6);

  /* divide by 2 for undirected graph */
  if (!directed || !igraph_is_directed(graph)) {
    for (j=0; j<igraph_vector_size(result); j++) {
      VECTOR(*result)[j] /= 2.0;
    }
  }
  
  return 0;
}

/**
 * \ingroup structural
 * \function igraph_pagerank
 * \brief Calculates the Google PageRank for the specified vertices.
 * 
 * Please note that the PageRank of a given vertex depends on the PageRank
 * of all other vertices, so even if you want to calculate the PageRank for
 * only some of the vertices, all of them must be calculated. Requesting
 * the PageRank for only some of the vertices does not result in any
 * performance increase at all.
 * </para>
 * <para>
 * Since the calculation is an iterative
 * process, the algorithm is stopped after a given count of iterations
 * or if the PageRank value differences between iterations are less than
 * a predefined value.
 * </para>
 * 
 * <para>
 * For the explanation of the PageRank algorithm, see the following
 * webpage:
 * http://www-db.stanford.edu/~backrub/google.html, or the
 * following reference:
 * </para>
 * 
 * <para>
 * Sergey Brin and Larry Page: The Anatomy of a Large-Scale Hypertextual
 * Web Search Engine. Proceedings of the 7th World-Wide Web Conference,
 * Brisbane, Australia, April 1998.
 * </para>
 * <para>
 * \param graph The graph object.
 * \param res The result vector containing the PageRank values for the
 * given nodes.
 * \param vids Vector with the vertex ids
 * \param directed Logical, if true directed paths will be considered
 *        for directed graphs. It is ignored for undirected graphs.
 * \param niter The maximum number of iterations to perform
 * \param eps The algorithm will consider the calculation as complete
 *        if the difference of PageRank values between iterations change
 *        less than this value for every node
 * \param damping The damping factor ("d" in the original paper)
 * \return Error code:
 *         \c IGRAPH_ENOMEM, not enough memory for
 *         temporary data. 
 *         \c IGRAPH_EINVVID, invalid vertex id in
 *         \p vids. 
 */

int igraph_pagerank(const igraph_t *graph, igraph_vector_t *res, 
		    const igraph_vs_t *vids, bool_t directed, integer_t niter, 
		    real_t eps, real_t damping) {
  long int no_of_nodes=igraph_vcount(graph);
  long int i, j, n, nodes_to_calc;
  real_t *prvec, *prvec_new, *prvec_aux, *prvec_scaled;
  igraph_vector_t *neis, outdegree;
  integer_t dirmode;
  igraph_i_adjlist_t allneis;
  real_t maxdiff=eps;
  igraph_vs_t myvids;
  const igraph_vector_t * myvidsv;

  if (niter<=0) IGRAPH_ERROR("Invalid iteration count", IGRAPH_EINVAL);
  if (eps<=0) IGRAPH_ERROR("Invalid epsilon value", IGRAPH_EINVAL);
  if (damping<=0 || damping>=1) IGRAPH_ERROR("Invalid damping factor", IGRAPH_EINVAL);

  IGRAPH_CHECK(igraph_vs_vectorview_it(graph, vids, &myvids));
  IGRAPH_FINALLY(igraph_vs_destroy, &myvids);
  myvidsv=igraph_vs_vector_getvector(graph, &myvids);
  nodes_to_calc=igraph_vector_size(myvidsv);    

  IGRAPH_CHECK(igraph_vector_resize(res, nodes_to_calc));
  igraph_vector_null(res);
  
  IGRAPH_VECTOR_INIT_FINALLY(&outdegree, no_of_nodes);
    
  prvec=Calloc(no_of_nodes, real_t);
  if (prvec==0) {
    IGRAPH_ERROR("pagerank failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, prvec);
  
  prvec_new=Calloc(no_of_nodes, real_t);
  if (prvec_new==0) {
    IGRAPH_ERROR("pagerank failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, prvec_new);
  
  prvec_scaled=Calloc(no_of_nodes, real_t);
  if (prvec_scaled==0) {
    IGRAPH_ERROR("pagerank failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, prvec_scaled);
  
  if (directed) { dirmode=IGRAPH_IN; } else { dirmode=IGRAPH_ALL; }  
  igraph_i_adjlist_init(graph, &allneis, dirmode);
  IGRAPH_FINALLY(igraph_i_adjlist_destroy, &allneis);

  /* Calculate outdegrees for every node */
  igraph_degree(graph, &outdegree, IGRAPH_VS_ALL(graph),
		directed?IGRAPH_OUT:IGRAPH_ALL, 0);
  /* Initialize PageRank values */
  for (i=0; i<no_of_nodes; i++) {
    prvec[i]=1-damping;
    /* The next line is necessary to avoid division by zero in the
     * calculation of prvec_scaled. This won't cause any problem,
     * since if a node doesn't have any outgoing links, its
     * prvec_scaled value won't be used anywhere */
    if (VECTOR(outdegree)[i]==0) VECTOR(outdegree)[i]=1;
  }
  
  /* We will always calculate the new PageRank values into prvec_new
   * based on the existing values from prvec. To avoid unnecessary
   * copying from prvec_new to prvec at the end of every iteration,
   * the pointers are swapped after every iteration */
  while (niter>0 && maxdiff >= eps) {
    niter--;
    maxdiff=0;

    /* Calculate the quotient of the actual PageRank value and the
     * outdegree for every node */
    for (i=0; i<no_of_nodes; i++) {
      prvec_scaled[i]=prvec[i]/VECTOR(outdegree)[i];
    }
    
    /* Calculate new PageRank values based on the old ones */
    for (i=0; i<no_of_nodes; i++) {
      prvec_new[i]=0;
      neis=igraph_i_adjlist_get(&allneis, i);
      n=igraph_vector_size(neis);
      for (j=0; j<n; j++) {
	long int neighbor=VECTOR(*neis)[j];
	prvec_new[i]+=prvec_scaled[neighbor];
      }
      prvec_new[i]*=damping;
      prvec_new[i]+=(1-damping);

      if (prvec_new[i]-prvec[i]>maxdiff) 
	maxdiff=prvec_new[i]-prvec[i];
      else if (prvec[i]-prvec_new[i]>maxdiff)
	maxdiff=prvec[i]-prvec_new[i];
    }
    
    /* Swap the vectors */
    prvec_aux=prvec_new;
    prvec_new=prvec;
    prvec=prvec_aux;
  }
  
  /* Copy results from prvec to res */
  for (i=0; i<nodes_to_calc; i++) {
    long int vid=VECTOR(*myvidsv)[i];
    VECTOR(*res)[i]=prvec[vid];
  }
  
  igraph_i_adjlist_destroy(&allneis);
  igraph_vs_destroy(&myvids);
  igraph_vector_destroy(&outdegree);
  Free(prvec);
  Free(prvec_new);  
  Free(prvec_scaled);
  
  IGRAPH_FINALLY_CLEAN(6);
  
  if (vids->shorthand) { igraph_vs_destroy((igraph_vs_t*)vids); }
  return 0;
}

/**
 * \ingroup structural
 * \function igraph_rewire
 * \brief Randomly rewires a graph while preserves the degree distribution.
 * 
 * This function generates a new graph based on the original one by randomly
 * rewiring edges while preserving the original graph's degree distribution.
 * Please note that the rewiring is done "in place", so no new graph will
 * be generated. If you would like to keep the original graph intact, use
 * \ref igraph_copy() before.
 * 
 * \param graph The graph object to be rewired.
 * \param n Number of rewiring trials to perform.
 * \param mode The rewiring algorithm to be used. It can be one of the following:
 *         \c IGRAPH_REWIRING_SIMPLE: simple rewiring algorithm which
 *         chooses two arbitrary edges in each step (namely (a,b) and (c,d))
 *         and substitutes them with (a,d) and (c,b) if they don't exist.
 *         Time complexity: TODO.
 * \return Error code:
 *         \clist
 *           \cli IGRAPH_EINVMODE
 *                Invalid rewiring mode.
 *           \cli IGRAPH_EINVAL
 *                Graph unsuitable for rewiring (e.g. it has
 *                less than 4 nodes in case of \c IGRAPH_REWIRING_SIMPLE)
 *           \cli IGRAPH_ENOMEM
 *                Not enough memory for temporary data.
 *         \endclist
 */

int igraph_rewire(igraph_t *graph, integer_t n, igraph_rewiring_t mode) {
  long int no_of_nodes=igraph_vcount(graph);
  long int i, a, b, c, d;
  igraph_i_adjlist_t allneis;
  igraph_vector_t *neis[2], edgevec;
  
  if (mode == IGRAPH_REWIRING_SIMPLE && no_of_nodes<4)
    IGRAPH_ERROR("graph unsuitable for rewiring", IGRAPH_EINVAL);
  
  RNG_BEGIN();
  
  igraph_i_adjlist_init(graph, &allneis, IGRAPH_OUT);
  IGRAPH_FINALLY(igraph_i_adjlist_destroy, &allneis);
  igraph_vector_init(&edgevec, 4);
  IGRAPH_FINALLY(igraph_vector_destroy, &edgevec);
  
  while (n>0) {
    switch (mode) {
    case IGRAPH_REWIRING_SIMPLE:
      a=RNG_INTEGER(0, no_of_nodes-1);
      do { c=RNG_INTEGER(0, no_of_nodes-1); } while (c==a);
      
      /* We don't want the algorithm to get stuck in an infinite loop when
       * it can't choose two edges satisfying the conditions. Instead of
       * this, we choose two arbitrary edges and if they have endpoints
       * in common, we just decrease the number of trials left and continue
       * (so unsuccessful rewirings still count as a trial)
       */
      neis[0]=igraph_i_adjlist_get(&allneis, a);
      i=igraph_vector_size(neis[0]);
      if (i==0) b=c;
      else b=VECTOR(*neis[0])[RNG_INTEGER(0, i-1)];
      
      neis[1]=igraph_i_adjlist_get(&allneis, c);
      i=igraph_vector_size(neis[1]);
      if (i==0) d=a;
      else d=VECTOR(*neis[1])[RNG_INTEGER(0, i-1)];
      
      /* Okay, we have two edges. Can they be rewired?
       * neis[0] mustn't contain d and neis[1] mustn't contain b
       */
      if (!igraph_vector_search(neis[0], 0, d, NULL) &&
	  !igraph_vector_search(neis[1], 0, b, NULL) &&
	  b!=c && a!=d && a!=b && c!=d) {
	/* TODO: maybe it would be more efficient to implement an
	 * igraph_replace_edges function?
	 */
	VECTOR(edgevec)[0]=a; VECTOR(edgevec)[1]=b;
	VECTOR(edgevec)[2]=c; VECTOR(edgevec)[3]=d;
	/*printf("Deleting: %d -> %d, %d -> %d\n", a, b, c, d);*/
	igraph_delete_edges(graph, &edgevec);
	VECTOR(edgevec)[0]=a; VECTOR(edgevec)[1]=d;
	VECTOR(edgevec)[2]=c; VECTOR(edgevec)[3]=b;
	/*printf("Adding: %d -> %d, %d -> %d\n", a, d, c, b);*/
	igraph_add_edges(graph, &edgevec);
	/* We have to adjust the adjacency list view as well.
	   It's a luck that we have the pointers in neis[0] and neis[1] */
	for (i=igraph_vector_size(neis[0])-1; i>=0; i--)
	  if (VECTOR(*neis[0])[i]==b) { VECTOR(*neis[0])[i]=d; break; }
	for (i=igraph_vector_size(neis[1])-1; i>=0; i--)
	  if (VECTOR(*neis[1])[i]==d) { VECTOR(*neis[1])[i]=b; break; }
      }
      break;
    default:
      RNG_END();
      IGRAPH_ERROR("unknown rewiring mode", IGRAPH_EINVMODE);
    }
    n--;
  }
  
  igraph_i_adjlist_destroy(&allneis);
  igraph_vector_destroy(&edgevec);
  IGRAPH_FINALLY_CLEAN(2);
  
  RNG_END();
  
  return 0;
}

/**
 * \ingroup structural
 * \function igraph_subgraph
 * \brief Creates a subgraph with the specified vertices.
 * 
 * This function collects the specified vertices and all edges between
 * them to a new graph.
 * As the vertex ids in a graph always start with one, this function
 * very likely needs to reassign ids to the vertices.
 * \param graph The graph object.
 * \param res The subgraph, another graph object will be stored here,
 *        do \em not initialize this object before calling this
 *        function, and call \ref igraph_destroy() on it if you don't need
 *        it any more.
 * \param vids Vector with the vertex ids to put in the subgraph.
 * \return Error code:
 *         \c IGRAPH_ENOMEM, not enough memory for
 *         temporary data. 
 *         \c IGRAPH_EINVVID, invalid vertex id in
 *         \p vids. 
 * 
 * Time complexity: O(|V|+|E|),
 * |V| and
 * |E| are the number of vertices and
 * edges in the original graph.
 *
 * \sa \ref igraph_delete_vertices() to delete the specified set of
 * vertices from a graph, the opposite of this function.
 */

int igraph_subgraph(const igraph_t *graph, igraph_t *res, 
		    const igraph_vs_t* vids) {
  
  long int no_of_nodes=igraph_vcount(graph);
  igraph_vector_t delete=IGRAPH_VECTOR_NULL;
  char *remain;
  long int i;
  igraph_vs_t myvids;
  const igraph_vector_t *myvidsv;
  
  IGRAPH_CHECK(igraph_vs_vectorview_it(graph, vids, &myvids));
  IGRAPH_FINALLY(igraph_vs_destroy, &myvids);
  myvidsv=igraph_vs_vector_getvector(graph, &myvids);

  if (!igraph_vector_isininterval(myvidsv, 0, no_of_nodes-1)) {
    IGRAPH_ERROR("subgraph failed", IGRAPH_EINVVID);
  }

  IGRAPH_VECTOR_INIT_FINALLY(&delete, 0);
  remain=Calloc(no_of_nodes, char);
  if (remain==0) {
    IGRAPH_ERROR("subgraph failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(free, remain);	/* TODO: hack */
  IGRAPH_CHECK(igraph_vector_reserve(&delete, no_of_nodes-igraph_vector_size(myvidsv)));
  
  for (i=0; i<igraph_vector_size(myvidsv); i++) {
    remain[ (long int) VECTOR(*myvidsv)[i] ] = 1;
  }

  for (i=0; i<no_of_nodes; i++) {
    if (remain[i] == 0) {
      IGRAPH_CHECK(igraph_vector_push_back(&delete, i));
    }
  }

  Free(remain);
  IGRAPH_FINALLY_CLEAN(1);
  IGRAPH_CHECK(igraph_copy(res, graph));
  IGRAPH_FINALLY(igraph_destroy, res);
  IGRAPH_CHECK(igraph_delete_vertices(res, IGRAPH_VS_VECTOR(res,&delete)));
  
  igraph_vector_destroy(&delete);
  igraph_vs_destroy(&myvids);
  IGRAPH_FINALLY_CLEAN(3);
  if (vids->shorthand) { igraph_vs_destroy((igraph_vs_t*)vids); }
  return 0;
}

/**
 * \ingroup structural
 * \function igraph_simplify
 * \brief Removes loop and/or multiple edges from the graph.
 * 
 * \param graph The graph object.
 * \param multiple Logical, if true, multiple edges will be removed. 
 * \param loops Logical, if true, loops (self edges) will be removed.
 * \return Error code:
 *    \c IGRAPH_ENOMEM if we are out of memory.
 *
 * Time complexity: O(|V|+|E|) for
 * removing the loops,
 * O(|V|d*log(d)+|E|) for removing
 * the multiple edges. |V| and
 * |E| are the number of vertices and
 * edges in the graph, d is the
 * highest out-degree in the graph.
 */

int igraph_simplify(igraph_t *graph, bool_t multiple, bool_t loops) {

  igraph_vector_t edges=IGRAPH_VECTOR_NULL;
  igraph_vector_t neis=IGRAPH_VECTOR_NULL;
  long int no_of_nodes=igraph_vcount(graph);
  long int i, j;

  IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);

  for (i=0; i<no_of_nodes; i++) {
    IGRAPH_CHECK(igraph_neighbors(graph, &neis, i, IGRAPH_OUT));

    if (loops) {
      for (j=0; j<igraph_vector_size(&neis); j++) {
	if (VECTOR(neis)[j]==i) {
	  IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
	  IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
	}
      }
    } /* if loops */
    
    if (multiple) {
      igraph_vector_sort(&neis);
      for (j=1; j<igraph_vector_size(&neis); j++) {
	if (VECTOR(neis)[j]==VECTOR(neis)[j-1]) {
	  IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
	  IGRAPH_CHECK(igraph_vector_push_back(&edges, VECTOR(neis)[j]));
	}
      }
    }
  }

  igraph_vector_destroy(&neis);
  IGRAPH_FINALLY_CLEAN(1);
  IGRAPH_CHECK(igraph_delete_edges(graph, &edges));
  igraph_vector_destroy(&edges);
  IGRAPH_FINALLY_CLEAN(1);

  return 0;
}

int igraph_transitivity_undirected(const igraph_t *graph, igraph_vector_t *res) {

  long int no_of_nodes=igraph_vcount(graph);
  real_t triples=0, triangles=0;
  long int node;
  long int *neis;
  long int deg;

  igraph_vs_t nit, nit2;
  
  neis=Calloc(no_of_nodes, long int);
  if (neis==0) {
    IGRAPH_ERROR("undirected transitivity failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(free, neis);	/* TODO: hack */
  IGRAPH_CHECK(igraph_vector_resize(res, 1));

  if (no_of_nodes != 0) {    
    IGRAPH_CHECK(igraph_vs_adj(graph, &nit, 0, IGRAPH_ALL));
    IGRAPH_CHECK(igraph_vs_adj(graph, &nit2, 0, IGRAPH_ALL));
  }
  for (node=0; node<no_of_nodes; node++) {
    igraph_vs_adj_set(graph, &nit, node, IGRAPH_ALL);
    deg=0;
    /* Mark the neighbors of 'node' */
    while (!igraph_vs_end(graph, &nit)) {
      neis[ (long int)igraph_vs_get(graph, &nit) ] = node+1;
      igraph_vs_next(graph, &nit);
      deg++;
    }
    triples += deg*(deg-1);
    
    /* Count the triangles and triples */
    igraph_vs_adj_set(graph, &nit, node, IGRAPH_ALL);
    while (!igraph_vs_end(graph, &nit)) {
      long int v=igraph_vs_get(graph, &nit);
      igraph_vs_adj_set(graph, &nit2, v, IGRAPH_ALL);
      while (!igraph_vs_end(graph, &nit2)) {
	long int v2=igraph_vs_get(graph, &nit2);
	if (neis[v2] == node+1) {
	  triangles += 1.0;
	}
	igraph_vs_next(graph, &nit2);
      }
      igraph_vs_next(graph, &nit);
    }
  }

  Free(neis);
  IGRAPH_FINALLY_CLEAN(1);

  VECTOR(*res)[0] = triangles/triples;

  return 0;
}

/**
 * \ingroup structural
 * \function igraph_transitivity
 * \brief Calculates the transitivity (clustering coefficient) of a graph.
 * 
 * The transitivity measures the probability that two neighbors of a
 * vertex are connected. See the \p type parameter for
 * different definitions (more to be expected soon).
 * \param graph The graph object.
 * \param res Pointer to an initialized vector, this will be resized
 *        to contain the result.
 * \param type It gives the type of the transitivity to
 *        calculate. Possible values:
 *        \clist
 *        \cli IGRAPH_TRANSITIVITY_UNDIRECTED
 *          the most common
 *          definition, the ratio of the triangles and connected
 *          triples in the graph, the result is a single real number
 *          or NaN (0/0) if there are no connected triples in the
 *          graph.  Directed graphs are considered as
 *          undirected ones. 
 *        \endclist
 * \return Error code:
 *         \c IGRAPH_EINVAL: unknown transitivity type.
 *         \c IGRAPH_ENOMEM: not enough memory for
 *         temporary data. 
 * 
 * Time complexity: O(|V|*d^2) for
 * IGRAPH_TRANSITIVITY_UNDIRECTED.
 * |V| is the number of vertices in
 * the graph, d is the highest node
 * degree. 
 */

int igraph_transitivity(const igraph_t *graph, igraph_vector_t *res, 
			igraph_transitivity_type_t type) {
  
  int retval;
  
  if (type == IGRAPH_TRANSITIVITY_UNDIRECTED) {
    retval=igraph_transitivity_undirected(graph, res);
  } else {
    IGRAPH_ERROR("unknown transitivity type", IGRAPH_EINVAL);
  }
  
  return retval;
}
