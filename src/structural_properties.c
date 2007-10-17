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
#include <limits.h>

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
 * \param pres Pointer to an integer, if not \c NULL then it will contain 
 *        the diameter (the actual distance).
 * \param pfrom Pointer to an integer, if not \c NULL it will be set to the 
 *        source vertex of the diameter path.
 * \param pto Pointer to an integer, if not \c NULL it will be set to the 
 *        target vertex of the diameter path.
 * \param path Pointer to an initialized vector. If not \c NULL the actual 
 *        longest geodesic path will be stored here. The vector will be 
 *        resized as needed.
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

int igraph_diameter(const igraph_t *graph, igraph_integer_t *pres, 
		    igraph_integer_t *pfrom, igraph_integer_t *pto, 
		    igraph_vector_t *path,
		    igraph_bool_t directed, igraph_bool_t unconn) {

  long int no_of_nodes=igraph_vcount(graph);
  long int i, j, n;
  long int *already_added;
  long int nodes_reached;
  long int from=0, to=0;
  long int res=0;

  igraph_dqueue_t q=IGRAPH_DQUEUE_NULL;
  igraph_vector_t *neis;
  igraph_integer_t dirmode;
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

    IGRAPH_PROGRESS("Diameter: ", 100.0*i/no_of_nodes, NULL);

    IGRAPH_ALLOW_INTERRUPTION();
    
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

  IGRAPH_PROGRESS("Diameter: ", 100.0, NULL);
  
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
					     igraph_vss_1(to), dirmode));
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

int igraph_average_path_length(const igraph_t *graph, igraph_real_t *res,
			       igraph_bool_t directed, igraph_bool_t unconn) {
  long int no_of_nodes=igraph_vcount(graph);
  long int i, j, n;
  long int *already_added;
  long int nodes_reached=0;
  igraph_real_t normfact=0.0;

  igraph_dqueue_t q=IGRAPH_DQUEUE_NULL;
  igraph_vector_t *neis;
  igraph_integer_t dirmode;
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

    IGRAPH_ALLOW_INTERRUPTION();
    
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

int igraph_path_length_hist(const igraph_t *graph, igraph_vector_t *res,
			    igraph_real_t *unconnected, igraph_bool_t directed) {

  long int no_of_nodes=igraph_vcount(graph);
  long int i,j,n;
  igraph_vector_long_t already_added;
  long int nodes_reached;
  
  igraph_dqueue_t q=IGRAPH_DQUEUE_NULL;
  igraph_vector_t *neis;
  igraph_integer_t dirmode;
  igraph_i_adjlist_t allneis;
  long int ressize;
  
  if (directed) { dirmode=IGRAPH_OUT; } else { dirmode=IGRAPH_ALL; }

  IGRAPH_CHECK(igraph_vector_long_init(&already_added, no_of_nodes));
  IGRAPH_FINALLY(igraph_vector_long_destroy, &already_added);
  IGRAPH_DQUEUE_INIT_FINALLY(&q, 100);
  IGRAPH_CHECK(igraph_i_adjlist_init(graph, &allneis, dirmode));
  IGRAPH_FINALLY(igraph_i_adjlist_destroy, &allneis);
  
  IGRAPH_CHECK(igraph_vector_resize(res, 0));
  ressize=0;
  
  *unconnected=0;

  for (i=0; i<no_of_nodes; i++) {
    nodes_reached=1;		/* itself */
    IGRAPH_CHECK(igraph_dqueue_push(&q, i));
    IGRAPH_CHECK(igraph_dqueue_push(&q, 0));
    VECTOR(already_added)[i]=i+1;
    
    igraph_progress("Path-hist: ", 100.0*i/no_of_nodes, NULL);

    IGRAPH_ALLOW_INTERRUPTION();
    
    while (!igraph_dqueue_empty(&q)) {
      long int actnode=igraph_dqueue_pop(&q);
      long int actdist=igraph_dqueue_pop(&q);
      
      neis=igraph_i_adjlist_get(&allneis, actnode);
      n=igraph_vector_size(neis);
      for (j=0; j<n; j++) {
	long int neighbor=VECTOR(*neis)[j];
	if (VECTOR(already_added)[neighbor] == i+1) { continue; }
	VECTOR(already_added)[neighbor] = i+1;
	nodes_reached++;
	if (actdist+1 > ressize) {
	  IGRAPH_CHECK(igraph_vector_resize(res, actdist+1));
	  for (; ressize<actdist+1; ressize++) {
	    VECTOR(*res)[ressize]=0;
	  }
	}
	VECTOR(*res)[actdist] += 1;

	IGRAPH_CHECK(igraph_dqueue_push(&q, neighbor));
	IGRAPH_CHECK(igraph_dqueue_push(&q, actdist+1));
      }
    } /* while !igraph_dqueue_empty */

    *unconnected += (no_of_nodes-nodes_reached);

  } /* for i<no_of_nodes */

  igraph_progress("Path-hist: ", 100.0, NULL);

  /* count every pair only once for an undirected graph */
  if (!directed || !igraph_is_directed(graph)) {
    for (i=1; i<ressize; i++) {
      VECTOR(*res)[i] /= 2;
    }
  }

  igraph_vector_long_destroy(&already_added);
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
 * </para><para>
 * If the graph has more minimum spanning trees (this is always the
 * case, except if it is a forest) this implementation returns only
 * the same one.
 * 
 * </para><para>
 * Directed graphs are considered as undirected for this computation.
 *
 * </para><para>
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
  long int no_of_edges=igraph_ecount(graph);
  char *already_added;
  char *added_edges;
  
  igraph_dqueue_t q=IGRAPH_DQUEUE_NULL;
  igraph_vector_t edges=IGRAPH_VECTOR_NULL;
  igraph_vector_t tmp=IGRAPH_VECTOR_NULL;
  long int i, j;

  added_edges=Calloc(no_of_edges, char);
  if (added_edges==0) {
    IGRAPH_ERROR("unweighted spanning tree failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, added_edges);
  IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
  already_added=Calloc(no_of_nodes, char);
  if (already_added==0) {
    IGRAPH_ERROR("unweighted spanning tree failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, already_added);
  IGRAPH_VECTOR_INIT_FINALLY(&tmp, 0);
  IGRAPH_DQUEUE_INIT_FINALLY(&q, 100);
  
  for (i=0; i<no_of_nodes; i++) {
    if (already_added[i]>0) { continue; }

    IGRAPH_ALLOW_INTERRUPTION();

    already_added[i]=1;
    IGRAPH_CHECK(igraph_dqueue_push(&q, i));
    while (! igraph_dqueue_empty(&q)) {
      long int act_node=igraph_dqueue_pop(&q);
      IGRAPH_CHECK(igraph_adjacent(graph, &tmp, act_node, IGRAPH_ALL));
      for (j=0; j<igraph_vector_size(&tmp); j++) {
	long int edge=VECTOR(tmp)[j];
	if (added_edges[edge]==0) {
	  igraph_integer_t from, to;
	  igraph_edge(graph, edge, &from, &to);
	  if (act_node==to) { to=from; }
	  if (already_added[(long int) to]==0) {
	    already_added[(long int) to]=1;
	    added_edges[edge]=1;
	    IGRAPH_CHECK(igraph_dqueue_push(&q, to));
	  }
	}
      }
    }
  }
  
  igraph_dqueue_destroy(&q);
  Free(already_added);
  igraph_vector_destroy(&tmp);
  IGRAPH_FINALLY_CLEAN(3);

  /* summarize the edges to delete */
  j=0;
  for (i=0; i<no_of_edges; i++) {
    if (added_edges[i]==0) {
      j++;
    }
  }
  IGRAPH_CHECK(igraph_vector_resize(&edges, j));
  j=0;
  for (i=0; i<no_of_edges; i++) {
    if (added_edges[i]==0) {
      VECTOR(edges)[j++]=i;
    }
  }
  
  IGRAPH_CHECK(igraph_copy(mst, graph));
  IGRAPH_FINALLY(igraph_destroy, mst);
  IGRAPH_CHECK(igraph_delete_edges(mst, igraph_ess_vector(&edges)));
  
  igraph_vector_destroy(&edges);
  Free(added_edges);
  IGRAPH_FINALLY_CLEAN(3);

  return 0;
}

/**
 * \ingroup structural
 * \function igraph_minimum_spanning_tree_prim
 * \brief Calculates one minimum spanning tree of a weighted graph.
 *
 * </para><para>
 * This function uses Prim's method for carrying out the computation,
 * see Prim, R.C.: Shortest connection networks and some
 * generalizations, Bell System Technical
 * Journal, Vol. 36, 
 * 1957, 1389--1401.
 * 
 * </para><para>
 * If the graph has more than one minimum spanning tree, the current
 * implementation returns always the same one.
 *
 * </para><para>
 * Directed graphs are considered as undirected for this computation. 
 * 
 * </para><para>
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
  long int no_of_edges=igraph_ecount(graph);
  char *already_added;
  char *added_edges;

  igraph_d_indheap_t heap=IGRAPH_D_INDHEAP_NULL;
  igraph_vector_t edges=IGRAPH_VECTOR_NULL;
  igraph_integer_t mode=IGRAPH_ALL;
  
  igraph_vector_t adj;

  long int i, j;

  if (igraph_vector_size(weights) != igraph_ecount(graph)) {
    IGRAPH_ERROR("Invalid weights length", IGRAPH_EINVAL);
  }

  added_edges=Calloc(no_of_edges, char);
  if (added_edges==0) {
    IGRAPH_ERROR("prim spanning tree failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_destroy, added_edges);
  IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
  already_added=Calloc(no_of_nodes, char);
  if (already_added == 0) {
    IGRAPH_ERROR("prim spanning tree failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, already_added);
  IGRAPH_CHECK(igraph_d_indheap_init(&heap, 0));
  IGRAPH_FINALLY(igraph_d_indheap_destroy, &heap);
  IGRAPH_VECTOR_INIT_FINALLY(&adj, 0);

  for (i=0; i<no_of_nodes; i++) {
    if (already_added[i]>0) { continue; }
    IGRAPH_ALLOW_INTERRUPTION();

    already_added[i]=1;
    /* add all edges of the first vertex */
    igraph_adjacent(graph, &adj, i, mode);
    for (j=0; j<igraph_vector_size(&adj); j++) {
      long int edgeno=VECTOR(adj)[j];
      igraph_integer_t edgefrom, edgeto;
      long int neighbor;
      igraph_edge(graph, edgeno, &edgefrom, &edgeto);
      neighbor= edgefrom != i ? edgefrom : edgeto;
      if (already_added[neighbor] == 0) {
	IGRAPH_CHECK(igraph_d_indheap_push(&heap, -VECTOR(*weights)[edgeno], i,
					   edgeno));
      }
    }

    while(! igraph_d_indheap_empty(&heap)) {
      /* Get minimal edge */
      long int from, edge;
      igraph_integer_t tmp, to;
      igraph_d_indheap_max_index(&heap, &from, &edge);
      igraph_edge(graph, edge, &tmp, &to);
      
      /* Erase it */
      igraph_d_indheap_delete_max(&heap);

      /* Is this edge already included? */
      if (added_edges[edge]==0) {
	if (from==to) { to=tmp; }
	/* Does it point to a visited node? */      
	if (already_added[(long int)to]==0) {
	  already_added[(long int)to]=1;
	  added_edges[edge]=1;
	  /* add all outgoing edges */
	  igraph_adjacent(graph, &adj, to, mode);
	  for (j=0; j<igraph_vector_size(&adj); j++) {
	    long int edgeno=VECTOR(adj)[j];
	    igraph_integer_t edgefrom, edgeto;
	    long int neighbor;
	    igraph_edge(graph, edgeno, &edgefrom, &edgeto);
	    neighbor= edgefrom != to ? edgefrom : edgeto;
	    if (already_added[neighbor] == 0) {
	      IGRAPH_CHECK(igraph_d_indheap_push(&heap, -VECTOR(*weights)[edgeno], to,
						 edgeno));
	    }
	  }
	} /* for */
      } /* if !already_added */
    } /* while in the same component */
  } /* for all nodes */

  igraph_d_indheap_destroy(&heap);
  Free(already_added);
  igraph_vector_destroy(&adj);
  IGRAPH_FINALLY_CLEAN(3);

  /* Ok, collect the edges to delete */
  j=0;
  for (i=0; i<no_of_edges; i++) {
    if (added_edges[i]==0) {
      j++;
    }
  }
  IGRAPH_CHECK(igraph_vector_resize(&edges, j));
  j=0;
  for (i=0; i<no_of_edges; i++) {
    if (added_edges[i]==0) {
      VECTOR(edges)[j++]=i;
    }
  }
  
  IGRAPH_CHECK(igraph_copy(mst, graph));
  IGRAPH_FINALLY(igraph_destroy, mst);
  IGRAPH_CHECK(igraph_delete_edges(mst, igraph_ess_vector(&edges)));
  
  igraph_vector_destroy(&edges);
  Free(added_edges);
  IGRAPH_FINALLY_CLEAN(3);
  
  return 0;
}

/**
 * \ingroup structural
 * \function igraph_closeness
 * \brief Closeness centrality calculations for some vertices.
 *
 * </para><para>
 * The closeness centrality of a vertex measures how easily other
 * vertices can be reached from it (or the other way: how easily it
 * can be reached from the other vertices). It is defined as the
 * number of the number of vertices minus one divided by the sum of the
 * lengths of all geodesics from/to the given vertex.
 *
 * </para><para>
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
		     const igraph_vs_t vids, 
		     igraph_neimode_t mode) {

  long int no_of_nodes=igraph_vcount(graph);
  igraph_vector_t already_counted;
  long int i, j;
  long int nodes_reached;

  igraph_dqueue_t q;
  
  long int nodes_to_calc;
  igraph_vector_t tmp;
  igraph_vit_t vit;

  IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
  IGRAPH_FINALLY(igraph_vit_destroy, &vit);

  nodes_to_calc=IGRAPH_VIT_SIZE(vit);
  
  if (mode != IGRAPH_OUT && mode != IGRAPH_IN && 
      mode != IGRAPH_ALL) {
    IGRAPH_ERROR("calculating closeness", IGRAPH_EINVMODE);
  }

  IGRAPH_VECTOR_INIT_FINALLY(&already_counted, no_of_nodes);
  IGRAPH_VECTOR_INIT_FINALLY(&tmp, 0);
  IGRAPH_DQUEUE_INIT_FINALLY(&q, 100);

  IGRAPH_CHECK(igraph_vector_resize(res, nodes_to_calc));
  igraph_vector_null(res);
  
  for (IGRAPH_VIT_RESET(vit), i=0; 
       !IGRAPH_VIT_END(vit); 
       IGRAPH_VIT_NEXT(vit), i++) {
    IGRAPH_CHECK(igraph_dqueue_push(&q, IGRAPH_VIT_GET(vit)));
    IGRAPH_CHECK(igraph_dqueue_push(&q, 0));
    nodes_reached=1;
    VECTOR(already_counted)[(long int)IGRAPH_VIT_GET(vit)]=i+1;

    IGRAPH_ALLOW_INTERRUPTION();
    
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
  igraph_vit_destroy(&vit);
  IGRAPH_FINALLY_CLEAN(4);
  
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
			  const igraph_vs_t from, igraph_neimode_t mode) {

  long int no_of_nodes=igraph_vcount(graph);
  long int no_of_from;
  long int *already_counted;
  
  igraph_dqueue_t q=IGRAPH_DQUEUE_NULL;

  long int i, j;
  igraph_vector_t tmp=IGRAPH_VECTOR_NULL;
  igraph_vit_t fromvit;

  IGRAPH_CHECK(igraph_vit_create(graph, from, &fromvit));
  IGRAPH_FINALLY(igraph_vit_destroy, &fromvit);

  no_of_from=IGRAPH_VIT_SIZE(fromvit);

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

  for (IGRAPH_VIT_RESET(fromvit), i=0; 
       !IGRAPH_VIT_END(fromvit); 
       IGRAPH_VIT_NEXT(fromvit), i++) {
    long int reached=1;
    IGRAPH_CHECK(igraph_dqueue_push(&q, IGRAPH_VIT_GET(fromvit)));
    IGRAPH_CHECK(igraph_dqueue_push(&q, 0));
    already_counted[ (long int) IGRAPH_VIT_GET(fromvit) ] = i+1;
    
    IGRAPH_ALLOW_INTERRUPTION();

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
      if (MATRIX(*res, i, j) == 0 && j != IGRAPH_VIT_GET(fromvit)) {
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
  igraph_vit_destroy(&fromvit);
  IGRAPH_FINALLY_CLEAN(4);

  return 0;
}

/**
 * \ingroup structural
 * \function igraph_get_shortest_paths
 * \brief Calculates the shortest paths from/to one vertex.
 * 
 * </para><para>
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
			      igraph_integer_t from, const igraph_vs_t to, 
			      igraph_neimode_t mode) {

  /* TODO: use adjlist_t if to is long (longer than 1?) */

  long int no_of_nodes=igraph_vcount(graph);
  long int *father;
  
  igraph_dqueue_t q=IGRAPH_DQUEUE_NULL;

  long int j;
  igraph_vector_t tmp=IGRAPH_VECTOR_NULL;

  igraph_vit_t vit;
  
  long int to_reach;
  long int reached=0;

  if (from<0 || from>=no_of_nodes) {
    IGRAPH_ERROR("cannot get shortest paths", IGRAPH_EINVVID);
  }
  if (mode != IGRAPH_OUT && mode != IGRAPH_IN && 
      mode != IGRAPH_ALL) {
    IGRAPH_ERROR("Invalid mode argument", IGRAPH_EINVMODE);
  }

  IGRAPH_CHECK(igraph_vit_create(graph, to, &vit));
  IGRAPH_FINALLY(igraph_vit_destroy, &vit);

  if (IGRAPH_VIT_SIZE(vit) != igraph_vector_ptr_size(res)) {
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
  to_reach=IGRAPH_VIT_SIZE(vit);
  for (IGRAPH_VIT_RESET(vit); !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit)) {
    if (father[ (long int) IGRAPH_VIT_GET(vit) ] == 0) {
      father[ (long int) IGRAPH_VIT_GET(vit) ] = -1;
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
  
  for (IGRAPH_VIT_RESET(vit), j=0; 
       !IGRAPH_VIT_END(vit);
       IGRAPH_VIT_NEXT(vit), j++) {
    long int node=IGRAPH_VIT_GET(vit);
    igraph_vector_t *vec=VECTOR(*res)[j];
    igraph_vector_clear(vec);

    IGRAPH_ALLOW_INTERRUPTION();

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
  igraph_vit_destroy(&vit);
  IGRAPH_FINALLY_CLEAN(4);

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
				  igraph_integer_t from, const igraph_vs_t to,
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

  igraph_vit_t vit;

  if (from<0 || from>=no_of_nodes) {
    IGRAPH_ERROR("cannot get shortest paths", IGRAPH_EINVVID);
  }
  if (mode != IGRAPH_OUT && mode != IGRAPH_IN && 
      mode != IGRAPH_ALL) {
    IGRAPH_ERROR("Invalid mode argument", IGRAPH_EINVMODE);
  }

  IGRAPH_CHECK(igraph_vit_create(graph, to, &vit));
  IGRAPH_FINALLY(igraph_vit_destroy, &vit);
  
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
    
    IGRAPH_ALLOW_INTERRUPTION();

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
  for (IGRAPH_VIT_RESET(vit); !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit)) {
    geodist[ (long int) IGRAPH_VIT_GET(vit) ] = 1;
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

    IGRAPH_ALLOW_INTERRUPTION();

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
  igraph_vit_destroy(&vit);
  IGRAPH_FINALLY_CLEAN(6);

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

int igraph_subcomponent(const igraph_t *graph, igraph_vector_t *res, igraph_real_t vertex, 
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

    IGRAPH_ALLOW_INTERRUPTION();

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
 * </para><para>
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
			const igraph_vs_t vids, 
			igraph_bool_t directed) {

  long int no_of_nodes=igraph_vcount(graph);
  igraph_dqueue_t q=IGRAPH_DQUEUE_NULL;
  long int *distance;
  long int *nrgeo;
  double *tmpscore;
  igraph_stack_t stack=IGRAPH_STACK_NULL;
  long int source;
  long int j, k;
  igraph_vector_t tmp=IGRAPH_VECTOR_NULL;
  igraph_integer_t modein, modeout;
  igraph_vit_t vit;

  IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
  IGRAPH_FINALLY(igraph_vit_destroy, &vit);

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
    
  IGRAPH_CHECK(igraph_vector_resize(res, IGRAPH_VIT_SIZE(vit)));
  igraph_vector_null(res);

  /* here we go */
  
  for (source=0; source<no_of_nodes; source++) {
    
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
    for (k=0, IGRAPH_VIT_RESET(vit); !IGRAPH_VIT_END(vit); 
	 IGRAPH_VIT_NEXT(vit), k++) {
      long int node=IGRAPH_VIT_GET(vit);
      VECTOR(*res)[k] += tmpscore[node];
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
  igraph_vit_destroy(&vit);
  IGRAPH_FINALLY_CLEAN(7);

  return 0;
}

/**
 * \ingroup structural
 * \function igraph_edge_betweenness
 * \brief Betweenness centrality of the edges.
 * 
 * </para><para>
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
			     igraph_bool_t directed) {
  
  long int no_of_nodes=igraph_vcount(graph);
  long int no_of_edges=igraph_ecount(graph);
  igraph_dqueue_t q=IGRAPH_DQUEUE_NULL;
  long int *distance;
  long int *nrgeo;
  double *tmpscore;
  igraph_stack_t stack=IGRAPH_STACK_NULL;
  long int source;
  long int j;

  igraph_i_adjedgelist_t elist_out, elist_in;
  igraph_i_adjedgelist_t *elist_out_p, *elist_in_p;
  igraph_vector_t *neip;
  long int neino;
  long int i;
  igraph_integer_t modein, modeout;

  directed=directed && igraph_is_directed(graph);
  if (directed) {
    modeout=IGRAPH_OUT;
    modein=IGRAPH_IN;
    IGRAPH_CHECK(igraph_i_adjedgelist_init(graph, &elist_out, IGRAPH_OUT));
    IGRAPH_FINALLY(igraph_i_adjedgelist_destroy, &elist_out);
    IGRAPH_CHECK(igraph_i_adjedgelist_init(graph, &elist_in, IGRAPH_IN));
    IGRAPH_FINALLY(igraph_i_adjedgelist_destroy, &elist_in);
    elist_out_p=&elist_out;
    elist_in_p=&elist_in;
  } else {
    modeout=modein=IGRAPH_ALL;
    IGRAPH_CHECK(igraph_i_adjedgelist_init(graph,&elist_out, IGRAPH_ALL));
    IGRAPH_FINALLY(igraph_i_adjedgelist_destroy, &elist_out);
    elist_out_p=elist_in_p=&elist_out;
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

  IGRAPH_CHECK(igraph_vector_resize(result, no_of_edges));

  igraph_vector_null(result);

  /* here we go */
  
  for (source=0; source<no_of_nodes; source++) {

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
	  VECTOR(*result)[edgeno] +=
	    (tmpscore[actnode]+1)*nrgeo[neighbor]/nrgeo[actnode];
	}
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
  IGRAPH_FINALLY_CLEAN(5);

  if (directed) {
    igraph_i_adjedgelist_destroy(&elist_out);
    igraph_i_adjedgelist_destroy(&elist_in);
    IGRAPH_FINALLY_CLEAN(2);
  } else {
    igraph_i_adjedgelist_destroy(&elist_out);
    IGRAPH_FINALLY_CLEAN(1);
  }

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
 * </para><para>
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
 * 
 * Time complexity: TODO.
 */

int igraph_pagerank(const igraph_t *graph, igraph_vector_t *res, 
		    const igraph_vs_t vids, igraph_bool_t directed,
            igraph_integer_t niter, igraph_real_t eps, igraph_real_t damping) {
  long int no_of_nodes=igraph_vcount(graph);
  long int i, j, n, nodes_to_calc;
  igraph_real_t *prvec, *prvec_new, *prvec_aux, *prvec_scaled;
  igraph_vector_t *neis, outdegree;
  igraph_integer_t dirmode;
  igraph_i_adjlist_t allneis;
  igraph_real_t maxdiff=eps;
  igraph_vit_t vit;

  if (niter<=0) IGRAPH_ERROR("Invalid iteration count", IGRAPH_EINVAL);
  if (eps<=0) IGRAPH_ERROR("Invalid epsilon value", IGRAPH_EINVAL);
  if (damping<=0 || damping>=1) IGRAPH_ERROR("Invalid damping factor", IGRAPH_EINVAL);

  IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
  IGRAPH_FINALLY(igraph_vit_destroy, &vit);
  nodes_to_calc=IGRAPH_VIT_SIZE(vit);

  IGRAPH_CHECK(igraph_vector_resize(res, nodes_to_calc));
  igraph_vector_null(res);
  
  IGRAPH_VECTOR_INIT_FINALLY(&outdegree, no_of_nodes);
    
  prvec=Calloc(no_of_nodes, igraph_real_t);
  if (prvec==0) {
    IGRAPH_ERROR("pagerank failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, prvec);
  
  prvec_new=Calloc(no_of_nodes, igraph_real_t);
  if (prvec_new==0) {
    IGRAPH_ERROR("pagerank failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, prvec_new);
  
  prvec_scaled=Calloc(no_of_nodes, igraph_real_t);
  if (prvec_scaled==0) {
    IGRAPH_ERROR("pagerank failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, prvec_scaled);
  
  if (directed) { dirmode=IGRAPH_IN; } else { dirmode=IGRAPH_ALL; }  
  igraph_i_adjlist_init(graph, &allneis, dirmode);
  IGRAPH_FINALLY(igraph_i_adjlist_destroy, &allneis);

  /* Calculate outdegrees for every node */
  igraph_degree(graph, &outdegree, igraph_vss_all(),
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
      
      IGRAPH_ALLOW_INTERRUPTION();

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
  for (IGRAPH_VIT_RESET(vit), i=0; 
       !IGRAPH_VIT_END(vit); 
       IGRAPH_VIT_NEXT(vit), i++) {
    long int vid=IGRAPH_VIT_GET(vit);
    VECTOR(*res)[i]=prvec[vid];
  }
  
  igraph_i_adjlist_destroy(&allneis);
  igraph_vit_destroy(&vit);
  igraph_vector_destroy(&outdegree);
  Free(prvec);
  Free(prvec_new);  
  Free(prvec_scaled);
  
  IGRAPH_FINALLY_CLEAN(6);
  
  return 0;
}

/**
 * \ingroup structural
 * \function igraph_rewire
 * \brief Randomly rewires a graph while preserving the degree distribution.
 * 
 * </para><para>
 * This function generates a new graph based on the original one by randomly
 * rewiring edges while preserving the original graph's degree distribution.
 * Please note that the rewiring is done "in place", so no new graph will
 * be allocated. If you would like to keep the original graph intact, use
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
 *
 * Time complexity: TODO.
 */

int igraph_rewire(igraph_t *graph, igraph_integer_t n, igraph_rewiring_t mode) {
  long int no_of_nodes=igraph_vcount(graph);
  long int i, a, b, c, d;
  igraph_i_adjlist_t allneis;
  igraph_vector_t *neis[2], edgevec;
  igraph_bool_t directed;
  igraph_es_t es;
  
  if (mode == IGRAPH_REWIRING_SIMPLE && no_of_nodes<4)
    IGRAPH_ERROR("graph unsuitable for rewiring", IGRAPH_EINVAL);
  
  directed = igraph_is_directed(graph);
  
  RNG_BEGIN();
  
  igraph_i_adjlist_init(graph, &allneis, IGRAPH_OUT);
  IGRAPH_FINALLY(igraph_i_adjlist_destroy, &allneis);
  igraph_vector_init(&edgevec, 4);
  IGRAPH_FINALLY(igraph_vector_destroy, &edgevec);

  while (n>0) {
    
    IGRAPH_ALLOW_INTERRUPTION();
    
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
	/* printf("Deleting: %d -> %d, %d -> %d\n", a, b, c, d); */
	IGRAPH_CHECK(igraph_es_pairs_small(&es, directed, a, b, c, d, -1));
	IGRAPH_FINALLY(igraph_es_destroy, &es);
	IGRAPH_CHECK(igraph_delete_edges(graph, es));	
	igraph_es_destroy(&es);
	IGRAPH_FINALLY_CLEAN(1);
	VECTOR(edgevec)[0]=a; VECTOR(edgevec)[1]=d;
	VECTOR(edgevec)[2]=c; VECTOR(edgevec)[3]=b;
	/* printf("Adding: %d -> %d, %d -> %d\n", a, d, c, b); */
	igraph_add_edges(graph, &edgevec, 0);
	/* We have to adjust the adjacency list view as well.
	   It's a luck that we have the pointers in neis[0] and neis[1] */
	for (i=igraph_vector_size(neis[0])-1; i>=0; i--)
	  if (VECTOR(*neis[0])[i]==b) { VECTOR(*neis[0])[i]=d; break; }
	for (i=igraph_vector_size(neis[1])-1; i>=0; i--)
	  if (VECTOR(*neis[1])[i]==d) { VECTOR(*neis[1])[i]=b; break; }
	/* In case of an undirected graph, we have to adjust the
	 * adjacency view of vertices b and d as well */
	if (!directed) {
	  neis[0] = igraph_i_adjlist_get(&allneis, b);
	  neis[1] = igraph_i_adjlist_get(&allneis, d);
	  for (i=igraph_vector_size(neis[0])-1; i>=0; i--)
	    if (VECTOR(*neis[0])[i]==a) { VECTOR(*neis[0])[i]=c; break; }
	  for (i=igraph_vector_size(neis[1])-1; i>=0; i--)
	    if (VECTOR(*neis[1])[i]==c) { VECTOR(*neis[1])[i]=a; break; }
	}
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
 * </para><para>
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
		    const igraph_vs_t vids) {
  
  long int no_of_nodes=igraph_vcount(graph);
  igraph_vector_t delete=IGRAPH_VECTOR_NULL;
  char *remain;
  long int i;
  igraph_vit_t vit;
  
  IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
  IGRAPH_FINALLY(igraph_vit_destroy, &vit);

  IGRAPH_VECTOR_INIT_FINALLY(&delete, 0);
  remain=Calloc(no_of_nodes, char);
  if (remain==0) {
    IGRAPH_ERROR("subgraph failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(free, remain);	/* TODO: hack */
  IGRAPH_CHECK(igraph_vector_reserve(&delete, no_of_nodes-IGRAPH_VIT_SIZE(vit)));
  
  for (IGRAPH_VIT_RESET(vit); !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit)) {
    remain[ (long int) IGRAPH_VIT_GET(vit) ] = 1;
  }

  for (i=0; i<no_of_nodes; i++) {

    IGRAPH_ALLOW_INTERRUPTION();

    if (remain[i] == 0) {
      IGRAPH_CHECK(igraph_vector_push_back(&delete, i));
    }
  }

  Free(remain);
  IGRAPH_FINALLY_CLEAN(1);
  
  /* must set res->attr to 0 before calling igraph_copy */
  res->attr=0;           /* Why is this needed? TODO */
  IGRAPH_CHECK(igraph_copy(res, graph));
  IGRAPH_FINALLY(igraph_destroy, res);
  IGRAPH_CHECK(igraph_delete_vertices(res, igraph_vss_vector(&delete)));
  
  igraph_vector_destroy(&delete);
  igraph_vit_destroy(&vit);
  IGRAPH_FINALLY_CLEAN(3);
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
 * Time complexity: O(|V|+|E|).
 */

int igraph_simplify(igraph_t *graph, igraph_bool_t multiple, igraph_bool_t loops) {

  igraph_vector_t edges=IGRAPH_VECTOR_NULL;
  igraph_vector_t neis=IGRAPH_VECTOR_NULL;
  long int no_of_nodes=igraph_vcount(graph);
  long int i, j;
  igraph_es_t es;
  igraph_bool_t directed=igraph_is_directed(graph);

  IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);

  if (directed) { 
    for (i=0; i<no_of_nodes; i++) {
      IGRAPH_CHECK(igraph_neighbors(graph, &neis, i, IGRAPH_OUT));
      
      IGRAPH_ALLOW_INTERRUPTION();
      
      if (loops) {
	for (j=0; j<igraph_vector_size(&neis); j++) {
	  if (VECTOR(neis)[j]==i) {
	    IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
	    IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
	  }
	}
      } /* if loops */
      
      if (multiple) {
	for (j=1; j<igraph_vector_size(&neis); j++) {
	  if (VECTOR(neis)[j]==VECTOR(neis)[j-1] && 
	      (!loops || VECTOR(neis)[j] != i) ) {
	    IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
	    IGRAPH_CHECK(igraph_vector_push_back(&edges, VECTOR(neis)[j]));
	  }
	}
      }
    }
  } else { 			/* not directed */
    for (i=0; i<no_of_nodes; i++) {
      int flip=0;
      IGRAPH_CHECK(igraph_neighbors(graph, &neis, i, IGRAPH_OUT));
      
      IGRAPH_ALLOW_INTERRUPTION();
      
      if (loops) {
	for (j=0; j<igraph_vector_size(&neis); j++) {
	  if (VECTOR(neis)[j]==i) {
	    flip=1-flip;
	    if (flip==0) {
	      IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
	      IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
	    }
	  }
	}
      } /* if loops */
      
      if (multiple) {
	for (j=1; j<igraph_vector_size(&neis); j++) {
	  if (VECTOR(neis)[j] > i && VECTOR(neis)[j]==VECTOR(neis)[j-1]) {
	    IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
	    IGRAPH_CHECK(igraph_vector_push_back(&edges, VECTOR(neis)[j]));
	  }
	  if (VECTOR(neis)[j]==i && VECTOR(neis)[j-1]==i && !loops) {
	    flip=1-flip;
	    if (flip==0) {
	      IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
	      IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
	    }
	  }
	}
      }
    }
  }
    
  igraph_vector_destroy(&neis);
  IGRAPH_FINALLY_CLEAN(1);
  IGRAPH_CHECK(igraph_es_multipairs(&es, &edges, IGRAPH_DIRECTED));
  IGRAPH_FINALLY(igraph_es_destroy, &es);
  IGRAPH_CHECK(igraph_delete_edges(graph, es));
  igraph_es_destroy(&es);
  IGRAPH_FINALLY_CLEAN(1);
  igraph_vector_destroy(&edges);
  IGRAPH_FINALLY_CLEAN(1);

  return 0;
}

/**
 * \function igraph_transitivity_avglocal_undirected
 * \brief Average local transitivity (clustering coefficient)
 * 
 * The transitivity measures the probability that two neighbors of a
 * vertex are connected. In case of the average local transitivity
 * this probability if calculated for each vertex and then the average
 * is taken for those vertices which have at least two neighbors. If
 * there are no such vertices then \c NaN is returned. 
 * \param graph The input graph, directed graphs are considered as 
 *    undirected ones.
 * \param res Pointer to a real variable, the result will be stored here.
 * \return Error code.
 * 
 * \sa \ref igraph_transitivity_undirected(), \ref
 * igraph_transitivity_local_undirected().
 * 
 * Time complexity: O(|V|*d^2), |V| is the number of vertices in the
 * graph and d is the average degree.
 */

int igraph_transitivity_avglocal_undirected(const igraph_t *graph,
					    igraph_real_t *res) {

  long int no_of_nodes=igraph_vcount(graph);
  igraph_real_t sum=0.0;
  igraph_integer_t count=0;
  long int node, i, j, nn;
  igraph_i_adjlist_t allneis;
  igraph_vector_t *neis1, *neis2;
  long int neilen1, neilen2;
  igraph_integer_t triples;
  long int *neis;
  long int maxdegree;

  igraph_vector_t order;
  igraph_vector_t rank;
  igraph_vector_t degree;
  igraph_vector_t triangles;

  IGRAPH_VECTOR_INIT_FINALLY(&order, no_of_nodes);
  IGRAPH_VECTOR_INIT_FINALLY(&degree, no_of_nodes);
  
  IGRAPH_CHECK(igraph_degree(graph, &degree, igraph_vss_all(), IGRAPH_ALL,
			     IGRAPH_LOOPS));
  maxdegree=igraph_vector_max(&degree)+1;
  igraph_vector_order1(&degree, &order, maxdegree);
  igraph_vector_destroy(&degree);
  IGRAPH_FINALLY_CLEAN(1);
  IGRAPH_VECTOR_INIT_FINALLY(&rank, no_of_nodes);
  for (i=0; i<no_of_nodes; i++) {
    VECTOR(rank)[ (long int) VECTOR(order)[i] ] = no_of_nodes-i-1;
  }
  
  IGRAPH_CHECK(igraph_i_adjlist_init(graph, &allneis, IGRAPH_ALL));
  IGRAPH_FINALLY(igraph_i_adjlist_destroy, &allneis);
  IGRAPH_CHECK(igraph_i_adjlist_simplify(&allneis));

  neis=Calloc(no_of_nodes, long int);
  if (neis==0) {
    IGRAPH_ERROR("undirected average local transitivity failed",
		 IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, neis);

  IGRAPH_VECTOR_INIT_FINALLY(&triangles, no_of_nodes);
  
  for (nn=no_of_nodes-1; nn >= 0; nn--) {
    node=VECTOR(order)[nn];

    IGRAPH_ALLOW_INTERRUPTION();
    
    neis1=igraph_i_adjlist_get(&allneis, node);
    neilen1=igraph_vector_size(neis1);
    triples = (double)neilen1 * (neilen1-1) / 2;
    /* Mark the neighbors of 'node' */
    for (i=0; i<neilen1; i++) {
      neis[ (long int)VECTOR(*neis1)[i] ] = node+1;
    }
    
    for (i=0; i<neilen1; i++) {
      long int nei=VECTOR(*neis1)[i];
      if (VECTOR(rank)[nei] > VECTOR(rank)[node]) {
	neis2=igraph_i_adjlist_get(&allneis, nei);
	neilen2=igraph_vector_size(neis2);
	for (j=0; j<neilen2; j++) {
	  long int nei2=VECTOR(*neis2)[j];
	  if (VECTOR(rank)[nei2] < VECTOR(rank)[nei]) {
	    continue;
	  }
	  if (neis[nei2] == node+1) {
	    VECTOR(triangles)[nei2] += 1;
	    VECTOR(triangles)[nei] += 1;
	    VECTOR(triangles)[node] += 1;
	  }
	}
      }
    }
    
    if (triples != 0) {
      sum += VECTOR(triangles)[node] / triples;
      count++;
    }
  }
  
  *res = sum/count;

  igraph_vector_destroy(&triangles);
  Free(neis);
  igraph_i_adjlist_destroy(&allneis);
  igraph_vector_destroy(&rank);
  igraph_vector_destroy(&order);
  IGRAPH_FINALLY_CLEAN(5);
  return 0;
}
  
int igraph_transitivity_local_undirected1(const igraph_t *graph, 
					  igraph_vector_t *res,
					  const igraph_vs_t vids) {
  
  long int no_of_nodes=igraph_vcount(graph);
  igraph_vit_t vit;
  long int nodes_to_calc;
  igraph_vector_t *neis1, *neis2;
  igraph_real_t triples, triangles;
  long int i, j, k;
  long int neilen1, neilen2;
  long int *neis;
  igraph_i_lazy_adjlist_t adjlist;

  IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
  IGRAPH_FINALLY(igraph_vit_destroy, &vit);
  nodes_to_calc=IGRAPH_VIT_SIZE(vit);

  neis=Calloc(no_of_nodes, long int);
  if (neis==0) {
    IGRAPH_ERROR("local undirected transitivity failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, neis);

  IGRAPH_CHECK(igraph_vector_resize(res, nodes_to_calc));

  igraph_i_lazy_adjlist_init(graph, &adjlist, IGRAPH_ALL, IGRAPH_I_SIMPLIFY);
  IGRAPH_FINALLY(igraph_i_lazy_adjlist_destroy, &adjlist);  

  for (i=0; !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit), i++) {
    long int node=IGRAPH_VIT_GET(vit);
    
    IGRAPH_ALLOW_INTERRUPTION();
    
    neis1=igraph_i_lazy_adjlist_get(&adjlist, node);
    neilen1=igraph_vector_size(neis1);
    for (j=0; j<neilen1; j++) {
      neis[ (long int)VECTOR(*neis1)[j] ] = i+1;
    }
    triples = (double)neilen1*(neilen1-1);
    triangles = 0;

    for (j=0; j<neilen1; j++) {
      long int v=VECTOR(*neis1)[j];
      neis2=igraph_i_lazy_adjlist_get(&adjlist, v);
      neilen2=igraph_vector_size(neis2);
      for (k=0; k<neilen2; k++) {
	long int v2=VECTOR(*neis2)[k];
	if (neis[v2] == i+1) {
	  triangles += 1.0;
	}
      }
    }
    VECTOR(*res)[i] = triangles/triples;
/*     fprintf(stderr, "%f %f\n", triangles, triples); */
  }

  igraph_i_lazy_adjlist_destroy(&adjlist);
  Free(neis);
  igraph_vit_destroy(&vit);
  IGRAPH_FINALLY_CLEAN(3);
  return 0;
}

int igraph_transitivity_local_undirected2(const igraph_t *graph, 
					  igraph_vector_t *res,
					  const igraph_vs_t vids) {
  
  long int no_of_nodes=igraph_vcount(graph);
  igraph_vit_t vit;
  long int nodes_to_calc, affected_nodes;
  long int maxdegree=0;
  long int i, j, k, nn;
  igraph_i_lazy_adjlist_t adjlist;
  igraph_vector_t index, avids, rank, order, triangles, degree;
  long int *neis;

  IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
  IGRAPH_FINALLY(igraph_vit_destroy, &vit);
  nodes_to_calc=IGRAPH_VIT_SIZE(vit);

  IGRAPH_CHECK(igraph_i_lazy_adjlist_init(graph, &adjlist, IGRAPH_ALL,
					  IGRAPH_I_SIMPLIFY));
  IGRAPH_FINALLY(igraph_i_lazy_adjlist_destroy, &adjlist);

  IGRAPH_VECTOR_INIT_FINALLY(&index, no_of_nodes);
  IGRAPH_VECTOR_INIT_FINALLY(&avids, 0);
  IGRAPH_CHECK(igraph_vector_reserve(&avids, nodes_to_calc));
  k=0;
  for (i=0; i<nodes_to_calc; IGRAPH_VIT_NEXT(vit), i++) {
    long int v=IGRAPH_VIT_GET(vit);
    igraph_vector_t *neis;
    long int neilen;
    if (VECTOR(index)[v]==0) {
      VECTOR(index)[v]=k+1; k++;
      IGRAPH_CHECK(igraph_vector_push_back(&avids, v));
    } 
    
    neis=igraph_i_lazy_adjlist_get(&adjlist, v);
    neilen=igraph_vector_size(neis);
    for (j=0; j<neilen; j++) {
      long int nei=VECTOR(*neis)[j];
      if (VECTOR(index)[nei]==0) {
	VECTOR(index)[nei]=k+1; k++;
	IGRAPH_CHECK(igraph_vector_push_back(&avids, nei));
      }
    }
  }

  /* Degree, ordering, ranking */
  affected_nodes=igraph_vector_size(&avids);
  IGRAPH_VECTOR_INIT_FINALLY(&order, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&degree, affected_nodes);
  for (i=0; i<affected_nodes; i++) {
    long int v=VECTOR(avids)[i];
    igraph_vector_t *neis;
    long int deg;
    neis=igraph_i_lazy_adjlist_get(&adjlist, v);
    VECTOR(degree)[i]=deg=igraph_vector_size(neis);
    if (deg > maxdegree) { maxdegree = deg; }
  }
  igraph_vector_order1(&degree, &order, maxdegree+1);
  igraph_vector_destroy(&degree);
  IGRAPH_FINALLY_CLEAN(1);
  IGRAPH_VECTOR_INIT_FINALLY(&rank, affected_nodes);
  for (i=0; i<affected_nodes; i++) {
    VECTOR(rank)[ (long int) VECTOR(order)[i] ] = affected_nodes-i-1;
  }
  
  neis=Calloc(no_of_nodes, long int);
  if (neis==0) {
    IGRAPH_ERROR("local transitivity calculation failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, neis);
  
  IGRAPH_VECTOR_INIT_FINALLY(&triangles, affected_nodes);
  for (nn=affected_nodes-1; nn>=0; nn--) {
    long int node=VECTOR(avids) [ (long int) VECTOR(order)[nn] ];
    igraph_vector_t *neis1, *neis2;
    long int neilen1, neilen2;
    long int nodeindex=VECTOR(index)[node];
    long int noderank=VECTOR(rank) [nodeindex-1];
    
/*     fprintf(stderr, "node %li (index %li, rank %li)\n", node, */
/* 	    (long int)VECTOR(index)[node]-1, noderank); */
    
    IGRAPH_ALLOW_INTERRUPTION();
    
    neis1=igraph_i_lazy_adjlist_get(&adjlist, node);
    neilen1=igraph_vector_size(neis1);
    for (i=0; i<neilen1; i++) {
      long int nei=VECTOR(*neis1)[i];
      neis[nei] = node+1;
    }
    for (i=0; i<neilen1; i++) {
      long int nei=VECTOR(*neis1)[i];
      long int neiindex=VECTOR(index)[nei];
      long int neirank=VECTOR(rank)[neiindex-1];

/*       fprintf(stderr, "  nei %li (index %li, rank %li)\n", nei, */
/* 	      neiindex, neirank); */
      if (neirank > noderank) {
	neis2=igraph_i_lazy_adjlist_get(&adjlist, nei);
	neilen2=igraph_vector_size(neis2);
	for (j=0; j<neilen2; j++) {	  
	  long int nei2=VECTOR(*neis2)[j];
	  long int nei2index=VECTOR(index)[nei2];
	  long int nei2rank=VECTOR(rank)[nei2index-1];
/* 	  fprintf(stderr, "    triple %li %li %li\n", node, nei, nei2); */
	  if (nei2rank < neirank) {
	    continue;
	  } 
	  if (neis[nei2] == node+1) {
/* 	    fprintf(stderr, "    triangle\n"); */
	    VECTOR(triangles) [ nei2index-1 ] += 1;
	    VECTOR(triangles) [ neiindex-1 ] += 1;
	    VECTOR(triangles) [ nodeindex-1 ] += 1;
	  }
	}
      }
    }    
  }
  
  /* Ok, for all affected vertices the number of triangles were counted */
  
  IGRAPH_CHECK(igraph_vector_resize(res, nodes_to_calc));
  IGRAPH_VIT_RESET(vit);
  for (i=0; i<nodes_to_calc; i++, IGRAPH_VIT_NEXT(vit)) {
    long int node=IGRAPH_VIT_GET(vit);
    long int idx=VECTOR(index)[node]-1;
    igraph_vector_t *neis=igraph_i_lazy_adjlist_get(&adjlist, node);
    long int deg=igraph_vector_size(neis);
    igraph_real_t triples=(double) deg * (deg-1) / 2;
    VECTOR(*res)[i] = VECTOR(triangles)[idx] / triples;
/*     fprintf(stderr, "%f %f\n", VECTOR(triangles)[idx], triples); */
  }
  
  igraph_vector_destroy(&triangles);
  igraph_free(neis);
  igraph_vector_destroy(&rank);
  igraph_vector_destroy(&order);
  igraph_vector_destroy(&avids);
  igraph_vector_destroy(&index);
  igraph_i_lazy_adjlist_destroy(&adjlist);
  igraph_vit_destroy(&vit);
  IGRAPH_FINALLY_CLEAN(8);

  return 0;
}

/* We don't use this, it is theoretically good, but practically not.
 */

/* int igraph_transitivity_local_undirected3(const igraph_t *graph, */
/* 				      igraph_vector_t *res, */
/* 				      const igraph_vs_t vids) { */

/*   igraph_vit_t vit; */
/*   long int nodes_to_calc; */
/*   igraph_i_lazy_adjlist_t adjlist; */
/*   long int i, j; */
  
/*   IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit)); */
/*   IGRAPH_FINALLY(igraph_vit_destroy, &vit); */
/*   nodes_to_calc=IGRAPH_VIT_SIZE(vit); */
  
/*   IGRAPH_CHECK(igraph_i_lazy_adjlist_init(graph, &adjlist, IGRAPH_ALL, */
/* 					  IGRAPH_I_SIMPLIFY)); */
/*   IGRAPH_FINALLY(igraph_i_lazy_adjlist_destroy, &adjlist); */
  
/*   IGRAPH_CHECK(igraph_vector_resize(res, nodes_to_calc)); */
/*   for (i=0, IGRAPH_VIT_RESET(vit); !IGRAPH_VIT_END(vit);  */
/*        i++, IGRAPH_VIT_NEXT(vit)) { */
/*     long int node=IGRAPH_VIT_GET(vit); */
/*     igraph_vector_t *neis=igraph_i_lazy_adjlist_get(&adjlist, node); */
/*     long int n1=igraph_vector_size(neis); */
/*     igraph_real_t triangles=0; */
/*     igraph_real_t triples=(double)n1*(n1-1); */
/*     IGRAPH_ALLOW_INTERRUPTION(); */
/*     for (j=0; j<n1; j++) { */
/*       long int node2=VECTOR(*neis)[j]; */
/*       igraph_vector_t *neis2=igraph_i_lazy_adjlist_get(&adjlist, node2); */
/*       long int n2=igraph_vector_size(neis2); */
/*       long int l1=0, l2=0; */
/*       while (l1 < n1 && l2 < n2) { */
/* 	long int nei1=VECTOR(*neis)[l1]; */
/* 	long int nei2=VECTOR(*neis2)[l2]; */
/* 	if (nei1 < nei2) {  */
/* 	  l1++; */
/* 	} else if (nei1 > nei2) { */
/* 	  l2++; */
/* 	} else { */
/* 	  triangles+=1; */
/* 	  l1++; l2++; */
/* 	} */
/*       } */
/*     } */
/*     /\* We're done with 'node' *\/ */
/*     VECTOR(*res)[i] = triangles / triples;   */
/*   } */

/*   igraph_i_lazy_adjlist_destroy(&adjlist); */
/*   igraph_vit_destroy(&vit); */
/*   IGRAPH_FINALLY_CLEAN(2); */

/*   return 0; */
/* } */

int igraph_transitivity_local_undirected4(const igraph_t *graph,
					  igraph_vector_t *res,
					  const igraph_vs_t vids) {

  long int no_of_nodes=igraph_vcount(graph);
  long int node, i, j, nn;
  igraph_i_adjlist_t allneis;
  igraph_vector_t *neis1, *neis2;
  long int neilen1, neilen2;
  igraph_integer_t triples;
  long int *neis;
  long int maxdegree;

  igraph_vector_t order;
  igraph_vector_t rank;
  igraph_vector_t degree;
  
  if (!igraph_vs_is_all((igraph_vs_t*)&vids)) {
    IGRAPH_ERROR("Internal error, wrong transitivity function called", 
		 IGRAPH_EINVAL);
  }

  IGRAPH_VECTOR_INIT_FINALLY(&order, no_of_nodes);
  IGRAPH_VECTOR_INIT_FINALLY(&degree, no_of_nodes);
  
  IGRAPH_CHECK(igraph_degree(graph, &degree, igraph_vss_all(), IGRAPH_ALL,
			     IGRAPH_LOOPS));
  maxdegree=igraph_vector_max(&degree)+1;
  igraph_vector_order1(&degree, &order, maxdegree);
  igraph_vector_destroy(&degree);
  IGRAPH_FINALLY_CLEAN(1);
  IGRAPH_VECTOR_INIT_FINALLY(&rank, no_of_nodes);
  for (i=0; i<no_of_nodes; i++) {
    VECTOR(rank)[ (long int)VECTOR(order)[i] ] = no_of_nodes-i-1;
  }
  
  IGRAPH_CHECK(igraph_i_adjlist_init(graph, &allneis, IGRAPH_ALL));
  IGRAPH_FINALLY(igraph_i_adjlist_destroy, &allneis);
  IGRAPH_CHECK(igraph_i_adjlist_simplify(&allneis));
  
  neis=Calloc(no_of_nodes, long int);
  if (neis==0) {
    IGRAPH_ERROR("undirected local transitivity failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, neis);
  
  IGRAPH_CHECK(igraph_vector_resize(res, no_of_nodes));
  igraph_vector_null(res);
  
  for (nn=no_of_nodes-1; nn>=0; nn--) {
    node=VECTOR(order)[nn];
    
    IGRAPH_ALLOW_INTERRUPTION();
    
    neis1=igraph_i_adjlist_get(&allneis, node);
    neilen1=igraph_vector_size(neis1);
    triples=(double)neilen1*(neilen1-1)/2;
    /* Mark the neighbors of the node */
    for (i=0; i<neilen1; i++) {
      neis[ (long int) VECTOR(*neis1)[i] ] = node+1;
    }
    
    for (i=0; i<neilen1; i++) {
      long int nei=VECTOR(*neis1)[i];
      if (VECTOR(rank)[nei] > VECTOR(rank)[node]) {
	neis2=igraph_i_adjlist_get(&allneis, nei);
	neilen2=igraph_vector_size(neis2);
	for (j=0; j<neilen2; j++) {
	  long int nei2=VECTOR(*neis2)[j];
	  if (VECTOR(rank)[nei2] < VECTOR(rank)[nei]) {
	    continue;
	  }
	  if (neis[nei2] == node+1) {
	    VECTOR(*res)[nei2] += 1;
	    VECTOR(*res)[nei] += 1;
	    VECTOR(*res)[node] += 1;
	  }
	}
      }
    }
    
    VECTOR(*res)[node] /= triples;
  }
  
  igraph_free(neis);
  igraph_i_adjlist_destroy(&allneis);
  igraph_vector_destroy(&rank);
  igraph_vector_destroy(&order);
  IGRAPH_FINALLY_CLEAN(4);
  
  return 0;
}

/**
 * \function igraph_transitivity_local_undirected
 * \brief Calculates the local transitivity (clustering coefficient)
 * of a graph
 * 
 * The transitivity measures the probability that two neighbors of a
 * vertex are connected. In case of the local transitivity, this
 * probability is calculated separately for each vertex.
 * \param graph The input graph, it can be directed but direction of
 *   the edges will be ignored.
 * \param res Pointer to an initialized vector, the result will be
 *   stored here. It will be resized as needed.
 * \param vids Vertex set, the vertices for which the local
 *   transitivity will be calculated.
 * \return Error code.
 * 
 * \sa \ref igraph_transitivity_undirected(), \ref
 * igraph_transitivity_avglocal_undirected().
 * 
 * Time complexity: O(n*d^2), n is the number of vertices for which
 * the transitivity is calculated, d is the average vertex degree.
 */

int igraph_transitivity_local_undirected(const igraph_t *graph,
					 igraph_vector_t *res,
					 const igraph_vs_t vids) {
  if (igraph_vs_is_all((igraph_vs_t*)&vids)) {
    return igraph_transitivity_local_undirected4(graph, res, vids);
  } else {
    igraph_vit_t vit;
    long int size;
    IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit);
    size=IGRAPH_VIT_SIZE(vit);
    igraph_vit_destroy(&vit);
    if (size < 100) {
      return igraph_transitivity_local_undirected1(graph, res, vids);
    } else {
      return igraph_transitivity_local_undirected2(graph, res, vids);
    }
  }
  
  return 0;
}

/**
 * \ingroup structural
 * \function igraph_transitivity_undirected
 * \brief Calculates the transitivity (clustering coefficient) of a graph.
 * 
 * </para><para>
 * The transitivity measures the probability that two neighbors of a
 * vertex are connected. More precisely this is the ratio of the
 * triangles and connected triples in the graph, the result is a
 * single real number or NaN (0/0) if there are no connected triples
 * in the graph.  Directed graphs are considered as undirected ones.
 * \param graph The graph object.  
 * \param res Pointer to a real variable, the result will be stored here.
 * \return Error code:
 *         \c IGRAPH_ENOMEM: not enough memory for
 *         temporary data. 
 *
 * \sa \ref igraph_transitivity_local_undirected(), 
 * \ref igraph_transitivity_avglocal_undirected().
 *
 * Time complexity: O(|V|*d^2), |V| is the number of vertices in 
 * the graph, d is the average node degree. 
 */


int igraph_transitivity_undirected(const igraph_t *graph,
				   igraph_real_t *res) {

  long int no_of_nodes=igraph_vcount(graph);
  igraph_real_t triples=0, triangles=0;
  long int node, nn;
  long int maxdegree;
  long int *neis;
  igraph_vector_t order;
  igraph_vector_t rank;
  igraph_vector_t degree;
  
  igraph_i_adjlist_t allneis;
  igraph_vector_t *neis1, *neis2;
  long int i, j, neilen1, neilen2;

  IGRAPH_VECTOR_INIT_FINALLY(&order, no_of_nodes);
  IGRAPH_VECTOR_INIT_FINALLY(&degree, no_of_nodes);

  IGRAPH_CHECK(igraph_degree(graph, &degree, igraph_vss_all(), IGRAPH_ALL,
			     IGRAPH_LOOPS));
  maxdegree=igraph_vector_max(&degree)+1;
  igraph_vector_order1(&degree, &order, maxdegree);
  igraph_vector_destroy(&degree);
  IGRAPH_FINALLY_CLEAN(1);
  IGRAPH_VECTOR_INIT_FINALLY(&rank, no_of_nodes);
  for (i=0; i<no_of_nodes; i++) {
    VECTOR(rank)[ (long int) VECTOR(order)[i] ]=no_of_nodes-i-1;
  }
  
  IGRAPH_CHECK(igraph_i_adjlist_init(graph, &allneis, IGRAPH_ALL));
  IGRAPH_FINALLY(igraph_i_adjlist_destroy, &allneis);
  IGRAPH_CHECK(igraph_i_adjlist_simplify(&allneis));

  neis=Calloc(no_of_nodes, long int);
  if (neis==0) {
    IGRAPH_ERROR("undirected transitivity failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, neis);
  
  for (nn=no_of_nodes-1; nn >=0; nn--) { 
    node=VECTOR(order)[nn];

    IGRAPH_ALLOW_INTERRUPTION();
    
    neis1=igraph_i_adjlist_get(&allneis, node);
    neilen1=igraph_vector_size(neis1);
    triples += (double)neilen1 * (neilen1-1);
    /* Mark the neighbors of 'node' */
    for (i=0; i<neilen1; i++) {
      long int nei=VECTOR(*neis1)[i];
      neis[nei] = node+1;
    }
    for (i=0; i<neilen1; i++) {
      long int nei=VECTOR(*neis1)[i];
      /* If 'nei' is not ready yet */      
      if (VECTOR(rank)[nei] > VECTOR(rank)[node]) {
	neis2=igraph_i_adjlist_get(&allneis, nei);
	neilen2=igraph_vector_size(neis2);
	for (j=0; j<neilen2; j++) {
	  long int nei2=VECTOR(*neis2)[j];
	  if (neis[nei2] == node+1) {
	    triangles += 1.0;
	  }
	}
      }
    }
  }
    
  Free(neis);
  igraph_i_adjlist_destroy(&allneis);
  igraph_vector_destroy(&rank);
  igraph_vector_destroy(&order);
  IGRAPH_FINALLY_CLEAN(4);
  
  *res = triangles / triples * 2.0;
  
  return 0;
}

/**
 * \ingroup structural
 * \function igraph_reciprocity
 * \brief Calculates the reciprocity of a directed graph.
 * 
 * </para><para>
 * A vertex pair (A, B) is said to be reciprocal if there are edges
 * between them in both directions. The reciprocity of a directed graph
 * is the proportion of all possible (A, B) pairs which are reciprocal,
 * provided there is at least one edge between A and B. The reciprocity
 * of an empty graph is undefined (results in an error code). Undirected
 * graphs always have a reciprocity of 1.0 unless they are empty.
 * 
 * \param graph The graph object.
 * \param res Pointer to an \c igraph_real_t which will contain the result.
 * \param ignore_loops Whether to ignore loop edges.
 * \return Error code:
 *         \c IGRAPH_EINVAL: graph has no edges
 *         \c IGRAPH_ENOMEM: not enough memory for
 *         temporary data. 
 * 
 * Time complexity: O(|V|+|E|), |V| is the number of vertices, 
 * |E| is the number of edges.
 */

int igraph_reciprocity(const igraph_t *graph, igraph_real_t *res, 
		       igraph_bool_t ignore_loops) {
  igraph_integer_t nonrec=0, rec=0;
  igraph_vector_t inneis, outneis;
  long int i;
  long int no_of_nodes=igraph_vcount(graph);
  
  /* THIS IS AN EXIT HERE !!!!!!!!!!!!!! */
  if (!igraph_is_directed(graph)) {
    *res=1.0;
    return 0;
  }

  IGRAPH_VECTOR_INIT_FINALLY(&inneis, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&outneis, 0);

  for (i=0; i<no_of_nodes; i++) {
    long int ip, op;
    igraph_neighbors(graph, &inneis, i, IGRAPH_IN);
    igraph_neighbors(graph, &outneis, i, IGRAPH_OUT);
    
    ip=op=0;
    while (ip < igraph_vector_size(&inneis) &&
	   op < igraph_vector_size(&outneis)) {
      if (VECTOR(inneis)[ip] < VECTOR(outneis)[op]) {
	nonrec += 1;
	ip++;
      } else if (VECTOR(inneis)[ip] > VECTOR(outneis)[op]) {
	nonrec += 1;
	op++;
      } else { 
	/* loop edge? */
	if (!ignore_loops || VECTOR(inneis)[ip] != i) {
	  rec += 1;
	}
	ip++;
	op++;
      }
    }
    nonrec += (igraph_vector_size(&inneis)-ip) + 
      (igraph_vector_size(&outneis)-op);
  }

  *res= rec/(rec+nonrec);
  
  igraph_vector_destroy(&inneis);
  igraph_vector_destroy(&outneis);
  IGRAPH_FINALLY_CLEAN(2);
  return 0;
}

/**
 * \function igraph_constraint
 * \brief Burt's constraint scores
 * 
 * </para><para>
 * This function calculates Burt's constraint scores for the given
 * vertices, also known as structural holes.
 * 
 * </para><para>
 * Burt's constraint is higher if ego has less, or mutually stronger
 * related (i.e. more redundant) contacts. Burt's measure of
 * constraint, C[i], of vertex i's ego network V[i], is defined for
 * directed and valued graphs,
 * <blockquote><para>
 * C[i] = sum( sum( (p[i,q] p[q,j])^2, q in V[i], q != i,j ), j in
 * V[], j != i)
 * </para></blockquote>
 * for a graph of order (ie. number od vertices) N, where proportional
 * tie strengths are defined as 
 * <blockquote><para>
 * p[i,j]=(a[i,j]+a[j,i]) / sum(a[i,k]+a[k,i], k in V[i], k != i),
 * </para></blockquote>
 * a[i,j] are elements of A and
 * the latter being the graph adjacency matrix. For isolated vertices,
 * constraint is undefined. 
 * 
 * </para><para>
 * Burt, R.S. (2004). Structural holes and good ideas. American
 * Journal of Sociology 110, 349-399.
 *
 * </para><para>
 * The first R version of this function was contributed by Jeroen
 * Bruggeman. 
 * \param graph A graph object.
 * \param res Pointer to an initialized vector, the result will be
 *        stored here. The vector will be resized to have the
 *        appropriate size for holding the result.
 * \param vids Vertex selector containing the vertices for which the
 *        constraint should be calculated.
 * \param weights Vector giving the weights of the edges. If it is
 *        \c NULL then each edge is supposed to have the same weight.
 * \return Error code.
 * 
 * Time complexity: O(|V|+E|+n*d^2), n is the number of vertices for
 * which the constraint is calculated and d is the average degree, |V|
 * is the number of vertices, |E| the number of edges in the
 * graph. If the weights argument is \c NULL then the time complexity
 * is O(|V|+n*d^2).
 */

int igraph_constraint(const igraph_t *graph, igraph_vector_t *res,
		      igraph_vs_t vids, const igraph_vector_t *weights) {
  
  long int no_of_nodes=igraph_vcount(graph);
  long int no_of_edges=igraph_ecount(graph);
  igraph_vit_t vit;
  long int nodes_to_calc;
  long int a, b, c, i, j, q;
  igraph_integer_t edge, from, to, edge2, from2, to2;
  
  igraph_vector_t contrib;
  igraph_vector_t degree;
  igraph_vector_t ineis_in, ineis_out, jneis_in, jneis_out;

  if (weights != 0 && igraph_vector_size(weights) != no_of_edges) {
    IGRAPH_ERROR("Invalid length of weight vector", IGRAPH_EINVAL);
  }

  IGRAPH_VECTOR_INIT_FINALLY(&contrib, no_of_nodes);
  IGRAPH_VECTOR_INIT_FINALLY(&degree, no_of_nodes);
  IGRAPH_VECTOR_INIT_FINALLY(&ineis_in, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&ineis_out, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&jneis_in, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&jneis_out, 0);
  
  IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
  IGRAPH_FINALLY(igraph_vit_destroy, &vit);
  nodes_to_calc=IGRAPH_VIT_SIZE(vit);

  if (weights==0) {
    IGRAPH_CHECK(igraph_degree(graph, &degree, igraph_vss_all(),
			       IGRAPH_ALL, IGRAPH_NO_LOOPS));
  } else {
    for (a=0; a<no_of_edges; a++) {
      igraph_edge(graph, a, &from, &to);
      if (from != to) {
	VECTOR(degree)[(long int) from] += VECTOR(*weights)[a];
	VECTOR(degree)[(long int) to  ] += VECTOR(*weights)[a];
      }
    }
  }
  
  IGRAPH_CHECK(igraph_vector_resize(res, nodes_to_calc));
  igraph_vector_null(res);
  
  for (a=0; a<nodes_to_calc; a++, IGRAPH_VIT_NEXT(vit)) {
    i=IGRAPH_VIT_GET(vit);
    
    /* get neighbors of i */
    IGRAPH_CHECK(igraph_adjacent(graph, &ineis_in, i, IGRAPH_IN));
    IGRAPH_CHECK(igraph_adjacent(graph, &ineis_out, i, IGRAPH_OUT));

    /* NaN for isolates */
    if (igraph_vector_size(&ineis_in) == 0 &&
	igraph_vector_size(&ineis_out) == 0) {
      VECTOR(*res)[a] = IGRAPH_NAN;
    }

    /* zero their contribution */
    for (b=0; b<igraph_vector_size(&ineis_in); b++) {
      edge=VECTOR(ineis_in)[b];
      igraph_edge(graph, edge, &from, &to);
      if (to==i) { to=from; }
      j=to;
      VECTOR(contrib)[j]=0.0;
    }
    for (b=0; b<igraph_vector_size(&ineis_out); b++) {
      edge=VECTOR(ineis_out)[b];
      igraph_edge(graph, edge, &from, &to);
      if (to==i) { to=from; }
      j=to;
      VECTOR(contrib)[j]=0.0;
    }

    /* add the direct contibutions, in-neighbors and out-neighbors */
    for (b=0; b<igraph_vector_size(&ineis_in); b++) {
      edge=VECTOR(ineis_in)[b];
      igraph_edge(graph, edge, &from, &to);
      if (to==i) { to=from; }
      j=to;
      if (i != j) {		/* excluding loops */
	if (weights) {
	  VECTOR(contrib)[j] += 
	    VECTOR(*weights)[(long int)edge]/VECTOR(degree)[i];
	} else {
	  VECTOR(contrib)[j] += 1.0/VECTOR(degree)[i];
	}
      }
    }
    if (igraph_is_directed(graph)) {
      for (b=0; b<igraph_vector_size(&ineis_out); b++) {
	edge=VECTOR(ineis_out)[b];
	igraph_edge(graph, edge, &from, &to);
	if (to==i) { to=from; }
	j=to;
	if (i != j) {
	  if (weights) {
	    VECTOR(contrib)[j] += 
	      VECTOR(*weights)[(long int)edge]/VECTOR(degree)[i];
	  } else {
	    VECTOR(contrib)[j] += 1.0/VECTOR(degree)[i];
	  }
	}
      }
    }

    /* add the indirect contributions, in-in, in-out, out-in, out-out */
    for (b=0; b<igraph_vector_size(&ineis_in); b++) {
      edge=VECTOR(ineis_in)[b];
      igraph_edge(graph, edge, &from, &to);
      if (to==i) { to=from; }
      j=to;
      if (i == j) { continue; }
      IGRAPH_CHECK(igraph_adjacent(graph, &jneis_in, j, IGRAPH_IN));
      IGRAPH_CHECK(igraph_adjacent(graph, &jneis_out, j, IGRAPH_OUT));
      for (c=0; c<igraph_vector_size(&jneis_in); c++) {
	edge2=VECTOR(jneis_in)[c];
	igraph_edge(graph, edge2, &from2, &to2);
	if (to2==j) { to2=from2; }
	q=to2;
	if (j != q) {
	  if (weights) {
	    VECTOR(contrib)[q] += 
	      VECTOR(*weights)[(long int)edge]*
	      VECTOR(*weights)[(long int)edge2]/
	      VECTOR(degree)[i]/VECTOR(degree)[j];
	  } else {
	    VECTOR(contrib)[q] += 1/VECTOR(degree)[i]/VECTOR(degree)[j];
	  }
	}
      }
      if (igraph_is_directed(graph)) {
	for (c=0; c<igraph_vector_size(&jneis_out); c++) {
	  edge2=VECTOR(jneis_out)[c];
	  igraph_edge(graph, edge2, &from2, &to2);
	  if (to2==j) { to2=from2; }
	  q=to2;
	  if (j != q) {
	    if (weights) {
	      VECTOR(contrib)[q] += 
		VECTOR(*weights)[(long int)edge]*
		VECTOR(*weights)[(long int)edge2]/
		VECTOR(degree)[i]/VECTOR(degree)[j];
	    } else {
	      VECTOR(contrib)[q] += 1/VECTOR(degree)[i]/VECTOR(degree)[j];
	    }
	  }
	}
      }
    }
    if (igraph_is_directed(graph)) {
      for (b=0; b<igraph_vector_size(&ineis_out); b++) {
	edge=VECTOR(ineis_out)[b];
	igraph_edge(graph, edge, &from, &to);
	if (to==i) { to=from; }
	j=to;
	if (i == j) { continue; }
	IGRAPH_CHECK(igraph_adjacent(graph, &jneis_in, j, IGRAPH_IN));
	IGRAPH_CHECK(igraph_adjacent(graph, &jneis_out, j, IGRAPH_OUT));
	for (c=0; c<igraph_vector_size(&jneis_in); c++) {
	  edge2=VECTOR(jneis_in)[c];
	  igraph_edge(graph, edge2, &from2, &to2);
	  if (to2==j) { to2=from2; }
	  q=to2;
	  if (j != q) {
	    if (weights) {
	      VECTOR(contrib)[q] += 
		VECTOR(*weights)[(long int)edge]*
		VECTOR(*weights)[(long int)edge2]/
		VECTOR(degree)[i]/VECTOR(degree)[j];
	    } else {
	      VECTOR(contrib)[q] += 1/VECTOR(degree)[i]/VECTOR(degree)[j];
	    }
	  }
	}
	for (c=0; c<igraph_vector_size(&jneis_out); c++) {
	  edge2=VECTOR(jneis_out)[c];
	  igraph_edge(graph, edge2, &from2, &to2);
	  if (to2==j) { to2=from2; }
	  q=to2;
	  if (j != q) {
	    if (weights) {
	      VECTOR(contrib)[q] += 
		VECTOR(*weights)[(long int)edge]*
		VECTOR(*weights)[(long int)edge2]/
		VECTOR(degree)[i]/VECTOR(degree)[j];
	    } else {
	      VECTOR(contrib)[q] += 1/VECTOR(degree)[i]/VECTOR(degree)[j];
	    }
	  }
	}
      }
    }
    
    /* squared sum of the contributions */
    for (b=0; b<igraph_vector_size(&ineis_in); b++) {
      edge=VECTOR(ineis_in)[b];
      igraph_edge(graph, edge, &from, &to);
      if (to==i) { to=from; }
      j=to;
      if (i == j) { continue; }
      VECTOR(*res)[a] += VECTOR(contrib)[j] * VECTOR(contrib)[j];
      VECTOR(contrib)[j]=0.0;
    }
    if (igraph_is_directed(graph)) {
      for (b=0; b<igraph_vector_size(&ineis_out); b++) {
	edge=VECTOR(ineis_out)[b];
	igraph_edge(graph, edge, &from, &to);
	if (to==i) { to=from; }
	j=to;
	if (i == j) { continue; }
	VECTOR(*res)[a] += VECTOR(contrib)[j] * VECTOR(contrib)[j];
	VECTOR(contrib)[j]=0.0;
      }
    }
  }

  igraph_vit_destroy(&vit);
  igraph_vector_destroy(&jneis_out);
  igraph_vector_destroy(&jneis_in);
  igraph_vector_destroy(&ineis_out);
  igraph_vector_destroy(&ineis_in);
  igraph_vector_destroy(&degree);
  igraph_vector_destroy(&contrib);
  IGRAPH_FINALLY_CLEAN(7);
  
  return 0;
}

/**
 * \function igraph_maxdegree
 * \brief Calculate the maximum degree in a graph (or set of vertices).
 * 
 * </para><para>
 * The largest in-, out- or total degree of the specified vertices is
 * calculated.
 * \param graph The input graph.
 * \param res Pointer to an integer (\c igraph_integer_t), the result
 *        will be stored here.
 * \param mode Defines the type of the degree.
 *        \c IGRAPH_OUT, out-degree,
 *        \c IGRAPH_IN, in-degree,
 *        \c IGRAPH_ALL, total degree (sum of the
 *        in- and out-degree). 
 *        This parameter is ignored for undirected graphs. 
 * \param loops Boolean, gives whether the self-loops should be
 *        counted.
 * \return Error code:
 *         \c IGRAPH_EINVVID: invalid vertex id.
 *         \c IGRAPH_EINVMODE: invalid mode argument.
 *
 * Time complexity: O(v) if
 * loops is 
 * TRUE, and
 * O(v*d)
 * otherwise. v is the number
 * vertices for which the degree will be calculated, and
 * d is their (average) degree. 
 */

int igraph_maxdegree(const igraph_t *graph, igraph_integer_t *res,
		     igraph_vs_t vids, igraph_neimode_t mode, 
		     igraph_bool_t loops) {

  igraph_vector_t tmp;
  
  IGRAPH_VECTOR_INIT_FINALLY(&tmp, 0);
  
  igraph_degree(graph, &tmp, vids, mode, loops);  
  *res=igraph_vector_max(&tmp);
  
  igraph_vector_destroy(&tmp);
  IGRAPH_FINALLY_CLEAN(1);
  return 0;
}

/**
 * \function igraph_density
 * Calculate the density of a graph.
 * 
 * </para><para>The density of a graph is simply the ratio number of
 * edges and the number of possible edges. Note that density is
 * ill-defined for graphs with multiple and/or loop edges, so consider
 * calling \ref igraph_simplify() on the graph if you know that it
 * contains multiple or loop edges. 
 * \param graph The input graph object.
 * \param res Pointer to a real number, the result will be stored
 *   here.
 * \param loops Logical constant, whether to include loops in the
 *   calculation. If this constant is TRUE then
 *   loop edges are thought to be possible in the graph (this does not
 *   neccessary means that the graph really contains any loops). If
 *   this FALSE then the result is only correct if the graph does not
 *   contain loops.
 * \return Error code.
 *
 * Time complexity: O(1).
 */

int igraph_density(const igraph_t *graph, igraph_real_t *res, 
		   igraph_bool_t loops) {

  igraph_integer_t no_of_nodes=igraph_vcount(graph);
  igraph_integer_t no_of_edges=igraph_ecount(graph);
  igraph_bool_t directed=igraph_is_directed(graph);
  
  if (!loops) {
    if (directed) {
      *res = no_of_edges / (no_of_nodes*(no_of_nodes-1));
    } else {
      *res = no_of_edges / (no_of_nodes*(no_of_nodes-1)/2);
    }
  } else {
    if (directed) {
      *res = no_of_edges / (no_of_nodes*no_of_nodes);
    } else {
      *res = no_of_edges / (no_of_nodes*no_of_nodes/2);
    }
  }
  
  return 0;
}

/**
 * \function igraph_neighborhood_size
 * \brief Calculates the size of the neighborhood of a given vertex
 * 
 * The neighborhood of a given order of a vertex includes all vertices
 * which are closer to the vertex than the order. Ie. order 0 is
 * always the vertex itself, order 1 is the vertex plus its immediate
 * neighbors, order 2 is order 1 plus the immediate neighbors of the
 * vertices in order 1, etc. 
 * 
 * </para><para> This function calculates the size of the neighborhood
 * of the given order for the given vertices. 
 * \param graph The input graph.
 * \param res Pointer to an initialized vector, the result will be
 *    stored here. It will be resized as needed.
 * \param vids The vertices for which the calculation is performed.
 * \param order Integer giving the order of the neighborhood.
 * \param mode Specifies how to use the direction of the edges if a
 *   directed graph is analyzed. For \c IGRAPH_OUT only the outgoing
 *   edges are followed, so all vertices reachable from the source
 *   vertex in at most \c order steps are counted. For \c IGRAPH_IN
 *   all vertices from which the source vertex is reachable in at most
 *   \c order steps are counted. \c IGRAPH_ALL ignores the direction
 *   of the edges. This argument is ignored for undirected graphs.
 * \return Error code.
 * 
 * \sa \ref igraph_neighborhood() for calculating the actual neighborhood, 
 * \ref igraph_neighborhood_graphs() for creating separate graphs from
 * the neighborhoods.
 * 
 * Time complexity: O(n*d*o), where n is the number vertices for which
 * the calculation is performed, d is the average degree, o is the order.
 */ 

int igraph_neighborhood_size(const igraph_t *graph, igraph_vector_t *res,
			     igraph_vs_t vids, igraph_integer_t order,
			     igraph_neimode_t mode) {

  long int no_of_nodes=igraph_vcount(graph);
  igraph_dqueue_t q;
  igraph_vit_t vit;
  long int i, j;
  long int *added;
  igraph_vector_t neis;
  
  if (order < 0) {
    IGRAPH_ERROR("Negative order in neighborhood size", IGRAPH_EINVAL);
  }

  added=Calloc(no_of_nodes, long int);
  if (added==0) {
    IGRAPH_ERROR("Cannot calculate neighborhood size", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, added);
  IGRAPH_DQUEUE_INIT_FINALLY(&q, 100);
  IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
  IGRAPH_FINALLY(igraph_vit_destroy, &vit);
  IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
  IGRAPH_CHECK(igraph_vector_resize(res, IGRAPH_VIT_SIZE(vit)));
  
  for (i=0; !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit), i++) {
    long int node=IGRAPH_VIT_GET(vit);
    long int size=1;
    added[node]=i+1;
    igraph_dqueue_clear(&q);
    if (order > 0) {
      igraph_dqueue_push(&q, node);
      igraph_dqueue_push(&q, 0);
    }
    
    while (!igraph_dqueue_empty(&q)) {
      long int actnode=igraph_dqueue_pop(&q);
      long int actdist=igraph_dqueue_pop(&q);
      long int n;
      igraph_neighbors(graph, &neis, actnode, mode);
      n=igraph_vector_size(&neis);

      if (actdist<order-1) {
	/* we add them to the q */
	for (j=0; j<n; j++) {
	  long int nei=VECTOR(neis)[j];
	  if (added[nei] != i+1) {
	    added[nei]=i+1;
	    IGRAPH_CHECK(igraph_dqueue_push(&q, nei));
	    IGRAPH_CHECK(igraph_dqueue_push(&q, actdist+1));
	    size++;
	  }
	}
      } else {
	/* we just count them, but don't add them */
	for (j=0; j<n; j++) {
	  long int nei=VECTOR(neis)[j];
	  if (added[nei] != i+1) {
	    added[nei]=i+1;
	    size++;
	  }
	}
      }

    } /* while q not empty */

    VECTOR(*res)[i]=size;
  } /* for VIT, i */

  igraph_vector_destroy(&neis);
  igraph_vit_destroy(&vit);
  igraph_dqueue_destroy(&q);
  Free(added);
  IGRAPH_FINALLY_CLEAN(4);

  return 0;
}

/** 
 * \function igraph_neighborhood
 * Calculate the neighborhood of vertices
 * 
 * The neighborhood of a given order of a vertex includes all vertices
 * which are closer to the vertex than the order. Ie. order 0 is
 * always the vertex itself, order 1 is the vertex plus its immediate
 * neighbors, order 2 is order 1 plus the immediate neighbors of the
 * vertices in order 1, etc. 
 * 
 * </para><para> This function calculates the vertices within the
 * neighborhood of the specified vertices.
 * \param graph The input graph.
 * \param res An initialized pointer vector. Note that the objects
 *    (pointers) in the vector will \em not be freed, but the pointer
 *    vector will be resized as needed. The result of the calculation
 *    will be stored here in \c vector_t objects.
 * \param vids The vertices for which the calculation is performed.
 * \param order Integer giving the order of the neighborhood.
 * \param mode Specifies how to use the direction of the edges if a
 *   directed graph is analyzed. For \c IGRAPH_OUT only the outgoing
 *   edges are followed, so all vertices reachable from the source
 *   vertex in at most \c order steps are included. For \c IGRAPH_IN
 *   all vertices from which the source vertex is reachable in at most
 *   \c order steps are included. \c IGRAPH_ALL ignores the direction
 *   of the edges. This argument is ignored for undirected graphs.
 * \return Error code.
 * 
 * \sa \ref igraph_neighborhood_size() to calculate the size of the
 * neighborhood, \ref igraph_neighborhood_graphs() for creating
 * graphs from the neighborhoods.
 * 
 * Time complexity: O(n*d*o), n is the number of vertices for which
 * the calculation is performed, d is the average degree, o is the
 * order.
 */

int igraph_neighborhood(const igraph_t *graph, igraph_vector_ptr_t *res,
			igraph_vs_t vids, igraph_integer_t order,
			igraph_neimode_t mode) {
  
  long int no_of_nodes=igraph_vcount(graph);
  igraph_dqueue_t q;
  igraph_vit_t vit;
  long int i, j;
  long int *added;
  igraph_vector_t neis;
  igraph_vector_t tmp;
  igraph_vector_t *newv;

  if (order < 0) {
    IGRAPH_ERROR("Negative order in neighborhood size", IGRAPH_EINVAL);
  }
  
  added=Calloc(no_of_nodes, long int);
  if (added==0) {
    IGRAPH_ERROR("Cannot calculate neighborhood size", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, added);
  IGRAPH_DQUEUE_INIT_FINALLY(&q, 100);
  IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
  IGRAPH_FINALLY(igraph_vit_destroy, &vit);
  IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&tmp, 0);
  IGRAPH_CHECK(igraph_vector_ptr_resize(res, IGRAPH_VIT_SIZE(vit)));
  
  for (i=0; !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit), i++) {
    long int node=IGRAPH_VIT_GET(vit);
    added[node]=i+1;
    igraph_vector_clear(&tmp);
    IGRAPH_CHECK(igraph_vector_push_back(&tmp, node));
    if (order > 0) {
      igraph_dqueue_push(&q, node);
      igraph_dqueue_push(&q, 0);
    }

    while (!igraph_dqueue_empty(&q)) {
      long int actnode=igraph_dqueue_pop(&q);
      long int actdist=igraph_dqueue_pop(&q);
      long int n;
      igraph_neighbors(graph, &neis, actnode, mode);
      n=igraph_vector_size(&neis);
      
      if (actdist<order-1) {
	/* we add them to the q */
	for (j=0; j<n; j++) {
	  long int nei=VECTOR(neis)[j];
	  if (added[nei] != i+1) {
	    added[nei]=i+1;
	    IGRAPH_CHECK(igraph_dqueue_push(&q, nei));
	    IGRAPH_CHECK(igraph_dqueue_push(&q, actdist+1));
	    IGRAPH_CHECK(igraph_vector_push_back(&tmp, nei));
	  }
	}
      } else {
	/* we just count them but don't add them to q */
	for (j=0; j<n; j++) {
	  long int nei=VECTOR(neis)[j];
	  if (added[nei] != i+1) {
	    added[nei]=i+1;
	    IGRAPH_CHECK(igraph_vector_push_back(&tmp, nei));
	  }
	}
      }

    } /* while q not empty */

    newv=Calloc(1, igraph_vector_t);
    if (newv==0) {
      IGRAPH_ERROR("Cannot calculate neighborhood", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, newv);
    IGRAPH_CHECK(igraph_vector_copy(newv, &tmp));
    VECTOR(*res)[i]=newv;
    IGRAPH_FINALLY_CLEAN(1);
  }

  igraph_vector_destroy(&tmp);
  igraph_vector_destroy(&neis);
  igraph_vit_destroy(&vit);
  igraph_dqueue_destroy(&q);
  Free(added);
  IGRAPH_FINALLY_CLEAN(5);

  return 0;
}

/**
 * \function igraph_neighborhood_graphs
 * Create graphs from the neighborhood(s) of some vertex/vertices
 * 
 * The neighborhood of a given order of a vertex includes all vertices
 * which are closer to the vertex than the order. Ie. order 0 is
 * always the vertex itself, order 1 is the vertex plus its immediate
 * neighbors, order 2 is order 1 plus the immediate neighbors of the
 * vertices in order 1, etc. 
 * 
 * </para><para> This function finds every vertex in the neighborhood
 * of a given parameter vertex and creates a graph from these
 * vertices. 
 * 
 * </para><para> The first version of this function was written by 
 * Vincent Matossian, thanks Vincent.
 * \param graph The input graph.
 * \param res Pointer to a pointer vector, the result will be stored
 *   here, ie. \c res will contain pointers to \c igraph_t
 *   objects. It will be resized if needed but note that the
 *   objects in the pointer vector will not be freed.
 * \param vids The vertices for which the calculation is performed.
 * \param order Integer giving the order of the neighborhood.
 * \param mode Specifies how to use the direction of the edges if a
 *   directed graph is analyzed. For \c IGRAPH_OUT only the outgoing
 *   edges are followed, so all vertices reachable from the source
 *   vertex in at most \c order steps are counted. For \c IGRAPH_IN
 *   all vertices from which the source vertex is reachable in at most
 *   \c order steps are counted. \c IGRAPH_ALL ignores the direction
 *   of the edges. This argument is ignored for undirected graphs.
 * \return Error code.
 * 
 * \sa \ref igraph_neighborhood_size() for calculating the neighborhood
 * sizes only, \ref igraph_neighborhood() for calculating the
 * neighborhoods (but not creating graphs).
 * 
 * Time complexity: O(n*(|V|+|E|)), where n is the number vertices for
 * which the calculation is performed, |V| and |E| are the number of
 * vertices and edges in the original input graph.
 */

int igraph_neighborhood_graphs(const igraph_t *graph, igraph_vector_ptr_t *res,
			       igraph_vs_t vids, igraph_integer_t order,
			       igraph_neimode_t mode) {
  long int no_of_nodes=igraph_vcount(graph);
  igraph_dqueue_t q;
  igraph_vit_t vit;
  long int i, j;
  long int *added;
  igraph_vector_t neis;
  igraph_vector_t tmp;
  igraph_t *newg;
  
  if (order < 0) {
    IGRAPH_ERROR("Negative order in neighborhood size", IGRAPH_EINVAL);
  }
  
  added=Calloc(no_of_nodes, long int);
  if (added==0) {
    IGRAPH_ERROR("Cannot calculate neighborhood size", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, added);
  IGRAPH_DQUEUE_INIT_FINALLY(&q, 100);
  IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
  IGRAPH_FINALLY(igraph_vit_destroy, &vit);
  IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&tmp, 0);
  IGRAPH_CHECK(igraph_vector_ptr_resize(res, IGRAPH_VIT_SIZE(vit)));
  
  for (i=0; !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit), i++) {
    long int node=IGRAPH_VIT_GET(vit);
    added[node]=i+1;
    igraph_vector_clear(&tmp);
    IGRAPH_CHECK(igraph_vector_push_back(&tmp, node));
    if (order > 0) {
      igraph_dqueue_push(&q, node);
      igraph_dqueue_push(&q, 0);
    }

    while (!igraph_dqueue_empty(&q)) {
      long int actnode=igraph_dqueue_pop(&q);
      long int actdist=igraph_dqueue_pop(&q);
      long int n;
      igraph_neighbors(graph, &neis, actnode, mode);
      n=igraph_vector_size(&neis);
      
      if (actdist<order-1) {
	/* we add them to the q */
	for (j=0; j<n; j++) {
	  long int nei=VECTOR(neis)[j];
	  if (added[nei] != i+1) {
	    added[nei]=i+1;
	    IGRAPH_CHECK(igraph_dqueue_push(&q, nei));
	    IGRAPH_CHECK(igraph_dqueue_push(&q, actdist+1));
	    IGRAPH_CHECK(igraph_vector_push_back(&tmp, nei));
	  }
	}
      } else {
	/* we just count them but don't add them to q */
	for (j=0; j<n; j++) {
	  long int nei=VECTOR(neis)[j];
	  if (added[nei] != i+1) {
	    added[nei]=i+1;
	    IGRAPH_CHECK(igraph_vector_push_back(&tmp, nei));
	  }
	}
      }

    } /* while q not empty */

    newg=Calloc(1, igraph_t);
    if (newg==0) {
      IGRAPH_ERROR("Cannot create neighborhood graph", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, newg);
    if (igraph_vector_size(&tmp) < no_of_nodes) {
      IGRAPH_CHECK(igraph_subgraph(graph, newg, igraph_vss_vector(&tmp)));
    } else {
      IGRAPH_CHECK(igraph_copy(newg, graph));
    }
    VECTOR(*res)[i]=newg;
    IGRAPH_FINALLY_CLEAN(1);
  }

  igraph_vector_destroy(&tmp);
  igraph_vector_destroy(&neis);
  igraph_vit_destroy(&vit);
  igraph_dqueue_destroy(&q);
  Free(added);
  IGRAPH_FINALLY_CLEAN(5);

  return 0;
}

/**
 * \function igraph_topological_sorting
 * Calculate a possible topological sorting of the graph
 *
 * </para><para>
 * A topological sorting of a directed acyclic graph is a linear ordering
 * of its nodes where each node comes before all nodes to which it has
 * edges. Every DAG has at least one topological sort, and may have many.
 * This function returns a possible topological sort among them. If the
 * graph is not acyclic (it has at least one cycle), a partial topological
 * sort is returned and a warning is issued.
 *
 * \param graph The input graph.
 * \param res Pointer to a vector, the result will be stored here.
 *   It will be resized if needed.
 * \param mode Specifies how to use the direction of the edges.
 *   For \c IGRAPH_OUT, the sorting order ensures that each node comes
 *   before all nodes to which it has edges, so nodes with no incoming
 *   edges go first. For \c IGRAPH_IN, it is quite the opposite: each
 *   node comes before all nodes from which it receives edges. Nodes 
 *   with no outgoing edges go first.
 * \return Error code.
 * 
 * Time complexity: O(|V|+|E|), where |V| and |E| are the number of
 * vertices and edges in the original input graph.
 */
int igraph_topological_sorting(const igraph_t* graph, igraph_vector_t *res,
			       igraph_neimode_t mode) {
  long int no_of_nodes=igraph_vcount(graph);
  igraph_vector_t degrees, neis;
  igraph_dqueue_t sources;
  igraph_neimode_t deg_mode;
  long int node, i, j;

  if (mode == IGRAPH_ALL || !igraph_is_directed(graph)) {
    IGRAPH_ERROR("topological sorting does not make sense for undirected graphs", IGRAPH_EINVAL);
  } else if (mode == IGRAPH_OUT) {
    deg_mode = IGRAPH_IN;
  } else if (mode == IGRAPH_IN) {
    deg_mode = IGRAPH_OUT;
  } else {
    IGRAPH_ERROR("invalid mode", IGRAPH_EINVAL);
  }

  IGRAPH_VECTOR_INIT_FINALLY(&degrees, no_of_nodes);
  IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
  IGRAPH_CHECK(igraph_dqueue_init(&sources, 0));
  IGRAPH_FINALLY(igraph_dqueue_destroy, &sources);
  IGRAPH_CHECK(igraph_degree(graph, &degrees, igraph_vss_all(), deg_mode, 0));

  igraph_vector_clear(res);

  /* Do we have nodes with no incoming vertices? */
  for (i=0; i<no_of_nodes; i++) {
    if (VECTOR(degrees)[i] == 0)
      IGRAPH_CHECK(igraph_dqueue_push(&sources, i));
  }

  /* Take all nodes with no incoming vertices and remove them */
  while (!igraph_dqueue_empty(&sources)) {
    node=(long)igraph_dqueue_pop(&sources);
    /* Add the node to the result vector */
    igraph_vector_push_back(res, node);
    /* Exclude the node from further source searches */
    VECTOR(degrees)[node]=-1;
    /* Get the neighbors and decrease their degrees by one */
    IGRAPH_CHECK(igraph_neighbors(graph, &neis, node, mode));
    j=igraph_vector_size(&neis);
    for (i=0; i<j; i++) {
      VECTOR(degrees)[(long)VECTOR(neis)[i]]--;
      if (VECTOR(degrees)[(long)VECTOR(neis)[i]] == 0)
	IGRAPH_CHECK(igraph_dqueue_push(&sources, VECTOR(neis)[i]));
    }
  }

  if (igraph_vector_size(res)<no_of_nodes)
    IGRAPH_WARNING("graph contains a cycle, partial result is returned");

  igraph_vector_destroy(&degrees);
  igraph_vector_destroy(&neis);
  igraph_dqueue_destroy(&sources);
  IGRAPH_FINALLY_CLEAN(3);

  return 0;
}

/**
 * \function igraph_is_simple
 * \brief Decides whether the input graph is a simple graph
 * 
 * </para><para>
 * A graph is a simple graph if it does not contain loop edges and 
 * multiple edges.
 * 
 * \param graph The input graph.
 * \param res Pointer to a boolean constant, the result
 *     is stored here.
 * \return Error code.
 * 
 * \sa \ref igraph_is_loop() and \ref igraph_is_multiple() to
 * find the loops and multiple edges and \ref igraph_simplify() to
 * get rid of them.
 * 
 * Time complexity: O(|V|+|E|). 
 */

int igraph_is_simple(const igraph_t *graph, igraph_bool_t *res) {
  long int vc=igraph_vcount(graph);
  long int ec=igraph_ecount(graph);
  
  if (vc==0 || ec==0) {
    *res=1;
  } else {
    igraph_vector_t neis;
    long int i, j, n;
    igraph_bool_t found=0;
    IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);    
    for (i=0; !found && i<vc; i++) {
      igraph_neighbors(graph, &neis, i, IGRAPH_OUT);
      n=igraph_vector_size(&neis);
      for (j=0; j<n; j++) {
	if (VECTOR(neis)[j]==i) { found=1; break; }
	if (j>0 && VECTOR(neis)[j-1]==VECTOR(neis)[j]) {
	  found=1; break;
	}
      }
    }
    *res=!found;
    igraph_vector_destroy(&neis);
    IGRAPH_FINALLY_CLEAN(1);
  }
  
  return 0;
}

/**
 * \function igraph_is_loop
 * \brief Find the loop edges in a graph
 * 
 * </para><para>
 * A loop edge is an edge from a vertex to itself.
 * \param graph The input graph.
 * \param res Pointer to an initialized boolean vector for storing the result, 
 *         it will be resized as needed.
 * \param es The edges to check, for all edges supply \ref igraph_ess_all() here.
 * \return Error code.
 * 
 * \sa \ref igraph_simplify() to get rid of loop edges.
 *
 * Time complexity: O(e), the number of edges to check.
 */

int igraph_is_loop(const igraph_t *graph, igraph_vector_bool_t *res, 
		   igraph_es_t es) {
  igraph_eit_t eit;  
  long int i;
  IGRAPH_CHECK(igraph_eit_create(graph, es, &eit));
  IGRAPH_FINALLY(igraph_eit_destroy, &eit);

  IGRAPH_CHECK(igraph_vector_bool_resize(res, IGRAPH_EIT_SIZE(eit)));

  for (i=0; !IGRAPH_EIT_END(eit); i++, IGRAPH_EIT_NEXT(eit)) {
    long int e=IGRAPH_EIT_GET(eit);
    VECTOR(*res)[i] = (IGRAPH_FROM(graph, e)==IGRAPH_TO(graph, e)) ? 1 : 0;
  }
  
  igraph_eit_destroy(&eit);
  IGRAPH_FINALLY_CLEAN(1);
  return 0;
}

/**
 * \function igraph_is_multiple
 * \brief Find the multiple edges in a graph
 * 
 * </para><para>
 * An edge is a multiple edge if there is another 
 * edge with the same head and tail vertices in the graph.
 * 
 * </para><para>
 * Note that this function returns true only for the second or more 
 * appereances of the multiple edges.
 * \param graph The input graph.
 * \param res Pointer to a boolean vector, the result will be stored 
 *        here. It will be resized as needed.
 * \param es The edges to check. Supply \ref igraph_ess_all() if you want
 *        to check all edges.
 * \return Error code.
 * 
 * \sa \ref igraph_count_multiple() and \ref igraph_simplify().
 * 
 * Time complexity: O(e*d), e is the number of edges to check and d is the 
 * average degree (out-degree in directed graphs) of the vertices at the 
 * tail of the edges.
 */

int igraph_is_multiple(const igraph_t *graph, igraph_vector_bool_t *res, 
		       igraph_es_t es) {
  igraph_eit_t eit;
  long int i;
  igraph_i_lazy_adjedgelist_t adjlist;
  
  IGRAPH_CHECK(igraph_eit_create(graph, es, &eit));
  IGRAPH_FINALLY(igraph_eit_destroy, &eit);
  IGRAPH_CHECK(igraph_i_lazy_adjedgelist_init(graph, &adjlist, IGRAPH_OUT));
  IGRAPH_FINALLY(igraph_i_lazy_adjedgelist_destroy, &adjlist);
  
  IGRAPH_CHECK(igraph_vector_bool_resize(res, IGRAPH_EIT_SIZE(eit)));
  
  for (i=0; !IGRAPH_EIT_END(eit); i++, IGRAPH_EIT_NEXT(eit)) {
    long int e=IGRAPH_EIT_GET(eit);
    long int from=IGRAPH_FROM(graph, e);
    long int to=IGRAPH_TO(graph, e);
    igraph_vector_t *neis=igraph_i_lazy_adjedgelist_get(&adjlist, from);
    long int j, n=igraph_vector_size(neis);
    VECTOR(*res)[i]=0;
    for (j=0; j<n; j++) {
      long int e2=VECTOR(*neis)[j];
      long int to2=IGRAPH_OTHER(graph,e2,from);
      if (to2==to && e2<e) {
	VECTOR(*res)[i]=1;
      }
    }
  }
  
  igraph_i_lazy_adjedgelist_destroy(&adjlist);
  igraph_eit_destroy(&eit);
  IGRAPH_FINALLY_CLEAN(2);
  return 0;
}

/**
 * \function igraph_count_multiple
 * \brief Count the number of appearance of the edges in a graph
 * 
 * </para><para>
 * If the graph has no multiple edges then the result vector will be 
 * filled with ones.
 * (An edge is a multiple edge if there is another 
 * edge with the same head and tail vertices in the graph.)
 * 
 * </para><para>
 * Note that this function returns true only for the second or more 
 * appereances of the multiple edges.
 * \param graph The input graph.
 * \param res Pointer to a vector, the result will be stored 
 *        here. It will be resized as needed.
 * \param es The edges to check. Supply \ref igraph_ess_all() if you want
 *        to check all edges.
 * \return Error code.
 * 
 * \sa \ref igraph_is_multiple() and \ref igraph_simplify().
 * 
 * Time complexity: O(e*d), e is the number of edges to check and d is the 
 * average degree (out-degree in directed graphs) of the vertices at the 
 * tail of the edges.
 */


int igraph_count_multiple(const igraph_t *graph, igraph_vector_t *res, igraph_es_t es) {
  igraph_eit_t eit;
  long int i;
  igraph_i_lazy_adjedgelist_t adjlist;
  
  IGRAPH_CHECK(igraph_eit_create(graph, es, &eit));
  IGRAPH_FINALLY(igraph_eit_destroy, &eit);
  IGRAPH_CHECK(igraph_i_lazy_adjedgelist_init(graph, &adjlist, IGRAPH_OUT));
  IGRAPH_FINALLY(igraph_i_lazy_adjedgelist_destroy, &adjlist);
  
  IGRAPH_CHECK(igraph_vector_resize(res, IGRAPH_EIT_SIZE(eit)));
  
  for (i=0; !IGRAPH_EIT_END(eit); i++, IGRAPH_EIT_NEXT(eit)) {
    long int e=IGRAPH_EIT_GET(eit);
    long int from=IGRAPH_FROM(graph, e);
    long int to=IGRAPH_TO(graph, e);
    igraph_vector_t *neis=igraph_i_lazy_adjedgelist_get(&adjlist, from);
    long int j, n=igraph_vector_size(neis);
    VECTOR(*res)[i]=1;
    for (j=0; j<n; j++) {
      long int e2=VECTOR(*neis)[j];
      long int to2=IGRAPH_OTHER(graph,e2,from);
      if (to2==to && e2!=e) {
	VECTOR(*res)[i] += 1;
      }
    }
  }
  
  igraph_i_lazy_adjedgelist_destroy(&adjlist);
  igraph_eit_destroy(&eit);
  IGRAPH_FINALLY_CLEAN(2);
  return 0;
}

/**
 * \function igraph_girth
 * \brief The girth of a graph is the length of the shortest circle in it.
 * 
 * </para><para>
 * The current implementation works for undirected graphs only, 
 * directed graphs are treated as undirected graphs. Loop edges and
 * multiple edges are ignored.
 * </para><para>
 * If the graph is a forest (ie. acyclic), then zero is returned.
 * </para><para>
 * This implementation is based on Alon Itai and Michael Rodeh:
 * Finding a minimum circuit in a graph
 * \emb Proceedings of the ninth annual ACM symposium on Theory of
 * computing \eme, 1-10, 1977. The first implementation of this
 * function was done by Keith Briggs, thanks Keith.
 * \param graph The input graph.
 * \param girth Pointer to an integer, if not \c NULL then the result
 *     will be stored here.
 * \param circle Pointer to an initialized vector, the vertex ids in
 *     the shortest circle will be stored here. If \c NULL then it is
 *     ignored. 
 * \return Error code.
 *
 * Time complexity: O((|V|+|E|)^2), |V| is the number of vertices, |E|
 * is the number of edges in the general case. If the graph has no
 * circles at all then the function needs O(|V|+|E|) time to realize
 * this and then it stops.
 */

int igraph_girth(const igraph_t *graph, igraph_integer_t *girth, 
		 igraph_vector_t *circle) {
  
  long int no_of_nodes=igraph_vcount(graph);
  igraph_dqueue_t q;
  igraph_i_lazy_adjlist_t adjlist;
  long int mincirc=LONG_MAX, minvertex=0;
  long int node;
  igraph_bool_t triangle=0;
  igraph_vector_t *neis;
  igraph_vector_long_t level;
  long int stoplevel=no_of_nodes+1;
  igraph_bool_t anycircle=0;
  long int t1=0, t2=0;
  
  IGRAPH_CHECK(igraph_i_lazy_adjlist_init(graph, &adjlist, IGRAPH_ALL, 
					  IGRAPH_I_SIMPLIFY));
  IGRAPH_FINALLY(igraph_i_lazy_adjlist_destroy, &adjlist);
  IGRAPH_DQUEUE_INIT_FINALLY(&q, 100);
  IGRAPH_CHECK(igraph_vector_long_init(&level, no_of_nodes));
  IGRAPH_FINALLY(igraph_vector_long_destroy, &level);
  
  for (node=0; !triangle && node<no_of_nodes; node++) {

    /* Are there circles in this graph at all? */
    if (node==1 && anycircle==0) { 
      igraph_bool_t conn;
      IGRAPH_CHECK(igraph_is_connected(graph, &conn, IGRAPH_WEAK));
      if (conn) {
	/* No, there are none */
	break;
      }
    }

    anycircle=0;
    igraph_dqueue_clear(&q);
    igraph_vector_long_null(&level);
    IGRAPH_CHECK(igraph_dqueue_push(&q, node));
    VECTOR(level)[node]=1;    
    
    IGRAPH_ALLOW_INTERRUPTION();
    
    while (!igraph_dqueue_empty(&q)) {
      long int actnode=igraph_dqueue_pop(&q);
      long int actlevel=VECTOR(level)[actnode];
      long int i, n;

      if (actlevel>=stoplevel) { break; }

      neis=igraph_i_lazy_adjlist_get(&adjlist, actnode);
      n=igraph_vector_size(neis);
      for (i=0; i<n; i++) {
	long int nei=VECTOR(*neis)[i];
	long int neilevel=VECTOR(level)[nei];
	if (neilevel != 0) {
	  if (neilevel==actlevel-1) {
	    continue;
	  } else {
	    /* found circle */
	    stoplevel=neilevel;
	    anycircle=1;
	    if (actlevel<mincirc) {
	      /* Is it a minimum circle? */
	      mincirc=actlevel+neilevel-1;
	      minvertex=node;
	      t1=actnode; t2=nei;
	      if (neilevel==2) {
		/* Is it a triangle? */
		triangle=1;
	      }	    
	    }
	    if (neilevel==actlevel) {
	      break;
	    }
	  }
	} else {
	  igraph_dqueue_push(&q, nei);
	  VECTOR(level)[nei]=actlevel+1;
	}
      }

    } /* while q !empty */
  } /* node */

  if (girth) {
    if (mincirc==LONG_MAX) {
      *girth=mincirc=0;
    } else {
      *girth=mincirc;
    }
  }

  /* Store the actual circle, if needed */
  if (circle) {
    IGRAPH_CHECK(igraph_vector_resize(circle, mincirc));
    if (mincirc != 0) {
      long int i, n, idx=0;
      igraph_dqueue_clear(&q);
      igraph_vector_long_null(&level); /* used for father pointers */
#define FATHER(x) (VECTOR(level)[(x)])
      IGRAPH_CHECK(igraph_dqueue_push(&q, minvertex));
      FATHER(minvertex)=minvertex;
      while (FATHER(t1)==0 || FATHER(t2)==0) {
	long int actnode=igraph_dqueue_pop(&q);
	neis=igraph_i_lazy_adjlist_get(&adjlist, actnode);
	n=igraph_vector_size(neis);
	for (i=0; i<n; i++) {
	  long int nei=VECTOR(*neis)[i];
	  if (FATHER(nei) == 0) {
	    FATHER(nei)=actnode+1;
	    igraph_dqueue_push(&q, nei);
	  }
	}
      }  /* while q !empty */
      /* Ok, now use FATHER to create the path */
      while (t1!=minvertex) {
	VECTOR(*circle)[idx++]=t1;
	t1=FATHER(t1)-1;
      }
      VECTOR(*circle)[idx]=minvertex;
      idx=mincirc-1;
      while (t2!=minvertex) {
	VECTOR(*circle)[idx--]=t2;
	t2=FATHER(t2)-1;
      }
    } /* anycircle */
  } /* circle */
#undef FATHER

  igraph_vector_long_destroy(&level);
  igraph_dqueue_destroy(&q);
  igraph_i_lazy_adjlist_destroy(&adjlist);
  IGRAPH_FINALLY_CLEAN(3);
  
  return 0;
}
    
int igraph_linegraph_undirected(const igraph_t *graph, igraph_t *linegraph) {
  
  long int no_of_edges=igraph_ecount(graph);  
  long int i, j, n;
  igraph_vector_t adjedges, adjedges2;
  igraph_vector_t edges;
  long int prev=-1;
  
  IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&adjedges, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&adjedges2, 0);
  
  for (i=0; i<no_of_edges; i++) {
    long int from=IGRAPH_FROM(graph, i);
    long int to=IGRAPH_TO(graph, i);
    
    IGRAPH_ALLOW_INTERRUPTION();

    if (from != prev) {
      IGRAPH_CHECK(igraph_adjacent(graph, &adjedges, from, IGRAPH_ALL));
    }
    n=igraph_vector_size(&adjedges);
    for (j=0; j<n; j++) {
      long int e=VECTOR(adjedges)[j];
      if (e!=i) {
	if (e<=i) { 
	  IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
	  IGRAPH_CHECK(igraph_vector_push_back(&edges, e));
	}
      }
    }
    
    IGRAPH_CHECK(igraph_adjacent(graph, &adjedges2, to, IGRAPH_ALL));
    n=igraph_vector_size(&adjedges2);
    for (j=0; j<n; j++) {
      long int e=VECTOR(adjedges2)[j];
      if (e!=i) {
	if (e<i) { 
	  IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
	  IGRAPH_CHECK(igraph_vector_push_back(&edges, e));
	}
      }
    }
    
    prev=from;
  }

  igraph_vector_destroy(&adjedges);
  igraph_vector_destroy(&adjedges2);
  IGRAPH_FINALLY_CLEAN(2);
  igraph_create(linegraph, &edges, no_of_edges, igraph_is_directed(graph));
  igraph_vector_destroy(&edges);
  IGRAPH_FINALLY_CLEAN(1);
  
  return 0;
}

int igraph_linegraph_directed(const igraph_t *graph, igraph_t *linegraph) {
  
  long int no_of_edges=igraph_ecount(graph);
  long int i, j, n;
  igraph_vector_t adjedges; 
  igraph_vector_t edges;
  long int prev=-1;
  
  IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&adjedges, 0);
  
  for (i=0; i<no_of_edges; i++) {
    long int from=IGRAPH_FROM(graph, i);

    IGRAPH_ALLOW_INTERRUPTION();
    
    if (from != prev) {
      IGRAPH_CHECK(igraph_adjacent(graph, &adjedges, from, IGRAPH_IN));
    }
    n=igraph_vector_size(&adjedges);
    for (j=0; j<n; j++) {
      long int e=VECTOR(adjedges)[j];
      IGRAPH_CHECK(igraph_vector_push_back(&edges, e));
      IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
    }
    
    prev=from;
  }
  
  igraph_vector_destroy(&adjedges);
  IGRAPH_FINALLY_CLEAN(1);
  igraph_create(linegraph, &edges, no_of_edges, igraph_is_directed(graph));
  igraph_vector_destroy(&edges);
  IGRAPH_FINALLY_CLEAN(1);
  
  return 0;
}

/**
 * \function igraph_linegraph
 * \brief Create the line graph of a graph
 *
 * The line graph L(G) of a G undirected graph is defined as follows. 
 * L(G) has one vertex for each edge in G and two vertices in L(G) are connected 
 * by an edge if their corresponding edges share an end point.
 *
 * </para><para>
 * The line graph L(G) of a G directed graph is slightly different, 
 * L(G) has one vertex for each edge in G and two vertices in L(G) are connected 
 * by a directed edge if the target of the first vertex's corresponding edge
 * is the same as the source of the second vertex's corresponding edge.
 * 
 * </para><para>
 * The first version of this function was contributed by Vincent Matossian, 
 * thanks.
 * \param graph The input graph, may be directed or undirected.
 * \param linegraph Pointer to an uninitialized graph object, the
 *        result is stored here.
 * \return Error code.
 * 
 * Time complexity: O(|V|+|E|), the number of edges plus the number of vertices.
 */

int igraph_linegraph(const igraph_t *graph, igraph_t *linegraph) {

  if (igraph_is_directed(graph)) {
    return igraph_linegraph_directed(graph, linegraph);
  } else {
    return igraph_linegraph_undirected(graph, linegraph);
  }
}

/**
 * \function igraph_add_edge
 * \brief Adds a single edge to a graph
 * 
 * </para><para>
 * For directed graphs the edge points from \p from to \p to.
 * 
 * </para><para>
 * Note that if you want to add many edges to a big graph, then it is
 * unefficient to add them one by one, it is better to collect them into
 * a vector and add all of them via a single \ref igraph_add_edges() call.
 * \param igraph The graph.
 * \param from The id of the first vertex of the edge.
 * \param to The id of the second vertex of the edge.
 * \return Error code.
 * 
 * \sa \ref igraph_add_edges() to add many edges, \ref
 * igraph_delete_edges() to remove edges and \ref
 * igraph_add_vertices() to add vertices.
 * 
 * Time complexity: O(|V|+|E|), the number of edges plus the number of
 * vertices.
 */

int igraph_add_edge(igraph_t *graph, igraph_integer_t from, igraph_integer_t to) {
  
  igraph_vector_t edges;
  int ret;
  
  IGRAPH_VECTOR_INIT_FINALLY(&edges, 2);

  VECTOR(edges)[0]=from;
  VECTOR(edges)[1]=to;
  IGRAPH_CHECK(ret=igraph_add_edges(graph, &edges, 0));
  
  igraph_vector_destroy(&edges);
  IGRAPH_FINALLY_CLEAN(1);
  return ret;
}
