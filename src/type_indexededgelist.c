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
#include "string.h"		/* memset & co. */

/* Internal functions */

int igraph_i_create_start(igraph_vector_t *res, igraph_vector_t *el, igraph_vector_t *index, 
			  igraph_integer_t nodes);

/**
 * \section about_basic_interface
 *
 * <para>This is the very minimal API in \a igraph. All the other
 * functions use this minimal set for creating and manipulating the
 * graphs.</para>
 * 
 * <para>This is a very important principle since it makes possible to
 * implement other data representations by implementing only this
 * minimal set.</para>
 */

/** 
 * \ingroup interface
 * \function igraph_empty
 * \brief Creates an empty graph with some vertices and no edges.
 *
 * The most basic constructor, all the other constructors should call
 * this to create a minimal graph object.
 * \param graph Pointer to a not-yet initialized graph object.
 * \param n The number of vertices in the graph, a non-negative
 *          integer number is expected.
 * \param directed Whether the graph is directed or not.
 * \return Error code:
 *         \c IGRAPH_EINVAL: invalid number of vertices.
 * 
 * Time complexity: O(|V|) for a graph with
 * |V| vertices (and no edges).
 */
int igraph_empty(igraph_t *graph, igraph_integer_t n, igraph_bool_t directed) {

  if (n<0) {
    IGRAPH_ERROR("cannot create empty graph with negative number of vertices",
		  IGRAPH_EINVAL);
  }

  graph->n=0;
  graph->directed=directed;
  IGRAPH_VECTOR_INIT_FINALLY(&graph->from, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&graph->to, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&graph->oi, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&graph->ii, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&graph->os, 1);
  IGRAPH_VECTOR_INIT_FINALLY(&graph->is, 1);

  VECTOR(graph->os)[0]=0;
  VECTOR(graph->is)[0]=0;

  /* init attributes */
  graph->attr=0;
  IGRAPH_CHECK(igraph_i_attribute_init(graph));

  /* add the vertices */
  IGRAPH_CHECK(igraph_add_vertices(graph, n));
  
  IGRAPH_FINALLY_CLEAN(6);
  return 0;
}

/**
 * \ingroup interface
 * \function igraph_destroy
 * \brief Frees the memory allocated for a graph object. 
 * 
 * This function should be called for every graph object exactly once.
 *
 * This function invalidates all iterators (of course), but the
 * iterators of are graph should be destroyed before the graph itself
 * anyway. 
 * \param graph Pointer to the graph to free.
 * \return Error code.
 * 
 * Time complexity: operating system specific.
 */
int igraph_destroy(igraph_t *graph) {

  IGRAPH_I_ATTRIBUTE_DESTROY(graph);

  igraph_vector_destroy(&graph->from);
  igraph_vector_destroy(&graph->to);
  igraph_vector_destroy(&graph->oi);
  igraph_vector_destroy(&graph->ii);
  igraph_vector_destroy(&graph->os);
  igraph_vector_destroy(&graph->is);
  
  return 0;
}

/**
 * \ingroup interface
 * \function igraph_copy
 * \brief Creates an exact (deep) copy of a graph.
 * 
 * This function deeply copies a graph object to create an exact
 * replica of it. The new replica should be destroyed by calling
 * \ref igraph_destroy() on it when not needed any more.
 * 
 * You can also create a shallow copy of a graph by simply using the
 * standard assignment operator, but be careful and do \em not
 * destroy a shallow replica. To avoid this mistake creating shallow
 * copies is not recommended.
 * \param to Pointer to an uninitialized graph object.
 * \param from Pointer to the graph object to copy.
 * \return Error code.
 *
 * Time complexity:  O(|V|+|E|) for a
 * graph with |V| vertices and
 * |E| edges.
 */

int igraph_copy(igraph_t *to, const igraph_t *from) {
  to->n=from->n;
  to->directed=from->directed;
  IGRAPH_CHECK(igraph_vector_copy(&to->from, &from->from));
  IGRAPH_FINALLY(igraph_vector_destroy, &to->from);
  IGRAPH_CHECK(igraph_vector_copy(&to->to, &from->to));
  IGRAPH_FINALLY(igraph_vector_destroy, &to->to);
  IGRAPH_CHECK(igraph_vector_copy(&to->oi, &from->oi));
  IGRAPH_FINALLY(igraph_vector_destroy, &to->oi);
  IGRAPH_CHECK(igraph_vector_copy(&to->ii, &from->ii));
  IGRAPH_FINALLY(igraph_vector_destroy, &to->ii);
  IGRAPH_CHECK(igraph_vector_copy(&to->os, &from->os));
  IGRAPH_FINALLY(igraph_vector_destroy, &to->os);
  IGRAPH_CHECK(igraph_vector_copy(&to->is, &from->is));
  IGRAPH_FINALLY(igraph_vector_destroy, &to->is);

  IGRAPH_I_ATTRIBUTE_COPY(to, from); /* does IGRAPH_CHECK */

  IGRAPH_FINALLY_CLEAN(6);
  return 0;
}

/**
 * \ingroup interface
 * \function igraph_add_edges
 * \brief Adds edges to a graph object. 
 * 
 * The edges are given in a vector, the
 * first two elements define the first edge (the order is
 * <code>from</code>, <code>to</code> for directed
 * graphs). The vector 
 * should contain even number of integer numbers between zero and the
 * number of vertices in the graph minus one (inclusive). If you also
 * want to add new vertices, call igraph_add_vertices() first.
 * \param graph The graph to which the edges will be added.
 * \param edges The edges themselves.
 * \return Error code:
 *    \c IGRAPH_EINVEVECTOR: invalid (odd)
 *    edges vector length, \c IGRAPH_EINVVID:
 *    invalid vertex id in edges vector. 
 *
 * This function invalidates all iterators.
 *
 * Time complexity: O(|V|+|E|) where
 * |V| is the number of vertices and
 * |E| is the number of
 * edges in the \em new, extended graph.
 */
int igraph_add_edges(igraph_t *graph, const igraph_vector_t *edges) {
  long int no_of_edges=igraph_vector_size(&graph->from);
  long int edges_to_add=igraph_vector_size(edges)/2;
  long int i=0;
  igraph_error_handler_t *oldhandler;
  int ret1, ret2;
  igraph_vector_t newoi, newii;

  if (igraph_vector_size(edges) % 2 != 0) {
    IGRAPH_ERROR("invalid (odd) length of edges vector", IGRAPH_EINVEVECTOR);
  }
  if (!igraph_vector_isininterval(edges, 0, igraph_vcount(graph)-1)) {
    IGRAPH_ERROR("cannot add edges", IGRAPH_EINVVID);
  }

  /* from & to */
  IGRAPH_CHECK(igraph_vector_reserve(&graph->from, no_of_edges+edges_to_add));
  IGRAPH_CHECK(igraph_vector_reserve(&graph->to  , no_of_edges+edges_to_add));

  while (i<edges_to_add*2) {
    igraph_vector_push_back(&graph->from, VECTOR(*edges)[i++]); /* reserved */
    igraph_vector_push_back(&graph->to,   VECTOR(*edges)[i++]); /* reserved */
  }

  /* disable the error handler temporarily */
  oldhandler=igraph_set_error_handler(igraph_error_handler_ignore);
    
  /* oi & ii */
  ret1=igraph_vector_init(&newoi, no_of_edges);
  ret2=igraph_vector_init(&newii, no_of_edges);
  if (ret1 != 0 || ret2 != 0) {
    igraph_vector_resize(&graph->from, no_of_edges); /* gets smaller */
    igraph_vector_resize(&graph->to, no_of_edges);   /* gets smaller */
    igraph_set_error_handler(oldhandler);
    IGRAPH_ERROR("cannot add edges", IGRAPH_ERROR_SELECT_2(ret1, ret2));
  }  
  ret1=igraph_vector_order(&graph->from, &newoi, graph->n);
  ret2=igraph_vector_order(&graph->to  , &newii, graph->n);
  if (ret1 != 0 || ret2 != 0) {
    igraph_vector_resize(&graph->from, no_of_edges);
    igraph_vector_resize(&graph->to, no_of_edges);
    igraph_vector_destroy(&newoi);
    igraph_vector_destroy(&newii);
    igraph_set_error_handler(oldhandler);
    IGRAPH_ERROR("cannot add edges", IGRAPH_ERROR_SELECT_2(ret1, ret2));
  }  

  /* Attributes */
  if (graph->attr) { 
    ret1=igraph_i_attribute_add_edges(graph, edges_to_add);
    if (ret1 != 0) {
      igraph_vector_resize(&graph->from, no_of_edges);
      igraph_vector_resize(&graph->to, no_of_edges);
      igraph_vector_destroy(&newoi);
      igraph_vector_destroy(&newii);
      igraph_set_error_handler(oldhandler);
      IGRAPH_ERROR("cannot add edges", ret1);
    }  
  }
  
  /* os & is, its length does not change, error safe */
  igraph_i_create_start(&graph->os, &graph->from, &newoi, graph->n);
  igraph_i_create_start(&graph->is, &graph->to  , &newii, graph->n);

  /* everything went fine  */
  igraph_vector_destroy(&graph->oi);
  igraph_vector_destroy(&graph->ii);
  graph->oi=newoi;
  graph->ii=newii;
  igraph_set_error_handler(oldhandler);
  
  return 0;
}

/**
 * \ingroup interface
 * \function igraph_add_vertices
 * \brief Adds vertices to a graph. 
 *
 * This function invalidates all iterators.
 *
 * \param graph The graph object to extend.
 * \param nv Non-negative integer giving the number of 
 *           vertices to add.
 * \return Error code: 
 *         \c IGRAPH_EINVAL: invalid number of new
 *         vertices. 
 *
 * Time complexity: O(|V|) where
 * |V| is 
 * the number of vertices in the \em new, extended graph.
 */
int igraph_add_vertices(igraph_t *graph, igraph_integer_t nv) {
  long int ec=igraph_ecount(graph);
  long int i;
  int ret;

  if (nv < 0) {
    IGRAPH_ERROR("cannot add negative number of vertices", IGRAPH_EINVAL);
  }

  if (graph->attr) {
    IGRAPH_CHECK(igraph_i_attribute_add_vertices(graph, nv));
  }

  IGRAPH_CHECK(igraph_vector_reserve(&graph->os, graph->n+nv+1));
  IGRAPH_CHECK(igraph_vector_reserve(&graph->is, graph->n+nv+1));
  
  igraph_vector_resize(&graph->os, graph->n+nv+1); /* reserved */
  igraph_vector_resize(&graph->is, graph->n+nv+1); /* reserved */
  for (i=graph->n+1; i<graph->n+nv+1; i++) {
    VECTOR(graph->os)[i]=ec;
    VECTOR(graph->is)[i]=ec;
  }
  
  graph->n += nv;   
  
  return 0;
}

/**
 * \ingroup interface
 * \function igraph_delete_edges
 * \brief Removes edges from a graph.
 *
 * The edges to remove are given in a vector, the first two numbers
 * define the first edge, etc. The vector should contain even number
 * of elements.
 * This function cannot remove vertices, they will be kept, even if
 * they lose all their edges.
 *
 * This function invalidates all iterators.
 * \param graph The graph to work on.
 * \param edges The edges to remove.
 * \return Error code:
 *   \c IGRAPH_EINVEVECTOR: invalid (odd) length of
 *   edges vector, \c IGRAPH_EINVVID: invalid vertex
 *   id in edges vector, \c IGRAPH_EINVAL: no such
 *   edge to delete. 
 *
 * Time complexity: O(|V|+|E|) where
 * |V| 
 * and |E| are the number of vertices
 * and edges in the \em original graph, respectively.
 */
int igraph_delete_edges(igraph_t *graph, const igraph_vector_t *edges) {

  int directed=graph->directed;
  long int no_of_edges=igraph_ecount(graph);
  long int edges_to_delete=igraph_vector_size(edges)/2;
  long int really_delete=0;
  long int i, p;
  igraph_vector_t newfrom=IGRAPH_VECTOR_NULL, newto=IGRAPH_VECTOR_NULL;  
  igraph_vector_t newoi=IGRAPH_VECTOR_NULL;
  long int idx=0;
  igraph_vector_t backup=IGRAPH_VECTOR_NULL;
  int ret1, ret2;
  igraph_error_handler_t *oldhandler;

  if (igraph_vector_size(edges) % 2 != 0) {
    IGRAPH_ERROR("invalid (odd) length of edges vector", IGRAPH_EINVEVECTOR);
  }
  if (!igraph_vector_isininterval(edges, 0, igraph_vcount(graph)-1)) {
    IGRAPH_ERROR("invalid vertex id in edges vector", IGRAPH_EINVVID);
  }

  /* backup copy */
  IGRAPH_CHECK(igraph_vector_copy(&backup, &graph->from));

  /* result */

  for (i=0; i<edges_to_delete; i++) {
    long int from=VECTOR(*edges)[2*i];
    long int to  =VECTOR(*edges)[2*i+1];
    long int j=VECTOR(graph->os)[from];
    long int d=-1;
    while (d==-1 && j<VECTOR(graph->os)[from+1]) {
      long int idx=VECTOR(graph->oi)[j];
      if ( VECTOR(graph->to)[idx] == to) {
	d=idx;
      }
      j++;
    }
    if (d!=-1) {
      VECTOR(graph->from)[d]=-1;
      VECTOR(graph->to)[d]=-1;
      really_delete++;
    }
    if (! directed && d==-1) {
      j=VECTOR(graph->is)[from];
      while(d==-1 && j<VECTOR(graph->is)[from+1]) {
	long int idx=VECTOR(graph->ii)[j];
	if( VECTOR(graph->from)[idx] == to) {
	  d=idx;
	}
	j++;
      }
      if (d!=-1) {
	VECTOR(graph->from)[d]=-1;
	VECTOR(graph->to)[d]=-1;
	really_delete++;
      }
    }
    if (d==-1) {
      igraph_vector_destroy(&graph->from);
      graph->from=backup;
      IGRAPH_ERROR("No such edge to delete", IGRAPH_EINVAL);
    }
  }

  /* OK, all edges to delete are marked with negative numbers */

  oldhandler=igraph_set_error_handler(igraph_error_handler_ignore);

  ret1=igraph_vector_init(&newfrom, no_of_edges-really_delete);
  ret2=igraph_vector_init(&newto  , no_of_edges-really_delete);
  if (ret1 != 0 || ret2 != 0) {
    igraph_vector_destroy(&newfrom);
    igraph_vector_destroy(&newto);
    igraph_vector_destroy(&graph->from);
    graph->from=backup;
    igraph_set_error_handler(oldhandler);
    IGRAPH_ERROR("cannot delete edges", IGRAPH_ERROR_SELECT_2(ret1, ret2));
  }
    
  for (i=0; idx<no_of_edges-really_delete; i++) {
    if (VECTOR(graph->from)[i] >= 0) {
      VECTOR(newfrom)[idx]=VECTOR(graph->from)[i];
      VECTOR(newto  )[idx]=VECTOR(graph->to  )[i];
      idx++;
    }
  }

  /* update indices */
  ret1=igraph_vector_init(&newoi, no_of_edges-really_delete);
  if (ret1 != 0) {
    igraph_vector_destroy(&newfrom);
    igraph_vector_destroy(&newto);
    igraph_vector_destroy(&graph->from);
    graph->from=backup;
    igraph_set_error_handler(oldhandler);
    IGRAPH_ERROR("cannot delete edges", ret1);
  }
  ret1=igraph_vector_order(&newfrom, &newoi, graph->n);
  if (ret1 != 0) {
    igraph_vector_destroy(&newfrom);
    igraph_vector_destroy(&newto);
    igraph_vector_destroy(&graph->from);
    graph->from=backup;
    igraph_vector_destroy(&newoi);
    igraph_set_error_handler(oldhandler);
    IGRAPH_ERROR("cannot delete edges", ret1);
  }
  ret1=igraph_vector_order(&newto  , &graph->ii, graph->n);
  if (ret1 != 0) {
    igraph_vector_destroy(&newfrom);
    igraph_vector_destroy(&newto);
    igraph_vector_destroy(&graph->from);
    graph->from=backup;
    igraph_vector_destroy(&newoi);
    igraph_set_error_handler(oldhandler);
    IGRAPH_ERROR("cannot delete edges", ret1);
  }
  igraph_vector_destroy(&graph->oi);
  graph->oi=newoi;
  igraph_vector_destroy(&backup);

  /* attributes */
  if (graph->attr) {
    long int i, j=1;
    for (i=0; i<igraph_vector_size(&graph->from); i++) {
      if (VECTOR(graph->from)[i] > 0) {
	VECTOR(graph->from)[i]=j++;
      } else {
	VECTOR(graph->from)[i]=0;
      }
    }
    igraph_i_attribute_delete_edges(graph, &graph->from);
  }

  igraph_set_error_handler(oldhandler);

  /* These are safe */
  igraph_i_create_start(&graph->os, &newfrom, &newoi, graph->n);
  igraph_i_create_start(&graph->is, &newto  , &graph->ii, graph->n);

  igraph_vector_destroy(&graph->from);
  igraph_vector_destroy(&graph->to);  
  graph->from=newfrom;
  graph->to=newto;  

  return 0;
}

/**
 * \ingroup interface
 * \function igraph_delete_vertices
 * \brief Removes vertices (with all their edges) from the graph.
 *
 * This function changes the ids of the vertices (except in some very
 * special cases, but these should not be relied on anyway).
 *
 * This function invalidates all iterators.
 * 
 * \param graph The graph to work on.
 * \param vertices The ids of the vertices to remove in a 
 *                 vector. The vector may contain the same id more
 *                 than once.
 * \return Error code:
 *         \c IGRAPH_EINVVID: invalid vertex id.
 *
 * Time complexity: O(|V|+|E|),
 * |V| and 
 * |E| are the number of vertices and
 * edges in the original graph.
 */
int igraph_delete_vertices(igraph_t *graph, const igraph_vs_t *vertices) {

  long int no_of_edges=igraph_ecount(graph);
  long int no_of_nodes=igraph_vcount(graph);
  long int vertices_to_delete;
  long int really_delete=0;
  long int edges_to_delete=0;
  igraph_vector_t index;
  long int i, j, p;
  long int idx2=0;
  igraph_t result;
  igraph_vs_t myvertices;
  const igraph_vector_t *myverticesv;

  IGRAPH_CHECK(igraph_copy(&result, graph));
  IGRAPH_FINALLY(igraph_destroy, &result);

  IGRAPH_CHECK(igraph_vs_vectorview_it(graph, vertices, &myvertices));
  IGRAPH_FINALLY(igraph_vs_destroy, &myvertices);
  myverticesv=igraph_vs_vector_getvector(graph, &myvertices);

  vertices_to_delete=igraph_vector_size(myverticesv);

  if (!igraph_vector_isininterval(myverticesv, 0, no_of_nodes-1)) {
    IGRAPH_ERROR("invalid vertex id", IGRAPH_EINVVID);
  }
  
  /* we use result.oi and result.os to mark vertices and edges to delete */
  igraph_vector_pop_back(&result.os);	/* adjust length a bit */
  for (i=0; i<vertices_to_delete; i++) {
    long int vid=VECTOR(*myverticesv)[i];
    if (VECTOR(graph->os)[vid] >= 0) {
      really_delete++;      
      for (j=VECTOR(graph->os)[vid]; j<VECTOR(graph->os)[vid+1]; j++) {
	long int idx=VECTOR(graph->oi)[j];
	if (VECTOR(result.oi)[idx]>=0) {
	  edges_to_delete++;
	  VECTOR(result.oi)[idx]=-1;
	}
      }
      VECTOR(result.os)[vid]=-1; /* mark as deleted */
      for (j=VECTOR(graph->is)[vid]; j<VECTOR(graph->is)[vid+1]; j++) {
	long int idx=VECTOR(graph->ii)[j];
	if (VECTOR(result.oi)[idx]>=0) {
	  edges_to_delete++;
	  VECTOR(result.oi)[idx]=-1;
	}
      }
    } /* if vertex !deleted yet */
  } /* for */
  
  /* Ok, deleted vertices and edges are marked */
  /* Create index for renaming vertices */

  IGRAPH_VECTOR_INIT_FINALLY(&index, no_of_nodes);
  j=1;
  for (i=0; i<no_of_nodes; i++) {
    if (VECTOR(result.os)[i]>=0) {
      VECTOR(index)[i]=j++;
    }
  }

  /* copy & rewrite edges */
  IGRAPH_CHECK(igraph_vector_resize(&result.from, no_of_edges-edges_to_delete));
  IGRAPH_CHECK(igraph_vector_resize(&result.to,   no_of_edges-edges_to_delete));

  for (i=0; idx2<no_of_edges-edges_to_delete; i++) {
    if (VECTOR(result.oi)[i] >= 0) {
      VECTOR(result.from)[idx2]=VECTOR(index)[ (long int) VECTOR(graph->from)[i] ]-1;
      VECTOR(result.to  )[idx2]=VECTOR(index)[ (long int) VECTOR(graph->to  )[i] ]-1;
      idx2++;
    }
  }

  result.n -= really_delete;

  /* Attributes */
  if (result.attr) {
    long int i, j=1;
    for (i=0; i<igraph_vector_size(&result.oi); i++) {
      if (VECTOR(result.oi)[i] >= 0) {
	VECTOR(result.oi)[i]=j++;
      } else {
	VECTOR(result.oi)[i]=0;
      }
    }
    IGRAPH_I_ATTRIBUTE_DELETE_VERTICES(&result, &result.oi, &index);
  }
    
  /* update indices */
  IGRAPH_CHECK(igraph_vector_order(&result.from, &result.oi, result.n));
  IGRAPH_CHECK(igraph_vector_order(&result.to,   &result.ii, result.n));
  
  IGRAPH_CHECK(igraph_i_create_start(&result.os, &result.from,
				     &result.oi, result.n));
  IGRAPH_CHECK(igraph_i_create_start(&result.is, &result.to, 
				     &result.ii, result.n));

  igraph_destroy(graph);
  igraph_vs_destroy(&myvertices);
  igraph_vector_destroy(&index);
  
  *graph=result;
  IGRAPH_FINALLY_CLEAN(3);

  if (vertices->shorthand) { igraph_vs_destroy((igraph_vs_t*)vertices); }
  return 0;
}

/**
 * \ingroup interface
 * \function igraph_vcount
 * \brief The number of vertices in a graph
 * 
 * \param graph The graph.
 * \return Number of vertices.
 *
 * Time complexity: O(1)
 */
igraph_integer_t igraph_vcount(const igraph_t *graph) {
  return graph->n;
}

/**
 * \ingroup interface
 * \function igraph_ecount
 * \brief The number of edges in a graph
 * 
 * \param graph The graph.
 * \return Number of edges.
 *
 * Time complexity: O(1)
 */
igraph_integer_t igraph_ecount(const igraph_t *graph) {
  return igraph_vector_size(&graph->from);
}

/**
 * \ingroup interface
 * \function igraph_neighbors
 * \brief Adjacent vertices to a vertex.
 *
 * \param graph The graph to work on.
 * \param neis This vector will contain the result. The vector should
 *        be initialized before and will be resized.
 * \param pnode The id of the node of which the adjacent vertices are
 *        searched. 
 * \param mode Defines the way adjacent vertices are searched for
 *        directed graphs. It can have the following values:
 *        \c IGRAPH_OUT, vertices reachable by an
 *        edge from the specified vertex are searched,
 *        \c IGRAPH_IN, vertices from which the
 *        specified vertex is reachable are searched.
 *        \c IGRAPH_ALL, both kind of vertices are
 *        searched. 
 *        This parameter is ignored for undirected graphs.
 * \return Error code:
 *         \c IGRAPH_EINVVID: invalid vertex id.
 *         \c IGRAPH_EINVMODE: invalid mode argument.
 *         \c IGRAPH_ENOMEM: not enough memory.
 * 
 * Time complexity: O(d),
 * d is the number
 * of adjacent vertices to the queried vertex.
 */
int igraph_neighbors(const igraph_t *graph, igraph_vector_t *neis, igraph_integer_t pnode, 
		     igraph_neimode_t mode) {

  long int length=0, idx=0;   
  long int no_of_edges;
  long int i, j;

  long int node=pnode;

  if (node<0 || node>igraph_vcount(graph)-1) {
    IGRAPH_ERROR("cannot get neighbors", IGRAPH_EINVVID);
  }
  if (mode != IGRAPH_OUT && mode != IGRAPH_IN && 
      mode != IGRAPH_ALL) {
    IGRAPH_ERROR("cannot get neighbors", IGRAPH_EINVMODE);
  }

  no_of_edges=igraph_vector_size(&graph->from);
  if (! graph->directed) {
    mode=IGRAPH_ALL;
  }

  /* Calculate needed space first & allocate it*/

  if (mode & IGRAPH_OUT) {
    length += (VECTOR(graph->os)[node+1] - VECTOR(graph->os)[node]);
  }
  if (mode & IGRAPH_IN) {
    length += (VECTOR(graph->is)[node+1] - VECTOR(graph->is)[node]);
  }
  
  IGRAPH_CHECK(igraph_vector_resize(neis, length));
  
  if (mode & IGRAPH_OUT) {
    j=VECTOR(graph->os)[node+1];
    for (i=VECTOR(graph->os)[node]; i<j; i++) {
      VECTOR(*neis)[idx++] = 
	VECTOR(graph->to)[ (long int)VECTOR(graph->oi)[i] ];
    }
  }
  if (mode & IGRAPH_IN) {
    j=VECTOR(graph->is)[node+1];
    for (i=VECTOR(graph->is)[node]; i<j; i++) {
      VECTOR(*neis)[idx++] =
	VECTOR(graph->from)[ (long int)VECTOR(graph->ii)[i] ];
    }
  }

  return 0;
}

/**
 * \ingroup internal
 * 
 */

int igraph_i_create_start(igraph_vector_t *res, igraph_vector_t *el, igraph_vector_t *index, 
			  igraph_integer_t nodes) {
  
# define EDGE(i) (VECTOR(*el)[ (long int) VECTOR(*index)[(i)] ])
  
  long int no_of_nodes;
  long int no_of_edges;
  long int i, j, idx;
  
  no_of_nodes=nodes;
  no_of_edges=igraph_vector_size(el);
  
  /* result */
  
  IGRAPH_CHECK(igraph_vector_resize(res, nodes+1));
  
  /* create the index */

  if (igraph_vector_size(el)==0) {
    /* empty graph */
    igraph_vector_null(res);
  } else {
    idx=-1;
    for (i=0; i<=EDGE(0); i++) {
      idx++; VECTOR(*res)[idx]=0;
    }
    for (i=1; i<no_of_edges; i++) {
      long int n=EDGE(i) - EDGE((long int)VECTOR(*res)[idx]);
      for (j=0; j<n; j++) {
	idx++; VECTOR(*res)[idx]=i;
      }
    }
    j=EDGE((long int)VECTOR(*res)[idx]);
    for (i=0; i<no_of_nodes-j; i++) {
      idx++; VECTOR(*res)[idx]=no_of_edges;
    }
  }

  /* clean */

# undef EDGE
  return 0;
}

/**
 * \ingroup interface
 * \function igraph_is_directed
 * \brief Is this a directed graph?
 *
 * \param graph The graph.
 * \return Logical value, <code>TRUE</code> if the graph is directed,
 * <code>FALSE</code> otherwise.
 *
 * Time complexity: O(1)
 */

igraph_bool_t igraph_is_directed(const igraph_t *graph) {
  return graph->directed;
}

/**
 * \ingroup interface
 * \function igraph_degree
 * \brief The degree of some vertices in a graph.
 *
 * This function calculates the in-, out- or total degree of the
 * specified vertices. 
 * \param graph The graph.
 * \param res Vector, this will contain the result. It should be
 *        initialized and will be resized to be the appropriate size.
 * \param vids Vector, giving the vertex ids of which the degree will
 *        be calculated.
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
int igraph_degree(const igraph_t *graph, igraph_vector_t *res, 
		  const igraph_vs_t *vids, 
		  igraph_neimode_t mode, igraph_bool_t loops) {

  long int nodes_to_calc;
  long int i, j;
  igraph_vs_t myvids;
  const igraph_vector_t *myvidsv;

  if (vids->shorthand) { 
    IGRAPH_FINALLY(igraph_vs_destroy, (igraph_vs_t*)vids); 
  }

  IGRAPH_CHECK(igraph_vs_vectorview_it(graph, vids, &myvids));
  IGRAPH_FINALLY(igraph_vs_destroy, &myvids);
  myvidsv=igraph_vs_vector_getvector(graph, &myvids);

  if (!igraph_vector_isininterval(myvidsv, 0, igraph_vcount(graph)-1)) {
    IGRAPH_ERROR("cannot count degree", IGRAPH_EINVVID);
  }
  if (mode != IGRAPH_OUT && mode != IGRAPH_IN && mode != IGRAPH_ALL) {
    IGRAPH_ERROR("degree calculation failed", IGRAPH_EINVMODE);
  }
  
  nodes_to_calc=igraph_vector_size(myvidsv);
  if (!igraph_is_directed(graph)) {
    mode=IGRAPH_ALL;
  }

  IGRAPH_CHECK(igraph_vector_resize(res, nodes_to_calc));
  igraph_vector_null(res);

  if (loops) {
    if (mode & IGRAPH_OUT) {
      for (i=0; i<nodes_to_calc; i++) {
	long int vid=VECTOR(*myvidsv)[i];
	VECTOR(*res)[i] += (VECTOR(graph->os)[vid+1]-VECTOR(graph->os)[vid]);
      }
    }
    if (mode & IGRAPH_IN) {
      for (i=0; i<nodes_to_calc; i++) {
	long int vid=VECTOR(*myvidsv)[i];
	VECTOR(*res)[i] += (VECTOR(graph->is)[vid+1]-VECTOR(graph->is)[vid]);
      }
    }
  } else { /* no loops */
    if (mode & IGRAPH_OUT) {
      for (i=0; i<nodes_to_calc; i++) {
	long int vid=VECTOR(*myvidsv)[i];
	VECTOR(*res)[i] += (VECTOR(graph->os)[vid+1]-VECTOR(graph->os)[vid]);
	for (j=VECTOR(graph->os)[vid]; j<VECTOR(graph->os)[vid+1]; j++) {
	  if (VECTOR(graph->to)[ (long int)VECTOR(graph->oi)[j] ]==vid) {
	    VECTOR(*res)[i] -= 1;
	  }
	}
      }
    }
    if (mode & IGRAPH_IN) {
      for (i=0; i<nodes_to_calc; i++) {
	long int vid=VECTOR(*myvidsv)[i];
	VECTOR(*res)[i] += (VECTOR(graph->is)[vid+1]-VECTOR(graph->is)[vid]);
	for (j=VECTOR(graph->is)[vid]; j<VECTOR(graph->is)[vid+1]; j++) {
	  if (VECTOR(graph->from)[ (long int)VECTOR(graph->ii)[j] ]==vid) {
	    VECTOR(*res)[i] -= 1;
	  }
	}
      }
    }
  }  /* loops */

  igraph_vs_destroy(&myvids);
  IGRAPH_FINALLY_CLEAN(1);

  if (vids->shorthand) { 
    igraph_vs_destroy((igraph_vs_t*)vids); 
    IGRAPH_FINALLY_CLEAN(1);
  }
  return 0;
}

/**
 * \function igraph_edge
 * \brief Gives the head and tail vertices of an edge.
 * 
 * \param graph The graph object.
 * \param eid The edge id. 
 * \param from Pointer to an \type igraph_integer_t. The tail of the edge
 * will be placed here. 
 * \param to Pointer to an \type igraph_integer_t. The head of the edge 
 * will be placed here.
 * \return Error code. The current implementation always returns with
 * success. 
 * 
 * Added in version 0.2.</para><para>
 * 
 * Time complexity: O(1).
 */

int igraph_edge(const igraph_t *graph, igraph_integer_t eid, 
		igraph_integer_t *from, igraph_integer_t *to) {
  
  *from = VECTOR(graph->from)[(long int)eid];
  *to   = VECTOR(graph->to  )[(long int)eid];

  return 0;
}

/**
 * \function igraph_adjacent
 * \brief Gives the adjacent edges of a vertex.
 * 
 * \param graph The graph object.
 * \param eids An initialized \type vector_t object. It will be resized
 * to hold the result.
 * \param pnode A vertex id.
 * \param mode Specifies what kind of edges to include for directed
 * graphs. \c IGRAPH_OUT means only outgoing edges, \c IGRAPH_IN only
 * incoming edges, \c IGRAPH_ALL both. This parameter is ignored for
 * undirected graphs.
 * \return Error code. \c IGRAPH_EINVVID: invalid \p pnode argument, 
 *   \c IGRAPH_EINVMODE: invalid \p mode argument.
 * 
 * Added in version 0.2.</para><para>
 * 
 * Time complexity: O(d), the number of adjacent edges to \p pnode.
 */

int igraph_adjacent(const igraph_t *graph, igraph_vector_t *eids, 
		    igraph_integer_t pnode, igraph_neimode_t mode) {
  
  long int length=0, idx=0;   
  long int no_of_edges;
  long int i, j;

  long int node=pnode;

  if (node<0 || node>igraph_vcount(graph)-1) {
    IGRAPH_ERROR("cannot get neighbors", IGRAPH_EINVVID);
  }
  if (mode != IGRAPH_OUT && mode != IGRAPH_IN && 
      mode != IGRAPH_ALL) {
    IGRAPH_ERROR("cannot get neighbors", IGRAPH_EINVMODE);
  }

  no_of_edges=igraph_vector_size(&graph->from);
  if (! graph->directed) {
    mode=IGRAPH_ALL;
  }

  /* Calculate needed space first & allocate it*/

  if (mode & IGRAPH_OUT) {
    length += (VECTOR(graph->os)[node+1] - VECTOR(graph->os)[node]);
  }
  if (mode & IGRAPH_IN) {
    length += (VECTOR(graph->is)[node+1] - VECTOR(graph->is)[node]);
  }
  
  IGRAPH_CHECK(igraph_vector_resize(eids, length));
  
  if (mode & IGRAPH_OUT) {
    j=VECTOR(graph->os)[node+1];
    for (i=VECTOR(graph->os)[node]; i<j; i++) {
      VECTOR(*eids)[idx++] = VECTOR(graph->oi)[i];
    }
  }
  if (mode & IGRAPH_IN) {
    j=VECTOR(graph->is)[node+1];
    for (i=VECTOR(graph->is)[node]; i<j; i++) {
      VECTOR(*eids)[idx++] = VECTOR(graph->ii)[i];
    }
  }

  return 0;
}
