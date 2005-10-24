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
#include "string.h"		/* memset & co. */

/* Internal functions */

int igraph_i_create_start(vector_t *res, vector_t *el, vector_t *index, 
			  integer_t nodes);

/** 
 * \ingroup interface
 * \brief Creates an empty graph with some vertices and no edges.
 *
 * The most basic constructor, all the other constructors should call
 * this to create a minimal graph object.
 * @param graph Pointer to a not-yet initialized graph object.
 * @param n The number of vertices in the graph, a non-negative
 *          integer number is expected.
 * @param directed Whether the graph is directed or not.
 * @return Error code.
 * 
 * Time complexity: <code>O(|V|)</code> for a graph with
 * <code>|V|</code> vertices (and no edges).
 */
int igraph_empty(igraph_t *graph, integer_t n, bool_t directed) {

  graph->n=0;
  graph->directed=directed;
  vector_init(&graph->from, 0);
  vector_init(&graph->to, 0);
  vector_init(&graph->oi, 0);
  vector_init(&graph->ii, 0);
  vector_init(&graph->os, 1);  VECTOR(graph->os)[0]=0;
  vector_init(&graph->is, 1);  VECTOR(graph->is)[0]=0;

  igraph_attribute_list_init(&graph->gal, 1);
  igraph_attribute_list_init(&graph->val, 0);

  /* add the vertices */
  igraph_add_vertices(graph, n);
  
  return 0;
}

/**
 * \ingroup interface
 * \brief Frees the memory allocated for a graph object. 
 * 
 * This function should be called for every graph object exactly once.
 *
 * This function invalidates all iterators (of course), but the
 * iterators of are graph should be destroyed before the graph itself
 * anyway. 
 * @param graph Pointer to the graph to free.
 * @return Error code.
 * 
 * Time complexity: operating system specific.
 */
int igraph_destroy(igraph_t *graph) {
  vector_destroy(&graph->from);
  vector_destroy(&graph->to);
  vector_destroy(&graph->oi);
  vector_destroy(&graph->ii);
  vector_destroy(&graph->os);
  vector_destroy(&graph->is);
  
  igraph_attribute_list_destroy(&graph->gal);
  igraph_attribute_list_destroy(&graph->val);

  return 0;
}

/**
 * \ingroup interface
 * \brief Creates an exact (deep) copy of a graph.
 * 
 * This function deeply copies a graph object to create an exact
 * replica of it. The new replica should be destroyed by calling
 * igraph_destroy() on it when not needed any more.
 * 
 * You can also create a shallow copy of a graph by simply using the
 * standard assignment operator, but be careful and do <em>not</em>
 * destroy a shallow replica. To avoid this mistake creating shallow
 * copies is not recommended.
 * @param to Pointer to an uninitialized graph object.
 * @param from Pointer to the graph object to copy.
 * @return Error code.
 *
 * Time complexity:  <code>O(|V|+|E|)</code> for a graph with
 * <code>|V|</code> vertices and <code>|E|</code> edges.
 */

int igraph_copy(igraph_t *to, igraph_t *from) {
  to->n=from->n;
  to->directed=from->directed;
  vector_copy(&to->from, &from->from);
  vector_copy(&to->to, &from->to);
  vector_copy(&to->oi, &from->oi);
  vector_copy(&to->ii, &from->ii);
  vector_copy(&to->os, &from->os);
  vector_copy(&to->is, &from->is);

  return 0;
}

/**
 * \ingroup interface
 * \brief Adds edges to a graph object. 
 * 
 * The edges are given in a vector, the
 * first two elements define the first edge (the order is
 * <code>from</code>, <code>to</code> for directed graphs). The vector
 * should contain even number of integer numbers between zero and the
 * number of vertices in the graph minus one (inclusive). If you also
 * want to add new vertices, call igraph_add_vertices() first.
 * @param graph The graph to which the edges will be added.
 * @param edges The edges themselves.
 * @return Error code.
 *
 * This function invalidates all iterators.
 *
 * Time complexity: <code>O(|V|+|E|)</code> where <code>|V|</code>
 * is the number of vertices and <code>|E|</code> is the number of
 * edges in the \em new, extended graph.
 */
int igraph_add_edges(igraph_t *graph, vector_t *edges) {
  long int no_of_edges=vector_size(&graph->from);
  long int edges_to_add=vector_size(edges)/2;
  long int i=0;
  real_t max=vector_max(edges);

  if (max > graph->n-1) {
    igraph_error("invalid vertex id\n");
  }

  /* from & to */
  vector_reserve(&graph->from, no_of_edges+edges_to_add);
  vector_reserve(&graph->to  , no_of_edges+edges_to_add);
  while (i<edges_to_add*2) {
    vector_push_back(&graph->from, VECTOR(*edges)[i++]);
    vector_push_back(&graph->to,   VECTOR(*edges)[i++]);
  }
  
  /* oi & ii */
  vector_order(&graph->from, &graph->oi, graph->n);
  vector_order(&graph->to  , &graph->ii, graph->n);
  
  /* os & is */
  igraph_i_create_start(&graph->os, &graph->from, &graph->oi, graph->n);
  igraph_i_create_start(&graph->is, &graph->to  , &graph->ii, graph->n);

  return 0;
}

/**
 * \ingroup interface
 * \brief Adds vertices to a graph. 
 *
 * This function invalidates all iterators.
 *
 * @param graph The graph object to extend.
 * @param nv Non-negative integer giving the number of 
 *           vertices to add.
 * @return Error code.
 *
 * Time complexity: <code>O(|V|)</code> where <code>|V|</code> is
 * the number of vertices in the \em new, extended graph.
 */
int igraph_add_vertices(igraph_t *graph, integer_t nv) {
  long int ec=igraph_ecount(graph);
  long int i;
  
  vector_resize(&graph->os, graph->n+nv+1);
  vector_resize(&graph->is, graph->n+nv+1);
  for (i=graph->n+1; i<graph->n+nv+1; i++) {
    VECTOR(graph->os)[i]=ec;
    VECTOR(graph->is)[i]=ec;
  }
  
  graph->n += nv;

  igraph_attribute_list_increase_length(&graph->val, nv);

  return 0;
}

/**
 * \ingroup interface
 * \brief Removes edges from a graph.
 *
 * The edges to remove are given in a vector, the first two numbers
 * define the first edge, etc. The vector should contain even number
 * of elements. If an edge to remove is not included in the graph,
 * that edge will be ignored. (Perhaps not the best approach, this
 * will likely change, as soon as there will be error handling.
 * This function cannot remove vertices, they will be kept, even if
 * they lose all their edges.
 *
 * This function invalidates all iterators.
 * @param graph The graph to work on.
 * @param edges The edges to remove.
 * @return Error code.
 *
 * Time complexity: <code>O(|V|+|E|)</code> where <code>|V|</code>
 * and <code>|E|</code> are the number of vertices and edges in the
 * \em original graph, respectively.
 */
int igraph_delete_edges(igraph_t *graph, vector_t *edges) {

  int directed=graph->directed;
  long int no_of_edges=igraph_ecount(graph);
  long int edges_to_delete=vector_size(edges)/2;
  long int really_delete=0;
  long int i;
  vector_t newfrom, newto;
  long int idx=0;

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
/*     if (d==-1) { */
/*       igraph_error("No such edge to delete"); */
/*     }        */
  }

  /* OK, all edges to delete are marked with negative numbers */
  vector_init(&newfrom, no_of_edges-really_delete);
  vector_init(&newto  , no_of_edges-really_delete);
  for (i=0; idx<no_of_edges-really_delete; i++) {
    if (VECTOR(graph->from)[i] >= 0) {
      VECTOR(newfrom)[idx]=VECTOR(graph->from)[i];
      VECTOR(newto  )[idx]=VECTOR(graph->to  )[i];
      idx++;
    }
  }
  vector_destroy(&graph->from);
  vector_destroy(&graph->to  );  
  graph->from=newfrom;
  graph->to  =newto  ;  

  /* update indices */
  vector_order(&graph->from, &graph->oi, graph->n);
  vector_order(&graph->to  , &graph->ii, graph->n);
  igraph_i_create_start(&graph->os, &graph->from, &graph->oi, graph->n);
  igraph_i_create_start(&graph->is, &graph->to  , &graph->ii, graph->n);
  
  return 0;
}

/**
 * \ingroup interface
 * \brief Removes vertices (with all their edges) from the graph.
 *
 * This function changes the ids of the vertices (except in some very
 * special cases, but these should not be relied on anyway).
 *
 * This function invalidates all iterators.
 * 
 * @param graph The graph to work on.
 * @param vertices The ids of the vertices to remove in a 
 *                 vector. The vector may contain the same id more
 *                 than once.
 * @return Error code.
 *
 * Time complexity: <code>O(|V|+|E|)</code>, <code>|V|</code> and
 * <code>|E|</code> are the number of vertices and edges in the
 * original graph.
 */
int igraph_delete_vertices(igraph_t *graph, vector_t *vertices) {

  long int no_of_edges=igraph_ecount(graph);
  long int no_of_nodes=igraph_vcount(graph);
  long int vertices_to_delete=vector_size(vertices);
  long int really_delete=0;
  long int edges_to_delete=0;
  long int i, j;
  long int *index;
  vector_t newfrom, newto;
  long int idx2=0;

  for (i=0; i<vertices_to_delete; i++) {
    long int vid=VECTOR(*vertices)[i];
    if (VECTOR(graph->os)[vid] >= 0) {
      really_delete++;      
      for (j=VECTOR(graph->os)[vid]; j<VECTOR(graph->os)[vid+1]; j++) {
	long int idx=VECTOR(graph->oi)[j];
	if (VECTOR(graph->from)[idx]>=0) {
	  edges_to_delete++;
	  VECTOR(graph->from)[idx]=-1;
	}
      }
      VECTOR(graph->os)[vid]=-1; /* mark as deleted */
      for (j=VECTOR(graph->is)[vid]; j<VECTOR(graph->is)[vid+1]; j++) {
	long int idx=VECTOR(graph->ii)[j];
	if (VECTOR(graph->from)[idx]>=0) {
	  edges_to_delete++;
	  VECTOR(graph->from)[idx]=-1;
	}
      }
    } /* if !deleted yet */
  } /* for */
  
  /* Ok, deletes vertices and edges are marked */
  /* Create index for renaming vertices */
  
  index=Calloc(no_of_nodes, long int);
  j=1;
  for (i=0; i<no_of_nodes; i++) {
    if (VECTOR(graph->os)[i]>=0) {
      index[i]=j++;
    }
  }

  /* copy & rewrite edges */
  vector_init(&newfrom, no_of_edges-edges_to_delete);
  vector_init(&newto  , no_of_edges-edges_to_delete);
  for (i=0; idx2<no_of_edges-edges_to_delete; i++) {
    if (VECTOR(graph->from)[i] >= 0) {
      VECTOR(newfrom)[idx2]=index[ (long int) VECTOR(graph->from)[i] ]-1;
      VECTOR(newto  )[idx2]=index[ (long int) VECTOR(graph->to  )[i] ]-1;
      idx2++;
    }
  }
  
  igraph_attribute_list_remove_elements(&graph->val, index, really_delete);

  Free(index);
  vector_destroy(&graph->from);
  vector_destroy(&graph->to  );  
  graph->from=newfrom;
  graph->to  =newto  ;
  
  /* update */
  graph->n -= really_delete;

  /* update indices */
  vector_order(&graph->from, &graph->oi, graph->n);
  vector_order(&graph->to  , &graph->ii, graph->n);
  igraph_i_create_start(&graph->os, &graph->from, &graph->oi, graph->n);
  igraph_i_create_start(&graph->is, &graph->to  , &graph->ii, graph->n);
  
  return 0;
}

/**
 * \ingroup interface
 * \brief The number of vertices in a graph
 * 
 * @param graph The graph.
 * @return Number of vertices.
 *
 * Time complexity: <code>O(1)</code>
 */
integer_t igraph_vcount(igraph_t *graph) {
  return graph->n;
}

/**
 * \ingroup interface
 * \brief The number of edges in a graph
 * 
 * @param graph The graph.
 * @return Number of edges.
 *
 * Time complexity: <code>O(1)</code>
 */
integer_t igraph_ecount(igraph_t *graph) {
  return vector_size(&graph->from);
}

/**
 * \ingroup interface
 * \brief Adjacent vertices to a vertex.
 *
 * @param graph The graph to work on.
 * @param neis This vector will contain the result. The vector should
 *        be initialized before and will be resized.
 * @param pnode The id of the node of which the adjacent vertices are
 *        searched. 
 * @param mode Defines the way adjacent vertices are searched for
 *        directed graphs. It can have the following values:
 *        - <b>IGRAPH_OUT</b>, vertices reachable by an edge from the specified
 *          vertex are searched,
 *        - <b>IGRAPH_IN</b>, vertices from which the specified 
 *          vertex is reachable are searched.
 *        - <b>IGRAPH_ALL</b>, both kind of vertices are searched.
 *        This parameter is ignored for undirected graphs.
 * @return Error code.
 * 
 * Time complexity: <code>O(d)</code>, <code>d</code> is the number
 * of adjacent vertices to the queried vertex.
 */
int igraph_neighbors(igraph_t *graph, vector_t *neis, integer_t pnode, 
		     igraph_neimode_t mode) {

  long int length=0, idx=0;   
  long int no_of_edges;
  long int i;

  long int node=pnode;

  no_of_edges=vector_size(&graph->from);
  if (! graph->directed) {
    mode=3;
  }

  /* Calculate needed space first & allocate it*/

  if (mode & IGRAPH_OUT) {
    length += (VECTOR(graph->os)[node+1] - VECTOR(graph->os)[node]);
  }
  if (mode & IGRAPH_IN) {
    length += (VECTOR(graph->is)[node+1] - VECTOR(graph->is)[node]);
  }
  
  vector_resize(neis, length);
  
  if (mode & IGRAPH_OUT) {
    for (i=VECTOR(graph->os)[node]; i<VECTOR(graph->os)[node+1]; i++) {
      VECTOR(*neis)[idx++] = 
	VECTOR(graph->to)[ (long int)VECTOR(graph->oi)[i] ];
    }
  }
  if (mode & IGRAPH_IN) {
    for (i=VECTOR(graph->is)[node]; i<VECTOR(graph->is)[node+1]; i++) {
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

int igraph_i_create_start(vector_t *res, vector_t *el, vector_t *index, 
			  integer_t nodes) {
  
# define EDGE(i) (VECTOR(*el)[ (long int) VECTOR(*index)[(i)] ])
  
  long int no_of_nodes;
  long int no_of_edges;
  long int i, j, idx;
  
  no_of_nodes=nodes;
  no_of_edges=vector_size(el);
  
  /* result */
  
  vector_resize(res, nodes+1);
  
  /* create the index */

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

  /* clean */

# undef EDGE
  return 0;
}

/**
 * \ingroup interface
 * \brief Is this a directed graph?
 *
 * @param graph The graph.
 * @return Logical value, <code>TRUE</code> if the graph is directed,
 * <code>FALSE</code> otherwise.
 *
 * Time complexity: <code>O(1)</code>
 */

bool_t igraph_is_directed(igraph_t *graph) {
  return graph->directed;
}

/**
 * \ingroup interface
 * \brief The degree of some vertices in a graph.
 *
 * This function calculates the in-, out- or total degree of the
 * specified vertices. 
 * @param graph The graph.
 * @param res Vector, this will contain the result. It should be
 *        initialized and will be resized to be the appropriate size.
 * @param vids Vector, giving the vertex ids of which the degree will
 *        be calculated.
 * @param mode Defines the type of the degree.
 *        - <b>IGRAPH_OUT</b>, out-degree,
 *        - <b>IGRAPH_IN</b>, in-degree,
 *        - <b>IGRAPH_ALL</b>, total degree (sum of the in- and out-degree).
 *        This parameter is ignored for undirected graphs. 
 * @param loops Boolean, gives whether the self-loops should be
 *        counted.
 *
 * Time complexity: <code>O(v)</code> if <code>loops</code> is
 * <code>TRUE</code>, and <code>O(v*d)</code>
 * otherwise. <code>v</code> is the number vertices for which the
 * degree will be calculated, and <code>d</code> is their (average)
 * degree. 
 */
int igraph_degree(igraph_t *graph, vector_t *res, vector_t *vids, 
		  igraph_neimode_t mode, bool_t loops) {

  long int nodes_to_calc;
  long int i, j;
  
  nodes_to_calc=vector_size(vids);
  if (!igraph_is_directed(graph)) {
    mode=IGRAPH_ALL;
  }

  vector_resize(res, nodes_to_calc);
  vector_null(res);

  if (loops) {
    if (mode & IGRAPH_OUT) {
      for (i=0; i<nodes_to_calc; i++) {
	long int vid=VECTOR(*vids)[i];
	VECTOR(*res)[i] += (VECTOR(graph->os)[vid+1]-VECTOR(graph->os)[vid]);
      }
    }
    if (mode & IGRAPH_IN) {
      for (i=0; i<nodes_to_calc; i++) {
	long int vid=VECTOR(*vids)[i];
	VECTOR(*res)[i] += (VECTOR(graph->is)[vid+1]-VECTOR(graph->is)[vid]);
      }
    }
  } else { /* no loops */
    if (mode & IGRAPH_OUT) {
      for (i=0; i<nodes_to_calc; i++) {
	long int vid=VECTOR(*vids)[i];
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
	long int vid=VECTOR(*vids)[i];
	VECTOR(*res)[i] += (VECTOR(graph->is)[vid+1]-VECTOR(graph->is)[vid]);
	for (j=VECTOR(graph->is)[vid]; j<VECTOR(graph->is)[vid+1]; j++) {
	  if (VECTOR(graph->from)[ (long int)VECTOR(graph->ii)[j] ]==vid) {
	    VECTOR(*res)[i] -= 1;
	  }
	}
      }
    }
  }  /* loops */

  return 0;
}

/**********************************************************
 * Testing purposes, indexed edgelist type                *
 *********************************************************/

/* #include <stdio.h> */

/* int print_vector(vector_t *v) { */
/*   long int i; */
/*   for (i=0; i<vector_size(v); i++) { */
/*     printf("%f ", VECTOR(*v)[i]); */
/*   } */
/*   printf("\n"); */
  
/*   return 0; */
/* } */

/* int print_igraph(igraph_t *graph) { */
/*   printf("Nodes: %li\n", (long int)graph->n); */
/*   printf("Directed: %i\n", (int)graph->directed); */
/*   printf("From:\n"); */
/*   print_vector(&graph->from); */
/*   printf("To:\n"); */
/*   print_vector(&graph->to); */
/*   printf("Oi:\n"); */
/*   print_vector(&graph->oi); */
/*   printf("Ii:\n"); */
/*   print_vector(&graph->ii); */
/*   printf("Os:\n"); */
/*   print_vector(&graph->os); */
/*   printf("Is:\n"); */
/*   print_vector(&graph->is); */
/*   printf("--------------------------------------\n"); */
/*   return 0; */
/* } */

/* int main() { */
  
/*   igraph_t g; */
/*   vector_t edges; */
/*   vector_t neis; */
/*   vector_t deledges; */
/*   vector_t delvert; */
/*   long int i; */
  
/*   /\* create empty graph with vertices *\/ */
/* /\*   igraph_empty(&g, 4, 1);	/\\* directed *\\/ *\/ */
/*   igraph_empty(&g, 4, 0);	/\* undirected *\/ */
/*   print_igraph(&g); */

/*   /\* add edges to it *\/ */
/*   vector_init(&edges, 8); */
/*   VECTOR(edges)[0]=0;  VECTOR(edges)[1]=1; */
/*   VECTOR(edges)[2]=1;  VECTOR(edges)[3]=2; */
/*   VECTOR(edges)[4]=2;  VECTOR(edges)[5]=2; */
/*   VECTOR(edges)[6]=2;  VECTOR(edges)[7]=3; */
/*   igraph_add_edges(&g, &edges); */
/*   print_igraph(&g); */

/*   /\* vcount, ecount *\/ */
/*   printf("Vertices: %f\n", igraph_vcount(&g)); */
/*   printf("Edges:    %f\n", igraph_ecount(&g)); */

/*   /\* neighbors *\/ */
/*   vector_init(&neis, 0); */
/*   printf("Out neighbors:\n"); */
/*   for (i=0; i<igraph_vcount(&g); i++) { */
/*     igraph_neighbors(&g, &neis, i, 1); */
/*     printf("%li: ", i); */
/*     print_vector(&neis); */
/*   } */
/*   printf("In neighbors:\n"); */
/*   for (i=0; i<igraph_vcount(&g); i++) { */
/*     igraph_neighbors(&g, &neis, i, 2); */
/*     printf("%li: ", i); */
/*     print_vector(&neis); */
/*   } */
/*   printf("All neighbors:\n"); */
/*   for (i=0; i<igraph_vcount(&g); i++) { */
/*     igraph_neighbors(&g, &neis, i, 3); */
/*     printf("%li: ", i); */
/*     print_vector(&neis); */
/*   } */

/*   /\* delete edges *\/ */
/*   vector_init(&deledges, 4); */
/*   VECTOR(deledges)[0]=0; VECTOR(deledges)[1]=1; */
/*   VECTOR(deledges)[2]=2; VECTOR(deledges)[3]=2; */
/*   igraph_delete_edges(&g, &deledges); */
/*   print_igraph(&g); */

/*   /\* delete vertices *\/ */
/*   igraph_destroy(&g); */
/*   igraph_empty(&g, 4, 1); */
/*   igraph_add_edges(&g, &edges); */
/*   print_igraph(&g); */
/*   vector_init(&delvert, 2); */
/*   VECTOR(delvert)[0]=2; */
/*   VECTOR(delvert)[1]=2; */
/*   igraph_delete_vertices(&g, &delvert); */
/*   print_igraph(&g); */

/*   /\* destroy *\/ */
/*   igraph_destroy(&g); */
/*   vector_destroy(&edges); */
/*   vector_destroy(&neis); */
/*   vector_destroy(&deledges); */
/*   vector_destroy(&delvert); */

/*   return 0; */
/* } */
