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
 * \defgroup iterators_generic Generic iterator functions
 * \ingroup iterators
 */

/**
 * \defgroup iterators_vertex Vertex iterators
 * \ingroup iterators
 */

/**
 * \defgroup iterators_edge Edge iterators
 * \ingroup iterators
 */
  
/* Vertex iterators */
#define IGRAPH_ITERATOR_VID           1
#define IGRAPH_ITERATOR_VNEIS         2

/* Edge iterators */
#define IGRAPH_ITERATOR_EID        1001
#define IGRAPH_ITERATOR_ENEIS      1002
#define IGRAPH_ITERATOR_EFROMORDER 1003

/* Constructors */

/**
 * \defgroup iterators_vid Simple vertex iterator
 * \ingroup iterators_vertex
 * \brief Iterates through the vertices of a graph in arbitrary
 * order.
 * 
 * Implemented generic operations: igraph_next(), igraph_prev(), 
 * igraph_end(), igraph_reset(), igraph_getvertex().
 * 
 * Specific operation: igraph_iterator_vid().
 * 
 * The time complexity of all generic operations is
 * <code>O(1)</code>. 
 */

/** 
 * \ingroup iterators_vid
 * \brief Initialization function for the simple vertex iterator.
 * 
 * This is the initialization function of the simple vertex iterator
 * which visits the vertices of a graph in arbitrary order. The
 * iterator <em>must</em> be destroyed by igraph_iterator_destroy() if
 * not needed any more to free the allocated memory.
 * @param graph A graph object.
 * @param it Pointer to an uninitialized iterator object.
 * @return Error code.
 * 
 * Time complexity: <code>O(1)</code>.
 */

int igraph_iterator_vid(igraph_t *graph, igraph_iterator_t *it) {  
  it->type=IGRAPH_ITERATOR_VID;
  it->data=Calloc(1, real_t);
  ((real_t*)it->data)[0]=0;
  it->next=igraph_next_vid;
  it->prev=igraph_prev_vid;
  it->end=igraph_end_vid;
  it->reset=igraph_reset_vid;
  it->getvertex=igraph_get_vertex_vid;
  it->getvertexfrom=it->getvertexto=it->getedge=0;
  it->getvertexnei=0;
  return 0;
}

/**
 * \defgroup iterators_eid Simple edge iterator.
 * \ingroup iterators_edge
 * \brief Iterates through the edges of a graph in arbitrary order.
 * 
 * Implemented generic operations: igraph_next(), igraph_prev(),
 * igraph_end(), igraph_reset(), igraph_get_vertex_from(),
 * igraph_get_vertex_to(), igraph_get_edge().
 * 
 * Specific operation: igraph_iterator_eid().
 *
 * Time complexity of all generic operations is <code>O(1)</code>.
 */

/**
 * \ingroup iterators_eid
 * \brief Initializes a simple edge iterator.
 * 
 * This iterator walks through the edges of a graph in arbitrary
 * order. Don't forget to call igraph_iterator_destroy() on it if it
 * is not needed any more.
 * @param graph Pointer to a graph.
 * @param it Pointer to an uninitialized iterator object.
 * @return Error code.
 * 
 * Time complexity: <code>O(1)</code>.
 */

int igraph_iterator_eid(igraph_t *graph, igraph_iterator_t *it) {  
  it->type=IGRAPH_ITERATOR_EID;
  it->data=Calloc(1, real_t);
  ((real_t*)it->data)[0]=0;
  it->next=igraph_next_eid;
  it->prev=igraph_prev_eid;
  it->end=igraph_end_eid;
  it->reset=igraph_reset_eid;
  it->getvertex=0;
  it->getvertexfrom=igraph_get_vertex_from_eid;
  it->getvertexto=igraph_get_vertex_to_eid;
  it->getedge=igraph_get_edge_eid;
  it->getvertexnei=0;
  return 0;
}

/**
 * \defgroup iterators_efromorder Ordered edge iterator for directed graphs
 * \ingroup iterators_edge
 * \brief Visits the edges in the order of their starting point.
 *
 * Implemented operations: igraph_next(), igraph_prev(), igraph_end(), 
 * igraph_reset(), igraph_get_vertex_from(), igraph_get_vertex_to(),
 * igraph_get_edge().
 * 
 * Specific operation: igraph_iterator_efromorder().
 * 
 * Time complexity of all generic operations: <code>O(1)</code>.
 */

/**
 * \ingroup iterators_efromorder
 * \brief Initializes a starting point ordered edge iterator.
 *
 * The order of the edges is determined by the id of the starting
 * vertex of the directed edges. Edges starting at the same vertex are
 * visited in arbitrary order.
 * @param graph Pointer to graph object.
 * @param it Pointer to an uninitialized iterator object.
 * @return Error code.
 * 
 * Time complexity: <code>O(1)</code>.
 */

int igraph_iterator_efromorder(igraph_t *graph, igraph_iterator_t *it) {
  it->type=IGRAPH_ITERATOR_EFROMORDER;
  it->data=Calloc(1, real_t);
  ((real_t*)it->data)[0]=0;
  it->next=igraph_next_efromorder;
  it->prev=igraph_prev_efromorder;
  it->end=igraph_end_efromorder;
  it->reset=igraph_reset_efromorder;
  it->getvertex=0;
  it->getvertexfrom=igraph_get_vertex_from_efromorder;
  it->getvertexto=igraph_get_vertex_to_efromorder;
  it->getedge=igraph_get_edge_efromorder;
  it->getvertexnei=0;
  return 0;
}

typedef struct igraph_iterator_eneis_data_t {
  integer_t vid;
  integer_t mode;
  integer_t oidx;
  integer_t iidx;
} igraph_iterator_eneis_data_t;

/**
 * \defgroup iterators_eneis Adjacent edges of a vertex.
 * \ingroup iterators_edge
 * \brief Iterates through the adjacenct edges of a vertex.
 * 
 * Implemented generic operations: igraph_next(), igraph_end(),
 * igraph_reset(), igraph_get_vertex_from(), igraph_get_vertex_to(),
 * igraph_get_edge(), igraph_get_vertex_nei().
 * 
 * Specific operations: igraph_iterator_eneis(),
 * igraph_iterator_eneis_set(). 
 * 
 * Time complexity of all implemented generic operation is
 * <code>O(1)</code>. 
 */

/**
 * \ingroup iterators_eneis
 * \brief Initializes an adjacent edge iterator.
 * 
 * The edges are visited in arbitrary order.
 * @param graph A graph object.
 * @param it Pointer to an uninitialized iterator object.
 * @param vid The id of the vertex of which the adjacent edges will be
 *        visited.
 * @param mode Defined which edges to visit: <b>IGRAPH_OUT</b> visits
 *        outgoing edges, <b>IGRAPH_IN</b> visits incoming edges,
 *        <b>IGRAPH_ALL</b> visits all adajcent edges. This argument
 *        is ignored for undirected graphs.
 * @return Error code.
 * 
 * Time complexity: <code>O(1)</code>.
 */

int igraph_iterator_eneis(igraph_t *graph, igraph_iterator_t *it, 
			  integer_t vid, igraph_neimode_t mode) {  

  igraph_iterator_eneis_data_t *data;
  if (!igraph_is_directed(graph)) {
    mode=IGRAPH_ALL;
  }
  
  it->type=IGRAPH_ITERATOR_ENEIS;
  it->next=igraph_next_eneis;
  it->prev=0;
  it->end=igraph_end_eneis;
  it->reset=igraph_reset_eneis;
  it->getvertex=0;
  it->getvertexfrom=igraph_get_vertex_from_eneis;
  it->getvertexto=igraph_get_vertex_to_eneis;
  it->getedge=igraph_get_edge_eneis;
  it->getvertexnei=igraph_get_vertex_nei_eneis;

  it->data=Calloc(1, igraph_iterator_eneis_data_t);
  data=it->data;
  data->vid=vid;
  data->mode=mode;
  if (mode & IGRAPH_OUT) {
    data->oidx=VECTOR(graph->os)[(long int)vid];
  } else {
    data->oidx=igraph_ecount(graph);
  }
  if (mode & IGRAPH_IN) {
    data->iidx=VECTOR(graph->is)[(long int)vid];
  } else {
    data->iidx=igraph_ecount(graph);
  }
  
  return 0;
}

typedef igraph_iterator_eneis_data_t igraph_iterator_vneis_data_t;

/**
 * \defgroup iterators_vneis Adjacent vertices of a vertex
 * \ingroup iterators_vertex
 * \brief Iterates through the adjacenct vertices of a vertex.
 * 
 * Implemented generic operations: igraph_next(), igraph_end(),
 * igraph_reset(), igraph_get_vertex().
 * 
 * Specific operation: igraph_iterator_vneis().
 * 
 * Time complexity of all implemented generic operations is
 * <code>O(1)</code>.
 */

/**
 * \ingroup iterators_vneis
 * \brief Initializes an adjacent vertex iterator.
 *
 * The adjacenct vertices are visited in arbitrary order.
 * @param graph Pointer to a graph object.
 * @param it Pointer to an uninitialized iterator object.
 * @param vid The id of the vertex of which the adjacenct vertices
 *        will be visited.
 * @param mode Defined which vertices to visit: <b>IGRAPH_OUT</b> visits
 *        vectices at outgoing edges, <b>IGRAPH_IN</b> checks incoming
 *        edges, <b>IGRAPH_ALL</b> visits all adajcent vertices. This
 *         argument is ignored for undirected graphs.
 * @return Error code.
 * 
 * Time complexity: <code>O(1)</code>.
 */

int igraph_iterator_vneis(igraph_t *graph, igraph_iterator_t *it, 
			  integer_t vid, igraph_neimode_t mode) {
  igraph_iterator_vneis_data_t *data;
  if (!igraph_is_directed(graph)) {
    mode=IGRAPH_ALL;
  }
  it->type=IGRAPH_ITERATOR_VNEIS;
  it->next=igraph_next_vneis;
  it->prev=0;
  it->end=igraph_end_vneis;
  it->reset=igraph_reset_vneis;
  it->getvertex=igraph_get_vertex_vneis;
  it->getvertexfrom=0;
  it->getvertexto=0;
  it->getedge=0;
  it->getvertexnei=0;

  it->data=Calloc(1, igraph_iterator_vneis_data_t);
  data=it->data;
  data->vid=vid;
  data->mode=mode;
  if (mode & IGRAPH_OUT) {
    data->oidx=VECTOR(graph->os)[(long int)vid];
  } else {
    data->oidx=igraph_ecount(graph);
  }
  if (mode & IGRAPH_IN) {
    data->iidx=VECTOR(graph->is)[(long int)vid];
  } else {
    data->iidx=igraph_ecount(graph);
  }
  
  return 0;  
}

/**
 * \ingroup iterators_generic
 * \brief Uninitialize an iterator and free the memory allocated for
 * it.
 *
 * All iterators must be properly destroyed by calling this function
 * on them if they are not needed any more. An iterator must be always
 * destroyed before its corresponding graph is destroyed.
 * @param graph The associated graph to the iterator.
 * @param it The iterator to destroy.
 * @return Error code.
 * 
 * Time complexity: iterator type dependent, ususally <code>O(1)</code>.
 */

int igraph_iterator_destroy(igraph_t *graph, igraph_iterator_t *it) {
  switch ((long int) it->type) {
  case IGRAPH_ITERATOR_VID:
    Free(it->data);
    break;
  case IGRAPH_ITERATOR_VNEIS:
    Free(it->data);
    break;
  case IGRAPH_ITERATOR_EID:
    free(it->data); 
    break;
  case IGRAPH_ITERATOR_ENEIS:
    Free(it->data);
    break;
  case IGRAPH_ITERATOR_EFROMORDER:
    Free(it->data);
    break;
  }
  return 0;
}

/** 
 * \ingroup iterators_generic
 * \brief Steps an iterator.
 *
 * This function is implemented by all vertex and edge iterators.
 * @param graph The associated graph of the iterator.
 * @param it The iterator object.
 * @return Error code.
 * 
 * Time complexity: iterator type dependent, ususally <code>O(1)</code>.
 */

int igraph_next(igraph_t *graph, igraph_iterator_t *it) {
  it->next(graph, it);
  return 0;
}

/** 
 * \ingroup iterators_generic
 * \brief Steps an iterator backward.
 *
 * This function is not always implemented for the iterators.
 * @param graph The associated graph of the iterator.
 * @param it The iterator object.
 * @return Error code.
 * 
 * Time complexity: iterator type dependent, ususally <code>O(1)</code>.
 */

int igraph_prev(igraph_t *graph, igraph_iterator_t *it) {
  it->prev(graph, it);
  return 0;
}

/** 
 * \ingroup iterators_generic
 * \brief Checks whether the iterator has already finished.
 *
 * @param graph The associated graph of the iterator.
 * @param it The iterator object.
 * @return Logical value, false if there are still more elements to
 * visit.
 * 
 * Time complexity: iterator type dependent, ususally <code>O(1)</code>.
 */

bool_t igraph_end(igraph_t *graph, igraph_iterator_t *it) {
  return it->end(graph, it);
}

/** 
 * \ingroup iterators_generic
 * \brief Gives the vertex at the end of the current edge.
 *
 * The exact semantics of this function depends on the type of the
 * iterator, but usually it returns the id of the ``other'' vertex of 
 * the current edge. For example if iterating through the adjacent
 * edges of a vertex, this function gives the other end point of the
 * current edge.
 * @param graph The associated graph of the iterator.
 * @param it The iterator object.
 * @return Error code.
 * 
 * Time complexity: iterator type dependent, ususally <code>O(1)</code>.
 */

integer_t igraph_get_vertex_nei(igraph_t *graph, igraph_iterator_t *it) {
  return it->getvertexnei(graph, it);
}

/** 
 * \ingroup iterators_generic
 * \brief Resets the iterator.
 *
 * The effect of this function is the same as destroying the iterator
 * and recreating it with the same parameters, but use it only if the
 * graph hasn't changed since the previous initialization.
 * @param graph The associated graph of the iterator.
 * @param it The iterator object.
 * @return Error code.
 * 
 * Time complexity: iterator type dependent, ususally <code>O(1)</code>.
 */

int igraph_reset(igraph_t *graph, igraph_iterator_t *it) {
  return it->reset(graph, it);
}

/* Vertices */

/** 
 * \ingroup iterators_generic
 * \brief Return the current vertex (mostly for vertex iterators).
 *
 * @param graph The associated graph of the iterator.
 * @param it The iterator object.
 * @return Error code.
 * 
 * Time complexity: iterator type dependent, ususally <code>O(1)</code>.
 */

integer_t igraph_get_vertex(igraph_t *graph, igraph_iterator_t *it) {
  return it->getvertex(graph, it);
}

/* Edges */

/** 
 * \ingroup iterators_generic
 * \brief Returns the first end point of the current edge (for edge
 * iterators). 
 *
 * For directed graphs this is always the <em>starting</em> end of the
 * edge. 
 * @param graph The associated graph of the iterator.
 * @param it The iterator object.
 * @return The id of the vertex.
 * 
 * Time complexity: iterator type dependent, ususally <code>O(1)</code>.
 */

integer_t igraph_get_vertex_from(igraph_t *graph, igraph_iterator_t *it) {
  return it->getvertexfrom(graph, it);
}

/** 
 * \ingroup iterators_generic
 * \brief The second end point of the current edge (for edge
 * iterators). 
 *
 * For directed graphs this is always the <em>end</em> of the edge. 
 * @param graph The associated graph of the iterator.
 * @param it The iterator object.
 * @return The id of the vertex.
 * 
 * Time complexity: iterator type dependent, ususally <code>O(1)</code>.
 */

integer_t igraph_get_vertex_to(igraph_t *graph, igraph_iterator_t *it) {
  return it->getvertexto(graph, it);
}

/** 
 * \ingroup iterators_generic
 * \brief Returns the id of the current edge (edge iterators).
 *
 * @param graph The associated graph of the iterator.
 * @param it The iterator object.
 * @return The id of the vertex.
 * 
 * Time complexity: iterator type dependent, ususally <code>O(1)</code>.
 */

integer_t igraph_get_edge(igraph_t *graph, igraph_iterator_t *it) {
  return it->getedge(graph, it);
}

/* Specifics, simple vertex iterator */
int igraph_next_vid(igraph_t *graph, igraph_iterator_t *it) {
  ((real_t*)(it->data))[0] += 1;
  return 0;
}

int igraph_prev_vid(igraph_t *graph, igraph_iterator_t *it) {
  ((real_t*)it->data)[0] -= 1;
  return 0;
}

bool_t igraph_end_vid(igraph_t *graph, igraph_iterator_t *it) {
  return ((real_t*)it->data)[0] >= igraph_vcount(graph);
}

integer_t igraph_get_vertex_vid(igraph_t *graph, igraph_iterator_t *it) {
  return ((real_t*)it->data)[0];
}

int igraph_reset_vid(igraph_t *graph, igraph_iterator_t *it) {
  ((real_t*)it->data)[0] = 0;
  return 0;
}

/* Specifics, simple edge iterator */
int igraph_next_eid(igraph_t *graph, igraph_iterator_t *it) {
  ((real_t*)it->data)[0] += 1;
  return 0;
}

int igraph_prev_eid(igraph_t *graph, igraph_iterator_t *it) {
  ((real_t*)it->data)[0] -= 1;
  return 0;
}

bool_t igraph_end_eid(igraph_t *graph, igraph_iterator_t *it) {
  return ((real_t*)it->data)[0] >= igraph_ecount(graph);
}

int igraph_reset_eid(igraph_t *graph, igraph_iterator_t *it) {
  ((real_t*)it->data)[0]=0;
  return 0;
}

integer_t igraph_get_vertex_from_eid(igraph_t *graph, igraph_iterator_t *it) {
  return VECTOR(graph->from)[ (long int) ((real_t*)it->data)[0] ];
}

integer_t igraph_get_vertex_to_eid(igraph_t *graph, igraph_iterator_t *it) {
  return VECTOR(graph->to)[ (long int) ((real_t*)it->data)[0] ];
}

integer_t igraph_get_edge_eid(igraph_t *graph, igraph_iterator_t *it) {
  return ((real_t*)it->data)[0];
}

int igraph_next_efromorder(igraph_t *graph, igraph_iterator_t *it) {
  ((real_t*)it->data)[0] += 1;
  return 0;
}

int igraph_prev_efromorder(igraph_t *graph, igraph_iterator_t *it) {
  ((real_t*)it->data)[0] -= 1;  
  return 0;
}

bool_t igraph_end_efromorder(igraph_t *graph, igraph_iterator_t *it) {
  return ((real_t*)it->data)[0] >= igraph_ecount(graph);
}

int igraph_reset_efromorder(igraph_t *graph, igraph_iterator_t *it) {
  ((real_t*)it->data)[0] = 0;
  return 0;
}

integer_t igraph_get_vertex_from_efromorder(igraph_t *graph, 
					    igraph_iterator_t *it) {
  long int idx=VECTOR(graph->oi)[ (long int) ((real_t*)it->data)[0] ];
  return VECTOR(graph->from)[ idx ];
}

integer_t igraph_get_vertex_to_efromorder(igraph_t *graph, 
					  igraph_iterator_t *it) {
  long int idx=VECTOR(graph->oi)[ (long int) ((real_t*)it->data)[0] ];
  return VECTOR(graph->to)[ idx ];
}

integer_t igraph_get_edge_efromorder(igraph_t *graph, igraph_iterator_t *it) {
  return VECTOR(graph->oi)[ (long int) ((real_t*)it->data)[0] ];
}

int igraph_next_vneis(igraph_t *graph, igraph_iterator_t *it) {
  igraph_iterator_vneis_data_t *data=it->data;
  data->oidx++;
  if (data->oidx > VECTOR(graph->os)[ (long int)data->vid+1 ]) {
    data->iidx ++;
  }
  return 0;
}

int igraph_end_vneis(igraph_t *graph, igraph_iterator_t *it) {
  igraph_iterator_vneis_data_t *data=it->data;
  return (data->oidx >= VECTOR(graph->os)[ (long int) data->vid+1 ] &&
	  data->iidx >= VECTOR(graph->is)[ (long int) data->vid+1 ]);  
}

int igraph_reset_vneis(igraph_t *graph, igraph_iterator_t *it) {
  igraph_iterator_vneis_data_t *data=it->data;
  if ((int)data->mode & IGRAPH_OUT) {
    data->oidx=VECTOR(graph->os)[(long int)data->vid];
  } else {
    data->oidx=igraph_ecount(graph);
  }
  if ((int)data->mode & IGRAPH_IN) {
    data->iidx=VECTOR(graph->is)[(long int)data->vid];
  } else {
    data->iidx=igraph_ecount(graph);
  }

  return 0;
}

integer_t igraph_get_vertex_vneis(igraph_t *graph, igraph_iterator_t *it) {
  igraph_iterator_vneis_data_t* data=it->data;
  if (data->oidx < VECTOR(graph->os)[ (long int)data->vid+1 ]) {
    long int idx=VECTOR(graph->oi)[(long int)data->oidx];
    return VECTOR(graph->to)[idx];
  } else {
    long int idx=VECTOR(graph->ii)[(long int)data->iidx];
    return VECTOR(graph->from)[idx];
  }  
}

int igraph_iterator_vneis_set(igraph_t *graph, igraph_iterator_t *it, 
			      integer_t vid, igraph_neimode_t mode) {

  igraph_iterator_vneis_data_t *data=it->data;
  if (!igraph_is_directed(graph)) {
    mode=IGRAPH_ALL;
  }
  data->vid=vid;
  data->mode=mode;
  if (mode & IGRAPH_OUT) {
    data->oidx=VECTOR(graph->os)[(long int)vid];
  } else {
    data->oidx=igraph_ecount(graph);
  }
  if (mode & IGRAPH_IN) {
    data->iidx=VECTOR(graph->is)[(long int)vid];
  } else {
    data->iidx=igraph_ecount(graph);
  }  
  
  return 0;
}

/* Iterates over the edges to and/or from a vertex */

/**
 * \ingroup iterators_eneis
 * \brief Resets the iterator but to a different vertex.
 * 
 * This function is similar to igraph_reset(), it reinitializes the
 * iterator but to a different vertex.
 * @param graph Pointer to the associated graph object.
 * @param it Pointer to the iterator object. 
 * @param vid Id of the vertex of which the adjacenct edges will be
 *        visited.
 * @param mode Type of the adjacenct edges to visit. Possible values:
 *        <b>IGRAPH_OUT</b>, <b>IGRAPH_IN</b>, <b>IGRAPH_ALL</b>.
 * @return Error code.
 */

int igraph_iterator_eneis_set(igraph_t *graph, igraph_iterator_t *it, 
			      integer_t vid, igraph_neimode_t mode) {
  igraph_iterator_eneis_data_t *data=it->data;
  if (!igraph_is_directed(graph)) {
    mode=IGRAPH_ALL;
  }
  data->vid=vid;
  data->mode=mode;
  if (mode & IGRAPH_OUT) {
    data->oidx=VECTOR(graph->os)[(long int)vid];
  } else {
    data->oidx=igraph_ecount(graph);
  }
  if (mode & IGRAPH_IN) {
    data->iidx=VECTOR(graph->is)[(long int)vid];
  } else {
    data->iidx=igraph_ecount(graph);
  }
  
  return 0;
}

int igraph_reset_eneis(igraph_t *graph, igraph_iterator_t *it) {
  igraph_iterator_eneis_data_t *data=it->data;
  if ((int)data->mode & IGRAPH_OUT) {
    data->oidx=VECTOR(graph->os)[(long int)data->vid];
  } else {
    data->oidx=igraph_ecount(graph);
  }
  if ((int)data->mode & IGRAPH_IN) {
    data->iidx=VECTOR(graph->is)[(long int)data->vid];
  } else {
    data->iidx=igraph_ecount(graph);
  }
  
  return 0;
}

int igraph_next_eneis(igraph_t *graph, igraph_iterator_t *it) {
  igraph_iterator_eneis_data_t* data=it->data;
  data->oidx ++; 
  if (data->oidx > VECTOR(graph->os)[ (long int)data->vid+1 ]) {
    data->iidx ++;
  }     
  return 0;
}

bool_t igraph_end_eneis(igraph_t *graph, igraph_iterator_t *it) {
  igraph_iterator_eneis_data_t* data=it->data;
  return (data->oidx >= VECTOR(graph->os)[ (long int) data->vid+1 ] &&
	  data->iidx >= VECTOR(graph->is)[ (long int) data->vid+1 ]);
}

integer_t igraph_get_vertex_from_eneis(igraph_t *graph, 
				       igraph_iterator_t *it) {
  igraph_iterator_eneis_data_t* data=it->data;
  if (data->oidx < VECTOR(graph->os)[ (long int)data->vid+1 ]) {
    long int idx=VECTOR(graph->oi)[(long int)data->oidx];
    return VECTOR(graph->from)[idx];
  } else {
    long int idx=VECTOR(graph->ii)[(long int)data->iidx];
    return VECTOR(graph->from)[idx];
  }  
}

integer_t igraph_get_vertex_to_eneis(igraph_t *graph, igraph_iterator_t *it) {
  igraph_iterator_eneis_data_t* data=it->data;
  if (data->oidx < VECTOR(graph->os)[ (long int)data->vid+1 ]) {
    long int idx=VECTOR(graph->oi)[(long int)data->oidx];
    return VECTOR(graph->to)[idx];
  } else {
    long int idx=VECTOR(graph->ii)[(long int)data->iidx];
    return VECTOR(graph->to)[idx];
  }  
}

integer_t igraph_get_vertex_nei_eneis(igraph_t *graph, igraph_iterator_t *it) {
  igraph_iterator_eneis_data_t* data=it->data;
  if (data->oidx < VECTOR(graph->os)[ (long int)data->vid+1 ]) {
    long int idx=VECTOR(graph->oi)[(long int)data->oidx];
    return VECTOR(graph->to)[idx];
  } else {
    long int idx=VECTOR(graph->ii)[(long int)data->iidx];
    return VECTOR(graph->from)[idx];
  }  
}  

integer_t igraph_get_edge_eneis(igraph_t *graph, igraph_iterator_t *it) {
  igraph_iterator_eneis_data_t* data=it->data;
  if (data->oidx < VECTOR(graph->os)[ (long int)data->vid+1 ]) {
    return VECTOR(graph->oi)[(long int)data->oidx];
  } else {
    return VECTOR(graph->ii)[(long int)data->iidx];
  }
}

/**********************************************************
 * Testing purposes, indexed edgelist type                *
 *********************************************************/

/* int main() { */
  
/*   igraph_t g; */
/*   vector_t edges; */
/*   igraph_iterator_t it; */

/*   vector_init(&edges, 10); */
/*   VECTOR(edges)[0]=0;  VECTOR(edges)[1]=1; */
/*   VECTOR(edges)[2]=0;  VECTOR(edges)[3]=2; */
/*   VECTOR(edges)[4]=3;  VECTOR(edges)[5]=0; */
/*   VECTOR(edges)[6]=0;  VECTOR(edges)[7]=4; */
/*   VECTOR(edges)[8]=3;  VECTOR(edges)[9]=4; */
/*   igraph_create(&g, &edges, 0, 1); */
  
/*   print_igraph(&g); */

/*   igraph_iterator_eneis(&g, &it, 3, 3); */
/*   while (! igraph_end(&g, &it)) { */
/*     printf("%f: %f -> %f\n", igraph_get_edge(&g, &it),  */
/* 	   igraph_get_vertex_from(&g, &it), */
/* 	   igraph_get_vertex_to(&g, &it)); */
/*     igraph_next(&g, &it); */
/*   } */

/*   printf("----------------------\n"); */
/*   igraph_iterator_eneis_set(&g, &it, 0, 3); */
/*   while (! igraph_end(&g, &it)) { */
/*     printf("%f: %f -> %f\n", igraph_get_edge(&g, &it),  */
/* 	   igraph_get_vertex_from(&g, &it), */
/* 	   igraph_get_vertex_to(&g, &it)); */
/*     igraph_next(&g, &it); */
/*   }   */
  
/*   igraph_iterator_destroy(&g, &it); */
/*   igraph_destroy(&g); */
/*   vector_destroy(&edges); */
  
/*   return 0; */
/* } */
	   
