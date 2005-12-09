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
#include "random.h"

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
 * @return Error code:
 *         - <b>IGRAPH_ENOMEM</b>, not enough memory.
 * 
 * Time complexity: <code>O(1)</code>.
 */

int igraph_iterator_vid(igraph_t *graph, igraph_iterator_t *it) {  
  it->type=IGRAPH_ITERATOR_VID;
  it->data=Calloc(1, real_t);
  if (it->data == 0) {
    IGRAPH_FERROR("cannot create iterator", IGRAPH_ENOMEM);
  }
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
 * @return Error code:
 *         - <b>IGRAPH_ENOMEM</b>: not enough memory.
 * 
 * Time complexity: <code>O(1)</code>.
 */

int igraph_iterator_eid(igraph_t *graph, igraph_iterator_t *it) {  
  it->type=IGRAPH_ITERATOR_EID;
  it->data=Calloc(1, real_t);
  if (it->data == 0) {
    IGRAPH_FERROR("cannot create iterator", IGRAPH_ENOMEM);
  }
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
 * @return Error code:
 *         - <b>IGRAPH_ENOMEM</b>, not enough memory.
 * 
 * Time complexity: <code>O(1)</code>.
 */

int igraph_iterator_efromorder(igraph_t *graph, igraph_iterator_t *it) {
  it->type=IGRAPH_ITERATOR_EFROMORDER;
  it->data=Calloc(1, real_t);
  if (it->data == 0) {
    IGRAPH_FERROR("cannot create iterator", IGRAPH_ENOMEM);
  }
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
 * @return Error code:
 *         - <b>IGRAPH_ENOMEM</b>, not enough memory.
 * 
 * Time complexity: <code>O(1)</code>.
 */

int igraph_iterator_eneis(igraph_t *graph, igraph_iterator_t *it, 
			  integer_t vid, igraph_neimode_t mode) {  

  real_t *data;
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

  it->data=Calloc(4, real_t);
  if (it->data == 0) {
    IGRAPH_FERROR("cannot create iterator", IGRAPH_ENOMEM);
  }
  data=it->data;
  data[0]=vid;
  data[1]=mode;
  if (mode & IGRAPH_OUT) {
    data[2]=VECTOR(graph->os)[(long int)vid];
  } else {
    data[2]=igraph_ecount(graph);
  }
  if (mode & IGRAPH_IN) {
    data[3]=VECTOR(graph->is)[(long int)vid];
  } else {
    data[3]=igraph_ecount(graph);
  }
  
  return 0;
}

/**
 * \defgroup iterators_vneis Adjacent vertices of a vertex
 * \ingroup iterators_vertex
 * \brief Iterates through the adjacenct vertices of a vertex.
 * 
 * Implemented generic operations: igraph_next(), igraph_end(),
 * igraph_reset(), igraph_get_vertex().
 * 
 * Specific operation: igraph_iterator_vneis(),
 * igraph_iterator_vneis_set(). 
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
 * @return Error code:
 *         - <b>IGRAPH_ENOMEM</b>, not enough memory.
 * 
 * Time complexity: <code>O(1)</code>.
 */

int igraph_iterator_vneis(igraph_t *graph, igraph_iterator_t *it, 
			  integer_t vid, igraph_neimode_t mode) {
  real_t *data;
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

  it->data=Calloc(4, real_t);
  if (it->data == 0) {
    IGRAPH_FERROR("cannot create iterator", IGRAPH_ENOMEM);
  }
  data=it->data;
  data[0]=vid;
  data[1]=mode;
  if (mode & IGRAPH_OUT) {
    data[2]=VECTOR(graph->os)[(long int)vid];
  } else {
    data[2]=igraph_ecount(graph);
  }
  if (mode & IGRAPH_IN) {
    data[3]=VECTOR(graph->is)[(long int)vid];
  } else {
    data[3]=igraph_ecount(graph);
  }
  
  return 0;  
}

/**
 * \defgroup iterators_randomwalk Random walker iterator.
 * \ingroup iterators_vertex
 * \brief Performs a random walk from a given staring vertex.
 * 
 * Implemented generic operations: igraph_next(), igraph_end(),
 * igraph_reset(), igraph_get_vertex().
 * 
 * Specific operation: igraph_iterator_randomwalk(),
 * igraph_iterator_randomwalk_length().
 * 
 * Time complexity of igraph_next() is <code>O(d)</code>, the other
 * generic operations are <code>O(1)</code>. <code>d</code> is the
 * number of neighboring vertices at the current vertex.
 * 
 * \sa iterators_randomwalk1
 */

/**
 * \ingroup iterators_randomwalk
 * \brief Creates a random walker iterator.
 * 
 * Initializes a random walker iterator. The random walker starts at
 * the given vertex and steps to a neighboring vertex if igraph_next()
 * is called. The random walker cannot leave the starting component of
 * the graph. igraph_end() returns TRUE only if there is no possible
 * legal steps to take.
 * @param graph The graph to walk on.
 * @param it Pointer to an uninitialized iterator.
 * @param vid The id of the vertex to start from.
 * @param mode Constant giving the type of the edges to use while
 *        walking. Possible values: <b>IGRAPH_OUT</b> follows the
 *        direction of the edges, <b>IGRAPH_IN</b> follows the
 *        opposite of the direction of the edges, <b>IGRAPH_ALL</b>
 *        ignores the direction of the edges. This argument is ignored
 *        for undirected graphs.
 * @return Error code:
 *         - <b>IGRAPH_ENOMEM</b>: not enough memory.
 * 
 * Time complexity: <code>O(1)</code>.
 */

int igraph_iterator_randomwalk(igraph_t *graph, igraph_iterator_t *it,
			       integer_t vid, igraph_neimode_t mode) {
  real_t *data;
  
  if (!igraph_is_directed(graph)) {
    mode=IGRAPH_ALL;
  }
  it->type=IGRAPH_ITERATOR_RANDOMWALK;
  it->next=igraph_next_randomwalk;
  it->prev=0;
  it->end=igraph_end_randomwalk;
  it->reset=igraph_reset_randomwalk;
  it->getvertex=igraph_get_vertex_randomwalk;
  it->getvertexfrom=0;
  it->getvertexto=0;
  it->getedge=0;
  it->getvertexnei=0;
  
  it->data=Calloc(4, real_t);
  if (it->data==0) {
    IGRAPH_ERROR("out of memory", IGRAPH_ENOMEM);
  }
  data=it->data;
  data[0]=vid;			/* actual vertex */
  data[1]=vid;			/* start vertex  */
  data[2]=mode;			/* mode */
  data[3]=0;			/* number of steps so far */

  return 0;
}

/**
 * \defgroup iterators_randomwalk1 Random walker with one step memory.
 * \ingroup iterators_vertex
 * \brief Perform a random walk from a given vertex, trying not to
 * step backwards.
 * 
 * Implemented generic operations: igraph_next(), igraph_end(),
 * igraph_reset(), igraph_get_vertex().
 * 
 * Specific operations: igraph_iterator_randomwalk1(),
 * igraph_iterator_randomwalk_length().
 * 
 * Time complexity of igraph_next() is <code>O(d)</code>, the other
 * generic operations are <code>O(1)</code>. <code>d</code> is the
 * number of neighboring vertices at the current vertex.
 * 
 * \sa iterators_randomwalk
 */

/**
 * \ingroup iterators_randomwalk1
 * \brief Creates a random walker with memory.
 *
 * This is the initializer of a ``smart'' random walker, which does
 * not step backwards whenever this is possible. If the only way is
 * backwards then it steps backwards.
 * @param graph The graph to walk on.
 * @param it Pointer to an uninitialized iterator.
 * @param vid The id of the vertex to start from.
 * @param mode Constant giving the type of the edges to use while
 *        walking. Possible values: <b>IGRAPH_OUT</b> follows the
 *        direction of the edges, <b>IGRAPH_IN</b> follows the
 *        opposite of the direction of the edges, <b>IGRAPH_ALL</b>
 *        ignores the direction of the edges. This argument is ignored
 *        for undirected graphs.
 * @return Error code:
 *         - <b>IGRAPH_ENOMEM</b>: not enough memory.
 * 
 * Time complexity: <code>O(1)</code>.
 */

int igraph_iterator_randomwalk1(igraph_t *graph, igraph_iterator_t *it,
			       integer_t vid, igraph_neimode_t mode) {
  real_t *data;
  
  if (!igraph_is_directed(graph)) {
    mode=IGRAPH_ALL;
  }
  it->type=IGRAPH_ITERATOR_RANDOMWALK1;
  it->next=igraph_next_randomwalk1;
  it->prev=0;
  it->end=igraph_end_randomwalk1;
  it->reset=igraph_reset_randomwalk1;
  it->getvertex=igraph_get_vertex_randomwalk1;
  it->getvertexfrom=0;
  it->getvertexto=0;
  it->getedge=0;
  it->getvertexnei=0;
  
  it->data=Calloc(5, real_t);
  if (it->data==0) {
    IGRAPH_ERROR("out of memory", IGRAPH_ENOMEM);
  }
  data=it->data;
  data[0]=vid;			/* actual vertex */
  data[1]=vid;			/* start vertex  */
  data[2]=mode;			/* mode */
  data[3]=0;			/* number of steps so far */
  data[4]=-1;                   /* the previous vertex */

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
    if (it->data != 0) { Free(it->data); }
    break;
  case IGRAPH_ITERATOR_VNEIS:
    if (it->data != 0) { Free(it->data); }
    break;
  case IGRAPH_ITERATOR_EID:
    if (it->data != 0) { Free(it->data); }
    break;
  case IGRAPH_ITERATOR_ENEIS:
    if (it->data != 0) { Free(it->data); }
    break;
  case IGRAPH_ITERATOR_EFROMORDER:
    if (it->data != 0) { Free(it->data); }
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
  return it->prev(graph, it);
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
  real_t *data=it->data;
  data[2] ++;
  if (data[2] > VECTOR(graph->os)[ (long int)data[0]+1 ]) {
    data[3] ++;
  }
  return 0;
}

int igraph_end_vneis(igraph_t *graph, igraph_iterator_t *it) {
  real_t *data=it->data;
  return (data[2] >= VECTOR(graph->os)[ (long int) data[0]+1 ] &&
	  data[3] >= VECTOR(graph->is)[ (long int) data[0]+1 ]);  
}

int igraph_reset_vneis(igraph_t *graph, igraph_iterator_t *it) {
  real_t *data=it->data;
  if ((int)data[1] & IGRAPH_OUT) {
    data[2]=VECTOR(graph->os)[(long int)data[0]];
  } else {
    data[2]=igraph_ecount(graph);
  }
  if ((int)data[1] & IGRAPH_IN) {
    data[3]=VECTOR(graph->is)[(long int)data[0]];
  } else {
    data[3]=igraph_ecount(graph);
  }

  return 0;
}

integer_t igraph_get_vertex_vneis(igraph_t *graph, igraph_iterator_t *it) {
  real_t* data=it->data;
  if (data[2] < VECTOR(graph->os)[ (long int)data[0]+1 ]) {
    long int idx=VECTOR(graph->oi)[(long int)data[2]];
    return VECTOR(graph->to)[idx];
  } else {
    long int idx=VECTOR(graph->ii)[(long int)data[3]];
    return VECTOR(graph->from)[idx];
  }  
}

/**
 * \ingroup iterators_vneis
 * \brief Reinitialize the iterator for another vertex and/or mode.
 * 
 * This is equivalent to destroying ad inializing again the operator,
 * but there is no need to free and allocate memory, so it is slightly
 * faster and cannot signal memory allocation errors.
 * @param graph The graph to work on.
 * @param it The already initialized iterator.
 * @param vid The vertex of which the neighbors will be visited.
 * @param mode Specifies the type of neighbors to visit. See
 *        igraph_iterators_vneis() for possible values.
 * @return Error code: none right now.
 * 
 * Time complexity: <code>O(1)</code>.
 */

int igraph_iterator_vneis_set(igraph_t *graph, igraph_iterator_t *it, 
			      integer_t vid, igraph_neimode_t mode) {

  real_t *data=it->data;
  if (!igraph_is_directed(graph)) {
    mode=IGRAPH_ALL;
  }
  data[0]=vid;
  data[1]=mode;
  if (mode & IGRAPH_OUT) {
    data[2]=VECTOR(graph->os)[(long int)vid];
  } else {
    data[2]=igraph_ecount(graph);
  }
  if (mode & IGRAPH_IN) {
    data[3]=VECTOR(graph->is)[(long int)vid];
  } else {
    data[3]=igraph_ecount(graph);
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
  real_t *data=it->data;
  if (!igraph_is_directed(graph)) {
    mode=IGRAPH_ALL;
  }
  data[0]=vid;
  data[1]=mode;
  if (mode & IGRAPH_OUT) {
    data[2]=VECTOR(graph->os)[(long int)vid];
  } else {
    data[2]=igraph_ecount(graph);
  }
  if (mode & IGRAPH_IN) {
    data[3]=VECTOR(graph->is)[(long int)vid];
  } else {
    data[3]=igraph_ecount(graph);
  }
  
  return 0;
}

int igraph_reset_eneis(igraph_t *graph, igraph_iterator_t *it) {
  real_t *data=it->data;
  if ((int)data[1] & IGRAPH_OUT) {
    data[2]=VECTOR(graph->os)[(long int)data[0]];
  } else {
    data[2]=igraph_ecount(graph);
  }
  if ((int)data[1] & IGRAPH_IN) {
    data[3]=VECTOR(graph->is)[(long int)data[0]];
  } else {
    data[3]=igraph_ecount(graph);
  }
  
  return 0;
}

int igraph_next_eneis(igraph_t *graph, igraph_iterator_t *it) {
  real_t* data=it->data;
  data[2] ++; 
  if (data[2] > VECTOR(graph->os)[ (long int)data[0]+1 ]) {
    data[3] ++;
  }     
  return 0;
}

bool_t igraph_end_eneis(igraph_t *graph, igraph_iterator_t *it) {
  real_t* data=it->data;
  return (data[2] >= VECTOR(graph->os)[ (long int) data[0]+1 ] &&
	  data[3] >= VECTOR(graph->is)[ (long int) data[0]+1 ]);
}

integer_t igraph_get_vertex_from_eneis(igraph_t *graph, 
				       igraph_iterator_t *it) {
  real_t* data=it->data;
  if (data[2] < VECTOR(graph->os)[ (long int)data[0]+1 ]) {
    long int idx=VECTOR(graph->oi)[(long int)data[2]];
    return VECTOR(graph->from)[idx];
  } else {
    long int idx=VECTOR(graph->ii)[(long int)data[3]];
    return VECTOR(graph->from)[idx];
  }  
}

integer_t igraph_get_vertex_to_eneis(igraph_t *graph, igraph_iterator_t *it) {
  real_t* data=it->data;
  if (data[2] < VECTOR(graph->os)[ (long int)data[0]+1 ]) {
    long int idx=VECTOR(graph->oi)[(long int)data[2]];
    return VECTOR(graph->to)[idx];
  } else {
    long int idx=VECTOR(graph->ii)[(long int)data[3]];
    return VECTOR(graph->to)[idx];
  }  
}

integer_t igraph_get_vertex_nei_eneis(igraph_t *graph, igraph_iterator_t *it) {
  real_t* data=it->data;
  if (data[2] < VECTOR(graph->os)[ (long int)data[0]+1 ]) {
    long int idx=VECTOR(graph->oi)[(long int)data[2]];
    return VECTOR(graph->to)[idx];
  } else {
    long int idx=VECTOR(graph->ii)[(long int)data[3]];
    return VECTOR(graph->from)[idx];
  }  
}  

integer_t igraph_get_edge_eneis(igraph_t *graph, igraph_iterator_t *it) {
  real_t* data=it->data;
  if (data[2] < VECTOR(graph->os)[ (long int)data[0]+1 ]) {
    return VECTOR(graph->oi)[(long int)data[2]];
  } else {
    return VECTOR(graph->ii)[(long int)data[3]];
  }
}

int igraph_next_randomwalk(igraph_t *graph, igraph_iterator_t *it) {
  real_t *data=it->data;
  long int vid=data[0], nvid;
  igraph_neimode_t mode=data[2];
  long int indegree=0, outdegree=0;
  
  if (mode & IGRAPH_OUT) {
    outdegree += (VECTOR(graph->os)[vid+1]-VECTOR(graph->os)[vid]);
  } 
  if (mode & IGRAPH_IN) {
    indegree += (VECTOR(graph->is)[vid+1]-VECTOR(graph->is)[vid]);
  }

  if (indegree+outdegree==0) {
    /* TODO: Nowhere to step, isolate vertex. What should we do? */
    return 0;
  }
  
  RNG_BEGIN();
  nvid=RNG_INTEGER(0, outdegree+indegree-1);  
  RNG_END();

  if (nvid < outdegree) {
    long int i=VECTOR(graph->os)[vid]+nvid;
    nvid=VECTOR(graph->to)[ (long int) VECTOR(graph->oi)[i] ];
  } else {
    long int i=VECTOR(graph->is)[vid]+nvid-outdegree;
    nvid=VECTOR(graph->from)[ (long int) VECTOR(graph->ii)[i] ];
  }
  
  data[0]=nvid;
  data[3] += 1.0;

  return 0;
}

bool_t igraph_end_randomwalk(igraph_t *graph, igraph_iterator_t *it) {
  real_t *data=it->data;
  long int vid=data[0];
  igraph_neimode_t mode=data[2];
  long int indegree=0, outdegree=0;
  
  if (mode & IGRAPH_OUT) {
    outdegree += (VECTOR(graph->os)[vid+1]-VECTOR(graph->is)[vid]);
  } 
  if (mode & IGRAPH_IN) {
    indegree += (VECTOR(graph->is)[vid+1]-VECTOR(graph->is)[vid]);
  }

  return indegree+outdegree == 0;
}

integer_t igraph_get_vertex_randomwalk(igraph_t *graph, 
				       igraph_iterator_t *it) {
  real_t *data=it->data;
  return data[0];
}

int igraph_reset_randomwalk(igraph_t *graph, igraph_iterator_t *it) {
  real_t *data=it->data;
  data[0]=data[1];
  data[3]=0;
  return 0;
}

/**
 * \ingroup iterators_randomwalk
 * \brief Reports the length of a random walk.
 * 
 * @param graph The graph to work on.
 * @param it The random walker iterator.
 * @return The length of the random walk.
 * 
 * Time complexity: <code>O(1)</code>.
 */

long int igraph_iterator_randomwalk_length(igraph_t *graph, 
					   igraph_iterator_t *it) {
  real_t *data=it->data;
  return data[3];
}

int igraph_next_randomwalk1(igraph_t *graph, igraph_iterator_t *it) {
  real_t *data=it->data;
  long int vid=data[0], nvid;
  igraph_neimode_t mode=data[2];
  long int indegree=0, outdegree=0;

  if (mode & IGRAPH_OUT) {
    outdegree += (VECTOR(graph->os)[vid+1]-VECTOR(graph->os)[vid]);
  } 
  if (mode & IGRAPH_IN) {
    indegree += (VECTOR(graph->is)[vid+1]-VECTOR(graph->is)[vid]);
  }

  if (indegree+outdegree==0) {
    /* TODO: Nowhere to step, isolate vertex. What should we do? */
    return 0;
  }
  
  if (indegree+outdegree==1) {
    /* There only one way, we may go back as well */
    nvid=0;
    
    if (nvid < outdegree) {
      long int i=VECTOR(graph->os)[vid]+nvid;
      nvid=VECTOR(graph->to)[ (long int) VECTOR(graph->oi)[i] ];
    } else {
      long int i=VECTOR(graph->is)[vid]+nvid-outdegree;
      nvid=VECTOR(graph->from)[ (long int) VECTOR(graph->ii)[i] ];
    }
    
  } else {
    /* There are other options */
    RNG_BEGIN();
    if (data[4] >=0) {
      nvid=RNG_INTEGER(0, outdegree+indegree-2);
    } else {
      nvid=RNG_INTEGER(0, outdegree+indegree-1);
    }      
    RNG_END();

    if (nvid < outdegree) {
      long int i=VECTOR(graph->os)[vid]+nvid;
      nvid=VECTOR(graph->to)[ (long int) VECTOR(graph->oi)[i] ];
    } else {
      long int i=VECTOR(graph->is)[vid]+nvid-outdegree;
      nvid=VECTOR(graph->from)[ (long int) VECTOR(graph->ii)[i] ];
    }
    
    /* In case we wanted to step back */
    /* TODO: it is still possible to step back but on a *different* edge, 
       this applies only to graphs with multiple edges of course */
    if (nvid==data[4]) {
      nvid=outdegree+indegree-1;

      if (nvid < outdegree) {
	long int i=VECTOR(graph->os)[vid]+nvid;
	nvid=VECTOR(graph->to)[ (long int) VECTOR(graph->oi)[i] ];
      } else {
	long int i=VECTOR(graph->is)[vid]+nvid-outdegree;
	nvid=VECTOR(graph->from)[ (long int) VECTOR(graph->ii)[i] ];
      }
    }
  }

  data[4]=data[0];
  data[0]=nvid;
  data[3]+=1.0;

  return 0;
}

bool_t igraph_end_randomwalk1(igraph_t *graph, igraph_iterator_t *it) {
  real_t *data=it->data;
  long int vid=data[0];
  igraph_neimode_t mode=data[2];
  long int indegree=0, outdegree=0;
  
  if (mode & IGRAPH_OUT) {
    outdegree += (VECTOR(graph->os)[vid+1]-VECTOR(graph->is)[vid]);
  } 
  if (mode & IGRAPH_IN) {
    indegree += (VECTOR(graph->is)[vid+1]-VECTOR(graph->is)[vid]);
  }

  return indegree+outdegree == 0;
}

integer_t igraph_get_vertex_randomwalk1(igraph_t *graph, 
				       igraph_iterator_t *it) {
  real_t *data=it->data;
  return data[0];
}

int igraph_reset_randomwalk1(igraph_t *graph, igraph_iterator_t *it) {
  real_t *data=it->data;
  data[0]=data[1];
  data[3]=0;
  data[5]=-1;
  return 0;
}

