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

#include <string.h>
#include <stdarg.h>

/**
 * \section about_iterators
 * 
 * <para>Vertex and edge sequences are probably the most central concepts of
 * \a igraph, They provide a method to step over a sequence of vertices
 * or edges. There are different types of vertex sequences (and this is
 * also true for edges sequences, like everything else in this section),
 * and they are all handled through a common abstract interface. </para>
 * 
 * <para>Vertex sequences are created by various vertex sequences
 * constructors, eg. \ref igraph_vs_all() creates a vertex sequence
 * containing all vertices of a graph, in increasing order.
 * We will discuss all vertex sequence constructors later.</para>
 * 
 * <para>Vertex sequences are usually used as iterators, but there is a
 * method unfold finite vertex sequences into a \type igraph_vector_t
 * object. (Infinite sequences cannot be unfolded of course.) 
 * This manual sometimes calles vertex (or edge) sequences as
 * iterators, these two names are used as synonims. </para>
 * 
 * <para>After creation the vertex sequence is set to \quote point\endquote
 * to the first vertex (if there are any in the sequence), the actual
 * vertex id can be queried by \ref igraph_vs_get(). \ref
 * igraph_vs_next() moves on to the next vertex in the sequence. \ref
 * igraph_vs_end() checks whether there are more vertices in the
 * sequence. \ref igraph_vs_destroy() should be called once the vertex
 * sequence is not needed any more, this deallocates the memory area
 * used by the iterator.</para>
 * 
 * <para>Every vertex sequence is related to an \type igraph_t
 * object, this is called the \em underlying graph. If the structure
 * of the underlying graph changes or it is destroyed by calling \ref
 * igraph_destroy() on it, the vertex sequence is said to be
 * invalidated. The result of all function calls except \ref
 * igraph_vs_destroy() on invalidated iterators is undefined.</para>
 * 
 * <para>Every vertex sequence type has a separate constructor, and some of
 * them also have special functions only meaningful if called with a
 * certain type of vertex sequence (unlike the general abstract
 * interface). These are discussed in the next section.</para>
 */

/**
 * \section about_iterator_shorthands
 * 
 * <para>There are a number of functions in the \a igraph library which take
 * vertex or edge sequences as parameters. For example the simple \ref
 * igraph_degree() function calculates the degree of the vertices
 * given as a vertex sequence parameter. Imagine that you want to
 * implement an algorithm which needs the degree of some vertex from
 * time to time for some purpose. Now what you would normally do is to
 * create a vertex sequence containing only one vertex every time you
 * need a vertex degree, call \ref igraph_degree() and destroy the
 * vertex sequence. This would be not only inefficient but also quite
 * annoying since you needed to write three statements (initialize the
 * vertex sequence, call igraph_degree(), destroy the vertex sequence)
 * instead of just one. </para>
 * 
 * <para>Vertex and edge sequence shorthands come over this. They allow you
 * to use a shorter (and usually also more efficient)
 * notation. Iterator shorthands can be used in the parameter list of
 * the function you want to call, they provide a compact
 * notation. Here is a little example for calculating the closeness
 * centraility for every vertex in graph <code>g</code> by
 * calling \ref igraph_closeness():
 * \verbatim ret = igraph_closeness(g, res, IGRAPH_VS_ALL(g), IGRAPH_ALL); \endverbatim </para>
 * 
 * <para>Please use shorthands only in the parameter list of a function, 
 * otherwise your program might contain memory leaks. (This is because
 * the memory allocated for the shorthand is freed by the function you
 * call.)</para>
 *
 * <para>Also use this notation only if the function you call handles
 * shorthands, as they require special traitment.</para>
 */

/**
 * \section iterator_examples
 * 
 * <para>TODO</para>
 */

/* -------------------------------------------------- */
/* Vertex iterator generics                           */
/* -------------------------------------------------- */

/**
 * \function igraph_vs_next
 * 
 * Steps to the next vertex in the sequence. 
 * You should call \ref igraph_vs_end() before calling this function 
 * to be sure that there is really a next vertex. 
 * \param graph The underlying graph object.
 * \param vs The vertex sequence.
 * 
 * Time complexity: usually O(1) 
 * but see also the specific types.
 */

void igraph_vs_next(const igraph_t *graph, igraph_vs_t *vs) {
  vs->table->next(graph, (struct igraph_vs_t*)vs);
}

/**
 * \function igraph_vs_end
 * 
 * Checks whether there are more vertices in the vertex sequence. 
 * Do not call \ref igraph_vs_next() or \ref igraph_vs_get() on vertex 
 * sets for which this function returns true.
 * \param graph The underlying graph object.
 * \param vs The vertex sequence. 
 * \return True if there are no more vertices in the sequence, ie. the 
 *   pointer points to the element \em after the last one.
 *   False otherwise. 
 * 
 * Time complexity: usually O(1) 
 * but see also the specific types.
 */

bool_t igraph_vs_end(const igraph_t *graph, const igraph_vs_t *vs) {
  return vs->table->end(graph, (struct igraph_vs_t*)vs);
}

/**
 * \function igraph_vs_reset
 * 
 * This function resets the iterator position of a vertex
 * sequence. This has the same effect as destroying and
 * reinitializing the iterator, only more efficient. Note that this
 * does \em not neccessarily means that the iterator touches the same
 * vertices, eg. the random walker iterator (\ref igraph_vs_rw()) will
 * simply start a new random walk from the same vertex, but will
 * usually produce a different sequence of vertices.
 * \param graph The underlying graph object.
 * \param vs The vertex sequence object.
 * 
 * Time complexity: usually O(1) 
 * but see also the specific types.
 */

void igraph_vs_reset(const igraph_t *graph, igraph_vs_t *vs) {
  vs->table->reset(graph, (struct igraph_vs_t*)vs);
}

/**
 * \function igraph_vs_get
 * 
 * This function returns the id of the current vertex an iterator. (Or
 * in other words the vertex the iterator \quote points\endquote to.
 * You should call \ref igraph_vs_end() on the iterator beforehand to
 * be sure that the iterator really points to a vertex.
 * \param graph The underlying graph object.
 * \param vs The vertex sequence object.
 * \return Id of the \quote current\endquote vertex.
 * 
 * Time complexity: usually O(1) 
 * but see also the specific types.
 */

integer_t igraph_vs_get(const igraph_t *graph, const igraph_vs_t *vs) {
  return vs->table->get(graph, (struct igraph_vs_t*)vs);
}

/**
 * \function igraph_vs_unfold
 * 
 * Steps over a vertex sequence and stores the visited vertices in a
 * vector. Note that it is a runtime error to try to unfold infinite
 * iterators, like random walkers.
 * \param graph The underlying graph object.
 * \param vs The vertex sequence object.
 * \param v An initialized \type igraph_vector_t object. These will be
 *   resized to hold the vertices.
 * \return Error code, \c IGRAPH_ENOMEM if there is
 *   not enough memory to resize \p v.
 * 
 * Time complexity: depends on the type of the iterator, but usually
 * it is at most O(n), the number of
 * visited vertices.
 */

int igraph_vs_unfold(const igraph_t *graph, const igraph_vs_t *vs, 
		     igraph_vector_t *v) {
  return vs->table->unfold(graph, (const struct igraph_vs_t*)vs, v);
}

/**
 * \function igraph_vs_destroy
 *
 * Destroys a vertex sequence, ie. frees all memory allocated for it. 
 * The iterator has to be reinitialized before using it again.
 * \param vs The vertex sequence object to destroy.
 * 
 * Time complexity: usually O(1) 
 * but see also the specific types.
 */

void igraph_vs_destroy(igraph_vs_t *vs) {
  vs->table->destroy((struct igraph_vs_t*)vs);
}

/* -------------------------------------------------- */
/* Simple vertex iterator                             */
/* -------------------------------------------------- */

void igraph_vs_next_all(const igraph_t *graph, igraph_vs_t *vs);
bool_t igraph_vs_end_all(const igraph_t *graph, const igraph_vs_t *vs);
void igraph_vs_reset_all(const igraph_t *graph, igraph_vs_t *vs);
integer_t igraph_vs_get_all(const igraph_t *graph, const igraph_vs_t *vs);
int igraph_vs_unfold_all(const igraph_t *graph, const igraph_vs_t *vs, 
			 igraph_vector_t *v);
void igraph_vs_destroy_all(igraph_vs_t *vs);

igraph_i_vstable_t igraph_i_vs_all_table = {
  igraph_vs_next_all, igraph_vs_end_all, igraph_vs_reset_all,
  igraph_vs_get_all, igraph_vs_unfold_all, igraph_vs_destroy_all
};

/**
 * \function igraph_vs_all 
 * 
 * Vertex sequence containing all vertices in increasing vertex id
 * order. This iterator doesn't have any special functions.
 * \param graph The underlying graph object.
 * \param vs Pointer to an uninitialized vertex set object.
 * \return Error code, the current implementation always returns with
 *   success. 
 * 
 * Time complexity: O(1).
 */

int igraph_vs_all(const igraph_t *graph, igraph_vs_t *vs) {
  vs->type=IGRAPH_ITERATOR_VS_ALL;
  vs->stdata[0]=0;
  vs->table=&igraph_i_vs_all_table;
  vs->shorthand=0;
  return 0;
}

/**
 * \function IGRAPH_VS_ALL
 * 
 * This is a shorthand for all vertices in a graph.
 * \param graph The underlying graph object.
 * \return Vertex sequence shorthand.
 *
 * Time complexity: O(1).
 */

const igraph_vs_t *IGRAPH_VS_ALL(const igraph_t *graph) {
  igraph_vs_t *vs=Calloc(1, igraph_vs_t);
  if (vs==0) {
    igraph_error("Cannot create iterator shorthand", __FILE__, __LINE__,
		 IGRAPH_ENOMEM);
    return 0;			/* TODO: how to sign error ??? */
  }
  igraph_vs_all(graph, vs);
  vs->shorthand=1;
  return vs;
}

void igraph_vs_next_all(const igraph_t *graph, igraph_vs_t *vs) {
  vs->stdata[0] ++;
}

bool_t igraph_vs_end_all(const igraph_t *graph, const igraph_vs_t *vs) {
  return vs->stdata[0] >= igraph_vcount(graph);
}

void igraph_vs_reset_all(const igraph_t *graph, igraph_vs_t *vs) {
  vs->stdata[0]=0;
}

integer_t igraph_vs_get_all(const igraph_t *graph, const igraph_vs_t *vs) {
  return vs->stdata[0];
}

int igraph_vs_unfold_all(const igraph_t *graph, const igraph_vs_t *vs, 
			 igraph_vector_t *v) {
  long int n;
  igraph_vector_t v2;
  n=igraph_vcount(graph);
  IGRAPH_CHECK(igraph_vector_init_seq(&v2, 0, n-1));
  igraph_vector_destroy(v);
  *v=v2;
  return 0;
}

void igraph_vs_destroy_all(igraph_vs_t *pvs) {
  igraph_vs_t *vs=(igraph_vs_t*)pvs;
  if (vs->shorthand) {
    Free(pvs);
  }
}

/* -------------------------------------------------- */
/* Adjacent vertices of a vertex                      */
/* -------------------------------------------------- */

void igraph_vs_next_adj(const igraph_t *graph, igraph_vs_t *vs);
bool_t igraph_vs_end_adj(const igraph_t *graph, const igraph_vs_t *vs);
void igraph_vs_reset_adj(const igraph_t *graph, igraph_vs_t *vs);
integer_t igraph_vs_get_adj(const igraph_t *graph, const igraph_vs_t *vs);
int igraph_vs_unfold_adj(const igraph_t *graph, const igraph_vs_t *vs, 
			 igraph_vector_t *v);
void igraph_vs_destroy_adj(igraph_vs_t *vs);

igraph_i_vstable_t igraph_i_vs_adj_table = {
  igraph_vs_next_adj, igraph_vs_end_adj, igraph_vs_reset_adj,
  igraph_vs_get_adj, igraph_vs_unfold_adj, igraph_vs_destroy_adj
};

/**
 * \function igraph_vs_adj
 *
 * A sequence containing all vertices connected by an edge to a
 * specified vertex.
 * \param graph The underlying graph object.
 * \param vs Pointer to an uninitialized vertex sequence object.
 * \param vid The vertex id of the vertex of which the neighboring
 *   vertices will be searched.
 * \param mode \c IGRAPH_OUT for vertices at the
 *   end of an edge originating from \p vid,
 *   \c IGRAPH_IN for vertices from which an edge
 *   points to \p vid. 
 *   \c IGRAPH_ALL for the union of the previous
 *   two. This parameter is ignored for undirected graphs.
 * \return Error code, the current implementation always returns with
 *   success.
 *
 * Time complexity: O(1).
 */

int igraph_vs_adj(const igraph_t *graph, igraph_vs_t *vs,
		   integer_t vid, igraph_neimode_t mode) {
  if (!igraph_is_directed(graph)) {
    mode=IGRAPH_ALL;
  }
  vs->type=IGRAPH_ITERATOR_VS_ADJ;
  vs->table=&igraph_i_vs_adj_table;
  vs->shorthand=0;
  
  vs->stdata[0]=vid;
  vs->stdata[1]=mode;
  if (mode & IGRAPH_OUT) {
    vs->stdata[2]=VECTOR(graph->os)[(long int)vid];
  } else {
    vs->stdata[2]=igraph_ecount(graph);
  }
  if (mode & IGRAPH_IN) {
    vs->stdata[3]=VECTOR(graph->is)[(long int)vid];
  } else {
    vs->stdata[3]=igraph_ecount(graph);
  }
  vs->stdata[4]=VECTOR(graph->os)[ (long int) vs->stdata[0]+1 ];
  vs->stdata[5]=VECTOR(graph->is)[ (long int) vs->stdata[0]+1 ];
  
  return 0;
}

void igraph_vs_next_adj(const igraph_t *graph, igraph_vs_t *vs) {
  vs->stdata[2] ++;
  if (vs->stdata[2] > vs->stdata[4]) {
    vs->stdata[3] ++;
  }
}

bool_t igraph_vs_end_adj(const igraph_t *graph, const igraph_vs_t *vs) {
  return (vs->stdata[2] >= vs->stdata[4] &&
	  vs->stdata[3] >= vs->stdata[5]);
}

void igraph_vs_reset_adj(const igraph_t *graph, igraph_vs_t *vs) {
  if ((int)vs->stdata[1] & IGRAPH_OUT) {
    vs->stdata[2]=VECTOR(graph->os)[(long int)vs->stdata[0]];
  } else {
    vs->stdata[2]=igraph_ecount(graph);
  }
  if ((int)vs->stdata[1] & IGRAPH_IN) {
    vs->stdata[3]=VECTOR(graph->is)[(long int)vs->stdata[0]];
  } else {
    vs->stdata[3]=igraph_ecount(graph);
  }
}

integer_t igraph_vs_get_adj(const igraph_t *graph, const igraph_vs_t *vs) {
  if (vs->stdata[2] < vs->stdata[4]) {
    long int idx=VECTOR(graph->oi)[(long int)vs->stdata[2]];
    return VECTOR(graph->to)[idx];
  } else {
    long int idx=VECTOR(graph->ii)[(long int)vs->stdata[3]];
    return VECTOR(graph->from)[idx];
  }
}

int igraph_vs_unfold_adj(const igraph_t *graph, const igraph_vs_t *vs,
			 igraph_vector_t *v) {
  IGRAPH_CHECK(igraph_neighbors(graph, v, vs->stdata[0], vs->stdata[1]));
  return 0;
}

/**
 * \function igraph_vs_adj_set
 * 
 * Reinitializes an iterator created by \ref igraph_vs_adj() for a
 * different vertex. This function is equivialent to destroying the
 * iterator and initializing it with a different
 * \p vid argument but it is more efficient.
 * 
 * The arguments are the same as for the \ref igraph_vs_adj()
 * function. 
 * 
 * Time complexity: O(1).
 */

void igraph_vs_adj_set(const igraph_t *graph, igraph_vs_t *vs,
			integer_t vid, igraph_neimode_t mode) {
  if (!igraph_is_directed(graph)) {
    mode=IGRAPH_ALL;
  }
  vs->stdata[0]=vid;
  vs->stdata[1]=mode;
  if (mode & IGRAPH_OUT) {
    vs->stdata[2]=VECTOR(graph->os)[(long int)vid];
  } else {
    vs->stdata[2]=igraph_ecount(graph);
  }
  if (mode & IGRAPH_IN) {
    vs->stdata[3]=VECTOR(graph->is)[(long int)vid];
  } else {
    vs->stdata[3]=igraph_ecount(graph);
  }
  vs->stdata[4]=VECTOR(graph->os)[ (long int) vs->stdata[0]+1 ];
  vs->stdata[5]=VECTOR(graph->is)[ (long int) vs->stdata[0]+1 ];
}

void igraph_vs_destroy_adj(igraph_vs_t *vs) {
}  

/* -------------------------------------------------- */
/* Random walker                                      */
/* -------------------------------------------------- */

void igraph_vs_next_rw(const igraph_t *graph, igraph_vs_t *vs);
bool_t igraph_vs_end_rw(const igraph_t *graph, const igraph_vs_t *vs);
void igraph_vs_reset_rw(const igraph_t *graph, igraph_vs_t *vs);
integer_t igraph_vs_get_rw(const igraph_t *graph, const igraph_vs_t *vs);
int igraph_vs_unfold_rw(const igraph_t *graph, const igraph_vs_t *vs, 
			igraph_vector_t *v);
void igraph_vs_destroy_rw(igraph_vs_t *vs);

igraph_i_vstable_t igraph_i_vs_rw_table = {
  igraph_vs_next_rw, igraph_vs_end_rw, igraph_vs_reset_rw,
  igraph_vs_get_rw, igraph_vs_unfold_rw, igraph_vs_destroy_rw
};

/**
 * \function igraph_vs_rw
 * 
 * A random walker vertex iterator without memory.
 * The random walker starts at the supplied vertex and steps to
 * another vertex along an edge randomly.
 * 
 * This is usually an infinite iterator, \ref igraph_vs_end() always
 * returns false, and it is a runtime error to attempt to unfold it.
 * \param graph The underlying graph object.
 * \param vs Pointer to an uninitialized vertex sequence object.
 * \param vid The id of the starting vertex.
 * \param mode Gives how to step along directed
 *   edges. \c IGRAPH_OUT moves along the direction
 *   of the edges, \c IGRAPH_IN the opposite
 *   way. \c IGRAPH_ALL ignored the direction of the
 *   edges. This argument is ignored for undirected graphs.
 * \return Error code, the current implementation always returns with
 *   success. 
 * 
 * Time complexity: O(1).
 */

int igraph_vs_rw(const igraph_t *graph, igraph_vs_t *vs,
		 integer_t vid, igraph_neimode_t mode) {
  if (!igraph_is_directed(graph)) {
    mode=IGRAPH_ALL;
  }
  vs->type=IGRAPH_ITERATOR_VS_RW;
  vs->shorthand=0;
  
  vs->stdata[0]=vid;		/* actual vertex */
  vs->stdata[1]=vid;		/* start vertex  */
  vs->stdata[2]=mode;		/* mode */
  vs->stdata[3]=0;		/* number of steps so far */

  vs->table=&igraph_i_vs_rw_table;
  return 0;
}

void igraph_vs_next_rw(const igraph_t *graph, igraph_vs_t *vs) {
  long int vid=vs->stdata[0], nvid;
  igraph_neimode_t mode=vs->stdata[2];
  long int indegree=0, outdegree=0;
  
  if (mode & IGRAPH_OUT) {
    outdegree += (VECTOR(graph->os)[vid+1]-VECTOR(graph->os)[vid]);
  }
  if (mode & IGRAPH_IN) {
    indegree += (VECTOR(graph->is)[vid+1]-VECTOR(graph->is)[vid]);
  }

  if (indegree+outdegree==0) {
    /* TODO: Nowhere to step, isolate vertex. What should we do? */
    return;
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
  
  vs->stdata[0]=nvid;
  vs->stdata[3] += 1.0;
}

bool_t igraph_vs_end_rw(const igraph_t *graph, const igraph_vs_t *vs) {
  long int vid=vs->stdata[0];
  igraph_neimode_t mode=vs->stdata[2];
  long int indegree=0, outdegree=0;
  
  if (mode & IGRAPH_OUT) {
    outdegree += (VECTOR(graph->os)[vid+1]-VECTOR(graph->is)[vid]);
  }

  if (mode & IGRAPH_IN) {
    indegree += (VECTOR(graph->is)[vid+1]-VECTOR(graph->is)[vid]);
  }

  return indegree+outdegree == 0;
}

void igraph_vs_reset_rw(const igraph_t *graph, igraph_vs_t *vs) {
  vs->stdata[0]=vs->stdata[1];
  vs->stdata[3]=0;
}

integer_t igraph_vs_get_rw(const igraph_t *graph, const igraph_vs_t *vs) {
  return vs->stdata[0];
}

int igraph_vs_unfold_rw(const igraph_t *graph, const igraph_vs_t *vs,
			igraph_vector_t *v) {
  IGRAPH_ERROR("attempt to unfold random walker", IGRAPH_EUNFOLDINF);
  return 0;			/* return unneccesary */
}

long int igraph_vs_rw_length(const igraph_t *graph, const igraph_vs_t *vs) {
  return vs->stdata[3];
}

void igraph_vs_destroy_rw(igraph_vs_t *pvs) {
  igraph_vs_t *vs=(igraph_vs_t*)pvs;
  if (vs->shorthand) {
    Free(pvs);
  }
}

/* -------------------------------------------------- */
/* Random walker with one unit memory                 */
/* -------------------------------------------------- */

void igraph_vs_next_rw1(const igraph_t *graph, igraph_vs_t *vs);
bool_t igraph_vs_end_rw1(const igraph_t *graph, const igraph_vs_t *vs);
void igraph_vs_reset_rw1(const igraph_t *graph, igraph_vs_t *vs);
integer_t igraph_vs_get_rw1(const igraph_t *graph, const igraph_vs_t *vs);
int igraph_vs_unfold_rw1(const igraph_t *graph, const igraph_vs_t *vs, 
			 igraph_vector_t *v);
void igraph_vs_destroy_rw1(igraph_vs_t *vs);

igraph_i_vstable_t igraph_i_vs_rw1_table = {
  igraph_vs_next_rw1, igraph_vs_end_rw1, igraph_vs_reset_rw1,
  igraph_vs_get_rw1, igraph_vs_unfold_rw1, igraph_vs_destroy_rw1
};

/**
 * \function igraph_vs_rw1
 * 
 * A random walker vertex iterator with one unit memory.
 * This random walker starts at the supplied vertex and steps to
 * another vertex along an edge randomly avoiding stepping back.
 * The iterators only steps backwards if this is the only way it can
 * go. 
 * 
 * This is usually an infinite iterator, \ref igraph_vs_end() always
 * returns false, and it is a runtime error to attempt to unfold it.
 * \param graph The underlying graph object.
 * \param vs Pointer to an uninitialized vertex sequence object.
 * \param vid The id of the starting vertex.
 * \param mode Gives how to step along directed
 *   edges. \c IGRAPH_OUT moves along the direction
 *   of the edges, \c IGRAPH_IN the opposite
 *   way. \c IGRAPH_ALL ignored the direction of the
 *   edges. This argument is ignored for undirected graphs.
 * \return Error code, the current implementation always returns with
 *   success. 
 * 
 * Time complexity: O(1).
 */

int igraph_vs_rw1(const igraph_t *graph, igraph_vs_t *vs,
		  integer_t vid, igraph_neimode_t mode) {
  if (!igraph_is_directed(graph)) {
    mode=IGRAPH_ALL;
  }
  vs->type=IGRAPH_ITERATOR_VS_RW1;
  vs->shorthand=0;
  
  vs->stdata[0]=vid;			/* actual vertex */
  vs->stdata[1]=vid;			/* start vertex  */
  vs->stdata[2]=mode;			/* mode */
  vs->stdata[3]=0;			/* number of steps so far */
  vs->stdata[4]=-1;		        /* the previous vertex */

  vs->table=&igraph_i_vs_rw1_table;
  return 0;
}

void igraph_vs_next_rw1(const igraph_t *graph, igraph_vs_t *vs) {

  long int vid=vs->stdata[0], nvid;
  igraph_neimode_t mode=vs->stdata[2];
  long int indegree=0, outdegree=0;

  if (mode & IGRAPH_OUT) {
    outdegree += (VECTOR(graph->os)[vid+1]-VECTOR(graph->os)[vid]);
  }
  if (mode & IGRAPH_IN) {
    indegree += (VECTOR(graph->is)[vid+1]-VECTOR(graph->is)[vid]);
  }

  if (indegree+outdegree==0) {
    /* TODO: Nowhere to step, isolate vertex. What should we do? */
    return;
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
    if (vs->stdata[4] >=0) {
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
    if (nvid==vs->stdata[4]) {
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

  vs->stdata[4]=vs->stdata[0];
  vs->stdata[0]=nvid;
  vs->stdata[3]+=1.0;
}

bool_t igraph_vs_end_rw1(const igraph_t *graph, const igraph_vs_t *vs) {
  long int vid=vs->stdata[0];
  igraph_neimode_t mode=vs->stdata[2];
  long int indegree=0, outdegree=0;
  
  if (mode & IGRAPH_OUT) {
    outdegree += (VECTOR(graph->os)[vid+1]-VECTOR(graph->is)[vid]);
  }
  if (mode & IGRAPH_IN) {
    indegree += (VECTOR(graph->is)[vid+1]-VECTOR(graph->is)[vid]);
  }

  return indegree+outdegree == 0;
}

integer_t igraph_vs_get_rw1(const igraph_t *graph, const igraph_vs_t *vs) {
  return vs->stdata[0];
}

void igraph_vs_reset_rw1(const igraph_t *graph, igraph_vs_t *vs) {
  vs->stdata[0]=vs->stdata[1];
  vs->stdata[3]=0;
  vs->stdata[5]=-1;
}

int igraph_vs_unfold_rw1(const igraph_t *graph, const igraph_vs_t *vs,
			 igraph_vector_t *v) {
  IGRAPH_ERROR("attempt to unfold random walker", IGRAPH_EUNFOLDINF);
  return 0;
}

long int igraph_vs_rw1_length(const igraph_t *graph, const igraph_vs_t *vs) {
  return vs->stdata[3];
}

void igraph_vs_destroy_rw1(igraph_vs_t *pvs) {
  if (pvs->shorthand) {
    Free(pvs);
  }
}

/* -------------------------------------------------- */
/* Empty vertex set                                   */
/* -------------------------------------------------- */

void igraph_vs_next_none(const igraph_t *graph, igraph_vs_t *vs);
bool_t igraph_vs_end_none(const igraph_t *graph, const igraph_vs_t *vs);
void igraph_vs_reset_none(const igraph_t *graph, igraph_vs_t *vs);
integer_t igraph_vs_get_none(const igraph_t *graph, const igraph_vs_t *vs);
int igraph_vs_unfold_none(const igraph_t *graph, const igraph_vs_t *vs,
			  igraph_vector_t *v);
void igraph_vs_destroy_none(igraph_vs_t *vs);

igraph_i_vstable_t igraph_i_vs_none_table = {
  igraph_vs_next_none, igraph_vs_end_none, igraph_vs_reset_none,
  igraph_vs_get_none, igraph_vs_unfold_none, igraph_vs_destroy_none
};

/**
 * \function igraph_vs_none
 * 
 * Creates an empty vertex sequence.
 * \param graph The underlying graph object.
 * \param vs Pointer to an uninitialized vertex sequence.
 * \return Error code, the current implementation always returns with
 *   success. 
 * 
 * Time complexity: O(1).
 */

int igraph_vs_none(const igraph_t *graph, igraph_vs_t *vs) {
  vs->type=IGRAPH_ITERATOR_VS_NONE;
  vs->table=&igraph_i_vs_none_table;
  vs->shorthand=0;
  return 0;
}

void igraph_vs_next_none(const igraph_t *graph, igraph_vs_t *vs) {
  /* nothing to do */
}

bool_t igraph_vs_end_none(const igraph_t *graph, const igraph_vs_t *vs) {
  return 1;
}

void igraph_vs_reset_none(const igraph_t *graph, igraph_vs_t *vs) {
  /* nothing to do */
}

integer_t igraph_vs_get_none(const igraph_t *graph, const igraph_vs_t *vs) {
  /* ooops this is an error, no way to signal it though... */
  return -1;
}

int igraph_vs_unfold_none(const igraph_t *graph, const igraph_vs_t *vs,
			  igraph_vector_t *v) {
  igraph_vector_clear(v);
  return 0;
}

void igraph_vs_destroy_none(igraph_vs_t *vs) {
  if (vs->shorthand) {
    Free(vs);
  }
}

/* -------------------------------------------------- */
/* vertex set with single vertex                      */
/* -------------------------------------------------- */

void igraph_vs_next_1(const igraph_t *graph, igraph_vs_t *vs);
bool_t igraph_vs_end_1(const igraph_t *graph, const igraph_vs_t *vs);
void igraph_vs_reset_1(const igraph_t *graph, igraph_vs_t *vs);
integer_t igraph_vs_get_1(const igraph_t *graph, const igraph_vs_t *vs);
int igraph_vs_unfold_1(const igraph_t *graph, const igraph_vs_t *vs, 
		       igraph_vector_t *v);
void igraph_vs_destroy_1(igraph_vs_t *vs);

igraph_i_vstable_t igraph_i_vs_1_table = {
  igraph_vs_next_1, igraph_vs_end_1, igraph_vs_reset_1,
  igraph_vs_get_1, igraph_vs_unfold_1, igraph_vs_destroy_1
};

/**
 * \function igraph_vs_1
 * 
 * Creates a vertex sequence with a single vertex.
 * \param igraph The underlying graph object.
 * \param vs Pointer to an uninitialized vertex sequence.
 * \param vid The id of the only vertex in the sequence. 
 * \return Error code, the current implementation always returns with
 *   success. 
 * 
 * Time complexity: O(1).
 */

int igraph_vs_1(const igraph_t *igraph, igraph_vs_t *vs, integer_t vid) {
  vs->type=IGRAPH_ITERATOR_VS_1;
  vs->stdata[0]=vid;
  vs->stdata[1]=0;		/* write 1 here if end */
  vs->table=&igraph_i_vs_1_table;
  vs->shorthand=0;
  return 0;
}

/**
 * \function IGRAPH_VS_1
 * 
 * Vertex sequence shorthand for a single vertex.
 * \param graph The underlying graph object.
 * \param vid The single vertex id to be included in the vertex
 * sequence. 
 * \return Vertex sequence shorthand.
 * 
 * Time complexity: O(1).
 */

const igraph_vs_t *IGRAPH_VS_1(const igraph_t *graph, integer_t vid) {
  igraph_vs_t *vs=Calloc(1, igraph_vs_t);
  if (vs==0) {
    igraph_error("Cannot create iterator shorthand", __FILE__, __LINE__,
		 IGRAPH_ENOMEM);
    return 0;			/* TODO: how to sign error ??? */
  }
  igraph_vs_1(graph, vs, vid);
  vs->shorthand=1;
  return vs;
}

void igraph_vs_next_1(const igraph_t *graph, igraph_vs_t *vs) {
  /* signal end */
  vs->stdata[1]=1;
}

bool_t igraph_vs_end_1(const igraph_t *graph, const igraph_vs_t *vs) {
  return (vs->stdata[1]==1);
}

void igraph_vs_reset_1(const igraph_t *graph, igraph_vs_t *vs) {
  vs->stdata[1]=0;
}

integer_t igraph_vs_get_1(const igraph_t *graph, const igraph_vs_t *vs) {
  return vs->stdata[0];
}

int igraph_vs_unfold_1(const igraph_t *graph, const igraph_vs_t *vs,
		       igraph_vector_t *v) {
  IGRAPH_CHECK(igraph_vector_resize(v, 1));
  VECTOR(*v)[0]=vs->stdata[0];
  return 0;
}

void igraph_vs_destroy_1(igraph_vs_t *pvs) {
  igraph_vs_t *vs=(igraph_vs_t*)pvs;
  if (vs->shorthand) {
    Free(pvs);
  }
}

/* -------------------------------------------------- */
/* Vertex set with sequence of vertices               */
/* -------------------------------------------------- */

void igraph_vs_next_seq(const igraph_t *graph, igraph_vs_t *vs);
bool_t igraph_vs_end_seq(const igraph_t *graph, const igraph_vs_t *vs);
void igraph_vs_reset_seq(const igraph_t *graph, igraph_vs_t *vs);
integer_t igraph_vs_get_seq(const igraph_t *graph, const igraph_vs_t *vs);
int igraph_vs_unfold_seq(const igraph_t *graph, const igraph_vs_t *vs,
			 igraph_vector_t *v);
void igraph_vs_destroy_seq(igraph_vs_t *vs);

igraph_i_vstable_t igraph_i_vs_seq_table = {
  igraph_vs_next_seq, igraph_vs_end_seq, igraph_vs_reset_seq,
  igraph_vs_get_seq, igraph_vs_unfold_seq, igraph_vs_destroy_seq
};

/**
 * \function igraph_vs_seq
 * 
 * Creates a regular sequence of vertices. 
 * \param igraph The underlying graph object. 
 * \param vs Pointer to an uninitialized vertex sequence object.
 * \param from The lower limit of the inverval (inclusive).
 * \param to The upper limit of the interval (inclusive).
 * \return Error code, the current implementation always returns with
 *   success. 
 * 
 * Time complexity: O(1).
 */

int igraph_vs_seq(const igraph_t *igraph, igraph_vs_t *vs, integer_t from,
		  integer_t to) {
  vs->type=IGRAPH_ITERATOR_VS_SEQ;
  vs->stdata[0]=from;
  vs->stdata[1]=from;
  vs->stdata[2]=to;
  vs->table=&igraph_i_vs_seq_table;
  vs->shorthand=0;
  return 0;
}

void igraph_vs_next_seq(const igraph_t *graph, igraph_vs_t *vs) {
  vs->stdata[0] ++;
}

bool_t igraph_vs_end_seq(const igraph_t *graph, const igraph_vs_t *vs) {
  return vs->stdata[0] > vs->stdata[2];
}

void igraph_vs_reset_seq(const igraph_t *graph, igraph_vs_t *vs) {
  vs->stdata[0]=vs->stdata[1];
}

integer_t igraph_vs_get_seq(const igraph_t *graph, const igraph_vs_t *vs) {
  return vs->stdata[0];
}

int igraph_vs_unfold_seq(const igraph_t *graph, const igraph_vs_t *vs,
			 igraph_vector_t *v) {
  igraph_vector_t v2;
  IGRAPH_CHECK(igraph_vector_init_seq(&v2, vs->stdata[1], vs->stdata[2]));
  igraph_vector_destroy(v);
  *v=v2;
  return 0;
}

void igraph_vs_destroy_seq(igraph_vs_t *pvs) {
  igraph_vs_t *vs=(igraph_vs_t*)pvs;
  if (vs->shorthand) {
    Free(pvs);
  }
}

/* -------------------------------------------------- */
/* Vertex ids in a vector,                            */
/*   this can be a view or a copy                     */
/* -------------------------------------------------- */

void igraph_vs_next_vectorview(const igraph_t *graph, igraph_vs_t *vs);
bool_t igraph_vs_end_vectorview(const igraph_t *graph, const igraph_vs_t *vs);
void igraph_vs_reset_vectorview(const igraph_t *graph, igraph_vs_t *vs);
integer_t igraph_vs_get_vectorview(const igraph_t *graph, 
				   const igraph_vs_t *vs);
int igraph_vs_unfold_vectorview(const igraph_t *graph, const igraph_vs_t *vs,
				igraph_vector_t *v);
void igraph_vs_destroy_vectorview(igraph_vs_t *vs);

igraph_i_vstable_t igraph_i_vs_vectorview_table = {
  igraph_vs_next_vectorview, igraph_vs_end_vectorview, 
  igraph_vs_reset_vectorview, igraph_vs_get_vectorview, 
  igraph_vs_unfold_vectorview, igraph_vs_destroy_vectorview
};

typedef struct igraph_i_vs_vectorview_pdata_t {
  const igraph_vector_t v;
  bool_t destroy;
} igraph_i_vs_vectorview_pdata_t;

/**
 * \function igraph_vs_vectorview
 * 
 * This iterator type allows to handle a \type igraph_vector_t object
 * as a vertex sequence. Note that this function does \em not make a
 * copy of the original vector, so be sure that you don't destroy the
 * underlying \type igraph_vector_t object before destroying the
 * vertex set. 
 *
 * The \ref igraph_vs_vector_getvector() specific function can be
 * called on iterators created by this function.
 * 
 * \param igraph The underlying graph object.
 * \param vs The vertex sequence object.
 * \param vids The underlying \type igraph_vector_t object.
 * \return Error code, \c IGRAPH_ENOMEM if there is
 *   not enough memory.
 * 
 * \sa \ref igraph_vs_vector() is a similar iterator but uses a
 * private copy of the vector, \ref igraph_vs_vectorview_it().
 *
 * Time complexity: O(1).
 */

int igraph_vs_vectorview(const igraph_t *igraph, igraph_vs_t *vs, 
			 const igraph_vector_t *vids) {
  igraph_i_vs_vectorview_pdata_t *data;
  igraph_vector_t *fakev;

  vs->type=IGRAPH_ITERATOR_VS_VECTOR;
  vs->stdata[0]=0;
  vs->table=&igraph_i_vs_vectorview_table;
  vs->shorthand=0;

  vs->pdata=Calloc(1, igraph_i_vs_vectorview_pdata_t);
  if (vs->pdata==0) {
    IGRAPH_ERROR("Cannot create vector iterator", IGRAPH_ENOMEM);
  }
  data=(igraph_i_vs_vectorview_pdata_t*)vs->pdata;
  fakev=(igraph_vector_t*) &data->v;
  *fakev = *vids;
  data->destroy=0;

  vs->stdata[1]=igraph_vector_size(&data->v);

  return 0;
}

/**
 * \function igraph_vs_vector
 * 
 * This function creates a vertex sequence object from a
 * \type igraph_vector_t containing vertex ids. Unlike \ref
 * igraph_vs_vectorview() this function makes a copy of the original
 * vector, which can be safely destroyed.
 * 
 * The \ref igraph_vs_vector_getvector() specific function can be
 * called on iterators created by this function.
 * 
 * \param igraph The underlying graph.
 * \param vs The vertex sequence object.
 * \param vids The original vector.
 * \return Error code, \c IGRAPH_ENOMEM if we don't
 * have enough memory.
 *
 * Time complexity: O(n),
 * n is the number of elements in the
 * original vector.
 * 
 * \sa \ref igraph_vs_vectorview(), \ref igraph_vs_vector_small(), 
 * \ref igraph_vs_vectorview_it()
 */

int igraph_vs_vector(const igraph_t *igraph, igraph_vs_t *vs,
		     const igraph_vector_t *vids) {
  igraph_i_vs_vectorview_pdata_t *data;
  igraph_vector_t *fakev;

  vs->type=IGRAPH_ITERATOR_VS_VECTOR;
  vs->stdata[0]=0;
  vs->table=&igraph_i_vs_vectorview_table;
  vs->shorthand=0;
  
  vs->pdata=Calloc(1, igraph_i_vs_vectorview_pdata_t);
  if (vs->pdata==0) {
    IGRAPH_ERROR("Cannot create vector iterator", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, vs->pdata);
  data=(igraph_i_vs_vectorview_pdata_t*)vs->pdata;
  fakev=(igraph_vector_t*)&data->v;
  IGRAPH_CHECK(igraph_vector_copy(fakev, vids));
  data->destroy=1;

  vs->stdata[1]=igraph_vector_size(&data->v);
  
  IGRAPH_FINALLY_CLEAN(1);
  return 0;
}

/**
 * \function igraph_vs_vectorview_it
 * 
 * This function creates an iterator based on another iterator. 
 * The new iterator is similar to those created by \ref
 * igraph_vs_vectorview() and \ref igraph_vs_vector(), because 
 * the \ref igraph_vs_vector_getvector() specific function can be
 * called on it.
 * 
 * This iterator is mainly used to create a vector \em view of another
 * iterator. 
 * 
 * Note that the new iterator is a \em view of the original iterator,
 * don't destroy the original iterator before destroying the newly
 * created view.
 * \param graph The underlying graph object.
 * \param vs The original iterator. 
 * \param newvs Pointer to an uninitialized vertex sequence, this 
 *   will be initialized as a view of \p vs.
 * \return Error code, \c IGRAPH_ENOMEM if there is
 *   not enough memory.
 * 
 * Time complexity: this is the function of the type of the original
 * iterator. At most O(n), the number
 * of visited vertices by the original iterator. 
 * 
 * \sa \ref igraph_vs_vectorview(), \ref igraph_vs_vector().
 */

int igraph_vs_vectorview_it(const igraph_t *graph, const igraph_vs_t *vs,
			    igraph_vs_t *newvs) {
  igraph_i_vs_vectorview_pdata_t *data, *vsdata;
  igraph_vector_t *fakev;

  newvs->type=IGRAPH_ITERATOR_VS_VECTOR;
  newvs->stdata[0]=0;
  newvs->table=&igraph_i_vs_vectorview_table;
  newvs->shorthand=0;
  
  newvs->pdata=Calloc(1, igraph_i_vs_vectorview_pdata_t);
  if (newvs->pdata==0) {
    IGRAPH_ERROR("Cannot create vector iterator", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, newvs->pdata);
  data=(igraph_i_vs_vectorview_pdata_t*)newvs->pdata;
  fakev=(igraph_vector_t*)&data->v;
  if (vs->type != IGRAPH_ITERATOR_VS_VECTOR) {
    IGRAPH_VECTOR_INIT_FINALLY((igraph_vector_t*) &data->v, 0);
    IGRAPH_CHECK(igraph_vs_unfold(graph, vs, fakev));
    data->destroy=1;
    IGRAPH_FINALLY_CLEAN(1);
  } else {
    /* If it is a vector we just create a view */
    vsdata=(igraph_i_vs_vectorview_pdata_t*)vs->pdata;    
    *fakev=vsdata->v;
    data->destroy=0;
  }

  newvs->stdata[1]=igraph_vector_size(&data->v);
  
  IGRAPH_FINALLY_CLEAN(1);
  return 0;
}

/**
 * \function igraph_vs_vector_small
 * 
 * This is a convenience function for creating small vertex
 * sequences. The vertices are simply given as function parameters,
 * the end of the vertices are marked by a -1 parameter. 
 * 
 * \param igraph The underlying graph object.
 * \param vs The vertex sequence object.
 * \param ... The ids of the vertices are given by additional
 *   parameters closed by a -1 parameter.
 * \return Error code, \c IGRAPH_ENOMEM if there is
 *   not enough memory to create the vertex sequence.
 * 
 * \sa \ref igraph_vs_vectorview() and \ref igraph_vs_vector().
 * 
 * Time complexity: usually O(n),
 * the number of vertices in the vertex sequence.
 */

int igraph_vs_vector_small(const igraph_t *igraph, igraph_vs_t *vs, ...) {
  va_list ap;
  igraph_i_vs_vectorview_pdata_t *data;
  igraph_vector_t *fakev;
  long int i, n=0;
  vs->type=IGRAPH_ITERATOR_VS_VECTOR;
  vs->stdata[0]=0;
  vs->table=&igraph_i_vs_vectorview_table;
  vs->shorthand=0;
  
  vs->pdata=Calloc(1, igraph_i_vs_vectorview_pdata_t);
  if (vs->pdata==0) {
    IGRAPH_ERROR("Cannot create vector iterator", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, vs->pdata);
  data=(igraph_i_vs_vectorview_pdata_t*)vs->pdata;
  fakev=(igraph_vector_t*)&data->v;
  data->destroy=1;

  va_start(ap, vs);
  while (1) {
    int num = va_arg(ap, int);
    if (num == -1) {
      break;
    }
    n++;
  }
  va_end(ap);

  IGRAPH_VECTOR_INIT_FINALLY(fakev, n);
  
  va_start(ap, vs);
  for (i=0; i<n; i++) {
    VECTOR(*fakev)[i]=(real_t) va_arg(ap, int);
  }
  va_end(ap);  
  
  vs->stdata[1]=igraph_vector_size(&data->v);
  
  IGRAPH_FINALLY_CLEAN(2);
  return 0;
}  

/**
 * \function IGRAPH_VS_VECTOR
 * 
 * Vertex sequence shorthand for a \type igraph_vector_t object. 
 * \param graph The underlying graph object.
 * \param vids An initialized \type igraph_vector_t object containing
 * vertex ids.
 * \return Vertex sequence shorthand.
 * 
 * Time complexity: O(1).
 */

const igraph_vs_t *IGRAPH_VS_VECTOR(const igraph_t *graph, 
				    const igraph_vector_t *vids) {
  igraph_vs_t *vs=Calloc(1, igraph_vs_t);
  if (vs==0) {
    igraph_error("Cannot create iterator shorthand", __FILE__, __LINE__,
		 IGRAPH_ENOMEM);
    return 0;			/* TODO: how to sign error ??? */
  }
  igraph_vs_vectorview(graph, vs, vids);
  vs->shorthand=1;
  return vs;
}

/**
 * \function IGRAPH_VS
 * 
 * Vertex sequence shorthand for the vertices given as arguments.
 * \param graph The underlying graph object.
 * \param ... The vertex ids, the last argument has to be -1.
 * \return Vertex sequence shorthand.
 * 
 * Time complexity: O(n), the number
 * of vertices in the sequence.
 */

const igraph_vs_t *IGRAPH_VS(const igraph_t *graph, ...) {
  igraph_vs_t *vs=Calloc(1, igraph_vs_t);
  va_list ap;
  igraph_i_vs_vectorview_pdata_t *data;
  igraph_vector_t *fakev;
  long int i, n=0;
  int ret;
  if (vs==0) {
    igraph_error("Cannot create iterator shorthand", __FILE__, __LINE__,
		 IGRAPH_ENOMEM);
    return 0;			/* TODO: how to sign error ??? */
  }

  vs->type=IGRAPH_ITERATOR_VS_VECTOR;
  vs->stdata[0]=0;
  vs->table=&igraph_i_vs_vectorview_table;
  vs->shorthand=0;
  
  vs->pdata=Calloc(1, igraph_i_vs_vectorview_pdata_t);
  if (vs->pdata==0) {
    igraph_error("Cannot create vector iterator", __FILE__, __LINE__,
		 IGRAPH_ENOMEM);
    return 0;
  }
  IGRAPH_FINALLY(igraph_free, vs->pdata);
  data=(igraph_i_vs_vectorview_pdata_t*)vs->pdata;
  fakev=(igraph_vector_t*)&data->v;
  data->destroy=1;

  va_start(ap, graph);
  while (1) {
    int num = va_arg(ap, int);
    if (num == -1) {
      break;
    }
    n++;
  }
  va_end(ap);

  ret=igraph_vector_init(fakev, n);
  if (ret != 0) {
    igraph_error("Cannot create vector for iterator shorthand", __FILE__, 
		 __LINE__, ret);
    return 0;
  }
  IGRAPH_FINALLY(igraph_vector_destroy, fakev);
  
  va_start(ap, graph);
  for (i=0; i<n; i++) {
    VECTOR(*fakev)[i]=(real_t) va_arg(ap, int);
  }
  va_end(ap);  
  
  vs->stdata[1]=igraph_vector_size(&data->v);
  
  IGRAPH_FINALLY_CLEAN(2);
  vs->shorthand=1;
  return vs;  
}

void igraph_vs_next_vectorview(const igraph_t *graph, igraph_vs_t *vs) {
  vs->stdata[0] ++;
}

bool_t igraph_vs_end_vectorview(const igraph_t *graph, const igraph_vs_t *vs) {
  return vs->stdata[0] >= vs->stdata[1];
}

void igraph_vs_reset_vectorview(const igraph_t *graph, igraph_vs_t *vs) {
  vs->stdata[0]=0;
}

integer_t igraph_vs_get_vectorview(const igraph_t *graph, 
				   const igraph_vs_t *vs) {
  igraph_i_vs_vectorview_pdata_t *data=
    (igraph_i_vs_vectorview_pdata_t*)vs->pdata;
  return VECTOR(data->v)[ (long int) (vs->stdata[0]) ]; 
}

int igraph_vs_unfold_vectorview(const igraph_t *graph, const igraph_vs_t *vs,
				igraph_vector_t *v) {
  igraph_vector_t v2;
  igraph_i_vs_vectorview_pdata_t *data=
    (igraph_i_vs_vectorview_pdata_t*)vs->pdata;
  IGRAPH_CHECK(igraph_vector_copy(&v2, &data->v));
  igraph_vector_destroy(v);
  *v=v2;
  return 0;
}

void igraph_vs_destroy_vectorview(igraph_vs_t *pvs) {
  igraph_vs_t *vs=(igraph_vs_t*)pvs;
  igraph_i_vs_vectorview_pdata_t *data=
    (igraph_i_vs_vectorview_pdata_t*)vs->pdata;
  if (data->destroy) {
    igraph_vector_destroy((igraph_vector_t*)&data->v);
  }
  Free(data);
  
  if (vs->shorthand) {
    Free(pvs);
  }
}

/**
 * \function igraph_vs_vector_getvector
 * 
 * This is a specific vertex sequence function, it can be called for
 * vector type iterators only (these are created by \ref
 * igraph_vs_vector(), \ref igraph_vs_vectorview(), \ref
 * igraph_vs_vectorview_it() or \ref igraph_vs_vector_small()).
 * It gives access to the vertex ids in the sequence as a
 * \type igraph_vector_t type.
 * 
 * The result is undefined if you call it with a different iterator
 * type. 
 * 
 * \param graph The underlying graph object.
 * \param vs The vertex sequence. 
 * \return Pointer to a \type igraph_vector_t object. This object should
 *   be considered as constant, don't change its elements or size.
 * 
 * Time complexity: O(1).
 */

const igraph_vector_t *igraph_vs_vector_getvector(const igraph_t *graph, 
					   const igraph_vs_t *vs) {
  igraph_i_vs_vectorview_pdata_t *data=
    (igraph_i_vs_vectorview_pdata_t*)vs->pdata;
  return &data->v;
}


/* -------------------------------------------------- */
/* Edge iterator generics                             */
/* -------------------------------------------------- */

/**
 * \function igraph_es_next
 *
 * Steps the iterator to the next edge. Usually you want to check that
 * there is a next edge in the edge sequence by calling \ref
 * igraph_es_end(). 
 * 
 * \param graph The underlying graph object.
 * \param es The edge sequence object.
 * 
 * Time complexity: usually O(1), but
 * see also the specific iterator types.
 */

void igraph_es_next(const igraph_t *graph, igraph_es_t *es) {
  es->table->next(graph, es);
}

/**
 * \function igraph_es_end
 * 
 * Checks whether there are more elements in an edge sequence. 
 * \param graph The underlying graph object.
 * \param es The edge sequence object.
 * \return True (positive integer) if there are no more edges in the
 *   sequence, ie. the pointer points to the element \em after the
 *   last one. Fast otherwise.
 *
 * Time complexity: usually O(1), but
 * see also the specific iterator types.
 */

bool_t igraph_es_end(const igraph_t *graph, const igraph_es_t *es) {
  return es->table->end(graph, es);
}

/**
 * \function igraph_es_reset
 * 
 * This function resets the iterator position of an edge
 * sequence. This has the same effect as destroying and reinitializing
 * the iterator, only more efficient. 
 * \param graph The underlying graph object.
 * \param es The edge sequence object.
 * 
 * Time complexity: usually O(1), but
 * see also the specific iterator types.
 */

void igraph_es_reset(const igraph_t *graph, igraph_es_t *es) {
  return es->table->reset(graph, es);
}

/**
 * \function igraph_es_from
 * 
 * Return the vertex at the starting point of the current edge. For
 * undirected graphs one of the two endpoints is returned and the
 * other one is provided by \ref igraph_es_to(). Be sure to call \ref
 * igraph_es_end() to see whether the iterator really points to an
 * edge. 
 * \param graph The underlying graph object.
 * \param es The edge sequence object.
 * \return Id of the vertex.
 * 
 * Time complexity: usually O(1), but
 * see also the specific iterator types.
 * 
 * \sa \ref igraph_es_to().
 */

integer_t igraph_es_from(const igraph_t *graph, const igraph_es_t *es) {
  return es->table->from(graph, es);
}

/**
 * \function igraph_es_to
 * 
 * Return the vertex at the end point of the current edge. For
 * undirected graphs one of the two endpoints is returned and the
 * other one is provided by \ref igraph_es_from(). Be sure to call \ref
 * igraph_es_end() to see whether the iterator really points to an
 * edge. 
 * \param graph The underlying graph object.
 * \param es The edge sequence object.
 * \return Id of the vertex.
 * 
 * Time complexity: usually O(1), but
 * see also the specific iterator types.
 * 
 * \sa \ref igraph_es_from().
 */

integer_t igraph_es_to(const igraph_t *graph, const igraph_es_t *es) {
  return es->table->to(graph, es);
}

/**
 * \function igraph_es_get
 * 
 * This function returns the id of the current edge an iterator. (Or
 * in other words the edge the iterator \quote points\endquote to.
 * You should call \ref igraph_es_end() on the iterator beforehand to
 * be sure that the iterator really points to an edge.
 * \param graph The underlying graph object.
 * \param es The vertex sequence object.
 * \return Id of the \quote current\endquote edge.
 * 
 * Time complexity: usually O(1) 
 * but see also the specific types.
 */

integer_t igraph_es_get(const igraph_t *graph, const igraph_es_t *es) {
  return es->table->get(graph, es);
}

/**
 * \function igraph_es_unfold
 * 
 * Steps over an edge sequence and stores the visited edges in a
 * vector. 
 * \param graph The underlying graph object.
 * \param es The edge sequence object.
 * \param v An initialized \type igraph_vector_t object. These will be
 *   resized to hold the edges.
 * \return Error code, \c IGRAPH_ENOMEM if there is
 *   not enough memory to resize \p v.
 * 
 * Time complexity: depends on the type of the iterator, but usually
 * it is at most O(n), the number of
 * visited vertices.
 */

int igraph_es_unfold(const igraph_t *graph, const igraph_es_t *es,
		     igraph_vector_t *v) {
  return es->table->unfold(graph, es, v);
}

/**
 * \function igraph_es_destroy
 *
 * Destroys an edge sequence, ie. frees all memory allocated for it. 
 * The iterator has to be reinitialized before using it again.
 * \param es The edge sequence object to destroy.
 * 
 * Time complexity: usually O(1) 
 * but see also the specific types.
 */

void igraph_es_destroy(igraph_es_t *es) {
  es->table->destroy(es);
}

/* -------------------------------------------------- */
/* Simple edge iterator                               */
/* -------------------------------------------------- */

void igraph_es_next_all(const igraph_t *graph, igraph_es_t *es);
bool_t igraph_es_end_all(const igraph_t *graph, const igraph_es_t *es);
void igraph_es_reset_all(const igraph_t *graph, igraph_es_t *es);
integer_t igraph_es_get_all(const igraph_t *graph, const igraph_es_t *es);
integer_t igraph_es_from_all(const igraph_t *graph, const igraph_es_t *es);
integer_t igraph_es_to_all(const igraph_t *graph, const igraph_es_t *es);
int igraph_es_unfold_all(const igraph_t *graph, const igraph_es_t *es,
			 igraph_vector_t *v);
void igraph_es_destroy_all(igraph_es_t *es);

igraph_i_estable_t igraph_i_es_all_table = {
  igraph_es_next_all, igraph_es_end_all, igraph_es_reset_all,
  igraph_es_get_all, igraph_es_from_all, igraph_es_to_all,
  igraph_es_unfold_all, igraph_es_destroy_all
};

/**
 * \function igraph_es_all 
 * 
 * Edge sequence containing all edges in arbitrary 
 * order. This iterator doesn't have any special functions.
 * \param graph The underlying graph object.
 * \param es Pointer to an uninitialized edge set object.
 * \return Error code, the current implementation always returns with
 *   success. 
 * 
 * Time complexity: O(1).
 * 
 * \sa \ref igraph_es_fromorder().
 */

int igraph_es_all(const igraph_t *graph, igraph_es_t *es) {
  es->type=IGRAPH_ITERATOR_ES_ALL;
  es->stdata[0]=0;
  es->table=&igraph_i_es_all_table;
  es->shorthand=0;
  return 0;
}

/**
 * \function IGRAPH_ES_ALL
 * 
 * Edge sequence shorthand for all edges in a graph.
 * \param graph The underlying graph object.
 * \return Edge sequence shorthand.
 *
 * Time complexity: O(1).
 */

const igraph_es_t *IGRAPH_ES_ALL(const igraph_t *graph) {
  igraph_es_t *es=Calloc(1, igraph_es_t);
  if (es==0) {
    igraph_error("Cannot create iterator shorthand", __FILE__, __LINE__,
		 IGRAPH_ENOMEM);
    return 0;			/* TODO: how to sign error ??? */
  }
  igraph_es_all(graph, es);
  es->shorthand=1;
  return es;
}

void igraph_es_next_all(const igraph_t *graph, igraph_es_t *es) {
  es->stdata[0] ++;
}

bool_t igraph_es_end_all(const igraph_t *graph, const igraph_es_t *es) {
  return es->stdata[0] >= igraph_ecount(graph);
}

void igraph_es_reset_all(const igraph_t *graph, igraph_es_t *es) {
  es->stdata[0]=0;
}

integer_t igraph_es_get_all(const igraph_t *graph, const igraph_es_t *es) {
  return es->stdata[0];
}

integer_t igraph_es_from_all(const igraph_t *graph, const igraph_es_t *es) {
  return VECTOR(graph->from)[ (long int) es->stdata[0] ];
}

integer_t igraph_es_to_all(const igraph_t *graph, const igraph_es_t *es) {
  return VECTOR(graph->to)[ (long int) es->stdata[0] ];
}

int igraph_es_unfold_all(const igraph_t *graph, const igraph_es_t *es, 
			 igraph_vector_t *v) {
  long int n;
  igraph_vector_t v2;
  n=igraph_ecount(graph);
  IGRAPH_CHECK(igraph_vector_init_seq(&v2, 0, n-1));
  igraph_vector_destroy(v);
  *v=v2;
  return 0;
}

void igraph_es_destroy_all(igraph_es_t *pvs) {
  igraph_es_t *es=(igraph_es_t*)pvs;
  if (es->shorthand) {
    Free(es);
  }
}

/* -------------------------------------------------- */
/* Ordered edge iterator                              */
/* -------------------------------------------------- */

void igraph_es_next_fromorder(const igraph_t *graph, igraph_es_t *es);
bool_t igraph_es_end_fromorder(const igraph_t *graph, const igraph_es_t *es);
void igraph_es_reset_fromorder(const igraph_t *graph, igraph_es_t *es);
integer_t igraph_es_get_fromorder(const igraph_t *graph, 
				  const igraph_es_t *es);
integer_t igraph_es_from_fromorder(const igraph_t *graph, 
				   const igraph_es_t *es);
integer_t igraph_es_to_fromorder(const igraph_t *graph, const igraph_es_t *es);
int igraph_es_unfold_fromorder(const igraph_t *graph, const igraph_es_t *es,
			       igraph_vector_t *v);
void igraph_es_destroy_fromorder(igraph_es_t *es);

igraph_i_estable_t igraph_i_es_fromorder_table = {
  igraph_es_next_fromorder, igraph_es_end_fromorder, igraph_es_reset_fromorder,
  igraph_es_get_fromorder, igraph_es_from_fromorder, igraph_es_to_fromorder,
  igraph_es_unfold_fromorder, igraph_es_destroy_fromorder
};
  
/**
 * \function igraph_es_fromorder
 * 
 * Edge sequence containing all edges in the order of increasing
 * starting vertex id. Ie. first come all edges \em from vertex 0 in
 * arbitrary order, then the ones \em from vertex 1 in arbitrary
 * order, etc. For undirected graph the edges are visited in arbitrary
 * order. 
 * \param graph The underlying graph object.
 * \param es Pointer to an uninitialized edge set object.
 * \return Error code, the current implementation always returns with
 *   success. 
 * 
 * Time complexity: O(1).
 * 
 * \sa \ref igraph_es_fromorder().
 */

int igraph_es_fromorder(const igraph_t *graph, igraph_es_t *es) {
  es->type=IGRAPH_ITERATOR_ES_FROMORDER;
  es->stdata[0]=0;
  es->table=&igraph_i_es_fromorder_table;
  es->shorthand=0;
  return 0;
}

void igraph_es_next_fromorder(const igraph_t *graph, igraph_es_t *es) {
  es->stdata[0] ++;
}

bool_t igraph_es_end_fromorder(const igraph_t *graph, const igraph_es_t *es) {
  return es->stdata[0] >= igraph_ecount(graph);
}

void igraph_es_reset_fromorder(const igraph_t *graph, igraph_es_t *es) {
  es->stdata[0]=0;
}

integer_t igraph_es_get_fromorder(const igraph_t *graph, 
				  const igraph_es_t *es) {
  return VECTOR(graph->oi)[ (long int) es->stdata[0] ];
}

integer_t igraph_es_from_fromorder(const igraph_t *graph, 
				   const igraph_es_t *es) {
  long int idx=VECTOR(graph->oi)[ (long int) es->stdata[0] ];
  return VECTOR(graph->from)[ idx ];
}

integer_t igraph_es_to_fromorder(const igraph_t *graph, 
				 const igraph_es_t *es) {
  long int idx=VECTOR(graph->oi)[ (long int) es->stdata[0] ];
  return VECTOR(graph->to)[ idx ];
}

int igraph_es_unfold_fromorder(const igraph_t *graph, const igraph_es_t *es,
			       igraph_vector_t *v) {
  long int i=0, n=igraph_ecount(graph);  
  IGRAPH_CHECK(igraph_vector_resize(v, n));
  for (i=0; i<n; i++) {
    VECTOR(*v)[i]=VECTOR(graph->oi)[i];
  }
  return 0;
}

void igraph_es_destroy_fromorder(igraph_es_t *pvs) {
  igraph_es_t *es=(igraph_es_t*)pvs;
  if (es->shorthand) {
    Free(es);
  }
}

/* -------------------------------------------------- */
/* Adjacent edges of a vertex                         */
/* -------------------------------------------------- */

void igraph_es_next_adj(const igraph_t *graph, igraph_es_t *es);
bool_t igraph_es_end_adj(const igraph_t *graph, const igraph_es_t *es);
void igraph_es_reset_adj(const igraph_t *graph, igraph_es_t *es);
integer_t igraph_es_get_adj(const igraph_t *graph, const igraph_es_t *es);
integer_t igraph_es_from_adj(const igraph_t *graph, const igraph_es_t *es);
integer_t igraph_es_to_adj(const igraph_t *graph, const igraph_es_t *es);
int igraph_es_unfold_adj(const igraph_t *graph, const igraph_es_t *es,
			 igraph_vector_t *v);
void igraph_es_destroy_adj(igraph_es_t *es);

igraph_i_estable_t igraph_i_es_adj_table = {
  igraph_es_next_adj, igraph_es_end_adj, igraph_es_reset_adj,
  igraph_es_get_adj, igraph_es_from_adj, igraph_es_to_adj,
  igraph_es_unfold_adj, igraph_es_destroy_adj
};
  
/**
 * \function igraph_es_adj
 * 
 * A sequence containing all adjacenct edges of a vertex.
 * 
 * \param graph The underlying graph object.
 * \param es The edge sequence object. 
 * \param vid The vertex of which the adjacent edges will be visited. 
 * \param mode Constant, specifies the type of adjacenct edges to
 *   visit. \c IGRAPH_OUT visits only outgoing,
 *   \c IGRAPH_IN only incoming edges,
 *   \c IGRAPH_ALL visits all adjacent edges.
 * \return Error code, the current implementation always returns with
 *   success. 
 * 
 * Time complexity: O(1).
 */

int igraph_es_adj(const igraph_t *graph, igraph_es_t *es,
		  integer_t vid, igraph_neimode_t mode) {
  if (!igraph_is_directed(graph)) {
    mode=IGRAPH_ALL;
  }
  
  es->type=IGRAPH_ITERATOR_ES_ADJ;
  es->shorthand=0;

  es->stdata[0]=vid;
  es->stdata[1]=mode;
  if (mode & IGRAPH_OUT) {
    es->stdata[2]=VECTOR(graph->os)[(long int)vid];
  } else {
    es->stdata[2]=igraph_ecount(graph);
  }
  if (mode & IGRAPH_IN) {
    es->stdata[3]=VECTOR(graph->is)[(long int)vid];
  } else {
    es->stdata[3]=igraph_ecount(graph);
  }

  es->table=&igraph_i_es_adj_table;
  return 0;
}

void igraph_es_next_adj(const igraph_t *graph, igraph_es_t *es) {
  es->stdata[2] ++;
  if (es->stdata[2] > VECTOR(graph->os)[ (long int)es->stdata[0]+1 ]) {
    es->stdata[3] ++;
  }
}

bool_t igraph_es_end_adj(const igraph_t *graph, const igraph_es_t *es) {
  return (es->stdata[2] >= VECTOR(graph->os)[ (long int) es->stdata[0]+1 ] &&
	  es->stdata[3] >= VECTOR(graph->is)[ (long int) es->stdata[0]+1 ]);
}

void igraph_es_reset_adj(const igraph_t *graph, igraph_es_t *es) {
  if ((int)es->stdata[1] & IGRAPH_OUT) {
    es->stdata[2]=VECTOR(graph->os)[(long int)es->stdata[0]];
  } else {
    es->stdata[2]=igraph_ecount(graph);
  }
  if ((int)es->stdata[1] & IGRAPH_IN) {
    es->stdata[3]=VECTOR(graph->is)[(long int)es->stdata[0]];
  } else {
    es->stdata[3]=igraph_ecount(graph);
  }
}

integer_t igraph_es_get_adj(const igraph_t *graph, const igraph_es_t *es) {
  if (es->stdata[2] < VECTOR(graph->os)[ (long int)es->stdata[0]+1 ]) {
    return VECTOR(graph->oi)[(long int)es->stdata[2]];
  } else {
    return VECTOR(graph->ii)[(long int)es->stdata[3]];
  }
}

integer_t igraph_es_from_adj(const igraph_t *graph, const igraph_es_t *es) {
  if (es->stdata[2] < VECTOR(graph->os)[ (long int)es->stdata[0]+1 ]) {
    long int idx=VECTOR(graph->oi)[(long int)es->stdata[2]];
    return VECTOR(graph->from)[idx];
  } else {
    long int idx=VECTOR(graph->ii)[(long int)es->stdata[3]];
    return VECTOR(graph->from)[idx];
  }
}

integer_t igraph_es_to_adj(const igraph_t *graph, const igraph_es_t *es) {
  if (es->stdata[2] < VECTOR(graph->os)[ (long int)es->stdata[0]+1 ]) {
    long int idx=VECTOR(graph->oi)[(long int)es->stdata[2]];
    return VECTOR(graph->to)[idx];
  } else {
    long int idx=VECTOR(graph->ii)[(long int)es->stdata[3]];
    return VECTOR(graph->to)[idx];
  }
}

int igraph_es_unfold_adj(const igraph_t *graph, const igraph_es_t *es,
			 igraph_vector_t *v) {
  long int length=0, idx=0;
  long int no_of_edges;
  long int i;
  
  long int node=es->stdata[0];
  int mode=es->stdata[1];
  
  no_of_edges=igraph_ecount(graph);

  /* Calculate needed space first & allocate it*/
  
  if (mode & IGRAPH_OUT) {
    length += (VECTOR(graph->os)[node+1] - VECTOR(graph->os)[node]);
  }
  if (mode & IGRAPH_IN) {
    length += (VECTOR(graph->is)[node+1] - VECTOR(graph->is)[node]);
  }
  
  IGRAPH_CHECK(igraph_vector_resize(v, length));
  
  if (mode & IGRAPH_OUT) {
    for (i=VECTOR(graph->os)[node]; i<VECTOR(graph->os)[node+1]; i++) {
      VECTOR(*v)[idx++] = VECTOR(graph->oi)[i];
    }
  }
  if (mode & IGRAPH_IN) {
    for (i=VECTOR(graph->is)[node]; i<VECTOR(graph->is)[node+1]; i++) {
      VECTOR(*v)[idx++] = VECTOR(graph->ii)[i];
    }
  }

  return 0;
}

void igraph_es_destroy_adj(igraph_es_t *pvs) {
  igraph_es_t *es=(igraph_es_t*)pvs;
  if (es->shorthand) {
    Free(es);
  }
}

/**
 * \function igraph_es_adj_set
 * 
 * Reinitialize an iterator created by \ref igraph_es_adj() for a
 * different vertex. This is the same as destroying the iterator and
 * initializing it with a different \p vid and/or
 * \p mode argument.
 * 
 * The arguments are the same as for the \ref igraph_es_adj()
 * function. 
 */

void igraph_es_adj_set(const igraph_t *graph, igraph_es_t *es,
		       integer_t vid, igraph_neimode_t mode) {
  if (!igraph_is_directed(graph)) {
    mode=IGRAPH_ALL;
  }
  es->stdata[0]=vid;
  es->stdata[1]=mode;
  if (mode & IGRAPH_OUT) {
    es->stdata[2]=VECTOR(graph->os)[(long int)vid];
  } else {
    es->stdata[2]=igraph_ecount(graph);
  }
  if (mode & IGRAPH_IN) {
    es->stdata[3]=VECTOR(graph->is)[(long int)vid];
  } else {
    es->stdata[3]=igraph_ecount(graph);
  }
}

/**
 * \function igraph_es_adj_vertex
 * 
 * Special function for \ref igraph_es_adj() iterators. Provides the
 * id of the vertex at the \em other end of the edge. Ie. if the
 * iterator visits the adjacent edges of vertex
 * \p vid then this function always returns the id
 * of the other vertex, not \p vid.
 * 
 * Time complexity: O(1).
 */

integer_t igraph_es_adj_vertex(const igraph_t *graph, const igraph_es_t *es) {
  if (es->stdata[2] < VECTOR(graph->os)[ (long int)es->stdata[0]+1 ]) {
    long int idx=VECTOR(graph->oi)[(long int)es->stdata[2]];
    return VECTOR(graph->to)[idx];
  } else {
    long int idx=VECTOR(graph->ii)[(long int)es->stdata[3]];
    return VECTOR(graph->from)[idx];
  }
}

/* -------------------------------------------------- */
/* Empty edge iterator                                */
/* -------------------------------------------------- */

void igraph_es_next_none(const igraph_t *graph, igraph_es_t *es);
bool_t igraph_es_end_none(const igraph_t *graph, const igraph_es_t *es);
void igraph_es_reset_none(const igraph_t *graph, igraph_es_t *es);
integer_t igraph_es_get_none(const igraph_t *graph, const igraph_es_t *es);
integer_t igraph_es_from_none(const igraph_t *graph, const igraph_es_t *es);
integer_t igraph_es_to_none(const igraph_t *graph, const igraph_es_t *es);
int igraph_es_unfold_none(const igraph_t *graph, const igraph_es_t *es, 
			  igraph_vector_t *v);
void igraph_es_destroy_none(igraph_es_t *es);

igraph_i_estable_t igraph_i_es_none_table = {
  igraph_es_next_none, igraph_es_end_none, igraph_es_reset_none,
  igraph_es_get_none, igraph_es_from_none, igraph_es_to_none,
  igraph_es_unfold_none, igraph_es_destroy_none
};

/**
 * \function igraph_es_none
 * 
 * Creates an empty edge sequence. 
 *
 * \param graph The underlying graph object.
 * \param es The edge sequence object.
 * \return Error code, the current implementation always returned with
 *   success. 
 * 
 * Time complexity: O(1).
 */

int igraph_es_none(const igraph_t *graph, igraph_es_t *es) {
  es->type=IGRAPH_ITERATOR_ES_NONE;
  es->shorthand=0;  
  es->table=&igraph_i_es_none_table;
  return 0;
}

void igraph_es_next_none(const igraph_t *graph, igraph_es_t *es) {
  /* nothing to do */
}

bool_t igraph_es_end_none(const igraph_t *graph, const igraph_es_t *es) {
  return 1;
}

void igraph_es_reset_none(const igraph_t *graph, igraph_es_t *es) {
  /* nothing to do */
}

integer_t igraph_es_get_none(const igraph_t *graph, const igraph_es_t *es) {
  /* ooops this is an error, no way to signal it though... */
  return -1;
}

integer_t igraph_es_from_none(const igraph_t *graph, const igraph_es_t *es) {
  /* error */
  return 0;
}

integer_t igraph_es_to_none(const igraph_t *graph, const igraph_es_t *es) {
  /* error */
  return 0;
}

int igraph_es_unfold_none(const igraph_t *graph, const igraph_es_t *es,
			  igraph_vector_t *v) {
  igraph_vector_clear(v);
  return 0;
}

void igraph_es_destroy_none(igraph_es_t *pvs) {
  igraph_es_t *es=(igraph_es_t*)pvs;
  if (es->shorthand) {
    Free(es);
  }
}

/* -------------------------------------------------- */
/* edge iterator, single edge                         */
/* -------------------------------------------------- */

void igraph_es_next_1(const igraph_t *graph, igraph_es_t *es);
bool_t igraph_es_end_1(const igraph_t *graph, const igraph_es_t *es);
void igraph_es_reset_1(const igraph_t *graph, igraph_es_t *es);
integer_t igraph_es_get_1(const igraph_t *graph, const igraph_es_t *es);
integer_t igraph_es_from_1(const igraph_t *graph, const igraph_es_t *es);
integer_t igraph_es_to_1(const igraph_t *graph, const igraph_es_t *es);
int igraph_es_unfold_1(const igraph_t *graph, const igraph_es_t *es, 
		       igraph_vector_t *v);
void igraph_es_destroy_1(igraph_es_t *es);

igraph_i_estable_t igraph_i_es_1_table = {
  igraph_es_next_1, igraph_es_end_1, igraph_es_reset_1,
  igraph_es_get_1, igraph_es_from_1, igraph_es_to_1,
  igraph_es_unfold_1, igraph_es_destroy_1
};

/**
 * \function igraph_es_1
 * 
 * Creates an edge sequence containing a single edge. 
 * 
 * \param igraph The underlying graph object.
 * \param es The edge sequence object.
 * \param eid The single edge id to be contained in the edge
 * sequence.
 * \return Error code, the current implementation always returns with
 * success. 
 */

int igraph_es_1(const igraph_t *igraph, igraph_es_t *es, integer_t eid) {
  es->type=IGRAPH_ITERATOR_ES_1;
  es->stdata[0]=eid;
  es->stdata[1]=0;		/* write 1 here if end */
  es->table=&igraph_i_es_1_table;
  es->shorthand=0;
  return 0;
}

/**
 * \function IGRAPH_ES_1
 * 
 * Edge sequence shorthand for a single edge.
 * \param graph The underlying graph object.
 * \param eid The id of the single edge in the sequence.
 * \return Edge sequence shorthand.
 * 
 * Time complexity: O(1).
 */

const igraph_es_t *IGRAPH_ES_1(const igraph_t *graph, integer_t eid) {
  igraph_es_t *es=Calloc(1, igraph_es_t);
  if (es==0) {
    igraph_error("Cannot create iterator shorthand", __FILE__, __LINE__,
		 IGRAPH_ENOMEM);
    return 0;			/* TODO: how to sign error ??? */
  }
  igraph_es_1(graph, es, eid);
  es->shorthand=1;
  return es;
}

void igraph_es_next_1(const igraph_t *graph, igraph_es_t *es) {
  /* signal end */
  es->stdata[1]=1;
}

bool_t igraph_es_end_1(const igraph_t *graph, const igraph_es_t *es) {
  return (es->stdata[1]==1);
}

void igraph_es_reset_1(const igraph_t *graph, igraph_es_t *es) {
  es->stdata[1]=0;
}

integer_t igraph_es_get_1(const igraph_t *graph, const igraph_es_t *es) {
  return es->stdata[0];
}

integer_t igraph_es_from_1(const igraph_t *graph, const igraph_es_t *es) {
  return VECTOR(graph->from)[ (long int) es->stdata[0] ];
}

integer_t igraph_es_to_1(const igraph_t *graph, const igraph_es_t *es) {
  return VECTOR(graph->to)[ (long int) es->stdata[0] ];
}

int igraph_es_unfold_1(const igraph_t *graph, const igraph_es_t *es,
		       igraph_vector_t *v) {
  IGRAPH_CHECK(igraph_vector_resize(v, 1));
  VECTOR(*v)[0]=es->stdata[0];
  return 0;
}

void igraph_es_destroy_1(igraph_es_t *pvs) {
  igraph_es_t *es=(igraph_es_t*)pvs;
  if (es->shorthand) {
    Free(es);
  }
}

/* -------------------------------------------------- */
/* Edge set with sequence of vertices               */
/* -------------------------------------------------- */

void igraph_es_next_seq(const igraph_t *graph, igraph_es_t *es);
bool_t igraph_es_end_seq(const igraph_t *graph, const igraph_es_t *es);
void igraph_es_reset_seq(const igraph_t *graph, igraph_es_t *es);
integer_t igraph_es_get_seq(const igraph_t *graph, const igraph_es_t *es);
integer_t igraph_es_from_seq(const igraph_t *graph, const igraph_es_t *es);
integer_t igraph_es_to_seq(const igraph_t *graph, const igraph_es_t *es);
int igraph_es_unfold_seq(const igraph_t *graph, const igraph_es_t *es,
			 igraph_vector_t *v);
void igraph_es_destroy_seq(igraph_es_t *es);

igraph_i_estable_t igraph_i_es_seq_table = {
  igraph_es_next_seq, igraph_es_end_seq, igraph_es_reset_seq,
  igraph_es_get_seq, igraph_es_from_seq, igraph_es_to_seq,
  igraph_es_unfold_seq, igraph_es_destroy_seq
};

/**
 * \function igraph_es_seq
 * 
 * Creates a regular sequence of edges. 
 * \param igraph The underlying graph object. 
 * \param es Pointer to an uninitialized edge sequence object.
 * \param from The lower limit of the inverval (inclusive).
 * \param to The upper limit of the interval (inclusive).
 * \return Error code, the current implementation always returns with
 *   success. 
 * 
 * Time complexity: O(1).
 */

int igraph_es_seq(const igraph_t *igraph, igraph_es_t *es, integer_t from,
		  integer_t to) {
  es->type=IGRAPH_ITERATOR_ES_SEQ;
  es->stdata[0]=from;
  es->stdata[1]=from;
  es->stdata[2]=to;
  es->table=&igraph_i_es_seq_table;
  es->shorthand=0;
  return 0;
}

void igraph_es_next_seq(const igraph_t *graph, igraph_es_t *es) {
  es->stdata[0] ++;
}

bool_t igraph_es_end_seq(const igraph_t *graph, const igraph_es_t *es) {
  return es->stdata[0] > es->stdata[2];
}

void igraph_es_reset_seq(const igraph_t *graph, igraph_es_t *es) {
  es->stdata[0]=es->stdata[1];
}

integer_t igraph_es_get_seq(const igraph_t *graph, const igraph_es_t *es) {
  return es->stdata[0];
}

integer_t igraph_es_from_seq(const igraph_t *graph, const igraph_es_t *es) {
  return VECTOR(graph->from)[ (long int) es->stdata[0] ];
}

integer_t igraph_es_to_seq(const igraph_t *graph, const igraph_es_t *es) {
  return VECTOR(graph->to)[ (long int) es->stdata[0] ];
}

int igraph_es_unfold_seq(const igraph_t *graph, const igraph_es_t *es,
			 igraph_vector_t *v) {
  igraph_vector_t v2;
  IGRAPH_CHECK(igraph_vector_init_seq(&v2, es->stdata[1], es->stdata[2]));
  igraph_vector_destroy(v);
  *v=v2;
  return 0;
}

void igraph_es_destroy_seq(igraph_es_t *pes) {
  igraph_es_t *es=(igraph_es_t*)pes;
  if (es->shorthand) {
    Free(pes);
  }
}

/* -------------------------------------------------- */
/* Edge ids in a vector, this is a view               */
/* -------------------------------------------------- */

void igraph_es_next_vectorview(const igraph_t *graph, igraph_es_t *es);
bool_t igraph_es_end_vectorview(const igraph_t *graph, const igraph_es_t *es);
void igraph_es_reset_vectorview(const igraph_t *graph, igraph_es_t *es);
integer_t igraph_es_get_vectorview(const igraph_t *graph, 
				   const igraph_es_t *es);
integer_t igraph_es_from_vectorview(const igraph_t *graph, 
				    const igraph_es_t *es);
integer_t igraph_es_to_vectorview(const igraph_t *graph, 
				  const igraph_es_t *es);
int igraph_es_unfold_vectorview(const igraph_t *graph, const igraph_es_t *es,
				igraph_vector_t *v);
void igraph_es_destroy_vectorview(igraph_es_t *es);

igraph_i_estable_t igraph_i_es_vectorview_table = {
  igraph_es_next_vectorview, igraph_es_end_vectorview, 
  igraph_es_reset_vectorview, igraph_es_get_vectorview, 
  igraph_es_from_vectorview, igraph_es_to_vectorview,
  igraph_es_unfold_vectorview, igraph_es_destroy_vectorview
};

typedef igraph_i_vs_vectorview_pdata_t igraph_i_es_vectorview_pdata_t;

/**
 * \function igraph_es_vectorview
 * 
 * This iterator type allows to handle a \type igraph_vector_t object
 * as an edge sequence. Note that this function does \em not make a
 * copy of the original vector, so be sure that you don't destroy the
 * underlying \type igraph_vector_t object before destroying the
 * edge set. 
 *
 * The \ref igraph_es_vector_getvector() specific function can be
 * called on iterators created by this function.
 * 
 * \param igraph The underlying graph object.
 * \param es The edge sequence object.
 * \param eids The underlying \type igraph_vector_t object.
 * \return Error code, \c IGRAPH_ENOMEM if there is
 *   not enough memory.
 * 
 * \sa \ref igraph_es_vector() is a similar iterator but uses a
 * private copy of the vector, \ref igraph_es_vectorview_it().
 *
 * Time complexity: O(1).
 */

int igraph_es_vectorview(const igraph_t *igraph, igraph_es_t *es, 
			 const igraph_vector_t *eids) {
  igraph_i_es_vectorview_pdata_t *data;
  igraph_vector_t *fakev;

  es->type=IGRAPH_ITERATOR_ES_VECTOR;
  es->stdata[0]=0;
  es->table=&igraph_i_es_vectorview_table;
  es->shorthand=0;
  
  es->pdata=Calloc(1, igraph_i_es_vectorview_pdata_t);
  if (es->pdata==0) {
    IGRAPH_ERROR("Cannot create vector iterator", IGRAPH_ENOMEM);
  }
  data=(igraph_i_es_vectorview_pdata_t*)es->pdata;
  fakev=(igraph_vector_t*) &data->v;
  *fakev = *eids;
  data->destroy=0;

  es->stdata[1]=igraph_vector_size(&data->v);

  return 0;
}

/**
 * \function igraph_es_vector
 * 
 * This function creates a edge sequence object from a
 * \type igraph_vector_t containing edge ids. Unlike \ref
 * igraph_es_vectorview() this function makes a copy of the original
 * vector, which can be safely destroyed.
 * 
 * The \ref igraph_es_vector_getvector() specific function can be
 * called on iterators created by this function.
 * 
 * \param igraph The underlying graph.
 * \param es The edge sequence object.
 * \param eids The original vector.
 * \return Error code, \c IGRAPH_ENOMEM if we don't
 * have enough memory.
 *
 * Time complexity: O(n),
 * n is the number of elements in the
 * original vector.
 * 
 * \sa \ref igraph_es_vectorview(), \ref igraph_es_vector_small(), 
 * \ref igraph_es_vectorview_it()
 */

int igraph_es_vector(const igraph_t *igraph, igraph_es_t *es,
		     const igraph_vector_t *eids) {
  igraph_i_es_vectorview_pdata_t *data;
  igraph_vector_t *fakev;

  es->type=IGRAPH_ITERATOR_ES_VECTOR;
  es->stdata[0]=0;
  es->table=&igraph_i_es_vectorview_table;
  es->shorthand=0;
  
  es->pdata=Calloc(1, igraph_i_es_vectorview_pdata_t);
  if (es->pdata==0) {
    IGRAPH_ERROR("Cannot create vector iterator", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, es->pdata);
  data=(igraph_i_es_vectorview_pdata_t*)es->pdata;
  fakev=(igraph_vector_t*)&data->v;
  IGRAPH_CHECK(igraph_vector_copy(fakev, eids));
  data->destroy=1;

  es->stdata[1]=igraph_vector_size(&data->v);
  
  IGRAPH_FINALLY_CLEAN(1);
  return 0;
}

/**
 * \function igraph_es_vectorview_it
 * 
 * This function creates an iterator based on another iterator. 
 * The new iterator is similar to those created by \ref
 * igraph_es_vectorview() and \ref igraph_es_vector(), because 
 * the \ref igraph_es_vector_getvector() specific function can be
 * called on it.
 * 
 * This iterator is mainly used to create a vector \em view of another
 * iterator. 
 * 
 * Note that the new iterator is a \em view of the original iterator,
 * don't destroy the original iterator before destroying the newly
 * created view.
 * \param graph The underlying graph object.
 * \param es The original iterator. 
 * \param newes Pointer to an uninitialized edge sequence, this 
 *   will be initialized as a view of \p es.
 * \return Error code, \c IGRAPH_ENOMEM if there is
 *   not enough memory.
 * 
 * Time complexity: this is the function of the type of the original
 * iterator. At most O(n), the number
 * of visited edges by the original iterator. 
 * 
 * \sa \ref igraph_es_vectorview(), \ref igraph_es_vector().
 */

int igraph_es_vectorview_it(const igraph_t *graph, const igraph_es_t *es,
			    igraph_es_t *newes) {
  igraph_i_es_vectorview_pdata_t *data, *esdata;
  igraph_vector_t *fakev;

  newes->type=IGRAPH_ITERATOR_ES_VECTOR;
  newes->stdata[0]=0;
  newes->table=&igraph_i_es_vectorview_table;
  newes->shorthand=0;
  
  newes->pdata=Calloc(1, igraph_i_es_vectorview_pdata_t);
  if (newes->pdata==0) {
    IGRAPH_ERROR("Cannot create vector iterator", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, newes->pdata);
  data=(igraph_i_es_vectorview_pdata_t*)newes->pdata;
  IGRAPH_VECTOR_INIT_FINALLY((igraph_vector_t*) &data->v, 0);
  fakev=(igraph_vector_t*)&data->v;
  if (es->type != IGRAPH_ITERATOR_ES_VECTOR) {
    igraph_vector_init(fakev, 0);
    IGRAPH_CHECK(igraph_es_unfold(graph, es, fakev));
    data->destroy=1;
  } else {
    /* If it is a vector we just create a view */
    esdata=(igraph_i_es_vectorview_pdata_t*)es->pdata;    
    *fakev=esdata->v;
    data->destroy=0;
  }

  newes->stdata[1]=igraph_vector_size(&data->v);
  
  IGRAPH_FINALLY_CLEAN(2);
  return 0;
}

/**
 * \function igraph_es_vector_small
 * 
 * This is a convenience function for creating small edge
 * sequences. The edges are simply given as function parameters,
 * the end of the edges are marked by a -1 parameter. 
 * 
 * \param igraph The underlying graph object.
 * \param es The edge sequence object.
 * \param ... The ids of the edges are given by additional
 *   parameters closed by a -1 parameter.
 * \return Error code, \c IGRAPH_ENOMEM if there is
 *   not enough memory to create the edge sequence.
 * 
 * \sa \ref igraph_es_vectorview() and \ref igraph_es_vector().
 * 
 * Time complexity: usually O(n),
 * the number of edges in the edge sequence.
 */

int igraph_es_vector_small(const igraph_t *igraph, igraph_es_t *es, ...) {
  va_list ap;
  igraph_i_es_vectorview_pdata_t *data;
  igraph_vector_t *fakev;
  long int i, n=0;
  es->type=IGRAPH_ITERATOR_ES_VECTOR;
  es->stdata[0]=0;
  es->table=&igraph_i_es_vectorview_table;
  es->shorthand=0;
  
  es->pdata=Calloc(1, igraph_i_es_vectorview_pdata_t);
  if (es->pdata==0) {
    IGRAPH_ERROR("Cannot create vector iterator", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, es->pdata);
  data=(igraph_i_es_vectorview_pdata_t*)es->pdata;
  fakev=(igraph_vector_t*)&data->v;
  data->destroy=1;

  va_start(ap, es);
  while (1) {
    int num = va_arg(ap, int);
    if (num == -1) {
      break;
    }
    n++;
  }
  va_end(ap);

  IGRAPH_VECTOR_INIT_FINALLY(fakev, n);
  
  va_start(ap, es);
  for (i=0; i<n; i++) {
    VECTOR(*fakev)[i]=(real_t) va_arg(ap, int);
  }
  va_end(ap);  
  
  es->stdata[1]=igraph_vector_size(&data->v);
  
  IGRAPH_FINALLY_CLEAN(2);
  return 0;
}

/** 
 * \function IGRAPH_ES_VECTOR
 * 
 * Edge sequence shorthand for a \type igraph_vector_t object. 
 * \param graph The underlying graph object.
 * \param eids The \type igraph_vector_t object containing the edge
 * ids. 
 * \return Edge sequence shorthand.
 *
 * Time complexity: O(1).
 */

const igraph_es_t *IGRAPH_ES_VECTOR(const igraph_t *graph, 
				    const igraph_vector_t *eids) {
  igraph_es_t *es=Calloc(1, igraph_es_t);
  if (es==0) {
    igraph_error("Cannot create iterator shorthand", __FILE__, __LINE__,
		 IGRAPH_ENOMEM);
    return 0;			/* TODO: how to sign error ??? */
  }
  igraph_es_vectorview(graph, es, eids);
  es->shorthand=1;
  return es;
}

/** 
 * \function IGRAPH_ES
 * 
 * Edge sequence shorthand for edge ids given as parameters.
 * \param graph The underlying graph object.
 * \param ... The edge ids for the sequence, the end of them is
 * denoted by -1.
 * \return Edge sequence shorthand.
 * 
 * Time complexity: O(n), the number
 * of edge ids.
 */

const igraph_es_t *IGRAPH_ES(const igraph_t *graph, ...) {
  igraph_es_t *es=Calloc(1, igraph_es_t);
  va_list ap;
  igraph_i_es_vectorview_pdata_t *data;
  igraph_vector_t *fakev;
  long int i, n=0;
  int ret;
  if (es==0) {
    igraph_error("Cannot create iterator shorthand", __FILE__, __LINE__,
		 IGRAPH_ENOMEM);
    return 0;			/* TODO: how to sign error ??? */
  }

  es->type=IGRAPH_ITERATOR_ES_VECTOR;
  es->stdata[0]=0;
  es->table=&igraph_i_es_vectorview_table;
  es->shorthand=0;
  
  es->pdata=Calloc(1, igraph_i_es_vectorview_pdata_t);
  if (es->pdata==0) {
    igraph_error("Cannot create vector iterator", __FILE__, __LINE__,
		 IGRAPH_ENOMEM);
    return 0;
  }
  IGRAPH_FINALLY(igraph_free, es->pdata);
  data=(igraph_i_es_vectorview_pdata_t*)es->pdata;
  fakev=(igraph_vector_t*)&data->v;
  data->destroy=1;

  va_start(ap, graph);
  while (1) {
    int num = va_arg(ap, int);
    if (num == -1) {
      break;
    }
    n++;
  }
  va_end(ap);

  ret=igraph_vector_init(fakev, n);
  if (ret != 0) {
    igraph_error("Cannot create vector for iterator shorthand", __FILE__, 
		 __LINE__, ret);
    return 0;
  }
  IGRAPH_FINALLY(igraph_vector_destroy, fakev);
  
  va_start(ap, graph);
  for (i=0; i<n; i++) {
    VECTOR(*fakev)[i]=(real_t) va_arg(ap, int);
  }
  va_end(ap);  
  
  es->stdata[1]=igraph_vector_size(&data->v);
  
  IGRAPH_FINALLY_CLEAN(2);
  es->shorthand=1;
  return es;  
}

void igraph_es_next_vectorview(const igraph_t *graph, igraph_es_t *es) {
  es->stdata[0] ++;
}

bool_t igraph_es_end_vectorview(const igraph_t *graph, const igraph_es_t *es) {
  return es->stdata[0] >= es->stdata[1];
}

void igraph_es_reset_vectorview(const igraph_t *graph, igraph_es_t *es) {
  es->stdata[0]=0;
}

integer_t igraph_es_get_vectorview(const igraph_t *graph, 
				   const igraph_es_t *es) {
  igraph_i_es_vectorview_pdata_t *data=
    (igraph_i_es_vectorview_pdata_t*)es->pdata;
  return VECTOR(data->v)[ (long int) (es->stdata[0]) ]; 
}

integer_t igraph_es_from_vectorview(const igraph_t *graph, 
				    const igraph_es_t *es) {
  igraph_i_es_vectorview_pdata_t *data=
    (igraph_i_es_vectorview_pdata_t*)es->pdata;
  long int id=VECTOR(data->v)[ (long int) (es->stdata[0]) ];
  return VECTOR(graph->from) [id];
}

integer_t igraph_es_to_vectorview(const igraph_t *graph, 
				  const igraph_es_t *es) {
  igraph_i_es_vectorview_pdata_t *data=
    (igraph_i_es_vectorview_pdata_t*)es->pdata;
  long int id=VECTOR(data->v)[ (long int) (es->stdata[0]) ];
  return VECTOR(graph->to) [id];
}

int igraph_es_unfold_vectorview(const igraph_t *graph, const igraph_es_t *es,
				igraph_vector_t *v) {
  igraph_vector_t v2;
  igraph_i_es_vectorview_pdata_t *data=
    (igraph_i_es_vectorview_pdata_t*)es->pdata;
  IGRAPH_CHECK(igraph_vector_copy(&v2, &data->v));
  igraph_vector_destroy(v);
  *v=v2;
  return 0;
}

void igraph_es_destroy_vectorview(igraph_es_t *pes) {
  igraph_es_t *es=(igraph_es_t*)pes;
  igraph_i_es_vectorview_pdata_t *data=
    (igraph_i_es_vectorview_pdata_t*)es->pdata;
  if (data->destroy) {
    igraph_vector_destroy((igraph_vector_t*)&data->v);
  }
  Free(data);
  
  if (es->shorthand) {
    Free(pes);
  }
}

/**
 * \function igraph_es_vector_getvector
 * 
 * This is a specific edge sequence function, it can be called for
 * vector type iterators only (these are created by \ref
 * igraph_es_vector(), \ref igraph_es_vectorview(), \ref
 * igraph_es_vectorview_it() or \ref igraph_es_vector_small()).
 * It gives access to the edge ids in the sequence as a
 * \type igraph_vector_t type.
 * 
 * The result is undefined if you call it with a different iterator
 * type. 
 * 
 * \param graph The underlying graph object.
 * \param es The edge sequence. 
 * \return Pointer to a \type igraph_vector_t object. This object should
 *   be considered as constant, don't change its elements or size.
 * 
 * Time complexity: O(1).
 */

const igraph_vector_t *igraph_es_vector_getvector(const igraph_t *graph, 
					   const igraph_es_t *es) {
  igraph_i_es_vectorview_pdata_t *data=
    (igraph_i_es_vectorview_pdata_t*)es->pdata;
  return &data->v;  
}

/* -------------------------------------------------- */
/* edge iterator, all edges between two vertex sets   */
/* -------------------------------------------------- */

/**
 * \function igraph_es_fromto
 * 
 * Edge sequence containing all edges between two vertex sets.
 * \ref igraph_es_vector_getvector() can be called on this edge
 * sequence. 
 * \param graph The underlying graph object.
 * \param es The uninitialized edge sequence objects.
 * \param from The first vertex sequence.
 * \param to The second vertex sequence. 
 * \param directed Logical, if true and \p graph
 *   is directed only edges \em from \p from \em to
 *   \p to will be included. If false or
 *   \p graph is undirected all edges between
 *   \p from and \p to will be
 *   included. 
 *
 * Time complexity: O(dn),
 * n is the number of vertices in
 * \p from and d
 * is the average (out- or total- depending on
 * \p directed) degree of the vertices in
 * \p from.
 */

int igraph_es_fromto(const igraph_t *graph, igraph_es_t *es, 
		     const igraph_vs_t *from, const igraph_vs_t *to, 
		     bool_t directed) {

  igraph_vs_t myfrom, myto;  
  long int i, j, lfrom;
  igraph_es_t edgeit;
  igraph_neimode_t mode;
  const igraph_vector_t *fromvect;
  igraph_vector_t tovect;
  
  igraph_i_es_vectorview_pdata_t *data;

  if (directed && igraph_is_directed(graph)) {
    mode=IGRAPH_OUT;
  } else {
    mode=IGRAPH_ALL;
  }

  IGRAPH_CHECK(igraph_vs_vectorview_it(graph, from, &myfrom));
  IGRAPH_FINALLY(igraph_vs_destroy, &myfrom);
  IGRAPH_CHECK(igraph_vs_vectorview_it(graph, to, &myto));
  IGRAPH_FINALLY(igraph_vs_destroy, &myto);  

  fromvect=igraph_vs_vector_getvector(graph, &myfrom);
  IGRAPH_CHECK(igraph_vector_copy(&tovect, igraph_vs_vector_getvector(graph, &myto)));
  IGRAPH_FINALLY(igraph_vector_destroy, &tovect);

  es->type=IGRAPH_ITERATOR_ES_VECTOR;
  es->shorthand=0;
  es->stdata[0]=0;
  es->table=&igraph_i_es_vectorview_table;    
  es->pdata=Calloc(1, igraph_i_es_vectorview_pdata_t);
  if (es->pdata == 0) {
    IGRAPH_ERROR("Cannot create iterator", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, es->pdata);
  data=(igraph_i_es_vectorview_pdata_t*)es->pdata;
  IGRAPH_VECTOR_INIT_FINALLY((igraph_vector_t*) &data->v, 0);
  data->destroy=1;

  lfrom=igraph_vector_size(fromvect);
  if (lfrom != 0 && igraph_vector_size(&tovect) != 0) {
    IGRAPH_CHECK(igraph_es_adj(graph, &edgeit, 0, mode));
    IGRAPH_FINALLY(igraph_es_destroy, &edgeit);
    igraph_vector_sort(&tovect);
    for (i=0; i<lfrom; i++) {
      long int vfrom=VECTOR(*fromvect)[i];
      igraph_es_adj_set(graph, &edgeit, vfrom, mode);
      while (!igraph_es_end(graph, &edgeit)) {
	long int vto=igraph_es_adj_vertex(graph, &edgeit);
	if (igraph_vector_binsearch(&tovect, vto, 0)) {
	  igraph_vector_push_back((igraph_vector_t*)&data->v, igraph_es_get(graph, &edgeit));
	}
	igraph_es_next(graph, &edgeit);
      }
    }
    IGRAPH_FINALLY_CLEAN(1);
  }

  es->stdata[1]=igraph_vector_size(&data->v);

  /* Clean */
  igraph_vs_destroy(&myfrom);
  igraph_vs_destroy(&myto);
  igraph_vector_destroy(&tovect);
  IGRAPH_FINALLY_CLEAN(5);

  if (from->shorthand) { igraph_vs_destroy((igraph_vs_t*) from); }
  if (to->shorthand) { igraph_vs_destroy((igraph_vs_t*) to); }

  return 0;
}

