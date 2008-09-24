/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2003, 2004, 2005  Gabor Csardi <csardi@rmki.kfki.hu>
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
#include "config.h"

#include <string.h>   /* memset */

/**
 * \section about_adjlists
 * <para>Sometimes it is easier to work with a graph which is in
 * adjacency list format: a list of vectors; each vector contains the
 * neighbor vertices or adjacent edges of a given vertex. Typically,
 * this representation is good if we need to iterate over the neigbors
 * of all vertices many times. E.g. when finding the shortest paths
 * between every pairs of vertices or calculating closeness centrality
 * for all the vertices.</para>
 * 
 * <para>The <type>igraph_adjlist_t</type> stores the adjacency lists
 * of a graph. After creation it is independent of the original graph,
 * it can be modified freely with the usual vector operations, the
 * graph is not affected. E.g. the adjacency list can be used to
 * rewire the edges of a graph efficiently. If one used the
 * straightforward \ref igraph_delete_edges() and \ref
 * igraph_add_edges() combination for this that needs O(|V|+|E|) time
 * for every single deletion and insertion operation, it is thus very
 * slow if many edges are rewired. Extracting the graph into an
 * adjacency list, do all the rewiring operations on the vectors of
 * the adjacency list and then creating a new graph needs (depending
 * on how exactly the rewiring is done) typically O(|V|+|E|) time for
 * the whole rewiring process.</para>
 * 
 * <para>Lazy adjacency lists are a bit different. When creating a
 * lazy adjacency list, the neighbors of the vertices are not queried,
 * only some memory is allocated for the vectors. When \ref
 * igraph_lazy_adjlist_get() is called for vertex v the first time,
 * the neighbors of v are queried and stored in a vector of the
 * adjacency list, so they don't need to be queried again. Lazy
 * adjacency lists are handy if you have an at least linear operation
 * (because initialization is generally linear in terms of number of
 * vertices), but you don't know how many vertices you will visit
 * during the computation.
 * </para>
 * 
 */

/**
 * \function igraph_adjlist_init
 * Initialize an adjacency list of vertices
 * 
 * Create a list of vectors containing the neighbors of all vertices
 * in a graph. The adjacency list is independent of the graph after
 * creation, e.g. the graph can be destroyed and modified, the
 * adjacency list contains the state of the graph at the time of its
 * initialization. 
 * \param graph The input graph. 
 * \param al Pointer to an uninitialized <type>igraph_adjlist_t</type> object.
 * \param mode Constant specifying whether outgoing
 *   (<code>IGRAPH_OUT</code>), incoming (<code>IGRAPH_IN</code>),
 *   or both (<code>IGRAPH_ALL</code>) types of neighbors to include
 *   in the adjacency list. It is ignored for undirected networks.
 * \return Error code.
 * 
 * Time complexity: O(|V|+|E|), linear in the number of vertices and
 * edges.
 */

int igraph_adjlist_init(const igraph_t *graph, igraph_adjlist_t *al, 
			  igraph_neimode_t mode) {
  long int i;

  if (mode != IGRAPH_IN && mode != IGRAPH_OUT && mode != IGRAPH_ALL) {
    IGRAPH_ERROR("Cannot create adjlist view", IGRAPH_EINVMODE);
  }

  if (!igraph_is_directed(graph)) { mode=IGRAPH_ALL; }

  al->length=igraph_vcount(graph);
  al->adjs=igraph_Calloc(al->length, igraph_vector_t);
  if (al->adjs == 0) {
    IGRAPH_ERROR("Cannot create adjlist view", IGRAPH_ENOMEM);
  }

  IGRAPH_FINALLY(igraph_adjlist_destroy, al);
  for (i=0; i<al->length; i++) {
    IGRAPH_ALLOW_INTERRUPTION();
    IGRAPH_CHECK(igraph_vector_init(&al->adjs[i], 0));
    IGRAPH_CHECK(igraph_neighbors(graph, &al->adjs[i], i, mode));
  }

  IGRAPH_FINALLY_CLEAN(1);
  return 0;
}

/**
 * \function igraph_adjlist_init_complementer
 * Adjacency lists for the complementer graph
 * 
 * This function creates adjacency lists for the complementer 
 * of the input graph. In the complementer graph all edges are present
 * which are not present in the original graph. Multiple edges in the
 * input graph are ignored.
 * \param graph The input graph.
 * \param al Pointer to a not yet initialized adjacency list.
 * \param mode Constant specifying whether outgoing
 *   (<code>IGRAPH_OUT</code>), incoming (<code>IGRAPH_IN</code>),
 *   or both (<code>IGRAPH_ALL</code>) types of neighbors (in the
 *   complementer graph) to include in the adjacency list. It is
 *   ignored for undirected networks.
 * \param loops Whether to consider loop edges.
 * \return Error code.
 * 
 * Time complexity: O(|V|^2+|E|), quadratic in the number of vertices.
 */

int igraph_adjlist_init_complementer(const igraph_t *graph,
				       igraph_adjlist_t *al, 
				       igraph_neimode_t mode,
				       igraph_bool_t loops) {
  long int i, j, k, n;
  igraph_bool_t* seen;
  igraph_vector_t vec;

  if (mode != IGRAPH_IN && mode != IGRAPH_OUT && mode != IGRAPH_ALL) {
    IGRAPH_ERROR("Cannot create complementer adjlist view", IGRAPH_EINVMODE);
  }

  if (!igraph_is_directed(graph)) { mode=IGRAPH_ALL; }

  al->length=igraph_vcount(graph);
  al->adjs=igraph_Calloc(al->length, igraph_vector_t);
  if (al->adjs == 0) {
    IGRAPH_ERROR("Cannot create complementer adjlist view", IGRAPH_ENOMEM);
  }

  IGRAPH_FINALLY(igraph_adjlist_destroy, al);

  n=al->length;
  seen=igraph_Calloc(n, igraph_bool_t);
  if (seen==0) {
    IGRAPH_ERROR("Cannot create complementer adjlist view", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, seen);

  IGRAPH_VECTOR_INIT_FINALLY(&vec, 0);

  for (i=0; i<al->length; i++) {
    IGRAPH_ALLOW_INTERRUPTION();
    igraph_neighbors(graph, &vec, i, mode);
    memset(seen, 0, sizeof(igraph_bool_t)*al->length);
    n=al->length;
    if (!loops) { seen[i] = 1; n--; }
    for (j=0; j<igraph_vector_size(&vec); j++) {
      if (! seen [ (long int) VECTOR(vec)[j] ] ) {
	n--;
	seen[ (long int) VECTOR(vec)[j] ] = 1;
      }
    }
    IGRAPH_CHECK(igraph_vector_init(&al->adjs[i], n));
    for (j=0, k=0; k<n; j++) {
      if (!seen[j]) {
	VECTOR(al->adjs[i])[k++] = j;
      }
    }
  }

  igraph_Free(seen);
  igraph_vector_destroy(&vec);
  IGRAPH_FINALLY_CLEAN(3);
  return 0;
}

/**
 * \function igraph_adjlist_destroy
 * Deallocate memory
 * 
 * Free all memory allocated for an adjacency list. 
 * \param al The adjacency list to destroy.
 * 
 * Time complexity: depends on memory management.
 */

void igraph_adjlist_destroy(igraph_adjlist_t *al) {
  long int i;
  for (i=0; i<al->length; i++) {
    if (&al->adjs[i]) { igraph_vector_destroy(&al->adjs[i]); }
  }
  igraph_Free(al->adjs);
}

/**
 * \function igraph_adjlist_size
 * Number of vertices in an adjacency list.
 * 
 * \param al The adjacency list.
 * \return The number of elements.
 * 
 * Time complexity: O(1).
 */

igraph_integer_t igraph_adjlist_size(const igraph_adjlist_t *al) {
  return al->length;
}

/* igraph_vector_t *igraph_adjlist_get(igraph_adjlist_t *al, igraph_integer_t no) { */
/*   return &al->adjs[(long int)no]; */
/* } */

/**
 * \function igraph_adjlist_sort
 * Sort each vector in an adjacency list.
 * 
 * Sorts every vector of the adjacency list.
 * \param al The adjacency list.
 * 
 * Time complexity: O(n log n), n is the total number of elements in
 * the adjacency list.
 */

void igraph_adjlist_sort(igraph_adjlist_t *al) {
  long int i;
  for (i=0; i<al->length; i++)
    igraph_vector_sort(&al->adjs[i]);
}

/**
 * \function igraph_adjlist_simplify
 * Simplify
 * 
 * Simplify an adjacency list, ie. remove loop and multiple edges.
 * \param al The adjacency list.
 * \return Error code.
 * 
 * Time complexity: O(|V|+|E|), linear in the number of edges and
 * vertices.
 */

int igraph_adjlist_simplify(igraph_adjlist_t *al) {
  long int i;
  long int n=al->length;
  igraph_vector_t mark;
  IGRAPH_VECTOR_INIT_FINALLY(&mark, n);
  for (i=0; i<n; i++) {
    igraph_vector_t *v=&al->adjs[i];
    long int j, l=igraph_vector_size(v);
    VECTOR(mark)[i] = i+1;
    for (j=0; j<l; /* nothing */) {
      long int e=VECTOR(*v)[j];
      if (VECTOR(mark)[e] != i+1) {
	VECTOR(mark)[e]=i+1;
	j++;
      } else {
	VECTOR(*v)[j] = igraph_vector_tail(v);
	igraph_vector_pop_back(v);
	l--;
      }
    }
  }
  
  igraph_vector_destroy(&mark);
  IGRAPH_FINALLY_CLEAN(1);
  return 0;
}

/**
 * \function igraph_adjedgelist_init
 * Initialize an adjacency list of edges
 * 
 * Create a list of vectors containing the adjacent edges for all
 * vertices. The adjacency list is independent of the graph after
 * creation, subsequent changes of the graph object do not update the
 * adjacency list, and changes to the adjacency list do no update the
 * graph.
 * \param graph The input graph.
 * \param ael Pointer to an uninitialized adjcency list.
 * \param mode Constant specifying whether incoming edges
 *   (<code>IGRAPH_IN</code>), outgoing edges (<code>IGRAPH_OUT</code>) or
 *   both (<code>IGRAPH_ALL</code>) to include in the adjacency lists
 *   of directed graphs. It is ignored for undirected graphs.
 * \return Error code.
 * 
 * Time complexity: O(|V|+|E|), linear in the number of vertices and
 * edges.
 */

int igraph_adjedgelist_init(const igraph_t *graph, 
			      igraph_adjedgelist_t *ael, 
			      igraph_neimode_t mode) {
  long int i;

  if (mode != IGRAPH_IN && mode != IGRAPH_OUT && mode != IGRAPH_ALL) {
    IGRAPH_ERROR("Cannot create adjedgelist view", IGRAPH_EINVMODE);
  }

  if (!igraph_is_directed(graph)) { mode=IGRAPH_ALL; }

  ael->length=igraph_vcount(graph);
  ael->adjs=igraph_Calloc(ael->length, igraph_vector_t);
  if (ael->adjs == 0) {
    IGRAPH_ERROR("Cannot create adjedgelist view", IGRAPH_ENOMEM);
  }

  IGRAPH_FINALLY(igraph_adjlist_destroy, ael);  
  for (i=0; i<ael->length; i++) {
    IGRAPH_ALLOW_INTERRUPTION();
    IGRAPH_CHECK(igraph_vector_init(&ael->adjs[i], 0));
    IGRAPH_CHECK(igraph_adjacent(graph, &ael->adjs[i], i, mode));
  }
  
  IGRAPH_FINALLY_CLEAN(1);
  return 0;
}

/**
 * \function igraph_adjedgelist_destroy
 * Destroy
 * 
 * Free all memory allocated for an adjacency list.
 * \param eal The adjcency list to destroy.
 * 
 * Time complexity: depends on memory management.
 */

void igraph_adjedgelist_destroy(igraph_adjedgelist_t *ael) {
  long int i;
  for (i=0; i<ael->length; i++) {
    /* This works if some igraph_vector_t's are 0, because igraph_vector_destroy can
       handle this. */
    igraph_vector_destroy(&ael->adjs[i]);
  }
  igraph_Free(ael->adjs);
}

/**
 * \function igraph_lazy_adjlist_init
 * Constructor
 *
 * Create a lazy adjacency list for vertices. This function only
 * allocates some memory for storing the vectors of an adjacency list,
 * but the neighbor vertices are not queried, only at the \ref
 * igraph_lazy_adjlist_get() calls. 
 * \param graph The input graph.
 * \param al Pointer to an uninitialized adjacency list object.
 * \param mode Constant, it gives whether incoming edges
 *   (<code>IGRAPH_IN</code>), outgoing edges
 *   (<code>IGRPAH_OUT</code>) or both types of edges
 *   (<code>IGRAPH_ALL</code>) are considered. It is ignored for
 *   undirected graphs.
 * \param simplify Constant, it gives whether to simplify the vectors
 *   in the adjacency list (<code>IGRAPH_SIMPLIFY</code>) ot not
 *   (<code>IGRAPH_DONT_SIMPLIFY</code>).
 * \return Error code.
 * 
 * Time complexity: O(|V|), the number of vertices, possibly, but
 * depends on the underlying memory management too.
 */

int igraph_lazy_adjlist_init(const igraph_t *graph,
			       igraph_lazy_adjlist_t *al,
			       igraph_neimode_t mode,
			       igraph_lazy_adlist_simplify_t simplify) {
  if (mode != IGRAPH_IN && mode != IGRAPH_OUT && mode != IGRAPH_ALL) {
    IGRAPH_ERROR("Cannor create adjlist view", IGRAPH_EINVMODE);
  }

  if (!igraph_is_directed(graph)) { mode=IGRAPH_ALL; }  
  al->mode=mode;
  al->simplify=simplify;
  al->graph=graph;
  
  al->length=igraph_vcount(graph);
  al->adjs=igraph_Calloc(al->length, igraph_vector_t*);
  if (al->adjs == 0) {
    IGRAPH_ERROR("Cannot create lazy adjlist view", IGRAPH_ENOMEM);
  }

  return 0;
}

/**
 * \function igraph_lazy_adjlist_destroy
 * Deallocate memory
 * 
 * Free all allocated memory for a lazy adjacency list.
 * \param al The adjacency list to deallocate.
 * 
 * Time complexity: depends on the memory management.
 */

void igraph_lazy_adjlist_destroy(igraph_lazy_adjlist_t *al) {
  long int i, n=al->length;
  for (i=0; i<n; i++) {
    if (al->adjs[i] != 0) {
      igraph_vector_destroy(al->adjs[i]);
      igraph_Free(al->adjs[i]);
    }
  }
  igraph_Free(al->adjs);
}

igraph_vector_t *igraph_lazy_adjlist_get_real(igraph_lazy_adjlist_t *al,
						igraph_integer_t pno) {
  long int no=pno;
  int ret;
  if (al->adjs[no] == 0) {
    al->adjs[no] = igraph_Calloc(1, igraph_vector_t);
    if (al->adjs[no] == 0) {
      igraph_error("Lazy adjlist failed", __FILE__, __LINE__, 
		   IGRAPH_ENOMEM);
    }
    ret=igraph_vector_init(al->adjs[no], 0);
    if (ret != 0) {
      igraph_error("", __FILE__, __LINE__, ret);
    }
    ret=igraph_neighbors(al->graph, al->adjs[no], no, al->mode);
    if (ret != 0) {
      igraph_error("", __FILE__, __LINE__, ret);
    }

    if (al->simplify == IGRAPH_SIMPLIFY) {
      igraph_vector_t *v=al->adjs[no];
      long int i, p=0, n=igraph_vector_size(v);
      for (i=0; i<n; i++) {
	if (VECTOR(*v)[i] != no && 
	    (i==n-1 || VECTOR(*v)[i+1] != VECTOR(*v)[i])) {
	  VECTOR(*v)[p]=VECTOR(*v)[i];
	  p++;
	}
      }
      igraph_vector_resize(v, p);
    }
  }
  
  return al->adjs[no];
}

/**
 * \function igraph_lazy_adjedgelist_init
 * Constructor
 * 
 * Create a lazy adjacency list for edges. This function only
 * allocates some memory for storing the vectors of an adjacency list,
 * but the adjacent edges are not queried, only when \ref
 * igraph_lazy_adjedgelist_get() is called.
 * \param graph The input graph.
 * \param al Pointer to an uninitialized adjacency list.
 * \param mode Constant, it gives whether incoming edges
 *   (<code>IGRAPH_IN</code>), outgoing edges
 *   (<code>IGRPAH_OUT</code>) or both types of edges
 *   (<code>IGRAPH_ALL</code>) are considered. It is ignored for
 *   undirected graphs.
 * \return Error code.
 * 
 * Time complexity: O(|V|), the number of vertices, possibly. But it
 * also depends on the underlying memory management too.
 */

int igraph_lazy_adjedgelist_init(const igraph_t *graph,
				   igraph_lazy_adjedgelist_t *al,
				   igraph_neimode_t mode) {

  if (mode != IGRAPH_IN && mode != IGRAPH_OUT && mode != IGRAPH_ALL) {
    IGRAPH_ERROR("Cannot create adjlist view", IGRAPH_EINVMODE);
  }
  
  if (!igraph_is_directed(graph)) { mode=IGRAPH_ALL; }
  
  al->mode=mode;
  al->graph=graph;
  
  al->length=igraph_vcount(graph);
  al->adjs=igraph_Calloc(al->length, igraph_vector_t*);
  if (al->adjs == 0) {
    IGRAPH_ERROR("Cannot create lazy adjedgelist view", IGRAPH_ENOMEM);
  }

  return 0;
  
}

/**
 * \function igraph_lazy_adjedgelist_destroy
 * Deallocate memory
 * 
 * Free all allocated memory for a lazy edge adjacency list.
 * \param al The adjacency list to deallocate.
 * 
 * Time complexity: depends on memory management.
 */

void igraph_lazy_adjedgelist_destroy(igraph_lazy_adjedgelist_t *al) {
  long int i, n=al->length;
  for (i=0; i<n; i++) {
    if (al->adjs[i] != 0) {
      igraph_vector_destroy(al->adjs[i]);
      igraph_Free(al->adjs[i]);
    }
  }
  igraph_Free(al->adjs);
}

igraph_vector_t *igraph_lazy_adjedgelist_get_real(igraph_lazy_adjedgelist_t *al,
						    igraph_integer_t pno) {
  long int no=pno;
  int ret;
  if (al->adjs[no] == 0) {
    al->adjs[no] = igraph_Calloc(1, igraph_vector_t);
    if (al->adjs[no] == 0) {
      igraph_error("Lazy adjedgelist failed", __FILE__, __LINE__, 
		   IGRAPH_ENOMEM);
    }
    ret=igraph_vector_init(al->adjs[no], 0);
    if (ret != 0) {
      igraph_error("", __FILE__, __LINE__, ret);
    }
    ret=igraph_adjacent(al->graph, al->adjs[no], no, al->mode);
    if (ret != 0) {
      igraph_error("", __FILE__, __LINE__, ret);
    }
  }
  return al->adjs[no];
}
