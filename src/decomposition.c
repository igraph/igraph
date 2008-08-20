/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2008  Gabor Csardi <csardi@rmki.kfki.hu>
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
#include "error.h"

/**
 * \function igraph_maximum_cardinality_search
 * Maximum cardinality search
 * 
 * This function implements the maximum cardinality search algorithm
 * discussed in 
 * 
 * Robert E Tarjan and Mihalis Yannakakis, Simple linear-time
 * algorithms to test chordality of graphs, test acyclicity of
 * hypergraphs, and selectively reduce acyclic hypergraphs.
 * SIAM Journal of Computation 13, 566--579, 1984.
 * 
 * \param graph The input graph. Can be directed, but the direction
 *   of the edges is ignored.
 * \param res Pointer to an initialized vector, the result is stored here. 
 *   It will be resized, as needed. Upon return it contains
 *   the rank of the each vertex.
 * \return Error code.
 * 
 * Time complexity: O(|V|+|E|), linear in terms of the number of
 * vertices and edges.  
 */

int igraph_maximum_cardinality_search(const igraph_t *graph,
				      igraph_vector_t *res) {
  
  long int no_of_nodes=igraph_vcount(graph);
  igraph_vector_long_t size;
  igraph_vector_long_t head, next, prev; /* doubly linked list with head */
  long int i;
  igraph_adjlist_t adjlist;
  
  /***************/
  /* local j, v; */
  /***************/
  
  long int j, v;

  IGRAPH_CHECK(igraph_vector_long_init(&size, no_of_nodes));
  IGRAPH_FINALLY(igraph_vector_long_destroy, &size);
  IGRAPH_CHECK(igraph_vector_long_init(&head, no_of_nodes));
  IGRAPH_FINALLY(igraph_vector_long_destroy, &head);
  IGRAPH_CHECK(igraph_vector_long_init(&next, no_of_nodes));
  IGRAPH_FINALLY(igraph_vector_long_destroy, &next);
  IGRAPH_CHECK(igraph_vector_long_init(&prev, no_of_nodes));
  IGRAPH_FINALLY(igraph_vector_long_destroy, &prev);

  IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist, IGRAPH_ALL));
  IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist);
  
  IGRAPH_CHECK(igraph_vector_resize(res, no_of_nodes));
  
  /***********************************************/
  /* for i in [0,n-1] -> set(i) := emptyset rof; */
  /***********************************************/

  /* nothing to do, 'head' contains all zeros */

  /*********************************************************/
  /* for v in vertices -> size(v):=0; add v to set(0) rof; */
  /*********************************************************/

  VECTOR(head)[0]=1;
  for (v=0; v<no_of_nodes; v++) {
    VECTOR(next)[v]=v+2;
    VECTOR(prev)[v]=v;
  }
  VECTOR(next)[no_of_nodes-1] = 0;
  /* size is allready all zero */

  /***************/
  /* i:=n; j:=0; */
  /***************/

  i=no_of_nodes; j=0;
  
  /**************/
  /* do i>=1 -> */
  /**************/

  while (i>=1) {
    long int x, k, len;
    igraph_vector_t *neis;

    /********************************/
    /* v :=  delete any from set(j) */
    /********************************/

    v=VECTOR(head)[j]-1;
    x=VECTOR(next)[v];
    VECTOR(head)[j]=x;
    if (x != 0) {
      VECTOR(prev)[x-1]=0;
    }
    
    /*************************************************/
    /* alpha(v) := i; alpham1(i) := v; size(v) := -1 */
    /*************************************************/

    VECTOR(*res)[v]=i-1; 
/*     VECTOR(alpham1)[i-1]=v;	/\* This is not needed *\/ */
    VECTOR(size)[v]=-1;
    
    /********************************************/
    /* for {v,w} in E such that size(w) >= 0 -> */
    /********************************************/
    
    neis=igraph_adjlist_get(&adjlist, v);
    len=igraph_vector_size(neis);
    for (k=0; k<len; k++) {
      long int w=VECTOR(*neis)[k];
      long int ws=VECTOR(size)[w];
      if (ws >= 0) {

	/******************************/
	/* delete w from set(size(w)) */
	/******************************/
	
	long int nw=VECTOR(next)[w];
	long int pw=VECTOR(prev)[w];
	if (nw != 0) {
	  VECTOR(prev)[nw-1] = pw;
	}
	if (pw != 0) {
	  VECTOR(next)[pw-1] = nw;
	} else {
	  VECTOR(head)[ws]=nw;
	}

	/******************************/
	/* size(w) := size(w)+1       */
	/******************************/

	VECTOR(size)[w] += 1;

	/******************************/
	/* add w to set(size(w))      */
	/******************************/

	ws=VECTOR(size)[w];
	nw=VECTOR(head)[ws];
	VECTOR(next)[w]=nw;
	VECTOR(prev)[w]=0;
	if (nw != 0) {
	  VECTOR(prev)[nw-1]=w+1;
	}
	VECTOR(head)[ws]=w+1;

      }
    }
    
    /***********************/
    /* i := i-1; j := j+1; */
    /***********************/
    
    i -= 1;
    j += 1;
    
    /*********************************************/
    /* do j>=0 and set(j)=emptyset -> j:=j-1; od */ 
    /*********************************************/
    
    while (j>=0 && VECTOR(head)[j]==0) j--;
    
  }
  
  igraph_adjlist_destroy(&adjlist);
  igraph_vector_long_destroy(&prev);
  igraph_vector_long_destroy(&next);
  igraph_vector_long_destroy(&head);
  igraph_vector_long_destroy(&size);
  IGRAPH_FINALLY_CLEAN(5);
  
  return 0;
}
