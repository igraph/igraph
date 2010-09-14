/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2010  Gabor Csardi <csardi.gabor@gmail.com>
   Rue de l'Industrie 5, Lausanne 1005, Switzerland
   
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

#include "igraph_cohesive_blocks.h"
#include "igraph_interface.h"
#include "igraph_memory.h"
#include "igraph_flow.h"
#include "igraph_separators.h"
#include "igraph_structural.h"
#include "igraph_components.h"
#include "igraph_dqueue.h"
#include "igraph_constructors.h"

void igraph_i_cohesive_blocks_free(igraph_vector_ptr_t *ptr) {
  long int i, n=igraph_vector_ptr_size(ptr);

  for (i=0; i<n; i++) {
    igraph_t *g=VECTOR(*ptr)[i];
    if (g) {
      igraph_destroy(g);
      igraph_free(g);
    }
  }
}

void igraph_i_cohesive_blocks_free2(igraph_vector_ptr_t *ptr) {
  long int i, n=igraph_vector_ptr_size(ptr);

  for (i=0; i<n; i++) {
    igraph_vector_long_t *v=VECTOR(*ptr)[i];
    if (v) {
      igraph_vector_long_destroy(v);
      igraph_free(v);
    }
  }
}

void igraph_i_cohesive_blocks_free3(igraph_vector_ptr_t *ptr) {
  long int i, n=igraph_vector_ptr_size(ptr);

  for (i=0; i<n; i++) {
    igraph_vector_t *v=VECTOR(*ptr)[i];
    if (v) {
      igraph_vector_destroy(v);
      igraph_free(v);
    }
  }
}

/* This is kind of a BFS to find the components of the graph, after
 * deleting the vertices marked in 'excluded'.
 * These vertices are not put in the BFS queue, but they are added to
 * all neighboring components.
 */

int igraph_i_cb_components(igraph_t *graph, 
			   const igraph_vector_bool_t *excluded,
			   igraph_vector_long_t *components,
			   long int *no,
			   /* working area follows */
			   igraph_vector_long_t *compid, 
			   igraph_dqueue_t *Q,
			   igraph_vector_t *neis) {

  long int no_of_nodes=igraph_vcount(graph);
  long int i;
  long int cno=0;

  igraph_vector_long_clear(components);
  igraph_dqueue_clear(Q);
  IGRAPH_CHECK(igraph_vector_long_resize(compid, no_of_nodes));
  igraph_vector_long_null(compid);
  
  for (i=0; i<no_of_nodes; i++) {    

    if (VECTOR(*compid)[i])   { continue; }
    if (VECTOR(*excluded)[i]) { continue; }
    
    IGRAPH_CHECK(igraph_dqueue_push(Q, i));
    IGRAPH_CHECK(igraph_vector_long_push_back(components, i));
    VECTOR(*compid)[i] = ++cno;
    
    while (!igraph_dqueue_empty(Q)) {
      long int node=igraph_dqueue_pop(Q);
      long int j, n;
      IGRAPH_CHECK(igraph_neighbors(graph, neis, node, IGRAPH_ALL));
      n=igraph_vector_size(neis);
      for (j=0; j<n; j++) {
	long int v=VECTOR(*neis)[j];
	if (VECTOR(*excluded)[v]) {
	  if (VECTOR(*compid)[v] != cno) {
	    VECTOR(*compid)[v] = cno;
	    IGRAPH_CHECK(igraph_vector_long_push_back(components, v));
	  }
	} else {
	  if (!VECTOR(*compid)[v]) {
	    VECTOR(*compid)[v] = cno; /* could be anything positive */
	    IGRAPH_CHECK(igraph_vector_long_push_back(components, v));
	    IGRAPH_CHECK(igraph_dqueue_push(Q, v));
	  }
	}
      }
    } /* while !igraph_dqueue_empty */

    IGRAPH_CHECK(igraph_vector_long_push_back(components, -1));
    
  } /* for i<no_of_nodes */
  
  *no=cno;

  return 0;
}

int igraph_cohesive_blocks(const igraph_t *graph,
			   igraph_vector_ptr_t *blocks,
			   igraph_vector_t *cohesion,
			   igraph_vector_t *parent,
			   igraph_t *block_tree) {

  /* Some implementation comments. Everything is relatively
     straightforward, except, that we need to follow the vertex ids
     of the various subgraphs, without having to store two-way
     mappings at each level. The subgraphs can overlap, this
     complicates things a bit.

     The 'Q' vector is used as a double ended queue and it contains
     the subgraphs to work on in the future. Some other vectors are
     associated with it. 'Qparent' gives the parent graph of a graph
     in Q. Qmapping gives the mapping of the vertices from the graph
     to the parent graph. Qcohesion is the vertex connectivity of the 
     graph. 

     Qptr is an integer and points to the next graph to work on.
  */
  
  igraph_vector_ptr_t Q;
  igraph_vector_ptr_t Qmapping;
  igraph_vector_long_t Qparent;
  igraph_vector_long_t Qcohesion;
  long int Qptr=0;
  igraph_integer_t conn;

  igraph_t *graph_copy;
  
  igraph_vector_ptr_t separators;
  igraph_vector_t compvertices;
  igraph_vector_long_t components;
  igraph_vector_bool_t marked;

  igraph_vector_long_t compid;
  igraph_dqueue_t bfsQ;
  igraph_vector_t neis;

  if (igraph_is_directed(graph)) {
    IGRAPH_ERROR("Cohesive blocking only works on undirected graphs",
		 IGRAPH_EINVAL);
  }

  if (blocks)   { igraph_vector_ptr_clear(blocks); }
  if (cohesion) { igraph_vector_clear(cohesion);   }
  if (parent)   { igraph_vector_clear(parent);     }

  IGRAPH_CHECK(igraph_vector_ptr_init(&Q, 1));
  IGRAPH_FINALLY(igraph_vector_ptr_destroy, &Q);
  IGRAPH_FINALLY(igraph_i_cohesive_blocks_free, &Q);

  IGRAPH_CHECK(igraph_vector_ptr_init(&Qmapping, 1));
  IGRAPH_FINALLY(igraph_vector_ptr_destroy, &Qmapping);
  IGRAPH_FINALLY(igraph_i_cohesive_blocks_free2, &Qmapping);

  IGRAPH_CHECK(igraph_vector_long_init(&Qparent, 1));
  IGRAPH_FINALLY(igraph_vector_long_destroy, &Qparent);

  IGRAPH_CHECK(igraph_vector_long_init(&Qcohesion, 1));
  IGRAPH_FINALLY(igraph_vector_long_destroy, &Qcohesion);
  
  IGRAPH_CHECK(igraph_vector_ptr_init(&separators, 0));
  IGRAPH_FINALLY(igraph_vector_ptr_destroy, &separators);

  IGRAPH_VECTOR_INIT_FINALLY(&compvertices, 0);
  IGRAPH_CHECK(igraph_vector_bool_init(&marked, 0));
  IGRAPH_FINALLY(igraph_vector_bool_destroy, &marked);
  IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
  IGRAPH_CHECK(igraph_dqueue_init(&bfsQ, 100));
  IGRAPH_FINALLY(igraph_dqueue_destroy, &bfsQ);
  IGRAPH_CHECK(igraph_vector_long_init(&compid, 0));
  IGRAPH_FINALLY(igraph_vector_long_destroy, &compid);
  IGRAPH_CHECK(igraph_vector_long_init(&components, 0));
  IGRAPH_FINALLY(igraph_vector_long_destroy, &components);
  
  /* Put the input graph in the queue */
  graph_copy=igraph_Calloc(1, igraph_t);
  if (!graph_copy) { 
    IGRAPH_ERROR("Cannot do cohesive blocking", IGRAPH_ENOMEM);
  }
  IGRAPH_CHECK(igraph_copy(graph_copy, graph));
  VECTOR(Q)[0] = graph_copy;
  VECTOR(Qmapping)[0] = 0;	/* Identity mapping */
  VECTOR(Qparent)[0] = -1;	/* Has no parent */
  IGRAPH_CHECK(igraph_vertex_connectivity(graph, &conn, /*checks=*/ 1));
  VECTOR(Qcohesion)[0] = conn;
  
  /* Then work until the queue is empty */
  while (Qptr < igraph_vector_ptr_size(&Q)) {
    igraph_t *mygraph=VECTOR(Q)[Qptr];
    long int mynodes=igraph_vcount(mygraph);
    long int i, nsep;
    long int no;
    long int cptr=0;
    long int nsepv=0;

    /* Get the separators */
    IGRAPH_CHECK(igraph_minimum_size_separators(mygraph, &separators));
    IGRAPH_FINALLY(igraph_i_cohesive_blocks_free3, &separators);

    /* Remove them from the graph, also mark them */    
    IGRAPH_CHECK(igraph_vector_bool_resize(&marked, mynodes));
    igraph_vector_bool_null(&marked);
    for (i=0, nsep=igraph_vector_ptr_size(&separators); i<nsep; i++) {
      igraph_vector_t *v=VECTOR(separators)[i];
      long int j, n=igraph_vector_size(v);
      for (j=0; j<n; j++) {
	long int vv=VECTOR(*v)[j];
	if (!VECTOR(marked)[vv]) {
	  nsepv++;
	  VECTOR(marked)[vv] = 1;
	}
      }
    }
    
    /* Find the connected components, omitting the separator vertices, 
       but including the neighboring separator vertices
     */
    IGRAPH_CHECK(igraph_i_cb_components(mygraph, &marked, 
					&components, &no,
					&compid, &bfsQ, &neis));

    /* Add the separator vertices themselves, as another component, 
       but only if there is at least one vertex not included in any
       separator. */
    if (nsepv != mynodes) {
      for (i=0; i<mynodes; i++) {
	if (VECTOR(marked)[i]) {
	  IGRAPH_CHECK(igraph_vector_long_push_back(&components, i));
	}
      }
      IGRAPH_CHECK(igraph_vector_long_push_back(&components, -1));
      no++;
    }
    
    for (i=0; i<no; i++) {
      igraph_vector_t *newmapping;
      igraph_t *newgraph;
      igraph_integer_t maxdeg;

      igraph_vector_clear(&compvertices);
      while (1) {
	long int v=VECTOR(components)[cptr++];
	if (v < 0) { break; }
	IGRAPH_CHECK(igraph_vector_push_back(&compvertices, v));
      }
      
      newmapping=igraph_Calloc(1, igraph_vector_t);
      if (!newmapping) {
	IGRAPH_ERROR("Cannot do cohesive blocking", IGRAPH_ENOMEM);
      }
      IGRAPH_FINALLY(igraph_free, newmapping);
      IGRAPH_VECTOR_INIT_FINALLY(newmapping, 0);
      newgraph=igraph_Calloc(1, igraph_t);
      if (!newgraph) { 
	IGRAPH_ERROR("Cannot do cohesive blocking", IGRAPH_ENOMEM);
      }
      IGRAPH_FINALLY(igraph_free, newgraph);
      IGRAPH_CHECK(igraph_induced_subgraph_map(mygraph, newgraph, 
					igraph_vss_vector(&compvertices),
					IGRAPH_SUBGRAPH_AUTO,
					/*map=*/ 0,
					/*invmap=*/ newmapping));
      IGRAPH_FINALLY(igraph_destroy, newgraph);

      IGRAPH_CHECK(igraph_maxdegree(newgraph, &maxdeg, igraph_vss_all(),
				    IGRAPH_ALL, IGRAPH_LOOPS));
      if (maxdeg > VECTOR(Qcohesion)[Qptr]) {
	igraph_integer_t newconn;
	IGRAPH_CHECK(igraph_vector_ptr_push_back(&Q, newgraph));
	IGRAPH_FINALLY_CLEAN(2);
	IGRAPH_CHECK(igraph_vector_ptr_push_back(&Qmapping, newmapping));
	IGRAPH_FINALLY_CLEAN(2);
	IGRAPH_CHECK(igraph_vertex_connectivity(newgraph, &newconn, 
						/*checks=*/ 1));
	IGRAPH_CHECK(igraph_vector_long_push_back(&Qcohesion, newconn));
	IGRAPH_CHECK(igraph_vector_long_push_back(&Qparent, Qptr));
      } else {
	igraph_destroy(newgraph);
	igraph_free(newgraph);
	igraph_vector_destroy(newmapping);
	igraph_free(newmapping);
	IGRAPH_FINALLY_CLEAN(4);
      }
    }
    
    igraph_destroy(mygraph);
    igraph_free(mygraph);
    VECTOR(Q)[Qptr] = 0;
    igraph_i_cohesive_blocks_free3(&separators);
    IGRAPH_FINALLY_CLEAN(1);

    Qptr++;
  }

  igraph_vector_long_destroy(&components);
  igraph_vector_long_destroy(&compid);
  igraph_dqueue_destroy(&bfsQ);
  igraph_vector_destroy(&neis);
  igraph_vector_bool_destroy(&marked);
  igraph_vector_destroy(&compvertices);
  igraph_vector_ptr_destroy(&separators);
  IGRAPH_FINALLY_CLEAN(7);

  if (blocks || cohesion || parent || block_tree) {
    long int noblocks=Qptr, badblocks=0;
    igraph_vector_bool_t removed;
    long int i, resptr=0;
    igraph_vector_long_t rewritemap;
    
    IGRAPH_CHECK(igraph_vector_bool_init(&removed, noblocks));
    IGRAPH_FINALLY(igraph_vector_bool_destroy, &removed);
    IGRAPH_CHECK(igraph_vector_long_init(&rewritemap, noblocks));
    IGRAPH_FINALLY(igraph_vector_long_destroy, &rewritemap);

    for (i=1; i<noblocks; i++) {
      long int p=VECTOR(Qparent)[i];
      if (VECTOR(Qcohesion)[p] >= VECTOR(Qcohesion)[i]) {
	VECTOR(removed)[i]=1;
	badblocks++;
      }
    }
    noblocks -= badblocks;

    if (blocks) { IGRAPH_CHECK(igraph_vector_ptr_resize(blocks, noblocks)); }
    if (cohesion) { IGRAPH_CHECK(igraph_vector_resize(cohesion, noblocks)); }
    if (parent) { IGRAPH_CHECK(igraph_vector_resize(parent, noblocks)); }

    /* Rewrite the mappings */
    for (i=1; i<Qptr; i++) {
      long int p=VECTOR(Qparent)[i];
      igraph_vector_t *mapping=VECTOR(Qmapping)[i];
      igraph_vector_t *pmapping=VECTOR(Qmapping)[p];
      long int j, n=igraph_vector_size(mapping);

      if (!pmapping) { continue; }
      for (j=0; j<n; j++) {
	long int v=VECTOR(*mapping)[j];
	VECTOR(*mapping)[j] = VECTOR(*pmapping)[v];
      }
    }

    for (i=0; i<Qptr; i++) {
      if (VECTOR(removed)[i]) { continue; }
      VECTOR(rewritemap)[i] = resptr;
      if (cohesion) { VECTOR(*cohesion)[resptr]=VECTOR(Qcohesion)[i]; }
      if (parent || block_tree) {
	long int p=VECTOR(Qparent)[i];
	while (p>=0 && VECTOR(removed)[p]) { p=VECTOR(Qparent)[p]; }
	if (p>=0) { p=VECTOR(rewritemap)[p]; }
	VECTOR(Qparent)[i]=p;
	if (parent) { VECTOR(*parent)[resptr]=p; }
      }
      if (blocks) {
	VECTOR(*blocks)[resptr]=VECTOR(Qmapping)[i];
	VECTOR(Qmapping)[i]=0;
      }
      resptr++;
    }

    /* Plus the original graph */
    if (blocks) {
      igraph_vector_t *orig=igraph_Calloc(1, igraph_vector_t);
      if (!orig) { 
	IGRAPH_ERROR("Cannot do cohesive blocking", IGRAPH_ENOMEM); 
      }
      IGRAPH_FINALLY(igraph_free, orig);
      IGRAPH_CHECK(igraph_vector_init_seq(orig, 0, igraph_vcount(graph)-1));
      VECTOR(*blocks)[0]=orig;
      IGRAPH_FINALLY_CLEAN(1);
    }

    if (block_tree) {
      igraph_vector_t edges;
      long int eptr=0;
      IGRAPH_VECTOR_INIT_FINALLY(&edges, noblocks*2-2);
      for (i=1; i<Qptr; i++) {
	if (VECTOR(removed)[i]) { continue; }
	VECTOR(edges)[eptr++] = VECTOR(Qparent)[i];
	VECTOR(edges)[eptr++] = VECTOR(rewritemap)[i];
      }
      
      IGRAPH_CHECK(igraph_create(block_tree, &edges, noblocks, 
				 IGRAPH_DIRECTED));
      igraph_vector_destroy(&edges);
      IGRAPH_FINALLY_CLEAN(1);
    }

    igraph_vector_long_destroy(&rewritemap);
    igraph_vector_bool_destroy(&removed);
    IGRAPH_FINALLY_CLEAN(2);

  }

  igraph_vector_long_destroy(&Qcohesion);
  igraph_vector_long_destroy(&Qparent);
  igraph_i_cohesive_blocks_free2(&Qmapping);
  IGRAPH_FINALLY_CLEAN(3);
  
  igraph_vector_ptr_destroy(&Qmapping);
  igraph_vector_ptr_destroy(&Q);
  IGRAPH_FINALLY_CLEAN(3);  	/* + the elements of Q, they were
				   already destroyed */
  return 0;
}
