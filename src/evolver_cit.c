/* -*- mode: C -*-  */
/* 
   IGraph R package.
   Copyright (C) 2007-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge, MA 02139 USA
   
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

#include "igraph_revolver.h"
#include "igraph_memory.h"
#include "igraph_random.h"
#include "igraph_constructors.h"
#include "igraph_psumtree.h"
#include "config.h"

#include <math.h>

/* This file contains the method for creating citation networks */

int igraph_i_create_outseq(igraph_vector_t *real_outseq,
			   igraph_integer_t nodes,
			   const igraph_vector_t *outseq,
			   const igraph_vector_t *outdist,
			   igraph_integer_t m,
			   igraph_integer_t *edges) {
  
  long int no_of_edges=0;

  if (outseq && nodes != igraph_vector_size(outseq)) {
    IGRAPH_ERROR("Invalid out-degree sequence length", IGRAPH_EINVAL);
  }
  if (!outseq && outdist && igraph_vector_size(outdist)==0) {
    IGRAPH_ERROR("Invalid out-degree distribution length", IGRAPH_EINVAL);
  }
  if (!outseq && !outdist && m<0) {
    IGRAPH_ERROR("Invalid constant out-degree", IGRAPH_EINVAL);
  }

  if (outseq) {
    igraph_vector_clear(real_outseq);
    igraph_vector_append(real_outseq, outseq);    
    no_of_edges=igraph_vector_sum(real_outseq)-VECTOR(*real_outseq)[0];
  } else if (outdist) {    
    igraph_vector_t cumsum;
    long int i, n=igraph_vector_size(outdist);   
    IGRAPH_VECTOR_INIT_FINALLY(&cumsum, n+1);
    IGRAPH_CHECK(igraph_vector_resize(real_outseq, nodes));
    VECTOR(cumsum)[0]=0;
    for (i=0; i<n; i++) {
      VECTOR(cumsum)[i+1] = VECTOR(cumsum)[i] + VECTOR(*outdist)[i];
    }
    RNG_BEGIN();
    VECTOR(*real_outseq)[0]=0;
    for (i=1; i<nodes; i++) {
      long int deg;
      igraph_vector_binsearch(&cumsum, RNG_UNIF(0, VECTOR(cumsum)[n]), &deg);
      VECTOR(*real_outseq)[0]=deg;
      no_of_edges += deg;
    }
    RNG_END();
    igraph_vector_destroy(&cumsum);
    IGRAPH_FINALLY_CLEAN(1);
  } else {
    long int i;
    for (i=0; i<nodes; i++) {
      VECTOR(*real_outseq)[i]=m;
    }
    no_of_edges=(nodes-1)*m;
  }
  
  if (edges) {
    *edges=no_of_edges;
  }
  
  return 0;
}

int igraph_evolver_d(igraph_t *graph,
		     igraph_integer_t nodes,
		     igraph_vector_t *kernel,
		     const igraph_vector_t *outseq,
		     const igraph_vector_t *outdist,
		     igraph_integer_t m,
		     igraph_bool_t directed) {

  igraph_vector_t real_outseq;
  long int no_of_nodes=nodes;
  igraph_integer_t no_of_edges;
  igraph_vector_t edges;
  long int kernel_size=igraph_vector_size(kernel);
  long int edgeptr=0;
  igraph_psumtree_t sumtree;
  igraph_vector_t degree;
  long int i, j;
  
  if (no_of_nodes<0) {
    IGRAPH_ERROR("Invalid number of vertices", IGRAPH_EINVAL);
  }
  if (kernel_size==0) {
    IGRAPH_ERROR("Zero length kernel", IGRAPH_EINVAL);
  }
  if (VECTOR(*kernel)[0]==0) {
    IGRAPH_ERROR("Zero attractivity for zero degree vertices no allowed", 
		 IGRAPH_EINVAL);
  }
  
  IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);

  IGRAPH_VECTOR_INIT_FINALLY(&real_outseq, no_of_nodes);
  IGRAPH_CHECK(igraph_i_create_outseq(&real_outseq, nodes, outseq, 
				      outdist, m, &no_of_edges));
  
  IGRAPH_CHECK(igraph_vector_resize(&edges, no_of_edges*2));
  IGRAPH_CHECK(igraph_psumtree_init(&sumtree, no_of_nodes));
  IGRAPH_FINALLY(igraph_psumtree_destroy, &sumtree);
  IGRAPH_VECTOR_INIT_FINALLY(&degree, no_of_nodes);
  
  RNG_BEGIN();
  
  /* first node */
  igraph_psumtree_update(&sumtree, 0, VECTOR(*kernel)[0]);
  
  for (i=1; i<no_of_nodes; i++) {
    igraph_real_t sum=igraph_psumtree_sum(&sumtree);
    long int no_of_neighbors=VECTOR(real_outseq)[i];
    long int edgeptr_save;
    long int to;
    
    /* Add edges */
    edgeptr_save=edgeptr;
    for (j=0; j<no_of_neighbors; j++) {
      igraph_psumtree_search(&sumtree, &to, RNG_UNIF(0, sum));
      VECTOR(degree)[to] += 1;
      VECTOR(edges)[edgeptr++] = i;
      VECTOR(edges)[edgeptr++] = to;
    }     
 
    /* Update probabilities */
    for (j=0; j<no_of_neighbors; j++) {
      long int to=VECTOR(edges)[edgeptr_save+j*2+1];
      long int deg=VECTOR(degree)[to];
      long int a= deg < kernel_size ? VECTOR(*kernel)[deg] : 
	VECTOR(*kernel)[kernel_size-1];
      igraph_psumtree_update(&sumtree, to, a);
    }
			     
    /* New vertex */
    igraph_psumtree_update(&sumtree, i, VECTOR(*kernel)[0]);
  }
  
  RNG_END();
  
  igraph_vector_destroy(&degree);
  igraph_psumtree_destroy(&sumtree);
  igraph_vector_destroy(&real_outseq);
  IGRAPH_FINALLY_CLEAN(3);
  
  IGRAPH_CHECK(igraph_create(graph, &edges, nodes, directed));
  igraph_vector_destroy(&edges);
  IGRAPH_FINALLY_CLEAN(1);
  
  return 0;
}
