/* -*- mode: C -*-  */
/* 
   IGraph R library.
   Copyright (C) 2003, 2004  Gabor Csardi <csardi@rmki.kfki.hu>
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
#include "random.h"
#include "memory.h"

int igraph_barabasi_game(igraph_t *graph, integer_t n, integer_t m, 
			 vector_t *outseq, bool_t outpref, bool_t directed) {

/* SEXP REST_ba_game(SEXP pn, SEXP pm, SEXP outseq, SEXP poutpref) { */
  
  long int no_of_nodes=n;
  long int no_of_neighbors;
  long int *bag;
  long int bagp;
  long int no_of_edges;
  vector_t edges;
  
  long int resp=0;

  long int i,j;

  if (vector_size(outseq) == 0) {
    no_of_neighbors=m;
    bag=Calloc(no_of_nodes * no_of_neighbors + no_of_nodes +
	       outpref * no_of_nodes * no_of_neighbors,
	       long int);
    no_of_edges=(no_of_nodes-1)*no_of_neighbors;
  } else {
    no_of_edges=0;
    for (i=1; i<vector_size(outseq); i++) {
      no_of_edges+=VECTOR(*outseq)[i];
    }
    bag=Calloc(no_of_nodes + no_of_edges + outpref * no_of_edges,
	       long int);
  }
  
  vector_init(&edges, no_of_edges*2);
  
  /* The first node */

  bagp=0;
  bag[bagp++]=0;
  
  RNG_BEGIN();

  /* and the others */
  
  for (i=1; i<no_of_nodes; i++) {
    /* draw edges */
    if (vector_size(outseq)!=0) { no_of_neighbors=VECTOR(*outseq)[i]; }
    for (j=0; j<no_of_neighbors; j++) {
      long int to=bag[RNG_INTEGER(0, bagp-1)];
      VECTOR(edges)[resp++] = i;
      VECTOR(edges)[resp++] = to;
    }
    /* update bag */
    bag[bagp++] = i;
    for (j=0; j<no_of_neighbors; j++) {
      bag[bagp++] = VECTOR(edges)[resp-2*j-1];
      if (outpref) {
	bag[bagp++] = i;
      }
    }
  }

  RNG_END();

  Free(bag);  
  igraph_create(graph, &edges, 0, directed);
  vector_destroy(&edges);

  return 0;
}

int igraph_erdos_renyi_game(igraph_t *graph, integer_t n, real_t p,
			    bool_t directed, bool_t loops) {
  /* TODO */
  return 0;
}

int igraph_degree_sequence_game(igraph_t *graph, vector_t *out_deg,
				vector_t *in_deg, integer_t method) {
  /* TODO */
  return 0;
}

int igraph_growing_random_game(igraph_t *graph, integer_t n, 
			       integer_t m, bool_t directed,
			       bool_t citation) {

  long int no_of_nodes=n;
  long int no_of_neighbors=m;
  long int no_of_edges;
  vector_t edges;
  
  long int resp=0;

  long int i,j;

  no_of_edges=(no_of_nodes-1) * no_of_neighbors;
    
  vector_init(&edges, no_of_edges*2);  

  RNG_BEGIN();

  for (i=1; i<no_of_nodes; i++) {
    for (j=0; j<no_of_neighbors; j++) {
      if (citation) {
	long int to=RNG_INTEGER(0, i-1);
	VECTOR(edges)[resp++] = i;
	VECTOR(edges)[resp++] = to;
      } else {
	long int from=RNG_INTEGER(0, i);
	long int to=RNG_INTEGER(1,i);
	VECTOR(edges)[resp++] = from;
	VECTOR(edges)[resp++] = to;
      }
    }
  }

  RNG_END();

  igraph_create(graph, &edges, n, directed);
  vector_destroy(&edges);

  return 0;
}

int igraph_aging_prefatt_game(igraph_t *graph, integer_t n, integer_t m,
			      integer_t aging_type, real_t aging_exp) {
  /* TODO */
  return 0;
}
