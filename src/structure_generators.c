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
 * \ingroup generators
 * \brief Creates a graph with the specified edges.
 * 
 * @param graph An uninitialized graph object.
 * @param edges The edges to add, the first two elements are the first
 *        edge, etc.
 * @param n The number of vertices in the graph, if smaller or equal
 *        to the highest vertex id in the <code>edges</code> vector it
 *        will be increased automatically. So it is safe to give 0
 *        here.
 * @param directed Boolean, whether to create a directed graph or
 *        not. If yes, then the first edge points from the first
 *        vertex id in <code>edges</code> to the second, etc.
 * @return Error code.
 *
 * Time complexity: <code>O(|V|+|E|)</code>, <code>|V|</code> is the
 * number of vertices, <code>|E|</code> the number of edges in the
 * graph.
 */
int igraph_create(igraph_t *graph, vector_t *edges, integer_t n, 
		  bool_t directed) {

  igraph_empty(graph, n, directed);
  if (vector_size(edges)>0) {
    real_t max=vector_max(edges)+1;
    integer_t vc=igraph_vcount(graph);
    if (vc < max) {
      igraph_add_vertices(graph, max-vc);
    }
    igraph_add_edges(graph, edges);
  }

  return 0;
}

int igraph_adjacency(igraph_t *graph, vector_t *adjmatrix,
		     integer_t dirmode) {
  /* TODO */
  return 0;
}

int igraph_star(igraph_t *graph, integer_t n, integer_t mode, 
		integer_t center, bool_t directed) {
  /* TODO */
  return 0;
}

int igraph_connect_neighborhood(igraph_t *graph, integer_t nei, 
				bool_t mutual) {
  /* TODO */
/* SEXP REST_connect_neighborhood(SEXP neis, SEXP pradius, SEXP pmutual) { */

/*   SEXP result; */
  
/*   long int no_of_nodes; */
/*   long int radius; */
/*   dqueue_t q; */
/*   long int *already_visited; */
/*   vector_t add; */
/*   long int i,j; */
/*   int mutual; */

/*   no_of_nodes=GET_LENGTH(neis); */
/*   radius=R(pradius); */
/*   mutual=LOGICAL(pmutual)[0]; */

/*   already_visited=(long int*) R_alloc(no_of_nodes, sizeof(long int)); */
/*   memset(already_visited, 0, no_of_nodes*sizeof(long int)); */

/*   dqueue_init(&q, 100); */
/*   vector_init(&add, 0); */
/*   vector_reserve(&add, radius*no_of_nodes); */
  
/*   for (i=1; i<=no_of_nodes; i++) { */
/*     dqueue_push(&q, i); */
/*     dqueue_push(&q, 0); */
/*     already_visited[i-1]=i; */
    
/*     while (!dqueue_empty(&q)) { */
/*       long int actnode=dqueue_pop(&q); */
/*       long int actdist=dqueue_pop(&q); */
/*       if (actdist >= 2 && (actnode > i || mutual)) { */
/* 	vector_push_back(&add, i); */
/* 	vector_push_back(&add, actnode); */
/*       } */
      
/*       if (actdist+1 <= radius) { */
/* 	for (j=0; j<GET_LENGTH(VECTOR_ELT(neis, actnode-1)); j++) { */
/* 	  long int neighbor=REAL(VECTOR_ELT(neis, actnode-1))[j]; */
/* 	  if (already_visited[neighbor-1] == i) { continue; } */
/* 	  already_visited[neighbor-1] = i; */
/* 	  dqueue_push(&q, neighbor); */
/* 	  dqueue_push(&q, actdist+1); */
/* 	} */
/*       } */
/*     } /\* while !dqueue_empty(q) *\/ */
/*     R_CheckUserInterrupt(); */
/*   } /\* for i<=no_of_nodes *\/ */

/*   dqueue_destroy(&q); */
  
/*   /\* Copy vector to result *\/ */
/*   j=vector_size(&add); */
/*   PROTECT(result=NEW_NUMERIC(j)); */
/*   for (i=0; i<j; i++) { */
/*     REAL(result)[i]=vector_e(&add, i); */
/*   } */

/*   /\* Clean and return *\/ */
/*   vector_destroy(&add); */
/*   UNPROTECT(1); */
  return 0;
}

/**
 * \ingroup generators 
 * \brief Creating all kinds of lattices.
 *
 * This function can create most kind of regular lattices.
 * @param graph An uninitialized graph object.
 * @param dimvector Vector giving the sizes of the lattice in each of
 *        its dimensions. IE. the dimension of the lattice will be the
 *        same as the length of this vector.
 * @param nei Integer value giving the distance (number of steps)
 *        within which two vertices will be connected. Not implemented
 *        yet. 
 * @param directed Boolean, whether to create a directed graph. The
 *        direction of the edges is determined by the generation
 *        algorithm and is unlikely to suit you, so this isn't a very
 *        useful option.
 * @param mutual Boolean, if the graph is directed this gives whether
 *        to create all connections as mutual.
 * @param circular Boolean, defines whether the generated lattice is
 *        periodic.
 * @return Error code.
 *
 * Time complexity: <code>O(|V|+|E|)</code> (as far as i remember),
 * <code>|V|</code> and <code>|E|</code> are the number of vertices
 * and edges in the generated graph.
 */
int igraph_lattice(igraph_t *graph, vector_t *dimvector, integer_t nei, 
		   bool_t directed, bool_t mutual, bool_t circular) {

  long int dims=vector_size(dimvector);
  long int no_of_nodes=vector_prod(dimvector);
  vector_t edges;
  long int *coords, *weights;
  long int i, j;
  long int resp=0;
  int carry, pos;

  /* init coords & weights */

  coords=Calloc(dims, long int);
  weights=Calloc(dims, long int);
  weights[0]=1;
  for (i=1; i<dims; i++) {
    weights[i]=weights[i-1]*VECTOR(*dimvector)[i-1];
  }
  
  vector_init(&edges, 0);
  vector_reserve(&edges, no_of_nodes*dims +
		 directed * no_of_nodes*dims);
  
  for (i=0; i<no_of_nodes; i++) {
    for (j=0; j<dims; j++) {
      if (circular || coords[j] != VECTOR(*dimvector)[j]-1) {
	long int new_nei;
	if (coords[j] != VECTOR(*dimvector)[j]-1) {
	  new_nei = i + weights[j] + 1;
	} else {
	  new_nei = i - (VECTOR(*dimvector)[j]-1) * weights[j] + 1;
	}
	if (new_nei != i+1 && (VECTOR(*dimvector)[j] != 2 || coords[j] != 1)) {
	  vector_push_back(&edges, i);
	  vector_push_back(&edges, new_nei-1);
	}
      } /* if circular || coords[j] */
      if (directed && (circular || coords[j] != 0)) {
	long int new_nei;
	if (coords[j]!=0) {
	  new_nei=i-weights[j]+1;
	} else {
	  new_nei=i+(VECTOR(*dimvector)[j]-1) * weights[j]+1;
	}
	if (VECTOR(*dimvector)[j] != 2 || coords[j] != 0) {
	  vector_push_back(&edges, i);
	  vector_push_back(&edges, new_nei-1);
	}
      } /* if circular || coords[0] */
    } /* for j<dims */
    
    /* increase coords */
    carry=1;
    pos=0;
    
    while (carry==1 && pos != dims) {
      if (coords[pos] != VECTOR(*dimvector)[pos]-1) {
	coords[pos]++;
	carry=0;
      } else {
	coords[pos]=0;
	pos++;
      }
    }
    
  } /* for i<no_of_nodes */

  igraph_create(graph, &edges, 0, directed);

  /* clean up */
  Free(coords);
  Free(weights);
  vector_destroy(&edges);

  return 0;
}

int igraph_ring(igraph_t *graph, integer_t n, bool_t directed, bool_t mutual) {
  /* TODO */
  return 0;
}

int igraph_tree(igraph_t *graph, integer_t n, integer_t children) {
  /* TODO */
  return 0;
}

/**********************************************************
 * Testing purposes                                       *
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
  
/*   vector_t v; */
/*   igraph_t g; */

/*   /\* igraph_create *\/ */
/*   vector_init(&v, 8); */
/*   VECTOR(v)[0]=0;  VECTOR(v)[1]=1; */
/*   VECTOR(v)[2]=1;  VECTOR(v)[3]=2; */
/*   VECTOR(v)[4]=2;  VECTOR(v)[5]=2; */
/*   VECTOR(v)[6]=2;  VECTOR(v)[7]=3; */
/*   igraph_create(&g, &v, 0, 1); */
  
/*   print_igraph(&g); */

/*   return 0; */
/* } */
