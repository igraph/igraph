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
  real_t max=vector_max(edges)+1;

  if (vector_size(edges) % 2 != 0) {
    IGRAPH_ERROR("Invalid (odd) edges vector", IGRAPH_EINVAL);
  }
  if (!vector_isininterval(edges, 0, max-1)) {
    IGRAPH_ERROR("Invalid (negative) vertex id", IGRAPH_EINVAL);
  }

  igraph_empty(graph, n, directed);
  if (vector_size(edges)>0) {
    integer_t vc=igraph_vcount(graph);
    if (vc < max) {
      igraph_add_vertices(graph, max-vc);
    }
    igraph_add_edges(graph, edges);
  }

  return 0;
}

int igraph_i_adjacency_directed(matrix_t *adjmatrix, vector_t *edges) {
  
  long int no_of_nodes=matrix_nrow(adjmatrix);
  long int i, j, k;
  
  for (i=0; i<no_of_nodes; i++) {
    for (j=0; j<no_of_nodes; j++) {
      long int M=MATRIX(*adjmatrix, i, j);
      for (k=0; k<M; k++) {
	vector_push_back(edges, i);
	vector_push_back(edges, j);
      }
    }
  }
  
  return 0;
}

int igraph_i_adjacency_max(matrix_t *adjmatrix, vector_t *edges) {
  
  long int no_of_nodes=matrix_nrow(adjmatrix);
  long int i, j, k;
  
  for (i=0; i<no_of_nodes; i++) {
    for (j=i; j<no_of_nodes; j++) {
      long int M1=MATRIX(*adjmatrix, i, j);
      long int M2=MATRIX(*adjmatrix, j, i);
      if (M1<M2) { M1=M2; }
      for (k=0; k<M1; k++) {
	vector_push_back(edges, i);
	vector_push_back(edges, j);
      }
    }
  }
  
  return 0;
}

int igraph_i_adjacency_upper(matrix_t *adjmatrix, vector_t *edges) {
  
  long int no_of_nodes=matrix_nrow(adjmatrix);
  long int i, j, k;
  
  for (i=0; i<no_of_nodes; i++) {
    for (j=i; j<no_of_nodes; j++) {
      long int M=MATRIX(*adjmatrix, i, j);
      for (k=0; k<M; k++) {
	vector_push_back(edges, i);
	vector_push_back(edges, j);
      }
    }
  }
  return 0;
}

int igraph_i_adjacency_lower(matrix_t *adjmatrix, vector_t *edges) {

  long int no_of_nodes=matrix_nrow(adjmatrix);
  long int i, j, k;
  
  for (i=0; i<no_of_nodes; i++) {
    for (j=0; j<=i; j++) {
      long int M=MATRIX(*adjmatrix, i, j);
      for (k=0; k<M; k++) {
	vector_push_back(edges, i);
	vector_push_back(edges, j);
      }
    }
  }
  return 0;
}

int igraph_i_adjacency_min(matrix_t *adjmatrix, vector_t *edges) {
  
  long int no_of_nodes=matrix_nrow(adjmatrix);
  long int i, j, k;
  
  for (i=0; i<no_of_nodes; i++) {
    for (j=i; j<no_of_nodes; j++) {
      long int M1=MATRIX(*adjmatrix, i, j);
      long int M2=MATRIX(*adjmatrix, j, i);
      if (M1>M2) { M1=M2; }
      for (k=0; k<M1; k++) {
	vector_push_back(edges, i);
	vector_push_back(edges, j);
      }
    }
  }
  
  return 0;
}

/**
 * \ingroup generators
 * \brief Creates a graph object from an adjacency matrix.
 * 
 * @param graph Pointer to an uninitialized graph object.
 * @param adjmatrix The adjacency matrix. How it is interpreted
 *        depends on the <code>mode</code> argument.
 * @param mode Constant to specify how the given matrix is interpreted
 *        as an adjacency matrix. Possible values (<code>A(i,j)</code>
 *        is the element in row <code>i</code> and column
 *        <code>j</code> in the adjacency matrix
 *        (<code>adjmatrix</code>): 
 *        - <b>IGRAPH_ADJ_DIRECTED</b>, the graph will be directed and
 *          an element gives the number of edges between two vertex.
 *        - <b>IGRAPH_ADJ_UNDIRECTED</b>, this is the same as
 *          <b>IGRAPH_ADJ_MAX</b>, for convenience.
 *        - <b>IGRAPH_ADJ_MAX</b>, undirected graph will be created
 *          and the number of edges between vertex <code>i</code> and
 *          <code>j</code> is <code>max(A(i,j), A(j,i))</code>.
 *        - <b>IGRAPH_ADJ_MIN</b>, undirected graph will be created
 *          with <code>min(A(i,j), A(j,i))</code> edges between vertex
 *          <code>i</code> and <code>j</code>.
 *        - <b>IGRAPH_ADJ_PLUS</b>, undirected graph will be created 
 *          with <code>A(i,j)+A(j,i)</code> edges between vertex
 *          <code>i</code> and <code>j</code>. 
 *        - <b>IGRAPH_ADJ_UPPER</b>, undirected graph will be created,
 *          only the upper right triangle (including the diagonal) is
 *          used for the number of edges.
 *        - <b>IGRAPH_ADJ_LOWER</b>, undirected graph will be created,
 *          only the lower left triangle (including the diagonal) is 
 *          used for creating the edges.
 * @return Error code.
 * 
 * Time complexity: <code>O(|V||V|+|E|)</code>, <code>|V|</code> and
 * <code>|E|</code> are number of vertices and edges in the graph.
 */

int igraph_adjacency(igraph_t *graph, matrix_t *adjmatrix,
		     igraph_adjacency_t mode) {

  vector_t edges;
  long int no_of_nodes;
  int retval=0;

  /* Some checks */
  if (matrix_nrow(adjmatrix) != matrix_ncol(adjmatrix)) {
    IGRAPH_ERROR("Non-square matrix", IGRAPH_EINVAL);
  }

  vector_init(&edges, 0);
  
  /* Collect the edges */
  no_of_nodes=matrix_nrow(adjmatrix);
  switch (mode) {
  case IGRAPH_ADJ_DIRECTED:
    retval=igraph_i_adjacency_directed(adjmatrix, &edges);
    break;
  case IGRAPH_ADJ_MAX:
    retval=igraph_i_adjacency_max(adjmatrix, &edges);
    break;
  case IGRAPH_ADJ_UPPER:
    retval=igraph_i_adjacency_upper(adjmatrix, &edges);
    break;
  case IGRAPH_ADJ_LOWER:
    retval=igraph_i_adjacency_lower(adjmatrix, &edges);
    break;
  case IGRAPH_ADJ_MIN:
    retval=igraph_i_adjacency_min(adjmatrix, &edges);
    break;
  case IGRAPH_ADJ_PLUS:
    retval=igraph_i_adjacency_directed(adjmatrix, &edges);
    break;
  }

  igraph_create(graph, &edges, no_of_nodes, (mode == IGRAPH_ADJ_DIRECTED));
  vector_destroy(&edges);
  
  return retval;
}

/**
 * \ingroup generators
 * \brief Created a \a star graph, every vertex connect to the center
 * only.
 *
 * @param graph Pointer to an uninitialized graph object, this will
 *        be the result.
 * @param n Integer constant, the number of vertices in the graph.
 * @param mode Contant, gives the type of the star graph to
 *        create. Possible values:
 *        - <b>IGRAPH_STAR_OUT</b>, directed star graph, edges point
 *          <em>from</em> the center to the other vertices.
 *        - <b>IGRAPH_STAR_IN</b>, directed star graph, edges point
 *          <em>to</em> the center from the other vertices.
 *        - <b>IGRAPH_STAR_UNDIRECTED</b>, an undirected star graph is
 *          created. 
 * @param center Id of the vertex which will be the center of the
 *          graph. 
 * @return Error code.
 *
 * Time complexity: <code>O(|V|)</code>, the number of vertices in the
 * graph.
 *
 * \sa igraph_lattice(), igraph_ring(), igraph_tree() for creating
 * other regular structures.
 */

int igraph_star(igraph_t *graph, integer_t n, igraph_star_mode_t mode, 
		integer_t center) {

  vector_t edges;
  long int i;

  if (n<0) { 
    IGRAPH_ERROR("Invalid number of vertices", IGRAPH_EINVAL);
  }
  if (center<0 || center >n-1) {
    IGRAPH_ERROR("Invalid center vertex", IGRAPH_EINVAL);
  }
  if (mode != IGRAPH_STAR_OUT && mode != IGRAPH_STAR_IN && 
      mode != IGRAPH_STAR_UNDIRECTED) {
    IGRAPH_ERROR("invalid mode", IGRAPH_EINVAL);
  }

  vector_init(&edges, (n-1)*2);
  if (mode == IGRAPH_STAR_OUT) {
    for (i=0; i<center; i++) {
      VECTOR(edges)[2*i]=center;
      VECTOR(edges)[2*i+1]=i;
    }
    for (i=center+1; i<n; i++) {
      VECTOR(edges)[2*(i-1)]=center;
      VECTOR(edges)[2*(i-1)+1]=i;
    }
  } else {
    for (i=0; i<center; i++) {
      VECTOR(edges)[2*i+1]=center;
      VECTOR(edges)[2*i]=i;
    }
    for (i=center+1; i<n; i++) {
      VECTOR(edges)[2*(i-1)+1]=center;
      VECTOR(edges)[2*(i-1)]=i;
    }
  }
  
  igraph_create(graph, &edges, 0, (mode != IGRAPH_STAR_UNDIRECTED));
  vector_destroy(&edges);
  
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
  int carry, pos;

  if (vector_any_smaller(dimvector, 0)) {
    IGRAPH_ERROR("Invalid dimension vector", IGRAPH_EINVAL);
  }

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

/**
 * \ingroup generators
 * \brief Creates a \a ring graph, a one dimensional lattice.
 * 
 * @param graph Pointer to an uninitialized graph object.
 * @param n The number of vertices in the ring.
 * @param directed Logical, whether to create a directed ring.
 * @param mutual Logical, whether to create mutual edges in a directed
 *        ring. It is ignored for undirected graphs.
 * @param circular Logical, if false, the ring will be open (this is
 *        not a real <em>ring</em> actually).
 * @return Error code.
 * 
 * Time complexity: <code>O(|V|)</code>, the number of vertices in the
 * graph.
 *
 * \sa igraph_lattice() for generating more general lattices.
 */

int igraph_ring(igraph_t *graph, integer_t n, bool_t directed, bool_t mutual,
		bool_t circular) {
  
  vector_t v;

  if (n<0) {
    IGRAPH_ERROR("negative number of vertices", IGRAPH_EINVAL);
  }
  
  vector_init(&v, 1);
  VECTOR(v)[0]=n;
  
  igraph_lattice(graph, &v, 1, directed, mutual, circular);
  vector_destroy(&v);		 
  
  return 0;
}

/**
 * \ingroup generators
 * \brief Creates a tree in which almost all vertices has the same
 * number of children.
 *
 * @param graph Pointer to an uninitialized graph object.
 * @param n Integer, the number of vertices in the graph.
 * @param children Integer, the number of children of a vertex in the
 *        tree. 
 * @param type Constant, gives whether to create a directed tree, and
 *        if this is the case, also its orientation. Possible values:
 *        - <b>IGRAPH_TREE_OUT</b>, directed tree, the edges point
 *          from the parents to their children,
 *        - <b>IGRAPH_TREE_IN</b>, directed tree, the edges point from
 *          the children to their parents.
 *        - <b>IGRAPH_TREE_UNDIRECTED</b>, undirected tree.
 * @return Error code.
 * 
 * Time complexity: <code>O(|V|+|E|)</code>, the number of vertices
 * plus the number of edges in the graph.
 * 
 * \sa igraph_lattice(), igraph_star() for creating other regular
 * structures. 
 */

int igraph_tree(igraph_t *graph, integer_t n, integer_t children, 
		igraph_tree_mode_t type) {
  
  vector_t edges;
  long int i, j;
  long int idx=0;
  long int to=1;

  if (n<0 || children<=0) {
    IGRAPH_ERROR("Invalid number of vertices or children", IGRAPH_EINVAL);
  }
  if (type != IGRAPH_TREE_OUT && type != IGRAPH_TREE_IN !=
      type != IGRAPH_TREE_UNDIRECTED) {
    IGRAPH_ERROR("Invalid mode argument", IGRAPH_EINVAL);
  }
  
  vector_init(&edges, 2*(n-1));
  
  i=0;
  if (type == IGRAPH_TREE_OUT) {
    while (idx<2*(n-1)) {
      for (j=0; j<children && idx<2*(n-1); j++) {
	VECTOR(edges)[idx++]=i;
	VECTOR(edges)[idx++]=to++;
      }
      i++;
    }
  } else {
    while (idx<2*(n-1)) {
      for (j=0; j<children && idx<2*(n-1); j++) {
	VECTOR(edges)[idx++]=to++;
	VECTOR(edges)[idx++]=i;
      }
      i++;
    }
  }
      
  igraph_create(graph, &edges, 0, type!=IGRAPH_TREE_UNDIRECTED);
  
  vector_destroy(&edges);
  return 0;
}

/**
 * \ingroup generators
 * \brief Creates a full graph (directed or undirected, with or
 * without loops). 
 * 
 * In a full graph every possible edge is present, every vertex is
 * connected to every other vertex. 
 * 
 * @param graph Pointer to an uninitialized graph object.
 * @param n Integer, the number of vertices in the graph.
 * @param directed Logical, whether to create a directed graph.
 * @param loops Logical, whether to include self-edges (loops).
 * @return Error code.
 * 
 * Time complexity: <code>O(|V|+|E|)</code>, <code>|V|</code> is the
 * number of vertices, <code>|E|</code> the number of edges in the
 * graph. Of course this is the same as <code>O(|E|)=O(|V||V|)</code>
 * here. 
 * 
 * \sa igraph.lattice(), igraph.star(), igraph.tree() for creating
 * other regular structures.
 */

int igraph_full(igraph_t *graph, integer_t n, bool_t directed, bool_t loops) {
  
  vector_t edges;
  long int i, j;

  if (n<0) {
    IGRAPH_ERROR("invalid number of vertices", IGRAPH_EINVAL);
  }

  vector_init(&edges, 0);
  if (directed && loops) {
    vector_reserve(&edges, n*n);
    for (i=0; i<n; i++) {
      for (j=0; j<n; j++) {
	vector_push_back(&edges, i);
	vector_push_back(&edges, j);
      }
    }
  } else if (directed && !loops) {
    vector_reserve(&edges, n*(n-1));
    for (i=0; i<n; i++) {
      for (j=0; j<i; j++) {
	vector_push_back(&edges, i);
	vector_push_back(&edges, j);
      }
      for (j=i+1; j<n; j++) {
	vector_push_back(&edges, i);
	vector_push_back(&edges, j);
      }
    }
  } else if (!directed && loops) {
    vector_reserve(&edges, n*(n+1)/2);
    for (i=0; i<n; i++) {
      for (j=i; j<n; j++) {
	vector_push_back(&edges, i);
	vector_push_back(&edges, j);
      }
    }
  } else {
    vector_reserve(&edges, n*(n-1)/2);
    for (i=0; i<n; i++) {
      for (j=i+1; j<n; j++) {
	vector_push_back(&edges, i);
	vector_push_back(&edges, j);
      }
    }
  }
  
  igraph_create(graph, &edges, n, directed);
  vector_destroy(&edges);
  
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
