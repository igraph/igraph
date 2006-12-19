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
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
   02110-1301 USA

*/

#include "igraph.h"
#include "memory.h"

#include <stdarg.h>

/** 
 * \section about_generators
 *
 * <para>Graph generators create graphs.</para>
 * 
 * <para>Almost all functions which create graph objects are documented
 * here. The exceptions are \ref igraph_subgraph() and alike, these 
 * create graphs based on another graph.</para>
 */


/**
 * \ingroup generators
 * \function igraph_create
 * \brief Creates a graph with the specified edges.
 * 
 * \param graph An uninitialized graph object.
 * \param edges The edges to add, the first two elements are the first
 *        edge, etc.
 * \param n The number of vertices in the graph, if smaller or equal
 *        to the highest vertex id in the \p edges vector it
 *        will be increased automatically. So it is safe to give 0
 *        here.
 * \param directed Boolean, whether to create a directed graph or
 *        not. If yes, then the first edge points from the first
 *        vertex id in \p edges to the second, etc.
 * \return Error code:
 *         \c IGRAPH_EINVEVECTOR: invalid edges
 *         vector (odd number of vertices).
 *         \c IGRAPH_EINVVID: invalid (negative)
 *         vertex id. 
 *
 * Time complexity: O(|V|+|E|),
 * |V| is the number of vertices,
 * |E| the number of edges in the
 * graph. 
 */
int igraph_create(igraph_t *graph, const igraph_vector_t *edges, igraph_integer_t n, 
		  igraph_bool_t directed) {
  igraph_real_t max=igraph_vector_max(edges)+1;

  if (igraph_vector_size(edges) % 2 != 0) {
    IGRAPH_ERROR("Invalid (odd) edges vector", IGRAPH_EINVEVECTOR);
  }
  if (!igraph_vector_isininterval(edges, 0, max-1)) {
    IGRAPH_ERROR("Invalid (negative) vertex id", IGRAPH_EINVVID);
  }

  IGRAPH_CHECK(igraph_empty(graph, n, directed));
  IGRAPH_FINALLY(igraph_destroy, graph);
  if (igraph_vector_size(edges)>0) {
    igraph_integer_t vc=igraph_vcount(graph);
    if (vc < max) {
      IGRAPH_CHECK(igraph_add_vertices(graph, max-vc, 0));
    }
    IGRAPH_CHECK(igraph_add_edges(graph, edges, 0));
  }

  IGRAPH_FINALLY_CLEAN(1);
  return 0;
}

int igraph_i_adjacency_directed(igraph_matrix_t *adjmatrix, igraph_vector_t *edges) {
  
  long int no_of_nodes=igraph_matrix_nrow(adjmatrix);
  long int i, j, k;
  
  for (i=0; i<no_of_nodes; i++) {
    for (j=0; j<no_of_nodes; j++) {
      long int M=MATRIX(*adjmatrix, i, j);
      for (k=0; k<M; k++) {
	IGRAPH_CHECK(igraph_vector_push_back(edges, i));
	IGRAPH_CHECK(igraph_vector_push_back(edges, j));
      }
    }
  }
  
  return 0;
}

int igraph_i_adjacency_max(igraph_matrix_t *adjmatrix, igraph_vector_t *edges) {
  
  long int no_of_nodes=igraph_matrix_nrow(adjmatrix);
  long int i, j, k;
  
  for (i=0; i<no_of_nodes; i++) {
    for (j=i; j<no_of_nodes; j++) {
      long int M1=MATRIX(*adjmatrix, i, j);
      long int M2=MATRIX(*adjmatrix, j, i);
      if (M1<M2) { M1=M2; }
      for (k=0; k<M1; k++) {
	IGRAPH_CHECK(igraph_vector_push_back(edges, i));
	IGRAPH_CHECK(igraph_vector_push_back(edges, j));
      }
    }
  }
  
  return 0;
}

int igraph_i_adjacency_upper(igraph_matrix_t *adjmatrix, igraph_vector_t *edges) {
  
  long int no_of_nodes=igraph_matrix_nrow(adjmatrix);
  long int i, j, k;
  
  for (i=0; i<no_of_nodes; i++) {
    for (j=i; j<no_of_nodes; j++) {
      long int M=MATRIX(*adjmatrix, i, j);
      for (k=0; k<M; k++) {
	IGRAPH_CHECK(igraph_vector_push_back(edges, i));
	IGRAPH_CHECK(igraph_vector_push_back(edges, j));
      }
    }
  }
  return 0;
}

int igraph_i_adjacency_lower(igraph_matrix_t *adjmatrix, igraph_vector_t *edges) {

  long int no_of_nodes=igraph_matrix_nrow(adjmatrix);
  long int i, j, k;
  
  for (i=0; i<no_of_nodes; i++) {
    for (j=0; j<=i; j++) {
      long int M=MATRIX(*adjmatrix, i, j);
      for (k=0; k<M; k++) {
	IGRAPH_CHECK(igraph_vector_push_back(edges, i));
	IGRAPH_CHECK(igraph_vector_push_back(edges, j));
      }
    }
  }
  return 0;
}

int igraph_i_adjacency_min(igraph_matrix_t *adjmatrix, igraph_vector_t *edges) {
  
  long int no_of_nodes=igraph_matrix_nrow(adjmatrix);
  long int i, j, k;
  
  for (i=0; i<no_of_nodes; i++) {
    for (j=i; j<no_of_nodes; j++) {
      long int M1=MATRIX(*adjmatrix, i, j);
      long int M2=MATRIX(*adjmatrix, j, i);
      if (M1>M2) { M1=M2; }
      for (k=0; k<M1; k++) {
	IGRAPH_CHECK(igraph_vector_push_back(edges, i));
	IGRAPH_CHECK(igraph_vector_push_back(edges, j));
      }
    }
  }
  
  return 0;
}

/**
 * \ingroup generators
 * \function igraph_adjacency
 * \brief Creates a graph object from an adjacency matrix.
 * 
 * \param graph Pointer to an uninitialized graph object.
 * \param adjmatrix The adjacency matrix. How it is interpreted
 *        depends on the \p mode argument.
 * \param mode Constant to specify how the given matrix is interpreted
 *        as an adjacency matrix. Possible values
 *        (A(i,j) 
 *        is the element in row i and column
 *        j in the adjacency matrix
 *        (\p adjmatrix): 
 *        \clist
 *        \cli IGRAPH_ADJ_DIRECTED
 *          the graph will be directed and
 *          an element gives the number of edges between two vertex.
 *        \cli IGRAPH_ADJ_UNDIRECTED
 *          this is the same as \c IGRAPH_ADJ_MAX,
 *          for convenience. 
 *        \cli IGRAPH_ADJ_MAX
 *          undirected graph will be created
 *          and the number of edges between vertex
 *          i and 
 *          j is
 *          max(A(i,j), A(j,i)). 
 *        \cli IGRAPH_ADJ_MIN
 *          undirected graph will be created
 *          with min(A(i,j), A(j,i))
 *          edges between vertex 
 *          i and
 *          j. 
 *        \cli IGRAPH_ADJ_PLUS 
 *          undirected graph will be created 
 *          with A(i,j)+A(j,i) edges
 *          between vertex 
 *          i and
 *          j.  
 *        \cli IGRAPH_ADJ_UPPER 
 *          undirected graph will be created,
 *          only the upper right triangle (including the diagonal) is
 *          used for the number of edges.
 *        \cli IGRAPH_ADJ_LOWER 
 *          undirected graph will be created,
 *          only the lower left triangle (including the diagonal) is 
 *          used for creating the edges.
 *       \endclist
 * \return Error code,
 *         \c IGRAPH_NONSQUARE: non-square matrix.
 * 
 * Time complexity: O(|V||V|+|E|),
 * |V| and 
 * |E| are number of vertices and
 * edges in the graph. 
 */

int igraph_adjacency(igraph_t *graph, igraph_matrix_t *adjmatrix,
		     igraph_adjacency_t mode) {

  igraph_vector_t edges=IGRAPH_VECTOR_NULL;
  long int no_of_nodes;

  /* Some checks */
  if (igraph_matrix_nrow(adjmatrix) != igraph_matrix_ncol(adjmatrix)) {
    IGRAPH_ERROR("Non-square matrix", IGRAPH_NONSQUARE);
  }

  IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
  
  /* Collect the edges */
  no_of_nodes=igraph_matrix_nrow(adjmatrix);
  switch (mode) {
  case IGRAPH_ADJ_DIRECTED:
    IGRAPH_CHECK(igraph_i_adjacency_directed(adjmatrix, &edges));
    break;
  case IGRAPH_ADJ_MAX:
    IGRAPH_CHECK(igraph_i_adjacency_max(adjmatrix, &edges));
    break;
  case IGRAPH_ADJ_UPPER:
    IGRAPH_CHECK(igraph_i_adjacency_upper(adjmatrix, &edges));
    break;
  case IGRAPH_ADJ_LOWER:
    IGRAPH_CHECK(igraph_i_adjacency_lower(adjmatrix, &edges));
    break;
  case IGRAPH_ADJ_MIN:
    IGRAPH_CHECK(igraph_i_adjacency_min(adjmatrix, &edges));
    break;
  case IGRAPH_ADJ_PLUS:
    IGRAPH_CHECK(igraph_i_adjacency_directed(adjmatrix, &edges));
    break;
  }

  IGRAPH_CHECK(igraph_create(graph, &edges, no_of_nodes, 
			     (mode == IGRAPH_ADJ_DIRECTED))); 
  igraph_vector_destroy(&edges);
  IGRAPH_FINALLY_CLEAN(1);
  
  return 0;
}

/**
 * \ingroup generators
 * \function igraph_star
 * \brief Creates a \em star graph, every vertex connects only to
 * the center.
 *
 * \param graph Pointer to an uninitialized graph object, this will
 *        be the result.
 * \param n Integer constant, the number of vertices in the graph.
 * \param mode Constant, gives the type of the star graph to
 *        create. Possible values:
 *        \clist
 *        \cli IGRAPH_STAR_OUT
 *          directed star graph, edges point
 *          \em from the center to the other vertices.
 *        \cli IGRAPH_STAR_IN
 *          directed star graph, edges point
 *          \em to the center from the other vertices.
 *        \cli IGRAPH_STAR_UNDIRECTED 
 *          an undirected star graph is
 *          created. 
 *        \endclist
 * \param center Id of the vertex which will be the center of the
 *          graph. 
 * \return Error code:
 *         \clist
 *         \cli IGRAPH_EINVVID 
 *           invalid number of vertices.
 *         \cli IGRAPH_EINVAL 
 *           invalid center vertex.
 *         \cli IGRAPH_EINVMODE 
 *           invalid mode argument.
 *         \endclist
 * 
 * Time complexity: O(|V|), the
 * number of vertices in the graph.
 *
 * \sa \ref igraph_lattice(), \ref igraph_ring(), \ref igraph_tree()
 * for creating other regular structures.
 */

int igraph_star(igraph_t *graph, igraph_integer_t n, igraph_star_mode_t mode, 
		igraph_integer_t center) {

  igraph_vector_t edges=IGRAPH_VECTOR_NULL;
  long int i;

  if (n<0) { 
    IGRAPH_ERROR("Invalid number of vertices", IGRAPH_EINVVID);
  }
  if (center<0 || center >n-1) {
    IGRAPH_ERROR("Invalid center vertex", IGRAPH_EINVAL);
  }
  if (mode != IGRAPH_STAR_OUT && mode != IGRAPH_STAR_IN && 
      mode != IGRAPH_STAR_UNDIRECTED) {
    IGRAPH_ERROR("invalid mode", IGRAPH_EINVMODE);
  }

  IGRAPH_VECTOR_INIT_FINALLY(&edges, (n-1)*2);
  
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
  
  IGRAPH_CHECK(igraph_create(graph, &edges, 0, 
			     (mode != IGRAPH_STAR_UNDIRECTED)));
  igraph_vector_destroy(&edges);
  IGRAPH_FINALLY_CLEAN(1);
  
  return 0;
}

/**
 * \ingroup generators 
 * \function igraph_lattice
 * \brief Creates most kind of lattices.
 *
 * \param graph An uninitialized graph object.
 * \param dimvector Vector giving the sizes of the lattice in each of
 *        its dimensions. Ie. the dimension of the lattice will be the
 *        same as the length of this vector.
 * \param nei Integer value giving the distance (number of steps)
 *        within which two vertices will be connected. Not implemented
 *        yet. 
 * \param directed Boolean, whether to create a directed graph. The
 *        direction of the edges is determined by the generation
 *        algorithm and is unlikely to suit you, so this isn't a very
 *        useful option.
 * \param mutual Boolean, if the graph is directed this gives whether
 *        to create all connections as mutual.
 * \param circular Boolean, defines whether the generated lattice is
 *        periodic.
 * \return Error code:
 *         \c IGRAPH_EINVAL: invalid (negative)
 *         dimension vector. 
 *
 * Time complexity: O(|V|+|E|) (as
 * far as i remember), |V| and
 * |E| are the number of vertices 
 * and edges in the generated graph.
 */
int igraph_lattice(igraph_t *graph, const igraph_vector_t *dimvector, igraph_integer_t nei, 
		   igraph_bool_t directed, igraph_bool_t mutual, igraph_bool_t circular) {

  long int dims=igraph_vector_size(dimvector);
  long int no_of_nodes=igraph_vector_prod(dimvector);
  igraph_vector_t edges=IGRAPH_VECTOR_NULL;
  long int *coords, *weights;
  long int i, j;
  int carry, pos;

  if (igraph_vector_any_smaller(dimvector, 0)) {
    IGRAPH_ERROR("Invalid dimension vector", IGRAPH_EINVAL);
  }

  /* init coords & weights */

  coords=Calloc(dims, long int);
  if (coords==0) {
    IGRAPH_ERROR("lattice failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(free, coords);	/* TODO: hack */
  weights=Calloc(dims, long int);
  if (weights == 0) {
    IGRAPH_ERROR("lattice failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(free, weights);
  weights[0]=1;
  for (i=1; i<dims; i++) {
    weights[i]=weights[i-1]*VECTOR(*dimvector)[i-1];
  }
  
  IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
  IGRAPH_CHECK(igraph_vector_reserve(&edges, no_of_nodes*dims +
				     mutual*directed * no_of_nodes*dims));

  for (i=0; i<no_of_nodes; i++) {
    IGRAPH_ALLOW_INTERRUPTION();
    for (j=0; j<dims; j++) {
      if (circular || coords[j] != VECTOR(*dimvector)[j]-1) {
	long int new_nei;
	if (coords[j] != VECTOR(*dimvector)[j]-1) {
	  new_nei = i + weights[j] + 1;
	} else {
	  new_nei = i - (VECTOR(*dimvector)[j]-1) * weights[j] + 1;
	}
	if (new_nei != i+1 && (VECTOR(*dimvector)[j] != 2 || coords[j] != 1)) {
	  igraph_vector_push_back(&edges, i); /* reserved */
	  igraph_vector_push_back(&edges, new_nei-1); /* reserved */
	}
      } /* if circular || coords[j] */
      if (mutual && directed && (circular || coords[j] != 0)) {
	long int new_nei;
	if (coords[j]!=0) {
	  new_nei=i-weights[j]+1;
	} else {
	  new_nei=i+(VECTOR(*dimvector)[j]-1) * weights[j]+1;
	}
	if (VECTOR(*dimvector)[j] != 2 || coords[j] != 0) {
	  igraph_vector_push_back(&edges, i); /* reserved */
	  igraph_vector_push_back(&edges, new_nei-1); /* reserved */
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

  IGRAPH_CHECK(igraph_create(graph, &edges, 0, directed));
  if (nei >= 2) {
    IGRAPH_CHECK(igraph_connect_neighborhood(graph, nei, IGRAPH_ALL));
  }

  /* clean up */
  Free(coords);
  Free(weights);
  igraph_vector_destroy(&edges);
  IGRAPH_FINALLY_CLEAN(3);

  return 0;
}

/**
 * \ingroup generators
 * \function igraph_ring
 * \brief Creates a \em ring graph, a one dimensional lattice.
 * 
 * \param graph Pointer to an uninitialized graph object.
 * \param n The number of vertices in the ring.
 * \param directed Logical, whether to create a directed ring.
 * \param mutual Logical, whether to create mutual edges in a directed
 *        ring. It is ignored for undirected graphs.
 * \param circular Logical, if false, the ring will be open (this is
 *        not a real \em ring actually).
 * \return Error code:
 *         \c IGRAPH_EINVAL: invalid number of vertices.
 * 
 * Time complexity: O(|V|), the
 * number of vertices in the graph.
 *
 * \sa \ref igraph_lattice() for generating more general lattices.
 */

int igraph_ring(igraph_t *graph, igraph_integer_t n, igraph_bool_t directed, igraph_bool_t mutual,
		igraph_bool_t circular) {
  
  igraph_vector_t v=IGRAPH_VECTOR_NULL;

  if (n<0) {
    IGRAPH_ERROR("negative number of vertices", IGRAPH_EINVAL);
  }
  
  IGRAPH_VECTOR_INIT_FINALLY(&v, 1);
  VECTOR(v)[0]=n;
  
  IGRAPH_CHECK(igraph_lattice(graph, &v, 1, directed, mutual, circular));
  igraph_vector_destroy(&v);		 
  
  IGRAPH_FINALLY_CLEAN(1);
  return 0;
}

/**
 * \ingroup generators
 * \function igraph_tree
 * \brief Creates a tree in which almost all vertices have the same
 * number of children.
 *
 * \param graph Pointer to an uninitialized graph object.
 * \param n Integer, the number of vertices in the graph.
 * \param children Integer, the number of children of a vertex in the
 *        tree. 
 * \param type Constant, gives whether to create a directed tree, and
 *        if this is the case, also its orientation. Possible values:
 *        \clist
 *        \cli IGRAPH_TREE_OUT 
 *          directed tree, the edges point
 *          from the parents to their children,
 *        \cli IGRAPH_TREE_IN 
 *          directed tree, the edges point from
 *          the children to their parents.
 *        \cli IGRAPH_TREE_UNDIRECTED
 *          undirected tree.
 *        \endclist
 * \return Error code:
 *         \c IGRAPH_EINVAL: invalid number of vertices.
 *         \c IGRAPH_INVMODE: invalid mode argument.
 * 
 * Time complexity: O(|V|+|E|), the
 * number of vertices plus the number of edges in the graph.
 * 
 * \sa \ref igraph_lattice(), \ref igraph_star() for creating other regular
 * structures. 
 */

int igraph_tree(igraph_t *graph, igraph_integer_t n, igraph_integer_t children, 
		igraph_tree_mode_t type) {
  
  igraph_vector_t edges=IGRAPH_VECTOR_NULL;
  long int i, j;
  long int idx=0;
  long int to=1;

  if (n<0 || children<=0) {
    IGRAPH_ERROR("Invalid number of vertices or children", IGRAPH_EINVAL);
  }
  if (type != IGRAPH_TREE_OUT && type != IGRAPH_TREE_IN &&
      type != IGRAPH_TREE_UNDIRECTED) {
    IGRAPH_ERROR("Invalid mode argument", IGRAPH_EINVMODE);
  }
  
  IGRAPH_VECTOR_INIT_FINALLY(&edges, 2*(n-1));
  
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
      
  IGRAPH_CHECK(igraph_create(graph, &edges, 0, type!=IGRAPH_TREE_UNDIRECTED));
  
  igraph_vector_destroy(&edges);
  IGRAPH_FINALLY_CLEAN(1);
  return 0;
}

/**
 * \ingroup generators
 * \function igraph_full
 * \brief Creates a full graph (directed or undirected, with or
 * without loops). 
 * 
 * </para><para>
 * In a full graph every possible edge is present, every vertex is
 * connected to every other vertex. 
 * 
 * \param graph Pointer to an uninitialized graph object.
 * \param n Integer, the number of vertices in the graph.
 * \param directed Logical, whether to create a directed graph.
 * \param loops Logical, whether to include self-edges (loops).
 * \return Error code:
 *         \c IGRAPH_EINVAL: invalid number of vertices.
 * 
 * Time complexity: O(|V|+|E|),
 * |V| is the number of vertices,
 * |E| the number of edges in the
 * graph. Of course this is the same as
 * O(|E|)=O(|V||V|) 
 * here. 
 * 
 * \sa \ref igraph_lattice(), \ref igraph_star(), \ref igraph_tree()
 * for creating other regular structures.
 */

int igraph_full(igraph_t *graph, igraph_integer_t n, igraph_bool_t directed, igraph_bool_t loops) {
  
  igraph_vector_t edges=IGRAPH_VECTOR_NULL;
  long int i, j;

  if (n<0) {
    IGRAPH_ERROR("invalid number of vertices", IGRAPH_EINVAL);
  }

  IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);

  if (directed && loops) {
    IGRAPH_CHECK(igraph_vector_reserve(&edges, n*n));
    for (i=0; i<n; i++) {
      for (j=0; j<n; j++) {
	igraph_vector_push_back(&edges, i); /* reserved */
	igraph_vector_push_back(&edges, j); /* reserved */
      }
    }
  } else if (directed && !loops) {
    IGRAPH_CHECK(igraph_vector_reserve(&edges, n*(n-1)));
    for (i=0; i<n; i++) {
      for (j=0; j<i; j++) {
	igraph_vector_push_back(&edges, i); /* reserved */
	igraph_vector_push_back(&edges, j); /* reserved */
      }
      for (j=i+1; j<n; j++) {
	igraph_vector_push_back(&edges, i); /* reserved */
	igraph_vector_push_back(&edges, j); /* reserved */
      }
    }
  } else if (!directed && loops) {
    IGRAPH_CHECK(igraph_vector_reserve(&edges, n*(n+1)/2));
    for (i=0; i<n; i++) {
      for (j=i; j<n; j++) {
	igraph_vector_push_back(&edges, i); /* reserved */
	igraph_vector_push_back(&edges, j); /* reserved */
      }
    }
  } else {
    IGRAPH_CHECK(igraph_vector_reserve(&edges, n*(n-1)/2));
    for (i=0; i<n; i++) {
      for (j=i+1; j<n; j++) {
	igraph_vector_push_back(&edges, i); /* reserved */
	igraph_vector_push_back(&edges, j); /* resreved */
      }
    }
  }
  
  IGRAPH_CHECK(igraph_create(graph, &edges, n, directed));
  igraph_vector_destroy(&edges);
  IGRAPH_FINALLY_CLEAN(1);
  
  return 0;
}

/**
 * \function igraph_small
 * \brief Shortand to create a short graph, giving the edges as agruments
 * 
 * </para><para>
 * This function is handy when a relatively small graph needs to be created. 
 * Instead giving the edges in vector, they are given simply as
 * arguments and a '-1' needs to be given after the last meaningful
 * edge argument. 
 * 
 * </para><para>Note that only graphs which have vertices less than
 * the highest value of the 'int' type can be created this way. If you
 * give larger values then the result is undefined.
 * 
 * \param graph Pointer to an uninitialized graph object, the result
 *        will be stored here.
 * \param n The number of vertices in the graph, an integer.
 * \param directed Logical constant, gives whether the graph should be
 *        directed. 
 * \param ... The additional arguments giving the edges of the
 *        graph. Don't forget to supply an additional '-1' after the last
 *        (meaningful) argument.
 * \return Error code.
 * 
 * Time complexity: O(|V|+|E|), the number of vertices plus the number
 * of edges in the graph to create.
 */

int igraph_small(igraph_t *graph, igraph_integer_t n, igraph_bool_t directed, 
		 ...) {
  igraph_vector_t edges;
  va_list ap;
  
  IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
  
  va_start(ap, directed);
  while (1) {
    int num = va_arg(ap, int);
    if (num == -1) {
      break;
    }
    igraph_vector_push_back(&edges, num);
  }

  IGRAPH_CHECK(igraph_create(graph, &edges, n, directed));
  
  igraph_vector_destroy(&edges);
  IGRAPH_FINALLY_CLEAN(1);
  return 0;
}

/**
 * \function igraph_extended_chordal_ring
 * Create an extended chordal ring
 * 
 * An extended chordal ring is regular graph, each node has the same
 * degree. It can be obtained from a simple ring by adding some extra
 * edges specified by a matrix. Let p denote the number of columns in
 * the <parameter>W</parameter> matrix. The extra edges of vertex i
 * are added according to column (i mod p) in
 * <parameter>W</parameter>. The number of extra edges is the number
 * of rows in <parameter>W</parameter>: for each row j an edge
 * i->i+w[ij] is added if i+w[ij] is less than the number of total
 * nodes. 
 * 
 * </para><para>
 * See also Kotsis, G: Interconnection Topologies for Parallel Processing
 * Systems, PARS Mitteilungen 11, 1-6, 1993.
 * 
 * \param graph Pointer to an uninitialized graph object, the result
 *   will be stored here. The result is always an undirected graph.
 * \param nodes Integer constant, the number of vertices in the
 *   graph. It must be at least 3.
 * \param W The matrix specifying the extra edges. The number of
 *   columns should divide the number of total vertices.
 * \return Error code.
 * 
 * \sa \ref igraph_ring().
 * 
 * Time complexity: O(|V|+|E|), the number of vertices plus the number
 * of edges.
 */

int igraph_extended_chordal_ring(igraph_t *graph, igraph_integer_t nodes, 
				 const igraph_matrix_t *W) {

  igraph_vector_t edges;
  long int period=igraph_matrix_ncol(W);
  long int degree=igraph_matrix_nrow(W)+2;
  long int i, j, mpos=0, epos=0;
  
  if (nodes<3) {
    IGRAPH_ERROR("An extended chordal ring has at least 3 nodes",
		 IGRAPH_EINVAL);
  }
  
  if ((long int)nodes % period != 0) {
    IGRAPH_ERROR("The period (number of columns in W) should divide the " 
		 "number of nodes", IGRAPH_EINVAL);
  }

  IGRAPH_VECTOR_INIT_FINALLY(&edges, nodes*degree);

  for (i=0; i<nodes-1; i++) {
    VECTOR(edges)[epos++] = i;
    VECTOR(edges)[epos++] = i+1;
  }
  VECTOR(edges)[epos++] = 0;
  VECTOR(edges)[epos++] = nodes-1;
  
  if (degree > 2) {
    for (i=0; i<nodes; i++) {
      for (j=0; j<degree-2; j++) {
	long int offset=MATRIX(*W, j, mpos);
	if (i+offset < nodes) {
	  VECTOR(edges)[epos++] = i;
	  VECTOR(edges)[epos++] = i+offset;
	}
      }
      mpos++; if (mpos==period) { mpos=0; }
    }
  }
  
  igraph_vector_resize(&edges, epos);
  IGRAPH_CHECK(igraph_create(graph, &edges, nodes, IGRAPH_UNDIRECTED));
  igraph_vector_destroy(&edges);
  IGRAPH_FINALLY_CLEAN(1);
  return 0;  
}

/**
 * \function igraph_connect_neighborhood
 * \brief Connects every vertex to its neighborhood
 */

int igraph_connect_neighborhood(igraph_t *graph, igraph_integer_t order,
				igraph_neimode_t mode) {
  
  long int no_of_nodes=igraph_vcount(graph);
  igraph_dqueue_t q;
  igraph_vector_t edges;
  long int i, j, in;
  long int *added;
  igraph_vector_t neis;
  
  if (order<0) {
    IGRAPH_ERROR("Negative order, cannot connect neighborhood", IGRAPH_EINVAL);
  }

  if (order<2) { 
    IGRAPH_WARNING("Order smaller than two, graph will be unchanged");
  }

  IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
  added=Calloc(no_of_nodes, long int);
  if (added==0) {
    IGRAPH_ERROR("Cannot connect neighborhood", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, added);
  IGRAPH_DQUEUE_INIT_FINALLY(&q, 100);
  IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
  
  for (i=0; i<no_of_nodes; i++) {
    added[i]=i+1;
    igraph_neighbors(graph, &neis, i, mode);
    in=igraph_vector_size(&neis);
    if (order > 1) {
      for (j=0; j<in; j++) {
	long int nei=VECTOR(neis)[j];
	added[nei]=i+1;
	igraph_dqueue_push(&q, nei);
	igraph_dqueue_push(&q, 1);
      }
    }
    
    while (!igraph_dqueue_empty(&q)) {
      long int actnode=igraph_dqueue_pop(&q);
      long int actdist=igraph_dqueue_pop(&q);
      long int n;
      igraph_neighbors(graph, &neis, actnode, mode);
      n=igraph_vector_size(&neis);
      
      if (actdist<order-1) {
	for (j=0; j<n; j++) {
	  long int nei=VECTOR(neis)[j];
	  if (added[nei] != i+1) {
	    added[nei]=i+1;
	    IGRAPH_CHECK(igraph_dqueue_push(&q, nei));
	    IGRAPH_CHECK(igraph_dqueue_push(&q, actdist+1));
	    if (mode == IGRAPH_OUT) {
	      IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
	      IGRAPH_CHECK(igraph_vector_push_back(&edges, nei));
	    } else {
	      IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
	      IGRAPH_CHECK(igraph_vector_push_back(&edges, nei));
	    }
	  }
	}
      } else { 
	for (j=0; j<n; j++) {
	  long int nei=VECTOR(neis)[j];
	  if (added[nei] != i+1) {
	    added[nei]=i+1;
	    if (mode == IGRAPH_OUT) {
	      IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
	      IGRAPH_CHECK(igraph_vector_push_back(&edges, nei));
	    } else {
	      IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
	      IGRAPH_CHECK(igraph_vector_push_back(&edges, nei));
	    }
	  }
	}
      }
      
    } /* while q not empty */
  } /* for i < no_of_nodes */
  
  igraph_vector_destroy(&neis);
  igraph_dqueue_destroy(&q);
  igraph_free(added);
  IGRAPH_FINALLY_CLEAN(3);
  
  IGRAPH_CHECK(igraph_add_edges(graph, &edges, 0));
  
  igraph_vector_destroy(&edges);
  IGRAPH_FINALLY_CLEAN(1);
  
  return 0;
}
