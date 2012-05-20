/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#include "igraph_operators.h"
#include "igraph_error.h"
#include "igraph_memory.h"
#include "igraph_interrupt_internal.h"
#include "igraph_interface.h"
#include "igraph_constructors.h"
#include "igraph_adjlist.h"
#include "igraph_attributes.h"
#include "config.h"

/**
 * \function igraph_disjoint_union
 * \brief Creates the union of two disjoint graphs
 *
 * </para><para>
 * First the vertices of the second graph will be relabeled with new
 * vertex ids to have two disjoint sets of vertex ids, then the union
 * of the two graphs will be formed.
 * If the two graphs have |V1| and |V2| vertices and |E1| and |E2|
 * edges respectively then the new graph will have |V1|+|V2| vertices
 * and |E1|+|E2| edges. 
 *
 * </para><para>
 * Both graphs need to have the same directedness, ie. either both
 * directed or both undirected.
 * 
 * </para><para>
 * The current version of this function cannot handle graph, vertex
 * and edge attributes, they will be lost.
 *
 * \param res  Pointer to an uninitialized graph object, the result
 *        will stored here.
 * \param left The first graph.
 * \param right The second graph.
 * \return Error code.
 * \sa \ref igraph_disjoint_union_many() for creating the disjoint union
 * of more than two graphs, \ref igraph_union() for non-disjoint
 * union.
 * 
 * Time complexity: O(|V1|+|V2|+|E1|+|E2|).
 * 
 * \example examples/simple/igraph_disjoint_union.c
 */

int igraph_disjoint_union(igraph_t *res, const igraph_t *left, 
			  const igraph_t *right) {

  long int no_of_nodes_left=igraph_vcount(left);
  long int no_of_nodes_right=igraph_vcount(right);
  long int no_of_edges_left=igraph_ecount(left);
  long int no_of_edges_right=igraph_ecount(right);
  igraph_vector_t edges;
  igraph_bool_t directed_left=igraph_is_directed(left);
  igraph_integer_t from, to;
  long int i;
  
  if (directed_left != igraph_is_directed(right)) {
    IGRAPH_ERROR("Cannot union directed and undirected graphs",
		 IGRAPH_EINVAL);
  }

  IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
  IGRAPH_CHECK(igraph_vector_reserve(&edges, 
				     2*(no_of_edges_left+no_of_edges_right)));
  for (i=0; i<no_of_edges_left; i++) {
    igraph_edge(left, i, &from, &to);
    igraph_vector_push_back(&edges, from);
    igraph_vector_push_back(&edges, to);
  }
  for (i=0; i<no_of_edges_right; i++) {
    igraph_edge(right, i, &from, &to);
    igraph_vector_push_back(&edges, from+no_of_nodes_left);
    igraph_vector_push_back(&edges, to+no_of_nodes_left);
  }
  
  IGRAPH_CHECK(igraph_create(res, &edges, no_of_nodes_left+no_of_nodes_right, 
			     directed_left));
  igraph_vector_destroy(&edges);
  IGRAPH_FINALLY_CLEAN(1);
  return 0;
}

/**
 * \function igraph_disjoint_union_many
 * \brief The disjint union of many graphs.
 *
 * </para><para>
 * First the vertices in the graphs will be relabeled with new vertex
 * ids to have pairwise disjoint vertex id sets and then the union of
 * the graphs is formed.
 * The number of vertices and edges in the result is the total number
 * of vertices and edges in the graphs.
 * 
 * </para><para>
 * Both graphs need to have the same directedness, ie. either both
 * directed or both undirected.
 * 
 * </para><para>
 * The current version of this function cannot handle graph, vertex
 * and edge attributes, they will be lost.
 *
 * \param res Pointer to an uninitialized graph object, the result of
 *        the operation will be stored here.
 * \param graphs Pointer vector, contains pointers to initialized
 *        graph objects.
 * \return Error code.
 * \sa \ref igraph_disjoint_union() for an easier syntax if you have
 * only two graphs, \ref igraph_union_many() for non-disjoint union.
 *
 * Time complexity: O(|V|+|E|), the number of vertices plus the number
 * of edges in the result.
 */

int igraph_disjoint_union_many(igraph_t *res, 
			       const igraph_vector_ptr_t *graphs) {
  long int no_of_graphs=igraph_vector_ptr_size(graphs);
  igraph_bool_t directed=1;
  igraph_vector_t edges;
  long int no_of_edges=0;
  long int shift=0;
  igraph_t *graph;
  long int i, j;
  igraph_integer_t from, to;
  
  if (no_of_graphs != 0) {
    graph=VECTOR(*graphs)[0];
    directed=igraph_is_directed(graph);
    for (i=0; i<no_of_graphs; i++) {      
      graph=VECTOR(*graphs)[i];
      no_of_edges += igraph_ecount(graph);
      if (directed != igraph_is_directed(graph)) {
	IGRAPH_ERROR("Cannot union directed and undirected graphs", 
		     IGRAPH_EINVAL);
      }
    }
  }
  
  IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
  IGRAPH_CHECK(igraph_vector_reserve(&edges, 2*no_of_edges));
  
  for (i=0; i<no_of_graphs; i++) {
    long int ec;
    graph=VECTOR(*graphs)[i];    
    ec=igraph_ecount(graph);
    for (j=0; j<ec; j++) {
      igraph_edge(graph, j, &from, &to);
      igraph_vector_push_back(&edges, from+shift);
      igraph_vector_push_back(&edges, to+shift);
    }
    shift += igraph_vcount(graph);
  }
  
  IGRAPH_CHECK(igraph_create(res, &edges, shift, directed));
  igraph_vector_destroy(&edges);
  IGRAPH_FINALLY_CLEAN(1);
  return 0;
}

/**
 * \function igraph_intersection
 * \brief Collect the common edges from two graphs.
 * 
 * </para><para>
 * The result graph contains only edges present both in the first and
 * the second graph. The number of vertices in the result graph is the
 * same as the larger from the two arguments.
 * 
 * \param res Pointer to an uninitialized graph object. This will
 * contain the result of the operation.
 * \param left The first operand, a graph object.
 * \param right The second operand, a graph object.
 * \return Error code.
 * \sa \ref igraph_intersection_many() to calculate the intersection
 * of many graphs at once, \ref igraph_union(), \ref
 * igraph_difference() for other operators.
 * 
 * Time complexity: O(|V|+|E|), |V| is the number of nodes, |E|
 * is the number of edges in the smaller graph of the two. (The one 
 * containing less vertices is considered smaller.)
 * 
 * \example examples/simple/igraph_intersection.c
 */

int igraph_intersection(igraph_t *res,
			const igraph_t *left, const igraph_t *right) {
  
  long int no_of_nodes_left=igraph_vcount(left);
  long int no_of_nodes_right=igraph_vcount(right);
  long int no_of_nodes;
  long int smaller_nodes;
  igraph_bool_t directed=igraph_is_directed(left);
  igraph_vector_t edges;
  igraph_vector_t nei1, nei2;
  long int i,j1,j2,n1,n2;
  igraph_integer_t v1, v2;

  if (directed != igraph_is_directed(right)) {
    IGRAPH_ERROR("Cannot intersect directed and undirected graph",
		 IGRAPH_EINVAL);
  }

  IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&nei1, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&nei2, 0);

  no_of_nodes= no_of_nodes_left > no_of_nodes_right ? 
    no_of_nodes_left : no_of_nodes_right;
  smaller_nodes= no_of_nodes_left < no_of_nodes_right ? 
    no_of_nodes_left : no_of_nodes_right;

  for (i=0; i<smaller_nodes; i++) {
    IGRAPH_ALLOW_INTERRUPTION();
    IGRAPH_CHECK(igraph_neighbors(left, &nei1, i, IGRAPH_OUT));
    IGRAPH_CHECK(igraph_neighbors(right, &nei2, i, IGRAPH_OUT));
    if (!directed) {
      igraph_vector_filter_smaller(&nei1, i);
      igraph_vector_filter_smaller(&nei2, i);
    }
    n1=igraph_vector_size(&nei1);
    n2=igraph_vector_size(&nei2);
    j1=j2=0;
    while (j1<n1 && j2<n2) {
      v1=VECTOR(nei1)[j1]; v2=VECTOR(nei2)[j2];
      if (v1 < v2) {
	j1++;
      } else if (v1 > v2) {
	j2++; 
      } else {
	IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
	IGRAPH_CHECK(igraph_vector_push_back(&edges, v1));
	j1++; 
	j2++;
      }
    }
  }

  IGRAPH_CHECK(igraph_create(res, &edges, no_of_nodes, directed));
  igraph_vector_destroy(&edges);
  igraph_vector_destroy(&nei1);
  igraph_vector_destroy(&nei2);
  IGRAPH_FINALLY_CLEAN(3);
  return 0;
}

void igraph_i_intersection_many_free(igraph_vector_ptr_t *v) {
  long int i, n=igraph_vector_ptr_size(v);
  for (i=0; i<n; i++) { 
    if (VECTOR(*v)[i] != 0) {
      igraph_vector_destroy(VECTOR(*v)[i]);
      igraph_Free(VECTOR(*v)[i]);
    }
  }
  igraph_vector_ptr_destroy(v);
}

/**
 * \function igraph_intersection_many
 * \brief The intersection of more than two graphs.
 * 
 * </para><para> 
 * This function calculates the intersection of the graphs stored in
 * the \c graphs argument. Only those edges will be included in the
 * result graph which are part of every graph in \c graphs.
 *
 * </para><para>
 * The number of vertices in the result graph will be the maximum
 * number of vertices in the argument graphs.
 *
 * \param res Pointer to an uninitialized graph object, the result of
 *        the operation will be stored here.
 * \param graphs Pointer vector, contains pointers to graphs objects,
 *        the operands of the intersection operator.
 * \return Error code.
 * \sa \ref igraph_intersection() for the intersection of two graphs, 
 * \ref igraph_union_many(), \ref igraph_union() and \ref
 * igraph_difference() for other operators.
 * 
 * Time complexity: O(|V|+|E|), |V| is the number of vertices,
 * |E| is the number of edges in the smallest graph (ie. the graph having
 * the less vertices).
 */

int igraph_intersection_many(igraph_t *res, 
			     const igraph_vector_ptr_t *graphs) {

  long int no_of_graphs=igraph_vector_ptr_size(graphs);
  long int no_of_nodes=0;
  long int smallest_nodes=0;
  igraph_bool_t directed=1;
  igraph_vector_t edges;
  igraph_vector_ptr_t neivects;
  igraph_vector_t neiptr;
  long int i, j;
  
  /* Check directedness */
  if (no_of_graphs != 0) {
    directed=igraph_is_directed(VECTOR(*graphs)[0]);
    no_of_nodes=smallest_nodes=igraph_vcount(VECTOR(*graphs)[0]);
  }
  for (i=1; i<no_of_graphs; i++) {
    if (directed != igraph_is_directed(VECTOR(*graphs)[i])) {
      IGRAPH_ERROR("Cannot intersect directed and undirected graphs",
		   IGRAPH_EINVAL);
    }
  }

  /* Calculate number of nodes */
  for (i=0; i<no_of_graphs; i++) {
    long int n=igraph_vcount(VECTOR(*graphs)[i]);
    if (n < smallest_nodes) {
      smallest_nodes=n;
    } else if (n > no_of_nodes) {
      no_of_nodes=n;
    }
  }

  IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&neiptr, no_of_graphs);

  /* Init neighbor vectors */
  if (no_of_graphs != 0) {
    IGRAPH_CHECK(igraph_vector_ptr_init(&neivects, no_of_graphs));
    IGRAPH_FINALLY(igraph_i_intersection_many_free, &neivects);
  }
  for (i=0; i<no_of_graphs; i++) {
    VECTOR(neivects)[i]=igraph_Calloc(1, igraph_vector_t);
    if (VECTOR(neivects)[i]==0) { 
      IGRAPH_ERROR("Cannot intersect graphs", IGRAPH_ENOMEM);
    }
    IGRAPH_CHECK(igraph_vector_init(VECTOR(neivects)[i], 0));
  }

  /* Main part */
  for (i=0; i<smallest_nodes; i++) {
    igraph_bool_t l;
    
    IGRAPH_ALLOW_INTERRUPTION();

    /* get neighbors */
    for (j=0; j<no_of_graphs; j++) {
      IGRAPH_CHECK(igraph_neighbors(VECTOR(*graphs)[j], VECTOR(neivects)[j], i,
				    IGRAPH_OUT));
      if (!directed) {
	igraph_vector_filter_smaller(VECTOR(neivects)[j], i);
      }
    }
    igraph_vector_null(&neiptr);

    /* check if there are more edges */
    l=1; j=0;
    while (l && j<no_of_graphs) {
      l = VECTOR(neiptr)[j] < igraph_vector_size(VECTOR(neivects)[j]);
      j++;
    }

    while (l) {
      
      /* get largest head element and decide whether all head elements
	 are the same */
      igraph_bool_t k=1;
      long int head=VECTOR(*(igraph_vector_t*)VECTOR(neivects)[0])
	[(long int)VECTOR(neiptr)[0]];
      long int maxhead=head;
      for (j=1; j<no_of_graphs; j++) {
	long int h=VECTOR(*(igraph_vector_t*)VECTOR(neivects)[j])
	  [(long int)VECTOR(neiptr)[j]];
	k= (k && head==h);
	if (h > maxhead) {
	  maxhead=h;
	}
      }
      
      /* add edge if common, remove head elements */
      if (k) {
	IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
	IGRAPH_CHECK(igraph_vector_push_back(&edges, head));
	for (j=0; j<no_of_graphs; j++) {
	  VECTOR(neiptr)[j] += 1;
	}
      } else {
	for (j=0; j<no_of_graphs; j++) {
	  while (VECTOR(neiptr)[j] < igraph_vector_size(VECTOR(neivects)[j]) &&
		 VECTOR(*(igraph_vector_t*)VECTOR(neivects)[j])
		 [ (long int) VECTOR(neiptr)[j]] < maxhead){
	    VECTOR(neiptr)[j] += 1;
	  }
	}
      }

      /* update l, check if there are more edges */
      l=1; j=0;
      while (l && j<no_of_graphs) {
	l = VECTOR(neiptr)[j] < igraph_vector_size(VECTOR(neivects)[j]);
	j++;
      }
      
    }
  }
      
  if (no_of_graphs != 0) {
    for (i=0; i<no_of_graphs; i++) {
      igraph_vector_destroy(VECTOR(neivects)[i]);
      igraph_Free(VECTOR(neivects)[i]);
    }
    igraph_vector_ptr_destroy(&neivects);
    IGRAPH_FINALLY_CLEAN(1);
  }

  IGRAPH_CHECK(igraph_create(res, &edges, no_of_nodes, directed));
  igraph_vector_destroy(&edges);
  igraph_vector_destroy(&neiptr);
  IGRAPH_FINALLY_CLEAN(2);
  return 0;
}

/**
 * \function igraph_union
 * \brief Calculates the union of two graphs.
 * 
 * </para><para>
 * The number of vertices in the result is that of the larger graph
 * from the two arguments. The result graph contains edges which are
 * present in at least one of the operand graphs.
 * 
 * \param res Pointer to an uninitialized graph object, the result
 *        will be stored here.
 * \param left The first graph.
 * \param right The second graph.
 * \return Error code.
 * \sa \ref igraph_union_many() for the union of many graphs, 
 * \ref igraph_intersection() and \ref igraph_difference() for other
 * operators. 
 * 
 * Time complexity: O(|V|+|E|), |V| is the number of
 * vertices, |E| the number of edges in the result graph.
 * 
 * \example examples/simple/igraph_union.c
 */

int igraph_union(igraph_t *res, 
		 const igraph_t *left, const igraph_t *right) {
  
  long int no_of_nodes_left=igraph_vcount(left);
  long int no_of_nodes_right=igraph_vcount(right);
  long int no_of_nodes;
  igraph_bool_t directed=igraph_is_directed(left);
  igraph_vector_t edges;
  igraph_vector_t nei1, nei2;
  long int i;

  if (directed != igraph_is_directed(right)) {
    IGRAPH_ERROR("Cannot make union of directed and undirected graph",
		 IGRAPH_EINVAL);
  }
  
  IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&nei1, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&nei2, 0);

  no_of_nodes=no_of_nodes_left > no_of_nodes_right ?
    no_of_nodes_left : no_of_nodes_right;
  
  for (i=0; i<no_of_nodes; i++) {
    long int n1=0, n2=0, p1=0, p2=0;
    IGRAPH_ALLOW_INTERRUPTION();
    if (i<no_of_nodes_left) {
      IGRAPH_CHECK(igraph_neighbors(left, &nei1, i, IGRAPH_OUT));
      if (!directed) {
	igraph_vector_filter_smaller(&nei1, i);
      }
      n1=igraph_vector_size(&nei1);
    }
    if (i<no_of_nodes_right) {
      IGRAPH_CHECK(igraph_neighbors(right, &nei2, i, IGRAPH_OUT));
      if (!directed) {
	igraph_vector_filter_smaller(&nei2, i);
      }
      n2=igraph_vector_size(&nei2);
    }
    
    while (p1<n1 || p2<n2) {
      if (p2>=n2 || (p1<n1 && VECTOR(nei1)[p1] < VECTOR(nei2)[p2])) {
	IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
	IGRAPH_CHECK(igraph_vector_push_back(&edges, VECTOR(nei1)[p1]));
	p1++;
      } else if (p1>=n1 || (p2<n2 && VECTOR(nei1)[p1] > VECTOR(nei2)[p2])) {
	IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
	IGRAPH_CHECK(igraph_vector_push_back(&edges, VECTOR(nei2)[p2]));
	p2++;
      } else {
	IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
	IGRAPH_CHECK(igraph_vector_push_back(&edges, VECTOR(nei1)[p1]));
	p1++;
	p2++;
      }
    }
  }
  
  IGRAPH_CHECK(igraph_create(res, &edges, no_of_nodes, directed));
  igraph_vector_destroy(&edges);
  igraph_vector_destroy(&nei1);
  igraph_vector_destroy(&nei2);
  IGRAPH_FINALLY_CLEAN(3);
  return 0;
}
	
void igraph_i_union_many_free(igraph_vector_ptr_t *v) {
  long int i, n=igraph_vector_ptr_size(v);
  for (i=0; i<n; i++) { 
    if (VECTOR(*v)[i] != 0) {
      igraph_vector_destroy(VECTOR(*v)[i]);
      igraph_Free(VECTOR(*v)[i]);
    }
  }
  igraph_vector_ptr_destroy(v);
}

/**
 * \function igraph_union_many
 * \brief Creates the union of many graphs.
 * 
 * </para><para> 
 * The result graph will contain as many vertices as the largest graph
 * among the arguments does, and an edge will be included in it if it
 * is part of at least one operand graph.
 * 
 * </para><para>
 * The directedness of the operand graphs must be the same.
 * 
 * \param res Pointer to an uninitialized graph object, this will
 *        contain the result.
 * \param graphs Pointer vector, contains pointers to the operands of
 *        the union operator, graph objects of course.
 * \return Error code.
 * \sa \ref igraph_union() for the union of two graphs, \ref
 * igraph_intersection_many(), \ref igraph_intersection() and \ref
 * igraph_difference for other operators.
 * 
 * 
 * Time complexity: O(|V|+|E|), |V| is the number of vertices
 * in largest graph and |E| is the number of edges in the result graph.
 * 
 * \example examples/simple/igraph_union.c
 */

int igraph_union_many(igraph_t *res, const igraph_vector_ptr_t *graphs) {
  
  long int no_of_graphs=igraph_vector_ptr_size(graphs);
  long int no_of_nodes=0;
  igraph_bool_t directed=1;
  igraph_vector_t edges;
  igraph_vector_ptr_t neivects;
  long int i, j;
  
  /* Check directedness */
  if (no_of_graphs != 0) {
    directed=igraph_is_directed(VECTOR(*graphs)[0]);
    no_of_nodes=igraph_vcount(VECTOR(*graphs)[0]);
  }
  for (i=1; i<no_of_graphs; i++) {
    if (directed != igraph_is_directed(VECTOR(*graphs)[i])) {
      IGRAPH_ERROR("Cannot union directed and undirected graphs",
		   IGRAPH_EINVAL);
    }
  }

  /* Calculate number of nodes */
  for (i=0; i<no_of_graphs; i++) {
    long int n=igraph_vcount(VECTOR(*graphs)[i]);
    if (n > no_of_nodes) {
      no_of_nodes=n;
    }
  }

  IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);

  /* Init neighbor vectors */
  if (no_of_graphs != 0) {
    IGRAPH_CHECK(igraph_vector_ptr_init(&neivects, no_of_graphs));
    IGRAPH_FINALLY(igraph_i_union_many_free, &neivects);
  }
  for (i=0; i<no_of_graphs; i++) {
    VECTOR(neivects)[i]=igraph_Calloc(1, igraph_vector_t);
    if (VECTOR(neivects)[i]==0) { 
      IGRAPH_ERROR("Cannot union graphs", IGRAPH_ENOMEM);
    }
    IGRAPH_CHECK(igraph_vector_init(VECTOR(neivects)[i], 0));
  }
  
  /* Main part */
  for (i=0; i<no_of_nodes; i++) {
    igraph_bool_t l=0;
    long int bigtail;
    
    IGRAPH_ALLOW_INTERRUPTION();

    /* get neighbors */
    for (j=0; j<no_of_graphs; j++) {
      if (i<igraph_vcount(VECTOR(*graphs)[j])) {
	IGRAPH_CHECK(igraph_neighbors(VECTOR(*graphs)[j], VECTOR(neivects)[j], 
				      i, IGRAPH_OUT));
	if (!directed) {
	  igraph_vector_filter_smaller(VECTOR(neivects)[j], i);
	}
	if (!igraph_vector_empty(VECTOR(neivects)[j])) {
	  l=1;
	}
      }
    }

    while (l) {

      /* Get the largest tail element */
      l=0; bigtail=0;
      for (j=0; j<no_of_graphs; j++) {
	if (!l && !igraph_vector_empty(VECTOR(neivects)[j])) {
	  l=1; 
	  bigtail=igraph_vector_tail(VECTOR(neivects)[j]);
	} else if (l && !igraph_vector_empty(VECTOR(neivects)[j])) {
	  long int tail=igraph_vector_tail(VECTOR(neivects)[j]);
	  if (tail > bigtail) {
	    bigtail=tail;
	  }
	}
      }
      
      /* add the edge */
      IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
      IGRAPH_CHECK(igraph_vector_push_back(&edges, bigtail));
      
      /* update neighbors */
      for (j=0; j<no_of_graphs; j++) {
	if (!igraph_vector_empty(VECTOR(neivects)[j]) &&
	    igraph_vector_tail(VECTOR(neivects)[j])==bigtail) {
	  igraph_vector_pop_back(VECTOR(neivects)[j]);
	}
      }
      
      /* Check if there are more edges */
      l=0; j=0;
      while (!l && j<no_of_graphs) {
	l=!igraph_vector_empty(VECTOR(neivects)[j]);
	j++;
      }
    }
  }
  
  if (no_of_graphs != 0) {
    for (i=0; i<no_of_graphs; i++) {
      igraph_vector_destroy(VECTOR(neivects)[i]);
      igraph_Free(VECTOR(neivects)[i]);
    }
    igraph_vector_ptr_destroy(&neivects);
    IGRAPH_FINALLY_CLEAN(1);
  }

  IGRAPH_CHECK(igraph_create(res, &edges, no_of_nodes, directed));
  igraph_vector_destroy(&edges);
  IGRAPH_FINALLY_CLEAN(1);
  return 0;
}

/**
 * \function igraph_difference
 * \brief Calculate the difference of two graphs
 * 
 * </para><para>
 * The number of vertices in the result is the number of vertices in
 * the original graph, ie. the left, first operand. In the results
 * graph only edges will be included from \c orig which are not
 * present in \c sub.
 * 
 * \param res Pointer to an uninitialized graph object, the result
 * will be stored here.
 * \param orig The left operand of the operator, a graph object.
 * \param sub The right operand of the operator, a graph object.
 * \return Error code.
 * \sa \ref igraph_intersection() and \ref igraph_union() for other
 * operators.
 * 
 * Time complexity: O(|V|+|E|), |V| is the number vertices in
 * the smaller graph, |E| is the
 * number of edges in the result graph.
 * 
 * \example examples/simple/igraph_difference.c
 */

int igraph_difference(igraph_t *res, 
		      const igraph_t *orig, const igraph_t *sub) {

  /* Quite nasty, but we will use that an edge adjacency list
     contains the vertices according to the order of the 
     vertex ids at the "other" end of the edge. */

  long int no_of_nodes_orig=igraph_vcount(orig);
  long int no_of_nodes_sub =igraph_vcount(sub);
  long int no_of_nodes=no_of_nodes_orig;
  long int smaller_nodes;
  igraph_bool_t directed=igraph_is_directed(orig);
  igraph_vector_t edges;
  igraph_vector_t edge_ids;
  igraph_vector_t *nei1, *nei2;
  igraph_inclist_t inc_orig, inc_sub;
  long int i;
  igraph_integer_t v1, v2;

  if (directed != igraph_is_directed(sub)) {
    IGRAPH_ERROR("Cannot subtract directed and undirected graphs",
		 IGRAPH_EINVAL);
  }
  
  IGRAPH_VECTOR_INIT_FINALLY(&edge_ids, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
  IGRAPH_CHECK(igraph_inclist_init(orig, &inc_orig, IGRAPH_OUT));
  IGRAPH_FINALLY(igraph_inclist_destroy, &inc_orig);
  IGRAPH_CHECK(igraph_inclist_init(sub, &inc_sub, IGRAPH_OUT));
  IGRAPH_FINALLY(igraph_inclist_destroy, &inc_sub);
  
  smaller_nodes=no_of_nodes_orig > no_of_nodes_sub ?
    no_of_nodes_sub : no_of_nodes_orig;
  
  for (i=0; i<smaller_nodes; i++) {
    long int n1, n2, e1, e2;
    IGRAPH_ALLOW_INTERRUPTION();
    nei1=igraph_inclist_get(&inc_orig, i);
    nei2=igraph_inclist_get(&inc_sub, i);
    n1=igraph_vector_size(nei1)-1;
    n2=igraph_vector_size(nei2)-1;
    while (n1>=0 && n2>=0) {
      e1=VECTOR(*nei1)[n1];
      e2=VECTOR(*nei2)[n2];
      v1=IGRAPH_OTHER(orig, e1, i);
      v2=IGRAPH_OTHER(sub, e2, i);
      
      if (!directed && v1<i) { 
	n1--;
      } else if (!directed && v2<i) {
	n2--;
      } else if (v1>v2) {
	IGRAPH_CHECK(igraph_vector_push_back(&edge_ids, e1));
	IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
	IGRAPH_CHECK(igraph_vector_push_back(&edges, v1));
	n1--;	
      } else if (v2>v1) {
	n2--;
      } else {
	n1--;
	n2--;
      }
    }
    
    /* Copy remaining edges */
    while (n1>=0) {
      e1=VECTOR(*nei1)[n1];
      v1=IGRAPH_OTHER(orig, e1, i);
      if (directed || v1 >= i) { 
	IGRAPH_CHECK(igraph_vector_push_back(&edge_ids, e1));
	IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
	IGRAPH_CHECK(igraph_vector_push_back(&edges, v1));
      }
      n1--;
    }
  }

  /* copy remaining edges, use the previous value of 'i' */
  for (; i<no_of_nodes_orig; i++) {
    long int n1, e1;
    nei1=igraph_inclist_get(&inc_orig, i);
    n1=igraph_vector_size(nei1)-1;
    while (n1>=0) {
      e1=VECTOR(*nei1)[n1];
      v1=IGRAPH_OTHER(orig, e1, i);
      if (directed || v1 >= i) { 
	IGRAPH_CHECK(igraph_vector_push_back(&edge_ids, e1));
	IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
	IGRAPH_CHECK(igraph_vector_push_back(&edges, v1));
      }
      n1--;
    }
  }

  igraph_inclist_destroy(&inc_sub);
  igraph_inclist_destroy(&inc_orig);
  IGRAPH_FINALLY_CLEAN(2);
  IGRAPH_CHECK(igraph_create(res, &edges, no_of_nodes, directed));
  igraph_vector_destroy(&edges);  
  IGRAPH_FINALLY_CLEAN(1);

  /* Attributes */
  if (orig->attr) {
    IGRAPH_I_ATTRIBUTE_DESTROY(res);
    IGRAPH_I_ATTRIBUTE_COPY(res, orig, /*graph=*/1, /*vertex=*/1, /*edge=*/0);
    IGRAPH_CHECK(igraph_i_attribute_permute_edges(orig, res, &edge_ids));
  }
  
  igraph_vector_destroy(&edge_ids);
  IGRAPH_FINALLY_CLEAN(1);
  
  return 0;
}

/**
 * \function igraph_complementer
 * \brief Create the complementer of a graph
 * 
 * </para><para>The complementer graph means that all edges which are
 * not part of the original graph will be included in the result.
 * 
 * \param res Pointer to an uninitialized graph object.
 * \param graph The original graph.
 * \param loops Whether to add loop edges to the complementer graph. 
 * \return Error code.
 * \sa \ref igraph_union(), \ref igraph_intersection() and \ref
 * igraph_difference().
 * 
 * Time complexity: O(|V|+|E1|+|E2|), |V| is the number of
 * vertices in the graph, |E1| is the number of edges in the original
 * and |E2| in the complementer graph.
 * 
 * \example examples/simple/igraph_complementer.c
 */

int igraph_complementer(igraph_t *res, const igraph_t *graph, 
			igraph_bool_t loops) {

  long int no_of_nodes=igraph_vcount(graph);
  igraph_vector_t edges;
  igraph_vector_t neis;
  long int i, j;
  long int zero=0, *limit;

  IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);

  if (igraph_is_directed(graph)) {
    limit=&zero;
  } else {
    limit=&i;
  }
  
  for (i=0; i<no_of_nodes; i++) {
    IGRAPH_ALLOW_INTERRUPTION();
    IGRAPH_CHECK(igraph_neighbors(graph, &neis, i, IGRAPH_OUT));
    if (loops) {
      for (j=no_of_nodes-1; j>=*limit; j--) {
	if (igraph_vector_empty(&neis) || j>igraph_vector_tail(&neis)) {
	  IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
	  IGRAPH_CHECK(igraph_vector_push_back(&edges, j));
	} else {
	  igraph_vector_pop_back(&neis);
	}
      }
    } else {
      for (j=no_of_nodes-1; j>=*limit; j--) {
	if (igraph_vector_empty(&neis) || j>igraph_vector_tail(&neis)) {
	  if (i!=j) {
	    IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
	    IGRAPH_CHECK(igraph_vector_push_back(&edges, j));
	  }
	} else {
	  igraph_vector_pop_back(&neis);
	}
      }
    }      
  }
  
  IGRAPH_CHECK(igraph_create(res, &edges, no_of_nodes, 
			     igraph_is_directed(graph)));  
  igraph_vector_destroy(&edges);
  igraph_vector_destroy(&neis);
  IGRAPH_I_ATTRIBUTE_DESTROY(res);
  IGRAPH_I_ATTRIBUTE_COPY(res, graph, /*graph=*/1, /*vertex=*/1, /*edge=*/0);
  IGRAPH_FINALLY_CLEAN(2);
  return 0;
}

/**
 * \function igraph_compose
 * \brief Calculates the composition of two graphs
 * 
 * The composition of graphs contains the same number of vertices as
 * the bigger graph of the two operands. It contains an (i,j) edge if
 * and only if there is a k vertex, such that the first graphs
 * contains an (i,k) edge and the second graph a (k,j) edge. 
 * 
 * </para><para>This is of course exactly the composition of two
 * binary relations.
 *
 * </para><para>Two two graphs must have the same directedness,
 * otherwise the function returns with an error message.
 * Note that for undirected graphs the two relations are by definition
 * symmetric. 
 * 
 * \param res Pointer to an uninitialized graph object, the result
 *        will be stored here.
 * \param g1 The firs operand, a graph object.
 * \param g2 The second operand, another graph object.
 * \return Error code.
 * 
 * Time complexity: O(|V|*d1*d2), |V| is the number of vertices in the
 * first graph, d1 and d2 the average degree in the first and second
 * graphs. 
 * 
 * \example examples/simple/igraph_compose.c
 */

int igraph_compose(igraph_t *res, const igraph_t *g1, const igraph_t *g2) {
  
  long int no_of_nodes_left=igraph_vcount(g1);
  long int no_of_nodes_right=igraph_vcount(g2);
  long int no_of_nodes;
  igraph_bool_t directed=igraph_is_directed(g1);
  igraph_vector_t edges;
  igraph_vector_t neis1, neis2;
  long int i;

  if (directed != igraph_is_directed(g2)) {
    IGRAPH_ERROR("Cannot compose directed and undirected graph",
		 IGRAPH_EINVAL);
  }

  no_of_nodes= no_of_nodes_left > no_of_nodes_right ? 
    no_of_nodes_left : no_of_nodes_right;

  IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&neis1, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&neis2, 0);
  
  for (i=0; i<no_of_nodes_left; i++) {
    IGRAPH_ALLOW_INTERRUPTION();
    IGRAPH_CHECK(igraph_neighbors(g1, &neis1, i, IGRAPH_OUT));
    while (!igraph_vector_empty(&neis1)) {
      long int con=igraph_vector_pop_back(&neis1);
      if (con<no_of_nodes_right) {
	IGRAPH_CHECK(igraph_neighbors(g2, &neis2, con, IGRAPH_OUT));
      }
      while (!igraph_vector_empty(&neis2)) {
	IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
	IGRAPH_CHECK(igraph_vector_push_back(&edges,
					     igraph_vector_pop_back(&neis2)));
      }
    }
  }

  igraph_vector_destroy(&neis1);
  igraph_vector_destroy(&neis2);
  IGRAPH_FINALLY_CLEAN(2);

  IGRAPH_CHECK(igraph_create(res, &edges, no_of_nodes, directed));

  igraph_vector_destroy(&edges);
  IGRAPH_FINALLY_CLEAN(1);
  return 0;
}
