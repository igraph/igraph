/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2013  Gabor Csardi <csardi.gabor@gmail.com>
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
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA 
   02110-1301 USA

*/

#include "igraph_scan.h"
#include "igraph_interface.h"
#include "igraph_adjlist.h"
#include "igraph_memory.h"
#include "igraph_interrupt_internal.h"
#include "igraph_arpack.h"
#include "igraph_eigen.h"
#include "igraph_centrality.h"
#include "igraph_operators.h"

int igraph_local_scan_0(const igraph_t *graph, igraph_vector_t *res,
			const igraph_vector_t *weights,
			igraph_neimode_t mode) {
  if (weights) {
    igraph_strength(graph, res, igraph_vss_all(), mode, /*loops=*/ 1,
		    weights);
  } else {
    igraph_degree(graph, res, igraph_vss_all(), mode, /*loops=*/ 1);
  }
  return 0;
}

/* From triangles.c */

int igraph_i_trans4_al_simplify(igraph_adjlist_t *al,
				const igraph_vector_int_t *rank);

/* This removes loop, multiple edges and edges that point
   "backwards" according to the rank vector. It works on
   edge lists */

int igraph_i_trans4_il_simplify(const igraph_t *graph, igraph_inclist_t *il,
				const igraph_vector_int_t *rank) {

  long int i;
  long int n=il->length;
  igraph_vector_int_t mark;
  igraph_vector_int_init(&mark, n);
  IGRAPH_FINALLY(igraph_vector_int_destroy, &mark);

  for (i=0; i<n; i++) {
    igraph_vector_int_t *v=&il->incs[i];
    int j, l=igraph_vector_int_size(v);
    int irank=VECTOR(*rank)[i];
    VECTOR(mark)[i] = i+1;
    for (j=0; j<l; /* nothing */) {
      long int edge=(long int) VECTOR(*v)[j];
      long int e=IGRAPH_OTHER(graph, edge, i);
      if (VECTOR(*rank)[e] > irank && VECTOR(mark)[e] != i+1) {
	VECTOR(mark)[e]=i+1;
	j++;
      } else {
	VECTOR(*v)[j] = igraph_vector_int_tail(v);
	igraph_vector_int_pop_back(v);
	l--;
      }
    }
  }

  igraph_vector_int_destroy(&mark);
  IGRAPH_FINALLY_CLEAN(1);
  return 0;

}

/* This one handles both weighted and unweighted cases */

int igraph_i_local_scan_1_directed(const igraph_t *graph,
				   igraph_vector_t *res,
				   const igraph_vector_t *weights,
				   igraph_neimode_t mode) {

  int no_of_nodes=igraph_vcount(graph);
  igraph_inclist_t incs;
  int i, node;

  igraph_vector_int_t neis;

  IGRAPH_CHECK(igraph_inclist_init(graph, &incs, mode));
  IGRAPH_FINALLY(igraph_inclist_destroy, &incs);

  igraph_vector_int_init(&neis, no_of_nodes);
  IGRAPH_FINALLY(igraph_vector_int_destroy, &neis);

  igraph_vector_resize(res, no_of_nodes);
  igraph_vector_null(res);

  for (node=0; node < no_of_nodes; node++) {
    igraph_vector_int_t *edges1=igraph_inclist_get(&incs, node);
    int edgeslen1=igraph_vector_int_size(edges1);

    IGRAPH_ALLOW_INTERRUPTION();

    /* Mark neighbors and self*/
    VECTOR(neis)[node] = node+1;
    for (i=0; i<edgeslen1; i++) {
      int e=VECTOR(*edges1)[i];
      int nei=IGRAPH_OTHER(graph, e, node);
      igraph_real_t w= weights ? VECTOR(*weights)[e] : 1;
      VECTOR(neis)[nei] = node+1;
      VECTOR(*res)[node] += w;
    }

    /* Crawl neighbors */
    for (i=0; i<edgeslen1; i++) {
      int e2=VECTOR(*edges1)[i];
      int nei=IGRAPH_OTHER(graph, e2, node);
      igraph_vector_int_t *edges2=igraph_inclist_get(&incs, nei);
      int j, edgeslen2=igraph_vector_int_size(edges2);
      for (j=0; j<edgeslen2; j++) {
	int e2=VECTOR(*edges2)[j];
	int nei2=IGRAPH_OTHER(graph, e2, nei);
	igraph_real_t w2= weights ? VECTOR(*weights)[e2] : 1;
	if (VECTOR(neis)[nei2] == node+1) {
	  VECTOR(*res)[node] += w2;
	}
      }
    }

  } /* node < no_of_nodes */

  igraph_vector_int_destroy(&neis);
  igraph_inclist_destroy(&incs);
  IGRAPH_FINALLY_CLEAN(2);

  return 0;
}

int igraph_i_local_scan_1_sumweights(const igraph_t *graph,
				     igraph_vector_t *res,
				     const igraph_vector_t *weights) {

  long int no_of_nodes=igraph_vcount(graph);
  long int node, i, j, nn;
  igraph_inclist_t allinc;
  igraph_vector_int_t *neis1, *neis2;
  long int neilen1, neilen2;
  long int *neis;
  long int maxdegree;

  igraph_vector_int_t order;
  igraph_vector_int_t rank;
  igraph_vector_t degree, *edge1=&degree; /* reuse degree as edge1 */

  if (igraph_vector_size(weights) != igraph_ecount(graph)) {
    IGRAPH_ERROR("Invalid weight vector length", IGRAPH_EINVAL);
  }

  igraph_vector_int_init(&order, no_of_nodes);
  IGRAPH_FINALLY(igraph_vector_int_destroy, &order);
  IGRAPH_VECTOR_INIT_FINALLY(&degree, no_of_nodes);

  IGRAPH_CHECK(igraph_degree(graph, &degree, igraph_vss_all(), IGRAPH_ALL,
			     IGRAPH_LOOPS));
  maxdegree=(long int) igraph_vector_max(&degree)+1;
  igraph_vector_order1_int(&degree, &order, maxdegree);
  igraph_vector_int_init(&rank, no_of_nodes);
  IGRAPH_FINALLY(igraph_vector_int_destroy, &rank);
  for (i=0; i<no_of_nodes; i++) {
    VECTOR(rank)[ VECTOR(order)[i] ] = no_of_nodes-i-1;
  }

  IGRAPH_CHECK(igraph_inclist_init(graph, &allinc, IGRAPH_ALL));
  IGRAPH_FINALLY(igraph_inclist_destroy, &allinc);
  IGRAPH_CHECK(igraph_i_trans4_il_simplify(graph, &allinc, &rank));

  neis=igraph_Calloc(no_of_nodes, long int);
  if (neis==0) {
    IGRAPH_ERROR("undirected local transitivity failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, neis);

  IGRAPH_CHECK(igraph_strength(graph, res, igraph_vss_all(), IGRAPH_ALL,
			       IGRAPH_LOOPS, weights));

  for (nn=no_of_nodes-1; nn>=0; nn--) {
    node=VECTOR(order)[nn];

    IGRAPH_ALLOW_INTERRUPTION();

    neis1=igraph_inclist_get(&allinc, node);
    neilen1=igraph_vector_int_size(neis1);

    /* Mark the neighbors of the node */
    for (i=0; i<neilen1; i++) {
      int edge = VECTOR(*neis1)[i];
      int nei = IGRAPH_OTHER(graph, edge, node);
      VECTOR(*edge1)[nei] = VECTOR(*weights)[edge];
      neis[nei] = node+1;
    }

    for (i=0; i<neilen1; i++) {
      long int edge=VECTOR(*neis1)[i];
      long int nei=IGRAPH_OTHER(graph, edge, node);
      igraph_real_t w=VECTOR(*weights)[edge];
      neis2=igraph_inclist_get(&allinc, nei);
      neilen2=igraph_vector_int_size(neis2);
      for (j=0; j<neilen2; j++) {
	long int edge2=VECTOR(*neis2)[j];
	long int nei2=IGRAPH_OTHER(graph, edge2, nei);
	igraph_real_t w2=VECTOR(*weights)[edge2];
	if (neis[nei2] == node+1) {
	  VECTOR(*res)[node] += w2;
	  VECTOR(*res)[nei2] += w;
	  VECTOR(*res)[nei] += VECTOR(*edge1)[nei2];
	}
      }
    }
  }

  igraph_free(neis);
  igraph_inclist_destroy(&allinc);
  igraph_vector_int_destroy(&rank);
  igraph_vector_destroy(&degree);
  igraph_vector_int_destroy(&order);
  IGRAPH_FINALLY_CLEAN(5);

  return 0;
}

int igraph_local_scan_1_ecount(const igraph_t *graph, igraph_vector_t *res,
			       const igraph_vector_t *weights,
			       igraph_neimode_t mode) {

  if ( (mode == IGRAPH_OUT || mode == IGRAPH_IN) &&
       igraph_is_directed(graph) ) {
    return igraph_i_local_scan_1_directed(graph, res, weights, mode);
  } else {
    if (weights) {
      return igraph_i_local_scan_1_sumweights(graph, res, weights);
    } else {

#define TRIEDGES
#include "triangles_template.h"
#undef TRIEDGES

    }
  }

  return 0;
}

int igraph_local_scan_1_ecount_approximate(const igraph_t *graph,
					   igraph_vector_t *res, int noevals,
					   igraph_arpack_options_t *options) {
  igraph_vector_t values;
  igraph_matrix_t vectors;
  igraph_eigen_which_t which;
  int no_of_nodes=igraph_vcount(graph);
  int i, j;

  if (noevals <= 0) {
    IGRAPH_ERROR("Number of eigenvalues should be positive",
                 IGRAPH_EINVAL);
  }
  if (noevals >= no_of_nodes) {
    IGRAPH_ERROR("Number of eigenvalues should be less than number "
                 "or vertices", IGRAPH_EINVAL);
  }

  IGRAPH_CHECK(igraph_vector_resize(res, no_of_nodes));
  igraph_vector_null(res);

  if (igraph_ecount(graph) == 0) { return 0; }

  igraph_vector_init(&values, 0);
  IGRAPH_FINALLY(igraph_vector_destroy, &values);
  igraph_matrix_init(&vectors, 0, 0);
  IGRAPH_FINALLY(igraph_matrix_destroy, &vectors);

  which.pos=IGRAPH_EIGEN_LM;
  which.howmany=noevals;

  igraph_eigen_adjacency(graph, IGRAPH_EIGEN_ARPACK, &which, options,
                         /*storage=*/ 0, &values, &vectors,
                         /*cmplxvalues=*/ 0, /*cmplxvectors=*/ 0);

  for (j=0; j<noevals; j++) {
    igraph_real_t v=VECTOR(values)[j];
    v = v * v * v;
    for (i=0; i<no_of_nodes; i++) {
      VECTOR(*res)[i] += v * MATRIX(vectors, i, j) * MATRIX(vectors, i, j);
    }
  }
  igraph_vector_scale(res, .5);

  igraph_matrix_destroy(&vectors);
  IGRAPH_FINALLY_CLEAN(1);

  /* Plus the degree */
  igraph_degree(graph, &values, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS);
  igraph_vector_add(res, &values);

  igraph_vector_destroy(&values);
  IGRAPH_FINALLY_CLEAN(1);

  return 0;
}

int igraph_local_scan_1_ecount_approximate_eigen(
				 const igraph_t *graph,
				 igraph_vector_t *res,
				 const igraph_vector_t *values,
				 const igraph_matrix_t *vectors) {

  igraph_vector_t degree;
  int no_of_nodes=igraph_vcount(graph);
  int noevals=igraph_vector_size(values);
  int i, j;

  if (igraph_matrix_nrow(vectors) != no_of_nodes) {
    IGRAPH_ERROR("Invalid eigenvector matrix, wrong number of rows",
                 IGRAPH_EINVAL);
  }
  if (igraph_matrix_ncol(vectors) != noevals) {
    IGRAPH_ERROR("Invalid eigenvector matrix, wrong number of columns",
                 IGRAPH_EINVAL);
  }

  IGRAPH_CHECK(igraph_vector_resize(res, no_of_nodes));
  igraph_vector_null(res);

  for (j=0; j<noevals; j++) {
    igraph_real_t v=VECTOR(*values)[j];
    v = v * v * v;
    for (i=0; i<no_of_nodes; i++) {
      VECTOR(*res)[i] += v * MATRIX(*vectors, i, j) *
        MATRIX(*vectors, i, j);
    }
  }
  igraph_vector_scale(res, .5);

  igraph_vector_init(&degree, no_of_nodes);
  IGRAPH_FINALLY(igraph_vector_destroy, &degree);

  igraph_degree(graph, &degree, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS);
  igraph_vector_add(res, &degree);

  igraph_vector_destroy(&degree);
  IGRAPH_FINALLY_CLEAN(1);

  return 0;
}

int igraph_i_local_scan_0_them_w(const igraph_t *us, const igraph_t *them,
			     igraph_vector_t *res,
			     const igraph_vector_t *weights_us,
			     const igraph_vector_t *weights_them,
			     igraph_neimode_t mode) {

  igraph_t is;
  igraph_vector_t map2;
  int i, m;

  if (!weights_us || !weights_them) {
    IGRAPH_ERROR("Edge weights not given for weighted scan-0",
		 IGRAPH_EINVAL);
  }
  if (igraph_vector_size(weights_us) != igraph_ecount(us)) {
    IGRAPH_ERROR("Invalid weights length (us) for scan-0", IGRAPH_EINVAL);
  }
  if (igraph_vector_size(weights_them) != igraph_ecount(them)) {
    IGRAPH_ERROR("Invalid weights length (them) for scan-0", IGRAPH_EINVAL);
  }

  IGRAPH_VECTOR_INIT_FINALLY(&map2, 0);
  igraph_intersection(&is, us, them, /*map1=*/ 0, &map2);
  IGRAPH_FINALLY(igraph_destroy, &is);

  /* Rewrite the map as edge weights */
  m=igraph_vector_size(&map2);
  for (i=0; i<m; i++) {
    VECTOR(map2)[i] = VECTOR(*weights_them)[ (int) VECTOR(map2)[i] ];
  }

  igraph_strength(&is, res, igraph_vss_all(), mode, IGRAPH_LOOPS,
		  /*weights=*/ &map2);

  igraph_destroy(&is);
  igraph_vector_destroy(&map2);
  IGRAPH_FINALLY_CLEAN(2);

  return 0;
}

int igraph_local_scan_0_them(const igraph_t *us, const igraph_t *them,
			     igraph_vector_t *res,
			     const igraph_vector_t *weights_us,
			     const igraph_vector_t *weights_them,
			     igraph_neimode_t mode) {

  igraph_t is;

  if (igraph_vcount(us) != igraph_vcount(them)) {
    IGRAPH_ERROR("Number of vertices don't match in scan-0", IGRAPH_EINVAL);
  }
  if (igraph_is_directed(us) != igraph_is_directed(them)) {
    IGRAPH_ERROR("Directedness don't match in scan-0", IGRAPH_EINVAL);
  }

  if (weights_us) {
    return igraph_i_local_scan_0_them_w(us, them, res, weights_us,
					weights_them, mode);
  }

  igraph_intersection(&is, us, them, /*edgemap1=*/ 0, /*edgemap2=*/ 0);
  IGRAPH_FINALLY(igraph_destroy, &is);

  igraph_degree(&is, res, igraph_vss_all(), mode, IGRAPH_LOOPS);

  igraph_destroy(&is);
  IGRAPH_FINALLY_CLEAN(1);

  return 0;
}

int igraph_local_scan_1_them(const igraph_t *us, const igraph_t *them,
			     igraph_vector_t *res, igraph_neimode_t mode) {

  int no_of_nodes=igraph_vcount(us);
  igraph_inclist_t incs_us, incs_them;
  igraph_vector_int_t neis;
  int node;

  if (igraph_vcount(them) != no_of_nodes) {
    IGRAPH_ERROR("Number of vertices must match in scan-1", IGRAPH_EINVAL);
  }
  if (igraph_is_directed(us) != igraph_is_directed(them)) {
    IGRAPH_ERROR("Directedness must match in scan-1", IGRAPH_EINVAL);
  }

  igraph_inclist_init(us, &incs_us, mode);
  IGRAPH_FINALLY(igraph_inclist_destroy, &incs_us);
  igraph_inclist_init(them, &incs_them, mode);
  IGRAPH_FINALLY(igraph_inclist_destroy, &incs_them);

  igraph_vector_int_init(&neis, no_of_nodes);
  IGRAPH_FINALLY(igraph_vector_int_destroy, &neis);

  igraph_vector_resize(res, no_of_nodes);
  igraph_vector_null(res);

  for (node=0; node < no_of_nodes; node++) {
    igraph_vector_int_t *edges1_us=igraph_inclist_get(&incs_us, node);
    igraph_vector_int_t *edges1_them=igraph_inclist_get(&incs_them, node);
    int len1_us=igraph_vector_int_size(edges1_us);
    int len1_them=igraph_vector_int_size(edges1_them);
    int i;

    IGRAPH_ALLOW_INTERRUPTION();

    /* Mark neighbors and self in us */
    VECTOR(neis)[node] = node+1;
    for (i = 0; i < len1_us; i++) {
      int e=VECTOR(*edges1_us)[i];
      int nei=IGRAPH_OTHER(us, e, node);
      VECTOR(neis)[nei] = node+1;
    }

    /* Crawl neighbors in them, first ego */
    for (i = 0; i < len1_them; i++) {
      int e=VECTOR(*edges1_them)[i];
      int nei=IGRAPH_OTHER(them, e, node);
      if (VECTOR(neis)[nei] == node+1) { VECTOR(*res)[node] += 1; }
    }
    /* Then the rest */
    for (i = 0; i < len1_us; i++) {
      int e=VECTOR(*edges1_us)[i];
      int nei=IGRAPH_OTHER(us, e, node);
      igraph_vector_int_t *edges2_them=igraph_inclist_get(&incs_them, nei);
      int j, len2_them=igraph_vector_int_size(edges2_them);
      for (j = 0; j < len2_them; j++) {
	int e2=VECTOR(*edges2_them)[j];
	int nei2=IGRAPH_OTHER(them, e2, nei);
	if (VECTOR(neis)[nei2] == node+1) {
	  VECTOR(*res)[node] += 1;
	}
      }
    }

    /* For undirected, it was double counted */
    if (mode == IGRAPH_ALL || ! igraph_is_directed(us)) {
      VECTOR(*res)[node] /= 2.0;
    }

  } /* node < no_of_nodes */

  igraph_vector_int_destroy(&neis);
  igraph_inclist_destroy(&incs_them);
  igraph_inclist_destroy(&incs_us);
  IGRAPH_FINALLY_CLEAN(3);

  return 0;
}
