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

int igraph_scan0(const igraph_t *graph, igraph_vector_t *res,
		 const igraph_vector_t *weights, igraph_neimode_t mode) {
  if (weights) {
    igraph_strength(graph, res, igraph_vss_all(), mode, /*loops=*/ 1,
		    weights);
  } else {
    igraph_degree(graph, res, igraph_vss_all(), mode, /*loops=*/ 1);
  }
  return 0;
}

int igraph_i_trans4_al_simplify(igraph_adjlist_t *al,
																const igraph_vector_int_t *rank);

int igraph_scan1_ecount(const igraph_t *graph, igraph_vector_t *res) {

#define TRIEDGES
#include "triangles_template.h"
#undef TRIEDGES

	return 0;
}

int igraph_scan1_ecount_approximate(const igraph_t *graph,
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

int igraph_scan1_ecount_approximate_eigen(const igraph_t *graph,
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
