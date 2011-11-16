/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2011  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard st, Cambridge, MA, 02138 USA
   
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

/*
 *  SCGlib : A C library for the spectral coarse graining of matrices
 *	as described in the paper: Shrinking Matrices while preserving their
 *	eigenpairs with Application to the Spectral Coarse Graining of Graphs.
 *	Preprint available at <http://people.epfl.ch/david.morton>
 *  
 *	Copyright (C) 2008 David Morton de Lachapelle <david.morton@a3.epfl.ch>
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
 *  02110-1301 USA
 *
 *  DESCRIPTION
 *	-----------
 *    The grouping function takes as argument 'nev' eigenvectors and 
 *	  and tries to minimize the eigenpair shifts induced by the coarse
 *	  graining (Section 5 of the above reference). The eigenvectors are
 *	  stored in a 'nev'x'n' matrix 'v'.
 *	  The 'algo' parameter can take the following values
 *		1  ->  Optimal method (sec. 5.3.1)
 *		2  ->  Intervals+k-means (sec. 5.3.3)
 *		3  ->  Intervals (sec. 5.3.2)
 *		4  ->  Exact SCG (sec. 5.4.1--last paragraph)
 *	  'nt' is a vector of length 'nev' giving either the size of the 
 *	  partitions (if algo = 1) or the number of intervals to cut the
 *	  eigenvectors if algo = 2 or algo = 3. When algo = 4 this parameter
 *	  is ignored. 'maxiter' fixes the maximum number of iterations of
 *	  the k-means algorithm, and is only considered when algo = 2.
 *	  All the algorithms try to find a minimizing partition of
 *	  ||v_i-Pv_i|| where P is a problem-specific projector and v_i denotes
 *	  the eigenvectors stored in v. The final partition is worked out
 *	  as decribed in Method 1 of Section 5.4.2.
 *	  'matrix' provides the type of SCG (i.e. the form of P). So far,
 *	  the options are those described in section 6, that is:
 *		1  ->  Symmetric (sec. 6.1)
 *		2  ->  Laplacian (sec. 6.2)
 *		3  ->  Stochastic (sec. 6.3)
 *	  In the stochastic case, a valid distribution probability 'p' must be
 *	  provided. In all other cases, 'p' is ignored and can be set to NULL.
 *	  The group labels in the final partition are given in 'gr' as positive
 *	  consecutive integers starting from 0.
 */

#include "igraph_scg.h"
#include "igraph_eigen.h"
#include "igraph_interface.h"
#include "igraph_structural.h"
#include "igraph_constructors.h"
#include "igraph_conversion.h"
#include "igraph_memory.h"

#include "scg_headers.h"

#include "math.h"

int igraph_scg_grouping(const igraph_matrix_t *V, 
			igraph_vector_t *groups,
			igraph_integer_t intervals,
			const igraph_vector_t *intervals_vector,
			igraph_scg_matrix_t matrix_type,
			igraph_scg_algorithm_t algorithm,
			const igraph_vector_t *p,
			igraph_integer_t maxiter) {

  int no_of_nodes=igraph_matrix_nrow(V);
  int nev=igraph_matrix_ncol(V);
  igraph_matrix_int_t gr_mat;
  int i;

  if (intervals_vector && igraph_vector_size(intervals_vector) != nev) {
    IGRAPH_ERROR("Invalid length for interval specification", IGRAPH_EINVAL);
  }

  if (!intervals_vector && algorithm != IGRAPH_SCG_EXACT) {
    if (intervals <= 1 || intervals >= no_of_nodes) {
      IGRAPH_ERROR("Invalid interval specification", IGRAPH_EINVAL);
    }
  } else if (algorithm != IGRAPH_SCG_EXACT) {
    igraph_real_t min, max;
    igraph_vector_minmax(intervals_vector, &min, &max);
    if (min <= 1 || max >= no_of_nodes) {
      IGRAPH_ERROR("Invalid interval specification", IGRAPH_EINVAL);
    }
  }

  if (matrix_type == IGRAPH_SCG_STOCHASTIC && !p) {
    IGRAPH_ERROR("`p' must be given for the stochastic matrix case", 
		 IGRAPH_EINVAL);
  }

  if (p && igraph_vector_size(p) != no_of_nodes) {
    IGRAPH_ERROR("Invalid `p' vector size", IGRAPH_EINVAL);
  }

  IGRAPH_CHECK(igraph_vector_resize(groups, no_of_nodes));

#define INVEC(i) (intervals_vector ? VECTOR(*intervals_vector)[i] : intervals)

  IGRAPH_CHECK(igraph_matrix_int_init(&gr_mat, no_of_nodes, nev));
  IGRAPH_FINALLY(igraph_matrix_int_destroy, &gr_mat);

  switch (algorithm) {
  case IGRAPH_SCG_OPTIMUM:
    for (i=0; i<nev; i++) {
      IGRAPH_CHECK(igraph_i_optimal_partition(&MATRIX(*V, 0, i), 
					      &MATRIX(gr_mat, 0, i), 
					      no_of_nodes, INVEC(i),
					      matrix_type, 
					      p ? VECTOR(*p) : 0));
    }
    break;
  case IGRAPH_SCG_INTERV_KM:
    for (i=0; i<nev; i++) {
      igraph_vector_t tmpv;
      igraph_vector_view(&tmpv, &MATRIX(*V, 0, i), no_of_nodes);
      IGRAPH_CHECK(igraph_i_intervals_plus_kmeans(&tmpv,
						  &MATRIX(gr_mat, 0, i),
						  no_of_nodes, INVEC(i),
						  maxiter));
    }
    break;
  case IGRAPH_SCG_INTERV:
    for (i=0; i<nev; i++) {
      igraph_vector_t tmpv;
      igraph_vector_view(&tmpv, &MATRIX(*V, 0, i), no_of_nodes);
      IGRAPH_CHECK(igraph_i_intervals_method(&tmpv,
					     &MATRIX(gr_mat, 0, i),
					     no_of_nodes, INVEC(i)));
    }
    break;
  case IGRAPH_SCG_EXACT:
    for (i=0; i<nev; i++) {
      IGRAPH_CHECK(igraph_i_exact_coarse_graining(&MATRIX(*V, 0, i),
						  &MATRIX(gr_mat, 0, i),
						  no_of_nodes));
    }
    break;
  }

#undef INVEC

  if (nev==1) {
    for (i=0; i<no_of_nodes; i++) {
      VECTOR(*groups)[i] = MATRIX(gr_mat, i, 0);
    }
  } else {
    igraph_i_scg_groups_t *g = igraph_Calloc(no_of_nodes, 
					     igraph_i_scg_groups_t);
    int gr_nb=0;

    IGRAPH_CHECK(igraph_matrix_int_transpose(&gr_mat));
    for(i=0; i<no_of_nodes; i++){
      g[i].ind = i;
      g[i].n = nev;
      g[i].gr = &MATRIX(gr_mat, 0, i);
    }
		
    qsort(g, no_of_nodes, sizeof(igraph_i_scg_groups_t), 
	  igraph_i_compare_groups);
    VECTOR(*groups)[g[0].ind] = gr_nb;
    for(i=1; i<no_of_nodes; i++){
      if(igraph_i_compare_groups(&g[i], &g[i-1]) != 0) gr_nb++;
      VECTOR(*groups)[g[i].ind] = gr_nb;
    }
    igraph_Free(g);
  }

  igraph_matrix_int_destroy(&gr_mat);
  IGRAPH_FINALLY_CLEAN(1);

  return 0;
}

int igraph_i_scg_semiprojectors_sym(const igraph_vector_t *groups,
				    igraph_matrix_t *L,
				    igraph_matrix_t *R,
				    igraph_sparsemat_t *Lsparse,
				    igraph_sparsemat_t *Rsparse,
				    int no_of_groups,
				    int no_of_nodes) {

  igraph_vector_t tab;
  int i;

  IGRAPH_VECTOR_INIT_FINALLY(&tab, no_of_groups);
  for (i=0; i<no_of_nodes; i++) {
    VECTOR(tab)[ (int) VECTOR(*groups)[i] ] += 1;
  }
  for (i=0; i<no_of_groups; i++) {
    VECTOR(tab)[i] = sqrt(VECTOR(tab)[i]);
  }
  
  if (L) { 
    IGRAPH_CHECK(igraph_matrix_resize(L, no_of_groups, no_of_nodes));
    igraph_matrix_null(L);
    for (i=0; i<no_of_nodes; i++) {
      int g=VECTOR(*groups)[i];
      MATRIX(*L, g, i) = 1/VECTOR(tab)[g];
    }
  }
  
  if (R) {
    if (L) { 
      IGRAPH_CHECK(igraph_matrix_update(R, L));
    } else {
      IGRAPH_CHECK(igraph_matrix_resize(R, no_of_groups, no_of_nodes));
      igraph_matrix_null(R);
      for (i=0; i<no_of_nodes; i++) {
	int g=VECTOR(*groups)[i];
	MATRIX(*R, g, i) = 1/VECTOR(tab)[g];
      }
    }
  }
  
  if (Lsparse) {
    IGRAPH_CHECK(igraph_sparsemat_init(Lsparse, no_of_groups, no_of_nodes, 
					 /* nzmax= */ no_of_nodes));
    for (i=0; i<no_of_nodes; i++) {
      int g=VECTOR(*groups)[i];
      IGRAPH_CHECK(igraph_sparsemat_entry(Lsparse, g, i, 1/VECTOR(tab)[g]));
    }
  }
  
  if (Rsparse) {
    IGRAPH_CHECK(igraph_sparsemat_init(Rsparse, no_of_groups, no_of_nodes,
					 /* nzmax= */ no_of_nodes));
    for (i=0; i<no_of_nodes; i++) {
      int g=VECTOR(*groups)[i];
      IGRAPH_CHECK(igraph_sparsemat_entry(Rsparse, g, i, 1/VECTOR(tab)[g]));
    }    
  }

  igraph_vector_destroy(&tab);
  IGRAPH_FINALLY_CLEAN(1);
  
  return 0;
}

int igraph_i_scg_semiprojectors_lap(const igraph_vector_t *groups,
				    igraph_matrix_t *L,
				    igraph_matrix_t *R,
				    igraph_sparsemat_t *Lsparse,
				    igraph_sparsemat_t *Rsparse,
				    int no_of_groups,
				    int no_of_nodes,
				    igraph_scg_norm_t norm) {
  
  igraph_vector_t tab;
  int i;

  IGRAPH_VECTOR_INIT_FINALLY(&tab, no_of_groups);
  for (i=0; i<no_of_nodes; i++) {
    VECTOR(tab)[ (int) VECTOR(*groups)[i] ] += 1;
  }
  for (i=0; i<no_of_groups; i++) {
    VECTOR(tab)[i] = VECTOR(tab)[i];
  }

  if (norm == IGRAPH_SCG_NORM_ROW) {
    if (L) {
      IGRAPH_CHECK(igraph_matrix_resize(L, no_of_groups, no_of_nodes));
      igraph_matrix_null(L);
      for (i=0; i<no_of_nodes; i++) {
	int g=VECTOR(*groups)[i];
	MATRIX(*L, g, i) = 1.0 / VECTOR(tab)[g];
      }
    }
    if (R) {
      IGRAPH_CHECK(igraph_matrix_resize(R, no_of_groups, no_of_nodes));
      igraph_matrix_null(R);
      for (i=0; i<no_of_nodes; i++) {
	int g=VECTOR(*groups)[i];
	MATRIX(*R, g, i) = 1.0;
      }
    }
    if (Lsparse) {
      IGRAPH_CHECK(igraph_sparsemat_init(Lsparse, no_of_groups, no_of_nodes,
					 /* nzmax= */ no_of_nodes));
      for (i=0; i<no_of_nodes; i++) {
	int g=VECTOR(*groups)[i];
	IGRAPH_CHECK(igraph_sparsemat_entry(Lsparse, g, i, 
					    1.0 / VECTOR(tab)[g]));
      }
    }
    if (Rsparse) {
      IGRAPH_CHECK(igraph_sparsemat_init(Rsparse, no_of_groups, no_of_nodes,
					   /* nzmax= */ no_of_nodes));
      for (i=0; i<no_of_nodes; i++) {
	int g=VECTOR(*groups)[i];
	IGRAPH_CHECK(igraph_sparsemat_entry(Rsparse, g, i, 1.0));
      }
    }
  } else {
    if (L) {
      IGRAPH_CHECK(igraph_matrix_resize(L, no_of_groups, no_of_nodes));
      igraph_matrix_null(L);
      for (i=0; i<no_of_nodes; i++) {
	int g=VECTOR(*groups)[i];
	MATRIX(*L, g, i) = 1.0;
      }
    }
    if (R) {
      IGRAPH_CHECK(igraph_matrix_resize(R, no_of_groups, no_of_nodes));
      igraph_matrix_null(R);
      for (i=0; i<no_of_nodes; i++) {
	int g=VECTOR(*groups)[i];
	MATRIX(*R, g, i) = 1.0 / VECTOR(tab)[g];
      }
    }
    if (Lsparse) {
      IGRAPH_CHECK(igraph_sparsemat_init(Lsparse, no_of_groups, no_of_nodes,
					 /* nzmax= */ no_of_nodes));
      for (i=0; i<no_of_nodes; i++) {
	int g=VECTOR(*groups)[i];
	IGRAPH_CHECK(igraph_sparsemat_entry(Lsparse, g, i, 1.0));
      }
    }
    if (Rsparse) {
      IGRAPH_CHECK(igraph_sparsemat_init(Rsparse, no_of_groups, no_of_nodes,
					 /* nzmax= */ no_of_nodes));
      for (i=0; i<no_of_nodes; i++) {
	int g=VECTOR(*groups)[i];
	IGRAPH_CHECK(igraph_sparsemat_entry(Rsparse, g, i, 
					    1.0 / VECTOR(tab)[g]));
      }
    }
    
  }

  igraph_vector_destroy(&tab);
  IGRAPH_FINALLY_CLEAN(1);

  return 0;
}

int igraph_i_scg_semiprojectors_sto(const igraph_vector_t *groups,
				    igraph_matrix_t *L,
				    igraph_matrix_t *R,
				    igraph_sparsemat_t *Lsparse,
				    igraph_sparsemat_t *Rsparse,
				    int no_of_groups,
				    int no_of_nodes,
				    const igraph_vector_t *p,
				    igraph_scg_norm_t norm) {
  
  igraph_vector_t pgr, pnormed;
  int i;

  IGRAPH_VECTOR_INIT_FINALLY(&pgr, no_of_groups);
  IGRAPH_VECTOR_INIT_FINALLY(&pnormed, no_of_nodes);
  for (i=0; i<no_of_nodes; i++) {
    int g=VECTOR(*groups)[i];
    VECTOR(pgr)[g] += VECTOR(*p)[i];
  }
  for (i=0; i<no_of_nodes; i++) {
    int g=VECTOR(*groups)[i];
    VECTOR(pnormed)[i] = VECTOR(*p)[i] / VECTOR(pgr)[g];
  }
  
  if (norm == IGRAPH_SCG_NORM_ROW) {
    if (L) {
      IGRAPH_CHECK(igraph_matrix_resize(L, no_of_groups, no_of_nodes));
      igraph_matrix_null(L);
      for (i=0; i<no_of_nodes; i++) {
	int g=VECTOR(*groups)[i];
	MATRIX(*L, g, i) = VECTOR(pnormed)[i];
      }
    }
    if (R) {
      IGRAPH_CHECK(igraph_matrix_resize(R, no_of_groups, no_of_nodes));
      igraph_matrix_null(R);
      for (i=0; i<no_of_nodes; i++) {
	int g=VECTOR(*groups)[i];
	MATRIX(*R, g, i) = 1.0;
      }      
    }
    if (Lsparse) {
      IGRAPH_CHECK(igraph_sparsemat_init(Lsparse, no_of_groups, no_of_nodes,
					 /* nzmax= */ no_of_nodes));
      for (i=0; i<no_of_nodes; i++) {
	int g=VECTOR(*groups)[i];
	IGRAPH_CHECK(igraph_sparsemat_entry(Lsparse, g, i,
					    VECTOR(pnormed)[i]));
      }
    }
    if (Rsparse) {
      IGRAPH_CHECK(igraph_sparsemat_init(Rsparse, no_of_groups, no_of_nodes,
					 /* nzmax= */ no_of_nodes));
      for (i=0; i<no_of_nodes; i++) {
	int g=VECTOR(*groups)[i];
	IGRAPH_CHECK(igraph_sparsemat_entry(Rsparse, g, i, 1.0));
      }
    }
  } else {
    if (L) {
      IGRAPH_CHECK(igraph_matrix_resize(L, no_of_groups, no_of_nodes));
      igraph_matrix_null(L);
      for (i=0; i<no_of_nodes; i++) {
	int g=VECTOR(*groups)[i];
	MATRIX(*L, g, i) = 1.0;
      }
    }
    if (R) {
      IGRAPH_CHECK(igraph_matrix_resize(R, no_of_groups, no_of_nodes));
      igraph_matrix_null(R);
      for (i=0; i<no_of_nodes; i++) {
	int g=VECTOR(*groups)[i];
	MATRIX(*R, g, i) = VECTOR(pnormed)[i];
      }      
    }
    if (Lsparse) {
      IGRAPH_CHECK(igraph_sparsemat_init(Lsparse, no_of_groups, no_of_nodes,
					   /* nzmax= */ no_of_nodes));
      for (i=0; i<no_of_nodes; i++) {
	int g=VECTOR(*groups)[i];
	IGRAPH_CHECK(igraph_sparsemat_entry(Lsparse, g, i, 1.0));
      }
    }
    if (Rsparse) {
      IGRAPH_CHECK(igraph_sparsemat_init(Rsparse, no_of_groups, no_of_nodes,
					   /* nzmax= */ no_of_nodes));
      for (i=0; i<no_of_nodes; i++) {
	int g=VECTOR(*groups)[i];
	IGRAPH_CHECK(igraph_sparsemat_entry(Rsparse, g, i,
					    VECTOR(pnormed)[i]));
      }
    }
  }
  
  
  igraph_vector_destroy(&pnormed);
  igraph_vector_destroy(&pgr);
  IGRAPH_FINALLY_CLEAN(2);

  return 0;
}

int igraph_scg_semiprojectors(const igraph_vector_t *groups,
			      igraph_scg_matrix_t matrix_type,
			      igraph_matrix_t *L,
			      igraph_matrix_t *R,
			      igraph_sparsemat_t *Lsparse,
			      igraph_sparsemat_t *Rsparse, 
			      const igraph_vector_t *p,
			      igraph_scg_norm_t norm) {
  
  int no_of_nodes=igraph_vector_size(groups);
  int no_of_groups;
  igraph_real_t min, max;

  igraph_vector_minmax(groups, &min, &max);
  no_of_groups=max+1;

  if (min < 0 || max >= no_of_nodes) {
    IGRAPH_ERROR("Invalid membership vector", IGRAPH_EINVAL);
  }

  if (matrix_type == IGRAPH_SCG_STOCHASTIC && !p) {
    IGRAPH_ERROR("`p' must be given for the stochastic matrix case", 
		 IGRAPH_EINVAL);
  }

  if (p && igraph_vector_size(p) != no_of_nodes) {
    IGRAPH_ERROR("Invalid `p' vector length, should match number of vertices",
		 IGRAPH_EINVAL);
  }

  switch (matrix_type) {
  case IGRAPH_SCG_SYMMETRIC:
    IGRAPH_CHECK(igraph_i_scg_semiprojectors_sym(groups, L, R, Lsparse, 
						 Rsparse, no_of_groups,
						 no_of_nodes));
    break;

  case IGRAPH_SCG_LAPLACIAN:
    IGRAPH_CHECK(igraph_i_scg_semiprojectors_lap(groups, L, R, Lsparse, 
						 Rsparse, no_of_groups,
						 no_of_nodes, norm));
    break;

  case IGRAPH_SCG_STOCHASTIC:
    IGRAPH_CHECK(igraph_i_scg_semiprojectors_sto(groups, L, R, Lsparse,
						 Rsparse, no_of_groups,
						 no_of_nodes, p, norm));
    break;    
  }
  
  return 0;
}

int igraph_scg_norm_eps(const igraph_matrix_t *V,
			const igraph_vector_t *groups,
			igraph_vector_t *eps,
			igraph_scg_matrix_t matrix_type,
			const igraph_vector_t *p,
			igraph_scg_norm_t norm) {

  int no_of_nodes=igraph_vector_size(groups);
  int no_of_groups;
  int no_of_vectors=igraph_matrix_ncol(V);
  igraph_real_t min, max;
  igraph_sparsemat_t Lsparse, Rsparse, Lsparse2, Rsparse2, Rsparse3, proj;
  igraph_vector_t x, res;
  int k, i;
  
  if (igraph_matrix_nrow(V) != no_of_nodes) {
    IGRAPH_ERROR("Eigenvector length and group vector length do not match",
		 IGRAPH_EINVAL);
  }

  igraph_vector_minmax(groups, &min, &max);
  no_of_groups=max+1;
  
  if (min < 0 || max >= no_of_nodes) {
    IGRAPH_ERROR("Invalid membership vector", IGRAPH_EINVAL);
  }

  if (matrix_type == IGRAPH_SCG_STOCHASTIC && !p) {
    IGRAPH_ERROR("`p' must be given for the stochastic matrix case", 
		 IGRAPH_EINVAL);
  }

  if (p && igraph_vector_size(p) != no_of_nodes) {
    IGRAPH_ERROR("Invalid `p' vector length, should match number of vertices",
		 IGRAPH_EINVAL);
  }

  IGRAPH_CHECK(igraph_sparsemat_init(&Lsparse, no_of_groups, no_of_nodes,
				     no_of_nodes));
  IGRAPH_FINALLY(igraph_sparsemat_destroy, &Lsparse);
  IGRAPH_CHECK(igraph_sparsemat_init(&Rsparse, no_of_groups, no_of_nodes,
				     no_of_nodes));
  IGRAPH_FINALLY(igraph_sparsemat_destroy, &Rsparse);
  
  IGRAPH_CHECK(igraph_scg_semiprojectors(groups, matrix_type, /* L= */ 0, 
					 /* R= */ 0, &Lsparse, &Rsparse, p,
					 norm));

  IGRAPH_CHECK(igraph_sparsemat_compress(&Lsparse, &Lsparse2));
  IGRAPH_FINALLY(igraph_sparsemat_destroy, &Lsparse2);
  IGRAPH_CHECK(igraph_sparsemat_compress(&Rsparse, &Rsparse2));
  IGRAPH_FINALLY(igraph_sparsemat_destroy, &Rsparse2);
  IGRAPH_CHECK(igraph_sparsemat_transpose(&Rsparse2, &Rsparse3, 
					  /*values=*/ 1));
  IGRAPH_FINALLY(igraph_sparsemat_destroy, &Rsparse3);
  
  IGRAPH_CHECK(igraph_sparsemat_multiply(&Rsparse3, &Lsparse2, &proj));
  IGRAPH_FINALLY(igraph_sparsemat_destroy, &proj);

  IGRAPH_VECTOR_INIT_FINALLY(&res, no_of_nodes);
  IGRAPH_CHECK(igraph_vector_resize(eps, no_of_vectors));
  
  for (k = 0; k < no_of_vectors; k++) {
    igraph_vector_view(&x, &MATRIX(*V, 0, k), no_of_nodes);
    igraph_vector_null(&res);
    IGRAPH_CHECK(igraph_sparsemat_gaxpy(&proj, &x, &res));
    VECTOR(*eps)[k] = 0.0;
    for (i = 0; i < no_of_nodes; i++) {
      igraph_real_t di=MATRIX(*V, i, k) - VECTOR(res)[i];
      VECTOR(*eps)[k] += di * di;
    }
    VECTOR(*eps)[k] = sqrt(VECTOR(*eps)[k]);
  }

  igraph_vector_destroy(&res);
  igraph_sparsemat_destroy(&proj);
  igraph_sparsemat_destroy(&Rsparse3);
  igraph_sparsemat_destroy(&Rsparse2);
  igraph_sparsemat_destroy(&Lsparse2);
  igraph_sparsemat_destroy(&Rsparse);
  igraph_sparsemat_destroy(&Lsparse);
  IGRAPH_FINALLY_CLEAN(7);

  return 0;
}

int igraph_get_stochastic(const igraph_t *graph, 
			  igraph_matrix_t *matrix,
			  igraph_scg_norm_t norm) {
  
  int no_of_nodes=igraph_vcount(graph);
  igraph_real_t sum;
  int i, j;
  
  IGRAPH_CHECK(igraph_get_adjacency(graph, matrix, 
				    IGRAPH_GET_ADJACENCY_BOTH, /*eids=*/ 0));
  
  switch (norm) {
  case IGRAPH_SCG_NORM_ROW:
    for (i=0; i<no_of_nodes; i++) {
      sum=0.0; 
      for (j=0; j<no_of_nodes; j++) { 
	sum += MATRIX(*matrix, i, j);
      }
      for (j=0; j<no_of_nodes; j++) {
	MATRIX(*matrix, i, j) /= sum;
      }
    }
    break;
  case IGRAPH_SCG_NORM_COL:
    for (i=0; i<no_of_nodes; i++) {
      sum=0.0; 
      for (j=0; j<no_of_nodes; j++) { 
	sum += MATRIX(*matrix, j, i);
      }
      for (j=0; j<no_of_nodes; j++) {
	MATRIX(*matrix, j, i) /= sum;
      }
    }
    break;
  }

  return 0;
}

int igraph_i_normalize_sparsemat(igraph_sparsemat_t *sparsemat, 
				 igraph_scg_norm_t norm) {
  igraph_vector_t sum;
  int no_of_nodes=igraph_sparsemat_nrow(sparsemat);
  int i;
  
  IGRAPH_VECTOR_INIT_FINALLY(&sum, no_of_nodes);

  switch (norm) {       
  case IGRAPH_SCG_NORM_ROW:
    IGRAPH_CHECK(igraph_sparsemat_rowsums(sparsemat, &sum));
    for (i=0; i<no_of_nodes; i++) {
      if (VECTOR(sum)[i] == 0.0) {
	IGRAPH_ERROR("Zero out-degree vertices not allowed", 
		     IGRAPH_EINVAL);
      }
      VECTOR(sum)[i] = 1.0 / VECTOR(sum)[i];
    }
    IGRAPH_CHECK(igraph_sparsemat_scale_rows(sparsemat, &sum));
    break;
  case IGRAPH_SCG_NORM_COL:
    IGRAPH_CHECK(igraph_sparsemat_colsums(sparsemat, &sum));
    for (i=0; i<no_of_nodes; i++) {
      if (VECTOR(sum)[i] == 0.0) {
	IGRAPH_ERROR("Zero out-degree vertices not allowed", 
		     IGRAPH_EINVAL);
      }
      VECTOR(sum)[i] = 1.0 / VECTOR(sum)[i];
    }
    IGRAPH_CHECK(igraph_sparsemat_scale_cols(sparsemat, &sum));
    break;
  }

  igraph_vector_destroy(&sum);
  IGRAPH_FINALLY_CLEAN(1);
  
  return 0;
}

int igraph_get_stochastic_sparsemat(const igraph_t *graph, 
				    igraph_sparsemat_t *sparsemat,
				    igraph_scg_norm_t norm) {

  IGRAPH_CHECK(igraph_get_sparsemat(graph, sparsemat));
  IGRAPH_FINALLY(igraph_sparsemat_destroy, sparsemat);
  IGRAPH_CHECK(igraph_i_normalize_sparsemat(sparsemat, norm));
  IGRAPH_FINALLY_CLEAN(1);
  
  return 0;
}

int igraph_i_matrix_laplacian(const igraph_matrix_t *matrix, 
			      igraph_matrix_t *mymatrix, 
			      igraph_scg_norm_t norm) {

  igraph_vector_t degree;
  int i, j, n=igraph_matrix_nrow(matrix);
  IGRAPH_CHECK(igraph_matrix_resize(mymatrix, n, n));

  IGRAPH_VECTOR_INIT_FINALLY(&degree, n);  

  if (norm==IGRAPH_SCG_NORM_ROW) { 
    IGRAPH_CHECK(igraph_matrix_rowsum(matrix, &degree));
  } else {
    IGRAPH_CHECK(igraph_matrix_colsum(matrix, &degree));
  }
  for (i=0; i<n; i++) {
    VECTOR(degree)[i] -= MATRIX(*matrix, i, i);
  }

  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) {
      MATRIX(*mymatrix, i, j) = - MATRIX(*matrix, i, j);
    }
    MATRIX(*mymatrix, i, i) = VECTOR(degree)[i];
  }
  
  igraph_vector_destroy(&degree);
  IGRAPH_FINALLY_CLEAN(1);

  return 0;
}

int igraph_i_matrix_stochastic(const igraph_matrix_t *matrix, 
			       igraph_matrix_t *mymatrix, 
			       igraph_scg_norm_t norm) {

  int i, j, n=igraph_matrix_nrow(matrix);
  IGRAPH_CHECK(igraph_matrix_copy(mymatrix, matrix));

  if (norm==IGRAPH_SCG_NORM_ROW){
    for (i=0; i<n; i++) {
      igraph_real_t sum=0.0;
      for (j=0; j<n; j++) {
	sum += MATRIX(*matrix, i, j);
      }
      if (sum == 0) { IGRAPH_WARNING("Zero degree vertices"); }
      for (j=0; j<n; j++) {
	MATRIX(*mymatrix, i, j) = MATRIX(*matrix, i, j) / sum;
      }
    }
  } else {
    for (i=0; i<n; i++) {
      igraph_real_t sum=0.0;
      for (j=0; j<n; j++) {
	sum += MATRIX(*matrix, j, i);
      }
      if (sum == 0) { IGRAPH_WARNING("Zero degree vertices"); }
      for (j=0; j<n; j++) {
	MATRIX(*mymatrix, j, i) = MATRIX(*matrix, j, i) / sum;
      }
    }    
  }	

  return 0;
}

int igraph_i_sparsemat_laplacian(const igraph_sparsemat_t *sparse,
				 igraph_sparsemat_t *mysparse,
				 igraph_scg_norm_t norm) {

  igraph_vector_t degree;
  int i, n=igraph_sparsemat_nrow(sparse);
  int nzmax=igraph_sparsemat_nzmax(sparse);
  igraph_sparsemat_iterator_t it;

  IGRAPH_CHECK(igraph_sparsemat_init(mysparse, n, n, nzmax+n));
  IGRAPH_FINALLY(igraph_sparsemat_destroy, mysparse);
  igraph_sparsemat_iterator_init(&it, (igraph_sparsemat_t *) sparse);

  IGRAPH_VECTOR_INIT_FINALLY(&degree, n);  
  for (igraph_sparsemat_iterator_reset(&it); 
       !igraph_sparsemat_iterator_end(&it);
       igraph_sparsemat_iterator_next(&it)) {
    int row=igraph_sparsemat_iterator_row(&it);
    int col=igraph_sparsemat_iterator_col(&it);
    if (row != col) { 
      igraph_real_t val=igraph_sparsemat_iterator_get(&it);
      if (norm == IGRAPH_SCG_NORM_ROW) { 
	VECTOR(degree)[row] += val;
      } else {
	VECTOR(degree)[col] += val;
      }
    }
  }

  /* Diagonal */
  for (i=0; i<n; i++) {
    igraph_real_t t=VECTOR(degree)[i] > 0 ? 1 : 0;
    igraph_sparsemat_entry(mysparse, i, i, t);
  }

  /* And the rest, filter out diagonal elements */
  for (igraph_sparsemat_iterator_reset(&it); 
       !igraph_sparsemat_iterator_end(&it);
       igraph_sparsemat_iterator_next(&it)) {
    int row=igraph_sparsemat_iterator_row(&it);
    int col=igraph_sparsemat_iterator_col(&it);
    if (row != col) { 
      igraph_real_t val=igraph_sparsemat_iterator_get(&it);
      igraph_sparsemat_entry(mysparse, row, col, -val);
    }
  }

  igraph_vector_destroy(&degree);
  IGRAPH_FINALLY_CLEAN(2); 	/* + mysparse */
  
  return 0;
}

int igraph_i_sparsemat_stochastic(const igraph_sparsemat_t *sparse,
				  igraph_sparsemat_t *mysparse,
				  igraph_scg_norm_t norm) {
  
  IGRAPH_CHECK(igraph_sparsemat_copy(mysparse, sparse));
  IGRAPH_FINALLY(igraph_sparsemat_destroy, mysparse);
  IGRAPH_CHECK(igraph_i_normalize_sparsemat(mysparse, norm));
  IGRAPH_FINALLY_CLEAN(1);
  
  return 0;
}

int igraph_i_scg_get_result(igraph_scg_matrix_t type,
			    const igraph_matrix_t *matrix,
			    const igraph_sparsemat_t *sparsemat,
			    const igraph_sparsemat_t *Lsparse,
			    const igraph_sparsemat_t *Rsparse_t,
			    igraph_t *scg_graph,
			    igraph_matrix_t *scg_matrix,
			    igraph_sparsemat_t *scg_sparsemat, 
			    igraph_bool_t directed) {

  /* We need to calculate either scg_matrix (if input is dense), or
     scg_sparsemat (if input is sparse). For the latter we might need
     to temporarily use another matrix. */
  
  
  if (matrix) {
    igraph_matrix_t *my_scg_matrix=scg_matrix, v_scg_matrix;
    igraph_matrix_t tmp;
    igraph_sparsemat_t *myLsparse=(igraph_sparsemat_t *) Lsparse, v_Lsparse;

    if (!scg_matrix) {
      my_scg_matrix=&v_scg_matrix;
      IGRAPH_CHECK(igraph_matrix_init(my_scg_matrix, 0, 0));
      IGRAPH_FINALLY(igraph_matrix_destroy, my_scg_matrix);
    }

    if (!igraph_sparsemat_is_cc(Lsparse)) {
      myLsparse=&v_Lsparse;
      IGRAPH_CHECK(igraph_sparsemat_compress(Lsparse, myLsparse));
      IGRAPH_FINALLY(igraph_sparsemat_destroy, myLsparse);
    }

    IGRAPH_CHECK(igraph_matrix_init(&tmp, 0, 0));
    IGRAPH_FINALLY(igraph_matrix_destroy, &tmp);
    IGRAPH_CHECK(igraph_sparsemat_dense_multiply(matrix, Rsparse_t, &tmp));
    IGRAPH_CHECK(igraph_sparsemat_multiply_by_dense(myLsparse, &tmp, 
						    my_scg_matrix));
    igraph_matrix_destroy(&tmp);
    IGRAPH_FINALLY_CLEAN(1);

    if (scg_sparsemat) {
      IGRAPH_CHECK(igraph_matrix_as_sparsemat(scg_sparsemat, my_scg_matrix,
					      /* tol= */ 0));
      IGRAPH_FINALLY(igraph_sparsemat_destroy, scg_sparsemat);
    }

    if (scg_graph) {
      if (type != IGRAPH_SCG_LAPLACIAN) {
	IGRAPH_CHECK(igraph_weighted_adjacency(scg_graph, my_scg_matrix, 
					       directed ? 
					       IGRAPH_ADJ_DIRECTED : 
					       IGRAPH_ADJ_UNDIRECTED, 
					       "weight", /*loops=*/ 1));
      } else {
	int i, j, n=igraph_matrix_nrow(my_scg_matrix);
	igraph_matrix_t tmp;
	IGRAPH_MATRIX_INIT_FINALLY(&tmp, n, n);
	for (i=0; i<n; i++) {
	  for (j=0; j<n; j++) {
	    MATRIX(tmp, i, j) = -MATRIX(*my_scg_matrix, i, j);
	  }
	  MATRIX(tmp, i, i) = 0;
	}
	IGRAPH_CHECK(igraph_weighted_adjacency(scg_graph, &tmp, directed ?
					       IGRAPH_ADJ_DIRECTED :
					       IGRAPH_ADJ_UNDIRECTED, 
					       "weight", /*loops=*/ 0));
	igraph_matrix_destroy(&tmp);
	IGRAPH_FINALLY_CLEAN(1);
      }
      IGRAPH_FINALLY(igraph_destroy, scg_graph);
    }

    if (scg_graph)     { IGRAPH_FINALLY_CLEAN(1); }
    if (scg_sparsemat) { IGRAPH_FINALLY_CLEAN(1); }

    if (!igraph_sparsemat_is_cc(Lsparse)) {
      igraph_sparsemat_destroy(myLsparse);
      IGRAPH_FINALLY_CLEAN(1);
    }

    if (!scg_matrix) {
      igraph_matrix_destroy(my_scg_matrix);
      IGRAPH_FINALLY_CLEAN(1);
    }

  } else { /* sparsemat */
    igraph_sparsemat_t *my_scg_sparsemat=scg_sparsemat, v_scg_sparsemat;
    igraph_sparsemat_t tmp, *mysparsemat=(igraph_sparsemat_t *) sparsemat,
      v_sparsemat, *myLsparse=(igraph_sparsemat_t *) Lsparse, v_Lsparse;
    if (!scg_sparsemat) {
      my_scg_sparsemat=&v_scg_sparsemat;
    }
    if (!igraph_sparsemat_is_cc(sparsemat)) {
      mysparsemat=&v_sparsemat;
      IGRAPH_CHECK(igraph_sparsemat_compress(sparsemat, mysparsemat));
      IGRAPH_FINALLY(igraph_sparsemat_destroy, mysparsemat);
    }
    if (!igraph_sparsemat_is_cc(Lsparse)) {
      myLsparse=&v_Lsparse;
      IGRAPH_CHECK(igraph_sparsemat_compress(Lsparse, myLsparse));
      IGRAPH_FINALLY(igraph_sparsemat_destroy, myLsparse);
    }
    IGRAPH_CHECK(igraph_sparsemat_multiply(mysparsemat, Rsparse_t, 
					   &tmp));
    IGRAPH_FINALLY(igraph_sparsemat_destroy, &tmp);
    IGRAPH_CHECK(igraph_sparsemat_multiply(myLsparse, &tmp,
					   my_scg_sparsemat));
    igraph_sparsemat_destroy(&tmp);
    IGRAPH_FINALLY_CLEAN(1);
    IGRAPH_FINALLY(igraph_sparsemat_destroy, my_scg_sparsemat);

    if (scg_matrix) {
      IGRAPH_CHECK(igraph_sparsemat_as_matrix(scg_matrix, my_scg_sparsemat));
    }
    if (scg_graph) {
      if (type != IGRAPH_SCG_LAPLACIAN) {
	IGRAPH_CHECK(igraph_weighted_sparsemat(scg_graph, my_scg_sparsemat,
					       directed, "weight",
					       /*loops=*/ 1));
      } else {
	igraph_sparsemat_t tmp;
	IGRAPH_CHECK(igraph_sparsemat_copy(&tmp, my_scg_sparsemat));
	IGRAPH_FINALLY(igraph_sparsemat_destroy, &tmp);
	IGRAPH_CHECK(igraph_sparsemat_neg(&tmp));
	IGRAPH_CHECK(igraph_weighted_sparsemat(scg_graph, &tmp, directed, 
					       "weight", /*loops=*/ 0));
	igraph_sparsemat_destroy(&tmp);
	IGRAPH_FINALLY_CLEAN(1);
      }
      IGRAPH_FINALLY(igraph_destroy, scg_graph);
    }

    if (scg_graph) { IGRAPH_FINALLY_CLEAN(1); }
    if (!scg_sparsemat) {
      igraph_sparsemat_destroy(my_scg_sparsemat);
    }
    IGRAPH_FINALLY_CLEAN(1);	/* my_scg_sparsemat */
    if (!igraph_sparsemat_is_cc(Lsparse)) {
      igraph_sparsemat_destroy(myLsparse);
      IGRAPH_FINALLY_CLEAN(1);
    }
    if (!igraph_sparsemat_is_cc(sparsemat)) {
      igraph_sparsemat_destroy(mysparsemat);
      IGRAPH_FINALLY_CLEAN(1);
    }
  }

  return 0;
}

int igraph_i_scg_common_checks(const igraph_t *graph, 
			       const igraph_matrix_t *matrix,
			       const igraph_sparsemat_t *sparsemat,
			       const igraph_vector_t *ev,
			       igraph_integer_t intervals,
			       const igraph_vector_t *intervals_vector,
			       const igraph_matrix_t *evec,
			       const igraph_matrix_complex_t *evec_cmplx,
			       const igraph_vector_t *groups,
			       const igraph_t *scg_graph,
			       const igraph_matrix_t *scg_matrix,
			       const igraph_sparsemat_t *scg_sparsemat,
			       const igraph_vector_t *p,
			       igraph_real_t *evmin, igraph_real_t *evmax) {

  int no_of_nodes=-1;
  igraph_real_t min, max;
  int no_of_ev=igraph_vector_size(ev);

  if ( (graph?1:0) + (matrix?1:0) + (sparsemat?1:0) != 1 ) {
    IGRAPH_ERROR("Give exactly one of `graph', `matrix' and `sparsemat'",
		 IGRAPH_EINVAL);
  }
  
  if (graph) { 
    no_of_nodes = igraph_vcount(graph);
  } else if (matrix) {
    no_of_nodes = igraph_matrix_nrow(matrix);
  } else if (sparsemat) {
    no_of_nodes = igraph_sparsemat_nrow(sparsemat);
  }
  
  if ((matrix && igraph_matrix_ncol(matrix) != no_of_nodes) ||
      (sparsemat && igraph_sparsemat_ncol(sparsemat) != no_of_nodes)) {
    IGRAPH_ERROR("Matrix must be square", IGRAPH_NONSQUARE);
  }

  igraph_vector_minmax(ev, evmin, evmax);
  if (*evmin < 0 || *evmax >= no_of_nodes) {
    IGRAPH_ERROR("Invalid eigenvectors given", IGRAPH_EINVAL);
  }

  if (!intervals_vector && (intervals <= 1 || intervals >= no_of_nodes)) {
    IGRAPH_ERROR("Invalid interval specification", IGRAPH_EINVAL);
  }
  
  if (intervals_vector) { 
    if (igraph_vector_size(intervals_vector) != no_of_ev) {
      IGRAPH_ERROR("Invalid length for interval specification", 
		   IGRAPH_EINVAL);
    }
    igraph_vector_minmax(intervals_vector, &min, &max);
    if (min <= 1 || max >= no_of_nodes) { 
      IGRAPH_ERROR("Invalid interval specification", IGRAPH_EINVAL);
    }
  }

  if (evec && igraph_matrix_size(evec) != 0 && 
      (igraph_matrix_ncol(evec) != no_of_ev ||
       igraph_matrix_nrow(evec) != no_of_nodes)) {
    IGRAPH_ERROR("Invalid eigenvector matrix size", IGRAPH_EINVAL);
  }

  if (evec_cmplx && igraph_matrix_complex_size(evec_cmplx) != 0 &&
      (igraph_matrix_complex_ncol(evec_cmplx) != no_of_ev || 
       igraph_matrix_complex_nrow(evec_cmplx) != no_of_nodes)) {
    IGRAPH_ERROR("Invalid eigenvector matrix size", IGRAPH_EINVAL);
  }

  if (groups && igraph_vector_size(groups) != 0 && 
      igraph_vector_size(groups) != no_of_nodes) {
    IGRAPH_ERROR("Invalid `groups' vector size", IGRAPH_EINVAL);
  }

  if ( (scg_graph!=0) + (scg_matrix!=0) + (scg_sparsemat!=0) == 0 ) {
    IGRAPH_ERROR("No output is requested, please give at least one of "
		 "`scg_graph', `scg_matrix' and `scg_sparsemat'", 
		 IGRAPH_EINVAL);
  }

  if (p && igraph_vector_size(p) != 0 &&
      igraph_vector_size(p) != no_of_nodes) {
    IGRAPH_ERROR("Invalid `p' vector size", IGRAPH_EINVAL);
  }

  return 0;
}
			       

int igraph_scg_adjacency(const igraph_t *graph,
			 const igraph_matrix_t *matrix,
			 const igraph_sparsemat_t *sparsemat,
			 const igraph_vector_t *ev,
			 igraph_integer_t intervals,
			 const igraph_vector_t *intervals_vector,
			 igraph_scg_algorithm_t algorithm,
			 igraph_vector_t *eval,
			 igraph_matrix_t *evec,
			 igraph_vector_t *groups,
			 igraph_bool_t use_arpack,
			 igraph_integer_t maxiter,
			 igraph_t *scg_graph,
			 igraph_matrix_t *scg_matrix,
			 igraph_sparsemat_t *scg_sparsemat,
			 igraph_matrix_t *L,
			 igraph_matrix_t *R,
			 igraph_sparsemat_t *Lsparse,
			 igraph_sparsemat_t *Rsparse) {

  igraph_sparsemat_t *mysparsemat=(igraph_sparsemat_t*) sparsemat, 
    real_sparsemat;  
  int no_of_ev=igraph_vector_size(ev);
  /* eigenvectors are calculated and returned */
  igraph_bool_t do_evec= evec && igraph_matrix_size(evec)==0;
  /* groups are calculated */
  igraph_bool_t do_groups= !groups || igraph_vector_size(groups)==0;
  /* eigenvectors are not returned but must be calculated for groups */
  igraph_bool_t tmp_evec= !do_evec && do_groups;
  /* need temporary vector for groups */
  igraph_bool_t tmp_groups= !groups;
  igraph_matrix_t myevec;
  igraph_vector_t mygroups;
  igraph_bool_t tmp_lsparse=!Lsparse, tmp_rsparse=!Rsparse;
  igraph_sparsemat_t myLsparse, myRsparse, tmpsparse, Rsparse_t;
  int no_of_nodes;
  igraph_real_t evmin, evmax;
  
  /* --------------------------------------------------------------------*/
  /* Argument checks */

  IGRAPH_CHECK(igraph_i_scg_common_checks(graph, matrix, sparsemat,
					  ev, intervals, intervals_vector,
					  evec, 0, groups, scg_graph, 
					  scg_matrix, scg_sparsemat, 
					  /*p=*/ 0, &evmin, &evmax));

  if (graph) { 
    no_of_nodes=igraph_vcount(graph); 
  } else if (matrix) {
    no_of_nodes=igraph_matrix_nrow(matrix);
  } else {
    no_of_nodes=igraph_sparsemat_nrow(sparsemat);
  }

  /* -------------------------------------------------------------------- */
  /* Convert graph, if needed */

  if (graph) {
    mysparsemat=&real_sparsemat;
    IGRAPH_CHECK(igraph_get_sparsemat(graph, mysparsemat));
    IGRAPH_FINALLY(igraph_sparsemat_destroy, mysparsemat);
  }   

  /* -------------------------------------------------------------------- */
  /* Compute eigenpairs, if needed */
  if (tmp_evec) {
    evec=&myevec;
    IGRAPH_MATRIX_INIT_FINALLY(evec, no_of_nodes, no_of_ev);
  }

  if (do_evec || tmp_evec) {
    igraph_arpack_options_t options;
    igraph_eigen_which_t which;
    igraph_matrix_t tmp;
    igraph_vector_t tmpev;
    int i;
    
    which.pos = IGRAPH_EIGEN_SELECT;
    which.il = no_of_nodes-evmax+1;
    which.iu = no_of_nodes-evmin+1;

    IGRAPH_CHECK(igraph_matrix_init(&tmp, no_of_nodes, 
				    which.iu-which.il+1));
    IGRAPH_FINALLY(igraph_matrix_destroy, &tmp);
    IGRAPH_CHECK(igraph_eigen_matrix_symmetric(matrix, mysparsemat,
					       /* fun= */ 0,
					       /* extra= */ 0, 
					       /* algorithm= */ 
					       use_arpack ? 
					       IGRAPH_EIGEN_ARPACK : 
					       IGRAPH_EIGEN_LAPACK, &which, 
					       &options, /* values= */ 0, 
					       &tmp));
    IGRAPH_VECTOR_INIT_FINALLY(&tmpev, no_of_ev);
    for (i=0; i<no_of_ev; i++) {
      VECTOR(tmpev)[i] = evmax - VECTOR(*ev)[i];
    }
    IGRAPH_CHECK(igraph_matrix_select_cols(&tmp, evec, &tmpev));
    igraph_vector_destroy(&tmpev);
    igraph_matrix_destroy(&tmp);
    IGRAPH_FINALLY_CLEAN(2);
  }

  /* -------------------------------------------------------------------- */
  /* Work out groups, if needed */
  if (tmp_groups) {
    groups=&mygroups;
    IGRAPH_VECTOR_INIT_FINALLY((igraph_vector_t*)groups, no_of_nodes);
  }
  if (do_groups) {
    IGRAPH_CHECK(igraph_scg_grouping(evec, (igraph_vector_t*)groups, 
				     intervals, intervals_vector, 
				     IGRAPH_SCG_SYMMETRIC, algorithm, 
				     /*p=*/ 0, maxiter));
  }
  
  /* -------------------------------------------------------------------- */
  /* Perform coarse graining */
  if (tmp_lsparse) {
    Lsparse=&myLsparse;
  }
  if (tmp_rsparse) {
    Rsparse=&myRsparse;
  }
  IGRAPH_CHECK(igraph_scg_semiprojectors(groups, IGRAPH_SCG_SYMMETRIC,
					 L, R, Lsparse, Rsparse, /*p=*/ 0,
					 IGRAPH_SCG_NORM_ROW));
  if (tmp_groups) {
    igraph_vector_destroy((igraph_vector_t*) groups);
    IGRAPH_FINALLY_CLEAN(1);
  }
  if (tmp_evec) {
    igraph_matrix_destroy(evec);
    IGRAPH_FINALLY_CLEAN(1);
  }
  if (Rsparse) { IGRAPH_FINALLY(igraph_sparsemat_destroy, Rsparse); }
  if (Lsparse) { IGRAPH_FINALLY(igraph_sparsemat_destroy, Lsparse); }

  /* -------------------------------------------------------------------- */
  /* Compute coarse grained matrix/graph/sparse matrix */
  IGRAPH_CHECK(igraph_sparsemat_compress(Rsparse, &tmpsparse));
  IGRAPH_FINALLY(igraph_sparsemat_destroy, &tmpsparse);
  IGRAPH_CHECK(igraph_sparsemat_transpose(&tmpsparse, &Rsparse_t,
					  /*values=*/ 1));
  igraph_sparsemat_destroy(&tmpsparse);
  IGRAPH_FINALLY_CLEAN(1);
  IGRAPH_FINALLY(igraph_sparsemat_destroy, &Rsparse_t);  

  IGRAPH_CHECK(igraph_i_scg_get_result(IGRAPH_SCG_SYMMETRIC,
				       matrix, mysparsemat,
				       Lsparse, &Rsparse_t,
				       scg_graph, scg_matrix, 
				       scg_sparsemat, /*directed=*/ 0));

  /* -------------------------------------------------------------------- */
  /* Clean up */
  
  igraph_sparsemat_destroy(&Rsparse_t);
  IGRAPH_FINALLY_CLEAN(1);
  if (Lsparse) { IGRAPH_FINALLY_CLEAN(1); }
  if (Rsparse) { IGRAPH_FINALLY_CLEAN(1); }

  if (graph) {
    igraph_sparsemat_destroy(mysparsemat);
    IGRAPH_FINALLY_CLEAN(1);
  }

  return 0;
}

int igraph_scg_stochastic(const igraph_t *graph,
			  const igraph_matrix_t *matrix,
			  const igraph_sparsemat_t *sparsemat,
			  const igraph_vector_t *ev,
			  igraph_integer_t intervals,
			  const igraph_vector_t *intervals_vector,
			  igraph_scg_algorithm_t algorithm,
			  igraph_scg_norm_t norm,
			  igraph_vector_complex_t *eval,
			  igraph_matrix_complex_t *evec,
			  igraph_vector_t *groups,
			  igraph_vector_t *p,
			  igraph_bool_t use_arpack,
			  igraph_integer_t maxiter,
			  igraph_t *scg_graph,
			  igraph_matrix_t *scg_matrix,
			  igraph_sparsemat_t *scg_sparsemat,
			  igraph_matrix_t *L,
			  igraph_matrix_t *R,
			  igraph_sparsemat_t *Lsparse,
			  igraph_sparsemat_t *Rsparse) {

  igraph_matrix_t *mymatrix=(igraph_matrix_t*) matrix, real_matrix;
  igraph_sparsemat_t *mysparsemat=(igraph_sparsemat_t*) sparsemat, 
    real_sparsemat;
  int no_of_nodes;
  igraph_real_t evmin, evmax;
  igraph_arpack_options_t options;
  igraph_eigen_which_t which;
  /* eigenvectors are calculated and returned */
  igraph_bool_t do_evec= evec && igraph_matrix_complex_size(evec)==0;
  /* groups are calculated */
  igraph_bool_t do_groups= !groups || igraph_vector_size(groups)==0;
  igraph_bool_t tmp_groups= !groups;
  /* eigenvectors are not returned but must be calculated for groups */
  igraph_bool_t tmp_evec= !do_evec && do_groups;
  igraph_matrix_complex_t myevec;
  igraph_vector_t mygroups;
  igraph_bool_t do_p= p && igraph_vector_size(p)==0;
  igraph_vector_t *myp=(igraph_vector_t *) p;
  int no_of_ev=igraph_vector_size(ev);
  igraph_bool_t tmp_lsparse=!Lsparse, tmp_rsparse=!Rsparse;
  igraph_sparsemat_t myLsparse, myRsparse, tmpsparse, Rsparse_t;

  /* --------------------------------------------------------------------*/
  /* Argument checks */

  IGRAPH_CHECK(igraph_i_scg_common_checks(graph, matrix, sparsemat,
					  ev, intervals, intervals_vector,
					  0, evec, groups, scg_graph, 
					  scg_matrix, scg_sparsemat, p,
					  &evmin, &evmax));
  
  if (graph) { 
    no_of_nodes=igraph_vcount(graph); 
  } else if (matrix) {
    no_of_nodes=igraph_matrix_nrow(matrix);
  } else {
    no_of_nodes=igraph_sparsemat_nrow(sparsemat);
  }

  /* -------------------------------------------------------------------- */
  /* Convert graph, if needed */

  if (graph) {
    mysparsemat=&real_sparsemat;
    IGRAPH_CHECK(igraph_get_stochastic_sparsemat(graph, mysparsemat, norm));
    IGRAPH_FINALLY(igraph_sparsemat_destroy, mysparsemat);
  } else if (matrix) {
    mymatrix=&real_matrix;
    IGRAPH_CHECK(igraph_i_matrix_stochastic(matrix, mymatrix, norm));
    IGRAPH_FINALLY(igraph_matrix_destroy, mymatrix);
  } else { /* sparsemat */
    mysparsemat=&real_sparsemat;
    IGRAPH_CHECK(igraph_i_sparsemat_stochastic(sparsemat, mysparsemat, norm));
    IGRAPH_FINALLY(igraph_sparsemat_destroy, mysparsemat);
  }

  /* -------------------------------------------------------------------- */
  /* Compute eigenpairs, if needed */

  if (tmp_evec) {
    evec=&myevec;
    IGRAPH_CHECK(igraph_matrix_complex_init(evec, no_of_nodes, no_of_ev));
    IGRAPH_FINALLY(igraph_matrix_complex_destroy, evec);
  }

  if (do_evec || tmp_evec) {
    igraph_matrix_complex_t tmp;
    igraph_vector_t tmpev;
    int i;

    which.pos = IGRAPH_EIGEN_SELECT;
    which.il = no_of_nodes-evmax+1;
    which.iu = no_of_nodes-evmin+1;

    IGRAPH_CHECK(igraph_matrix_complex_init(&tmp, no_of_nodes, 
					    which.iu-which.il+1));
    IGRAPH_FINALLY(igraph_matrix_complex_destroy, &tmp);
    IGRAPH_CHECK(igraph_eigen_matrix(mymatrix, mysparsemat, /*fun=*/ 0,
				     /*extra=*/ 0, use_arpack ? 
				     IGRAPH_EIGEN_ARPACK : 
				     IGRAPH_EIGEN_LAPACK, &which, &options, 
				     /*values=*/ 0, &tmp));
    
    IGRAPH_VECTOR_INIT_FINALLY(&tmpev, no_of_ev);
    for (i=0; i<no_of_ev; i++) {
      VECTOR(tmpev)[i] = evmax - VECTOR(*ev)[i];
    }
    IGRAPH_CHECK(igraph_matrix_complex_select_cols(&tmp, evec, &tmpev));
    igraph_vector_destroy(&tmpev);
    igraph_matrix_complex_destroy(&tmp);
    IGRAPH_FINALLY_CLEAN(2);
  }

  /* Compute p if not supplied */
  if (do_p) {
    igraph_eigen_which_t w;
    igraph_matrix_complex_t tmp;
    igraph_arpack_options_t o;
    igraph_matrix_t trans, *mytrans=&trans;
    igraph_sparsemat_t sparse_trans, *mysparse_trans=&sparse_trans;
    int i;
    igraph_arpack_options_init(&o);
    IGRAPH_CHECK(igraph_matrix_complex_init(&tmp, 0, 0));
    IGRAPH_FINALLY(igraph_matrix_complex_destroy, &tmp);
    w.pos=IGRAPH_EIGEN_LR;
    w.howmany=1;

    if (mymatrix) { 
      IGRAPH_CHECK(igraph_matrix_copy(&trans, mymatrix));
      IGRAPH_FINALLY(igraph_matrix_destroy, &trans);
      IGRAPH_CHECK(igraph_matrix_transpose(&trans));
      mysparse_trans=0;
    } else { 
      IGRAPH_CHECK(igraph_sparsemat_transpose(mysparsemat, &sparse_trans, 
					      /*values=*/ 1));
      IGRAPH_FINALLY(igraph_sparsemat_destroy, mysparse_trans);
      mytrans=0;
    }

    IGRAPH_CHECK(igraph_eigen_matrix(mytrans, mysparse_trans, /*fun=*/ 0, 
				     /*extra=*/ 0, /*algorith=*/ 
				     use_arpack ? 
				     IGRAPH_EIGEN_ARPACK : 
				     IGRAPH_EIGEN_LAPACK, &w, &o, 
				     /*values=*/ 0, &tmp));

    if (mymatrix) {
      igraph_matrix_destroy(&trans);
      IGRAPH_FINALLY_CLEAN(1);
    } else {
      igraph_sparsemat_destroy(mysparse_trans);
      IGRAPH_FINALLY_CLEAN(1);
    }

    IGRAPH_CHECK(igraph_vector_resize(myp, no_of_nodes));
    for (i=0; i<no_of_nodes; i++) {
      VECTOR(*myp)[i] = fabs(IGRAPH_REAL(MATRIX(tmp, i, 0)));
    }
    igraph_matrix_complex_destroy(&tmp);
    IGRAPH_FINALLY_CLEAN(1);
  }    

  /* -------------------------------------------------------------------- */
  /* Work out groups, if needed */
  /* TODO: use complex part as well */
  if (tmp_groups) {
    groups=&mygroups;
    IGRAPH_VECTOR_INIT_FINALLY((igraph_vector_t*)groups, no_of_nodes);
  }
  if (do_groups) {
    igraph_matrix_t tmp;
    IGRAPH_MATRIX_INIT_FINALLY(&tmp, 0, 0);
    IGRAPH_CHECK(igraph_matrix_complex_real(evec, &tmp));
    IGRAPH_CHECK(igraph_scg_grouping(&tmp, (igraph_vector_t*)groups, 
				     intervals, intervals_vector, 
				     IGRAPH_SCG_STOCHASTIC, algorithm, 
				     myp, maxiter));
    igraph_matrix_destroy(&tmp);
    IGRAPH_FINALLY_CLEAN(1);
  }
  
  /* -------------------------------------------------------------------- */
  /* Perform coarse graining */
  if (tmp_lsparse) {
    Lsparse=&myLsparse;
  }
  if (tmp_rsparse) {
    Rsparse=&myRsparse;
  }
  IGRAPH_CHECK(igraph_scg_semiprojectors(groups, IGRAPH_SCG_STOCHASTIC,
					 L, R, Lsparse, Rsparse, myp, norm));
  if (tmp_groups) {
    igraph_vector_destroy((igraph_vector_t*) groups);
    IGRAPH_FINALLY_CLEAN(1);
  }
  if (tmp_evec) {
    igraph_matrix_complex_destroy(evec);
    IGRAPH_FINALLY_CLEAN(1);
  }
  if (Rsparse) { IGRAPH_FINALLY(igraph_sparsemat_destroy, Rsparse); }
  if (Lsparse) { IGRAPH_FINALLY(igraph_sparsemat_destroy, Lsparse); }

  /* -------------------------------------------------------------------- */
  /* Compute coarse grained matrix/graph/sparse matrix */
  IGRAPH_CHECK(igraph_sparsemat_compress(Rsparse, &tmpsparse));
  IGRAPH_FINALLY(igraph_sparsemat_destroy, &tmpsparse);
  IGRAPH_CHECK(igraph_sparsemat_transpose(&tmpsparse, &Rsparse_t,
					  /*values=*/ 1));
  igraph_sparsemat_destroy(&tmpsparse);
  IGRAPH_FINALLY_CLEAN(1);
  IGRAPH_FINALLY(igraph_sparsemat_destroy, &Rsparse_t);  

  IGRAPH_CHECK(igraph_i_scg_get_result(IGRAPH_SCG_STOCHASTIC,
				       mymatrix, mysparsemat,
				       Lsparse, &Rsparse_t,
				       scg_graph, scg_matrix, 
				       scg_sparsemat, /*directed=*/ 0));

  /* -------------------------------------------------------------------- */
  /* Clean up */

  igraph_sparsemat_destroy(&Rsparse_t);
  IGRAPH_FINALLY_CLEAN(1);
  if (Lsparse) { IGRAPH_FINALLY_CLEAN(1); }
  if (Rsparse) { IGRAPH_FINALLY_CLEAN(1); }  

  if (graph) {
    igraph_sparsemat_destroy(mysparsemat);
    IGRAPH_FINALLY_CLEAN(1);
  } else if (matrix) { 
    igraph_matrix_destroy(mymatrix);
    IGRAPH_FINALLY_CLEAN(1);
  } else {
    igraph_sparsemat_destroy(mysparsemat);
    IGRAPH_FINALLY_CLEAN(1);
  }

  return 0;
}

int igraph_scg_laplacian(const igraph_t *graph,
			 const igraph_matrix_t *matrix,
			 const igraph_sparsemat_t *sparsemat,
			 const igraph_vector_t *ev,
			 igraph_integer_t intervals,
			 const igraph_vector_t *intervals_vector,
			 igraph_scg_algorithm_t algorithm,
			 igraph_scg_norm_t norm,
			 igraph_scg_direction_t direction,
			 igraph_vector_complex_t *eval,
			 igraph_matrix_complex_t *evec,
			 igraph_vector_t *groups,
			 igraph_bool_t use_arpack,
			 igraph_integer_t maxiter,
			 igraph_t *scg_graph,
			 igraph_matrix_t *scg_matrix,
			 igraph_sparsemat_t *scg_sparsemat,
			 igraph_matrix_t *L,
			 igraph_matrix_t *R,
			 igraph_sparsemat_t *Lsparse,
			 igraph_sparsemat_t *Rsparse) {

  igraph_matrix_t *mymatrix=(igraph_matrix_t*) matrix, real_matrix;
  igraph_sparsemat_t *mysparsemat=(igraph_sparsemat_t*) sparsemat, 
    real_sparsemat;
  int no_of_nodes;
  igraph_real_t evmin, evmax;
  igraph_arpack_options_t options;
  igraph_eigen_which_t which;
  /* eigenvectors are calculated and returned */
  igraph_bool_t do_evec= evec && igraph_matrix_complex_size(evec)==0;
  /* groups are calculated */
  igraph_bool_t do_groups= !groups || igraph_vector_size(groups)==0;
  igraph_bool_t tmp_groups= !groups;
  /* eigenvectors are not returned but must be calculated for groups */
  igraph_bool_t tmp_evec= !do_evec && do_groups;
  igraph_matrix_complex_t myevec;
  igraph_vector_t mygroups;
  int no_of_ev=igraph_vector_size(ev);
  igraph_bool_t tmp_lsparse=!Lsparse, tmp_rsparse=!Rsparse;
  igraph_sparsemat_t myLsparse, myRsparse, tmpsparse, Rsparse_t;

  /* --------------------------------------------------------------------*/
  /* Argument checks */

  IGRAPH_CHECK(igraph_i_scg_common_checks(graph, matrix, sparsemat,
					  ev, intervals, intervals_vector,
					  0, evec, groups, scg_graph, 
					  scg_matrix, scg_sparsemat, /*p=*/ 0,
					  &evmin, &evmax));
  
  if (graph) { 
    no_of_nodes=igraph_vcount(graph); 
  } else if (matrix) {
    no_of_nodes=igraph_matrix_nrow(matrix);
  } else {
    no_of_nodes=igraph_sparsemat_nrow(sparsemat);
  }

  /* -------------------------------------------------------------------- */
  /* Convert graph, if needed, get Laplacian matrix */

  if (graph) {
    mysparsemat=&real_sparsemat;
    IGRAPH_CHECK(igraph_laplacian(graph, 0, mysparsemat, /*normalized=*/ 0, 
				  /*weights=*/ 0));
    IGRAPH_FINALLY(igraph_sparsemat_destroy, mysparsemat);
  } else if (matrix) { 
    mymatrix=&real_matrix;
    IGRAPH_CHECK(igraph_i_matrix_laplacian(matrix, mymatrix, norm));
    IGRAPH_FINALLY(igraph_matrix_destroy, mymatrix);
  } else { /* sparsemat */
    mysparsemat=&real_sparsemat;
    IGRAPH_CHECK(igraph_i_sparsemat_laplacian(sparsemat, mysparsemat, norm));
    IGRAPH_FINALLY(igraph_sparsemat_destroy, mysparsemat);
  }

  /* -------------------------------------------------------------------- */
  /* Compute eigenpairs, if needed */

  if (tmp_evec) {
    evec=&myevec;
    IGRAPH_CHECK(igraph_matrix_complex_init(evec, no_of_nodes, no_of_ev));
    IGRAPH_FINALLY(igraph_matrix_complex_destroy, evec);
  }

  if (do_evec || tmp_evec) {
    igraph_matrix_complex_t tmp;
    igraph_vector_t tmpev;
    int i;

    which.pos = IGRAPH_EIGEN_SELECT;
    which.il = no_of_nodes-evmax+1;
    which.iu = no_of_nodes-evmin+1;

    IGRAPH_CHECK(igraph_matrix_complex_init(&tmp, no_of_nodes, 
					    which.iu-which.il+1));
    IGRAPH_FINALLY(igraph_matrix_complex_destroy, &tmp);
    IGRAPH_CHECK(igraph_eigen_matrix(mymatrix, mysparsemat, /*fun=*/ 0,
				     /*extra=*/ 0, use_arpack ? 
				     IGRAPH_EIGEN_ARPACK : 
				     IGRAPH_EIGEN_LAPACK, &which, &options, 
				     /*values=*/ 0, &tmp));
    
    IGRAPH_VECTOR_INIT_FINALLY(&tmpev, no_of_ev);
    for (i=0; i<no_of_ev; i++) {
      VECTOR(tmpev)[i] = evmax - VECTOR(*ev)[i];
    }
    IGRAPH_CHECK(igraph_matrix_complex_select_cols(&tmp, evec, &tmpev));
    igraph_vector_destroy(&tmpev);
    igraph_matrix_complex_destroy(&tmp);
    IGRAPH_FINALLY_CLEAN(2);
  }

  /* -------------------------------------------------------------------- */
  /* Work out groups, if needed */
  /* TODO: use complex part as well */
  if (tmp_groups) {
    groups=&mygroups;
    IGRAPH_VECTOR_INIT_FINALLY((igraph_vector_t*)groups, no_of_nodes);
  }
  if (do_groups) {
    igraph_matrix_t tmp;
    IGRAPH_MATRIX_INIT_FINALLY(&tmp, 0, 0);
    IGRAPH_CHECK(igraph_matrix_complex_real(evec, &tmp));
    IGRAPH_CHECK(igraph_scg_grouping(&tmp, (igraph_vector_t*)groups, 
				     intervals, intervals_vector, 
				     IGRAPH_SCG_LAPLACIAN, algorithm, 
				     /*p=*/ 0, maxiter));
    igraph_matrix_destroy(&tmp);
    IGRAPH_FINALLY_CLEAN(1);
  }
  
  /* -------------------------------------------------------------------- */
  /* Perform coarse graining */
  if (tmp_lsparse) {
    Lsparse=&myLsparse;
  }
  if (tmp_rsparse) {
    Rsparse=&myRsparse;
  }
  IGRAPH_CHECK(igraph_scg_semiprojectors(groups, IGRAPH_SCG_LAPLACIAN,
					 L, R, Lsparse, Rsparse, /*p=*/ 0, 
					 norm));
  if (tmp_groups) {
    igraph_vector_destroy((igraph_vector_t*) groups);
    IGRAPH_FINALLY_CLEAN(1);
  }
  if (tmp_evec) {
    igraph_matrix_complex_destroy(evec);
    IGRAPH_FINALLY_CLEAN(1);
  }
  if (Rsparse) { IGRAPH_FINALLY(igraph_sparsemat_destroy, Rsparse); }
  if (Lsparse) { IGRAPH_FINALLY(igraph_sparsemat_destroy, Lsparse); }

  /* -------------------------------------------------------------------- */
  /* Compute coarse grained matrix/graph/sparse matrix */
  IGRAPH_CHECK(igraph_sparsemat_compress(Rsparse, &tmpsparse));
  IGRAPH_FINALLY(igraph_sparsemat_destroy, &tmpsparse);
  IGRAPH_CHECK(igraph_sparsemat_transpose(&tmpsparse, &Rsparse_t,
					  /*values=*/ 1));
  igraph_sparsemat_destroy(&tmpsparse);
  IGRAPH_FINALLY_CLEAN(1);
  IGRAPH_FINALLY(igraph_sparsemat_destroy, &Rsparse_t);  

  IGRAPH_CHECK(igraph_i_scg_get_result(IGRAPH_SCG_LAPLACIAN,
				       mymatrix, mysparsemat,
				       Lsparse, &Rsparse_t,
				       scg_graph, scg_matrix, 
				       scg_sparsemat, /*directed=*/ 1));

  /* -------------------------------------------------------------------- */
  /* Clean up */

  igraph_sparsemat_destroy(&Rsparse_t);
  IGRAPH_FINALLY_CLEAN(1);
  if (Lsparse) { IGRAPH_FINALLY_CLEAN(1); }
  if (Rsparse) { IGRAPH_FINALLY_CLEAN(1); }  

  if (graph) {
    igraph_sparsemat_destroy(mysparsemat);
    IGRAPH_FINALLY_CLEAN(1);
  } else if (matrix) {
    igraph_matrix_destroy(mymatrix);
    IGRAPH_FINALLY_CLEAN(1);
  } else {
    igraph_sparsemat_destroy(mysparsemat);
    IGRAPH_FINALLY_CLEAN(1);
  }

  return 0;
}
