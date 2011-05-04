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

  if (!intervals_vector) {
    if (intervals <= 1 || intervals >= no_of_nodes) {
      IGRAPH_ERROR("Invalid interval specification", IGRAPH_EINVAL);
    }
  } else {
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
      IGRAPH_CHECK(igraph_i_intervals_plus_kmeans(&MATRIX(*V, 0, i),
						  &MATRIX(gr_mat, 0, i),
						  no_of_nodes, INVEC(i),
						  maxiter));
    }
    break;
  case IGRAPH_SCG_INTERV:
    for (i=0; i<nev; i++) {
      IGRAPH_CHECK(igraph_i_intervals_method(&MATRIX(*V, 0, i),
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
    GROUPS *g = (GROUPS*)CALLOC(no_of_nodes, sizeof(GROUPS));
    IGRAPH_CHECK(igraph_matrix_int_transpose(&gr_mat));
    for(i=0; i<no_of_nodes; i++){
      g[i].ind = i;
      g[i].n = nev;
      g[i].gr = &MATRIX(gr_mat, 0, i);
    }
		
    qsort(g, no_of_nodes, sizeof(GROUPS), igraph_i_compare_groups);
    UINT gr_nb = FIRST_GROUP_NB;
    VECTOR(*groups)[g[0].ind] = gr_nb;
    for(i=1; i<no_of_nodes; i++){
      if(igraph_i_compare_groups(&g[i], &g[i-1]) != 0) gr_nb++;
      VECTOR(*groups)[g[i].ind] = gr_nb;
    }
    FREE(g);
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
    VECTOR(tab)[i] = sqrt(VECTOR(tab)[i]);
  }

  if (norm == IGRAPH_SCG_NORM_ROW) {
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
  } else {
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
