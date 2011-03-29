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

/* void grouping(REAL **v, UINT *gr, const UINT n, const UINT *nt, const UINT nev, */
/* 			const UINT matrix, const REAL *p, const UINT algo, const UINT maxiter) */
/* { */
/* 	UINT i,j; */
/* 	UINT **gr_mat = uint_matrix(nev,n); */
/* 	UINT **gr_mat_t = uint_matrix(n,nev); */
	
/* 	switch (algo) */
/* 	{ */
/* 		case 1: */
/* 			for(i=0; i<nev; i++) */
/* 				optimal_partition(v[i], gr_mat[i], n, nt[i], matrix, p);               */
/* 			break; */
			
/* 		case 2: */
/* 			for(i=0; i<nev; i++){ */
/* 				if(!intervals_plus_kmeans(v[i], gr_mat[i], n, nt[i], maxiter)) */
/* 					warning("kmeans did not converge"); */
/* 			} */
/* 			break; */
			
/* 		case 3: */
/* 			for(i=0; i<nev; i++) */
/* 				intervals_method(v[i], gr_mat[i], n, nt[i]); */
/* 			break; */
			
/* 		case 4: */
/* 			for(i=0; i<nev; i++) */
/* 				exact_coarse_graining(v[i],gr_mat[i],n); */
/* 			break;	 */
	 
/* 		default: */
/* 			free_uint_matrix(gr_mat, nev); */
/* 			free_uint_matrix(gr_mat_t, n); */
/* 			error("Choose a grouping method: 1-Optimal, 2-Fixed_size intervals+kmeans\ */
/* 					3-Fixed_size intervals, 4-Exact coarse graining"); */
/* 	} */
	
/* 	//If only one vector copy the groups and jump out */
/*  	if(nev==1){ */
/* 		for(i=0; i<n; i++) */
/* 			gr[i] = gr_mat[0][i]; */
			
/* 		free_uint_matrix(gr_mat, nev); */
/* 		free_uint_matrix(gr_mat_t, n); */
		
/* 		return; */
/* 	} */
	
/* 	//Otherwise works out the final groups as decribed in section 5.4.2 */
/* 		//First, works with the tranpose of gr_mat */
/* 	for(i=0; i<n; i++) */
/* 		for(j=0; j<nev; j++) */
/* 			gr_mat_t[i][j] = gr_mat[j][i]; */
/* 	free_uint_matrix(gr_mat, nev); */
	
/* 		//Then computes the final groups. Use qsort for speed */
/* 	GROUPS *g = (GROUPS*)CALLOC(n, sizeof(GROUPS)); */
/* 	for(i=0; i<n; i++){ */
/* 		g[i].ind = i; */
/* 		g[i].n = nev; */
/* 		g[i].gr = gr_mat_t[i]; */
/* 	} */
		
/* 	qsort(g, n, sizeof(GROUPS), compare_groups); */
/* 	UINT gr_nb = FIRST_GROUP_NB; */
/* 	gr[g[0].ind] = gr_nb; */
/* 	for(i=1; i<n; i++){ */
/* 		if(compare_groups(&g[i], &g[i-1]) != 0) gr_nb++; */
/* 		gr[g[i].ind] = gr_nb; */
/* 	} */
/* 	FREE(g); */
/* 	free_uint_matrix(gr_mat_t,n); */
/* } */



