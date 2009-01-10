/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2008  Gabor Csardi <Gabor.Csardi@unil.ch>
   University of Lausanne, Rue de Bugnon 27, CH-1005 Lausanne, Switzerland
   
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

#include "igraph.h"

#include <math.h>

int igraph_scg(const igraph_t *graph, 
	       igraph_t *res_graph,
	       igraph_matrix_t *res_matrix,
	       igraph_matrix_t *L, 
	       igraph_matrix_t *R,
	       const igraph_vector_t *ev,
	       const igraph_vector_t *nt,
	       igraph_scg_matrix_t matrix_type, 
	       igraph_scg_algorithm_t algo,
	       igraph_scg_norm_t norm,
	       igraph_scg_direction_t direction,
	       const igraph_matrix_t *evec,
	       const igraph_vector_t *markovp,
	       igraph_vector_t *group,
	       igraph_bool_t recalculate_group,
	       igraph_bool_t use_arpack,
	       igraph_integer_t maxiter) {
  
  /* TODO */
  
  return 0;
}

/* Handle the complex eigenvectors. Basically we just remove 
   imaginary columns that are all zero. We also rewrite oldnt into 
   newnt accordingly. */

int igraph_i_scg_handle_complex(igraph_matrix_t *vectors, 
				const igraph_vector_t *oldnt, 
				igraph_vector_t *newnt) {
  igraph_vector_t keep_cols;
  long int ncol=igraph_matrix_ncol(vectors);
  long int nrow=igraph_matrix_nrow(vectors);
  long int ntlen=igraph_vector_size(oldnt);
  igraph_matrix_t tmp;
  long int i, j;
  
  if (ncol != ntlen * 2) { 
    IGRAPH_ERROR("Internal error, this functions should be called "
		 "for the non-symmetric case only", IGRAPH_EINVAL);
  }

  IGRAPH_VECTOR_INIT_FINALLY(&keep_cols, 0);
  IGRAPH_CHECK(igraph_vector_reserve(&keep_cols, ncol));
  IGRAPH_MATRIX_INIT_FINALLY(&tmp, 0, 0);
  igraph_vector_clear(newnt);

  for (i=0; i<ntlen; i++) {
    /* real part first, always needed */
    igraph_vector_push_back(&keep_cols, 2*i); /* reserved already */
    IGRAPH_CHECK(igraph_vector_push_back(newnt, VECTOR(*oldnt)[i]));
    /* imaginary part, only if not all zero */
    j=0;
    while (j<nrow && MATRIX(*vectors, j, 2*i+1) == 0) {
      j++;
    }
    if (j!=nrow) {
      /* Not all zero */
      igraph_vector_push_back(&keep_cols, 2*i+1); /* reserved already */
      IGRAPH_CHECK(igraph_vector_push_back(newnt, VECTOR(*oldnt)[i]));
    }
  }
  
  IGRAPH_CHECK(igraph_matrix_select_cols(vectors, &tmp, &keep_cols));
  igraph_matrix_destroy(vectors);  *vectors = tmp;
  igraph_vector_destroy(&keep_cols);
  IGRAPH_FINALLY_CLEAN(2);

  return 0;
}

int igraph_scg_matrix(const igraph_matrix_t *matrix, 
		      const igraph_vector_long_t *ev,
		      const igraph_vector_t *nt,
		      igraph_matrix_t *res_matrix,
		      igraph_matrix_t *L, 
		      igraph_matrix_t *R,
		      igraph_matrix_t *evals,
		      igraph_matrix_t *evec,
		      igraph_bool_t recalculate_epairs,
		      igraph_vector_t *group,
		      igraph_bool_t recalculate_group,
		      igraph_scg_matrix_t matrix_type,
		      igraph_scg_algorithm_t algo,
		      igraph_scg_norm_t norm_type,
		      igraph_scg_direction_t direction,
		      const igraph_vector_t *markovp,
		      igraph_bool_t use_arpack,
		      igraph_integer_t maxiter) {

  long int n=igraph_matrix_nrow(matrix);
  long int nnt=igraph_vector_size(nt);
  long int nev=igraph_vector_long_size(ev);
  igraph_vector_t *mynt=(igraph_vector_t*)nt, mynt_v;
  igraph_vector_t *mynt2, mynt2_v;
  igraph_bool_t sym, isTransposedX;
  igraph_matrix_t *X=(igraph_matrix_t*)matrix, Xm;

  igraph_vector_t *mygroup=group, mygroup_v;
  igraph_matrix_t *myL=L, myL_v, *myR=R, myR_v;

  igraph_matrix_t *myevec=evec, my_evec_v;

  /************** Check arguments ***********/
  
  if (n != igraph_matrix_ncol(matrix)) {
    IGRAPH_ERROR("The input matrix must be square", IGRAPH_NONSQUARE);
  }
  
  if (nev==0) {
    IGRAPH_ERROR("Empty `ev' vector", IGRAPH_EINVAL);
  }

  if (igraph_vector_long_max(ev) > n || igraph_vector_long_min(ev) <= 0) { 
    IGRAPH_ERROR("`ev' must contain integers between 1 and the order of the "
		 "graph/matrix", IGRAPH_EINVAL);
  }
  
  if (nnt != nev && nnt != 1) {
    IGRAPH_ERROR("Invalid `nt' length, see docs", IGRAPH_EINVAL);
  }
  
  if (igraph_vector_max(nt) >= n || igraph_vector_min(nt) <= 1) {
    IGRAPH_ERROR("`nt' must be smaller than the order of the graph/matrix "
		 "and greater than 1", IGRAPH_EINVAL);
  }
  
  if (evec && !recalculate_epairs &&
      (igraph_matrix_nrow(evec) != n || igraph_matrix_ncol(evec) != nev)) {
    IGRAPH_ERROR("Invalid `evec' size, must be n x nev", IGRAPH_EINVAL);
  }

  if (matrix_type == IGRAPH_SCG_MATRIX_STOCHASTIC && markovp && 
      igraph_vector_size(markovp) != n) {
    IGRAPH_ERROR("Stationary Markov probability vector (`markovp') must be "
		 "NULL or of length `n', the order of the matrix/graph", 
		 IGRAPH_EINVAL);
  }

  if (group && !recalculate_group && igraph_vector_size(group) != n) {
    IGRAPH_ERROR("`group' must be either NULL or of length `n'", IGRAPH_EINVAL);
  }
  
  /*******************************************/
  
  sym = igraph_matrix_is_symmetric(matrix);

  isTransposedX = (direction == IGRAPH_SCG_DIR_LEFT);
  isTransposedX = isTransposedX || 
    (direction==IGRAPH_SCG_DIR_DEFAULT && 
     (matrix_type==IGRAPH_SCG_MATRIX_LAPLACIAN || 
      matrix_type==IGRAPH_SCG_MATRIX_STOCHASTIC) && norm_type==IGRAPH_SCG_NORM_COL);
  if (isTransposedX) {
    X=&Xm;
    IGRAPH_CHECK(igraph_matrix_copy(X, matrix));
    IGRAPH_FINALLY(igraph_matrix_destroy, X);
    IGRAPH_CHECK(igraph_matrix_transpose(X));
  }  
  
  if (nnt==1) { 
    mynt=&mynt_v;
    IGRAPH_VECTOR_INIT_FINALLY(mynt, nev);
    igraph_vector_fill(mynt, VECTOR(*nt)[0]);
  }

  
  /**************** Compute eigenpairs if not supplied *************/
  
  if (recalculate_group) {

    if (recalculate_epairs) {
      
      if (!evec) {
	myevec=&my_evec_v;
	IGRAPH_MATRIX_INIT_FINALLY(myevec, 0, 0);
      }
      
      if (use_arpack) {
	igraph_arpack_options_t arpack_opts;
	igraph_arpack_options_init(&arpack_opts);
	if (sym) {
	  igraph_arpack_eigen_symmetric_matrix(matrix, ev, evals, myevec, 
					       &arpack_opts);
	} else {
	  igraph_arpack_eigen_matrix(matrix, ev, evals, myevec, 
				     &arpack_opts);
	}
	
      } else {			/* ! use_arpack */
	IGRAPH_ERROR("Currently only arpack is implemented", 
		     IGRAPH_UNIMPLEMENTED);
      }
      
    } /* evec */
      
  } /* group */

  /* ------handle complex eigenvectors if any------------ */
  /* ------for now real and immaginary parts------------- */
  /* ------are partitioned the same way------------------ */

  if (!sym) {
    mynt2=&mynt2_v;
    IGRAPH_VECTOR_INIT_FINALLY(mynt2, 0);
    igraph_i_scg_handle_complex(myevec, mynt, mynt2);
  } else {
    mynt2=mynt;
  }

  /* ------work out the groups if not supplied----------- */
  
  if (!group) {
    mygroup=&mygroup_v;
    IGRAPH_VECTOR_INIT_FINALLY(mygroup, 0);
    recalculate_group=1;
  }
  
  if (recalculate_group) {
    igraph_scg_grouping(myevec, mygroup, mynt2, matrix_type, markovp, 
			algo, maxiter);
  }

  /* ------perform the coarse graining------------------- */

  if (!L) {
    myL=&myL_v;
    IGRAPH_MATRIX_INIT_FINALLY(myL, 0, 0);
  }
  if (!R) {
    myR=&myR_v;
    IGRAPH_MATRIX_INIT_FINALLY(myR, 0, 0); 
  }
  IGRAPH_CHECK(igraph_scg_semi_projectors(mygroup, myL, myR, matrix_type,
					  norm_type, markovp));

  /* ------computes a coarse-grained matrix-------------- */

  if (res_matrix) {
    
    igraph_matrix_t tmp;
    IGRAPH_MATRIX_INIT_FINALLY(&tmp, 0, 0);
    IGRAPH_CHECK(igraph_matrix_tcrossprod(X, myR, &tmp));
    IGRAPH_CHECK(igraph_matrix_mprod(myL, &tmp, res_matrix));
    igraph_matrix_destroy(&tmp);
    IGRAPH_FINALLY_CLEAN(1);

  } /* res_matrix */

  /*****************************************************************/

  if (!R) {
    igraph_matrix_destroy(myR);
    IGRAPH_FINALLY_CLEAN(1);
  }
  
  if (!L) {
    igraph_matrix_destroy(myL);
    IGRAPH_FINALLY_CLEAN(1);
  }

  if (!group) {
    igraph_vector_destroy(mygroup);
    IGRAPH_FINALLY_CLEAN(1);
  }

  if (mynt2 != mynt) {
    igraph_vector_destroy(mynt2);
    IGRAPH_FINALLY_CLEAN(1);
  }

  if (!evec) {
    igraph_matrix_destroy(myevec);
    IGRAPH_FINALLY_CLEAN(1);
  }

  if (mynt != nt) {
    igraph_vector_destroy(mynt);
    IGRAPH_FINALLY_CLEAN(1);
  }

  if (isTransposedX) {
    igraph_matrix_destroy(&Xm);
    IGRAPH_FINALLY_CLEAN(1);
  }

  return 0;
}

int igraph_i_scg_freq(const igraph_vector_t *v,
		      igraph_vector_t *freq,
		      long int max) {
  
  long int i, n=igraph_vector_size(v);
  IGRAPH_CHECK(igraph_vector_resize(freq, max+1));
  for (i=0; i<n; i++) {
    long int j=VECTOR(*v)[i];
    VECTOR(*freq)[j] += 1;
  }

  return 0;
}

int igraph_i_scg_sum(const igraph_vector_t *group, 
		     const igraph_vector_t *p, 
		     igraph_vector_t *sum, 
		     long int max) {
  
  long int i, n=igraph_vector_size(p); /* same as size(group) */
  IGRAPH_CHECK(igraph_vector_resize(sum, max+1));
  igraph_vector_null(sum);
  for (i=0; i<n; i++) {
    long int idx=VECTOR(*group)[i];
    VECTOR(*sum)[idx] += VECTOR(*p)[i];
  }

  return 0;
}

int igraph_scg_semi_projectors(const igraph_vector_t *group,
			       igraph_matrix_t *L,
			       igraph_matrix_t *R,
			       igraph_scg_matrix_t matrix_type,
			       igraph_scg_norm_t norm_type, 
			       const igraph_vector_t *markovp) {

  long int n=igraph_vector_size(group);
  long int m=igraph_vector_max(group);
  long int i, j;
  igraph_vector_t freq, sum, myp;
  igraph_matrix_t *tonorm;

  if (matrix_type == IGRAPH_SCG_MATRIX_STOCHASTIC && 
      (!markovp || igraph_vector_size(markovp) != n)) {
    IGRAPH_ERROR("no stationary probabibities (markovp) and stochastic "
		 "matrix", IGRAPH_EINVAL);
  }  
  
  IGRAPH_CHECK(igraph_matrix_resize(L, m, n));
  IGRAPH_CHECK(igraph_matrix_resize(R, m, n));
  igraph_matrix_null(L);
  
  switch(matrix_type) {

  case IGRAPH_SCG_MATRIX_SYMMETRIC:

    IGRAPH_VECTOR_INIT_FINALLY(&freq, m+1);
    IGRAPH_CHECK(igraph_i_scg_freq(group, &freq, m));

    for (i=0; i<n; i++) {
      long int j=VECTOR(*group)[i]-1;
      igraph_real_t div=sqrt(VECTOR(freq)[j+1]);
      MATRIX(*L, j, i) = 1.0 / div;
    }

    igraph_vector_destroy(&freq);
    IGRAPH_FINALLY_CLEAN(1);

    IGRAPH_CHECK(igraph_matrix_update(R, L));

    break;

  case IGRAPH_SCG_MATRIX_LAPLACIAN:

    for (i=0; i<n; i++) {
      long int j=VECTOR(*group)[i]-1;
      MATRIX(*L, j, i) = 1;
    }
  
    IGRAPH_CHECK(igraph_matrix_update(R, L));
    
    IGRAPH_VECTOR_INIT_FINALLY(&freq, m+1);
    IGRAPH_CHECK(igraph_i_scg_freq(group, &freq, m));
        
    tonorm = norm_type == IGRAPH_SCG_NORM_ROW ? L : R;
    for (i=0; i<n; i++) {
      long int j=VECTOR(*group)[i]-1;
      igraph_real_t div=VECTOR(freq)[j+1];
      MATRIX(*tonorm, j, i) /= div;
    }

    igraph_vector_destroy(&freq);
    IGRAPH_FINALLY_CLEAN(1);
    break;

  case IGRAPH_SCG_MATRIX_STOCHASTIC:
    
    for (i=0; i<n; i++) {
      long int j=VECTOR(*group)[i]-1;
      MATRIX(*L, j, i) = 1;
    }
  
    IGRAPH_CHECK(igraph_matrix_update(R, L));

    IGRAPH_VECTOR_INIT_FINALLY(&myp, n);
    IGRAPH_VECTOR_INIT_FINALLY(&sum, m+1);
    IGRAPH_CHECK(igraph_i_scg_sum(group, markovp, &sum, m));
    for (i=0; i<m; i++) {
      VECTOR(myp)[i] = VECTOR(*markovp)[i] / VECTOR(sum)[i+1];
    }
    
    tonorm = norm_type == IGRAPH_SCG_NORM_ROW ? L : R;
    for (i=0; i<n; i++) {
      igraph_real_t mul=VECTOR(myp)[i];
      MATRIX(*tonorm, j, i) *= mul;
    }

    igraph_vector_destroy(&sum);
    igraph_vector_destroy(&myp);
    IGRAPH_FINALLY_CLEAN(2);
    break;
    
  }

  return 0;
}
