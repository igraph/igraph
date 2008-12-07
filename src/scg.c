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

typedef struct {
  igraph_matrix_t *mat;
} igraph_i_scg_t;

int igraph_i_scg(igraph_real_t *to, const igraph_real_t *from,
		 long int n, void *extra) {

  igraph_matrix_t *mat= ((igraph_i_scg_t*)extra) -> mat;
  long int i, j;
  igraph_real_t sum;

  for (i=0; i<n; i++) {
    sum=0.0;
    for (j=0; j<n; j++) {
      sum += MATRIX(*mat,i,j) * from[j];
    }
    to[i] = sum;
  }
  
  return 0;
}

int igraph_scg(const igraph_t *graph, 
	       igraph_t *res_graph,
	       igraph_matrix_t *res_matrix,
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
	       igraph_integer_t maxiter
	       /* igraph_?_t *semproj, */
	       /* igraph_?_t *epairs, */
	       /* igraph_?_t *c_markovp, */) {
  
  
  
  return 0;
}

int igraph_i_scg_reorder_arpack(igraph_matrix_t *vectors, 
				const igraph_vector_long_t *order, 
				const igraph_vector_t *ev1) {
  /* TODO */
  return 0;
}


int igraph_scg_matrix(const igraph_matrix_t *matrix, 
		      igraph_t *res_graph,
		      igraph_matrix_t *res_matrix,
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
		      igraph_integer_t maxiter
		      /* igraph_?_t *semproj, */
		      /* igraph_?_t *epairs, */
		      /* igraph_?_t *c_markovp, */) {

  long int n=igraph_matrix_nrow(matrix);
  long int nnt=igraph_vector_size(nt);
  long int nev=igraph_vector_size(ev);
  igraph_vector_t *mynt=(igraph_vector_t*)nt, mynt_v;
  igraph_bool_t sym, isTransposedX;
  igraph_matrix_t *X=(igraph_matrix_t*)matrix, Xm;

  igraph_vector_t values;
  igraph_matrix_t values2, vectors;

  igraph_vector_t *mygroup=group, mygroup_v;
  
  /************** Check arguments ***********/
  
  if (n != igraph_matrix_ncol(matrix)) {
    IGRAPH_ERROR("The input matrix must be square", IGRAPH_NONSQUARE);
  }
  
  if (igraph_vector_size(ev)==0) {
    IGRAPH_ERROR("Empty `ev' vector", IGRAPH_EINVAL);
  }

  if (igraph_vector_max(ev) > n || igraph_vector_min(ev) <= 0) { 
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
  
  if (evec && 
      (igraph_matrix_nrow(evec) != n || igraph_matrix_ncol(evec) != nev)) {
    IGRAPH_ERROR("Invalid `evec' size, must be n x nev", IGRAPH_EINVAL);
  }

  if (matrix_type == IGRAPH_SCG_MATRIX_STOCHASTIC && markovp && 
      igraph_vector_size(markovp) != n) {
    IGRAPH_ERROR("Stationary Markov probability vector (`markovp') must be "
		 "NULL or of length `n', the order of the matrix/graph", 
		 IGRAPH_EINVAL);
  }

  if (group && igraph_vector_size(group) != n) {
    IGRAPH_ERROR("`group' must be either NULL or of length `n'", IGRAPH_EINVAL);
  }
  
  /*******************************************/
  
  sym = igraph_matrix_is_symmetric(matrix);

  isTransposedX = (direction == IGRAPH_SCG_DIR_LEFT);
  isTransposedX = isTransposedX || 
    (direction==IGRAPH_SCG_DIR_DEFAULT && 
     (matrix_type==IGRAPH_SCG_MATRIX_LAPLACIAN || 
      matrix_type==IGRAPH_SCG_MATRIX_STOCHASTIC) && norm==IGRAPH_SCG_NORM_COL);
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
  
  if (!group) {

    if (!evec) {
      
      IGRAPH_MATRIX_INIT_FINALLY(&vectors, 0, 0);
      
      if (use_arpack) {
	igraph_vector_t ev1, ev2;
	long int i;
	igraph_arpack_options_t arpack_opts;
	igraph_i_scg_t extra;
	extra.mat=X;

	IGRAPH_VECTOR_INIT_FINALLY(&ev1, 0);
	IGRAPH_VECTOR_INIT_FINALLY(&ev2, 0);
	for (i=0; i<nev; i++) {
	  long int actev=VECTOR(*ev)[i];
	    if (actev > n/2) {
	      IGRAPH_CHECK(igraph_vector_push_back(&ev1, actev));
	    } else {
	      IGRAPH_CHECK(igraph_vector_push_back(&ev2, actev));
	    }
	}

	if (sym) {
	  IGRAPH_VECTOR_INIT_FINALLY(&values, 0);
	} else { 
	  IGRAPH_MATRIX_INIT_FINALLY(&values2, 0, 0);
	}

	if (igraph_vector_size(&ev1)>0) {
	  igraph_arpack_options_init(&arpack_opts);
	  arpack_opts.n=n;
	  arpack_opts.nev=igraph_vector_max(&ev1);
	  if (sym) {
	    igraph_vector_long_t order;
	    arpack_opts.which[0]='L'; arpack_opts.which[1]='A';
	    arpack_opts.ncv=2*igraph_vector_max(&ev1)+1;
	    igraph_arpack_rssolve(igraph_i_scg, &extra,
				  &arpack_opts, /*storage=*/ 0, 
				  &values, &vectors);
	    IGRAPH_CHECK(igraph_vector_long_init(&order, n));
	    igraph_vector_sort_order(&values, &order, /*reverse=*/ 1);
	    igraph_i_scg_reorder_arpack(&vectors, &order, &ev1);
	    igraph_vector_long_destroy(&order);
	    IGRAPH_FINALLY_CLEAN(1);
	  } else {
	    arpack_opts.which[0]='L'; arpack_opts.which[1]='M';
	    arpack_opts.ncv=2*igraph_vector_max(&ev1)+2;
	    igraph_arpack_rnsolve(igraph_i_scg, &extra, 
				  &arpack_opts, /*storage=*/ 0,
				  &values2, &vectors);
	    igraph_arpack_unpack_complex(&vectors, &values2,
					 igraph_vector_max(&ev1));
	    
	    /* TODO */
	  }
	  
	  
	}

	if (igraph_vector_size(&ev2)>0) {
	  igraph_arpack_options_init(&arpack_opts);
	  
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

  /* TODO */

  /* ------work out the groups if not supplied----------- */
  
  /* TODO: complex case */
  if (!group) {
    mygroup=&mygroup_v;
    IGRAPH_VECTOR_INIT_FINALLY(mygroup, 0);
    recalculate_group=1;
  }
  
  if (recalculate_group) {
    igraph_scg_grouping(&vectors, group, mynt, matrix_type, markovp, 
			algo, maxiter);
  }

  /* ------perform the coarse graining------------------- */
  
  /* TODO */

  /* ------computes a coarse-grained matrix-------------- */
  
  /* TODO */
  
  /*****************************************************************/

  if (!group) {
    igraph_vector_destroy(mygroup);
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

