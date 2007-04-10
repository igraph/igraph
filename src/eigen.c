/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2007  Gabor Csardi <csardi@rmki.kfki.hu>
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
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA 
   02110-1301 USA

*/

#include "igraph.h"

void *tred2_();
void *tql2_();
void *tred1_();
void *tqlrat_();

int igraph_eigen_tred2(const igraph_matrix_t *A,
		       igraph_vector_t *D,
		       igraph_vector_t *E,
		       igraph_matrix_t *Z) {
  
  long int Arows=igraph_matrix_nrow(A);
  long int Acols=igraph_matrix_ncol(A);
  int nodes=Arows;

  if (Arows != Acols) {
    IGRAPH_ERROR("Invalid matrix", IGRAPH_NONSQUARE);
  }
  if (Arows == Acols && nodes != Arows) {
    IGRAPH_ERROR("Matrix too big, cannot represent its size as an `int'",
		 IGRAPH_EINVAL);
  }

  IGRAPH_CHECK(igraph_vector_resize(D, nodes));
  IGRAPH_CHECK(igraph_vector_resize(E, nodes));
  IGRAPH_CHECK(igraph_matrix_resize(Z, nodes, nodes));
  
  tred2_(&nodes, &nodes, &MATRIX(*A, 0, 0), VECTOR(*D), 
	 VECTOR(*E), &MATRIX(*Z, 0, 0));
  
  return 0;
}

int igraph_eigen_tql2(igraph_vector_t *D,
		      igraph_vector_t *E,
		      igraph_matrix_t *Z) {

  long int Dsize=igraph_vector_size(D);
  long int Esize=igraph_vector_size(E);
  long int Zrows=igraph_matrix_nrow(Z);
  long int Zcols=igraph_matrix_ncol(Z);
  int nodes=Dsize;
  int ierr;
  
  if (Dsize != Esize) {
    IGRAPH_ERROR("Sizes of `D' and `E' don't match", IGRAPH_EINVAL);
  }
  if (Zrows != Zcols) {
    IGRAPH_ERROR("Invalid matrix", IGRAPH_NONSQUARE);
  }
  if (Esize != Zrows) {
    IGRAPH_ERROR("Sizes of `E' and `Z' don't match", IGRAPH_EINVAL);
  }
  if (nodes != Dsize) {
    IGRAPH_ERROR("Matrix too big, cannot represent its size as an `int'",
		 IGRAPH_EINVAL);
  }
  
  tql2_(&nodes, &nodes, VECTOR(*D), VECTOR(*E), &MATRIX(*Z, 0, 0), &ierr);
  
  if (ierr != 0) {
    IGRAPH_ERROR("No convergence after 30 steps", IGRAPH_DIVERGED);
  }
  
  return 0;
}

int igraph_eigen_tred1(const igraph_matrix_t *A,
		       igraph_vector_t *D,
		       igraph_vector_t *E2) {
  
  igraph_vector_t E;
  long int Arows=igraph_matrix_nrow(A);
  long int Acols=igraph_matrix_ncol(A);
  int nodes=Arows;
  
  if (Arows != Acols) {
    IGRAPH_ERROR("Invalid matrix", IGRAPH_NONSQUARE);
  }
  if (nodes != Arows) {
    IGRAPH_ERROR("Matrix too big to represent its size in an `int'",
		 IGRAPH_EINVAL);
  }

  IGRAPH_VECTOR_INIT_FINALLY(&E, nodes);
  
  IGRAPH_CHECK(igraph_vector_resize(D, nodes));
  IGRAPH_CHECK(igraph_vector_resize(E2, nodes));
  
  tred1_(&nodes, &nodes, &MATRIX(*A, 0, 0), VECTOR(*D), 
	 VECTOR(E), VECTOR(*E2));
  
  igraph_vector_destroy(&E);
  IGRAPH_FINALLY_CLEAN(1);
  
  return 0;
}

int igraph_eigen_tqlrat(igraph_vector_t *D,
			igraph_vector_t *E2) {

  long int Dsize=igraph_vector_size(D);
  long int E2size=igraph_vector_size(E2);
  int nodes=Dsize;
  int ierr;
  
  if (Dsize != E2size) {
    IGRAPH_ERROR("Sizes of vectors differ", IGRAPH_EINVAL);
  }
  if (nodes != Dsize) {
    IGRAPH_ERROR("Vector size too big to represent it in an `int'",
		 IGRAPH_EINVAL);
  }
  
  tqlrat_(&nodes, VECTOR(*D), VECTOR(*E2), &ierr);
  
  if (ierr != 0) {
    IGRAPH_ERROR("No convergence after 30 steps", IGRAPH_DIVERGED);
  }
  
  return 0;
}

int igraph_eigen_rs(const igraph_matrix_t *A,
		    igraph_vector_t *values,
		    igraph_matrix_t *vectors) {
  
  long int Arows=igraph_matrix_nrow(A);
  long int Acols=igraph_matrix_ncol(A);
  int nodes=Arows;
  
  if (Arows != Acols) {
    IGRAPH_ERROR("Invalid matrix", IGRAPH_NONSQUARE);
  }
  if (nodes != Arows) {
    IGRAPH_ERROR("Matrix too big, cannot represent its size as an `int'",
		 IGRAPH_EINVAL);
  }
  
  if (!vectors) {
    /* Eigenvalues only */
    igraph_vector_t fv1;
    IGRAPH_VECTOR_INIT_FINALLY(&fv1, nodes);
    igraph_eigen_tred1(A, values, &fv1);
    igraph_eigen_tqlrat(values, &fv1);
    igraph_vector_destroy(&fv1);
    IGRAPH_FINALLY_CLEAN(1);
  } else {
    /* Eigenvectors too */
    igraph_vector_t fv1;
    IGRAPH_VECTOR_INIT_FINALLY(&fv1, nodes);    
    igraph_eigen_tred2(A, values, &fv1, vectors);
    igraph_eigen_tql2(values, &fv1, vectors);
    igraph_vector_destroy(&fv1);
    IGRAPH_FINALLY_CLEAN(1);
  }

  return 0;
}
