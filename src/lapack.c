/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2010  Gabor Csardi <csardi.gabor@gmail.com>
   Rue de l'Industrie 5, Lausanne 1005, Switzerland
   
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

#include "igraph_lapack.h"
#include "igraph_lapack_internal.h"

int igraph_lapack_dgetrf(igraph_matrix_t *a, igraph_vector_int_t *ipiv, 
			 int *info) {
  int m=igraph_matrix_nrow(a);
  int n=igraph_matrix_ncol(a);
  int lda=m > 0 ? m : 1;
  igraph_vector_int_t *myipiv=ipiv, vipiv;

  if (!ipiv) {
    IGRAPH_CHECK(igraph_vector_int_init(&vipiv, m<n ? m : n));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &vipiv);
    myipiv=&vipiv;
  }

  igraphdgetrf_(&m, &n, VECTOR(a->data), &lda, VECTOR(*myipiv), info);

  if (*info > 0) {
    IGRAPH_WARNING("LU: factor is exactly singular");
  } else if (*info < 0) {
    switch(*info) { 
    case -1:
      IGRAPH_ERROR("Invalid number of rows", IGRAPH_ELAPACK);
      break;
    case -2:
      IGRAPH_ERROR("Invalid number of columns", IGRAPH_ELAPACK);
      break;
    case -3:
      IGRAPH_ERROR("Invalid input matrix", IGRAPH_ELAPACK);
      break;
    case -4:
      IGRAPH_ERROR("Invalid LDA parameter", IGRAPH_ELAPACK);
      break;
    case -5:
      IGRAPH_ERROR("Invalid pivot vector", IGRAPH_ELAPACK);
      break;
    case -6:
      IGRAPH_ERROR("Invalid info argument", IGRAPH_ELAPACK);
      break;
    default:
      IGRAPH_ERROR("Unknown LAPACK error", IGRAPH_ELAPACK);
      break;
    }
  }

  if (!ipiv) {
    igraph_vector_int_destroy(&vipiv);
    IGRAPH_FINALLY_CLEAN(1);
  }
  
  return 0;
}

int igraph_lapack_dgetrs(igraph_bool_t transpose, const igraph_matrix_t *a,
			 igraph_vector_int_t *ipiv, igraph_matrix_t *b, 
			 int *info) {
  char trans = transpose ? 'T' : 'N';
  int n=igraph_matrix_nrow(a);
  int nrhs=igraph_matrix_ncol(b);
  int lda= n > 0 ? n : 1;
  int ldb= n > 0 ? n : 1;
  igraph_vector_int_t *myipiv=ipiv, vipiv;

  if (n != igraph_matrix_ncol(a)) {
    IGRAPH_ERROR("Cannot LU solve matrix", IGRAPH_NONSQUARE);
  }
  if (n != igraph_matrix_nrow(b)) {
    IGRAPH_ERROR("Cannot LU solve matrix, RHS of wrong size", IGRAPH_EINVAL);
  }

  if (!ipiv) {
    IGRAPH_CHECK(igraph_vector_int_init(&vipiv, n));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &vipiv);
    myipiv=&vipiv;
  }
  
  igraphdgetrs_(&trans, &n, &nrhs, VECTOR(a->data), &lda, VECTOR(*myipiv),
		VECTOR(b->data), &ldb, info);

  if (*info < 0) {
    switch(*info) { 
    case -1:
      IGRAPH_ERROR("Invalid transpose argument", IGRAPH_ELAPACK);
      break;
    case -2:
      IGRAPH_ERROR("Invalid number of rows/columns", IGRAPH_ELAPACK);
      break;
    case -3:
      IGRAPH_ERROR("Invalid number of RHS vectors", IGRAPH_ELAPACK);
      break;
    case -4:
      IGRAPH_ERROR("Invalid LU matrix", IGRAPH_ELAPACK);
      break;
    case -5: 
      IGRAPH_ERROR("Invalid LDA parameter", IGRAPH_ELAPACK);
      break;
    case -6:
      IGRAPH_ERROR("Invalid pivot vector", IGRAPH_ELAPACK);
      break;
    case -7:
      IGRAPH_ERROR("Invalid RHS matrix", IGRAPH_ELAPACK);
      break;
    case -8:
      IGRAPH_ERROR("Invalid LDB parameter", IGRAPH_ELAPACK);
      break;
    case -9:
      IGRAPH_ERROR("Invalid info argument", IGRAPH_ELAPACK);
      break;
    default:
      IGRAPH_ERROR("Unknown LAPACK error", IGRAPH_ELAPACK);
      break;
    }
  }
		
  if (!ipiv) {
    igraph_vector_int_destroy(&vipiv);
    IGRAPH_FINALLY_CLEAN(1);
  }

  return 0;
}

int igraph_lapack_dgesv(igraph_matrix_t *a, igraph_vector_int_t *ipiv,
			igraph_matrix_t *b, int *info) {

  int n=igraph_matrix_nrow(a);
  int nrhs=igraph_matrix_ncol(b);
  int lda= n > 0 ? n : 1;
  int ldb= n > 0 ? n : 1;
  igraph_vector_int_t *myipiv=ipiv, vipiv;

  if (n != igraph_matrix_ncol(a)) {
    IGRAPH_ERROR("Cannot LU solve matrix", IGRAPH_NONSQUARE);
  }
  if (n != igraph_matrix_nrow(b)) {
    IGRAPH_ERROR("Cannot LU solve matrix, RHS of wrong size", IGRAPH_EINVAL);
  }

  if (!ipiv) {
    IGRAPH_CHECK(igraph_vector_int_init(&vipiv, n));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &vipiv);
    myipiv=&vipiv;
  }
  
  igraphdgesv_(&n, &nrhs, VECTOR(a->data), &lda, VECTOR(*myipiv),
	       VECTOR(b->data), &ldb, info);

  if (*info > 0) {
    IGRAPH_WARNING("LU: factor is exactly singular");
  } else if (*info < 0) {
    switch(*info) { 
    case -1:
      IGRAPH_ERROR("Invalid number of rows/column", IGRAPH_ELAPACK);
      break;
    case -2:
      IGRAPH_ERROR("Invalid number of RHS vectors", IGRAPH_ELAPACK);
      break;
    case -3:
      IGRAPH_ERROR("Invalid input matrix", IGRAPH_ELAPACK);
      break;
    case -4:
      IGRAPH_ERROR("Invalid LDA parameter", IGRAPH_ELAPACK);
      break;
    case -5:
      IGRAPH_ERROR("Invalid pivot vector", IGRAPH_ELAPACK);
      break;
    case -6:
      IGRAPH_ERROR("Invalid RHS matrix", IGRAPH_ELAPACK);
      break;
    case -7:
      IGRAPH_ERROR("Invalid LDB parameter", IGRAPH_ELAPACK);
      break;
    case -8:
      IGRAPH_ERROR("Invalid info argument", IGRAPH_ELAPACK);
      break;
    default:
      IGRAPH_ERROR("Unknown LAPACK error", IGRAPH_ELAPACK);
      break;
    }
  }
		
  if (!ipiv) {
    igraph_vector_int_destroy(&vipiv);
    IGRAPH_FINALLY_CLEAN(1);
  }
  
  return 0;
}
