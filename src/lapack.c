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
			 igraph_vector_int_t *ipiv, igraph_matrix_t *b) {
  char trans = transpose ? 'T' : 'N';
  int n=igraph_matrix_nrow(a);
  int nrhs=igraph_matrix_ncol(b);
  int lda= n > 0 ? n : 1;
  int ldb= n > 0 ? n : 1;
  int info;

  if (n != igraph_matrix_ncol(a)) {
    IGRAPH_ERROR("Cannot LU solve matrix", IGRAPH_NONSQUARE);
  }
  if (n != igraph_matrix_nrow(b)) {
    IGRAPH_ERROR("Cannot LU solve matrix, RHS of wrong size", IGRAPH_EINVAL);
  }

  igraphdgetrs_(&trans, &n, &nrhs, VECTOR(a->data), &lda, VECTOR(*ipiv),
		VECTOR(b->data), &ldb, &info);

  if (info < 0) {
    switch(info) { 
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

int igraph_lapack_dsyevr(const igraph_matrix_t *A, 
			 igraph_lapack_dsyev_which_t which,
			 igraph_real_t vl, igraph_real_t vu, int vestimate, 
			 int il, int iu, igraph_real_t abstol,
			 igraph_vector_t *values, igraph_matrix_t *vectors,
			 igraph_vector_int_t *support) {

  igraph_matrix_t Acopy;
  char jobz = vectors ? 'V' : 'N', range, uplo='U';
  int n=igraph_matrix_nrow(A), lda=n, ldz=n;
  int m, info; 
  igraph_vector_t *myvalues=values, vvalues;
  igraph_vector_int_t *mysupport=support, vsupport;
  igraph_vector_t work;
  igraph_vector_int_t iwork;
  int lwork=-1, liwork=-1;

  if (n != igraph_matrix_ncol(A)) {
    IGRAPH_ERROR("Cannot find eigenvalues/vectors", IGRAPH_NONSQUARE);
  }
  if (which==IGRAPH_LAPACK_DSYEV_INTERVAL && 
      (vestimate < 1 || vestimate > n)) {
    IGRAPH_ERROR("Estimated (upper bound) number of eigenvalues must be "
		 "between 1 and n", IGRAPH_EINVAL);
  }
  if (which==IGRAPH_LAPACK_DSYEV_SELECT && iu-il < 0) {
    IGRAPH_ERROR("Invalied 'il' and/or 'iu' values", IGRAPH_EINVAL);
  }

  IGRAPH_CHECK(igraph_matrix_copy(&Acopy, A));
  IGRAPH_FINALLY(igraph_matrix_destroy, &Acopy);

  IGRAPH_VECTOR_INIT_FINALLY(&work, 1);
  IGRAPH_CHECK(igraph_vector_int_init(&iwork, 1));
  IGRAPH_FINALLY(igraph_vector_int_destroy, &iwork);

  if (!values) {
    IGRAPH_VECTOR_INIT_FINALLY(&vvalues, 0);
    myvalues=&vvalues;
  }
  if (!support) {
    IGRAPH_CHECK(igraph_vector_int_init(&vsupport, 0));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &vsupport);
    mysupport=&vsupport;
  }
  
  switch (which) {
  case IGRAPH_LAPACK_DSYEV_ALL:
    range = 'A';
    IGRAPH_CHECK(igraph_vector_resize(myvalues, n));
    IGRAPH_CHECK(igraph_vector_int_resize(mysupport, 2*n));
    if (vectors) { IGRAPH_CHECK(igraph_matrix_resize(vectors, n, n)); }
    break;
  case IGRAPH_LAPACK_DSYEV_INTERVAL:
    range = 'V';
    IGRAPH_CHECK(igraph_vector_resize(myvalues, vestimate));
    IGRAPH_CHECK(igraph_vector_int_resize(mysupport, 2*vestimate));
    if (vectors) { IGRAPH_CHECK(igraph_matrix_resize(vectors,n, vestimate)); }
   break;
  case IGRAPH_LAPACK_DSYEV_SELECT:
    range = 'I';
    IGRAPH_CHECK(igraph_vector_resize(myvalues, iu-il+1));
    IGRAPH_CHECK(igraph_vector_int_resize(mysupport, 2*(iu-il+1)));
    if (vectors) { IGRAPH_CHECK(igraph_matrix_resize(vectors, n, iu-il+1)); }
    break;
  }
  
  igraphdsyevr_(&jobz, &range, &uplo, &n, &MATRIX(Acopy,0,0), &lda,
		&vl, &vu, &il, &iu, &abstol, &m, VECTOR(*myvalues), 
		vectors ? &MATRIX(*vectors,0,0) : 0, &ldz, VECTOR(*mysupport),
		VECTOR(work), &lwork, VECTOR(iwork), &liwork, &info);
  
  lwork=VECTOR(work)[0];
  liwork=VECTOR(iwork)[0];
  IGRAPH_CHECK(igraph_vector_resize(&work, lwork));
  IGRAPH_CHECK(igraph_vector_int_resize(&iwork, liwork));

  igraphdsyevr_(&jobz, &range, &uplo, &n, &MATRIX(Acopy,0,0), &lda,
		&vl, &vu, &il, &iu, &abstol, &m, VECTOR(*myvalues), 
		vectors ? &MATRIX(*vectors,0,0) : 0, &ldz, VECTOR(*mysupport),
		VECTOR(work), &lwork, VECTOR(iwork), &liwork, &info);

  if (values) { 
    IGRAPH_CHECK(igraph_vector_resize(values, m));
  }
  if (vectors) { 
    IGRAPH_CHECK(igraph_matrix_resize(vectors, n, m));
  }
  if (support) {
    IGRAPH_CHECK(igraph_vector_int_resize(support, m));
  }

  if (!support) {
    igraph_vector_int_destroy(&vsupport);
    IGRAPH_FINALLY_CLEAN(1);
  }
  if (!values) {
    igraph_vector_destroy(&vvalues);
    IGRAPH_FINALLY_CLEAN(1);
  }

  igraph_vector_int_destroy(&iwork);
  igraph_vector_destroy(&work);
  igraph_matrix_destroy(&Acopy);
  IGRAPH_FINALLY_CLEAN(3);
  
  return 0;
}

int igraph_lapack_dgeev(const igraph_matrix_t *A, 
			igraph_vector_t *valuesreal,
			igraph_vector_t *valuesimag, 
			igraph_matrix_t *vectorsleft,
			igraph_matrix_t *vectorsright, 
			int *info) {

  char jobvl= vectorsleft  ? 'V' : 'N';
  char jobvr= vectorsright ? 'V' : 'N';
  int n=igraph_matrix_nrow(A);
  int lda=n, ldvl=n, ldvr=n, lwork=-1;
  igraph_vector_t work;
  igraph_vector_t *myreal=valuesreal, *myimag=valuesimag, vreal, vimag;
  igraph_matrix_t Acopy;
  int error=*info;

  if (igraph_matrix_ncol(A) != n) { 
    IGRAPH_ERROR("Cannot calculate eigenvalues (dgeev)", IGRAPH_NONSQUARE);
  }
  
  IGRAPH_CHECK(igraph_matrix_copy(&Acopy, A));
  IGRAPH_FINALLY(igraph_matrix_destroy, &Acopy);
  
  IGRAPH_VECTOR_INIT_FINALLY(&work, 1);
  
  if (!valuesreal) {
    IGRAPH_VECTOR_INIT_FINALLY(&vreal, n);
    myreal=&vreal;
  } else {
    IGRAPH_CHECK(igraph_vector_resize(myreal, n));
  }
  if (!valuesimag) {
    IGRAPH_VECTOR_INIT_FINALLY(&vimag, n);
    myimag=&vimag;
  } else {
    IGRAPH_CHECK(igraph_vector_resize(myimag, n));
  }
  if (vectorsleft) { 
    IGRAPH_CHECK(igraph_matrix_resize(vectorsleft, n, n));
  }
  if (vectorsright) {
    IGRAPH_CHECK(igraph_matrix_resize(vectorsright, n, n));
  }

  igraphdgeev_(&jobvl, &jobvr, &n, &MATRIX(Acopy,0,0), &lda, 
	       VECTOR(*myreal), VECTOR(*myimag), 
	       vectorsleft  ? &MATRIX(*vectorsleft ,0,0) : 0, &ldvl,
	       vectorsright ? &MATRIX(*vectorsright,0,0) : 0, &ldvr,
	       VECTOR(work), &lwork, info);

  lwork=VECTOR(work)[0];
  IGRAPH_CHECK(igraph_vector_resize(&work, lwork));
  
  igraphdgeev_(&jobvl, &jobvr, &n, &MATRIX(Acopy,0,0), &lda, 
	       VECTOR(*myreal), VECTOR(*myimag), 
	       vectorsleft  ? &MATRIX(*vectorsleft ,0,0) : 0, &ldvl,
	       vectorsright ? &MATRIX(*vectorsright,0,0) : 0, &ldvr,
	       VECTOR(work), &lwork, info);  

  if (*info < 0) {
      IGRAPH_ERROR("Cannot calculate eigenvalues (dgeev)", IGRAPH_ELAPACK);
  } else if (*info > 0) {    
    if (error) {
      IGRAPH_ERROR("Cannot calculate eigenvalues (dgeev)", IGRAPH_ELAPACK);
    } else {
      IGRAPH_WARNING("Cannot calculate eigenvalues (dgeev)");
    }
  }

  igraph_vector_destroy(&work);
  igraph_matrix_destroy(&Acopy);
  IGRAPH_FINALLY_CLEAN(2);

  if (!valuesimag) {
    igraph_vector_destroy(&vimag);
    IGRAPH_FINALLY_CLEAN(1);
  }
  if (!valuesreal) { 
    igraph_vector_destroy(&vreal);
    IGRAPH_FINALLY_CLEAN(1);
  }
  
  return 0;
}
