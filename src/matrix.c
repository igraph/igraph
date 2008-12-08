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

#include "types.h"

#define BASE_IGRAPH_REAL
#include "igraph_pmt.h"
#include "matrix.pmt"
#include "igraph_pmt_off.h"
#undef BASE_IGRAPH_REAL

#define BASE_LONG
#include "igraph_pmt.h"
#include "matrix.pmt"
#include "igraph_pmt_off.h"
#undef BASE_LONG

#define BASE_CHAR
#include "igraph_pmt.h"
#include "matrix.pmt"
#include "igraph_pmt_off.h"
#undef BASE_CHAR

#define BASE_BOOL
#include "igraph_pmt.h"
#include "matrix.pmt"
#include "igraph_pmt_off.h"
#undef BASE_BOOL

/*--------------------------------------------------------------*/
/* Product of two matrices, these are written only for 'double' */
/*--------------------------------------------------------------*/

int igraphdgemm_(char *transa, char *transb, long int *m, long int *n,
		 long int *k, igraph_real_t *alpha, igraph_real_t *a, 
		 long int *lda, igraph_real_t *b, long int *ldb, igraph_real_t *beta,
		 igraph_real_t *c, long int *ldc);

int igraph_matrix_dgemm(const igraph_matrix_t *m1,
			const igraph_matrix_t *m2,
			igraph_matrix_t *res, 
			igraph_real_t alpha,
			igraph_real_t beta,
			igraph_bool_t transpose_m1,
			igraph_bool_t transpose_m2) {
  
  long int nrow1=igraph_matrix_nrow(m1);
  long int ncol1=igraph_matrix_ncol(m1);
  long int nrow2=igraph_matrix_nrow(m2);
  long int ncol2=igraph_matrix_ncol(m2);
  int ret;
  char t1 = transpose_m1 ? 't' : 'n';
  char t2 = transpose_m2 ? 't' : 'n';

  long int m=transpose_m1 ? ncol1 : nrow1;
  long int n=transpose_m2 ? nrow2 : ncol2;
  long int k=transpose_m1 ? nrow1 : ncol1;
  long int k2=transpose_m2 ? ncol2 : nrow2;
  long int lda=nrow1;
  long int ldb=nrow2;
  long int ldc=m;

  if (m1==res || m2==res) {
    IGRAPH_ERROR("The input and output matrices must be different",
		 IGRAPH_EINVAL);
  }

  if (k != k2) {
    IGRAPH_ERROR("Invalid matrix sizes for multiplication", IGRAPH_EINVAL);
  }
  
  if (beta != 0.0 && 
      (igraph_matrix_nrow(res) != m ||
       igraph_matrix_ncol(res) != n)) {
    IGRAPH_ERROR("Non-zero `beta' and bad `res' matrix size, possible mistake",
		 IGRAPH_EINVAL);
  }
  
  IGRAPH_CHECK(igraph_matrix_resize(res, m, n));
  
  ret= igraphdgemm_(&t1, &t2, &m, &n, &k, &alpha,
		    &MATRIX(*m1,0,0), &lda, &MATRIX(*m2,0,0), &ldb, 
		    &beta, &MATRIX(*res, 0, 0), &ldc);
  
  if (ret) {
    IGRAPH_ERROR("Could not perform matrix multiplication", IGRAPH_EINVAL);
  }
  
  return ret;
}			

int igraph_matrix_mprod(const igraph_matrix_t *m1,
			const igraph_matrix_t *m2,
			igraph_matrix_t *res) {

  return igraph_matrix_dgemm(m1, m2, res, 1.0, 0.0, 
			     /* transpose_m1= */ 0, 
			     /* transpose_m2= */ 0);
}

int igraph_matrix_crossprod(const igraph_matrix_t *m1,
			    const igraph_matrix_t *m2,
			    igraph_matrix_t *res) {
  
  return igraph_matrix_dgemm(m1, m2, res, 1.0, 0.0, 
			     /* transpose_m1= */ 1, 
			     /* transpose_m2= */ 0);
}

int igraph_matrix_tcrossprod(const igraph_matrix_t *m1,
			     const igraph_matrix_t *m2,
			     igraph_matrix_t *res) {

  return igraph_matrix_dgemm(m1, m2, res, 1.0, 0.0, 
			     /* transpose_m1= */ 0, 
			     /* transpose_m2= */ 1);
}
