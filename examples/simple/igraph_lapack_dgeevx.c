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

#include <igraph.h>
#include <time.h>
#include <stdio.h>

#define DIM 10

igraph_bool_t check_ev(const igraph_matrix_t *A, 
		       const igraph_vector_t *values_real,
		       const igraph_vector_t *values_imag,
		       const igraph_matrix_t *vectors_left, 
		       const igraph_matrix_t *vectors_right, 
		       igraph_real_t tol) {

  int n=igraph_matrix_nrow(A);
  int i, j;
  
  if (igraph_matrix_ncol(A)             != n) { return 1; }
  if (igraph_vector_size(values_real)   != n) { return 1; }
  if (igraph_vector_size(values_imag)   != n) { return 1; }
  if (igraph_matrix_nrow(vectors_left)  != n) { return 1; }
  if (igraph_matrix_ncol(vectors_left)  != n) { return 1; }
  if (igraph_matrix_nrow(vectors_right) != n) { return 1; }
  if (igraph_matrix_ncol(vectors_right) != n) { return 1; }

  /* TODO, right now we just compare to the output that I 
     checked in R */
  
  return 0;
}

int main() {
  
  igraph_matrix_t A;
  igraph_matrix_t vectors_left, vectors_right;
  igraph_vector_t values_real, values_imag;
  int i, j;
  int info=1;
  int ilo, ihi;
  igraph_real_t abnrm;
  
  igraph_rng_seed(&igraph_rng_default, 42);
  
  igraph_matrix_init(&A, DIM, DIM);
  igraph_matrix_init(&vectors_left, 0, 0);
  igraph_matrix_init(&vectors_right, 0, 0);
  igraph_vector_init(&values_real, 0);
  igraph_vector_init(&values_imag, 0);

  for (i=0; i<DIM; i++) {
    for (j=0; j<DIM; j++) {
      MATRIX(A, i, j) = igraph_rng_get_integer(&igraph_rng_default, 1, 10);
    }
  }
  
  igraph_lapack_dgeevx(IGRAPH_LAPACK_DGEEVX_BALANCE_BOTH,
		       &A, &values_real, &values_imag, 
		       &vectors_left, &vectors_right, &ilo, &ihi,
		       /*scale=*/ 0, &abnrm, /*rconde=*/ 0, 
		       /*rcondv=*/ 0, &info);

  if (check_ev(&A, &values_real, &values_imag, 
	       &vectors_left, &vectors_right, /*tol=*/ 1e-8)) {
    return 1;
  }
  
  igraph_matrix_print(&A);
  igraph_vector_print(&values_real);
  igraph_vector_print(&values_imag);
  igraph_matrix_print(&vectors_left);
  igraph_matrix_print(&vectors_right);
  
  igraph_vector_destroy(&values_imag);
  igraph_vector_destroy(&values_real);
  igraph_matrix_destroy(&vectors_right);
  igraph_matrix_destroy(&vectors_left);
  igraph_matrix_destroy(&A);

  return 0;
}
