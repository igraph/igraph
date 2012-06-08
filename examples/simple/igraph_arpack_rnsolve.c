/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard st, Cambridge MA, 02139 USA
   
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

typedef struct cb2_data_t {
  igraph_matrix_t *A;
} cb2_data_t;

int cb2(igraph_real_t *to, const igraph_real_t *from, int n, void *extra) {
  cb2_data_t *data=(cb2_data_t*) extra;
  igraph_blas_dgemv_array(/*transpose=*/ 0, /*alpha=*/ 1.0, 
			  data->A, from, /*beta=*/ 0.0, to);
  return 0;
}

#define DIM 10

int main() {
  igraph_matrix_t A;
  igraph_matrix_t values, vectors;
  igraph_arpack_options_t options;
  cb2_data_t data = { &A };  
  int i, j;

  igraph_rng_seed(igraph_rng_default(), 42 * 42);

  igraph_matrix_init(&A, DIM, DIM);

  for (i=0; i<DIM; i++) {
    for (j=0; j<DIM; j++) {
      MATRIX(A, i, j) = igraph_rng_get_integer(igraph_rng_default(), -10, 10);
    }
  }
  
  igraph_arpack_options_init(&options);
  options.n=DIM;
  options.start=0;
  options.nev=4;
  options.ncv=9;
  options.which[0]='L' ; options.which[1]='M';

  igraph_matrix_init(&values, 0, 0);
  igraph_matrix_init(&vectors, options.n, 1);

  igraph_arpack_rnsolve(cb2, /*extra=*/ &data, &options, /*storage=*/ 0, 
			&values, &vectors);

  if (MATRIX(values, 2, 1) > 0) {
    MATRIX(values, 2, 1) = -MATRIX(values, 2, 1);
    MATRIX(values, 3, 1) = -MATRIX(values, 3, 1);    
  }

  igraph_matrix_print(&values);
  printf("---\n");
  igraph_matrix_print(&vectors);
  printf("---\n");

  /* -------------- */

  options.nev=3;
  options.which[0]='L' ; options.which[1]='M';

  igraph_arpack_rnsolve(cb2, /*extra=*/ &data, &options, /*storage=*/ 0, 
			&values, &vectors);

  if (MATRIX(values, 2, 1) > 0) {
    MATRIX(values, 2, 1) = -MATRIX(values, 2, 1);
  }

  igraph_matrix_print(&values);
  printf("---\n");
  igraph_matrix_print(&vectors);
  printf("---\n");

  /* -------------- */

  options.nev=3;
  options.which[0]='S' ; options.which[1]='R';

  igraph_arpack_rnsolve(cb2, /*extra=*/ &data, &options, /*storage=*/ 0, 
			&values, &vectors);

  igraph_matrix_print(&values);
  printf("---\n");
  igraph_matrix_print(&vectors);
  printf("---\n");

  /* -------------- */

  options.nev=3;
  options.which[0]='L' ; options.which[1]='I';

  igraph_arpack_rnsolve(cb2, /*extra=*/ &data, &options, /*storage=*/ 0, 
			&values, &vectors);

  igraph_matrix_print(&values);
  printf("---\n");
  igraph_matrix_print(&vectors);
  printf("---\n");

  /* -------------- */

  igraph_matrix_destroy(&values);
  igraph_matrix_destroy(&vectors);
  igraph_matrix_destroy(&A);
  
  return 0;
}
			
