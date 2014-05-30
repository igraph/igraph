/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2014  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge, MA 02139 USA
   
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
#include <stdio.h>

#define N  10
#define M  20
#define NZ 50

int main() {
  
  int i;
  igraph_sparsemat_t A, A2;
  igraph_vector_t vec;

  igraph_rng_seed(igraph_rng_default(), 42);

  /* Triplet diagonal matrix */

  igraph_vector_init(&vec, N);
  for (i = 0; i < N; i++) { VECTOR(vec)[i] = i; }
  igraph_sparsemat_diag(&A, /*nzmax=*/ N, /*values=*/ &vec,
			/*compress=*/ 0);

  igraph_vector_null(&vec);
  igraph_sparsemat_rowmins(&A, &vec);
  
  for (i = 0; i < N; i++) { if (VECTOR(vec)[i] != i) return 1; }
  
  igraph_vector_null(&vec);
  igraph_sparsemat_colmins(&A, &vec);
  for (i = 0; i < N; i++) { if (VECTOR(vec)[i] != i) return 2; }

  igraph_vector_destroy(&vec);
  igraph_sparsemat_destroy(&A);

  /* Compressed diagonal matrix */
  
  igraph_vector_init(&vec, N);
  for (i = 0; i < N; i++) { VECTOR(vec)[i] = i; }
  igraph_sparsemat_diag(&A, /*nzmax=*/ N, /*values=*/ &vec,
			/*compress=*/ 1);

  igraph_vector_null(&vec);
  igraph_sparsemat_rowmins(&A, &vec);
  for (i = 0; i < N; i++) { if (VECTOR(vec)[i] != i) return 3; }
  
  igraph_vector_null(&vec);
  igraph_sparsemat_colmins(&A, &vec);
  for (i = 0; i < N; i++) { if (VECTOR(vec)[i] != i) return 4; }

  igraph_vector_destroy(&vec);
  igraph_sparsemat_destroy(&A);  
  

  /* Random triplet matrix */
  
  igraph_sparsemat_init(&A, /*rows=*/ N, /*cols=*/ M, /*nzmax=*/ NZ+5);
  for (i = 0; i < NZ; i++) {
    int r = igraph_rng_get_integer(igraph_rng_default(), 0, N - 1);
    int c = igraph_rng_get_integer(igraph_rng_default(), 0, M - 1);
    igraph_real_t x = igraph_rng_get_integer(igraph_rng_default(),
					     -10, 10);
    igraph_sparsemat_entry(&A, r, c, x);
  }

  igraph_vector_init(&vec, 0);
  igraph_sparsemat_colmins(&A, &vec);
  igraph_vector_print(&vec);
  
  igraph_vector_null(&vec);
  igraph_sparsemat_rowmins(&A, &vec);
  igraph_vector_print(&vec);

  /* Random compresssed matrix */
  
  igraph_sparsemat_compress(&A, &A2);
  
  igraph_vector_null(&vec);
  igraph_sparsemat_colmins(&A2, &vec);
  igraph_vector_print(&vec);
  
  igraph_vector_null(&vec);
  igraph_sparsemat_rowmins(&A2, &vec);
  igraph_vector_print(&vec);

  igraph_vector_destroy(&vec);
  igraph_sparsemat_destroy(&A);
  igraph_sparsemat_destroy(&A2);

  /* Matrix with zero rows, triplet */
  
  igraph_sparsemat_init(&A, /*rows=*/ 0, /*cols=*/ M, /*nzmax=*/ NZ);
  
  igraph_vector_init(&vec, 5);
  igraph_sparsemat_rowmins(&A, &vec);
  if (igraph_vector_size(&vec) != 0) { return 5; }

  igraph_vector_null(&vec);
  igraph_sparsemat_colmins(&A, &vec);
  igraph_vector_print(&vec);

  /* Matrix with zero rows, compressed */

  igraph_sparsemat_compress(&A, &A2);
  
  igraph_vector_null(&vec);
  igraph_sparsemat_rowmins(&A, &vec);
  if (igraph_vector_size(&vec) != 0) { return 6; }

  igraph_vector_null(&vec);
  igraph_sparsemat_colmins(&A, &vec);
  igraph_vector_print(&vec);
  
  igraph_vector_destroy(&vec);
  igraph_sparsemat_destroy(&A);
  igraph_sparsemat_destroy(&A2);
  
  /* Matrix with zero columns, triplet */
  
  igraph_sparsemat_init(&A, /*rows=*/ N, /*cols=*/ 0, /*nzmax=*/ NZ);
  
  igraph_vector_init(&vec, 5);
  igraph_sparsemat_colmins(&A, &vec);
  if (igraph_vector_size(&vec) != 0) { return 5; }

  igraph_vector_null(&vec);
  igraph_sparsemat_rowmins(&A, &vec);
  igraph_vector_print(&vec);

  /* Matrix with zero columns, compressed */

  igraph_sparsemat_compress(&A, &A2);
  
  igraph_vector_null(&vec);
  igraph_sparsemat_colmins(&A, &vec);
  if (igraph_vector_size(&vec) != 0) { return 6; }

  igraph_vector_null(&vec);
  igraph_sparsemat_rowmins(&A, &vec);
  igraph_vector_print(&vec);
  
  igraph_vector_destroy(&vec);
  igraph_sparsemat_destroy(&A);
  igraph_sparsemat_destroy(&A2);

  return 0;
}
