/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2009  Gabor Csardi <csardi@rmki.kfki.hu>
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
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
   02110-1301 USA

*/

#include <igraph.h>
#include <igraph_sparsemat.h>

#define RNG_INTEGER(l, h) ((long int)((rand())/((double)RAND_MAX+1)*((h)-(l)+1)+(l)))
#define EPS 1e-13

int main() {

  igraph_sparsemat_t A, B;
  igraph_matrix_t vectors, values2;
  igraph_vector_t values;
  long int i;
  igraph_arpack_options_t options;
  igraph_real_t min, max;
  igraph_t g1, g2, g3;

  /* Identity matrix */
#define DIM 10
  igraph_sparsemat_init(&A, DIM, DIM, DIM);
  for (i=0; i<DIM; i++) {
    igraph_sparsemat_entry(&A, i, i, 1.0);
  }
  igraph_sparsemat_compress(&A, &B);
  igraph_sparsemat_destroy(&A);

  igraph_vector_init(&values, 0);
  igraph_arpack_options_init(&options);
  
  igraph_sparsemat_arpack_rssolve(&B, &options, /*storage=*/ 0, 
				  &values, /*vectors=*/ 0);

  if (VECTOR(values)[0] != 1.0) { return 1; }

  igraph_vector_destroy(&values);
  igraph_sparsemat_destroy(&B);  

#undef DIM

  /* Diagonal matrix */
#define DIM 10
  igraph_sparsemat_init(&A, DIM, DIM, DIM);
  for (i=0; i<DIM; i++) {
    igraph_sparsemat_entry(&A, i, i, i+1.0);
  }
  igraph_sparsemat_compress(&A, &B);
  igraph_sparsemat_destroy(&A);
  
  igraph_vector_init(&values, 0);
  igraph_matrix_init(&vectors, 0, 0);
  
  igraph_sparsemat_arpack_rssolve(&B, &options, /*storage=*/ 0,
				  &values, /*vectors=*/ &vectors);

  if ( fabs(VECTOR(values)[0] - DIM) > EPS ) {
    printf("VECTOR(values)[0] numerical precision is only %g, should be %g",
			fabs((double)VECTOR(values)[0]-DIM), EPS);
	return 2;
  }

  if ( fabs(fabs(MATRIX(vectors, DIM-1, 0)) - 1.0) > EPS) { return 3; }
  MATRIX(vectors, DIM-1, 0) = 0.0;
  igraph_matrix_minmax(&vectors, &min, &max);
  if (fabs(min) > EPS) { return 3; }
  if (fabs(max) > EPS) { return 3; }
  
  igraph_vector_destroy(&values);
  igraph_matrix_destroy(&vectors);
  igraph_sparsemat_destroy(&B);
#undef DIM

  /* A tree, plus a ring */
#define DIM 10
  igraph_tree(&g1, DIM, /*children=*/ 2, IGRAPH_TREE_UNDIRECTED);
  igraph_ring(&g2, DIM, IGRAPH_UNDIRECTED, /*mutual=*/ 0, /*circular=*/ 1);
  igraph_union(&g3, &g1, &g2);
  igraph_destroy(&g1);
  igraph_destroy(&g2);

  igraph_get_sparsemat(&g3, &A);
  igraph_destroy(&g3);
  igraph_sparsemat_compress(&A, &B);
  igraph_sparsemat_destroy(&A);
  
  igraph_vector_init(&values, 0);
  igraph_matrix_init(&vectors, 0, 0);
  
  igraph_sparsemat_arpack_rssolve(&B, &options, /*storage=*/ 0,
				  &values, &vectors);

  if (MATRIX(vectors, 0, 0) < 0.0) { 
    igraph_matrix_scale(&vectors, -1.0);
  }
  
  igraph_vector_print(&values);
  igraph_matrix_print(&vectors);
  
  igraph_vector_destroy(&values);
  igraph_matrix_destroy(&vectors);
  igraph_sparsemat_destroy(&B);
#undef DIM

  printf("--\n");

  /* A directed tree and a directed, mutual ring */
#define DIM 10
  igraph_tree(&g1, DIM, /*children=*/ 2, IGRAPH_TREE_OUT);
  igraph_ring(&g2, DIM, IGRAPH_DIRECTED, /*mutual=*/ 1, /*circular=*/ 1);
  igraph_union(&g3, &g1, &g2);
  igraph_destroy(&g1);
  igraph_destroy(&g2);

  igraph_get_sparsemat(&g3, &A);
  igraph_destroy(&g3);
  igraph_sparsemat_compress(&A, &B);
  igraph_sparsemat_destroy(&A);
  
  igraph_matrix_init(&values2, 0, 0);
  igraph_matrix_init(&vectors, 0, 0);
  
  igraph_sparsemat_arpack_rnsolve(&B, &options, /*storage=*/ 0,
				  &values2, &vectors);

  if (MATRIX(vectors, 0, 0) < 0.0) { 
    igraph_matrix_scale(&vectors, -1.0);
  }
  
  igraph_matrix_print(&values2);
  igraph_matrix_print(&vectors);
  
  igraph_matrix_destroy(&values2);
  igraph_matrix_destroy(&vectors);
  igraph_sparsemat_destroy(&B);
#undef DIM


  return 0;
}
