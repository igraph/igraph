/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
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
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA 
   02110-1301 USA

*/

#include <igraph.h>

int main() {
  igraph_real_t coords_array[][2] =
   {{3, 2}, {5, 1}, {4, 4}, {6, 4}, {4, 3},
    {2, 5}, {1, 3}, {2, 4}, {6, 3}, {9, 2}
   };
    
  igraph_matrix_t coords, resmat;
  igraph_vector_t result;
  long i;
  
  igraph_matrix_init(&coords, 10, 2);
  for (i=0; i<20; i++) MATRIX(coords, i/2, i%2) = coords_array[i/2][i%2];
  
  /* Testing with index output mode */
  igraph_vector_init(&result, 1);
  if (igraph_convex_hull(&coords, &result, 0))
    return 1;

  for (i=0; i<igraph_vector_size(&result); i++)
    printf("%ld ", (long)VECTOR(result)[i]);
  printf("\n");
  igraph_vector_destroy(&result);

  /* Testing with coordinate output mode */
  igraph_matrix_init(&resmat, 0, 0);
  if (igraph_convex_hull(&coords, 0, &resmat))
    return 1;

  for (i=0; i<igraph_matrix_nrow(&resmat); i++)
    printf("%ld %ld ", (long)MATRIX(resmat, i, 0), (long)MATRIX(resmat, i, 1));
  printf("\n");
  
  igraph_matrix_destroy(&resmat);
  igraph_matrix_destroy(&coords);
  
  return 0;
}
