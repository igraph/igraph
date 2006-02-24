
#include <igraph.h>

int main() {
  real_t coords_array[][2] =
   {{3, 2}, {5, 1}, {4, 4}, {6, 4}, {4, 3},
    {2, 5}, {1, 3}, {2, 4}, {6, 3}, {9, 2}
   };
    
  igraph_matrix_t coords;
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
  
  /* Testing with coordinate output mode */
  igraph_vector_init(&result, 1);
  if (igraph_convex_hull(&coords, &result, 1))
    return 1;

  for (i=0; i<igraph_vector_size(&result); i++)
    printf("%ld ", (long)VECTOR(result)[i]);
  printf("\n");
  
  igraph_vector_destroy(&result);
  igraph_matrix_destroy(&coords);
  
  return 0;
}
