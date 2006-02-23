
#include <igraph.h>

int main() {
  
  igraph_t g;
  integer_t result;
  igraph_vector_t edges, res;
  igraph_vs_t vids;
  long i, n;
  
  igraph_tree(&g, 10, 3, IGRAPH_TREE_OUT);
  igraph_rewire(&g, 1000, IGRAPH_REWIRING_SIMPLE);
  
  n=igraph_vcount(&g);
  igraph_vector_init(&res, 0);
  igraph_degree(&g, &res, IGRAPH_VS_ALL(&g), IGRAPH_IN, 0);
  for (i=0; i<n; i++)
    printf("%ld ", (long)igraph_vector_e(&res, i));
  printf("\n");
  
  igraph_degree(&g, &res, IGRAPH_VS_ALL(&g), IGRAPH_OUT, 0);
  for (i=0; i<n; i++)
    printf("%ld ", (long)igraph_vector_e(&res, i));
  printf("\n");
    
  igraph_destroy(&g);
  igraph_vector_destroy(&res);
  
  return 0;
}
