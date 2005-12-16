
#include <igraph.h>

int main() {
  
  igraph_t g;
  integer_t result;
  
  igraph_barabasi_game(&g, 30, 30, 0, 0, IGRAPH_DIRECTED);
  igraph_average_path_length(&g, &result, IGRAPH_UNDIRECTED, 1);
  
/*   printf("Length of the average shortest paths: %f\n", (float) result); */
  
  igraph_destroy(&g);
  return 0;
}
