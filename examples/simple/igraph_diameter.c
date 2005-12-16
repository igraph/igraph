
#include <igraph.h>

int main() {
  
  igraph_t g;
  integer_t result;
  
  igraph_barabasi_game(&g, 30, 30, 0, 0, IGRAPH_DIRECTED);
  igraph_diameter(&g, &result, IGRAPH_UNDIRECTED, 1);
  
/*   printf("Diameter: %li\n", (long int) result); */
  
  igraph_destroy(&g);
  return 0;
}
