
#include <igraph.h>

int main() {

  igraph_t g;
  int i, ret;
  
  /* G(n,p) */
  igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNP, 10, 0.0, 
			  IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
  if (igraph_ecount(&g) != 0) {
    return 1;
  }
  if (igraph_is_directed(&g)) {
    return 2;
  }
  igraph_destroy(&g);
  
  igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNP, 10, 1.0,
			  IGRAPH_DIRECTED, IGRAPH_NO_LOOPS);
  if (igraph_ecount(&g) != 10*9) {
    return 3;
  }
  if (!igraph_is_directed(&g)) {
    return 4;
  }
  igraph_destroy(&g);

  igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNM, 10, 0.5,
			  IGRAPH_DIRECTED, IGRAPH_NO_LOOPS);
  igraph_destroy(&g);

  /* G(n,m) */
  

  return 0;
}
