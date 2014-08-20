#include "igraph.h"

// TODO this is a placeholder. Since this algorithm is stochastic, need to figure out how to 
// make consistent output.
int main(){
  igraph_t graph;
  igraph_bipartite_game(&graph, NULL, IGRAPH_ERDOS_RENYI_GNM, 100, 100, 1, 500, 0, IGRAPH_ALL);
  printf("%d\n", igraph_vcount(&graph));
  igraph_vector_t membership;
  igraph_vector_init(&membership, 0);
  igraph_community_bipartite_sbm(&graph, &membership, 3, 3, 10);
  return 0;
}
