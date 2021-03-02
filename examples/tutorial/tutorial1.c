#include <igraph.h>

int main() {
  igraph_real_t diameter;
  igraph_t graph;

  igraph_rng_seed(igraph_rng_default(), 42);

  igraph_erdos_renyi_game(&graph, IGRAPH_ERDOS_RENYI_GNM, 1000, 3000,
                          IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);

  igraph_diameter(&graph, &diameter, 0, 0, 0, IGRAPH_UNDIRECTED, 1);
  printf("Diameter of a random graph with average degree %g: %g\n",
          2.0 * igraph_ecount(&graph) / igraph_vcount(&graph), 
          (double) diameter);

  igraph_destroy(&graph);

  return 0;
}

