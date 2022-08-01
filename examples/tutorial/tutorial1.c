#include <igraph.h>

int main() {
  igraph_integer_t num_vertices = 1000;
  igraph_integer_t num_edges = 1000;
  igraph_real_t diameter;
  igraph_t graph;

  igraph_rng_seed(igraph_rng_default(), 42);

  igraph_erdos_renyi_game(&graph, IGRAPH_ERDOS_RENYI_GNM,
                          num_vertices, num_edges,
                          IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);

  igraph_diameter(&graph, &diameter, NULL, NULL, NULL, NULL, IGRAPH_UNDIRECTED, /* unconn= */ true);
  printf("Diameter of a random graph with average degree %g: %g\n",
          2.0 * igraph_ecount(&graph) / igraph_vcount(&graph),
          (double) diameter);

  igraph_destroy(&graph);

  return 0;
}
