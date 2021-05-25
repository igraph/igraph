#include <igraph.h>

int main() {
  igraph_t graph;
  igraph_vector_t dimvector;
  igraph_vector_t edges;
  igraph_real_t avg_path_len;
  int i;

  igraph_vector_init(&dimvector, 2);
  VECTOR(dimvector)[0]=30;
  VECTOR(dimvector)[1]=30;
  igraph_lattice(&graph, &dimvector, 0, IGRAPH_UNDIRECTED, 0, 1);

  igraph_average_path_length(&graph, &avg_path_len, NULL, IGRAPH_UNDIRECTED, 1);
  printf("Average path length (lattice):            %g\n", (double) avg_path_len);

  igraph_rng_seed(igraph_rng_default(), 42);
  igraph_vector_init(&edges, 20);
  for (i=0; i < igraph_vector_size(&edges); i++) {
    VECTOR(edges)[i] = RNG_INTEGER(0, igraph_vcount(&graph) - 1);
  }

  igraph_add_edges(&graph, &edges, 0);
  igraph_average_path_length(&graph, &avg_path_len, NULL, IGRAPH_UNDIRECTED, 1);
  printf("Average path length (randomized lattice): %g\n", (double) avg_path_len);

  igraph_vector_destroy(&dimvector);
  igraph_vector_destroy(&edges);
  igraph_destroy(&graph);

  return 0;
}
