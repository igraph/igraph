#include <igraph.h>

int main(void) {
  igraph_t graph;
  igraph_vector_int_t dimvector;
  igraph_vector_int_t edges;
  igraph_vector_bool_t periodic;
  igraph_real_t avg_path_len;

  igraph_vector_int_init(&dimvector, 2);
  VECTOR(dimvector)[0]=30;
  VECTOR(dimvector)[1]=30;

  igraph_vector_bool_init(&periodic, 2);
  igraph_vector_bool_fill(&periodic, true);
  igraph_square_lattice(&graph, &dimvector, 0, IGRAPH_UNDIRECTED, /* mutual= */ false, &periodic);

  igraph_average_path_length(&graph, &avg_path_len, NULL, IGRAPH_UNDIRECTED, /* unconn= */ true);
  printf("Average path length (lattice):            %g\n", (double) avg_path_len);

  igraph_rng_seed(igraph_rng_default(), 42); /* seed RNG before first use */
  igraph_vector_int_init(&edges, 20);
  for (igraph_integer_t i=0; i < igraph_vector_int_size(&edges); i++) {
    VECTOR(edges)[i] = RNG_INTEGER(0, igraph_vcount(&graph) - 1);
  }

  igraph_add_edges(&graph, &edges, NULL);
  igraph_average_path_length(&graph, &avg_path_len, NULL, IGRAPH_UNDIRECTED, /* unconn= */ true);
  printf("Average path length (randomized lattice): %g\n", (double) avg_path_len);

  igraph_vector_bool_destroy(&periodic);
  igraph_vector_int_destroy(&dimvector);
  igraph_vector_int_destroy(&edges);
  igraph_destroy(&graph);

  return 0;
}
