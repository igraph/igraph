
#include <igraph.h>

void print_vector(igraph_vector_t *v) {
  long int i;
  real_t sum=igraph_vector_sum(v);
  for (i=0; i<igraph_vector_size(v); i++) {
    printf("%2.2f ", VECTOR(*v)[i]/sum);
  }
  printf("\n");
}

int main() {

  igraph_t g;
  igraph_vector_t hist;
  igraph_vector_t cp;

  igraph_vector_init_real(&cp, 8, 0.0, 0.0, 0.0, 0.0);

  igraph_ring(&g, 1000, IGRAPH_DIRECTED, 0, 1);
  igraph_vector_init(&hist, 0);
  igraph_motifs_randesu(&g, &hist, 3, &cp);
  print_vector(&hist);
  igraph_destroy(&g);
  igraph_vector_destroy(&hist);
  igraph_vector_destroy(&cp);
  
  return 0;
}
