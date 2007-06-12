
#include <igraph.h>

#include <stdlib.h>

int print_vector(igraph_vector_t *v) {
  long int i, l=igraph_vector_size(v);
  for (i=0; i<l; i++) {
    printf(" %li", (long int) VECTOR(*v)[i]);
  }
  printf("\n");
}

int main() {

  igraph_t g;
  igraph_vector_ptr_t vecs;
  long int i;
  igraph_vs_t vs;

  igraph_ring(&g, 10, IGRAPH_DIRECTED, 0, 1);
  
  igraph_vector_ptr_init(&vecs, 5);
  for (i=0; i<igraph_vector_ptr_size(&vecs); i++) {
    VECTOR(vecs)[i] = calloc(1, sizeof(igraph_vector_t));
    igraph_vector_init(VECTOR(vecs)[i], 0);
  }
  igraph_vs_vector_small(&vs, 1, 3, 5, 2, 1,  -1);
  
  igraph_get_shortest_paths(&g, &vecs, 0, vs, IGRAPH_OUT);
  
  for (i=0; i<igraph_vector_ptr_size(&vecs); i++) {
    print_vector(VECTOR(vecs)[i]);
    igraph_vector_destroy(VECTOR(vecs)[i]);
    free(VECTOR(vecs)[i]);
  }
  igraph_vector_ptr_destroy(&vecs);
  igraph_vs_destroy(&vs);
  igraph_destroy(&g);
  
  return 0;
}
