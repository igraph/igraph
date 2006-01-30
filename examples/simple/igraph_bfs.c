
#include <igraph.h>

void vector_print(igraph_vector_t *v) {
  long int i;
  for (i=0; i<igraph_vector_size(v); i++) {
    printf(" %li", (long int) VECTOR(*v)[i]);
  }
  printf("\n");
}

int main() {

  igraph_t g;
  igraph_vector_t vids, layers, parents;

  igraph_ring(&g, 10, IGRAPH_UNDIRECTED, 0, 0);
  igraph_vector_init(&vids, 0);
  igraph_vector_init(&layers, 0);
  igraph_vector_init(&parents, 0);
  igraph_bfs(&g, 0, IGRAPH_ALL, &vids, &layers, &parents);
  vector_print(&vids);
  vector_print(&layers);
  vector_print(&parents);
  igraph_destroy(&g);  

  igraph_tree(&g, 20, 2, IGRAPH_TREE_UNDIRECTED);
  igraph_bfs(&g, 0, IGRAPH_ALL, &vids, &layers, &parents);
  vector_print(&vids);
  vector_print(&layers);
  vector_print(&parents);
  igraph_destroy(&g);  
  
  igraph_vector_destroy(&vids);
  igraph_vector_destroy(&layers);
  igraph_vector_destroy(&parents);

  return 0;
}
