
#include <igraph.h>

int main() {
  
  igraph_t g;
  igraph_vector_t edges;
  igraph_vector_t vids;
  int class;
  
  igraph_vector_init_int_end(&edges, -1, 
			     0,1, 1,3, 1,4, 1,6, 3,1,
			     4,1, 4,2, 6,4, 6,5, 7,8,
			     8,7, 7,9, 9,7, 8,9, 9,8,
			     -1);
  igraph_create(&g, &edges, 0, IGRAPH_DIRECTED);
  igraph_vector_destroy(&edges);
  
  igraph_vector_init_int_end(&vids, -1, 1,4,6, -1);
  igraph_isoclass_subgraph(&g, &vids, &class);
  printf("class: %i\n", class);
  igraph_vector_destroy(&vids);

  igraph_vector_init_int_end(&vids, -1, 0,1,3, -1);
  igraph_isoclass_subgraph(&g, &vids, &class);
  printf("class: %i\n", class);
  igraph_vector_destroy(&vids);

  igraph_vector_init_int_end(&vids, -1, 7,8,9, -1);
  igraph_isoclass_subgraph(&g, &vids, &class);
  printf("class: %i\n", class);
  igraph_vector_destroy(&vids);

  igraph_vector_init_int_end(&vids, -1, 0,2,5, -1);
  igraph_isoclass_subgraph(&g, &vids, &class);
  printf("class: %i\n", class);
  igraph_vector_destroy(&vids);
  
  igraph_destroy(&g);
  return 0;
}
