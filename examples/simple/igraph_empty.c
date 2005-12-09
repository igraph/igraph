
#include <igraph.h>

int main() {
  
  igraph_t g;
  int ret;

  /* empty directed graph, zero vertices */
  igraph_empty(&g, 0, 1);
  if (igraph_vcount(&g) != 0) {
    return 1;
  }
  if (igraph_ecount(&g) != 0) {
    return 2;
  }
  igraph_destroy(&g);
  
  /* empty undirected graph, zero vertices */
  igraph_empty(&g, 0, 0);
  if (igraph_vcount(&g) != 0) {
    return 3;
  }
  if (igraph_ecount(&g) != 0) {
    return 4;
  }
  igraph_destroy(&g);

  /* empty directed graph, 20 vertices */
  igraph_empty(&g, 20, 1);
  if (igraph_vcount(&g) != 20) {
    return 5;
  }
  if (igraph_ecount(&g) != 0) {
    return 6;
  }
  igraph_destroy(&g);
  
  /* empty undirected graph, 30 vertices */
  igraph_empty(&g, 30, 0);
  if (igraph_vcount(&g) != 30) {
    return 7;
  }
  if (igraph_ecount(&g) != 0) {
    return 8;
  }
  igraph_destroy(&g);

  /* error: negative number of vertices */
  igraph_set_error_handler(igraph_error_handler_ignore);
  ret=igraph_empty(&g, -1, 0);
  if (ret != IGRAPH_EINVAL) {
    return 9;
  }

  return 0;
}
