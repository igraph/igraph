
#include <igraph.h>

int main() {

  igraph_t g;
  int ret;

  igraph_atlas(&g, 45);
  igraph_write_graph_edgelist(&g, stdout);
  printf("\n");
  igraph_destroy(&g);

  igraph_atlas(&g, 0);
  igraph_write_graph_edgelist(&g, stdout);
  printf("\n");
  igraph_destroy(&g);

  igraph_atlas(&g, 1252);
  igraph_write_graph_edgelist(&g, stdout);
  printf("\n");
  igraph_destroy(&g);

  igraph_set_error_handler(igraph_error_handler_ignore);
  ret=igraph_atlas(&g, -1);
  if (ret != IGRAPH_EINVAL) {
    return 1;
  }

  ret=igraph_atlas(&g, 1253);
  if (ret != IGRAPH_EINVAL) {
    return 2;
  }

  return 0;
}
