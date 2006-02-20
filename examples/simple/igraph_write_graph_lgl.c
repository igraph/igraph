#include <igraph.h>
#include <stdio.h>

int main(int argc, char **argv) {

  igraph_t g;
  igraph_error_handler_t *oldhandler;
  FILE *ofile;
  int ret;

  /* Testing error handling */
  igraph_barabasi_game(&g, 10, 1, 0, 0, IGRAPH_DIRECTED);
  oldhandler=igraph_set_error_handler(igraph_error_handler_ignore);
  ofile=fopen("test.txt", "w");
  ret=igraph_write_graph_lgl(&g, ofile, "names", "weights", 1);
  if (ret != IGRAPH_EINVAL) {
    return 1;
  }
  fclose(ofile);
  igraph_destroy(&g);
  igraph_set_error_handler(oldhandler);

  return 0;
}
