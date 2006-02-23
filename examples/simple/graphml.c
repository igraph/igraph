#include <igraph.h>
#include <stdio.h>

int main(int argc, char **argv) {
  igraph_t g;
  igraph_error_handler_t* oldhandler;
  int result;
  FILE *ifile, *ofile;

  /* GraphML */
  ifile=fopen("test.gxl", "r");
  if (ifile==0) {
    return 10;
  }
  
  oldhandler=igraph_set_error_handler(igraph_error_handler_ignore);
  if (result=igraph_read_graph_graphml(&g, ifile, 1, 0)) {
    // maybe it is simply disabled at compile-time
    if (result == IGRAPH_UNIMPLEMENTED) return 77;
    return 1;
  }
  igraph_set_error_handler(oldhandler);
    
  fclose(ifile);

  /* Write it back */
  ofile=fopen("test2.gxl", "w");
  /* If we can't create the test file, just skip the test */
  if (ofile) {
    if (result=igraph_write_graph_graphml(&g, ofile)) {
      return 1;
    }
    fclose(ofile);
    unlink("test2.gxl");
  }
  
  printf("The directed graph:\n");
  printf("Vertices: %li\n", (long int) igraph_vcount(&g));
  printf("Edges: %li\n", (long int) igraph_ecount(&g));
  printf("Directed: %i\n", (int) igraph_is_directed(&g));
  igraph_write_graph_edgelist(&g, stdout);
  igraph_destroy(&g);
 
  /* The same with undirected graph */
  ifile=fopen("test.gxl", "r");
  if (result=igraph_read_graph_graphml(&g, ifile, 0, 0)) {
    return 1;
  }
  fclose(ifile);
  printf("The undirected graph:\n");
  printf("Vertices: %li\n", (long int) igraph_vcount(&g));
  printf("Edges: %li\n", (long int) igraph_ecount(&g));
  printf("Directed: %i\n", (int) igraph_is_directed(&g));
  igraph_write_graph_edgelist(&g, stdout);
  igraph_destroy(&g);
  
  return 0;
}
