
#include <igraph.h>
#include <stdio.h>

int main(int argc, char **argv) {

  igraph_t g;
  FILE *ifile;

  /* PAJEK */
  ifile=fopen("LINKS.NET", "r");
  if (ifile==0) {
    return 10;
  }
  igraph_read_graph_pajek(&g, ifile);
  fclose(ifile);
  printf("The graph:\n");
  printf("Vertices: %li\n", (long int) igraph_vcount(&g));
  printf("Edges: %li\n", (long int) igraph_ecount(&g));
  printf("Directed: %i\n", (int) igraph_is_directed(&g));
  igraph_write_graph_edgelist(&g, stdout);
  igraph_destroy(&g);
  return 0;
}
