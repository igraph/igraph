
#include <igraph.h>

int main() {
  
  igraph_t g1, g2;
  
  /* complementer of the empty graph */
  igraph_empty(&g1, 5, IGRAPH_DIRECTED);
  igraph_complementer(&g2, &g1, IGRAPH_LOOPS);
  igraph_write_graph_edgelist(&g2, stdout);
  igraph_destroy(&g1);
  igraph_destroy(&g2);

  printf("---\n");

  /* the same without loops */
  igraph_empty(&g1, 5, IGRAPH_DIRECTED);
  igraph_complementer(&g2, &g1, IGRAPH_NO_LOOPS);
  igraph_write_graph_edgelist(&g2, stdout);
  igraph_destroy(&g1);
  igraph_destroy(&g2);

  printf("---\n");
  
  /* complementer of the full graph */
  igraph_full(&g1, 5, IGRAPH_DIRECTED, IGRAPH_LOOPS);
  igraph_complementer(&g2, &g1, IGRAPH_LOOPS);
  if (igraph_ecount(&g2) != 0) {
    return 1;
  }
  igraph_destroy(&g1);
  igraph_destroy(&g2);

  printf("---\n");

  /* complementer of the full graph, results loops only */
  igraph_full(&g1, 5, IGRAPH_DIRECTED, IGRAPH_NO_LOOPS);
  igraph_complementer(&g2, &g1, IGRAPH_LOOPS);
  igraph_write_graph_edgelist(&g2, stdout);
  igraph_destroy(&g1);
  igraph_destroy(&g2);

  printf("---\n");

  /**************
   * undirected *
   *************/

  /* complementer of the empty graph */
  igraph_empty(&g1, 5, IGRAPH_UNDIRECTED);
  igraph_complementer(&g2, &g1, IGRAPH_LOOPS);
  igraph_write_graph_edgelist(&g2, stdout);
  igraph_destroy(&g1);
  igraph_destroy(&g2);

  printf("---\n");

  /* the same without loops */
  igraph_empty(&g1, 5, IGRAPH_UNDIRECTED);
  igraph_complementer(&g2, &g1, IGRAPH_NO_LOOPS);
  igraph_write_graph_edgelist(&g2, stdout);
  igraph_destroy(&g1);
  igraph_destroy(&g2);

  printf("---\n");
  
  /* complementer of the full graph */
  igraph_full(&g1, 5, IGRAPH_UNDIRECTED, IGRAPH_LOOPS);
  igraph_complementer(&g2, &g1, IGRAPH_LOOPS);
  if (igraph_ecount(&g2) != 0) {
    return 1;
  }
  igraph_destroy(&g1);
  igraph_destroy(&g2);

  printf("---\n");

  /* complementer of the full graph, results loops only */
  igraph_full(&g1, 5, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
  igraph_complementer(&g2, &g1, IGRAPH_LOOPS);
  igraph_write_graph_edgelist(&g2, stdout);
  igraph_destroy(&g1);
  igraph_destroy(&g2);
  
  return 0;
}
