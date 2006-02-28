
#include <igraph.h>

int main() {

  igraph_t g1, g2, res;
  igraph_vector_t v;

  /* composition with the empty graph */
  igraph_empty(&g1, 5, IGRAPH_DIRECTED);
  igraph_full(&g2, 5, IGRAPH_DIRECTED, IGRAPH_NO_LOOPS);
  igraph_compose(&res, &g1, &g2);
  if (igraph_ecount(&res) != 0) { 
    return 1;
  }
  igraph_destroy(&res);
  igraph_compose(&res, &g2, &g1);
  if (igraph_ecount(&res) != 0) { 
    return 2;
  }
  igraph_destroy(&res);
  igraph_destroy(&g1);
  igraph_destroy(&g2);

  /* same but undirected */
  igraph_empty(&g1, 5, IGRAPH_UNDIRECTED);
  igraph_full(&g2, 5, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
  igraph_compose(&res, &g1, &g2);
  if (igraph_ecount(&res) != 0) { 
    return 1;
  }
  igraph_destroy(&res);
  igraph_compose(&res, &g2, &g1);
  if (igraph_ecount(&res) != 0) { 
    return 2;
  }
  igraph_destroy(&res);
  igraph_destroy(&g1);
  igraph_destroy(&g2);

  /* proper directed graph */
  igraph_vector_init_int_end(&v, -1, 0,1, 1,2, 5,6, -1);
  igraph_create(&g1, &v, 0, IGRAPH_DIRECTED);
  igraph_vector_destroy(&v);

  igraph_vector_init_int_end(&v, -1, 0,1, 2,4, 5,6, -1);
  igraph_create(&g2, &v, 0, IGRAPH_DIRECTED);
  igraph_vector_destroy(&v);
  
  igraph_compose(&res, &g1, &g2);
  igraph_write_graph_edgelist(&res, stdout);
  igraph_destroy(&res);
  igraph_destroy(&g1);
  igraph_destroy(&g2);

  /* undirected graph */
  igraph_vector_init_int_end(&v, -1, 0,1, 1,2, 5,6, -1);
  igraph_create(&g1, &v, 0, IGRAPH_UNDIRECTED);
  igraph_vector_destroy(&v);

  igraph_vector_init_int_end(&v, -1, 0,1, 0,4, 5,6, -1);
  igraph_create(&g2, &v, 0, IGRAPH_UNDIRECTED);
  igraph_vector_destroy(&v);
  
  igraph_compose(&res, &g1, &g2);
  igraph_write_graph_edgelist(&res, stdout);
  igraph_destroy(&res);
  igraph_destroy(&g1);
  igraph_destroy(&g2);

  return 0;
}
