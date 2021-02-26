#include <igraph.h>

int main() {
  igraph_t graph;
  igraph_vector_t v;
  igraph_vector_t result;
  igraph_real_t edges[] = { 0, 1, 0, 2, 0, 3, 0, 4, 0, 5, 0, 6, 0, 7, 0, 8,
                            0,10, 0,11, 0,12, 0,13, 0,17, 0,19, 0,21, 0,31,
                            1, 2, 1, 3, 1, 7, 1,13, 1,17, 1,19, 1,21, 1,30,
                            2, 3, 2, 7, 2,27, 2,28, 2,32, 2, 9, 2, 8, 2,13,
                            3, 7, 3,12, 3,13, 4, 6, 4,10, 5, 6, 5,10, 5,16,
                            6,16, 8,30, 8,32, 8,33, 9,33,13,33,14,32,14,33,
                           15,32,15,33,18,32,18,33,19,33,20,32,20,33,
                           22,32,22,33,23,25,23,27,23,32,23,33,23,29,
                           24,25,24,27,24,31,25,31,26,29,26,33,27,33,
                           28,31,28,33,29,32,29,33,30,32,30,33,31,32,31,33,
                           32,33
  };

  igraph_vector_view(&v, edges, sizeof(edges) / sizeof(double));
  igraph_create(&graph, &v, 0, IGRAPH_UNDIRECTED);

  igraph_vector_init(&result, 0);

  igraph_degree(&graph, &result, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS);
  printf("Maximum degree is      %10i, vertex %2i.\n",
         (int) igraph_vector_max(&result), (int) igraph_vector_which_max(&result));

  igraph_closeness(&graph, &result, NULL, NULL, igraph_vss_all(), IGRAPH_ALL,
                   /*weights=*/ NULL, /*normalized=*/ 0);
  printf("Maximum closeness is   %10g, vertex %2i.\n",
          (double) igraph_vector_max(&result), (int) igraph_vector_which_max(&result));

  igraph_betweenness(&graph, &result, igraph_vss_all(),
                     IGRAPH_UNDIRECTED, /*weights=*/ NULL);
  printf("Maximum betweenness is %10g, vertex %2i.\n",
          (double) igraph_vector_max(&result), (int) igraph_vector_which_max(&result));

  igraph_vector_destroy(&result);
  igraph_destroy(&graph);

  return 0;
}

