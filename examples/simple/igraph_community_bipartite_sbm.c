#include "igraph.h"

/* This only tests a case that really *should* work.
   TODO: add tests for any edge cases and tricky situations. */
void test_easy() {
  igraph_t g;
  igraph_real_t score;
  igraph_vector_t membership;
  igraph_vector_init(&membership, 0);

  /* The "Davis southern club women" dataset from 
  http://vlado.fmf.uni-lj.si/pub/networks/data/ucinet/ucidata.htm 

  Originally published here:
  Allison Davis, Burleigh B. Gardner, and Mary R. Gardner. 
    Deep South; a Social Anthropological Study of Caste and Class. 
    The University of Chicago Press, Chicago, 1941.

  This is a bipartite graph representing 18 women and the 14
  events that they might have attended. Women are indexed 0-17,
  events 18-31. There is an edge between woman w and a event e if
  w attended e.*/
  printf("Testing Southern Women dataset...\n");
  igraph_small(&g, 0, IGRAPH_UNDIRECTED, 
               0,18, 0,19, 0,20, 0,21, 0,22, 0,23, 0,25, 0,26,
               1,18, 1,19, 1,20, 1,22, 1,23, 1,24, 1,25,
               2,19, 2,20, 2,21, 2,22, 2,23, 2,24, 2,25, 2,26,
               3,18, 3,20, 3,21, 3,22, 3,23, 3,24, 3,25,
               4,20, 4,21, 4,22, 4,24,
               5,20, 5,22, 5,23, 5,25,
               6,22, 6,23, 6,24, 6,25,
               7,23, 7,25, 7,26,
               8,22, 8,24, 8,25, 8,26,
               9,24, 9,25, 9,26, 9,29,
               10,25, 10,26, 10,27, 10,29,
               11,25, 11,26, 11,27, 11,29, 11,30, 11,31,
               12,24, 12,25, 12,26, 12,27, 12,29, 12,30, 12,31,
               13,23, 13,24, 13,26, 13,27, 13,28, 13,29, 13,30, 13,31,
               14,24, 14,25, 14,27, 14,28, 14,29, 14,30, 14,31,
               15,25, 15,26, 15,27, 15,29,
               16,26, 16,28,
               17,26, 17,28,
               -1);

  printf("%d vertices\n", igraph_vcount(&g));
  printf("%d edges\n", igraph_ecount(&g));
  //igraph_rng_seed(igraph_rng_default(), 42);
  igraph_community_bipartite_sbm(&g, &membership, 2, 3, 0, &score);
  printf("score: %f\n", score);
  printf("community assignments: ");
  igraph_vector_print(&membership);
  igraph_vector_destroy(&membership);
  igraph_destroy(&g);
}

int main() {
  igraph_rng_seed(igraph_rng_default(), 42);
  test_easy();
  return 0;
}
