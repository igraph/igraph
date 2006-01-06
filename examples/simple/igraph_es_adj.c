
#include <igraph.h>

void vector_print(const vector_t *v) {
  long int i;
  for (i=0; i<vector_size(v); i++) {
    printf("%li ", (long int)VECTOR(*v)[i]);
  }
  printf("\n");
}

int main() {
  
  igraph_t g;
  const vector_t v;
  real_t edges1[] = { 0,1, 1,2, 2,2, 2,3, 2,4, 3,4 };
  vector_t was;
  igraph_es_t it;
  long int i;
  
  vector_view(&v, edges1, sizeof(edges1)/sizeof(real_t));
  vector_init(&was, 0);

  /******************************************/
  /* Directed graph                         */
  /******************************************/
 
  igraph_create(&g, &v, 0, IGRAPH_DIRECTED);
  
  /* Simple test, all neighbors */
  for (i=0; i<=vector_max(&v); i++) {
    vector_clear(&was);
    igraph_es_adj(&g, &it, i, IGRAPH_ALL);
    while (!igraph_es_end(&g, &it)) {
      vector_push_back(&was, igraph_es_adj_vertex(&g, &it));
      igraph_es_next(&g, &it);
    }
    igraph_es_destroy(&it);
    vector_sort(&was);
    vector_print(&was);
  }

  /* Simple test, outgoing neighbors */
  for (i=0; i<=vector_max(&v); i++) {
    vector_clear(&was);
    igraph_es_adj(&g, &it, i, IGRAPH_OUT);
    while (!igraph_es_end(&g, &it)) {
      vector_push_back(&was, igraph_es_adj_vertex(&g, &it));
      igraph_es_next(&g, &it);
    }
    igraph_es_destroy(&it);
    vector_sort(&was);
    vector_print(&was);
  }


  /* Simple test, incoming neighbors */
  for (i=0; i<=vector_max(&v); i++) {
    vector_clear(&was);
    igraph_es_adj(&g, &it, i, IGRAPH_IN);
    while (!igraph_es_end(&g, &it)) {
      vector_push_back(&was, igraph_es_adj_vertex(&g, &it));
      igraph_es_next(&g, &it);
    }
    igraph_es_destroy(&it);
    vector_sort(&was);
    vector_print(&was);
  }
		       
  igraph_destroy(&g);

  /******************************************/
  /* Undirected graph                       */
  /******************************************/

  igraph_create(&g, &v, 0, IGRAPH_UNDIRECTED);

  /* Simple test, all neighbors */
  for (i=0; i<=vector_max(&v); i++) {
    vector_clear(&was);
    igraph_es_adj(&g, &it, i, IGRAPH_ALL);
    while (!igraph_es_end(&g, &it)) {
      vector_push_back(&was, igraph_es_adj_vertex(&g, &it));
      igraph_es_next(&g, &it);
    }
    igraph_es_destroy(&it);
    vector_sort(&was);
    vector_print(&was);
  }

  igraph_destroy(&g);
  vector_destroy(&was);
  
  return 0;
}
