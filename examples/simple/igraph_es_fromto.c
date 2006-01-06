
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
  vector_t from, to;  
  igraph_es_t it;
  long int i;

  vector_view(&v, edges1, sizeof(edges1)/sizeof(real_t));

  /******************************************/
  /* Directed graph                         */
  /******************************************/
  
  igraph_create(&g, &v, 0, IGRAPH_DIRECTED);
  
  /* {0,1} -> {2,3}, result should be { 1->2 } */
  vector_init(&from, 2); VECTOR(from)[0]=0; VECTOR(from)[1]=1;
  vector_init(&to, 2);   VECTOR(to)  [0]=2; VECTOR(to)  [1]=3;
  igraph_es_fromto(&g, &it, IGRAPH_VS_VECTOR(&g, &from), 
		   IGRAPH_VS_VECTOR(&g, &to), IGRAPH_DIRECTED);
  vector_clear(&from); vector_clear(&to);
  while (!igraph_es_end(&g, &it)) {
    vector_push_back(&from, igraph_es_from(&g, &it));
    vector_push_back(&to, igraph_es_to(&g, &it));
    igraph_es_next(&g, &it);
  }
  vector_sort(&from); vector_sort(&to);
  vector_print(&from); vector_print(&to);

  igraph_es_destroy(&it);

  vector_destroy(&from);
  vector_destroy(&to);
  
  return 0;
}
