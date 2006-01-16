
#include <igraph.h>

void igraph_vector_print(const igraph_vector_t *v) {
  long int i;
  for (i=0; i<igraph_vector_size(v); i++) {
    printf("%li ", (long int)VECTOR(*v)[i]);
  }
  printf("\n");
}

int main() {

  igraph_t g, g2;
  const igraph_vector_t v;
  real_t edges[] = { 0,1, 1,2, 2,2, 2,3, 2,4, 3,4 };
  igraph_vector_t v2, v3;
  long int i;
  igraph_vs_t vs;

  igraph_vector_view(&v, edges, sizeof(edges)/sizeof(real_t));
  igraph_create(&g, &v, 0, IGRAPH_DIRECTED);

  /* Create iterator based on a vector (view) */
  igraph_vector_init(&v2, 6);
  VECTOR(v2)[0]=0;   VECTOR(v2)[1]=2;
  VECTOR(v2)[2]=4;   VECTOR(v2)[3]=0;
  VECTOR(v2)[4]=2;   VECTOR(v2)[5]=4;

  igraph_vs_vectorview(&g, &vs, &v2);
  i=0;
  while (!igraph_vs_end(&g, &vs)) {
    if (igraph_vs_get(&g, &vs) != VECTOR(v2)[i]) {
      return 1;
    }
    igraph_vs_next(&g, &vs); 
    i++;
  }
  if (i != igraph_vector_size(&v2)) {
    return 2;
  }

  igraph_vs_destroy(&vs);
  igraph_vector_destroy(&v2);

  /* Create small vector iterator */

  igraph_vs_vector_small(&g, &vs, 0, 2, 4, 0, 2, 4, 2, -1);
  igraph_vector_print(igraph_vs_vector_getvector(&g, &vs));

  igraph_vs_destroy(&vs);

  /* Using a shorthand */  

  igraph_subgraph(&g, &g2, IGRAPH_VS(&g, 0,1,2,3,4,-1));
  igraph_vector_init(&v2, 0); igraph_get_edgelist(&g, &v2, 0); igraph_vector_sort(&v2);
  igraph_vector_init(&v3, 0); igraph_get_edgelist(&g, &v3, 0); igraph_vector_sort(&v3);
  if (!igraph_vector_is_equal(&v2, &v3)) {
    return 3;
  }
  igraph_vector_destroy(&v2);
  igraph_vector_destroy(&v3);

  /* Clean up */

  igraph_destroy(&g);
  igraph_destroy(&g2);

  return 0;
}
