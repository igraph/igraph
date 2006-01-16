
#include <igraph.h>

void print_vector(igraph_vector_t *v, FILE *f) {
  long int i;
  for (i=0; i<igraph_vector_size(v); i++) {
    fprintf(f, " %li", (long int) VECTOR(*v)[i]);
  }
  fprintf(f, "\n");
}

int main() {

  igraph_t g;
  igraph_vector_t v, seq;
  int ret;

  /* Create graph */
  igraph_vector_init(&v, 8);
  VECTOR(v)[0]=0; VECTOR(v)[1]=1;
  VECTOR(v)[2]=1; VECTOR(v)[3]=2;
  VECTOR(v)[4]=2; VECTOR(v)[5]=3;
  VECTOR(v)[6]=2; VECTOR(v)[7]=2;
  igraph_create(&g, &v, 0, IGRAPH_DIRECTED);
  
  igraph_degree(&g, &v, IGRAPH_VS_ALL(&g), IGRAPH_OUT, IGRAPH_NO_LOOPS);
  print_vector(&v, stdout);

  igraph_degree(&g, &v, IGRAPH_VS_ALL(&g), IGRAPH_OUT, IGRAPH_LOOPS);
  print_vector(&v, stdout);
  
  igraph_degree(&g, &v, IGRAPH_VS_ALL(&g), IGRAPH_IN, IGRAPH_NO_LOOPS);
  print_vector(&v, stdout);

  igraph_degree(&g, &v, IGRAPH_VS_ALL(&g), IGRAPH_IN, IGRAPH_LOOPS);
  print_vector(&v, stdout);
  
  igraph_degree(&g, &v, IGRAPH_VS_ALL(&g), IGRAPH_ALL, IGRAPH_NO_LOOPS);
  print_vector(&v, stdout);

  igraph_degree(&g, &v, IGRAPH_VS_ALL(&g), IGRAPH_ALL, IGRAPH_LOOPS);
  print_vector(&v, stdout);
  
  igraph_destroy(&g);
  
  igraph_vector_resize(&v, 8);
  VECTOR(v)[0]=0; VECTOR(v)[1]=1;
  VECTOR(v)[2]=1; VECTOR(v)[3]=2;
  VECTOR(v)[4]=2; VECTOR(v)[5]=3;
  VECTOR(v)[6]=2; VECTOR(v)[7]=2;
  igraph_create(&g, &v, 0, IGRAPH_UNDIRECTED);

  igraph_degree(&g, &v, IGRAPH_VS_ALL(&g), IGRAPH_OUT, IGRAPH_NO_LOOPS);
  print_vector(&v, stdout);

  igraph_degree(&g, &v, IGRAPH_VS_ALL(&g), IGRAPH_OUT, IGRAPH_LOOPS);
  print_vector(&v, stdout);
  
  igraph_degree(&g, &v, IGRAPH_VS_ALL(&g), IGRAPH_IN, IGRAPH_NO_LOOPS);
  print_vector(&v, stdout);

  igraph_degree(&g, &v, IGRAPH_VS_ALL(&g), IGRAPH_IN, IGRAPH_LOOPS);
  print_vector(&v, stdout);
  
  igraph_degree(&g, &v, IGRAPH_VS_ALL(&g), IGRAPH_ALL, IGRAPH_NO_LOOPS);
  print_vector(&v, stdout);

  igraph_degree(&g, &v, IGRAPH_VS_ALL(&g), IGRAPH_ALL, IGRAPH_LOOPS);
  print_vector(&v, stdout);

  /* Degree of the same vertex multiple times */
  
  igraph_vector_init(&seq, 3);
  VECTOR(seq)[0]=2; VECTOR(seq)[1]=0; VECTOR(seq)[2]=2;
  igraph_degree(&g, &v, IGRAPH_VS_VECTOR(&g, &seq), IGRAPH_ALL, IGRAPH_LOOPS);
  print_vector(&v, stdout);

  /* Errors */
  igraph_set_error_handler(igraph_error_handler_ignore);
  ret=igraph_degree(&g, &v, IGRAPH_VS_VECTOR(&g, &seq), 0, IGRAPH_LOOPS);
  if (ret != IGRAPH_EINVMODE) {
    return 1;
  }

  VECTOR(seq)[0]=4;
  ret=igraph_degree(&g, &v, IGRAPH_VS_VECTOR(&g, &seq), IGRAPH_ALL, IGRAPH_LOOPS);
  if (ret != IGRAPH_EINVVID) {
    return 2;
  }  

  igraph_destroy(&g);
  igraph_vector_destroy(&v);
  igraph_vector_destroy(&seq);

  return 0;
}
