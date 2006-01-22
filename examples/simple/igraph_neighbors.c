
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
  igraph_vector_t v;
  int ret;

  igraph_vector_init(&v, 8);
  VECTOR(v)[0]=0; VECTOR(v)[1]=1;
  VECTOR(v)[2]=1; VECTOR(v)[3]=2;
  VECTOR(v)[4]=2; VECTOR(v)[5]=3;
  VECTOR(v)[6]=2; VECTOR(v)[7]=2;
  igraph_create(&g, &v, 0, 1);

  igraph_neighbors(&g, &v, 2, IGRAPH_OUT);
  igraph_vector_sort(&v);
  print_vector(&v, stdout);
  
  igraph_neighbors(&g, &v, 2, IGRAPH_IN);
  igraph_vector_sort(&v);
  print_vector(&v, stdout);

  igraph_neighbors(&g, &v, 2, IGRAPH_ALL);
  igraph_vector_sort(&v);
  print_vector(&v, stdout);
  
  /* Errors */
  igraph_set_error_handler(igraph_error_handler_ignore);
  ret=igraph_neighbors(&g, &v, 2, (igraph_neimode_t)0); /* conv for c++ */
  if (ret != IGRAPH_EINVMODE) {
    return 1;
  }
  
  ret=igraph_neighbors(&g, &v, 4, IGRAPH_ALL);
  if (ret != IGRAPH_EINVVID) {
    return 2;
  }

  igraph_vector_destroy(&v);
  igraph_destroy(&g);
  return 0;
}
