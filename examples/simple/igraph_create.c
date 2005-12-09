
#include <igraph.h>

int main() {
  
  igraph_t g;
  vector_t v1, v2;
  int ret;
  
  /* simple use */
  vector_init(&v1, 8);
  VECTOR(v1)[0]=0; VECTOR(v1)[1]=1;
  VECTOR(v1)[2]=1; VECTOR(v1)[3]=2;
  VECTOR(v1)[4]=2; VECTOR(v1)[5]=3;
  VECTOR(v1)[6]=2; VECTOR(v1)[7]=2;
  igraph_create(&g, &v1, 0, 0);
  if (igraph_vcount(&g) != 4) {
    return 1;
  }
  vector_init(&v2, 0);
  igraph_get_edgelist(&g, &v2, 0);
  vector_sort(&v1);
  vector_sort(&v2);
  if (!vector_is_equal(&v1, &v2)) {
    return 2;
  }
  igraph_destroy(&g);
  
  /* higher number of vertices */
  igraph_create(&g, &v1, 10, 0);
  if (igraph_vcount(&g) != 10) {
    return 1;
  }
  igraph_get_edgelist(&g, &v2, 0);
  vector_sort(&v1);
  vector_sort(&v2);
  if (!vector_is_equal(&v1, &v2)) {
    return 3;
  }
  igraph_destroy(&g);

  /* error: IGRAPH_EINVEVECTOR */
  igraph_set_error_handler(igraph_error_handler_ignore);
  vector_resize(&v1, 9);
  VECTOR(v1)[8]=0;
  ret=igraph_create(&g, &v1, 0, 0);
  if (ret != IGRAPH_EINVEVECTOR) {
    return 4;
  }
  
  /* error: IGRAPH_EINVVID */
  vector_resize(&v1, 8);
  VECTOR(v1)[7]=-1;
  ret=igraph_create(&g, &v1, 10, 1);
  if (ret != IGRAPH_EINVVID) {
    return 5;
  }
  vector_destroy(&v1);
  vector_destroy(&v2);

  return 0;
}
