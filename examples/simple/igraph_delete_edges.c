
#include <igraph.h>

int main() {

  igraph_t g;
  vector_t v;
  int ret;

  vector_init(&v, 8);
  VECTOR(v)[0]=0; VECTOR(v)[1]=1;
  VECTOR(v)[2]=1; VECTOR(v)[3]=2;
  VECTOR(v)[4]=2; VECTOR(v)[5]=3;
  VECTOR(v)[6]=2; VECTOR(v)[7]=2;
  igraph_create(&g, &v, 0, 0);
  
  /* resize vector */
  vector_resize(&v, 2);
  VECTOR(v)[0]=3; VECTOR(v)[1]=2;
  igraph_delete_edges(&g, &v);
  if (igraph_ecount(&g) != 3) {
    return 1;
  }

  /* error test, no such edge to delete */
  igraph_set_error_handler(igraph_error_handler_ignore);
  VECTOR(v)[0]=3; VECTOR(v)[1]=2;
  ret=igraph_delete_edges(&g, &v);
  if (ret != IGRAPH_EINVAL) {
    return 2;
  } 
  if (igraph_ecount(&g) != 3) {
    return 3;
  }

  /* error test, invalid vertex id */
  VECTOR(v)[0]=10;
  ret=igraph_delete_edges(&g, &v);
  if (ret != IGRAPH_EINVVID) {
    return 4;
  } 
  if (igraph_ecount(&g) != 3) {
    return 5;
  }
  
  /* error test, invalid vertex id */
  vector_resize(&v, 3);
  VECTOR(v)[0]=0; VECTOR(v)[1]=1;
  VECTOR(v)[2]=2;
  ret=igraph_delete_edges(&g, &v);
  if (ret != IGRAPH_EINVEVECTOR) {
    return 6;
  } 
  if (igraph_ecount(&g) != 3) {
    return 7;
  }  

  vector_destroy(&v);
  igraph_destroy(&g);
  
  return 0;
}
