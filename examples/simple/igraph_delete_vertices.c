
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
  vector_destroy(&v);

  /* resize vector */
  igraph_delete_vertices(&g, IGRAPH_VS_1(2));
  if (igraph_vcount(&g) != 3) {
    return 1;
  }
  if (igraph_ecount(&g) != 1) {
    return 2;
  }

  /* error test */
  igraph_set_error_handler(igraph_error_handler_ignore);
  ret=igraph_delete_vertices(&g, IGRAPH_VS_1(3));
  if (ret != IGRAPH_EINVVID) {
    return 3;
  }
  
  igraph_destroy(&g);

  return 0;
}
