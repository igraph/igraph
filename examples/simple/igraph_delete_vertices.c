
#include <igraph.h>

int main() {

  igraph_t g;
  igraph_vector_t v;
  int ret;

  /* without edges */
  igraph_empty(&g, 5, IGRAPH_DIRECTED);
  igraph_add_vertices(&g, 2);
  igraph_add_vertices(&g, 3);
  igraph_add_vertices(&g, 1);
  igraph_add_vertices(&g, 4);
  if (igraph_vcount(&g) != 15)  {
    return 1;
  }
  igraph_delete_vertices(&g, IGRAPH_VS_1(&g, 2));
  if (igraph_vcount(&g) != 14)  {
    return 1;
  }
  igraph_destroy(&g);
   
  igraph_vector_init(&v, 8);
  VECTOR(v)[0]=0; VECTOR(v)[1]=1;
  VECTOR(v)[2]=1; VECTOR(v)[3]=2;
  VECTOR(v)[4]=2; VECTOR(v)[5]=3;
  VECTOR(v)[6]=2; VECTOR(v)[7]=2;
  igraph_create(&g, &v, 0, 0);
  igraph_vector_destroy(&v);

  /* resize vector */
  igraph_delete_vertices(&g, IGRAPH_VS_1(&g, 2));
  if (igraph_vcount(&g) != 3) {
    return 1;
  }
  if (igraph_ecount(&g) != 1) {
    return 2;
  }

  /* error test */
  igraph_set_error_handler(igraph_error_handler_ignore);
  ret=igraph_delete_vertices(&g, IGRAPH_VS_1(&g, 3));
  if (ret != IGRAPH_EINVVID) {
    return 3;
  }
  
  igraph_destroy(&g);

  return 0;
}
