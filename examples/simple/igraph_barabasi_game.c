
#include <igraph.h>

int main() {
  
  igraph_t g;
  vector_t v, v2;
  int i, ret;
  
  igraph_barabasi_game(&g, 10, 2, 0, 0, 1);
  if (igraph_ecount(&g) != 18) {
    return 1;
  }
  if (igraph_vcount(&g) != 10) {
    return 2;
  }
  if (!igraph_is_directed(&g)) {
    return 3;
  }

  vector_init(&v, 0);
  igraph_get_edgelist(&g, &v, 0);
  for (i=0; i<igraph_ecount(&g); i++) {
    if (VECTOR(v)[2*i] <= VECTOR(v)[2*i+1]) {
      return 4;
    }
  }
  igraph_destroy(&g);
  
  /* out degree sequence */
  vector_resize(&v, 10);
  VECTOR(v)[0]=0; VECTOR(v)[1]=1;
  VECTOR(v)[2]=3; VECTOR(v)[3]=3;
  VECTOR(v)[4]=4; VECTOR(v)[5]=5;
  VECTOR(v)[6]=6; VECTOR(v)[7]=7;
  VECTOR(v)[8]=8; VECTOR(v)[9]=9;
  
  igraph_barabasi_game(&g, 10, 0, &v, 0, 1);
  if (igraph_ecount(&g) != vector_sum(&v)) {
    return 5;
  }
  vector_init(&v2, 0);
  igraph_degree(&g, &v2, IGRAPH_VS_ALL, IGRAPH_OUT, 1);
  for (i=0; i<igraph_vcount(&g); i++) {
    if (VECTOR(v)[i] != VECTOR(v2)[i]) {
      return 6;
    }
  }
  vector_destroy(&v);
  vector_destroy(&v2);
  igraph_destroy(&g);
  
  /* outpref, we cannot really test this quantitatively,
     would need to set random seed */
  igraph_barabasi_game(&g, 10, 2, 0, 1, 0);
  vector_init(&v, 0);
  igraph_get_edgelist(&g, &v, 0);
  for (i=0; i<igraph_ecount(&g); i++) {
    if (VECTOR(v)[2*i] <= VECTOR(v)[2*i+1]) {
      return 7;
    }
  }
  if (igraph_is_directed(&g)) {
    return 8;
  }
  igraph_destroy(&g);

  /* Error tests */
  igraph_set_error_handler(igraph_error_handler_ignore);
  ret=igraph_barabasi_game(&g, -10, 1, 0, 0, 0);
  if (ret != IGRAPH_EINVAL) {
    return 9;
  }
  ret=igraph_barabasi_game(&g, 10, -2, 0, 0, 0);
  if (ret != IGRAPH_EINVAL) {
    return 10;
  }
  vector_init(&v, 9);
  ret=igraph_barabasi_game(&g, 10, 0, &v, 0, 0);
  if (ret != IGRAPH_EINVAL) {
    return 11;
  }
  vector_destroy(&v);
  
  return 0;
}
