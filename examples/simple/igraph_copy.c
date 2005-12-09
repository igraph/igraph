
#include <igraph.h>

int main() {

  igraph_t g1, g2;
  vector_t v1, v2;

  vector_init(&v1, 8);
  VECTOR(v1)[0]=0; VECTOR(v1)[1]=1;
  VECTOR(v1)[2]=1; VECTOR(v1)[3]=2;
  VECTOR(v1)[4]=2; VECTOR(v1)[5]=3;
  VECTOR(v1)[6]=2; VECTOR(v1)[7]=2;

  igraph_create(&g1, &v1, 0, 0);
  igraph_copy(&g2, &g1);

  vector_init(&v2, 0);
  igraph_get_edgelist(&g2, &v2, 0);
  if (!vector_is_equal(&v1, &v2)) {
    return 1;
  }

  vector_destroy(&v1);
  vector_destroy(&v2);
  igraph_destroy(&g1);
  igraph_destroy(&g2);

  return 0;
}
