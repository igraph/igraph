
#include <igraph.h>

int main() {

  igraph_t g;
  
  igraph_empty(&g, 0, 0);
  if (igraph_is_directed(&g)) {
    return 1;
  }
  igraph_destroy(&g);

  igraph_empty(&g, 0, 1);
  if (!igraph_is_directed(&g)) {
    return 2;
  }
  igraph_destroy(&g);
  
  return 0;
}
