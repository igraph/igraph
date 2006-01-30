
#include <igraph.h>

int main() {

  igraph_t g;
  igraph_matrix_t coords;
  
  igraph_ring(&g, 100, IGRAPH_UNDIRECTED, 0, 1);
  igraph_matrix_init(&coords, 0, 0);
  igraph_layout_lgl(&g, &coords, 
		    /* maxiter */    150, 
		    /* maxdelta */   100,
		    /* area */       1000,
		    /* coolexp */    1.5,
		    /* repulserad */ 1000,
		    /* cellsize */   10);

  igraph_matrix_destroy(&coords);
  igraph_destroy(&g);
  return 0;
}
