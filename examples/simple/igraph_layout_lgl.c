
#include <igraph.h>
#include <math.h>

int main() {

  igraph_t g;
  igraph_matrix_t coords;
  real_t vc;
  
  igraph_tree(&g, 100, 3, IGRAPH_TREE_UNDIRECTED);
/*   igraph_barabasi_game(&g, 1000, 1, 0, 0, IGRAPH_UNDIRECTED); */
  igraph_matrix_init(&coords, 0, 0);
  vc=igraph_vcount(&g);
  igraph_layout_lgl(&g, &coords, 
		    /* maxiter */    150, 
		    /* maxdelta */   vc,
		    /* area */       vc*vc,
		    /* coolexp */    1.5,
		    /* repulserad */ vc*vc*vc,
		    /* cellsize */   sqrt(sqrt(vc)),
		    /* root */       0);
  
  igraph_matrix_destroy(&coords);
  igraph_destroy(&g);
  return 0;
}
