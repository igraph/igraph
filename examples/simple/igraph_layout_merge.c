
#include <igraph.h>
#include <stdlib.h>

int igraph_i_layout_merge_dla(igraph_i_layout_mergegrid_t *grid, 
			      long int actg, real_t *x, real_t *y, real_t r,
			      real_t cx, real_t cy, real_t startr, 
			      real_t killr);

int main() {
  
  /*******************/
  /* Testing the DLA */
  /*******************/
  long int nodes=10;
  igraph_i_layout_mergegrid_t grid;
  igraph_vector_t x, y, r;
  long int i;

  srand(time(0));
  
  igraph_vector_init(&x, nodes);
  igraph_vector_init(&y, nodes);
  igraph_vector_init(&r, nodes);  
  igraph_i_layout_mergegrid_init(&grid, -5, 5, 100, -5, 5, 100);

  /* radius */
  for (i=0; i<nodes; i++) {
    VECTOR(r)[i]=rand()/(double)RAND_MAX;
  }
  igraph_vector_sort(&r);

  /* place */
  VECTOR(x)[0]=0;
  VECTOR(y)[0]=0;
  igraph_i_layout_merge_place_sphere(&grid, 0, 0, VECTOR(r)[nodes-1], 0);

  for (i=1; i<nodes; i++) {
/*     fprintf(stderr, "%li ", i); */
    igraph_i_layout_merge_dla(&grid, i, 
			      igraph_vector_e_ptr(&x, i),
			      igraph_vector_e_ptr(&y, i),
			      VECTOR(r)[nodes-i-1], 0, 0, 4, 7);
    igraph_i_layout_merge_place_sphere(&grid, VECTOR(x)[i], VECTOR(y)[i], 
				       VECTOR(r)[nodes-i-1], i);
  }

/*   for (i=0; i<nodes; i++) {  */
/*     printf("%f %f\n", VECTOR(x)[i], VECTOR(y)[i]); */
/*   } */

/*   print_grid(&grid, stdout); */
  
  igraph_vector_destroy(&x);
  igraph_vector_destroy(&y);
  igraph_vector_destroy(&r);
  igraph_i_layout_mergegrid_destroy(&grid);
  return 0;
}
