
#include <igraph.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>

int igraph_i_layout_sphere_2d(igraph_matrix_t *coords, real_t *x, real_t *y,
			      real_t *r);
int igraph_i_layout_sphere_3d(igraph_matrix_t *coords, real_t *x, real_t *y,
			      real_t *z, real_t *r);

int main () {
  long int i; 
  igraph_matrix_t m;
  real_t x, y, z, r;

  srand(time(0));

  /* 2D */
  igraph_matrix_init(&m, 1000, 2);
  for (i=0; i<igraph_matrix_nrow(&m); i++) {
    MATRIX(m,i,0)=rand()/(double)RAND_MAX;
    MATRIX(m,i,1)=rand()/(double)RAND_MAX;
  }
  igraph_i_layout_sphere_2d(&m, &x, &y, &r);
  
  for (i=0; i<igraph_matrix_nrow(&m); i++) {
    real_t dist=sqrt((MATRIX(m,i,0)-x)*(MATRIX(m,i,0)-x) + 
		     (MATRIX(m,i,1)-y)*(MATRIX(m,i,1)-y));
    if (dist > r) {
      printf("x: %f y: %f r: %f\n", x, y, r);
      printf("x: %f y: %f dist: %f (%li)\n", 
	     MATRIX(m,i,0), MATRIX(m,i,1), dist, i);
      return 1;
    }
  }
  igraph_matrix_destroy(&m);

  /* 3D */
  igraph_matrix_init(&m, 1000, 3);
  for (i=0; i<igraph_matrix_nrow(&m); i++) {
    MATRIX(m,i,0)=rand()/(double)RAND_MAX;
    MATRIX(m,i,1)=rand()/(double)RAND_MAX;
    MATRIX(m,i,2)=rand()/(double)RAND_MAX;
  }
  igraph_i_layout_sphere_3d(&m, &x, &y, &z, &r);
  
  for (i=0; i<igraph_matrix_nrow(&m); i++) {
    real_t dist=sqrt((MATRIX(m,i,0)-x)*(MATRIX(m,i,0)-x) + 
		     (MATRIX(m,i,1)-y)*(MATRIX(m,i,1)-y) +
		     (MATRIX(m,i,2)-z)*(MATRIX(m,i,2)-z));
    if (dist > r) {
      printf("x: %f y: %f z: %f r: %f\n", x, y, z, r);
      printf("x: %f y: %f z: %f dist: %f (%li)\n", 
	     MATRIX(m,i,0), MATRIX(m,i,1), MATRIX(m,i,2), dist, i);
      return 1;
    }
  }
  igraph_matrix_destroy(&m);

  

  return 0;
}
