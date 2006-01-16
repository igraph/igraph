
#include <igraph.h>

void print_matrix(igraph_matrix_t *m, FILE *f) {
  long int i, j;
  for (i=0; i<igraph_matrix_nrow(m); i++) {
    for (j=0; j<igraph_matrix_ncol(m); j++) {
      fprintf(f, " %li", (long int)MATRIX(*m, i, j));
    }
    fprintf(f, "\n");
  }  
}

int main() {
  igraph_matrix_t m, m1;
  long int i, j;
  
  /* igraph_matrix_init, igraph_matrix_destroy */
  igraph_matrix_init(&m, 10, 10);
  igraph_matrix_destroy(&m);
  
  igraph_matrix_init(&m, 0, 0);
  igraph_matrix_destroy(&m);
  
  /* igraph_matrix_ncol, igraph_matrix_nrow */
  igraph_matrix_init(&m, 10, 5);
  if (igraph_matrix_nrow(&m) != 10) {
    return 1;
  }
  if (igraph_matrix_ncol(&m) != 5) {
    return 2;
  }

  /* igraph_matrix_size, igraph_matrix_resize */
  igraph_matrix_resize(&m, 6, 5);
  if (igraph_matrix_size(&m) != 30) {
    return 3;
  }
  if (igraph_matrix_nrow(&m) != 6) {
    return 4;
  }
  if (igraph_matrix_ncol(&m) != 5) {
    return 5;
  }
  igraph_matrix_resize(&m, 2, 4);
  if (igraph_matrix_nrow(&m) != 2) {
    return 6;
  }
  if (igraph_matrix_ncol(&m) != 4) {
    return 7;
  }
  igraph_matrix_destroy(&m);
  
  /* MATRIX, igraph_matrix_null */
  igraph_matrix_init(&m, 3, 4);
  for (i=0; i<igraph_matrix_nrow(&m); i++) {
    for (j=0; j<igraph_matrix_ncol(&m); j++) {
      MATRIX(m, i, j)= i+1;
    }
  }
  print_matrix(&m, stdout);
  igraph_matrix_null(&m);
  print_matrix(&m, stdout);
  igraph_matrix_destroy(&m);
  
  /* igraph_matrix_add_cols, igraph_matrix_add_rows */
  igraph_matrix_init(&m, 4, 3);
  for (i=0; i<igraph_matrix_nrow(&m); i++) {
    for (j=0; j<igraph_matrix_ncol(&m); j++) {
      MATRIX(m, i, j)= (i+1)*(j+1);
    }
  }
  igraph_matrix_add_cols(&m, 2);
  igraph_matrix_add_rows(&m, 2);
  if (igraph_matrix_ncol(&m) != 5) {
    return 8;
  }
  if (igraph_matrix_nrow(&m) != 6) {
    return 9;
  }
  igraph_matrix_destroy(&m);

  /* igraph_matrix_remove_col */
  igraph_matrix_init(&m, 5, 3);
  for (i=0; i<igraph_matrix_nrow(&m); i++) {
    for (j=0; j<igraph_matrix_ncol(&m); j++) {
      MATRIX(m, i, j)= (i+1)*(j+1);
    }
  }
  igraph_matrix_remove_col(&m, 0);
  print_matrix(&m, stdout);
  igraph_matrix_remove_col(&m, 1);
  print_matrix(&m, stdout);
  igraph_matrix_destroy(&m);

  /* TODO: igraph_matrix_permdelete_rows */
  /* TODO: igraph_matrix_delete_rows_neg */

  /* igraph_matrix_copy */
  igraph_matrix_init(&m, 2, 3);
  for (i=0; i<igraph_matrix_nrow(&m); i++) {
    for (j=0; j<igraph_matrix_ncol(&m); j++) {
      MATRIX(m, i, j)= (i+1)*(j+1);
    }
  }
  igraph_matrix_copy(&m1, &m);
  print_matrix(&m1, stdout);
  igraph_matrix_destroy(&m);
  igraph_matrix_destroy(&m1);

  return 0;
}
