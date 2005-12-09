
#include <igraph.h>

void print_matrix(matrix_t *m, FILE *f) {
  long int i, j;
  for (i=0; i<matrix_nrow(m); i++) {
    for (j=0; j<matrix_ncol(m); j++) {
      fprintf(f, " %li", (long int)MATRIX(*m, i, j));
    }
    fprintf(f, "\n");
  }  
}

int main() {
  matrix_t m, m1;
  long int i, j;
  
  /* matrix_init, matrix_destroy */
  matrix_init(&m, 10, 10);
  matrix_destroy(&m);
  
  matrix_init(&m, 0, 0);
  matrix_destroy(&m);
  
  /* matrix_ncol, matrix_nrow */
  matrix_init(&m, 10, 5);
  if (matrix_nrow(&m) != 10) {
    return 1;
  }
  if (matrix_ncol(&m) != 5) {
    return 2;
  }

  /* matrix_size, matrix_resize */
  matrix_resize(&m, 6, 5);
  if (matrix_size(&m) != 30) {
    return 3;
  }
  if (matrix_nrow(&m) != 6) {
    return 4;
  }
  if (matrix_ncol(&m) != 5) {
    return 5;
  }
  matrix_resize(&m, 2, 4);
  if (matrix_nrow(&m) != 2) {
    return 6;
  }
  if (matrix_ncol(&m) != 4) {
    return 7;
  }
  matrix_destroy(&m);
  
  /* MATRIX, matrix_null */
  matrix_init(&m, 3, 4);
  for (i=0; i<matrix_nrow(&m); i++) {
    for (j=0; j<matrix_ncol(&m); j++) {
      MATRIX(m, i, j)= i+1;
    }
  }
  print_matrix(&m, stdout);
  matrix_null(&m);
  print_matrix(&m, stdout);
  matrix_destroy(&m);
  
  /* matrix_add_cols, matrix_add_rows */
  matrix_init(&m, 4, 3);
  for (i=0; i<matrix_nrow(&m); i++) {
    for (j=0; j<matrix_ncol(&m); j++) {
      MATRIX(m, i, j)= (i+1)*(j+1);
    }
  }
  matrix_add_cols(&m, 2);
  matrix_add_rows(&m, 2);
  if (matrix_ncol(&m) != 5) {
    return 8;
  }
  if (matrix_nrow(&m) != 6) {
    return 9;
  }
  matrix_destroy(&m);

  /* matrix_remove_col */
  matrix_init(&m, 5, 3);
  for (i=0; i<matrix_nrow(&m); i++) {
    for (j=0; j<matrix_ncol(&m); j++) {
      MATRIX(m, i, j)= (i+1)*(j+1);
    }
  }
  matrix_remove_col(&m, 0);
  print_matrix(&m, stdout);
  matrix_remove_col(&m, 1);
  print_matrix(&m, stdout);
  matrix_destroy(&m);

  /* TODO: matrix_permdelete_rows */
  /* TODO: matrix_delete_rows_neg */

  /* matrix_copy */
  matrix_init(&m, 2, 3);
  for (i=0; i<matrix_nrow(&m); i++) {
    for (j=0; j<matrix_ncol(&m); j++) {
      MATRIX(m, i, j)= (i+1)*(j+1);
    }
  }
  matrix_copy(&m1, &m);
  print_matrix(&m1, stdout);
  matrix_destroy(&m);
  matrix_destroy(&m1);

  return 0;
}
