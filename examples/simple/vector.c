
#include <igraph.h>

void print_vector(vector_t *v, FILE *f) {
  long int i;
  for (i=0; i<vector_size(v); i++) {
    fprintf(f, " %li", (long int) VECTOR(*v)[i]);
  }
  fprintf(f, "\n");
}

int main() {
  
  vector_t v, v2;
  int i;
  real_t *ptr;

  /* simple init */
  vector_init(&v, 0);
  vector_destroy(&v);

  /* VECTOR(), vector_size */
  vector_init(&v, 10);
  for (i=0; i<vector_size(&v); i++) {
    VECTOR(v)[i]=10-i;
  }
  print_vector(&v, stdout);
  vector_destroy(&v);

  /* vector_reserve, vector_push_back */
  vector_init(&v, 0);
  vector_reserve(&v, 10);
  for (i=0; i<10; i++) {
    vector_push_back(&v, i);
  }

  /* vector_empty, vector_clear */
  if (vector_empty(&v)) {
    return 1;
  }
  vector_clear(&v);
  if (!vector_empty(&v)) {
    return 2;
  }
  vector_destroy(&v);

  /* vector_e, vector_e_ptr */
  vector_init(&v, 5);
  for (i=0; i<vector_size(&v); i++) {
    *vector_e_ptr(&v, i) = 100*i;
  }
  for (i=0; i<vector_size(&v); i++) {
    fprintf(stdout, " %li", (long int)vector_e(&v, i));
  }
  fprintf(stdout, "\n");
  vector_destroy(&v);

  /* vector_set */
  vector_init(&v, 5);
  for (i=0; i<vector_size(&v); i++) {
    vector_set(&v, i, 20*i);
  }
  print_vector(&v, stdout);
  vector_destroy(&v);

  /* vector_null */
  vector_init(&v, 0);
  vector_null(&v);
  vector_destroy(&v);
  vector_init(&v, 10);
  for (i=0; i<vector_size(&v); i++) {
    VECTOR(v)[i]=i+1;
  }
  vector_null(&v);
  print_vector(&v, stdout);
  vector_destroy(&v);

  /* vector_tail, vector_pop_back */
  vector_init(&v, 10);
  for (i=0; i<vector_size(&v); i++) {
    VECTOR(v)[i]=i+1;
  }
  while (!vector_empty(&v)) {
    fprintf(stdout, " %li", (long int)vector_tail(&v));
    fprintf(stdout, " %li", (long int)vector_pop_back(&v));
  }
  fprintf(stdout, "\n");
  vector_destroy(&v);

  /* vector_init_seq, vector_order */
  vector_init_seq(&v, 1, 10);
  vector_init(&v2, 0);
  vector_order(&v, &v2, 10);
  print_vector(&v2, stdout);
  vector_destroy(&v2);
  vector_destroy(&v);

  /* vector_resize, vector_sort */
  vector_init(&v, 20);  
  for (i=0; i<10; i++) {
    VECTOR(v)[i]=10-i;
  }
  vector_resize(&v, 10);
  vector_sort(&v);
  print_vector(&v, stdout);
  vector_destroy(&v);

  /* vector_as_vector, vector_max, vector_init_copy */
  ptr=(real_t*) malloc(10* sizeof(real_t));
  v=vector_as_vector(ptr, 10);
  for (i=0; i<vector_size(&v); i++) {
    VECTOR(v)[i]=100-i;    
  }
  for (i=0; i<10; i++) {
    fprintf(stdout, " %li", (long int)ptr[i]);
  }
  fprintf(stdout, "\n");
  fprintf(stdout, " %li\n", (long int)vector_max(&v));

  vector_init_copy(&v, ptr, 10);
  free(ptr);
  print_vector(&v, stdout);
  vector_destroy(&v);
  
  /* vector_copy_to */
  ptr=(real_t*) malloc(10* sizeof(real_t));
  vector_init_seq(&v, 11, 20);
  vector_copy_to(&v, ptr);
  for (i=0; i<10; i++) {
    fprintf(stdout, " %li", (long int)ptr[i]);
  }
  fprintf(stdout, "\n");
  free(ptr);
  vector_destroy(&v);

  /* vector_init_seq, vector_sum, vector_prod */
  vector_init_seq(&v, 1, 5);
  fprintf(stdout, " %li", (long int)vector_sum(&v));
  fprintf(stdout, " %li\n", (long int)vector_prod(&v));
  
  /* vector_remove_section */
  vector_remove_section(&v, 2, 4);
  fprintf(stdout, " %li", (long int)vector_sum(&v));
  fprintf(stdout, " %li\n", (long int)vector_prod(&v));
  vector_destroy(&v);

  /* vector_remove */
  vector_init_seq(&v, 1, 10);
  vector_remove(&v, 9);
  vector_remove(&v, 0);
  vector_remove(&v, 4);
  fprintf(stdout, " %li\n", (long int)vector_sum(&v));
  vector_destroy(&v);

  /* vector_move_interval */
  vector_init_seq(&v, 0, 9);
  vector_move_interval(&v, 5, 10, 0);
  if (vector_sum(&v) != 70) {
    return 3;
  }
  vector_destroy(&v);
    
  /* vector_isininterval */
  vector_init_seq(&v, 1, 10);
  if (!vector_isininterval(&v, 1, 10)) {
    return 4;
  }
  if (vector_isininterval(&v, 2, 10)) {
    return 5;
  }
  if (vector_isininterval(&v, 1, 9)) {
    return 6;
  }

  /* vector_any_smaller */
  if (vector_any_smaller(&v, 1)) {
    return 7;
  }
  if (!vector_any_smaller(&v, 2)) {
    return 8;
  }

  /* vector_permdelete */
  /* vector_remove_negidx */

  return 0;
}
  
