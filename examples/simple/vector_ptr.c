
#include <igraph.h>
#include <stdlib.h>

int main() {
  
  vector_ptr_t v1, v2;
  const vector_ptr_t v3;
  int i;
  void ** ptr;
  int d1=1, d2=2, d3=3, d4=4, d5=5;

  /* vector_ptr_init, vector_ptr_destroy */
  vector_ptr_init(&v1, 10);
  vector_ptr_destroy(&v1);
  vector_ptr_init(&v1, 0);
  vector_ptr_destroy(&v1);

  /* vector_ptr_free_all, vector_ptr_destroy_all */
  vector_ptr_init(&v1, 5);
  for (i=0; i<vector_ptr_size(&v1); i++) {
    VECTOR(v1)[i]=(void*)malloc(i*10);
  }
  vector_ptr_free_all(&v1);
  for (i=0; i<vector_ptr_size(&v1); i++) {
    VECTOR(v1)[i]=(void*)malloc(i*10);
  }
  vector_ptr_destroy_all(&v1);     
  
  /* vector_ptr_reserve */
  vector_ptr_init(&v1, 0);
  vector_ptr_reserve(&v1, 5);
  vector_ptr_reserve(&v1, 15);
  vector_ptr_reserve(&v1, 1);
  vector_ptr_reserve(&v1, 0);
  vector_ptr_destroy(&v1);

  /* vector_ptr_empty, vector_ptr_clear */
  vector_ptr_init(&v1, 10);
  if (vector_ptr_empty(&v1)) {
    return 1;
  }
  vector_ptr_clear(&v1);
  if (!vector_ptr_empty(&v1)) {
    return 2;
  }

  /* vector_ptr_size */
  if (vector_ptr_size(&v1) != 0) {
    return 3;
  }
  vector_ptr_resize(&v1, 10);
  if (vector_ptr_size(&v1) != 10) {
    return 4;
  }
  vector_ptr_destroy(&v1);

  /* vector_ptr_push_back */
  vector_ptr_init(&v1, 0);
  for (i=0; i<10; i++) {
    vector_ptr_push_back(&v1, (void*)malloc(i*10));
  }
  vector_ptr_destroy_all(&v1);
  
  /* vector_ptr_e */
  vector_ptr_init(&v1, 5);
  VECTOR(v1)[0]=&d1;
  VECTOR(v1)[1]=&d2;
  VECTOR(v1)[2]=&d3;
  VECTOR(v1)[3]=&d4;
  VECTOR(v1)[4]=&d5;
  if (vector_ptr_e(&v1, 0) != &d1) {
    return 5;
  }
  if (vector_ptr_e(&v1, 1) != &d2) {
    return 6;
  }
  if (vector_ptr_e(&v1, 2) != &d3) {
    return 7;
  }
  if (vector_ptr_e(&v1, 3) != &d4) {
    return 8;
  }
  if (vector_ptr_e(&v1, 4) != &d5) {
    return 9;
  }
  vector_ptr_destroy(&v1);

  /* vector_ptr_set */
  vector_ptr_init(&v1, 5);
  vector_ptr_set(&v1, 0, &d1);
  vector_ptr_set(&v1, 1, &d2);
  vector_ptr_set(&v1, 2, &d3);
  vector_ptr_set(&v1, 3, &d4);
  vector_ptr_set(&v1, 4, &d5);
  if (vector_ptr_e(&v1, 0) != &d1) {
    return 5;
  }
  if (vector_ptr_e(&v1, 1) != &d2) {
    return 6;
  }
  if (vector_ptr_e(&v1, 2) != &d3) {
    return 7;
  }
  if (vector_ptr_e(&v1, 3) != &d4) {
    return 8;
  }
  if (vector_ptr_e(&v1, 4) != &d5) {
    return 9;
  }
  vector_ptr_destroy(&v1);

  /* vector_ptr_null */
  vector_ptr_init(&v1, 5);
  vector_ptr_set(&v1, 0, &d1);
  vector_ptr_set(&v1, 1, &d2);
  vector_ptr_set(&v1, 2, &d3);
  vector_ptr_set(&v1, 3, &d4);
  vector_ptr_set(&v1, 4, &d5);
  vector_ptr_null(&v1);
  for (i=0; i<vector_ptr_size(&v1); i++) {
    if (VECTOR(v1)[i] != 0) {
      return 10;
    }
  }
  vector_ptr_destroy(&v1);

  /* vector_ptr_resize */
  vector_ptr_init(&v1, 10);
  vector_ptr_set(&v1, 0, &d1);
  vector_ptr_set(&v1, 1, &d2);
  vector_ptr_set(&v1, 2, &d3);
  vector_ptr_set(&v1, 3, &d4);
  vector_ptr_set(&v1, 4, &d5);
  vector_ptr_resize(&v1, 10);
  vector_ptr_resize(&v1, 15);
  vector_ptr_resize(&v1, 5);
  if (vector_ptr_size(&v1) != 5) {
    return 11;
  }
  if (vector_ptr_e(&v1, 0) != &d1) {
    return 12;
  }
  if (vector_ptr_e(&v1, 1) != &d2) {
    return 13;
  }
  if (vector_ptr_e(&v1, 2) != &d3) {
    return 14;
  }
  if (vector_ptr_e(&v1, 3) != &d4) {
    return 15;
  }
  if (vector_ptr_e(&v1, 4) != &d5) {
    return 16;
  }

  /* vector_ptr_view */
  ptr=(void**) malloc(5 * sizeof(void*));
  vector_ptr_view(&v3, ptr, 5);
  ptr[0]=&d1; ptr[1]=&d2; ptr[2]=&d3; ptr[3]=&d4; ptr[4]=&d5;
  for (i=0; i<vector_ptr_size(&v3); i++) {
    if ( *((int*)VECTOR(v3)[i]) != i+1) {
      return 17;
    }
  }
  
  /* vector_ptr_init_copy */
  vector_ptr_init_copy(&v1, ptr, 5);
  for (i=0; i<vector_ptr_size(&v1); i++) {
    if ( *((int*)VECTOR(v1)[i]) != i+1) {
      return 18;
    }
  }

  /* vector_ptr_copy_to */
  vector_ptr_copy_to(&v1, ptr);
  for (i=0; i<vector_ptr_size(&v1); i++) {
    if ( *((int*)ptr[i]) != i+1) {
      return 19;
    }
  }
  free(ptr);
  vector_ptr_destroy(&v1);

  /* vector_ptr_copy */
  vector_ptr_init(&v1, 5);
  vector_ptr_set(&v1, 0, &d1);
  vector_ptr_set(&v1, 1, &d2);
  vector_ptr_set(&v1, 2, &d3);
  vector_ptr_set(&v1, 3, &d4);
  vector_ptr_set(&v1, 4, &d5);
  vector_ptr_copy(&v2, &v1);
  vector_ptr_destroy(&v1);
  for (i=0; i<vector_ptr_size(&v2); i++) {
    if ( *((int*)VECTOR(v2)[i]) != i+1) {
      return 20;
    }
  }

  /* vector_ptr_remove */
  vector_ptr_remove(&v2, 0);
  vector_ptr_remove(&v2, 3);
  if ( *((int*)VECTOR(v2)[0]) != 2) {
      return 21;
  }
  if ( *((int*)VECTOR(v2)[1]) != 3) {
      return 22;
  }
  if ( *((int*)VECTOR(v2)[2]) != 4) {
      return 23;
  }

  return 0;
}
