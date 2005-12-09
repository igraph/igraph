
#include <igraph.h>

int main() {
  
  d_indheap_t h;
  int i;
  long int idx1, idx2;

  /* d_indheap_init, d_indheap_destroy */
  d_indheap_init(&h, 0);
  d_indheap_destroy(&h);
  d_indheap_init(&h, 100);
  d_indheap_destroy(&h);  

  /* d_indheap_empty, d_indheap_size */
  d_indheap_init(&h, 10);
  if (!d_indheap_empty(&h)) {
    return 1;
  }
  if (d_indheap_size(&h) != 0) {
    return 2;
  }
  d_indheap_push(&h, 10, 0, 0);
  if (d_indheap_empty(&h)) {
    return 3;
  }
  if (d_indheap_size(&h) != 1) {
    return 4;
  }
  
  /* d_indheap_push */
  d_indheap_push(&h, 9, 9, 8);  
  d_indheap_push(&h, 3, 3, 2);  
  d_indheap_push(&h, 2, 2, 1);  
  d_indheap_push(&h, 1, 1, 0);  
  d_indheap_push(&h, 7, 7, 6);  
  d_indheap_push(&h, 4, 4, 3);  
  d_indheap_push(&h, 0, 0, 1);  
  d_indheap_push(&h, 6, 6, 5);  
  d_indheap_push(&h, 5, 5, 4);  
  d_indheap_push(&h, 8, 8, 7);  

  /* d_indheap_max, d_indheap_delete_max */
  while (!d_indheap_empty(&h)) {
    printf("% li", (long int)d_indheap_max(&h));
    printf("% li\n", (long int)d_indheap_delete_max(&h));
  }

  /* d_indheap_reserve */
  d_indheap_reserve(&h, 5);
  d_indheap_reserve(&h, 20);
  d_indheap_reserve(&h, 0);
  d_indheap_reserve(&h, 3);
  
  /* d_indheap_max_index */
  d_indheap_push(&h, 0, 0, 1);  
  d_indheap_push(&h, 8, 8, 7);  
  d_indheap_push(&h, 2, 2, 1);  
  d_indheap_push(&h, 7, 7, 6);  
  d_indheap_push(&h, 9, 9, 8);  
  d_indheap_push(&h, 4, 4, 3);  
  d_indheap_push(&h, 3, 3, 2);  
  d_indheap_push(&h, 5, 5, 4);  
  d_indheap_push(&h, 1, 1, 0);  
  d_indheap_push(&h, 6, 6, 5);  
  while (!d_indheap_empty(&h)) {
    d_indheap_max_index(&h, &idx1, &idx2);
    printf(" %li %li", idx1, idx2);
    d_indheap_delete_max(&h);
  }
  printf("\n");
  d_indheap_destroy(&h);
  
  return 0;
}
