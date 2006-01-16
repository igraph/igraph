
#include <igraph.h>

int main() {
  
  igraph_d_indheap_t h;
  int i;
  long int idx1, idx2;

  /* igraph_d_indheap_init, igraph_d_indheap_destroy */
  igraph_d_indheap_init(&h, 0);
  igraph_d_indheap_destroy(&h);
  igraph_d_indheap_init(&h, 100);
  igraph_d_indheap_destroy(&h);  

  /* igraph_d_indheap_empty, igraph_d_indheap_size */
  igraph_d_indheap_init(&h, 10);
  if (!igraph_d_indheap_empty(&h)) {
    return 1;
  }
  if (igraph_d_indheap_size(&h) != 0) {
    return 2;
  }
  igraph_d_indheap_push(&h, 10, 0, 0);
  if (igraph_d_indheap_empty(&h)) {
    return 3;
  }
  if (igraph_d_indheap_size(&h) != 1) {
    return 4;
  }
  
  /* igraph_d_indheap_push */
  igraph_d_indheap_push(&h, 9, 9, 8);  
  igraph_d_indheap_push(&h, 3, 3, 2);  
  igraph_d_indheap_push(&h, 2, 2, 1);  
  igraph_d_indheap_push(&h, 1, 1, 0);  
  igraph_d_indheap_push(&h, 7, 7, 6);  
  igraph_d_indheap_push(&h, 4, 4, 3);  
  igraph_d_indheap_push(&h, 0, 0, 1);  
  igraph_d_indheap_push(&h, 6, 6, 5);  
  igraph_d_indheap_push(&h, 5, 5, 4);  
  igraph_d_indheap_push(&h, 8, 8, 7);  

  /* igraph_d_indheap_max, igraph_d_indheap_delete_max */
  while (!igraph_d_indheap_empty(&h)) {
    printf("% li", (long int)igraph_d_indheap_max(&h));
    printf("% li\n", (long int)igraph_d_indheap_delete_max(&h));
  }

  /* igraph_d_indheap_reserve */
  igraph_d_indheap_reserve(&h, 5);
  igraph_d_indheap_reserve(&h, 20);
  igraph_d_indheap_reserve(&h, 0);
  igraph_d_indheap_reserve(&h, 3);
  
  /* igraph_d_indheap_max_index */
  igraph_d_indheap_push(&h, 0, 0, 1);  
  igraph_d_indheap_push(&h, 8, 8, 7);  
  igraph_d_indheap_push(&h, 2, 2, 1);  
  igraph_d_indheap_push(&h, 7, 7, 6);  
  igraph_d_indheap_push(&h, 9, 9, 8);  
  igraph_d_indheap_push(&h, 4, 4, 3);  
  igraph_d_indheap_push(&h, 3, 3, 2);  
  igraph_d_indheap_push(&h, 5, 5, 4);  
  igraph_d_indheap_push(&h, 1, 1, 0);  
  igraph_d_indheap_push(&h, 6, 6, 5);  
  while (!igraph_d_indheap_empty(&h)) {
    igraph_d_indheap_max_index(&h, &idx1, &idx2);
    printf(" %li %li", idx1, idx2);
    igraph_d_indheap_delete_max(&h);
  }
  printf("\n");
  igraph_d_indheap_destroy(&h);
  
  return 0;
}
