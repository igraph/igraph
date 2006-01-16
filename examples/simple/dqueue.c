
#include <igraph.h>

int main() {
  
  igraph_dqueue_t q;
  int i;

  /* igraph_dqueue_init, igraph_dqueue_destroy, igraph_dqueue_empty */
  igraph_dqueue_init(&q, 5);
  if (!igraph_dqueue_empty(&q)) {
    return 1;
  }
  igraph_dqueue_destroy(&q);

  /* igraph_dqueue_push, igraph_dqueue_pop */
  igraph_dqueue_init(&q, 4);
  igraph_dqueue_push(&q, 1);
  igraph_dqueue_push(&q, 2);
  igraph_dqueue_push(&q, 3);
  igraph_dqueue_push(&q, 4);
  if (igraph_dqueue_pop(&q) != 1) {
    return 2;
  }
  if (igraph_dqueue_pop(&q) != 2) {
    return 3;
  }
  if (igraph_dqueue_pop(&q) != 3) {
    return 4;
  }
  if (igraph_dqueue_pop(&q) != 4) {
    return 5;
  }
  igraph_dqueue_destroy(&q);

  /* igraph_dqueue_clear, igraph_dqueue_size */
  igraph_dqueue_init(&q, 0);
  if (igraph_dqueue_size(&q) != 0) {
    return 6;
  }
  igraph_dqueue_clear(&q);
  if (igraph_dqueue_size(&q) != 0) {
    return 7;
  }
  for (i=0; i<10; i++) {
    igraph_dqueue_push(&q, i);
  }
  igraph_dqueue_clear(&q);
  if (igraph_dqueue_size(&q) != 0) {
    return 8;
  }
  igraph_dqueue_destroy(&q);

  /* TODO: igraph_dqueue_full */

  /* igraph_dqueue_head, igraph_dqueue_back, igraph_dqueue_pop_back */
  igraph_dqueue_init(&q, 0);
  for (i=0; i<10; i++) {
    igraph_dqueue_push(&q, i);
  }
  for (i=0; i<10; i++) {
    if (igraph_dqueue_head(&q) != 0) {
      return 9;
    }
    if (igraph_dqueue_back(&q) != 9-i) {
      return 10;
    }
    if (igraph_dqueue_pop_back(&q) != 9-i) {
      return 11;
    }
  }
  igraph_dqueue_destroy(&q);

  return 0;
}
