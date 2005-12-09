
#include <igraph.h>

int main() {
  
  dqueue_t q;
  int i;

  /* dqueue_init, dqueue_destroy, dqueue_empty */
  dqueue_init(&q, 5);
  if (!dqueue_empty(&q)) {
    return 1;
  }
  dqueue_destroy(&q);

  /* dqueue_push, dqueue_pop */
  dqueue_init(&q, 4);
  dqueue_push(&q, 1);
  dqueue_push(&q, 2);
  dqueue_push(&q, 3);
  dqueue_push(&q, 4);
  if (dqueue_pop(&q) != 1) {
    return 2;
  }
  if (dqueue_pop(&q) != 2) {
    return 3;
  }
  if (dqueue_pop(&q) != 3) {
    return 4;
  }
  if (dqueue_pop(&q) != 4) {
    return 5;
  }
  dqueue_destroy(&q);

  /* dqueue_clear, dqueue_size */
  dqueue_init(&q, 0);
  if (dqueue_size(&q) != 0) {
    return 6;
  }
  dqueue_clear(&q);
  if (dqueue_size(&q) != 0) {
    return 7;
  }
  for (i=0; i<10; i++) {
    dqueue_push(&q, i);
  }
  dqueue_clear(&q);
  if (dqueue_size(&q) != 0) {
    return 8;
  }
  dqueue_destroy(&q);

  /* TODO: dqueue_full */

  /* dqueue_head, dqueue_back, dqueue_pop_back */
  dqueue_init(&q, 0);
  for (i=0; i<10; i++) {
    dqueue_push(&q, i);
  }
  for (i=0; i<10; i++) {
    if (dqueue_head(&q) != 0) {
      return 9;
    }
    if (dqueue_back(&q) != 9-i) {
      return 10;
    }
    if (dqueue_pop_back(&q) != 9-i) {
      return 11;
    }
  }
  dqueue_destroy(&q);

  return 0;
}
