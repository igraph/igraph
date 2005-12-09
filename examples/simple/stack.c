
#include <igraph.h>

int main() {

  igraph_stack_t st;
  int i;

  /* igraph_stack_init, igraph_stack_destroy */
  igraph_stack_init(&st, 0);
  igraph_stack_destroy(&st);
  igraph_stack_init(&st, 10);
  igraph_stack_destroy(&st);

  /* igraph_stack_reserve */
  igraph_stack_init(&st, 0);
  igraph_stack_reserve(&st, 10);
  igraph_stack_reserve(&st, 5);

  /* igraph_stack_empty */
  if (!igraph_stack_empty(&st)) {
    return 1;
  }
  igraph_stack_push(&st, 1);
  if (igraph_stack_empty(&st)) {
    return 2;
  }

  /* igraph_stack_size */
  if (igraph_stack_size(&st) != 1) {
    return 3;
  }
  for (i=0; i<10; i++) {
    igraph_stack_push(&st, i);
  }
  if (igraph_stack_size(&st) != 11) {
    return 4;
  }

  /* igraph_stack_clear */
  igraph_stack_clear(&st);
  if (!igraph_stack_empty(&st)) {
    return 5;
  }
  igraph_stack_push(&st, 100);
  if (igraph_stack_pop(&st) != 100) {
    return 6;
  }
  igraph_stack_clear(&st);
  igraph_stack_clear(&st);

  /* igraph_stack_push, igraph_stack_pop */
  for (i=0; i<100; i++) {
    igraph_stack_push(&st, 100-i);
  }
  for (i=0; i<100; i++) {
    if (igraph_stack_pop(&st) != i+1) {
      return 7;
    }
  }
  if (!igraph_stack_empty(&st)) {
    return 8;
  }
  
  return 0;
}
