
#include <stdio.h>
#include <igraph.h>

int main() {
  
  igraph_trie_t trie;
  long int id;
  int i;
  char *str;

  /* init */
  igraph_trie_init(&trie, 0);

  /* add and get values */
  igraph_trie_get(&trie, "hello", &id);  printf("hello: %li\n", id);
  igraph_trie_get(&trie, "hepp", &id);   printf("hepp:  %li\n", id);
  igraph_trie_get(&trie, "alma", &id);   printf("alma:  %li\n", id);
  igraph_trie_get(&trie, "also", &id);   printf("also:  %li\n", id);

  igraph_trie_get(&trie, "hello", &id);  printf("hello: %li\n", id);
  igraph_trie_get(&trie, "hepp", &id);   printf("hepp:  %li\n", id);
  igraph_trie_get(&trie, "alma", &id);   printf("alma:  %li\n", id);
  igraph_trie_get(&trie, "also", &id);   printf("also:  %li\n", id);

  igraph_trie_get(&trie, "a", &id);      printf("a:     %li\n", id);
  igraph_trie_get(&trie, "axon", &id);   printf("axon:  %li\n", id);

  igraph_trie_get(&trie, "hello", &id);  printf("hello: %li\n", id);
  igraph_trie_get(&trie, "hepp", &id);   printf("hepp:  %li\n", id);
  igraph_trie_get(&trie, "alma", &id);   printf("alma:  %li\n", id);
  igraph_trie_get(&trie, "also", &id);   printf("also:  %li\n", id);
 
  /* destroy */
  igraph_trie_destroy(&trie);

  /* the same with index */
  igraph_trie_init(&trie, 1);

  igraph_trie_get(&trie, "hello", &id);  printf("hello: %li\n", id);
  igraph_trie_get(&trie, "hepp", &id);   printf("hepp:  %li\n", id);
  igraph_trie_get(&trie, "alma", &id);   printf("alma:  %li\n", id);
  igraph_trie_get(&trie, "also", &id);   printf("also:  %li\n", id);

  igraph_trie_get(&trie, "hello", &id);  printf("hello: %li\n", id);
  igraph_trie_get(&trie, "hepp", &id);   printf("hepp:  %li\n", id);
  igraph_trie_get(&trie, "alma", &id);   printf("alma:  %li\n", id);
  igraph_trie_get(&trie, "also", &id);   printf("also:  %li\n", id);

  igraph_trie_get(&trie, "a", &id);      printf("a:     %li\n", id);
  igraph_trie_get(&trie, "axon", &id);   printf("axon:  %li\n", id);

  igraph_trie_get(&trie, "hello", &id);  printf("hello: %li\n", id);
  igraph_trie_get(&trie, "hepp", &id);   printf("hepp:  %li\n", id);
  igraph_trie_get(&trie, "alma", &id);   printf("alma:  %li\n", id);
  igraph_trie_get(&trie, "also", &id);   printf("also:  %li\n", id);

  for (i=0; i<igraph_trie_size(&trie); i++) {
    igraph_trie_idx(&trie, i, &str);
    printf("%d: %s\n", i, str);
  }
  igraph_trie_destroy(&trie);
  
  return 0;
}
