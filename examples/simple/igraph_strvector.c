
#include <igraph.h>

int main() {

  igraph_strvector_t sv1, sv2;
  char *str1, *str2;
  int i;

  /* igraph_strvector_init, igraph_strvector_destroy */
  igraph_strvector_init(&sv1, 10);
  igraph_strvector_destroy(&sv1);
  igraph_strvector_init(&sv1, 0);
  igraph_strvector_destroy(&sv1);

  /* igraph_strvector_get, igraph_strvector_set */
  igraph_strvector_init(&sv1, 5);
  for (i=0; i<igraph_strvector_size(&sv1); i++) {
    igraph_strvector_get(&sv1, i, &str1);
    printf("---%s---\n", str1);
  }
  igraph_strvector_set(&sv1, 0, "zero");
  igraph_strvector_set(&sv1, 1, "one");
  igraph_strvector_set(&sv1, 2, "two");
  igraph_strvector_set(&sv1, 3, "three");
  igraph_strvector_set(&sv1, 4, "four");
  for (i=0; i<igraph_strvector_size(&sv1); i++) {
    igraph_strvector_get(&sv1, i, &str1);
    printf("---%s---\n", str1);
  }

  /* igraph_strvector_remove_section, igraph_strvector_remove, 
     igraph_strvector_resize, igraph_strvector_size */
  igraph_strvector_remove_section(&sv1, 0, 5);
  if (igraph_strvector_size(&sv1) != 0) {
    return 1;
  }
  igraph_strvector_resize(&sv1, 10);
  igraph_strvector_set(&sv1, 0, "zero");
  igraph_strvector_set(&sv1, 1, "one");
  igraph_strvector_set(&sv1, 2, "two");
  igraph_strvector_set(&sv1, 3, "three");
  igraph_strvector_set(&sv1, 4, "four");
  igraph_strvector_resize(&sv1, 5);
  for (i=0; i<igraph_strvector_size(&sv1); i++) {
    igraph_strvector_get(&sv1, i, &str1);
    printf("---%s---\n", str1);
  }

  /* igraph_strvector_move_interval */
  igraph_strvector_move_interval(&sv1, 3, 5, 0);
  for (i=0; i<igraph_strvector_size(&sv1); i++) {
    igraph_strvector_get(&sv1, i, &str1);
    printf("---%s---\n", str1);
  }

  /* igraph_strvector_copy */
  igraph_strvector_copy(&sv2, &sv1);
  for (i=0; i<igraph_strvector_size(&sv2); i++) {
    igraph_strvector_get(&sv2, i, &str1);
    printf("---%s---\n", str1);
  }
  igraph_strvector_resize(&sv1, 0);
  igraph_strvector_destroy(&sv2);
  igraph_strvector_copy(&sv2, &sv1);
  if (igraph_strvector_size(&sv2) != 0) {
    return 2;
  }
  igraph_strvector_destroy(&sv2);

  /* igraph_strvector_add */
  igraph_strvector_add(&sv1, "zeroth");
  igraph_strvector_add(&sv1, "first");
  igraph_strvector_add(&sv1, "second");
  igraph_strvector_add(&sv1, "third");
  igraph_strvector_add(&sv1, "fourth");
  for (i=0; i<igraph_strvector_size(&sv1); i++) {
    igraph_strvector_get(&sv1, i, &str1);
    printf("---%s---\n", str1);
  }

  /* TODO: igraph_strvector_permdelete */
  /* TODO: igraph_strvector_remove_negidx */
  
  return 0;
}
