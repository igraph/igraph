#include <igraph.h>

int main(void)
{
   igraph_t graph;
   igraph_vector_ptr_t seps;
   igraph_integer_t i;

   igraph_famous(&graph, "tutte");
   igraph_vector_ptr_init(&seps, 0);
   igraph_minimum_size_separators(&graph, &seps);

   for (i=0; i<igraph_vector_ptr_size(&seps); i++) {
     igraph_vector_int_t *v=VECTOR(seps)[i];
     igraph_vector_int_print(v);
     igraph_vector_int_destroy(v);
     igraph_free(v);
   }

   igraph_vector_ptr_destroy(&seps);
   igraph_destroy(&graph);
   return 0;
}
