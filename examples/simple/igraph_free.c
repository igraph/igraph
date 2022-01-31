#include <igraph.h>

int main(void)
{
   igraph_t graph;
   igraph_vector_ptr_t blocks;
   igraph_integer_t i;

   igraph_famous(&graph, "tutte");
   igraph_vector_ptr_init(&blocks, 0);
   igraph_cohesive_blocks(&graph, &blocks, 0, 0, 0);

   for (i = 0; i < igraph_vector_ptr_size(&blocks); i++) {
     igraph_vector_int_t *v = VECTOR(blocks)[i];
     igraph_vector_int_print(v);
     igraph_vector_int_destroy(v);
     igraph_free(v);
   }

   igraph_vector_ptr_destroy(&blocks);
   igraph_destroy(&graph);
   return 0;
}
