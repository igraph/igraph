/* -*- mode: C -*-  */

#include <igraph.h>
#include <stdio.h>

int main(void)
{
    igraph_t graph;
    igraph_vector_int_t dimvector;

    /* triangular triangulated mesh */
    igraph_vector_int_init(&dimvector, 1);
    VECTOR(dimvector)[0] = 5;

    igraph_triangulated_mesh(&graph, &dimvector, true, false);
    igraph_integer_t size = igraph_vector_int_size(&graph.from);

    igraph_vector_int_print(&graph.from);
    igraph_vector_int_print(&graph.to);

    igraph_destroy(&graph);

    /* rectangular triangulated mesh */
    igraph_vector_int_init(&dimvector, 2);
    VECTOR(dimvector)[0] = 4;
    VECTOR(dimvector)[1] = 5;

    igraph_triangulated_mesh(&graph, &dimvector, true, false);
    size = igraph_vector_int_size(&graph.from);

    printf("\n");
    igraph_vector_int_print(&graph.from);
    igraph_vector_int_print(&graph.to);

    igraph_destroy(&graph);

    /* hexagonal triangulated mesh */
    igraph_vector_int_init(&dimvector, 3);
    VECTOR(dimvector)[0] = 3;
    VECTOR(dimvector)[1] = 4;
    VECTOR(dimvector)[2] = 5;

    igraph_triangulated_mesh(&graph, &dimvector, true, false);
    size = igraph_vector_int_size(&graph.from);

    printf("\n");
    igraph_vector_int_print(&graph.from);
    igraph_vector_int_print(&graph.to);

    igraph_destroy(&graph);

    igraph_vector_int_destroy(&dimvector);

    return 0;
}
