#include <igraph.h>
#include <stdio.h>

int main()
{
    igraph_t graph;
    igraph_vector_int_t dimvector;

    igraph_vector_int_init(&dimvector, 1);
    VECTOR(dimvector)[0] = 5;


    igraph_triangulated_mesh(&graph, &dimvector, true, true);
    igraph_integer_t i, size = igraph_vector_int_size(&graph.from);

    printf("Triangle-like triangulated mesh with dims=(5)):\n");
    for (i = 0; i < size; i++) {
        printf("(%d, %d), ", VECTOR(graph.from)[i], VECTOR(graph.to)[i]);
    }
    printf("\n\n");

    igraph_destroy(&graph);

    igraph_vector_int_init(&dimvector, 2);
    VECTOR(dimvector)[0] = 4;
    VECTOR(dimvector)[1] = 5;

    igraph_triangulated_mesh(&graph, &dimvector, true, true);
    size = igraph_vector_int_size(&graph.from);

    printf("Rectangle-like triangulated mesh with dims=(4, 5)):\n");
    for (i = 0; i < size; i++) {
        printf("(%d, %d), ", VECTOR(graph.from)[i], VECTOR(graph.to)[i]);
    }
    printf("\n\n");

    igraph_destroy(&graph);

    igraph_vector_int_init(&dimvector, 3);
    VECTOR(dimvector)[0] = 3;
    VECTOR(dimvector)[1] = 4;
    VECTOR(dimvector)[2] = 5;

    igraph_triangulated_mesh(&graph, &dimvector, true, true);
    size = igraph_vector_int_size(&graph.from);

    printf("Hexagon-like triangulated mesh with dims=(3, 4, 5)):\n");
    for (i = 0; i < size; i++) {
        printf("(%d, %d), ", VECTOR(graph.from)[i], VECTOR(graph.to)[i]);
    }

    igraph_destroy(&graph);

    igraph_vector_int_destroy(&dimvector);

    return 0;
}
