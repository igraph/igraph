
#include <igraph.h>

#include "test_utilities.h"

void print_check_destroy(igraph_t *graph, igraph_vector_int_list_t *result) {
    igraph_integer_t i, rank = igraph_vector_int_list_size(result);
    igraph_integer_t ecount = igraph_ecount(graph), vcount = igraph_vcount(graph);
    igraph_integer_t ccount;

    for (i=0; i < rank; ++i) {
        print_vector_int(igraph_vector_int_list_get_ptr(result, i));
    }
    igraph_connected_components(graph, NULL, NULL, &ccount, IGRAPH_WEAK);
    IGRAPH_ASSERT(rank == ecount - vcount + ccount);
    igraph_destroy(graph);
}

int main() {
    igraph_t graph;
    igraph_vector_int_list_t result;
    igraph_integer_t rank;

    igraph_vector_int_list_init(&result, 0);

    printf("Null graph\n");
    igraph_empty(&graph, 0, IGRAPH_UNDIRECTED);
    igraph_minimum_cycle_basis(&graph, /* cutoff */ -1, /* complete */ 1, &result);
    print_check_destroy(&graph, &result);

    printf("\nSingleton graph\n");
    igraph_empty(&graph, 1, IGRAPH_UNDIRECTED);
    igraph_minimum_cycle_basis(&graph, /* cutoff */ -1, /* complete */ 1, &result);
    print_check_destroy(&graph, &result);

    printf("\nSingle vertex with loop\n");
    igraph_small(&graph, 1, IGRAPH_UNDIRECTED,
                 0,0,
                 -1);
    igraph_minimum_cycle_basis(&graph, /* cutoff */ -1, /* complete */ 1, &result);
    print_check_destroy(&graph, &result);

    printf("\nTree\n");
    igraph_kary_tree(&graph, 3, 2, IGRAPH_TREE_UNDIRECTED);
    igraph_minimum_cycle_basis(&graph, /* cutoff */ -1, /* complete */ 1, &result);
    print_check_destroy(&graph, &result);

    printf("\n2-cycle\n");
    igraph_small(&graph, 0, IGRAPH_UNDIRECTED,
                 0,1, 0,1,
                 -1);
    igraph_minimum_cycle_basis(&graph, /* cutoff */ -1, /* complete */ 1, &result);
    print_check_destroy(&graph, &result);

    printf("\nDisconnected\n");
    igraph_small(&graph, 0, IGRAPH_UNDIRECTED,
                 1,2, 2,3, 3,1,
                 4,5, 5,4, 4,5,
                 6,7, 7,8, 8,9, 9,6, 6,8,
                 10,10, 10,11,
                 12,12,
                 -1);
    igraph_minimum_cycle_basis(&graph, /* cutoff */ -1, /* complete */ 1, &result);
    print_check_destroy(&graph, &result);

    printf("\nPeriodic (5,6)-grid\n");


    {
        igraph_vector_int_t dimvec;
        igraph_vector_bool_t periodic;

        igraph_vector_int_init_int_end(&dimvec, -1,
                                       1, 5, 6, -1);

        igraph_vector_bool_init(&periodic, igraph_vector_int_size(&dimvec));
        igraph_vector_bool_fill(&periodic, 1);

        igraph_square_lattice(&graph, &dimvec, 1, IGRAPH_ADJ_UNDIRECTED, 0, &periodic);

        igraph_vector_bool_destroy(&periodic);
        igraph_vector_int_destroy(&dimvec);
    }

    igraph_minimum_cycle_basis(&graph, /* cutoff */ -1, /* complete */ 1, &result);

    rank = igraph_vector_int_list_size(&result);

    /* In a periodic grid graph, all elements in the minimum basis have size 4,
     * except two, which have size equal to the grid dimensions. */

    IGRAPH_ASSERT(igraph_vector_int_size(&VECTOR(result)[0]) == 4);
    IGRAPH_ASSERT(igraph_vector_int_size(&VECTOR(result)[rank-1]) == 6);
    IGRAPH_ASSERT(igraph_vector_int_size(&VECTOR(result)[rank-2]) == 5);
    IGRAPH_ASSERT(igraph_vector_int_size(&VECTOR(result)[rank-3]) == 4);

    print_check_destroy(&graph, &result);

#if 0
    /* This is for benchmarking */
    igraph_rng_seed(igraph_rng_default(), 42);
    igraph_watts_strogatz_game(&graph, 3, 10, 2, 0.1, 0, 0);
    igraph_minimum_cycle_basis(&graph, /* cutoff */ -1, /* complete */ 1, &result);
    print_check_destroy(&graph, &result);
#endif

    igraph_vector_int_list_destroy(&result);

    VERIFY_FINALLY_STACK();

    return 0;
}
