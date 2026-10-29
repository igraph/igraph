
#include <igraph.h>

#include "test_utilities.h"

int main(void) {

    igraph_t graph;
    igraph_vector_int_t membership, generators;

    igraph_rng_seed(igraph_rng_default(), 42);

    igraph_vector_int_init(&membership, 0);
    igraph_vector_int_init(&generators, 0);

    printf("Null graph:\n");

    igraph_empty(&graph, 0, IGRAPH_UNDIRECTED);

    igraph_community_voronoi(&graph, &membership, &generators, NULL, NULL, NULL, IGRAPH_ALL, -1);

    printf("Membership: "); print_vector_int(&membership);
    printf("Generators: "); print_vector_int(&generators);

    igraph_destroy(&graph);

    printf("\nSingleton graph:\n");

    igraph_empty(&graph, 1, IGRAPH_UNDIRECTED);

    igraph_community_voronoi(&graph, &membership, &generators, NULL, NULL, NULL, IGRAPH_ALL, -1);

    printf("Membership: "); print_vector_int(&membership);
    printf("Generators: "); print_vector_int(&generators);

    igraph_destroy(&graph);

    printf("\nTwo isolated nodes:\n");

    igraph_empty(&graph, 2, IGRAPH_UNDIRECTED);

    igraph_community_voronoi(&graph, &membership, &generators, NULL, NULL, NULL, IGRAPH_ALL, -1);

    printf("Membership: "); print_vector_int(&membership);
    printf("Generators: "); print_vector_int(&generators);

    {
        igraph_real_t modularity = 0;
        igraph_community_voronoi(&graph, NULL, NULL, &modularity, NULL, NULL, IGRAPH_ALL, -1);
        IGRAPH_ASSERT(isnan(modularity));
    }

    igraph_destroy(&graph);

    printf("\nZachary:\n");

    igraph_famous(&graph, "Zachary");

    igraph_community_voronoi(&graph, &membership, &generators, NULL, NULL, NULL, IGRAPH_ALL, -1);

    printf("Membership: "); print_vector_int(&membership);
    printf("Generators: "); print_vector_int(&generators);

    {
        igraph_vector_t betw, ibetw;

        igraph_vector_init(&betw, 0);

        igraph_edge_betweenness(&graph, NULL, &betw, igraph_ess_all(IGRAPH_EDGEORDER_ID), IGRAPH_UNDIRECTED, false);

        igraph_int_t n = igraph_vector_size(&betw);
        igraph_vector_init_copy(&ibetw, &betw);
        for (igraph_int_t i=0; i < n; i++) {
            VECTOR(ibetw)[i] = 1.0 / VECTOR(betw)[i];
        }

        printf("\nZachary, betweenness weighted:\n");

        igraph_community_voronoi(&graph, &membership, &generators, NULL, /*lengths=*/ &betw, /*weights=*/ &ibetw, IGRAPH_ALL, -1);

        printf("Membership: "); print_vector_int(&membership);
        printf("Generators: "); print_vector_int(&generators);

        igraph_vector_destroy(&ibetw);
        igraph_vector_destroy(&betw);
    }

    igraph_destroy(&graph);

    {
        igraph_t g1, g2;

        igraph_erdos_renyi_game_gnm(&g1, 10, 30, IGRAPH_DIRECTED, IGRAPH_SIMPLE_SW, IGRAPH_EDGE_UNLABELED);
        igraph_erdos_renyi_game_gnm(&g2, 10, 30, IGRAPH_DIRECTED, IGRAPH_SIMPLE_SW, IGRAPH_EDGE_UNLABELED);

        igraph_disjoint_union(&graph, &g1, &g2);
        igraph_add_edge(&graph, 9, 10);
        igraph_add_edge(&graph, 10, 9);

        igraph_destroy(&g1);
        igraph_destroy(&g2);

        printf("\nDirected two-block, in:\n");

        igraph_community_voronoi(&graph, &membership, &generators, NULL, NULL, NULL, IGRAPH_IN, -1);

        printf("Membership: "); print_vector_int(&membership);
        printf("Generators: "); print_vector_int(&generators);

        printf("\nDirected two-block, out:\n");

        igraph_community_voronoi(&graph, &membership, &generators, NULL, NULL, NULL, IGRAPH_OUT, -1);

        printf("Membership: "); print_vector_int(&membership);
        printf("Generators: "); print_vector_int(&generators);

        igraph_destroy(&graph);
    }

    /* Modularity with a fixed radius must not depend on whether membership is requested. */
    {
        igraph_real_t mod1, mod2;

        igraph_famous(&graph, "Zachary");

        igraph_rng_seed(igraph_rng_default(), 137);
        igraph_community_voronoi(&graph, &membership, NULL, &mod1, NULL, NULL, IGRAPH_ALL, 0.5);

        igraph_rng_seed(igraph_rng_default(), 137);
        igraph_community_voronoi(&graph, NULL, NULL, &mod2, NULL, NULL, IGRAPH_ALL, 0.5);

        IGRAPH_ASSERT(igraph_vector_int_max(&membership) > 0);
        IGRAPH_ASSERT(mod1 > 0);
        IGRAPH_ASSERT(mod1 == mod2);

        igraph_destroy(&graph);
    }

    igraph_vector_int_destroy(&generators);
    igraph_vector_int_destroy(&membership);

    VERIFY_FINALLY_STACK();

    return 0;
}
