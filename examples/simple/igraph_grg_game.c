
#include <igraph.h>
#include <math.h>

int main() {
    igraph_t graph;
    igraph_vector_t x, y;
    igraph_vector_t weights;
    igraph_eit_t eit;
    igraph_real_t avg_dist;

    /* Set random seed for reproducible results */

    igraph_rng_seed(igraph_rng_default(), 42);

    /* Create a random geometric graph and retrieve vertex coordinates */

    igraph_vector_init(&x, 0);
    igraph_vector_init(&y, 0);

    igraph_grg_game(&graph, 200, 0.1, /* torus */ 0, &x, &y);

    /* Compute edge weights as geometric distance */

    igraph_vector_init(&weights, igraph_ecount(&graph));
    igraph_eit_create(&graph, igraph_ess_all(IGRAPH_EDGEORDER_ID), &eit);
    for (; ! IGRAPH_EIT_END(eit); IGRAPH_EIT_NEXT(eit)) {
        long int e = IGRAPH_EIT_GET(eit);
        long int u = IGRAPH_FROM(&graph, e);
        long int v = IGRAPH_TO(&graph, e);

        VECTOR(weights)[e] = hypot(VECTOR(x)[u] - VECTOR(x)[v], VECTOR(y)[u] - VECTOR(y)[v]);
    }
    igraph_eit_destroy(&eit);

    /* Compute average path length */

    igraph_average_path_length_dijkstra(&graph, &avg_dist, NULL, &weights, IGRAPH_UNDIRECTED, /* unconn */ 1);

    printf("Average distance in the geometric graph: %g.\n", avg_dist);

    /* Destroy data structures when no longer needed */

    igraph_vector_destroy(&weights);
    igraph_destroy(&graph);
    igraph_vector_destroy(&x);
    igraph_vector_destroy(&y);

    return 0;
}
