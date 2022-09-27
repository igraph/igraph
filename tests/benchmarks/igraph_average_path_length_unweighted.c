
#include <igraph.h>

#include "bench.h"

int main(void) {
    igraph_t graph;
    igraph_real_t avglen;
    igraph_matrix_t mat;

    BENCH_INIT();
    igraph_rng_seed(igraph_rng_default(), 42);

    igraph_matrix_init(&mat, 0, 0);

    igraph_kautz(&graph, 4, 5);
    igraph_matrix_resize(&mat, igraph_vcount(&graph), igraph_vcount(&graph)); /* preallocate matrix */

    BENCH(" 1 Kautz(4, 5) average_path_length directed",
          igraph_average_path_length(&graph, &avglen, NULL, IGRAPH_DIRECTED, 1);
    );
    BENCH(" 2 Kautz(4, 5) distances directed",
          igraph_distances(&graph, &mat, igraph_vss_all(), igraph_vss_all(), IGRAPH_OUT);
    );
    BENCH(" 3 Kautz(4, 5) average_path_length undirected",
          igraph_average_path_length(&graph, &avglen, NULL, IGRAPH_UNDIRECTED, 1);
    );
    BENCH(" 4 Kautz(4, 5) distances undirected",
          igraph_distances(&graph, &mat, igraph_vss_all(), igraph_vss_all(), IGRAPH_ALL);
    );

    igraph_destroy(&graph);

    {
        igraph_vector_int_t dims;
        igraph_vector_bool_t periodic;
        igraph_vector_int_init_int(&dims, 3, 15, 15, 15);
        igraph_vector_bool_init_int(&periodic, 3, 1, 1, 1);
        igraph_square_lattice(&graph, &dims, 1, IGRAPH_UNDIRECTED, 0, &periodic);
        igraph_vector_int_destroy(&dims);
        igraph_vector_bool_destroy(&periodic);
        igraph_rewire(&graph, 100, IGRAPH_REWIRING_SIMPLE);
        igraph_matrix_resize(&mat, igraph_vcount(&graph), igraph_vcount(&graph)); /* preallocate matrix */
    }

    BENCH(" 5 Rewired 15x15x15 lattice average_path_length",
          igraph_average_path_length(&graph, &avglen, NULL, IGRAPH_UNDIRECTED, 1);
    );
    BENCH(" 6 Rewired 15x15x15 lattice distances undirected",
          igraph_distances(&graph, &mat, igraph_vss_all(), igraph_vss_all(), IGRAPH_ALL);
    );

    igraph_destroy(&graph);

    igraph_erdos_renyi_game_gnm(&graph, 10000, 12000, IGRAPH_DIRECTED, IGRAPH_LOOPS);
    igraph_matrix_resize(&mat, igraph_vcount(&graph), igraph_vcount(&graph)); /* preallocate matrix */

    BENCH(" 7 Erdos-Renyi n=10000 m=12000 average_path_length directed",
          igraph_average_path_length(&graph, &avglen, NULL, IGRAPH_DIRECTED, 1);
    );
    BENCH(" 8 Erdos-Renyi n=10000 m=12000 distances directed",
          igraph_distances(&graph, &mat, igraph_vss_all(), igraph_vss_all(), IGRAPH_OUT);
    );

    /* The undirected computation will be much slower on this graph, as the largest weakly connected
     * component is much larger. */
    BENCH(" 9 Erdos-Renyi n=10000 m=12000 average_path_length undirected",
          igraph_average_path_length(&graph, &avglen, NULL, IGRAPH_UNDIRECTED, 1);
    );
    BENCH("10 Erdos-Renyi n=10000 m=12000 distances undirected",
          igraph_distances(&graph, &mat, igraph_vss_all(), igraph_vss_all(), IGRAPH_ALL);
    );

    igraph_destroy(&graph);
    igraph_matrix_destroy(&mat);

    return 0;
}
