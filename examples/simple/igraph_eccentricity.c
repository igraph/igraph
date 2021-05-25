
#include <igraph.h>

int main() {

    igraph_t g;
    igraph_vector_t ecc;

    igraph_vector_init(&ecc, 0);

    igraph_star(&g, 10, IGRAPH_STAR_UNDIRECTED, 0);
    igraph_eccentricity(&g, &ecc, igraph_vss_all(), IGRAPH_OUT);
    igraph_vector_print(&ecc);
    igraph_destroy(&g);

    igraph_star(&g, 10, IGRAPH_STAR_OUT, 0);
    igraph_eccentricity(&g, &ecc, igraph_vss_all(), IGRAPH_ALL);
    igraph_vector_print(&ecc);
    igraph_destroy(&g);

    igraph_star(&g, 10, IGRAPH_STAR_OUT, 0);
    igraph_eccentricity(&g, &ecc, igraph_vss_all(), IGRAPH_OUT);
    igraph_vector_print(&ecc);
    igraph_destroy(&g);

    igraph_vector_destroy(&ecc);

    return 0;
}
