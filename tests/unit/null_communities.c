
#include <igraph.h>

#include "test_utilities.inc"

/* Testing community detection on null graph:
 *
 *  - The modularity should be NaN.
 *  - In hierarchical methods, the modularity vector should have size 1.
 */

int main() {
    igraph_t g;
    igraph_vector_t modularity, membership;
    igraph_matrix_t merges;
    igraph_real_t m;
    igraph_integer_t nb_communities;
    igraph_arpack_options_t ao;

    igraph_empty(&g, 0, IGRAPH_UNDIRECTED);

    igraph_vector_init(&modularity, 2);
    igraph_vector_init(&membership, 1);
    igraph_matrix_init(&merges, 1, 1);

    /* Edge betweenness */

    igraph_community_edge_betweenness(&g, NULL, NULL, &merges, NULL, &modularity, &membership, 0, NULL);

    IGRAPH_ASSERT(igraph_matrix_nrow(&merges) == 0);
    IGRAPH_ASSERT(igraph_matrix_ncol(&merges) == 2);
    IGRAPH_ASSERT(igraph_vector_size(&membership) == 0);
    IGRAPH_ASSERT(igraph_vector_size(&modularity) == 1);
    IGRAPH_ASSERT(igraph_is_nan(VECTOR(modularity)[0]));

    /* Fast greedy */

    igraph_vector_resize(&modularity, 2);
    igraph_vector_resize(&membership, 1);
    igraph_matrix_resize(&merges, 1, 1);

    igraph_community_fastgreedy(&g, NULL, &merges, &modularity, &membership);

    IGRAPH_ASSERT(igraph_matrix_nrow(&merges) == 0);
    IGRAPH_ASSERT(igraph_matrix_ncol(&merges) == 2);
    IGRAPH_ASSERT(igraph_vector_size(&membership) == 0);
    IGRAPH_ASSERT(igraph_vector_size(&modularity) == 1);
    IGRAPH_ASSERT(igraph_is_nan(VECTOR(modularity)[0]));

    /* Fluid communities */

    m = 2;
    igraph_vector_resize(&membership, 1);

    igraph_community_fluid_communities(&g, 0, &membership, &m);

    IGRAPH_ASSERT(igraph_vector_size(&membership) == 0);
    IGRAPH_ASSERT(igraph_is_nan(m));

    /* InfoMAP */

    igraph_vector_resize(&membership, 1);

    igraph_community_infomap(&g, NULL, NULL, 3, &membership, &m);

    IGRAPH_ASSERT(igraph_vector_size(&membership) == 0);

    /* Label propagation */

    m = 2;
    igraph_vector_resize(&membership, 1);

    igraph_community_label_propagation(&g, &membership, NULL, NULL, NULL, &m);

    IGRAPH_ASSERT(igraph_vector_size(&membership) == 0);
    IGRAPH_ASSERT(igraph_is_nan(m));

    /* Leading eigenvector */

    m = 2;
    igraph_vector_resize(&membership, 1);
    igraph_matrix_resize(&merges, 1, 1);

    igraph_arpack_options_init(&ao);
    igraph_community_leading_eigenvector(&g, NULL, &merges, &membership, 1, &ao, &m, 0, NULL, NULL, NULL, NULL, NULL);

    IGRAPH_ASSERT(igraph_matrix_nrow(&merges) == 0);
    IGRAPH_ASSERT(igraph_matrix_ncol(&merges) == 2);
    IGRAPH_ASSERT(igraph_vector_size(&membership) == 0);
    IGRAPH_ASSERT(igraph_is_nan(m));

    /* Leiden */

    m = 2;
    igraph_vector_resize(&membership, 1);

    igraph_community_leiden(&g, NULL, NULL, 1, 0.01, 0, &membership, &nb_communities, &m);

    IGRAPH_ASSERT(igraph_vector_size(&membership) == 0);
    IGRAPH_ASSERT(igraph_is_nan(m));

    /* Multilevel */

    igraph_vector_resize(&membership, 1);
    igraph_vector_resize(&modularity, 2);

    igraph_community_multilevel(&g, NULL, 1, &membership, NULL, &modularity);

    IGRAPH_ASSERT(igraph_vector_size(&membership) == 0);
    IGRAPH_ASSERT(igraph_vector_size(&modularity) == 1);
    IGRAPH_ASSERT(igraph_is_nan(VECTOR(modularity)[0]));

    /* Optimal modularity */

    m = 2;
    igraph_vector_resize(&membership, 1);

    igraph_community_optimal_modularity(&g, &m, &membership, NULL);

    IGRAPH_ASSERT(igraph_vector_size(&membership) == 0);
    IGRAPH_ASSERT(igraph_is_nan(m));

    /* Spinglass */

    m = 2;
    igraph_vector_resize(&membership, 1);

    igraph_community_spinglass(&g, NULL, &m, NULL, &membership, NULL, 5, 0, 1, 0.01, 0.99, IGRAPH_SPINCOMM_UPDATE_SIMPLE, 1, IGRAPH_SPINCOMM_IMP_ORIG, 1);

    IGRAPH_ASSERT(igraph_vector_size(&membership) == 0);
    IGRAPH_ASSERT(igraph_is_nan(m));

    m = 2;
    igraph_vector_resize(&membership, 1);

    igraph_community_spinglass(&g, NULL, &m, NULL, &membership, NULL, 5, 0, 1, 0.01, 0.99, IGRAPH_SPINCOMM_UPDATE_SIMPLE, 1, IGRAPH_SPINCOMM_IMP_NEG, 1);

    IGRAPH_ASSERT(igraph_vector_size(&membership) == 0);
    IGRAPH_ASSERT(igraph_is_nan(m));

    /* Walktrap */

    igraph_vector_resize(&modularity, 2);
    igraph_vector_resize(&membership, 1);
    igraph_matrix_resize(&merges, 1, 1);

    igraph_community_walktrap(&g, NULL, 4, &merges, &modularity, &membership);

    IGRAPH_ASSERT(igraph_matrix_nrow(&merges) == 0);
    IGRAPH_ASSERT(igraph_matrix_ncol(&merges) == 2);
    IGRAPH_ASSERT(igraph_vector_size(&membership) == 0);
    IGRAPH_ASSERT(igraph_vector_size(&modularity) == 1);
    IGRAPH_ASSERT(igraph_is_nan(VECTOR(modularity)[0]));

    /* Cleanup */

    igraph_matrix_destroy(&merges);
    igraph_vector_destroy(&membership);
    igraph_vector_destroy(&modularity);

    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();

    return 0;
}
