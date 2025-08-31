/*
   igraph library.
   Copyright (C) 2022  The igraph development team <igraph@igraph.org>

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include <igraph.h>

#include "test_utilities.h"

/* This test verifies that the membership vectors returns by community detection
 * functions are properly indexed, i.e. there are no negative indices and there
 * are no empty communities. */

void check(const igraph_vector_int_t *m) {
    igraph_vector_int_t m2;
    igraph_vector_int_init_copy(&m2, m);
    igraph_reindex_membership(&m2, NULL, NULL);
    IGRAPH_ASSERT(igraph_vector_int_min(m) == 0);
    IGRAPH_ASSERT(igraph_vector_int_max(m) == igraph_vector_int_max(&m2));
    igraph_vector_int_destroy(&m2);
}

int main(void) {
    igraph_t graph;
    igraph_vector_int_t membership;
    igraph_real_t m;
    igraph_error_handler_t *handler;
    igraph_error_t ret;

    igraph_rng_seed(igraph_rng_default(), 42);

    igraph_vector_int_init(&membership, 0);

    igraph_grg_game(&graph, 100, 0.2, false, NULL, NULL);

    igraph_community_fastgreedy(&graph, NULL, NULL, NULL, &membership);
    check(&membership);

    igraph_community_label_propagation(&graph, &membership, IGRAPH_ALL, NULL, NULL, NULL, IGRAPH_LPA_DOMINANCE);
    check(&membership);

    igraph_community_walktrap(&graph, NULL, 4, NULL, NULL, &membership);
    check(&membership);

    igraph_community_edge_betweenness(&graph, NULL, NULL, NULL, NULL, NULL, &membership, IGRAPH_UNDIRECTED, NULL, NULL);
    check(&membership);

    igraph_community_leading_eigenvector(&graph, NULL, NULL, &membership, igraph_vcount(&graph), NULL, NULL, false, NULL, NULL, NULL, NULL, NULL);
    check(&membership);

    igraph_community_leiden(&graph, NULL, NULL, NULL, 1, 0.01, 1, false, &membership, NULL, NULL);
    check(&membership);

    igraph_community_multilevel(&graph, NULL, 1, &membership, NULL, NULL);
    check(&membership);

    igraph_community_fluid_communities(&graph, 10, &membership);
    check(&membership);

    igraph_community_voronoi(&graph, &membership, NULL, NULL, NULL, NULL, IGRAPH_ALL, -1);
    check(&membership);

    igraph_community_spinglass(&graph, NULL, &m, NULL, &membership, NULL, 5, false, 1.0, 0.01, 0.99, IGRAPH_SPINCOMM_UPDATE_SIMPLE, 1, IGRAPH_SPINCOMM_IMP_ORIG, 1);
    check(&membership);

    igraph_community_spinglass(&graph, NULL, &m, NULL, &membership, NULL, 5, false, 1.0, 0.01, 0.99, IGRAPH_SPINCOMM_UPDATE_SIMPLE, 1, IGRAPH_SPINCOMM_IMP_NEG, 1);
    check(&membership);

    igraph_community_infomap(&graph, NULL, NULL, 1, false, 0, &membership, NULL);
    check(&membership);

    igraph_destroy(&graph);

    igraph_grg_game(&graph, 20, 0.5, false, NULL, NULL);

    handler = igraph_set_error_handler(&igraph_error_handler_ignore);
    ret = igraph_community_optimal_modularity(&graph, NULL, 1, NULL, &membership);
    igraph_set_error_handler(handler);
    if (ret != IGRAPH_UNIMPLEMENTED) { /* Test only when GLPK is available */
        IGRAPH_ASSERT(ret == IGRAPH_SUCCESS);
        check(&membership);
    }

    igraph_destroy(&graph);

    igraph_vector_int_destroy(&membership);

    VERIFY_FINALLY_STACK();

    return 0;
}
