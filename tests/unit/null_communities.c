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

/* Testing community detection on the null graph:
 *
 *  - The modularity should be NaN.
 *  - In hierarchical methods, the modularity vector should have size 1.
 */

int main(void) {
    igraph_t g;
    igraph_vector_t modularity;
    igraph_vector_int_t membership;
    igraph_matrix_int_t merges;
    igraph_real_t m;
    igraph_arpack_options_t ao;

    igraph_empty(&g, 0, IGRAPH_UNDIRECTED);

    igraph_vector_init(&modularity, 2);
    igraph_vector_int_init(&membership, 1);
    igraph_matrix_int_init(&merges, 1, 1);

    /* Edge betweenness */

    igraph_community_edge_betweenness(&g, NULL, NULL, &merges, NULL, &modularity, &membership, IGRAPH_UNDIRECTED, NULL, NULL);

    IGRAPH_ASSERT(igraph_matrix_int_nrow(&merges) == 0);
    IGRAPH_ASSERT(igraph_matrix_int_ncol(&merges) == 2);
    IGRAPH_ASSERT(igraph_vector_int_size(&membership) == 0);
    IGRAPH_ASSERT(igraph_vector_size(&modularity) == 1);
    IGRAPH_ASSERT(isnan(VECTOR(modularity)[0]));

    /* Fast greedy */

    igraph_vector_resize(&modularity, 2);
    igraph_vector_int_resize(&membership, 1);
    igraph_matrix_int_resize(&merges, 1, 1);

    igraph_community_fastgreedy(&g, NULL, &merges, &modularity, &membership);

    IGRAPH_ASSERT(igraph_matrix_int_nrow(&merges) == 0);
    IGRAPH_ASSERT(igraph_matrix_int_ncol(&merges) == 2);
    IGRAPH_ASSERT(igraph_vector_int_size(&membership) == 0);
    IGRAPH_ASSERT(igraph_vector_size(&modularity) == 1);
    IGRAPH_ASSERT(isnan(VECTOR(modularity)[0]));

    /* Fluid communities */

    igraph_vector_int_resize(&membership, 1);

    igraph_community_fluid_communities(&g, 0, &membership);

    IGRAPH_ASSERT(igraph_vector_int_size(&membership) == 0);

    /* InfoMAP */

    m = 2;
    igraph_vector_int_resize(&membership, 1);

    igraph_community_infomap(&g, NULL, NULL, 3, false, 0, &membership, &m);

    IGRAPH_ASSERT(igraph_vector_int_size(&membership) == 0);

    /* Label propagation */

    igraph_vector_int_resize(&membership, 1);

    igraph_lpa_variant_t variants[3] = {IGRAPH_LPA_DOMINANCE, IGRAPH_LPA_RETENTION, IGRAPH_LPA_FAST};
    for (igraph_int_t i = 0; i < 3; i++) {
        igraph_community_label_propagation(&g, &membership, IGRAPH_ALL, NULL, NULL, NULL, variants[i]);
    }

    IGRAPH_ASSERT(igraph_vector_int_size(&membership) == 0);

    /* Leading eigenvector */

    m = 2;
    igraph_vector_int_resize(&membership, 1);
    igraph_matrix_int_resize(&merges, 1, 1);

    igraph_arpack_options_init(&ao);
    igraph_community_leading_eigenvector(&g, NULL, &merges, &membership, 1, &ao, &m, 0, NULL, NULL, NULL, NULL, NULL);

    IGRAPH_ASSERT(igraph_matrix_int_nrow(&merges) == 0);
    IGRAPH_ASSERT(igraph_matrix_int_ncol(&merges) == 2);
    IGRAPH_ASSERT(igraph_vector_int_size(&membership) == 0);
    IGRAPH_ASSERT(isnan(m));

    /* Leiden */

    m = 2;
    igraph_vector_int_resize(&membership, 1);

    igraph_community_leiden(&g, NULL, NULL, NULL, 1, 0.01, 0, 1, &membership, NULL, &m);

    IGRAPH_ASSERT(igraph_vector_int_size(&membership) == 0);
    IGRAPH_ASSERT(isnan(m));

    /* Multilevel */

    igraph_vector_int_resize(&membership, 1);
    igraph_vector_resize(&modularity, 2);

    igraph_community_multilevel(&g, NULL, 1, &membership, NULL, &modularity);

    IGRAPH_ASSERT(igraph_vector_int_size(&membership) == 0);
    IGRAPH_ASSERT(igraph_vector_size(&modularity) == 1);
    IGRAPH_ASSERT(isnan(VECTOR(modularity)[0]));

    /* Optimal modularity */
    /* Test only when GLPK is available */

    {
        igraph_error_t ret;
        igraph_error_handler_t *handler;

        m = 2;
        igraph_vector_int_resize(&membership, 1);

        handler = igraph_set_error_handler(igraph_error_handler_ignore);
        ret = igraph_community_optimal_modularity(&g, NULL, 1, &m, &membership);
        igraph_set_error_handler(handler);

        if (ret != IGRAPH_UNIMPLEMENTED) {
            IGRAPH_ASSERT(ret == IGRAPH_SUCCESS);
            IGRAPH_ASSERT(igraph_vector_int_size(&membership) == 0);
            IGRAPH_ASSERT(isnan(m));
        }
    }

    /* Spinglass */

    m = 2;
    igraph_vector_int_resize(&membership, 1);

    igraph_community_spinglass(&g, NULL, &m, NULL, &membership, NULL, 5, 0, 1, 0.01, 0.99, IGRAPH_SPINCOMM_UPDATE_SIMPLE, 1, IGRAPH_SPINCOMM_IMP_ORIG, 1);

    IGRAPH_ASSERT(igraph_vector_int_size(&membership) == 0);
    IGRAPH_ASSERT(isnan(m));

    m = 2;
    igraph_vector_int_resize(&membership, 1);

    igraph_community_spinglass(&g, NULL, &m, NULL, &membership, NULL, 5, false, 1.0, 0.01, 0.99, IGRAPH_SPINCOMM_UPDATE_SIMPLE, 1, IGRAPH_SPINCOMM_IMP_NEG, 1);

    IGRAPH_ASSERT(igraph_vector_int_size(&membership) == 0);
    IGRAPH_ASSERT(isnan(m));

    /* Walktrap */

    igraph_vector_resize(&modularity, 2);
    igraph_vector_int_resize(&membership, 1);
    igraph_matrix_int_resize(&merges, 1, 1);

    igraph_community_walktrap(&g, NULL, 4, &merges, &modularity, &membership);

    IGRAPH_ASSERT(igraph_matrix_int_nrow(&merges) == 0);
    IGRAPH_ASSERT(igraph_matrix_int_ncol(&merges) == 2);
    IGRAPH_ASSERT(igraph_vector_int_size(&membership) == 0);
    IGRAPH_ASSERT(igraph_vector_size(&modularity) == 1);
    IGRAPH_ASSERT(isnan(VECTOR(modularity)[0]));

    /* Cleanup */

    igraph_matrix_int_destroy(&merges);
    igraph_vector_int_destroy(&membership);
    igraph_vector_destroy(&modularity);

    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();

    return 0;
}
