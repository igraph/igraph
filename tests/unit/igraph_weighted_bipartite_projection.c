/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2008-2023  The igraph development team

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include <igraph.h>
#include <stdlib.h>
#include "test_utilities.h"

int check_projection_sizes(
    const char* test_name, const igraph_t *graph, const igraph_vector_bool_t *types,
    const igraph_vector_t *weights, const igraph_t *proj1, const igraph_t *proj2
) {
    igraph_integer_t vcount1, ecount1, vcount2, ecount2;
    igraph_bipartite_projection_size(graph, types, &vcount1, &ecount1,
                                     &vcount2, &ecount2);

    /* igraph_bipartite_projection_size() is for the unweighted projection;
     * the weighted projection has mutual edges */
    ecount1 *= 2;
    ecount2 *= 2;

    if (proj1 && !igraph_is_directed(proj1)) {
        printf("%s: first projection is not directed", test_name);
        exit(14);
    }
    if (proj1 && igraph_vcount(proj1) != vcount1) {
        printf(
            "%s: first projection is expected to have %" IGRAPH_PRId
            " vertices, got %" IGRAPH_PRId,
            test_name, vcount1, igraph_vcount(proj1)
        );
        exit(10);
    }
    if (proj1 && igraph_ecount(proj1) != ecount1) {
        printf(
            "%s: first projection is expected to have %" IGRAPH_PRId
            " edges, got %" IGRAPH_PRId,
            test_name, ecount1, igraph_ecount(proj1)
        );
        exit(11);
    }

    if (proj2 && !igraph_is_directed(proj2)) {
        printf("%s: second projection is not directed", test_name);
        exit(15);
    }
    if (proj2 && igraph_vcount(proj2) != vcount2) {
        printf(
            "%s: second projection is expected to have %" IGRAPH_PRId
            " vertices, got %" IGRAPH_PRId,
            test_name, vcount2, igraph_vcount(proj2)
        );
        exit(12);
    }
    if (proj2 && igraph_ecount(proj2) != ecount2) {
        printf(
            "%s: second projection is expected to have %" IGRAPH_PRId
            " edges, got %" IGRAPH_PRId,
            test_name, ecount2, igraph_ecount(proj2)
        );
        exit(13);
    }

    return 0;
}

int main(void) {

    igraph_t g, p1, p2, full, ring;
    igraph_vector_bool_t types;
    igraph_bool_t iso;
    igraph_integer_t i;
    igraph_vector_t proj_weights1, proj_weights2;
    igraph_vector_t weights;

    /*******************************************************/
    /* Full bipartite graph -> full graphs                 */
    /*******************************************************/

    igraph_vector_init(&weights, 15);
    igraph_vector_fill(&weights, 1);
    igraph_vector_bool_init(&types, 0);

    igraph_full_bipartite(&g, &types, 5, 3, /*directed=*/ 0,
                          /*mode=*/ IGRAPH_ALL);

    /* Get both projections */
    igraph_weighted_bipartite_projection(&g, &types, &weights, &p1, &p2, 0, 0, /*probe1=*/ -1);
    check_projection_sizes("Full bipartite graph", &g, &types, &weights, &p1, &p2);

    /* Check first projection */
    igraph_full(&full, igraph_vcount(&p1), IGRAPH_DIRECTED, /*loops=*/0);
    igraph_isomorphic_bliss(&p1, &full, 0, 0, &iso, 0, 0,
                            IGRAPH_BLISS_FM, 0, 0);
    if (!iso) {
        return 1;
    }
    igraph_destroy(&full);

    /* Check second projection */
    igraph_full(&full, igraph_vcount(&p2), IGRAPH_DIRECTED, /*loops=*/0);
    igraph_isomorphic_bliss(&p2, &full, 0, 0, &iso, 0, 0,
                            IGRAPH_BLISS_FM, 0, 0);
    if (!iso) {
        return 2;
    }
    igraph_destroy(&full);

    igraph_destroy(&p1);
    igraph_destroy(&p2);
    igraph_destroy(&g);
    igraph_vector_bool_destroy(&types);
    igraph_vector_destroy(&weights);

    /*******************************************************/
    /* More sophisticated test                             */
    /*******************************************************/

    igraph_vector_init(&weights, 200);
    igraph_vector_fill(&weights, 1);

    igraph_ring(&g, 100, /*directed=*/ 1, /*mutual=*/ 1,
                /*circular=*/ 1);
    igraph_vector_bool_init(&types, igraph_vcount(&g));
    for (i = 0; i < igraph_vector_bool_size(&types); i++) {
        VECTOR(types)[i] = i % 2 ? 0 : 1;
    }

    /* Get both projections */
    igraph_weighted_bipartite_projection(
        &g, &types, &weights, &p1, &p2, 0, 0, /*probe1=*/ -1
    );
    check_projection_sizes("Ring graph", &g, &types, &weights, &p1, &p2);

    /* Check first projection */
    igraph_ring(&ring, igraph_vcount(&g) / 2, IGRAPH_DIRECTED, /* mutual= */ 1, /*circular=*/ 1);
    igraph_isomorphic_bliss(&p1, &ring, 0, 0, &iso, 0, 0,
                            IGRAPH_BLISS_FM, 0, 0);
    if (!iso) {
        printf("Ring graph: first projection is not isomorphic to the expected result\n");
        printf("Observed result was:\n\n");
        print_graph_canon(&p1);
        return 1;
    }

    /* Check second projection */
    igraph_isomorphic_bliss(&p2, &ring, 0, 0, &iso, 0, 0,
                            IGRAPH_BLISS_FM, 0, 0);
    if (!iso) {
        printf("Ring graph: second projection is not isomorphic to the expected result\n");
        printf("Observed result was:\n\n");
        print_graph_canon(&p2);
        return 2;
    }
    igraph_destroy(&ring);

    igraph_destroy(&p1);
    igraph_destroy(&p2);
    igraph_destroy(&g);
    igraph_vector_bool_destroy(&types);
    igraph_vector_destroy(&weights);

    /************************************/
    /* Examples from Tore Opsahl's page */
    /************************************/

    /* Source: https://toreopsahl.com/tnet/two-mode-networks/projection/ */
    /* Vertices A-F are encoded with vertex IDs 0-5.
     * Orange vertices are encoded with vertex IDs 6-11.
     */
    igraph_vector_init_int_end(
        &weights, -1,
        4, 2, 2, 1, 4, 3, 2, 5, 6, 2, 4, 1, 1,
        -1
    );

    igraph_small(
        &g, 12, IGRAPH_UNDIRECTED,
        0, 6, 0, 7,
        1, 6, 1, 7, 1, 8, 1, 9, 1, 10,
        2, 7,
        3, 8,
        4, 9, 4, 10, 4, 11,
        5, 11,
        -1
    );
    igraph_vector_bool_init_int_end(
        &types, -1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, -1
    );

    /* First we will check the original weights */

    igraph_vector_init(&proj_weights1, 0);
    igraph_vector_init(&proj_weights2, 0);
    igraph_weighted_bipartite_projection(
        &g, &types, &weights, &p1, &p2, &proj_weights1, &proj_weights2, /*probe=*/ -1
    );
    check_projection_sizes("Small example graph", &g, &types, &weights, &p1, &p2);

    printf("Small example graph first projection:\n");
    print_weighted_graph(&p1, &proj_weights1);
    printf("\n");
    printf("Small example graph second projection:\n");
    print_weighted_graph(&p2, &proj_weights2);
    printf("\n");

    igraph_vector_destroy(&proj_weights1);
    igraph_vector_destroy(&proj_weights2);
    igraph_destroy(&p1);
    igraph_destroy(&p2);

    igraph_destroy(&g);
    igraph_vector_bool_destroy(&types);
    igraph_vector_destroy(&weights);

    VERIFY_FINALLY_STACK();
    return 0;
}
