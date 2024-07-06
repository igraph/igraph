/*
   IGraph library.
   Copyright (C) 2024  The igraph development team

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
#include <cstdlib>

extern "C" int LLVMFuzzerTestOneInput(const uint8_t *Data, size_t Size) {
    igraph_t graph;
    igraph_vector_int_t edges;
    igraph_vector_t weights;

    igraph_set_warning_handler(igraph_warning_handler_ignore);

    if (Size % 3 == 0 || Size % 3 == 2 || Size > (3 * 256) + 1 || Size < 1) {
        return 0;
    }

    igraph_vector_int_init(&edges, ((Size-1) / 3) * 2);
    igraph_vector_init(&weights, (Size-1) / 3);
    for (size_t i=0; i < ((Size-1) / 3); ++i) {
        VECTOR(edges)[i * 2] = Data[i * 3 + 1];
        VECTOR(edges)[i * 2 + 1] = Data[i * 3 + 2];
        // We keep the weights strictly positive, as this is required by some algorithms.
        VECTOR(weights)[i] = ((double) Data[i * 3 + 3] + 1.0) / 105.0;
    }

    // Turn on attribute handling.
    igraph_set_attribute_table(&igraph_cattribute_table);

    igraph_rng_seed(igraph_rng_default(), 42);

    if (igraph_create(&graph, &edges, Data[0], IGRAPH_DIRECTED) == IGRAPH_SUCCESS) {
        igraph_matrix_int_t merges, im;
        igraph_matrix_t mat;
        igraph_vector_int_t membership, membership2, iv, iv2;
        igraph_vector_t mv, v;
        igraph_real_t m, r;
        igraph_integer_t i;
        igraph_bool_t b;

        /* Limit graph size for the sake of performance. */
        if (igraph_vcount(&graph) <= 64) {
            igraph_attribute_combination_t comb;

            igraph_matrix_int_init(&merges, 0, 0);
            igraph_matrix_int_init(&im, 0, 0);
            igraph_matrix_init(&mat, 0, 0);
            igraph_vector_int_init(&membership, 0);
            igraph_vector_int_init(&membership2, 0);
            igraph_vector_int_init(&iv, 0);
            igraph_vector_int_init(&iv2, 0);
            igraph_vector_init(&mv, 0);
            igraph_vector_init(&v, 0);

            SETEANV(&graph, "weight", &weights);

            // Currently "strength" is used only used as a vertex attribute to exercise the
            // attribute handling code durig vertex contraction, which is convenient to do
            // using the output from community detection. It is not used as input to any
            // community detection methods.
            igraph_strength(&graph, &v, igraph_vss_all(), IGRAPH_OUT, IGRAPH_LOOPS, &weights);
            SETVANV(&graph, "strength", &v);

            igraph_attribute_combination(&comb,
                                         "weight", IGRAPH_ATTRIBUTE_COMBINE_SUM,
                                         "strength", IGRAPH_ATTRIBUTE_COMBINE_MAX,
                                         IGRAPH_NO_MORE_ATTRIBUTES);

            // Ignore edge directions in label propagation for now, as on some weighted graphs the
            // algorithm seems to never complete. See https://github.com/igraph/igraph/issues/2561
            igraph_community_label_propagation(&graph, &membership, IGRAPH_ALL, &weights, NULL, NULL);

            igraph_community_walktrap(&graph, &weights, 3, &merges, &mv, &membership);
            igraph_community_edge_betweenness(&graph, &iv, &v, &merges, &iv2, &mv, &membership2, IGRAPH_DIRECTED, &weights);

            // Take the opportunity to run functions that can use the output of community detection,
            // potentially with weights.

            {
                igraph_community_comparison_t method[] = {
                    IGRAPH_COMMCMP_VI,
                    IGRAPH_COMMCMP_NMI,
                    IGRAPH_COMMCMP_SPLIT_JOIN,
                    IGRAPH_COMMCMP_RAND,
                    IGRAPH_COMMCMP_ADJUSTED_RAND
                };
                for (size_t i=0; i < sizeof(method) / sizeof(method[0]); i++) {
                    if ((method[i] == IGRAPH_COMMCMP_RAND ||
                         method[i] == IGRAPH_COMMCMP_ADJUSTED_RAND) &&
                        igraph_vcount(&graph) < 2) {
                        continue;
                    }
                    igraph_compare_communities(&membership, &membership2, &r, method[i]);
                }
            }

            {
                igraph_t graph2;
                igraph_copy(&graph2, &graph);
                igraph_contract_vertices(&graph2, &membership, &comb);
                igraph_destroy(&graph2);
            }

            igraph_modularity(&graph, &membership, &weights, 1.5, IGRAPH_UNDIRECTED, &m);
            igraph_modularity_matrix(&graph, &weights, 0.75, &mat, IGRAPH_DIRECTED);

            // NOTE: Currently Infomap dominates this fuzz target due to being
            // slower and more complex than the other algorithms.
            // TODO: Currently disabled because of extremely bad performance.
            // igraph_community_infomap(&graph, &weights, NULL, 1, &membership, &r);

            igraph_simplify(&graph, true, true, &comb);
            EANV(&graph, "weight", &weights);

            // Compute 'length' vector for community_voronoi()
            igraph_vector_update(&v, &weights);
            for (igraph_integer_t i=0; i < igraph_vector_size(&v); i++) {
                VECTOR(v)[i] = 1 / (1 + VECTOR(v)[i]);
            }
            igraph_community_voronoi(&graph, &membership, &iv, &m, &v, &weights, IGRAPH_OUT, 1.0);

            igraph_to_undirected(&graph, IGRAPH_TO_UNDIRECTED_COLLAPSE, &comb);
            EANV(&graph, "weight", &weights);

            igraph_community_fastgreedy(&graph, &weights, &merges, &mv, &membership);
            igraph_community_leiden(&graph, &weights, NULL, 1.5, 0.01, false, 2, &membership, &i, &r);
            igraph_community_multilevel(&graph, &weights, 0.8, &membership, &im, &mv);

            // community_spinglass() only works on connected graphs
            igraph_is_connected(&graph, &b, IGRAPH_WEAK);
            if (b) {
                igraph_community_spinglass(&graph, &weights, &r, NULL, &membership, NULL, 10, false, 1.0, 0.01, 0.99, IGRAPH_SPINCOMM_UPDATE_CONFIG, 1.0, IGRAPH_SPINCOMM_IMP_ORIG, 1.0);
            }

            DELALL(&graph);

            igraph_attribute_combination_destroy(&comb);

            igraph_vector_destroy(&v);
            igraph_vector_destroy(&mv);
            igraph_vector_int_destroy(&iv2);
            igraph_vector_int_destroy(&iv);
            igraph_vector_int_destroy(&membership2);
            igraph_vector_int_destroy(&membership);
            igraph_matrix_destroy(&mat);
            igraph_matrix_int_destroy(&im);
            igraph_matrix_int_destroy(&merges);
        }

        igraph_destroy(&graph);
    }

    igraph_vector_int_destroy(&edges);
    igraph_vector_destroy(&weights);

    IGRAPH_ASSERT(IGRAPH_FINALLY_STACK_EMPTY);

    return 0;  // Non-zero return values are reserved for future use.
}
