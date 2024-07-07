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
        // Error at src/centrality/betweenness.c:437 : Weight vector must be positive. - Invalid value.
        VECTOR(weights)[i] = ((double) Data[i * 3 + 3] + 1.0) / 105.0;
    }

    // Turn on attribute handling. Weights will be stored as an edge attribute
    // in order to allow graph simplification while retainingweights.
    igraph_set_attribute_table(&igraph_cattribute_table);

    igraph_rng_seed(igraph_rng_default(), 42);

    if (igraph_create(&graph, &edges, Data[0], IGRAPH_DIRECTED) == IGRAPH_SUCCESS) {
        igraph_vector_t v;
        igraph_vector_int_t iv;
        igraph_bool_t b;
        igraph_real_t r;

        /* Limit graph size for the sake of performance. */
        if (igraph_vcount(&graph) <= 64) {
            igraph_vector_init(&v, 0);
            igraph_vector_int_init(&iv, 0);

            igraph_betweenness_cutoff(&graph, &v, igraph_vss_all(), IGRAPH_ALL, &weights, 4);
            igraph_betweenness_cutoff(&graph, &v, igraph_vss_all(), IGRAPH_IN, &weights, 5);
            igraph_edge_betweenness_cutoff(&graph, &v, IGRAPH_DIRECTED, &weights, 4);
            igraph_edge_betweenness_cutoff(&graph, &v, IGRAPH_UNDIRECTED, &weights, 3);
            if (igraph_vcount(&graph) >= 10) {
                igraph_betweenness_subset(&graph, &v, igraph_vss_all(), IGRAPH_DIRECTED, igraph_vss_range(0,5), igraph_vss_range(5,10), &weights);
                igraph_edge_betweenness_subset(&graph, &v, igraph_ess_all(IGRAPH_EDGEORDER_ID), IGRAPH_DIRECTED, igraph_vss_range(0,10), igraph_vss_range(0,10), &weights);
            }
            igraph_closeness_cutoff(&graph, &v, &iv, &b, igraph_vss_all(), IGRAPH_ALL, &weights, true, 3);
            igraph_closeness_cutoff(&graph, &v, &iv, &b, igraph_vss_all(), IGRAPH_OUT, &weights, true, 4);
            igraph_harmonic_centrality_cutoff(&graph, &v, igraph_vss_all(), IGRAPH_ALL, &weights, true, 3);
            igraph_harmonic_centrality_cutoff(&graph, &v, igraph_vss_all(), IGRAPH_IN, &weights, true, 4);
            igraph_global_efficiency(&graph, &r, &weights, IGRAPH_DIRECTED);
            igraph_local_efficiency(&graph, &v, igraph_vss_all(), &weights, IGRAPH_DIRECTED, IGRAPH_OUT);
            igraph_pagerank(&graph, IGRAPH_PAGERANK_ALGO_PRPACK, &v, &r, igraph_vss_all(), IGRAPH_DIRECTED, 0.6, &weights, NULL);
            igraph_constraint(&graph, &v, igraph_vss_all(), &weights);
            igraph_spanner(&graph, &iv, 2.34, &weights);

            {
                igraph_attribute_combination_t comb;
                SETEANV(&graph, "weight", &weights);
                igraph_attribute_combination(&comb,
                                             "weight", IGRAPH_ATTRIBUTE_COMBINE_SUM,
                                             IGRAPH_NO_MORE_ATTRIBUTES);
                // This operation would be simpler if we use IGRAPH_TO_UNDIRECTED_EACH,
                // as collapsing edges happens in the simplification step anyway.
                // We use IGRAPH_TO_UNDIRECTED_COLLAPSE to exercise more code.
                igraph_to_undirected(&graph, IGRAPH_TO_UNDIRECTED_COLLAPSE, &comb);
                igraph_simplify(&graph, true, true, &comb);
                igraph_attribute_combination_destroy(&comb);
                EANV(&graph, "weight", &weights);
                DELEA(&graph, "weight");
            }

            igraph_diversity(&graph, &weights, &v, igraph_vss_all());

            igraph_vector_int_destroy(&iv);
            igraph_vector_destroy(&v);
        }

        igraph_destroy(&graph);
    }

    igraph_vector_int_destroy(&edges);
    igraph_vector_destroy(&weights);

    IGRAPH_ASSERT(IGRAPH_FINALLY_STACK_EMPTY);

    return 0;  // Non-zero return values are reserved for future use.
}
