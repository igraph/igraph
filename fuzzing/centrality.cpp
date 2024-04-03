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

    igraph_set_warning_handler(igraph_warning_handler_ignore);

    if (Size % 2 == 0 || Size > 512+1 || Size < 1) {
        return 0;
    }

    igraph_vector_int_init(&edges, Size-1);
    for (size_t i=0; i < Size-1; ++i) {
        VECTOR(edges)[i] = Data[i+1];
    }

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

            igraph_betweenness_cutoff(&graph, &v, igraph_vss_all(), IGRAPH_ALL, NULL, 4);
            igraph_betweenness_cutoff(&graph, &v, igraph_vss_all(), IGRAPH_IN, NULL, 5);
            igraph_edge_betweenness_cutoff(&graph, &v, IGRAPH_DIRECTED, NULL, 4);
            igraph_edge_betweenness_cutoff(&graph, &v, IGRAPH_UNDIRECTED, NULL, 3);
            if (igraph_vcount(&graph) >= 10) {
                igraph_betweenness_subset(&graph, &v, igraph_vss_all(), IGRAPH_DIRECTED, igraph_vss_range(0,5), igraph_vss_range(5,10), NULL);
                igraph_edge_betweenness_subset(&graph, &v, igraph_ess_all(IGRAPH_EDGEORDER_ID), IGRAPH_DIRECTED, igraph_vss_range(0,10), igraph_vss_range(0,10), NULL);
            }
            igraph_closeness_cutoff(&graph, &v, &iv, &b, igraph_vss_all(), IGRAPH_ALL, NULL, true, 3);
            igraph_closeness_cutoff(&graph, &v, &iv, &b, igraph_vss_all(), IGRAPH_OUT, NULL, true, 4);
            igraph_harmonic_centrality_cutoff(&graph, &v, igraph_vss_all(), IGRAPH_ALL, NULL, true, 3);
            igraph_harmonic_centrality_cutoff(&graph, &v, igraph_vss_all(), IGRAPH_IN, NULL, true, 4);
            igraph_global_efficiency(&graph, &r, NULL, IGRAPH_DIRECTED);
            igraph_local_efficiency(&graph, &v, igraph_vss_all(), NULL, IGRAPH_DIRECTED, IGRAPH_OUT);
            igraph_transitivity_undirected(&graph, &r, IGRAPH_TRANSITIVITY_NAN);
            igraph_transitivity_local_undirected(&graph, &v, igraph_vss_all(), IGRAPH_TRANSITIVITY_NAN);
            igraph_transitivity_avglocal_undirected(&graph, &r, IGRAPH_TRANSITIVITY_ZERO);
            igraph_adjacent_triangles(&graph, &v, igraph_vss_all());
            igraph_pagerank(&graph, IGRAPH_PAGERANK_ALGO_PRPACK, &v, &r, igraph_vss_all(), IGRAPH_DIRECTED, 0.6, NULL, NULL);
            igraph_constraint(&graph, &v, igraph_vss_all(), NULL);
            igraph_spanner(&graph, &iv, 2.34, NULL);

            igraph_to_undirected(&graph, IGRAPH_TO_UNDIRECTED_COLLAPSE, NULL);
            igraph_simplify(&graph, /* multiple */ true, /* loops */ false, NULL);
            igraph_trussness(&graph, &iv);

            igraph_vector_int_destroy(&iv);
            igraph_vector_destroy(&v);
        }

        igraph_destroy(&graph);
    }

    igraph_vector_int_destroy(&edges);

    IGRAPH_ASSERT(IGRAPH_FINALLY_STACK_EMPTY);

    return 0;  // Non-zero return values are reserved for future use.
}
