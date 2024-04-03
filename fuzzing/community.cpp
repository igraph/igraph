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
        igraph_matrix_int_t merges, im;
        igraph_matrix_t mat;
        igraph_vector_int_t membership, membership2, iv, iv2;
        igraph_vector_t mv, v;
        igraph_real_t m, r;
        igraph_integer_t i;
        igraph_bool_t b;

        /* Limit graph size for the sake of performance. */
        if (igraph_vcount(&graph) <= 64) {
            igraph_matrix_int_init(&merges, 0, 0);
            igraph_matrix_int_init(&im, 0, 0);
            igraph_matrix_init(&mat, 0, 0);
            igraph_vector_int_init(&membership, 0);
            igraph_vector_int_init(&membership2, 0);
            igraph_vector_int_init(&iv, 0);
            igraph_vector_int_init(&iv2, 0);
            igraph_vector_init(&mv, 0);
            igraph_vector_init(&v, 0);

            igraph_community_label_propagation(&graph, &membership, IGRAPH_OUT, NULL, NULL, NULL);
            igraph_community_walktrap(&graph, NULL, 3, &merges, &mv, &membership);
            igraph_community_edge_betweenness(&graph, &iv, &v, &merges, &iv2, &mv, &membership2, IGRAPH_DIRECTED, NULL);

            // Take the opportunity to run functions that can use the output of community detection.

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
                igraph_contract_vertices(&graph2, &membership, NULL);
                igraph_destroy(&graph2);
            }

            igraph_modularity(&graph, &membership, NULL, 1.5, IGRAPH_UNDIRECTED, &m);
            igraph_modularity_matrix(&graph, NULL, 0.75, &mat, IGRAPH_DIRECTED);
            igraph_assortativity_nominal(&graph, &membership, &r, IGRAPH_DIRECTED, true);

            igraph_simplify(&graph, true, true, NULL);
            igraph_community_voronoi(&graph, &membership, &iv, &m, NULL, NULL, IGRAPH_OUT, 1.0);

            igraph_to_undirected(&graph, IGRAPH_TO_UNDIRECTED_COLLAPSE, NULL);
            igraph_community_fastgreedy(&graph, NULL, &merges, &mv, &membership);
            igraph_community_leiden(&graph, NULL, NULL, 1.5, 0.01, false, 2, &membership, &i, &r);
            igraph_community_multilevel(&graph, NULL, 0.8, &membership, &im, &mv);

            igraph_is_connected(&graph, &b, IGRAPH_WEAK);
            if (b) {
                igraph_integer_t no_comm = 4;
                if (no_comm > igraph_vcount(&graph)) {
                    no_comm = igraph_vcount(&graph);
                }
                igraph_community_fluid_communities(&graph, no_comm, &membership);
            }

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

    IGRAPH_ASSERT(IGRAPH_FINALLY_STACK_EMPTY);

    return 0;  // Non-zero return values are reserved for future use.
}
