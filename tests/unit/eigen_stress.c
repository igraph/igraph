/*
   igraph library.
   Copyright (C) 2024-2026  The igraph development team <igraph@igraph.org>

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
#include <stdio.h>

#include "test_utilities.h"

#define CHECK_PRINT(expr) \
    do { \
        igraph_error_t err = expr; \
        if (err != IGRAPH_SUCCESS) { \
            fprintf(stderr, "Failed computation on a %s graph with %" IGRAPH_PRId \
                    " vertices and %" IGRAPH_PRId " edges. Edge list:\n", \
                    igraph_is_directed(&g) ? "directed" : "undirected", \
                    igraph_vcount(&g), igraph_ecount(&g)); \
            igraph_write_graph_edgelist(&g, stderr); \
            IGRAPH_FATALF("'%s' failed with error %d (%s).", \
                          #expr, (int) err, igraph_strerror(err)); \
        } \
    } while (0)

int main(void) {
    igraph_vector_t vec;
    igraph_real_t val;

    igraph_rng_seed(igraph_rng_default(), 987);
    igraph_set_warning_handler(&igraph_warning_handler_ignore);
    igraph_set_error_handler(&igraph_error_handler_printignore);

    igraph_vector_init(&vec, 0);

    for (igraph_int_t i = 0; i < 1252; i++) {
        igraph_t g;
        igraph_atlas(&g, i);

        CHECK_PRINT(igraph_eigenvector_centrality(&g, &vec, &val, IGRAPH_ALL, NULL, NULL));
        CHECK_PRINT(igraph_hub_and_authority_scores(&g, &vec, NULL, &val, NULL, NULL));

        igraph_to_directed(&g, IGRAPH_TO_DIRECTED_RANDOM);
        CHECK_PRINT(igraph_eigenvector_centrality(&g, &vec, &val, IGRAPH_OUT, NULL, NULL));
        CHECK_PRINT(igraph_eigenvector_centrality(&g, &vec, &val, IGRAPH_IN, NULL, NULL));
        CHECK_PRINT(igraph_hub_and_authority_scores(&g, &vec, NULL, &val, NULL, NULL));
        CHECK_PRINT(igraph_pagerank(&g, NULL, &vec, &val, 0.85, IGRAPH_DIRECTED,
                                    igraph_vss_all(), IGRAPH_PAGERANK_ALGO_ARPACK, NULL));
        CHECK_PRINT(igraph_pagerank(&g, NULL, &vec, &val, 0.99, IGRAPH_DIRECTED,
                                    igraph_vss_all(), IGRAPH_PAGERANK_ALGO_ARPACK, NULL));

        igraph_destroy(&g);
    }

    for (igraph_int_t i = 3; i < 100; i++) {
        igraph_t g;
        igraph_ring(&g, i, IGRAPH_DIRECTED, /*mutual=*/ false, /*circular=*/ true);
        CHECK_PRINT(igraph_eigenvector_centrality(&g, &vec, NULL, IGRAPH_OUT, NULL, NULL));
        igraph_destroy(&g);
    }

    {
        /* A <--> C, D <--> E, B */
        igraph_t g;
        igraph_small(&g, 5, IGRAPH_DIRECTED, 0, 2, 2, 0, 3, 4, 4, 2, -1);

        CHECK_PRINT(igraph_eigenvector_centrality(&g, &vec, &val, IGRAPH_OUT, NULL, NULL));
        CHECK_PRINT(igraph_eigenvector_centrality(&g, &vec, &val, IGRAPH_IN, NULL, NULL));
        CHECK_PRINT(igraph_hub_and_authority_scores(&g, &vec, NULL, &val, NULL, NULL));
        CHECK_PRINT(igraph_pagerank(&g, NULL, &vec, &val, 0.85, IGRAPH_DIRECTED,
                                    igraph_vss_all(), IGRAPH_PAGERANK_ALGO_ARPACK, NULL));
        CHECK_PRINT(igraph_pagerank(&g, NULL, &vec, &val, 0.99, IGRAPH_DIRECTED,
                                    igraph_vss_all(), IGRAPH_PAGERANK_ALGO_ARPACK, NULL));

        igraph_destroy(&g);
    }

    for (igraph_int_t n = 4; n <= 7; n++) {
        for (igraph_int_t i = 0; i < 100; i++) {
            igraph_t g;

            igraph_erdos_renyi_game_gnp(&g, n, 0.5, IGRAPH_DIRECTED, IGRAPH_LOOPS_SW, IGRAPH_EDGE_UNLABELED);

            CHECK_PRINT(igraph_eigenvector_centrality(&g, &vec, &val, IGRAPH_OUT, NULL, NULL));
            CHECK_PRINT(igraph_eigenvector_centrality(&g, &vec, &val, IGRAPH_IN, NULL, NULL));
            CHECK_PRINT(igraph_hub_and_authority_scores(&g, &vec, NULL, &val, NULL, NULL));
            CHECK_PRINT(igraph_pagerank(&g, NULL, &vec, &val, 0.85, IGRAPH_DIRECTED,
                                        igraph_vss_all(), IGRAPH_PAGERANK_ALGO_ARPACK, NULL));
            CHECK_PRINT(igraph_pagerank(&g, NULL, &vec, &val, 0.99, IGRAPH_DIRECTED,
                                        igraph_vss_all(), IGRAPH_PAGERANK_ALGO_ARPACK, NULL));

            igraph_destroy(&g);
        }
    }

    igraph_vector_destroy(&vec);

    VERIFY_FINALLY_STACK();

    return 0;
}
