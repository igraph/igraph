/*
  igraph library.
  Constructing realizations of degree sequences and bi-degree sequences.
  Copyright (C) 2024  The igraph development team <igraph@igraph.org>

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

void check_degrees(
    const igraph_t *graph,
    const igraph_vector_int_t *ds1, const igraph_vector_int_t *ds2
) {
    const igraph_int_t vcount = igraph_vcount(graph);
    const igraph_int_t n1 = igraph_vector_int_size(ds1);
    const igraph_int_t n2 = igraph_vector_int_size(ds2);
    igraph_vector_int_t degrees, expected_degrees;

    if (vcount != n1+n2) {
        printf("Bad graph size %" IGRAPH_PRId
               " for bidegree sequence size (%" IGRAPH_PRId
               ", %" IGRAPH_PRId ")",
               vcount, n1, n2);
        return;
    }

    igraph_vector_int_init(&degrees, vcount);
    igraph_degree(graph, &degrees, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS);

    igraph_vector_int_init_copy(&expected_degrees, ds1);
    igraph_vector_int_append(&expected_degrees, ds2);
    if (! igraph_vector_int_all_e(&degrees, &expected_degrees)) {
        printf("Degrees not as expected!");
        printf("Expected: "); print_vector_int(&expected_degrees);
        printf("Actual:   "); print_vector_int(&degrees);
    }

    igraph_vector_int_destroy(&expected_degrees);
    igraph_vector_int_destroy(&degrees);
}

int main(void) {
    igraph_vector_int_t ds1, ds2;
    igraph_t g;
    igraph_bool_t is_connected;
    igraph_bool_t is_bipartite;
    igraph_bool_t is_simple;

    igraph_realize_degseq_t methods[] = {
        IGRAPH_REALIZE_DEGSEQ_SMALLEST,
        IGRAPH_REALIZE_DEGSEQ_LARGEST,
        IGRAPH_REALIZE_DEGSEQ_INDEX
    };
    char *method_names[] = {
        "smallest", "largest", "index"
    };

    // Edge cases
    printf("===EDGE CASES===\n");

    // Tests that should fail with IGRAPH_EINVAL
    printf("--ds1 empty--\n");
    igraph_vector_int_init(&ds1, 0);
    igraph_vector_int_init_int(&ds2, 2,
                               1, 1);
    CHECK_ERROR(
        igraph_realize_bipartite_degree_sequence(&g, &ds1, &ds2, IGRAPH_SIMPLE_SW, IGRAPH_REALIZE_DEGSEQ_SMALLEST),
        IGRAPH_EINVAL
    );

    printf("--ds2 empty--\n");
    CHECK_ERROR(
        igraph_realize_bipartite_degree_sequence(&g, &ds2, &ds1, IGRAPH_SIMPLE_SW, IGRAPH_REALIZE_DEGSEQ_SMALLEST),
        IGRAPH_EINVAL
    );
    igraph_vector_int_destroy(&ds1);
    igraph_vector_int_destroy(&ds2);

    printf("\n===Empty degree sequences===\n");

    igraph_vector_int_init(&ds1, 0);
    igraph_vector_int_init(&ds2, 0);
    for (size_t i=0; i < sizeof(methods) / sizeof(methods[0]); i++) {
        printf("Method: %s\n", method_names[i]);
        igraph_realize_bipartite_degree_sequence(
            &g, &ds1, &ds2, IGRAPH_SIMPLE_SW, methods[i]);
        print_graph(&g);
        printf("\n");
        igraph_destroy(&g);
    }
    igraph_vector_int_destroy(&ds1);
    igraph_vector_int_destroy(&ds2);

    printf("\n===All 0 degree sequences===\n");

    igraph_vector_int_init(&ds1, 3);
    igraph_vector_int_init(&ds2, 4);
    for (size_t i=0; i < sizeof(methods) / sizeof(methods[0]); i++) {
        printf("Method: %s\n", method_names[i]);
        igraph_realize_bipartite_degree_sequence(
            &g, &ds1, &ds2, IGRAPH_SIMPLE_SW, methods[i]);
        print_graph(&g);
        printf("\n");
        igraph_destroy(&g);
    }
    igraph_vector_int_destroy(&ds1);
    igraph_vector_int_destroy(&ds2);

    printf("\n===TESTING SMALLEST METHOD===\n");
    // Simple undirected bipartite graph.
    printf("\n===Simple, connected graph===\n");
    igraph_int_t deg1[] = {2, 3, 2, 1};
    igraph_int_t deg2[] = {3, 1, 2, 1, 1};

    igraph_vector_int_init_array(&ds1, deg1, 4);
    igraph_vector_int_init_array(&ds2, deg2, 5);

    igraph_realize_bipartite_degree_sequence(&g, &ds1, &ds2, IGRAPH_SIMPLE_SW, IGRAPH_REALIZE_DEGSEQ_SMALLEST);

    igraph_is_simple(&g, &is_simple, IGRAPH_DIRECTED);
    igraph_is_connected(&g, &is_connected, IGRAPH_STRONG);
    igraph_is_bipartite(&g, &is_bipartite, NULL);

    print_graph(&g);
    check_degrees(&g, &ds1, &ds2);

    printf("Simple: %d\n", is_simple);
    printf("Connected: %d\n", is_connected);
    printf("Bipartite: %d\n", is_bipartite);

    igraph_vector_int_destroy(&ds1);
    igraph_vector_int_destroy(&ds2);
    igraph_destroy(&g);

    // undirected bipartite multigraph.
    igraph_bool_t has_multi;
    printf("\n===Connected multigraph===\n");

    igraph_int_t deg1m[] = {2, 3, 1};
    igraph_int_t deg2m[] = {4, 2};

    igraph_vector_int_init_array(&ds1, deg1m, 3);
    igraph_vector_int_init_array(&ds2, deg2m, 2);

    igraph_realize_bipartite_degree_sequence(&g, &ds1, &ds2, IGRAPH_MULTI_SW, IGRAPH_REALIZE_DEGSEQ_SMALLEST);

    igraph_has_multiple(&g, &has_multi);
    igraph_is_connected(&g, &is_connected, IGRAPH_STRONG);
    igraph_is_bipartite(&g, &is_bipartite, NULL);

    print_graph(&g);
    check_degrees(&g, &ds1, &ds2);

    printf("vcount: %" IGRAPH_PRId "\n", igraph_vcount(&g));
    printf("Has multiedges: %d\n", has_multi);
    printf("Connected: %d\n", is_connected);
    printf("Bipartite: %d\n", is_bipartite);

    igraph_vector_int_destroy(&ds1);
    igraph_vector_int_destroy(&ds2);
    igraph_destroy(&g);

    printf("\n===TESTING LARGEST METHOD===\n");
    // Simple undirected bipartite graph.
    printf("\n===Simple graph===\n");
    igraph_int_t deg1l[] = {2, 3, 2, 1};
    igraph_int_t deg2l[] = {3, 1, 2, 1, 1};

    igraph_vector_int_init_array(&ds1, deg1l, 4);
    igraph_vector_int_init_array(&ds2, deg2l, 5);

    igraph_realize_bipartite_degree_sequence(&g, &ds1, &ds2, IGRAPH_SIMPLE_SW, IGRAPH_REALIZE_DEGSEQ_LARGEST);

    igraph_is_simple(&g, &is_simple, IGRAPH_DIRECTED);
    // For this method, it does not have to be connected
    igraph_is_connected(&g, &is_connected, IGRAPH_WEAK);
    igraph_is_bipartite(&g, &is_bipartite, NULL);

    print_graph(&g);
    check_degrees(&g, &ds1, &ds2);

    printf("Simple: %d\n", is_simple);
    printf("Connected: %d\n", is_connected);
    printf("Bipartite: %d\n", is_bipartite);

    igraph_vector_int_destroy(&ds1);
    igraph_vector_int_destroy(&ds2);
    igraph_destroy(&g);

    // undirected bipartite multigraph.
    printf("\n===Multigraph===\n");

    igraph_int_t deg1l_m[] = {2, 3, 1};
    igraph_int_t deg2l_m[]= {4, 2};

    igraph_vector_int_init_array(&ds1, deg1l_m, 3);
    igraph_vector_int_init_array(&ds2, deg2l_m, 2);

    igraph_realize_bipartite_degree_sequence(&g, &ds1, &ds2, IGRAPH_MULTI_SW, IGRAPH_REALIZE_DEGSEQ_LARGEST);

    igraph_has_multiple(&g, &has_multi);
    // For this method, it does not need to be connected
    igraph_is_connected(&g, &is_connected, IGRAPH_WEAK);
    igraph_is_bipartite(&g, &is_bipartite, NULL);

    print_graph(&g);
    check_degrees(&g, &ds1, &ds2);

    printf("Has multiedges: %d\n", has_multi);
    printf("Connected: %d\n", is_connected);
    printf("Bipartite: %d\n", is_bipartite);

    igraph_vector_int_destroy(&ds1);
    igraph_vector_int_destroy(&ds2);
    igraph_destroy(&g);

    printf("\n===TESTING INDEX METHOD===\n");
    // Simple, undirected bipartite.
    printf("\n===Simple graph===\n");
    igraph_int_t deg1i[] = {1, 4, 3, 2};
    igraph_int_t deg2i[] = {2, 1, 1, 4, 2};

    igraph_vector_int_init_array(&ds1, deg1i, 4);
    igraph_vector_int_init_array(&ds2, deg2i, 5);

    igraph_realize_bipartite_degree_sequence(&g, &ds1, &ds2, IGRAPH_SIMPLE_SW, IGRAPH_REALIZE_DEGSEQ_INDEX);

    igraph_is_simple(&g, &is_simple, IGRAPH_DIRECTED);
    // For this method, it does not have to be connected
    igraph_is_connected(&g, &is_connected, IGRAPH_STRONG);
    igraph_is_bipartite(&g, &is_bipartite, NULL);

    print_graph(&g);
    check_degrees(&g, &ds1, &ds2);

    printf("Simple: %d\n", is_simple);
    printf("Connected: %d\n", is_connected);
    printf("Bipartite: %d\n", is_bipartite);

    igraph_vector_int_destroy(&ds1);
    igraph_vector_int_destroy(&ds2);
    igraph_destroy(&g);

    // undirected bipartite multigraph.
    printf("\n===Multigraph===\n");

    igraph_int_t deg1i_m[] = {2, 3, 1};
    igraph_int_t deg2i_m[]= {4, 2};

    igraph_vector_int_init_array(&ds1, deg1i_m, 3);
    igraph_vector_int_init_array(&ds2, deg2i_m, 2);

    igraph_realize_bipartite_degree_sequence(&g, &ds1, &ds2, IGRAPH_MULTI_SW, IGRAPH_REALIZE_DEGSEQ_INDEX);

    print_graph(&g);
    check_degrees(&g, &ds1, &ds2);

    igraph_has_multiple(&g, &has_multi);
    // For this method, it does not need to be connected
    igraph_is_connected(&g, &is_connected, IGRAPH_STRONG);
    igraph_is_bipartite(&g, &is_bipartite, NULL);

    printf("Has multiedges: %d\n", has_multi);
    printf("Connected: %d\n", is_connected);
    printf("Bipartite: %d\n", is_bipartite);

    igraph_vector_int_destroy(&ds1);
    igraph_vector_int_destroy(&ds2);
    igraph_destroy(&g);

    // A bidegree sequence that is not potentially connected
    igraph_vector_int_init_int(&ds1, 3,
                               1, 1, 1);
    for (size_t i=0; i < sizeof(methods) / sizeof(methods[0]); i++) {
        igraph_realize_bipartite_degree_sequence(
            &g, &ds1, &ds1, IGRAPH_SIMPLE_SW, methods[i]);
        igraph_is_bipartite(&g, &is_bipartite, NULL);
        IGRAPH_ASSERT(is_bipartite);
        igraph_destroy(&g);
    }
    igraph_vector_int_destroy(&ds1);

    // A bidegree sequence that can only be realized as a multigraph,
    // but not a simple graph
    igraph_vector_int_init_int(&ds1, 3,
                               1, 3, 1);
    igraph_vector_int_init_int(&ds2, 2,
                               2, 3);
    for (size_t i=0; i < sizeof(methods) / sizeof(methods[0]); i++) {
        CHECK_ERROR(
            igraph_realize_bipartite_degree_sequence(
            &g, &ds1, &ds2, IGRAPH_SIMPLE_SW, methods[i]),
            IGRAPH_EINVAL
        );

        igraph_realize_bipartite_degree_sequence(
            &g, &ds1, &ds2, IGRAPH_MULTI_SW, methods[i]);

        igraph_is_bipartite(&g, &is_bipartite, NULL);
        IGRAPH_ASSERT(is_bipartite);
        igraph_destroy(&g);
    }
    igraph_vector_int_destroy(&ds1);
    igraph_vector_int_destroy(&ds2);

    VERIFY_FINALLY_STACK();
    return 0;
}
