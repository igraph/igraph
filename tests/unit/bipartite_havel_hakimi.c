//
// Created by Lara Holm on 1.11.2023.
//
#include <igraph.h>
#include <stdio.h>

#include "test_utilities.h"
#include "math/safe_intop.h"

int main(void) {
    igraph_vector_int_t ds1, ds2;
    igraph_t g;
    igraph_bool_t is_connected;
    igraph_bool_t is_bipartite;
    igraph_bool_t is_simple;


    // Simple undirected bipartite graph.
    igraph_integer_t deg1[] = {2, 3, 2, 1};
    igraph_integer_t deg2[] = {3, 1, 2, 1, 1};

    igraph_vector_int_init_array(&ds1, deg1, 4);
    igraph_vector_int_init_array(&ds2, deg2, 5);

//    igraph_vector_int_print(&ds1);
//    igraph_vector_int_print(&ds2);

    igraph_realize_bipartite_degree_sequence(&g, &ds1, &ds2, false);

    igraph_is_simple(&g, &is_simple);
    igraph_is_connected(&g, &is_connected, IGRAPH_STRONG);
    igraph_is_bipartite(&g, &is_bipartite, NULL);
    // TODO: Check the degrees of the graph - are they the same
    // Verify that there is a graph, and that it is simple, bipartite and connected if it was potentially connected
    // igraph_has_multiple - returns a bool if it has a multiedge or igraph_is_simple for both self loops and multiedges
    // helper function to check if the initial dss are potentially connected for the tests
    // should we get simple graph
    // random dss

    printf("Simple: %d\n", is_simple);
    printf("Connected: %d\n", is_connected);
    printf("Bipartite: %d\n", is_bipartite);

    igraph_vector_int_destroy(&ds1);
    igraph_vector_int_destroy(&ds2);
    igraph_destroy(&g);


    // undirected bipartite multigraph.
    igraph_bool_t has_multi;

    igraph_integer_t deg1m[] = {2, 3, 1};
    igraph_integer_t deg2m[]= {4, 2};

    igraph_vector_int_init_array(&ds1, deg1m, 3);
    igraph_vector_int_init_array(&ds2, deg2m, 2);

    igraph_realize_bipartite_degree_sequence(&g, &ds1, &ds2, true);

    igraph_vector_int_t edges;
    igraph_vector_int_init(&edges, 0);

    igraph_get_edgelist(&g, &edges, false);

    igraph_vector_int_print(&edges);

    igraph_has_multiple(&g, &has_multi);
    // Unconnected?
    igraph_is_connected(&g, &is_connected, IGRAPH_STRONG);
    igraph_is_bipartite(&g, &is_bipartite, NULL);

    printf("Has multiedges: %d\n", has_multi);
    printf("Connected: %d\n", is_connected);
    printf("Bipartite: %d\n", is_bipartite);

    igraph_vector_int_destroy(&edges);
    igraph_vector_int_destroy(&ds1);
    igraph_vector_int_destroy(&ds2);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();
    return 0;
}
