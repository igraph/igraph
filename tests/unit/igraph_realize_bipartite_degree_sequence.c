//
// Created by Lara Holm on 1.11.2023.
//
#include <igraph.h>
#include <stdio.h>

#include "test_utilities.h"
#include "math/safe_intop.h"

int main(void) {
    // TODO: Edge cases - empty dss, all 0 degrees,
    igraph_vector_int_t ds1, ds2;
    igraph_t g;
    igraph_bool_t is_connected;
    igraph_bool_t is_bipartite;
    igraph_bool_t is_simple;
    igraph_vector_int_t degrees;
    igraph_vector_int_t edges;

    printf("===TESTING SMALLEST METHOD===\n");
    // Simple undirected bipartite graph.
    printf("\n===Simple, connected graph===\n");
    igraph_integer_t deg1[] = {2, 3, 2, 1};
    igraph_integer_t deg2[] = {3, 1, 2, 1, 1};

    igraph_vector_int_init_array(&ds1, deg1, 4);
    igraph_vector_int_init_array(&ds2, deg2, 5);

    igraph_realize_bipartite_degree_sequence(&g, &ds1, &ds2, IGRAPH_SIMPLE_SW, IGRAPH_REALIZE_DEGSEQ_SMALLEST);

    igraph_is_simple(&g, &is_simple);
    igraph_is_connected(&g, &is_connected, IGRAPH_STRONG);
    igraph_is_bipartite(&g, &is_bipartite, NULL);

    printf("Edges: ");
    igraph_vector_int_init(&edges, 0);
    igraph_get_edgelist(&g, &edges, false);
    igraph_vector_int_print(&edges);

    printf("Degrees: ");
    igraph_vector_int_init(&degrees, 0);
    igraph_degree(&g, &degrees, igraph_vss_all(), IGRAPH_ALL, false);
    igraph_vector_int_print(&degrees);

    printf("vcount: %d\n", (int)igraph_vcount(&g));
    printf("Simple: %d\n", is_simple);
    printf("Connected: %d\n", is_connected);
    printf("Bipartite: %d\n", is_bipartite);

    igraph_vector_int_destroy(&degrees);
    igraph_vector_int_destroy(&edges);
    igraph_vector_int_destroy(&ds1);
    igraph_vector_int_destroy(&ds2);
    igraph_destroy(&g);

    // undirected bipartite multigraph.
    igraph_bool_t has_multi;
    printf("===Connected multigraph===\n");

    igraph_integer_t deg1m[] = {2, 3, 1};
    igraph_integer_t deg2m[]= {4, 2};

    igraph_vector_int_init_array(&ds1, deg1m, 3);
    igraph_vector_int_init_array(&ds2, deg2m, 2);

    igraph_realize_bipartite_degree_sequence(&g, &ds1, &ds2, IGRAPH_MULTI_SW, IGRAPH_REALIZE_DEGSEQ_SMALLEST);

    printf("Edges: ");
    igraph_vector_int_init(&edges, 0);
    igraph_get_edgelist(&g, &edges, false);
    igraph_vector_int_print(&edges);

    printf("Degrees: ");
    igraph_vector_int_init(&degrees, 0);
    igraph_degree(&g, &degrees, igraph_vss_all(), IGRAPH_ALL, false);
    igraph_vector_int_print(&degrees);

    igraph_has_multiple(&g, &has_multi);
    igraph_is_connected(&g, &is_connected, IGRAPH_STRONG);
    igraph_is_bipartite(&g, &is_bipartite, NULL);

    printf("vcount: %d\n", (int)igraph_vcount(&g));
    printf("Has multiedges: %d\n", has_multi);
    printf("Connected: %d\n", is_connected);
    printf("Bipartite: %d\n", is_bipartite);

    igraph_vector_int_destroy(&degrees);
    igraph_vector_int_destroy(&edges);
    igraph_vector_int_destroy(&ds1);
    igraph_vector_int_destroy(&ds2);
    igraph_destroy(&g);

    printf("\n===TESTING LARGEST METHOD===\n");
    // Simple undirected bipartite graph.
    printf("\n===Simple graph===\n");
    igraph_integer_t deg1l[] = {2, 3, 2, 1};
    igraph_integer_t deg2l[] = {3, 1, 2, 1, 1};

    igraph_vector_int_init_array(&ds1, deg1l, 4);
    igraph_vector_int_init_array(&ds2, deg2l, 5);

    igraph_realize_bipartite_degree_sequence(&g, &ds1, &ds2, IGRAPH_SIMPLE_SW, IGRAPH_REALIZE_DEGSEQ_LARGEST);

    igraph_is_simple(&g, &is_simple);
    // For this method, it does not have to be connected
    igraph_is_connected(&g, &is_connected, IGRAPH_STRONG);
    igraph_is_bipartite(&g, &is_bipartite, NULL);

    printf("Edges: ");
    igraph_vector_int_init(&edges, 0);
    igraph_get_edgelist(&g, &edges, false);
    igraph_vector_int_print(&edges);

    printf("Degrees: ");
    igraph_vector_int_init(&degrees, 0);
    igraph_degree(&g, &degrees, igraph_vss_all(), IGRAPH_ALL, false);
    igraph_vector_int_print(&degrees);

    printf("vcount: %d\n", (int)igraph_vcount(&g));
    printf("Simple: %d\n", is_simple);
    printf("Connected: %d\n", is_connected);
    printf("Bipartite: %d\n", is_bipartite);

    igraph_vector_int_destroy(&degrees);
    igraph_vector_int_destroy(&edges);
    igraph_vector_int_destroy(&ds1);
    igraph_vector_int_destroy(&ds2);
    igraph_destroy(&g);

    // undirected bipartite multigraph.
//    igraph_bool_t has_multi;
    printf("===Multigraph===\n");

    igraph_integer_t deg1l_m[] = {2, 3, 1};
    igraph_integer_t deg2l_m[]= {4, 2};

    igraph_vector_int_init_array(&ds1, deg1l_m, 3);
    igraph_vector_int_init_array(&ds2, deg2l_m, 2);

    igraph_realize_bipartite_degree_sequence(&g, &ds1, &ds2, IGRAPH_MULTI_SW, IGRAPH_REALIZE_DEGSEQ_LARGEST);

    printf("Edges: ");
    igraph_vector_int_init(&edges, 0);
    igraph_get_edgelist(&g, &edges, false);
    igraph_vector_int_print(&edges);

    printf("Degrees: ");
    igraph_vector_int_init(&degrees, 0);
    igraph_degree(&g, &degrees, igraph_vss_all(), IGRAPH_ALL, false);
    igraph_vector_int_print(&degrees);

    igraph_has_multiple(&g, &has_multi);
    // For this method, it does not need to be connected
    igraph_is_connected(&g, &is_connected, IGRAPH_STRONG);
    igraph_is_bipartite(&g, &is_bipartite, NULL);

    printf("vcount: %d\n", (int)igraph_vcount(&g));
    printf("Has multiedges: %d\n", has_multi);
    printf("Connected: %d\n", is_connected);
    printf("Bipartite: %d\n", is_bipartite);

    igraph_vector_int_destroy(&degrees);
    igraph_vector_int_destroy(&edges);
    igraph_vector_int_destroy(&ds1);
    igraph_vector_int_destroy(&ds2);
    igraph_destroy(&g);


//    printf("===Random===\n");
//    igraph_t rg;
//    //igraph_vector_int_t dss1, dss2;
//    igraph_integer_t n1 = 5;
//    igraph_integer_t n2 = 5;
//    igraph_integer_t success = 0;
//    igraph_integer_t fail = 0;
//    igraph_integer_t total = 100;
//    igraph_bool_t potentially_connected;
//    igraph_vector_int_t dss1;
//    igraph_vector_int_t dss2;

        // Testing simple
//    printf("MULTI\n");
//    for (igraph_integer_t i=0;i < total;i++) {
//        // This ensures that we don't get degree sequences which have a 0.
//        while (true) {
//            igraph_bipartite_game_gnp(&rg, NULL, n1, n2, 0.5, false, IGRAPH_ALL);
//            igraph_rewire(&rg, 5, IGRAPH_REWIRING_SIMPLE);
//
//            igraph_vector_int_init(&ds1, 0);
//            igraph_degree(&rg, &ds1, igraph_vss_range(0, n1), IGRAPH_ALL, false);
//
//            igraph_vector_int_init(&ds2, 0);
//            igraph_degree(&rg, &ds2, igraph_vss_range(n1, n1 + n2), IGRAPH_ALL, false);
//
//            if ((!igraph_vector_int_contains(&ds1, 0)) && (!igraph_vector_int_contains(&ds2, 0))) {
//                break;
//            } else {
//                igraph_destroy(&rg);
//                igraph_vector_int_destroy(&ds1);
//                igraph_vector_int_destroy(&ds2);
//            }
//        }
//
//        // Check if potentially connected
//        if (igraph_ecount(&rg) >= igraph_vcount(&rg) - 1) {
//            potentially_connected = true;
//        } else {
//            potentially_connected = false;
//        }
//
//        igraph_realize_bipartite_degree_sequence(&g, &ds1, &ds2, IGRAPH_MULTI_SW, IGRAPH_REALIZE_DEGSEQ_LARGEST);
//
//        igraph_has_multiple(&g, &has_multi);
//        igraph_is_connected(&g, &is_connected, IGRAPH_STRONG);
//        igraph_is_bipartite(&g, &is_bipartite, NULL);
//
////        printf("vcount: %d\n", (int)igraph_vcount(&g));
////        printf("Has multiedges: %d\n", has_multi);
////        printf("Is simple: %d\n", is_simple);
////        printf("Connected: %d\n", is_connected);
////        printf("Bipartite: %d\n", is_bipartite);
//
//        // Generated graph degree sequences
//        //printf("DS1: ");
//        igraph_vector_int_init(&dss1, 0);
//        igraph_degree(&g, &dss1, igraph_vss_range(0, n1), IGRAPH_ALL, false);
//        //igraph_vector_int_print(&dss1);
//        //printf("DS2: ");
//        igraph_vector_int_init(&dss2, 0);
//        igraph_degree(&g, &dss2, igraph_vss_range(n1, n1 + n2), IGRAPH_ALL, false);
//        //igraph_vector_int_print(&dss2);
//
//        // Success criteria:
//        // If potentially connected and the graph is connected
//        // If not potentially connected, and the graph is bipartite and simple
//
//        // Prelims: Are degrees the same? Node count?
//        // Degree sequences for the nodes?
//        if (!is_bipartite) {
//            printf("not bipartite\n");
//            fail++;
//        } else if (igraph_ecount(&rg) != igraph_ecount(&g)) {
//            printf("Edge counts don't match.\n");
//            fail++;
//        } else if (igraph_vcount(&rg) != igraph_vcount(&g)) {
//            printf("Vertex counts don't match.\n");
//            fail++;
//        } else if ((!igraph_vector_int_all_e(&ds1, &dss1)) || (!igraph_vector_int_all_e(&ds2, &dss2))) {
//            printf("Degree sequences doesn't match.\n");
//            fail++;
//        } else if (is_bipartite) {
//            success++;
//        }
//        else {
//            printf("Else reached with dss:\n");
//            igraph_vector_int_print(&ds1);
//            igraph_vector_int_print(&ds2);
//            printf("vcount: %d\n", (int)igraph_vcount(&g));
//            printf("Has multiedges: %d\n", has_multi);
//            printf("Is simple: %d\n", is_simple);
//            printf("Potentially connected: %d\n", potentially_connected);
//            printf("Connected: %d\n", is_connected);
//            printf("Bipartite: %d\n", is_bipartite);
//            fail++;
//        }
//
//        igraph_vector_int_destroy(&ds1);
//        igraph_vector_int_destroy(&ds2);
//        igraph_vector_int_destroy(&dss1);
//        igraph_vector_int_destroy(&dss2);
//        igraph_vector_int_destroy(&degrees);
//        igraph_vector_int_destroy(&edges);
//        igraph_destroy(&rg);
//        igraph_destroy(&g);
//    }
//    // igraph_bipartate_game_gnp, use p
//    // rewire with low prob, should remove the bipartite
//    // then test that the function fails - it should if the dss are not bipartite
//    //igraph_vector_bool_destroy(&types);
//
//    printf("Total: %d\n", (int)total);
//    printf("Success: %d\n", (int)success);
//    printf("Fail: %d\n", (int)fail);

//    // Testing simple
//    for (igraph_integer_t i=0;i < total;i++) {
//        // This ensures that we don't get degree sequences which have a 0.
//        while (true) {
//            igraph_bipartite_game_gnp(&rg, NULL, n1, n2, 0.5, false, IGRAPH_ALL);
//            igraph_rewire(&rg, 5, IGRAPH_REWIRING_SIMPLE);
//
//            igraph_vector_int_init(&ds1, 0);
//            igraph_degree(&rg, &ds1, igraph_vss_range(0, n1), IGRAPH_ALL, false);
//
//            igraph_vector_int_init(&ds2, 0);
//            igraph_degree(&rg, &ds2, igraph_vss_range(n1, n1 + n2), IGRAPH_ALL, false);
//
//            if ((!igraph_vector_int_contains(&ds1, 0)) && (!igraph_vector_int_contains(&ds2, 0))) {
//                break;
//            } else {
//                igraph_destroy(&rg);
//                igraph_vector_int_destroy(&ds1);
//                igraph_vector_int_destroy(&ds2);
//            }
//        }
//
//        // Check if potentially connected
//        if (igraph_ecount(&rg) >= igraph_vcount(&rg) - 1) {
//            potentially_connected = true;
//        } else {
//            potentially_connected = false;
//        }
//
//        igraph_realize_bipartite_degree_sequence(&g, &ds1, &ds2, IGRAPH_SIMPLE_SW, IGRAPH_REALIZE_DEGSEQ_SMALLEST);
//
//        igraph_has_multiple(&g, &has_multi);
//        igraph_is_connected(&g, &is_connected, IGRAPH_STRONG);
//        igraph_is_bipartite(&g, &is_bipartite, NULL);
//
//        printf("vcount: %d\n", (int)igraph_vcount(&g));
//        printf("Has multiedges: %d\n", has_multi);
//        printf("Is simple: %d\n", is_simple);
//        printf("Connected: %d\n", is_connected);
//        printf("Bipartite: %d\n", is_bipartite);
//
//        // Generated graph degree sequences
//        //printf("DS1: ");
//        igraph_vector_int_init(&dss1, 0);
//        igraph_degree(&g, &dss1, igraph_vss_range(0, n1), IGRAPH_ALL, false);
//        //igraph_vector_int_print(&dss1);
//        //printf("DS2: ");
//        igraph_vector_int_init(&dss2, 0);
//        igraph_degree(&g, &dss2, igraph_vss_range(n1, n1 + n2), IGRAPH_ALL, false);
//        //igraph_vector_int_print(&dss2);
//
//        // Success criteria:
//        // If potentially connected and the graph is connected
//        // If not potentially connected, and the graph is bipartite and simple
//
//        // Prelims: Are degrees the same? Node count?
//        // Degree sequences for the nodes?
//        if (!is_bipartite) {
//            printf("not bipartite\n");
//            fail++;
//        } else if (igraph_ecount(&rg) != igraph_ecount(&g)) {
//            printf("Edge counts don't match.\n");
//            fail++;
//        } else if (igraph_vcount(&rg) != igraph_vcount(&g)) {
//            printf("Vertex counts don't match.\n");
//            fail++;
//        } else if ((!igraph_vector_int_all_e(&ds1, &dss1)) || (!igraph_vector_int_all_e(&ds2, &dss2))) {
//            printf("Degree sequences doesn't match.\n");
//            fail++;
//        } else if (potentially_connected & !is_connected) {
//            printf("Potentially connected didn't end up connected. Dss:\n");
//            igraph_vector_int_print(&ds1);
//            igraph_vector_int_print(&ds2);
//            fail++;
//        } else if (is_bipartite & potentially_connected & is_connected & is_simple) {
//            success++;
//        } else if (is_bipartite & !potentially_connected & !is_connected & is_simple) {
//            success++;
//        }
//        else {
//            printf("Else reached with dss:\n");
//            igraph_vector_int_print(&ds1);
//            igraph_vector_int_print(&ds2);
//            printf("vcount: %d\n", (int)igraph_vcount(&g));
//            printf("Has multiedges: %d\n", has_multi);
//            printf("Is simple: %d\n", is_simple);
//            printf("Potentially connected: %d\n", potentially_connected);
//            printf("Connected: %d\n", is_connected);
//            printf("Bipartite: %d\n", is_bipartite);
//            fail++;
//        }
//
//        igraph_vector_int_destroy(&ds1);
//        igraph_vector_int_destroy(&ds2);
//        igraph_vector_int_destroy(&dss1);
//        igraph_vector_int_destroy(&dss2);
//        igraph_vector_int_destroy(&degrees);
//        igraph_vector_int_destroy(&edges);
//        igraph_destroy(&rg);
//        igraph_destroy(&g);
//    }
//    // igraph_bipartate_game_gnp, use p
//    // rewire with low prob, should remove the bipartite
//    // then test that the function fails - it should if the dss are not bipartite
//    //igraph_vector_bool_destroy(&types);
//
//    printf("Total: %d\n", (int)total);
//    printf("Success: %d\n", (int)success);
//    printf("Fail: %d\n", (int)fail);

    VERIFY_FINALLY_STACK();
    return 0;
}
