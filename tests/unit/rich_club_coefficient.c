/*
   igraph library.
   Copyright (C) 2026  The igraph development team <igraph@igraph.org>

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

#include <stdio.h>
#include <igraph.h>
#include "test_utilities.h"

void test_null_graph(void) {
    igraph_t graph;
    igraph_vector_t result;

    igraph_vector_init(&result, 0);
    igraph_empty(&graph, 0, IGRAPH_UNDIRECTED);

    printf("Test 1: null graph\n");
    igraph_rich_club_coefficient(&graph, NULL, &result, IGRAPH_ALL, false, false);
    print_vector(&result);

    igraph_vector_destroy(&result);
    igraph_destroy(&graph);
}

void test_singleton_graph(void) {
    igraph_t graph;
    igraph_vector_t result;

    igraph_vector_init(&result, 0);

    /* Singleton no loop */
    printf("Test 2a: singleton without loop\n");
    igraph_empty(&graph, 1, IGRAPH_UNDIRECTED);
    igraph_rich_club_coefficient(&graph, NULL, &result, IGRAPH_ALL, false, false);
    print_vector(&result);
    igraph_destroy(&graph);

    /* Singleton with self loop */
    printf("Test 2b: singleton with loop (loops=true)\n");
    igraph_empty(&graph, 1, IGRAPH_UNDIRECTED);
    igraph_add_edge(&graph, 0, 0);
    igraph_rich_club_coefficient(&graph, NULL, &result, IGRAPH_ALL, true, false);
    print_vector(&result);
    igraph_destroy(&graph);

    igraph_vector_destroy(&result);
}

void test_complete_graph(void) {
    igraph_t graph;
    igraph_vector_t result;

    igraph_vector_init(&result, 0);
    igraph_full(&graph, 5, IGRAPH_UNDIRECTED, false);

    printf("Test 3: complete graph K5\n");
    igraph_rich_club_coefficient(&graph, NULL, &result, IGRAPH_ALL, false, false);
    print_vector(&result);

    igraph_vector_destroy(&result);
    igraph_destroy(&graph);
}

void test_cycle_graph(void) {
    igraph_t graph;
    igraph_vector_t result;

    igraph_vector_init(&result, 0);
    igraph_ring(&graph, 5, IGRAPH_UNDIRECTED, false, true);

    printf("Test 4: cycle graph C5\n");
    igraph_rich_club_coefficient(&graph, NULL, &result, IGRAPH_ALL, false, false);
    print_vector(&result);

    igraph_vector_destroy(&result);
    igraph_destroy(&graph);
}

void test_star_graph(void) {
    igraph_t graph;
    igraph_vector_t result;

    igraph_vector_init(&result, 0);
    igraph_star(&graph, 5, IGRAPH_STAR_UNDIRECTED, 0);

    printf("Test 5: star graph S5\n");
    igraph_rich_club_coefficient(&graph, NULL, &result, IGRAPH_ALL, false, false);
    print_vector(&result);

    igraph_vector_destroy(&result);
    igraph_destroy(&graph);
}

void test_directed_reciprocal(void) {
    igraph_t graph;
    igraph_vector_t result;

    igraph_vector_init(&result, 0);
    igraph_empty(&graph, 2, IGRAPH_DIRECTED);
    igraph_add_edge(&graph, 0, 1);
    igraph_add_edge(&graph, 1, 0);

    printf("Test 6a: directed reciprocal (directed=true)\n");
    igraph_rich_club_coefficient(&graph, NULL, &result, IGRAPH_OUT, false, true);
    print_vector(&result);

    printf("Test 6b: directed reciprocal (directed=false)\n");
    igraph_rich_club_coefficient(&graph, NULL, &result, IGRAPH_ALL, false, false);
    print_vector(&result);

    igraph_vector_destroy(&result);
    igraph_destroy(&graph);
}

void test_weighted(void) {
    igraph_t graph;
    igraph_vector_t result;
    igraph_vector_t weights;

    igraph_vector_init(&result, 0);
    igraph_full(&graph, 3, IGRAPH_UNDIRECTED, false);
    igraph_vector_init_real(&weights, 3, 2.0, 3.0, 5.0);

    printf("Test 7: weighted triangle\n");
    igraph_rich_club_coefficient(&graph, &weights, &result, IGRAPH_ALL, false, false);
    print_vector(&result);

    igraph_vector_destroy(&weights);
    igraph_vector_destroy(&result);
    igraph_destroy(&graph);
}

void test_error_checks(void) {
    igraph_t graph;
    igraph_vector_t result;
    igraph_vector_t weights;

    igraph_vector_init(&result, 0);
    igraph_empty(&graph, 4, IGRAPH_DIRECTED);
    igraph_add_edge(&graph, 0, 1);

    printf("Test 8: error handling\n");
    igraph_vector_init(&weights, 5); // wrong weights size (!= 1)
    CHECK_ERROR(igraph_rich_club_coefficient(&graph, &weights, &result, IGRAPH_ALL, false, false),
                IGRAPH_EINVAL);

    CHECK_ERROR(igraph_rich_club_coefficient(&graph, NULL, &result, (igraph_neimode_t) 999, false, false),
                IGRAPH_EINVMODE);

    igraph_vector_destroy(&weights);
    igraph_vector_destroy(&result);
    igraph_destroy(&graph);
}

int main(void) {
    test_null_graph();
    printf("\n");

    test_singleton_graph();
    printf("\n");

    test_complete_graph();
    printf("\n");

    test_cycle_graph();
    printf("\n");

    test_star_graph();
    printf("\n");

    test_directed_reciprocal();
    printf("\n");

    test_weighted();
    printf("\n");

    test_error_checks();

    VERIFY_FINALLY_STACK();
    return 0;
}
