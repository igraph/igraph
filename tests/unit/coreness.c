/*
   igraph library.
   Copyright (C) 2022  The igraph development team <igraph@igraph.org>

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

#include "test_utilities.h"

void validate_coreness(
    const igraph_t* graph, const igraph_vector_int_t* coreness,
    igraph_neimode_t mode
) {
    igraph_int_t i, j, min_coreness, max_coreness;
    igraph_int_t nv = igraph_vcount(graph);
    igraph_t subgraph;
    igraph_vs_t vs;
    igraph_vector_int_t vids, degree;

    IGRAPH_ASSERT(igraph_vector_int_size(coreness) == nv);
    if (igraph_vcount(graph) < 1) {
        return;
    }

    min_coreness = igraph_vector_int_min(coreness);
    max_coreness = igraph_vector_int_max(coreness);

    IGRAPH_ASSERT(min_coreness >= 0);

    igraph_vector_int_init(&vids, 0);
    igraph_vector_int_init(&degree, 0);

    for (i = max_coreness; i >= 0; i--) {
        igraph_vector_int_clear(&vids);
        for (j = 0; j < nv; j++) {
            if (VECTOR(*coreness)[j] >= i) {
                igraph_vector_int_push_back(&vids, j);
            }
        }

        igraph_vs_vector(&vs, &vids);
        igraph_induced_subgraph(graph, &subgraph, vs, IGRAPH_SUBGRAPH_AUTO);
        igraph_vs_destroy(&vs);

        igraph_degree(&subgraph, &degree, igraph_vss_all(), mode, IGRAPH_LOOPS_TWICE);
        for (j = 0; j < igraph_vcount(&subgraph); j++) {
            IGRAPH_ASSERT(VECTOR(degree)[j] >= i);
        }

        igraph_destroy(&subgraph);
    }

    igraph_vector_int_destroy(&degree);
    igraph_vector_int_destroy(&vids);
}

void test_graph(const igraph_t* graph, igraph_bool_t print) {
    igraph_vector_int_t coreness;

    igraph_vector_int_init(&coreness, 0);

    if (igraph_is_directed(graph)) {
        igraph_coreness(graph, &coreness, IGRAPH_ALL);
        validate_coreness(graph, &coreness, IGRAPH_ALL);
        if (print) {
            printf("mode = ALL: ");
            print_vector_int(&coreness);
        }

        igraph_coreness(graph, &coreness, IGRAPH_OUT);
        validate_coreness(graph, &coreness, IGRAPH_OUT);
        if (print) {
            printf("mode = OUT: ");
            print_vector_int(&coreness);
        }

        igraph_coreness(graph, &coreness, IGRAPH_IN);
        validate_coreness(graph, &coreness, IGRAPH_IN);
        if (print) {
            printf("mode = IN: ");
            print_vector_int(&coreness);
        }
    } else {
        igraph_coreness(graph, &coreness, IGRAPH_ALL);
        validate_coreness(graph, &coreness, IGRAPH_ALL);
        if (print) {
            print_vector_int(&coreness);
        }
    }

    igraph_vector_int_destroy(&coreness);
}

void add_loop_and_multiple_edges(igraph_t* graph, igraph_real_t loop_prob, igraph_real_t multi_prob) {
    igraph_int_t i, n, from, to;
    igraph_vector_int_t extra_edges;

    igraph_vector_int_init(&extra_edges, 0);

    n = igraph_vcount(graph);
    for (i = 0; i < n; i++) {
        if (RNG_UNIF01() < loop_prob) {
            igraph_vector_int_push_back(&extra_edges, i);
            igraph_vector_int_push_back(&extra_edges, i);
        }
    }

    n = igraph_ecount(graph);
    for (i = 0; i < n; i++) {
        if (RNG_UNIF01() < multi_prob) {
            igraph_edge(graph, i, &from, &to);
            igraph_vector_int_push_back(&extra_edges, from);
            igraph_vector_int_push_back(&extra_edges, to);
        }
    }

    igraph_add_edges(graph, &extra_edges, 0);

    igraph_vector_int_destroy(&extra_edges);
}

void remove_some_edges(igraph_t* graph, igraph_real_t prob) {
    igraph_int_t i, n;
    igraph_vector_int_t to_remove;

    igraph_vector_int_init(&to_remove, 0);

    n = igraph_ecount(graph);
    for (i = 0; i < n; i++) {
        if (igraph_rng_get_unif01(igraph_rng_default()) < prob) {
            igraph_vector_int_push_back(&to_remove, i);
        }
    }

    igraph_delete_edges(graph, igraph_ess_vector(&to_remove));

    igraph_vector_int_destroy(&to_remove);
}

int main(void) {
    igraph_t g;

    igraph_rng_seed(igraph_rng_default(), 137);

    /* Empty and singleton graph */
    printf("Empty graph.\n");
    igraph_empty(&g, 0, IGRAPH_UNDIRECTED);
    test_graph(&g, /* print = */ 1);
    igraph_destroy(&g);
    printf("Singleton graph.\n");
    igraph_empty(&g, 1, IGRAPH_UNDIRECTED);
    test_graph(&g, /* print = */ 1);
    igraph_destroy(&g);

    /* Simple full graph */
    printf("Full graph.\n");
    igraph_full(&g, 5, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    test_graph(&g, /* print = */ 1);
    igraph_destroy(&g);

    /* Full graph with loops */
    printf("Full graph with loops.\n");
    igraph_full(&g, 5, IGRAPH_UNDIRECTED, IGRAPH_LOOPS);
    test_graph(&g, /* print = */ 1);
    igraph_destroy(&g);

    /* Full directed graph */
    printf("Full directed graph.\n");
    igraph_full(&g, 5, IGRAPH_DIRECTED, IGRAPH_NO_LOOPS);
    test_graph(&g, /* print = */ 1);
    igraph_destroy(&g);

    /* Full directed graph */
    printf("Full directed graph with loops.\n");
    igraph_full(&g, 5, IGRAPH_DIRECTED, IGRAPH_LOOPS);
    test_graph(&g, /* print = */ 1);
    igraph_destroy(&g);

    /* Zachary karate club */
    printf("Zachary karate club.\n");
    igraph_famous(&g, "zachary");
    test_graph(&g, /* print = */ 1);
    igraph_destroy(&g);

    /* Zachary karate club, randomly directed edges */
    printf("Zachary karate club, directed edges.\n");
    igraph_famous(&g, "zachary");
    igraph_to_directed(&g, IGRAPH_TO_DIRECTED_MUTUAL);
    remove_some_edges(&g, 0.2);
    test_graph(&g, /* print = */ 1);
    igraph_destroy(&g);

    /* Zachary karate club with random loops and multi-edges */
    printf("Zachary karate club with loops and multi-edges.\n");
    for (int i = 0; i < 20; i++) {
        igraph_famous(&g, "zachary");
        add_loop_and_multiple_edges(&g, 0.5, 0.2);
        test_graph(&g, /* print = */ 0);
        igraph_destroy(&g);
    }

    /* Geometric random graph */
    printf("Geometric random graph.\n");
    igraph_grg_game(&g, 100, 0.2, /* torus = */ 0, /* x = */ 0, /* y = */ 0);
    test_graph(&g, /* print = */ 0);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();

    return 0;
}
