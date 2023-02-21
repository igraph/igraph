/*
   IGraph library.
   Copyright (C) 2021  The igraph development team <igraph@igraph.org>

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

/* Direct test for undirected forests:
 * Decompose into connected components and check that each one is a tree */
igraph_bool_t is_forest(const igraph_t *graph) {
    igraph_bool_t res = 1;
    igraph_graph_list_t components;

    igraph_graph_list_init(&components, 0);

    igraph_decompose(graph, &components, IGRAPH_WEAK, -1, 0);

    igraph_integer_t n = igraph_graph_list_size(&components);
    for (igraph_integer_t i=0; i < n; ++i) {
        igraph_is_tree(igraph_graph_list_get_ptr(&components, i), &res, NULL, IGRAPH_ALL);

        if (! res) {
            break;
        }
    }

    igraph_graph_list_destroy(&components);

    return res;
}

/* This test generates 'trials' random graphs, allowing self-loops and
 * multi-edges, and exercises the forest detection function on each. */
int main(void) {
    const igraph_integer_t n = 100; /* vertex count */
    const igraph_integer_t m = n / 2; /* edge count */
    const igraph_integer_t trials = 300;
    igraph_integer_t true_count;
    igraph_bool_t res1, res2;
    igraph_vector_int_t edges;

    igraph_rng_seed(igraph_rng_default(), 847532);

    igraph_vector_int_init(&edges, 2 * m);

    RNG_BEGIN();

    true_count = 0;
    for (igraph_integer_t k = 0; k < trials; ++k) {
        igraph_t graph;

        for (igraph_integer_t i = 0; i < 2 * m; ++i) {
            VECTOR(edges)[i] = RNG_INTEGER(0, n - 1);
        }

        igraph_create(&graph, &edges, n, IGRAPH_UNDIRECTED);

        igraph_is_forest(&graph, &res1, NULL, IGRAPH_ALL);
        res2 = is_forest(&graph);

        if (res1 != res2) {
            printf("Invalid result for the following graph.\nExpected result: %s\nActual result: %s\n\n",
                   res2 ? "true" : "false",
                   res1 ? "true" : "false");
            print_graph(&graph);
        }

        igraph_destroy(&graph);

        if (res1 != res2) {
            break;
        }

        if (res1) {
            true_count += 1;
        }
    }

    RNG_END();

    igraph_vector_int_destroy(&edges);

    /* ensure that test fails if we got an unexpected result */
    IGRAPH_ASSERT(res1 == res2);

    /* Check that the test is covering all cases: what fraction of random graphs was a tree?
     * With random multigraphs, this should be around 20-30% when the edge count is half
     * that of the vertex count/ */
    printf("Fraction of random multigraphs which were a forest: %g\n", ((double) true_count) / trials);

    VERIFY_FINALLY_STACK();

    return 0;
}
