/*
   IGraph library.
   Copyright (C) 2013-2021  The igraph development team <igraph@igraph.org>

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

#include "bench.h"

igraph_integer_t test_adj(igraph_t *g, igraph_adjlist_t *adj)
{
    igraph_integer_t dummy = 0;

    for (int i = 0; i < igraph_vcount(g); i++) {
        igraph_vector_int_t *neis = igraph_adjlist_get(adj, i);
        igraph_integer_t nneis = igraph_vector_int_size(neis);
        for (int j = 0; j < nneis; j++) {
            igraph_integer_t neighbor = VECTOR(*neis)[j];
            dummy += neighbor;
        }
    }
    return dummy;
}

igraph_integer_t test_inc_adj(igraph_t *g, igraph_inclist_t *inc, igraph_adjlist_t *adj)
{
    igraph_integer_t dummy = 0;

    for (int i = 0; i < igraph_vcount(g); i++) {
        igraph_vector_int_t *adjs = igraph_adjlist_get(adj, i);
        igraph_vector_int_t *incs = igraph_inclist_get(inc, i);
        igraph_integer_t nneis = igraph_vector_int_size(adjs);
        for (int j = 0; j < nneis; j++) {
            igraph_integer_t edge = VECTOR(*incs)[j];
            igraph_integer_t neighbor = VECTOR(*adjs)[j];
            dummy += neighbor + edge;
        }
    }
    return dummy;
}

igraph_integer_t test_inc_other(igraph_t *g, igraph_inclist_t *inc)
{
    igraph_integer_t dummy = 0;
    for (int i = 0; i < igraph_vcount(g); i++) {
        igraph_vector_int_t *neis = igraph_inclist_get(inc, i);
        igraph_integer_t nneis = igraph_vector_int_size(neis);
        for (int j = 0; j < nneis; j++) {
            igraph_integer_t edge = VECTOR(*neis)[j];
            igraph_integer_t neighbor = IGRAPH_OTHER(g, edge, i);
            dummy += neighbor;
        }
    }
    return dummy;
}

igraph_integer_t test_inc_to(igraph_t *g, igraph_inclist_t *inc)
{
    igraph_integer_t dummy = 0;
    for (int i = 0; i < igraph_vcount(g); i++) {
        igraph_vector_int_t *neis = igraph_inclist_get(inc, i);
        igraph_integer_t nneis = igraph_vector_int_size(neis);
        for (int j = 0; j < nneis; j++) {
            igraph_integer_t edge = VECTOR(*neis)[j];
            igraph_integer_t neighbor = IGRAPH_TO(g, edge);
            dummy += neighbor;
        }
    }
    return dummy;
}

igraph_integer_t test_inc_nop(igraph_t *g, igraph_inclist_t *inc)
{
    igraph_integer_t dummy = 0;
    for (int i = 0; i < igraph_vcount(g); i++) {
        igraph_vector_int_t *neis = igraph_inclist_get(inc, i);
        igraph_integer_t nneis = igraph_vector_int_size(neis);
        for (int j = 0; j < nneis; j++) {
            igraph_integer_t edge = VECTOR(*neis)[j];
            dummy += edge;
        }
    }
    return dummy;
}

void do_benchmarks(char *name, igraph_t *g, int repeat) {
    igraph_adjlist_t adj;
    igraph_inclist_t inc;
    volatile igraph_integer_t result;

    printf("%s", name);
    BENCH("1 initialize adjlist.",
            REPEAT(igraph_adjlist_init(g, &adj, IGRAPH_ALL, IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE), repeat);
         );

    printf("%s", name);
    BENCH("2 initialize inclist.",
            REPEAT(igraph_inclist_init(g, &inc, IGRAPH_ALL, IGRAPH_NO_LOOPS), repeat);
         );

    printf("%s", name);
    BENCH("3 go over vertices (multiple times) using adjlist.",
            REPEAT(result = test_adj(g, &adj), repeat);
         );

    printf("%s", name);
    BENCH("4 go over vertices (multiple times) using inclist, IGRAPH_OTHER.",
            REPEAT(result = test_inc_other(g, &inc), repeat);
         );

    printf("%s", name);
    BENCH("5 go over vertices (multiple times) using inclist, IGRAPH_TO.",
            REPEAT(result = test_inc_to(g, &inc), repeat);
         );

    printf("%s", name);
    BENCH("6 go over edges using inclist, don't retrieve vertex.",
            REPEAT(result = test_inc_nop(g, &inc), repeat);
         );

    printf("%s", name);
    BENCH("7 go over edges and vertices using adjlist and inclist.",
            REPEAT(result = test_inc_adj(g, &inc, &adj), repeat);
         );

    igraph_adjlist_destroy(&adj);
    igraph_inclist_destroy(&inc);

    BENCH("8 initialize adjlist, include loops and multiple (which aren't present).",
            REPEAT(igraph_adjlist_init(g, &adj, IGRAPH_ALL, IGRAPH_LOOPS, IGRAPH_MULTIPLE), repeat);
         );

    BENCH("9 initialize inclist, include loops (which aren't present).",
            REPEAT(igraph_inclist_init(g, &inc, IGRAPH_ALL, IGRAPH_LOOPS), repeat);
         );

    igraph_adjlist_destroy(&adj);
    igraph_inclist_destroy(&inc);
}

int main() {
    igraph_t g;

    igraph_rng_seed(igraph_rng_default(), 42);
    BENCH_INIT();

    printf("Full graph tests:\n");
    igraph_full(&g, 10000, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    do_benchmarks(" fg - ", &g, 1);
    igraph_destroy(&g);

    printf("\nRandom graph tests:\n");
    igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNM, 10000, 49994999, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    do_benchmarks(" rg - ", &g, 1);
    igraph_destroy(&g);

    printf("\nSmall graph tests:\n");
    igraph_full(&g, 1000, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    do_benchmarks(" sg - ", &g, 100);
    igraph_destroy(&g);


    return 0;
}
