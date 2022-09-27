/*
   IGraph library.
   Copyright (C) 2013-2022  The igraph development team <igraph@igraph.org>

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

typedef struct igraph_incadjlist_inter_t {
    igraph_integer_t length;
    igraph_vector_int_t *incadjs;
} igraph_incadjlist_inter_t;

#define igraph_incadjlist_inter_get(il,no) (&(il)->incadjs[(igraph_integer_t)(no)])
#define igraph_incadjlist_sep_get_inc(il,no) (&(il)->incs[(igraph_integer_t)(no)])
#define igraph_incadjlist_sep_get_adj(il,no) (&(il)->adjs[(igraph_integer_t)(no)])

void igraph_incadjlist_inter_destroy(igraph_incadjlist_inter_t *il) {
    igraph_integer_t i;
    for (i = 0; i < il->length; i++) {
        /* This works if some igraph_vector_int_t's contain NULL,
           because igraph_vector_int_destroy can handle this. */
        igraph_vector_int_destroy(&il->incadjs[i]);
    }
    IGRAPH_FREE(il->incadjs);
}

igraph_error_t igraph_incadjlist_inter_init(const igraph_t *graph,
                        igraph_incadjlist_inter_t *il,
                        igraph_neimode_t mode) {
    igraph_integer_t i, j, n;
    igraph_vector_int_t tmp;

    if (mode != IGRAPH_IN && mode != IGRAPH_OUT && mode != IGRAPH_ALL) {
        IGRAPH_ERROR("Cannot create incadjlist.", IGRAPH_EINVMODE);
    }

    igraph_vector_int_init(&tmp, 0);
    IGRAPH_FINALLY(igraph_vector_int_destroy, &tmp);

    if (!igraph_is_directed(graph)) {
        mode = IGRAPH_ALL;
    }

    il->length = igraph_vcount(graph);
    il->incadjs = IGRAPH_CALLOC(il->length, igraph_vector_int_t);
    if (il->incadjs == 0) {
        IGRAPH_ERROR("Cannot create incadjlist.", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
    }
    IGRAPH_FINALLY(igraph_incadjlist_inter_destroy, il);

    for (i = 0; i < il->length; i++) {
        //IGRAPH_ALLOW_INTERRUPTION();

        IGRAPH_CHECK(igraph_incident(graph, &tmp, i, mode));

        n = igraph_vector_int_size(&tmp);
        IGRAPH_CHECK(igraph_vector_int_init(&il->incadjs[i], n * 2));

        for (j = 0; j < n; j++) {
            VECTOR(il->incadjs[i])[j * 2] = VECTOR(tmp)[j];
            VECTOR(il->incadjs[i])[j * 2 + 1] =
                IGRAPH_OTHER(graph, VECTOR(tmp)[j], i);
        }
    }

    igraph_vector_int_destroy(&tmp);
    IGRAPH_FINALLY_CLEAN(2); /* + igraph_incadjlist_inter_destroy */

    return IGRAPH_SUCCESS;
}


typedef struct igraph_incadjlist_sep_t {
    igraph_integer_t length;
    igraph_vector_int_t *incs;
    igraph_vector_int_t *adjs;
} igraph_incadjlist_sep_t;

void igraph_incadjlist_sep_destroy(igraph_incadjlist_sep_t *il) {
    igraph_integer_t i;
    for (i = 0; i < il->length; i++) {
        /* This works if some igraph_vector_int_t's contain NULL,
           because igraph_vector_int_destroy can handle this. */
        igraph_vector_int_destroy(&il->incs[i]);
        igraph_vector_int_destroy(&il->adjs[i]);
    }
    IGRAPH_FREE(il->incs);
    IGRAPH_FREE(il->adjs);
}

igraph_error_t igraph_incadjlist_sep_init(const igraph_t *graph,
                        igraph_incadjlist_sep_t *il,
                        igraph_neimode_t mode) {
    igraph_integer_t i, j, n;
    igraph_vector_int_t tmp;

    if (mode != IGRAPH_IN && mode != IGRAPH_OUT && mode != IGRAPH_ALL) {
        IGRAPH_ERROR("Cannot create incadjlist.", IGRAPH_EINVMODE);
    }

    igraph_vector_int_init(&tmp, 0);
    IGRAPH_FINALLY(igraph_vector_int_destroy, &tmp);

    if (!igraph_is_directed(graph)) {
        mode = IGRAPH_ALL;
    }

    il->length = igraph_vcount(graph);
    il->incs = IGRAPH_CALLOC(il->length, igraph_vector_int_t);
    if (il->incs == 0) {
        IGRAPH_ERROR("Cannot create incadjlist.", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
    }

    il->adjs = IGRAPH_CALLOC(il->length, igraph_vector_int_t);
    if (il->adjs == 0) {
        IGRAPH_FREE(il->incs);
        IGRAPH_ERROR("Cannot create incadjlist.", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
    }

    IGRAPH_FINALLY(igraph_incadjlist_sep_destroy, il);

    for (i = 0; i < il->length; i++) {
        //IGRAPH_ALLOW_INTERRUPTION();

        IGRAPH_CHECK(igraph_incident(graph, &tmp, i, mode));

        n = igraph_vector_int_size(&tmp);
        IGRAPH_CHECK(igraph_vector_int_init(&il->incs[i], n));
        IGRAPH_CHECK(igraph_vector_int_init(&il->adjs[i], n));

        for (j = 0; j < n; j++) {
            VECTOR(il->incs[i])[j] = VECTOR(tmp)[j];
            VECTOR(il->adjs[i])[j] =
                IGRAPH_OTHER(graph, VECTOR(tmp)[j], i);
        }

    }

    igraph_vector_int_destroy(&tmp);
    IGRAPH_FINALLY_CLEAN(2); /* + igraph_incadjlist_sep_destroy */

    return IGRAPH_SUCCESS;
}

/* In the below tests, the 'dummy' variable is used to prevent the compilers
 * from optimizing away the entire side-effects-free function. We use XOR
 * instead of e.g. addition in order to prevent triggering undefined behaviour. */

igraph_integer_t test_direct(igraph_t *g)
{
    igraph_integer_t dummy = 0;
    igraph_integer_t vcount = igraph_vcount(g);
    igraph_integer_t ecount = igraph_ecount(g);

    igraph_integer_t ei = 0;
    igraph_integer_t eo = 0;
    for (igraph_integer_t i = 0; i < vcount; i++) {
        while (ei < ecount && VECTOR(g->from)[VECTOR(g->oi)[ei]] == i) {
            igraph_integer_t neighbor = VECTOR(g->to)[VECTOR(g->oi)[ei]];
            dummy ^= neighbor ^ VECTOR(g->oi)[ei];
            ei++;
        }
        while (eo < ecount && VECTOR(g->to)[VECTOR(g->ii)[eo]] == i) {
            igraph_integer_t neighbor = VECTOR(g->from)[VECTOR(g->ii)[eo]];
            dummy ^= neighbor ^ VECTOR(g->ii)[eo] ;
            eo++;
        }
    }
    return dummy;
}

igraph_integer_t test_adj(igraph_t *g, igraph_adjlist_t *adj)
{
    igraph_integer_t dummy = 0;
    igraph_integer_t vcount = igraph_vcount(g);

    for (igraph_integer_t i = 0; i < vcount; i++) {
        igraph_vector_int_t *neis = igraph_adjlist_get(adj, i);
        igraph_integer_t nneis = igraph_vector_int_size(neis);
        for (int j = 0; j < nneis; j++) {
            igraph_integer_t neighbor = VECTOR(*neis)[j];
            dummy ^= neighbor;
        }
    }
    return dummy;
}

igraph_integer_t test_inc_adj(igraph_t *g, igraph_inclist_t *inc, igraph_adjlist_t *adj)
{
    igraph_integer_t dummy = 0;
    igraph_integer_t vcount = igraph_vcount(g);

    for (igraph_integer_t i = 0; i < vcount; i++) {
        igraph_vector_int_t *adjs = igraph_adjlist_get(adj, i);
        igraph_vector_int_t *incs = igraph_inclist_get(inc, i);
        igraph_integer_t nneis = igraph_vector_int_size(adjs);
        for (int j = 0; j < nneis; j++) {
            igraph_integer_t edge = VECTOR(*incs)[j];
            igraph_integer_t neighbor = VECTOR(*adjs)[j];
            dummy ^= neighbor ^ edge;
        }
    }
    return dummy;
}

igraph_integer_t test_inc_other(igraph_t *g, igraph_inclist_t *inc)
{
    igraph_integer_t dummy = 0;
    igraph_integer_t vcount = igraph_vcount(g);

    for (igraph_integer_t i = 0; i < vcount; i++) {
        igraph_vector_int_t *neis = igraph_inclist_get(inc, i);
        igraph_integer_t nneis = igraph_vector_int_size(neis);
        for (igraph_integer_t j = 0; j < nneis; j++) {
            igraph_integer_t edge = VECTOR(*neis)[j];
            igraph_integer_t neighbor = IGRAPH_OTHER(g, edge, i);
            dummy ^= neighbor ^ edge;
        }
    }
    return dummy;
}

igraph_integer_t test_incadj_sep(igraph_t *g, igraph_incadjlist_sep_t *inc)
{
    igraph_integer_t dummy = 0;
    igraph_integer_t vcount = igraph_vcount(g);

    for (igraph_integer_t i = 0; i < vcount; i++) {
        igraph_vector_int_t *incs = igraph_incadjlist_sep_get_inc(inc, i);
        igraph_vector_int_t *adjs = igraph_incadjlist_sep_get_adj(inc, i);
        igraph_integer_t nneis = igraph_vector_int_size(incs);
        for (igraph_integer_t j = 0; j < nneis; j++) {
            igraph_integer_t edge = VECTOR(*incs)[j];
            igraph_integer_t neighbor = VECTOR(*adjs)[j];
            dummy ^= neighbor ^ edge;
        }
    }
    return dummy;
}

igraph_integer_t test_incadj_inter(igraph_t *g, igraph_incadjlist_inter_t *inc)
{
    igraph_integer_t dummy = 0;
    igraph_integer_t vcount = igraph_vcount(g);

    for (igraph_integer_t i = 0; i < vcount; i++) {
        igraph_vector_int_t *ias = igraph_incadjlist_inter_get(inc, i);
        igraph_integer_t nneis = igraph_vector_int_size(ias) / 2;
        for (igraph_integer_t j = 0; j < nneis; j++) {
            igraph_integer_t edge = VECTOR(*ias)[j * 2];
            igraph_integer_t neighbor = VECTOR(*ias)[j * 2 + 1];
            dummy ^= neighbor ^ edge;
        }
    }
    return dummy;
}

igraph_integer_t test_inc_to(igraph_t *g, igraph_inclist_t *inc)
{
    igraph_integer_t dummy = 0;
    igraph_integer_t vcount = igraph_vcount(g);

    for (igraph_integer_t i = 0; i < vcount; i++) {
        igraph_vector_int_t *neis = igraph_inclist_get(inc, i);
        igraph_integer_t nneis = igraph_vector_int_size(neis);
        for (igraph_integer_t j = 0; j < nneis; j++) {
            igraph_integer_t edge = VECTOR(*neis)[j];
            igraph_integer_t neighbor = IGRAPH_TO(g, edge);
            dummy ^= neighbor ^ edge;
        }
    }
    return dummy;
}

igraph_integer_t test_inc_nop(igraph_t *g, igraph_inclist_t *inc)
{
    igraph_integer_t dummy = 0;
    igraph_integer_t vcount = igraph_vcount(g);

    for (igraph_integer_t i = 0; i < vcount; i++) {
        igraph_vector_int_t *neis = igraph_inclist_get(inc, i);
        igraph_integer_t nneis = igraph_vector_int_size(neis);
        for (igraph_integer_t j = 0; j < nneis; j++) {
            igraph_integer_t edge = VECTOR(*neis)[j];
            dummy ^= edge;
        }
    }
    return dummy;
}

/* Used to prevent optimizing away the result.
 * This must be a global variable to prevent "variable set but not used"
 * warnings from some compilers. */
volatile igraph_integer_t result;

void do_benchmarks(char *name, igraph_t *g, int repeat) {
    igraph_adjlist_t adj;
    igraph_inclist_t inc;
    igraph_incadjlist_inter_t incadj_inter;
    igraph_incadjlist_sep_t incadj_sep;

    /* The init() call is in a REPEAT loop, so we include destroy() as well to avoid a memory leak.
     * This is representative of real use cases, where init/destroy should always come in pairs. */
    printf("%s", name);
    BENCH("1 init/destroy adjlist.",
            REPEAT(
              do {
                  igraph_adjlist_init(g, &adj, IGRAPH_ALL, IGRAPH_LOOPS_TWICE, IGRAPH_MULTIPLE);
                  igraph_adjlist_destroy(&adj);
              } while (0),
              repeat);
         );

    printf("%s", name);
    BENCH("2 init/destroy inclist.",
            REPEAT(
              do {
                  igraph_inclist_init(g, &inc, IGRAPH_ALL, IGRAPH_LOOPS_TWICE);
                  igraph_inclist_destroy(&inc);
              } while (0),
              repeat);
         );

    printf("%s", name);
    BENCH("3 init/destroy adjlist, remove loops and multiple (which aren't present).",
            REPEAT(
              do {
                  igraph_adjlist_init(g, &adj, IGRAPH_ALL, IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE);
                  igraph_adjlist_destroy(&adj);
              } while (0),
              repeat);
         );

    printf("%s", name);
    BENCH("4 init/destroy inclist, remove loops (which aren't present).",
            REPEAT(
              do {
                  igraph_inclist_init(g, &inc, IGRAPH_ALL, IGRAPH_NO_LOOPS);
                  igraph_inclist_destroy(&inc);
              } while (0),
              repeat);
         );

    /* Initialize adjlist  / inclist for the following benchmarks. */
    igraph_adjlist_init(g, &adj, IGRAPH_ALL, IGRAPH_LOOPS_TWICE, IGRAPH_MULTIPLE);
    igraph_inclist_init(g, &inc, IGRAPH_ALL, IGRAPH_LOOPS_TWICE);

    printf("%s", name);
    BENCH("5 go over vertices (multiple times) using adjlist.",
            REPEAT(result = test_adj(g, &adj), repeat);
         );

    printf("%s", name);
    BENCH("6 go over vertices (multiple times) using inclist, IGRAPH_OTHER.",
            REPEAT(result = test_inc_other(g, &inc), repeat);
         );

    printf("%s", name);
    BENCH("7 go over vertices (multiple times) using inclist, IGRAPH_TO.",
            REPEAT(result = test_inc_to(g, &inc), repeat);
         );

    printf("%s", name);
    BENCH("8 go over edges using inclist, don't retrieve vertex.",
            REPEAT(result = test_inc_nop(g, &inc), repeat);
         );

    printf("%s", name);
    BENCH("9 go over edges and vertices using adjlist and inclist.",
            REPEAT(result = test_inc_adj(g, &inc, &adj), repeat);
         );

    igraph_adjlist_destroy(&adj);
    igraph_inclist_destroy(&inc);

    printf("%s", name);
    BENCH("10 go over edges and vertices using graph internals directly.",
            REPEAT(result = test_direct(g), repeat);
         );

    printf("%s", name);
    BENCH("11 init/destroy interleaved incadjlist.",
            REPEAT(
              do {
                  igraph_incadjlist_inter_init(g, &incadj_inter, IGRAPH_ALL);
                  igraph_incadjlist_inter_destroy(&incadj_inter);
              } while (0),
              repeat);
         );

    igraph_incadjlist_inter_init(g, &incadj_inter, IGRAPH_ALL);

    printf("%s", name);
    BENCH("12 go over edges and vertices using interleaved incadjlist.",
            REPEAT(result = test_incadj_inter(g, &incadj_inter), repeat);
         );

    igraph_incadjlist_inter_destroy(&incadj_inter);

    printf("%s", name);
    BENCH("13 init/destroy incadjlist with two vectors.",
            REPEAT(
              do {
                  igraph_incadjlist_sep_init(g, &incadj_sep, IGRAPH_ALL);
                  igraph_incadjlist_sep_destroy(&incadj_sep);
              } while (0), repeat);
         );

    igraph_incadjlist_sep_init(g, &incadj_sep, IGRAPH_ALL);

    printf("%s", name);
    BENCH("14 go over edges and vertices using incadjlist.",
            REPEAT(result = test_incadj_sep(g, &incadj_sep), repeat);
         );

    igraph_incadjlist_sep_destroy(&incadj_sep);
}

int main(void) {
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
