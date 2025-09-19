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

void print_params(const igraph_t *graph, igraph_neimode_t mode, igraph_loops_t loops) {
    printf(igraph_is_directed(graph) ? "directed;   " : "undirected; ");
    switch (mode) {
        case IGRAPH_ALL: printf("ALL; "); break;
        case IGRAPH_OUT: printf("OUT; "); break;
        case IGRAPH_IN : printf("IN;  "); break;
    }
    switch (loops) {
        case IGRAPH_NO_LOOPS:    printf("NO_LOOPS;    "); break;
        case IGRAPH_LOOPS_ONCE:  printf("LOOPS_ONCE;  "); break;
        case IGRAPH_LOOPS_TWICE: printf("LOOPS_TWICE; "); break;
    }
}

void print_multiple(igraph_bool_t multiple) {
    printf(multiple ? "MULTIPLE;    " : "NO_MULTIPLE; ");
}

/* Print adjacency list based on igraph_neighbors() */
void print_adj1(const igraph_t *graph, igraph_neimode_t mode, igraph_loops_t loops, igraph_bool_t multiple) {
    igraph_int_t vcount = igraph_vcount(graph);
    igraph_vector_int_t neis;

    print_params(graph, mode, loops);
    print_multiple(multiple);
    printf("\n");

    igraph_vector_int_init(&neis, 0);
    for (igraph_int_t v=0; v < vcount; v++) {
        printf("%3" IGRAPH_PRId ": ", v);
        igraph_neighbors(graph, &neis, v, mode, loops, multiple);
        print_vector_int(&neis);
    }
    igraph_vector_int_destroy(&neis);
}

/* Print adjacency list based on igraph_adjlist_init() */
void print_adj2(const igraph_t *graph, igraph_neimode_t mode, igraph_loops_t loops, igraph_bool_t multiple) {
    igraph_int_t vcount = igraph_vcount(graph);
    igraph_adjlist_t al;

    print_params(graph, mode, loops);
    print_multiple(multiple);
    printf("\n");

    igraph_adjlist_init(graph, &al, mode, loops, multiple);
    for (igraph_int_t v=0; v < vcount; v++) {
        printf("%3" IGRAPH_PRId ": ", v);
        print_vector_int(igraph_adjlist_get(&al, v));
    }
    igraph_adjlist_destroy(&al);
}

/* Print incidence list based on igraph_incident() */
void print_inc1(const igraph_t *graph, igraph_neimode_t mode, igraph_loops_t loops) {
    igraph_int_t vcount = igraph_vcount(graph);
    igraph_vector_int_t incs;

    print_params(graph, mode, loops);
    printf("\n");

    igraph_vector_int_init(&incs, 0);
    for (igraph_int_t v=0; v < vcount; v++) {
        printf("%3" IGRAPH_PRId ": ", v);
        igraph_incident(graph, &incs, v, mode, loops);
        print_vector_int(&incs);
    }
    igraph_vector_int_destroy(&incs);
}

/* Print incidence list based on igraph_inclist_init() */
void print_inc2(const igraph_t *graph, igraph_neimode_t mode, igraph_loops_t loops) {
    igraph_int_t vcount = igraph_vcount(graph);
    igraph_inclist_t il;

    print_params(graph, mode, loops);
    printf("\n");

    igraph_inclist_init(graph, &il, mode, loops);
    for (igraph_int_t v=0; v < vcount; v++) {
        printf("%3" IGRAPH_PRId ": ", v);
        print_vector_int(igraph_inclist_get(&il, v));
    }
    igraph_inclist_destroy(&il);
}

/* Verifies that igraph_neighbors() and igraph_incident() produce results in the same order. */
void verify_ordering1(const igraph_t *graph, igraph_neimode_t mode, igraph_loops_t loops) {
    igraph_vector_int_t neis, incs;
    igraph_int_t vcount = igraph_vcount(graph);

    if (! igraph_is_directed(graph)) {
        mode = IGRAPH_ALL;
    }

    igraph_vector_int_init(&neis, 0);
    igraph_vector_int_init(&incs, 0);

    for (igraph_int_t v=0; v < vcount; v++) {
        igraph_neighbors(graph, &neis, v, mode, loops, IGRAPH_MULTIPLE);
        igraph_incident(graph, &incs, v, mode, loops);

        igraph_int_t n = igraph_vector_int_size(&neis);
        IGRAPH_ASSERT(igraph_vector_int_size(&incs) == n);

        for (igraph_int_t i=0; i < n; i++) {
            igraph_int_t u1, u2;

            u1 = VECTOR(neis)[i];
            switch (mode) {
            case IGRAPH_ALL: u2 = IGRAPH_OTHER(graph, VECTOR(incs)[i], v); break;
            case IGRAPH_OUT: u2 = IGRAPH_TO(graph, VECTOR(incs)[i]); break;
            case IGRAPH_IN:  u2 = IGRAPH_FROM(graph, VECTOR(incs)[i]); break;
            }

            if (u1 != u2) {
                IGRAPH_FATALF(
                    "For vertex %d, the %dth neighbor vertex between adj (%d) and inc (%d) differs.",
                    (int) v, (int) i, (int) u1, (int) u2
                );
            }
        }
    }

    igraph_vector_int_destroy(&incs);
    igraph_vector_int_destroy(&neis);
}

/* Verifies that igraph_adjlist_init() and igraph_inclist_init() produce results in the same order. */
void verify_ordering2(const igraph_t *graph, igraph_neimode_t mode, igraph_loops_t loops) {
    igraph_adjlist_t al;
    igraph_inclist_t il;

    igraph_int_t vcount = igraph_vcount(graph);

    if (! igraph_is_directed(graph)) {
        mode = IGRAPH_ALL;
    }

    igraph_adjlist_init(graph, &al, mode, loops, IGRAPH_MULTIPLE);
    igraph_inclist_init(graph, &il, mode, loops);

    for (igraph_int_t v=0; v < vcount; v++) {
        igraph_vector_int_t *neis = igraph_adjlist_get(&al, v);
        igraph_vector_int_t *incs = igraph_inclist_get(&il, v);

        igraph_int_t n = igraph_vector_int_size(neis);
        IGRAPH_ASSERT(igraph_vector_int_size(incs) == n);

        for (igraph_int_t i=0; i < n; i++) {
            igraph_int_t u1, u2;

            u1 = VECTOR(*neis)[i];
            switch (mode) {
            case IGRAPH_ALL: u2 = IGRAPH_OTHER(graph, VECTOR(*incs)[i], v); break;
            case IGRAPH_OUT: u2 = IGRAPH_TO(graph, VECTOR(*incs)[i]); break;
            case IGRAPH_IN:  u2 = IGRAPH_FROM(graph, VECTOR(*incs)[i]); break;
            }

            if (u1 != u2) {
                IGRAPH_FATALF(
                    "For vertex %d, the %dth neighbor vertex between adj (%d) and inc (%d) differs.",
                    (int) v, (int) i, (int) u1, (int) u2
                );
            }
        }
    }

    igraph_adjlist_destroy(&al);
    igraph_inclist_destroy(&il);
}

int main(void) {
    igraph_t dg, ug;

    igraph_small(&dg, 7, IGRAPH_DIRECTED,
              // 0    1    2    3    4    5    6    7    8    9    10   11   12
                 4,4, 0,3, 0,3, 2,3, 1,4, 2,4, 3,4, 3,2, 4,4, 2,3, 1,1, 6,6, 6,6,
                 -1);

    igraph_copy(&ug, &dg);
    igraph_to_undirected(&ug, IGRAPH_TO_UNDIRECTED_EACH, NULL);

    printf("UNDIRECTED\n\n");

    printf("igraph_neighbors()\n");
    print_adj1(&ug, IGRAPH_ALL, IGRAPH_NO_LOOPS, IGRAPH_MULTIPLE);
    print_adj1(&ug, IGRAPH_ALL, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE);
    print_adj1(&ug, IGRAPH_ALL, IGRAPH_LOOPS_TWICE, IGRAPH_MULTIPLE);

    printf("\nigraph_adjlist_init()\n");
    print_adj2(&ug, IGRAPH_ALL, IGRAPH_NO_LOOPS, IGRAPH_MULTIPLE);
    print_adj2(&ug, IGRAPH_ALL, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE);
    print_adj2(&ug, IGRAPH_ALL, IGRAPH_LOOPS_TWICE, IGRAPH_MULTIPLE);

    printf("\nigraph_incident()\n");
    print_inc1(&ug, IGRAPH_ALL, IGRAPH_NO_LOOPS);
    print_inc1(&ug, IGRAPH_ALL, IGRAPH_LOOPS_ONCE);
    print_inc1(&ug, IGRAPH_ALL, IGRAPH_LOOPS_TWICE);

    printf("\nigraph_inclist_init()\n");
    print_inc2(&ug, IGRAPH_ALL, IGRAPH_NO_LOOPS);
    print_inc2(&ug, IGRAPH_ALL, IGRAPH_LOOPS_ONCE);
    print_inc2(&ug, IGRAPH_ALL, IGRAPH_LOOPS_TWICE);

    printf("\nDIRECTED ALL\n\n");

    printf("igraph_neighbors()\n");
    print_adj1(&dg, IGRAPH_ALL, IGRAPH_NO_LOOPS, IGRAPH_MULTIPLE);
    print_adj1(&dg, IGRAPH_ALL, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE);
    print_adj1(&dg, IGRAPH_ALL, IGRAPH_LOOPS_TWICE, IGRAPH_MULTIPLE);

    printf("\nigraph_adjlist_init()\n");
    print_adj2(&dg, IGRAPH_ALL, IGRAPH_NO_LOOPS, IGRAPH_MULTIPLE);
    print_adj2(&dg, IGRAPH_ALL, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE);
    print_adj2(&dg, IGRAPH_ALL, IGRAPH_LOOPS_TWICE, IGRAPH_MULTIPLE);

    printf("\nigraph_incident()\n");
    print_inc1(&dg, IGRAPH_ALL, IGRAPH_NO_LOOPS);
    print_inc1(&dg, IGRAPH_ALL, IGRAPH_LOOPS_ONCE);
    print_inc1(&dg, IGRAPH_ALL, IGRAPH_LOOPS_TWICE);

    printf("\nigraph_inclist_init()\n");
    print_inc2(&dg, IGRAPH_ALL, IGRAPH_NO_LOOPS);
    print_inc2(&dg, IGRAPH_ALL, IGRAPH_LOOPS_ONCE);
    print_inc2(&dg, IGRAPH_ALL, IGRAPH_LOOPS_TWICE);

    printf("\nDIRECTED OUT\n\n");

    printf("igraph_neighbors()\n");
    print_adj1(&dg, IGRAPH_OUT, IGRAPH_NO_LOOPS, IGRAPH_MULTIPLE);
    print_adj1(&dg, IGRAPH_OUT, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE);

    printf("\nigraph_adjlist_init()\n");
    print_adj2(&dg, IGRAPH_OUT, IGRAPH_NO_LOOPS, IGRAPH_MULTIPLE);
    print_adj2(&dg, IGRAPH_OUT, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE);

    printf("\nigraph_incident()\n");
    print_inc1(&dg, IGRAPH_OUT, IGRAPH_NO_LOOPS);
    print_inc1(&dg, IGRAPH_OUT, IGRAPH_LOOPS_ONCE);

    printf("\nigraph_inclist_init()\n");
    print_inc2(&dg, IGRAPH_OUT, IGRAPH_NO_LOOPS);
    print_inc2(&dg, IGRAPH_OUT, IGRAPH_LOOPS_ONCE);

    /* Verify that neighbour lists (adjlist) and incident edge lists (inclist)
     * have the same ordering with all parameter combinations. */

    verify_ordering1(&ug, IGRAPH_ALL, IGRAPH_NO_LOOPS);
    verify_ordering1(&ug, IGRAPH_ALL, IGRAPH_LOOPS_ONCE);
    verify_ordering1(&ug, IGRAPH_ALL, IGRAPH_LOOPS_TWICE);

    verify_ordering1(&dg, IGRAPH_ALL, IGRAPH_NO_LOOPS);
    verify_ordering1(&dg, IGRAPH_ALL, IGRAPH_LOOPS_ONCE);
    verify_ordering1(&dg, IGRAPH_ALL, IGRAPH_LOOPS_TWICE);

    verify_ordering1(&dg, IGRAPH_OUT, IGRAPH_NO_LOOPS);
    verify_ordering1(&dg, IGRAPH_OUT, IGRAPH_LOOPS_ONCE);

    verify_ordering1(&dg, IGRAPH_IN, IGRAPH_NO_LOOPS);
    verify_ordering1(&dg, IGRAPH_IN, IGRAPH_LOOPS_ONCE);

    verify_ordering2(&ug, IGRAPH_ALL, IGRAPH_NO_LOOPS);
    verify_ordering2(&ug, IGRAPH_ALL, IGRAPH_LOOPS_ONCE);
    verify_ordering2(&ug, IGRAPH_ALL, IGRAPH_LOOPS_TWICE);

    verify_ordering2(&dg, IGRAPH_ALL, IGRAPH_NO_LOOPS);
    verify_ordering2(&dg, IGRAPH_ALL, IGRAPH_LOOPS_ONCE);
    verify_ordering2(&dg, IGRAPH_ALL, IGRAPH_LOOPS_TWICE);

    verify_ordering2(&dg, IGRAPH_OUT, IGRAPH_NO_LOOPS);
    verify_ordering2(&dg, IGRAPH_OUT, IGRAPH_LOOPS_ONCE);

    verify_ordering2(&dg, IGRAPH_IN, IGRAPH_NO_LOOPS);
    verify_ordering2(&dg, IGRAPH_IN, IGRAPH_LOOPS_ONCE);

    igraph_destroy(&ug);
    igraph_destroy(&dg);

    VERIFY_FINALLY_STACK();

    return 0;
}
