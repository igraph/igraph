/*
   IGraph library.
   Copyright (C) 2011-2022  The igraph development team <igraph@igraph.org>

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

#define N 1000

igraph_bool_t has_loops(const igraph_t *graph) {
    igraph_bool_t res;
    igraph_has_loop(graph, &res);
    return res;
}

igraph_bool_t has_multiple(const igraph_t *graph) {
    igraph_bool_t res;
    igraph_has_multiple(graph, &res);
    return res;
}

#define ERR() do {              \
        printf("Seed: %" IGRAPH_PRIu "\n", seed);           \
        print_graph(&ws); \
    } while (0)

#define SEED() do {                         \
        seed=igraph_rng_get_integer(igraph_rng_default(), 1, 10000);        \
        igraph_rng_seed(igraph_rng_default(), seed);                \
    } while (0)

int main(void) {

    igraph_t ws;
    igraph_bool_t sim, seen_loops, seen_multiple;
    igraph_integer_t i;
    igraph_uint_t seed = 1305473657;

    igraph_rng_seed(igraph_rng_default(), seed);

    /* No loops, no multiple edges */
    for (i = 0; i < N; i++) {
        SEED();
        igraph_watts_strogatz_game(&ws, /*dim=*/ 1, /*size=*/ 5, /*nei=*/ 1,
                                   /*p=*/ 0.5, IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE);
        igraph_is_simple(&ws, &sim);
        if (!sim) {
            ERR();
            return 1;
        }
        if (has_loops(&ws)) {
            ERR();
            return 1;
        }
        if (has_multiple(&ws)) {
            ERR();
            return 2;
        }
        igraph_destroy(&ws);
    }

    /* No loops, multiple edges possible */
    seen_multiple = 0;
    for (i = 0; i < N; i++) {
        SEED();
        igraph_watts_strogatz_game(&ws, /*dim=*/ 1, /*size=*/ 5, /*nei=*/ 1,
                                   /*p=*/ 0.5, IGRAPH_NO_LOOPS, IGRAPH_MULTIPLE);
        if (has_loops(&ws)) {
            ERR();
            return 3;
        }
        seen_multiple = seen_multiple || has_multiple(&ws);
        igraph_destroy(&ws);
    }
    /* This might actually happen */
    /* if (!seen_multiple) { return 4; } */

    /* Loops possible, no multiple edges */
    seen_loops = 0;
    for (i = 0; i < N; i++) {
        SEED();
        igraph_watts_strogatz_game(&ws, /*dim=*/ 1, /*size=*/ 5, /*nei=*/ 1,
                                   /*p=*/ 0.5, IGRAPH_LOOPS, IGRAPH_NO_MULTIPLE);
        if (has_multiple(&ws)) {
            return 5;
        }
        seen_loops = seen_loops || has_loops(&ws);
        igraph_destroy(&ws);
    }
    /* This might actually happen */
    /* if (!seen_loops) { return 6; } */

    /* Both loops and multiple edges are possible */
    for (i = 0; i < N; i++) {
        SEED();
        igraph_watts_strogatz_game(&ws, /*dim=*/ 1, /*size=*/ 5, /*nei=*/ 1,
                                   /*p=*/ 0.5, IGRAPH_LOOPS, IGRAPH_MULTIPLE);
        seen_loops = seen_loops || has_loops(&ws);
        seen_multiple = seen_multiple || has_multiple(&ws);
        igraph_destroy(&ws);
    }
    /* This might actually happen */
    /* if (!seen_loops) { return 7; } */
    /* if (!seen_multiple) { return 8; }   */

    VERIFY_FINALLY_STACK();

    return 0;
}
