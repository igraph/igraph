/*
   igraph library.
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

#include "bench.h"

void match(const igraph_t *graph,
           const igraph_t patt[], igraph_int_t n,
           igraph_vector_int_list_t *maps) {

    igraph_vector_int_list_clear(maps);
    for (igraph_int_t i=0; i < n; i++) {
        igraph_subisomorphic_lad(&patt[i], graph, NULL, NULL, NULL, maps, false);
    }
}

#define NP 8

int main(void) {
    igraph_t graph;
    igraph_t patt[NP];
    igraph_vector_int_list_t maps;

    BENCH_INIT();

    igraph_vector_int_list_init(&maps, 0);

    for (igraph_int_t i=0; i < NP; i++) {
        igraph_ring(&patt[i], i+1, IGRAPH_DIRECTED, false, true);
    }

    igraph_kautz(&graph, 3, 3);
    BENCH("1 Kautz(3,3) 10x", REPEAT(match(&graph, patt, NP, &maps), 10));
    igraph_destroy(&graph);

    igraph_kautz(&graph, 3, 4);
    BENCH("2 Kautz(3,4) 3x", REPEAT(match(&graph, patt, NP, &maps), 3));
    igraph_destroy(&graph);

    igraph_kautz(&graph, 4, 3);
    BENCH("3 Kautz(4,3) 1x", REPEAT(match(&graph, patt, NP, &maps), 1));
    igraph_destroy(&graph);

    for (igraph_int_t i=0; i < NP; i++) {
        igraph_destroy(&patt[i]);
    }

    igraph_vector_int_list_destroy(&maps);

    return 0;
}
