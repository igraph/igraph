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

igraph_integer_t test_adj(igraph_t *g)
{
    igraph_integer_t dummy = 0;
    igraph_adjlist_t adj;

    igraph_adjlist_init(g, &adj, IGRAPH_ALL, IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE);
    for (int i = 0; i < igraph_vcount(g); i++) {
        igraph_vector_int_t *neis = igraph_adjlist_get(&adj, i);
        igraph_integer_t nneis = igraph_vector_int_size(neis);
        for (int j = 0; j < nneis; j++) {
            igraph_integer_t neighbor = VECTOR(*neis)[j];
            dummy += neighbor;
        }
    }
    igraph_adjlist_destroy(&adj);
    return dummy;
}


igraph_integer_t test_inc(igraph_t *g)
{
    igraph_integer_t dummy = 0;
    igraph_inclist_t inc;

    igraph_inclist_init(g, &inc, IGRAPH_ALL, IGRAPH_NO_LOOPS);
    for (int i = 0; i < igraph_vcount(g); i++) {
        igraph_vector_int_t *neis = igraph_inclist_get(&inc, i);
        igraph_integer_t nneis = igraph_vector_int_size(neis);
        for (int j = 0; j < nneis; j++) {
            igraph_integer_t edge = VECTOR(*neis)[j];
            igraph_integer_t neighbor = IGRAPH_OTHER(g, edge, i);
            dummy += neighbor;
        }
    }
    igraph_inclist_destroy(&inc);
    return dummy;
}


int main() {
    igraph_t g;

    igraph_full(&g, 10000, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    BENCH_INIT();

    BENCH(" 1 go over vertices (multiple times) using adjlist",
            test_adj(&g);
         );

    BENCH(" 2 go over vertices (multiple times) using inclist",
            test_inc(&g);
         );
    return 0;
}
