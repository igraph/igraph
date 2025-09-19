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

#define TOSTR1(x) #x
#define TOSTR(x) TOSTR1(x)

void tree_game(igraph_int_t n, igraph_bool_t directed, igraph_random_tree_t method) {
    igraph_t g;
    igraph_tree_game(&g, n, directed, method);
    igraph_destroy(&g);
}

int main(void) {

    igraph_rng_seed(igraph_rng_default(), 137);
    BENCH_INIT();

#define VCOUNT 100
#define REP 100000

    BENCH(" 1 vcount=" TOSTR(VCOUNT) ", Prufer, " TOSTR(REP) "x",
          REPEAT(tree_game(VCOUNT, IGRAPH_UNDIRECTED, IGRAPH_RANDOM_TREE_PRUFER), REP);
    );

    BENCH(" 2 vcount=" TOSTR(VCOUNT) ", LERW, " TOSTR(REP) "x",
          REPEAT(tree_game(VCOUNT, IGRAPH_UNDIRECTED, IGRAPH_RANDOM_TREE_LERW), REP);
    );

#undef VCOUNT
#undef DENS
#undef REP

    printf("\n");

#define VCOUNT 1000
#define REP 10000

    BENCH(" 1 vcount=" TOSTR(VCOUNT) ", Prufer, " TOSTR(REP) "x",
          REPEAT(tree_game(VCOUNT, IGRAPH_UNDIRECTED, IGRAPH_RANDOM_TREE_PRUFER), REP);
    );

    BENCH(" 2 vcount=" TOSTR(VCOUNT) ", LERW, " TOSTR(REP) "x",
          REPEAT(tree_game(VCOUNT, IGRAPH_UNDIRECTED, IGRAPH_RANDOM_TREE_LERW), REP);
    );

#undef VCOUNT
#undef DENS
#undef REP

    printf("\n");

#define VCOUNT 10000
#define REP 1000

    BENCH(" 3 vcount=" TOSTR(VCOUNT) ", Prufer, " TOSTR(REP) "x",
          REPEAT(tree_game(VCOUNT, IGRAPH_UNDIRECTED, IGRAPH_RANDOM_TREE_PRUFER), REP);
    );

    BENCH(" 4 vcount=" TOSTR(VCOUNT) ", LERW, " TOSTR(REP) "x",
          REPEAT(tree_game(VCOUNT, IGRAPH_UNDIRECTED, IGRAPH_RANDOM_TREE_LERW), REP);
    );

#undef VCOUNT
#undef DENS
#undef REP

    printf("\n");

#define VCOUNT 100000
#define REP 100

    BENCH(" 3 vcount=" TOSTR(VCOUNT) ", Prufer, " TOSTR(REP) "x",
          REPEAT(tree_game(VCOUNT, IGRAPH_UNDIRECTED, IGRAPH_RANDOM_TREE_PRUFER), REP);
    );

    BENCH(" 4 vcount=" TOSTR(VCOUNT) ", LERW, " TOSTR(REP) "x",
          REPEAT(tree_game(VCOUNT, IGRAPH_UNDIRECTED, IGRAPH_RANDOM_TREE_LERW), REP);
    );

#undef VCOUNT
#undef DENS
#undef REP

    return 0;
}
