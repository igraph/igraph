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

#define BENCH_BLOCK() \
    BENCH(" 1 vc=" TOSTR(VCOUNT) ", ec=" TOSTR(ECOUNT) ", k=3, all, " TOSTR(REP) "x", \
          REPEAT(igraph_ecc(&g, &ecc, igraph_ess_all(IGRAPH_EDGEORDER_ID), 3, false, true), REP); \
    ); \
    BENCH(" 2 vc=" TOSTR(VCOUNT) ", ec=" TOSTR(ECOUNT) ", k=3, subset: all, " TOSTR(REP) "x", \
          REPEAT(igraph_ecc(&g, &ecc, igraph_ess_range(0, igraph_ecount(&g)), 3, false, true), REP); \
    ); \
    \
    igraph_random_sample(&eids, 0, igraph_ecount(&g)-1, SS); \
    \
    BENCH(" 3 vc=" TOSTR(VCOUNT) ", ec=" TOSTR(ECOUNT) ", k=3, subset: " TOSTR(SS) ", " TOSTR(SREP) "x", \
          REPEAT(igraph_ecc(&g, &ecc, igraph_ess_vector(&eids), 3, false, true), SREP); \
    ); \
    \
    BENCH(" 4 vc=" TOSTR(VCOUNT) ", ec=" TOSTR(ECOUNT) ", k=4, all, " TOSTR(REP) "x", \
          REPEAT(igraph_ecc(&g, &ecc, igraph_ess_all(IGRAPH_EDGEORDER_ID), 4, false, true), REP); \
    ); \
    \
    BENCH(" 4 vc=" TOSTR(VCOUNT) ", ec=" TOSTR(ECOUNT) ", k=4, subset: all, " TOSTR(REP) "x", \
          REPEAT(igraph_ecc(&g, &ecc, igraph_ess_range(0, igraph_ecount(&g)), 4, false, true), REP); \
    ); \
    BENCH(" 5 vc=" TOSTR(VCOUNT) ", ec=" TOSTR(ECOUNT) ", k=4, subset: " TOSTR(SS) ", " TOSTR(SREP) "x", \
          REPEAT(igraph_ecc(&g, &ecc, igraph_ess_vector(&eids), 4, false, true), SREP); \
    );

int main(void) {
    igraph_t g;
    igraph_vector_t ecc;
    igraph_vector_int_t eids;

    igraph_rng_seed(igraph_rng_default(), 137);
    BENCH_INIT();

    igraph_vector_init(&ecc, 0);
    igraph_vector_int_init(&eids, 0);

    printf("Erdos-Renyi GNM:\n\n");

#define VCOUNT 100
#define ECOUNT 5000
#define REP 100
#define SS 100
#define SREP 1000

    igraph_erdos_renyi_game_gnm(&g, VCOUNT, ECOUNT, IGRAPH_DIRECTED, IGRAPH_SIMPLE_SW, IGRAPH_EDGE_UNLABELED);
    igraph_vector_resize(&ecc, igraph_ecount(&g));

    BENCH_BLOCK()

    igraph_destroy(&g);

#undef VCOUNT
#undef ECOUNT
#undef REP
#undef SS
#undef SREP

    printf("\n");

#define VCOUNT 1000
#define ECOUNT 5000
#define REP 100
#define SS 100
#define SREP 1000

    igraph_erdos_renyi_game_gnm(&g, VCOUNT, ECOUNT, IGRAPH_DIRECTED, IGRAPH_SIMPLE_SW, IGRAPH_EDGE_UNLABELED);
    igraph_vector_resize(&ecc, igraph_ecount(&g));

    BENCH_BLOCK()

    igraph_destroy(&g);

#undef VCOUNT
#undef ECOUNT
#undef REP
#undef SS
#undef SREP

    printf("\n");

#define VCOUNT 10000
#define ECOUNT 10000
#define REP 100
#define SS 100
#define SREP 1000

    igraph_erdos_renyi_game_gnm(&g, VCOUNT, ECOUNT, IGRAPH_DIRECTED, IGRAPH_SIMPLE_SW, IGRAPH_EDGE_UNLABELED);
    igraph_vector_resize(&ecc, igraph_ecount(&g));

    BENCH_BLOCK()

    igraph_destroy(&g);

#undef VCOUNT
#undef ECOUNT
#undef REP
#undef SS
#undef SREP

    printf("\n");

#define VCOUNT 10000
#define ECOUNT 50000
#define REP 100
#define SS 100
#define SREP 1000

    igraph_erdos_renyi_game_gnm(&g, VCOUNT, ECOUNT, IGRAPH_DIRECTED, IGRAPH_SIMPLE_SW, IGRAPH_EDGE_UNLABELED);
    igraph_vector_resize(&ecc, igraph_ecount(&g));

    BENCH_BLOCK()

    igraph_destroy(&g);

#undef VCOUNT
#undef ECOUNT
#undef REP
#undef SS
#undef SREP

    printf("\n");

#define VCOUNT 10000
#define ECOUNT 50000
#define REP 100
#define SS 100
#define SREP 1000

    printf("Barabasi:\n\n");

    igraph_barabasi_game(&g,
                         VCOUNT, 1, ECOUNT / VCOUNT, NULL, true, 0,
                         IGRAPH_UNDIRECTED, IGRAPH_BARABASI_PSUMTREE_MULTIPLE, NULL);
    igraph_vector_resize(&ecc, igraph_ecount(&g));

    BENCH_BLOCK()

    igraph_destroy(&g);

#undef VCOUNT
#undef ECOUNT
#undef REP
#undef SS
#undef SREP

    igraph_vector_int_destroy(&eids);
    igraph_vector_destroy(&ecc);

    return 0;
}
