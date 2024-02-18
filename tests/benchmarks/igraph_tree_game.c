
#include <igraph.h>

#include "bench.h"

#define TOSTR1(x) #x
#define TOSTR(x) TOSTR1(x)

int main(void) {
    igraph_t g;

    igraph_rng_seed(igraph_rng_default(), 137);
    BENCH_INIT();

#define VCOUNT 100
#define REP 100000

    BENCH(" 1 vcount=" TOSTR(VCOUNT) ", Prufer, " TOSTR(REP) "x",
          REPEAT(igraph_tree_game(&g, VCOUNT, IGRAPH_UNDIRECTED, IGRAPH_RANDOM_TREE_PRUFER), REP);
          );
    igraph_destroy(&g);

    BENCH(" 2 vcount=" TOSTR(VCOUNT) ", LERW, " TOSTR(REP) "x",
          REPEAT(igraph_tree_game(&g, VCOUNT, IGRAPH_UNDIRECTED, IGRAPH_RANDOM_TREE_LERW), REP);
          );
    igraph_destroy(&g);

#undef VCOUNT
#undef DENS
#undef REP

    printf("\n");

#define VCOUNT 1000
#define REP 10000

    BENCH(" 1 vcount=" TOSTR(VCOUNT) ", Prufer, " TOSTR(REP) "x",
          REPEAT(igraph_tree_game(&g, VCOUNT, IGRAPH_UNDIRECTED, IGRAPH_RANDOM_TREE_PRUFER), REP);
    );
    igraph_destroy(&g);

    BENCH(" 2 vcount=" TOSTR(VCOUNT) ", LERW, " TOSTR(REP) "x",
          REPEAT(igraph_tree_game(&g, VCOUNT, IGRAPH_UNDIRECTED, IGRAPH_RANDOM_TREE_LERW), REP);
    );
    igraph_destroy(&g);

#undef VCOUNT
#undef DENS
#undef REP

    printf("\n");

#define VCOUNT 10000
#define REP 1000

    BENCH(" 3 vcount=" TOSTR(VCOUNT) ", Prufer, " TOSTR(REP) "x",
          REPEAT(igraph_tree_game(&g, VCOUNT, IGRAPH_UNDIRECTED, IGRAPH_RANDOM_TREE_PRUFER), REP);
          );
    igraph_destroy(&g);

    BENCH(" 4 vcount=" TOSTR(VCOUNT) ", LERW, " TOSTR(REP) "x",
          REPEAT(igraph_tree_game(&g, VCOUNT, IGRAPH_UNDIRECTED, IGRAPH_RANDOM_TREE_LERW), REP);
          );
    igraph_destroy(&g);

#undef VCOUNT
#undef DENS
#undef REP

    printf("\n");

#define VCOUNT 100000
#define REP 100

    BENCH(" 3 vcount=" TOSTR(VCOUNT) ", Prufer, " TOSTR(REP) "x",
          REPEAT(igraph_tree_game(&g, VCOUNT, IGRAPH_UNDIRECTED, IGRAPH_RANDOM_TREE_PRUFER), REP);
          );
    igraph_destroy(&g);

    BENCH(" 4 vcount=" TOSTR(VCOUNT) ", LERW, " TOSTR(REP) "x",
          REPEAT(igraph_tree_game(&g, VCOUNT, IGRAPH_UNDIRECTED, IGRAPH_RANDOM_TREE_LERW), REP);
          );
    igraph_destroy(&g);

#undef VCOUNT
#undef DENS
#undef REP

    return 0;
}
