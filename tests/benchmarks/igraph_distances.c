
#include <igraph.h>

#include "bench.h"

#define TOSTR1(x) #x
#define TOSTR(x) TOSTR1(x)

int main(void) {
    igraph_t g;
    igraph_matrix_t res;
    igraph_vector_t weights;

    igraph_rng_seed(igraph_rng_default(), 137);
    BENCH_INIT();

    igraph_matrix_init(&res, 0, 0);
    igraph_vector_init(&weights, 0);

#define VCOUNT 100
#define DENS 0.5
#define REP 200

    igraph_erdos_renyi_game_gnp(&g, VCOUNT, DENS, IGRAPH_DIRECTED, IGRAPH_NO_LOOPS);
    igraph_matrix_resize(&res, igraph_vcount(&g), igraph_vcount(&g));
    igraph_vector_resize(&weights, igraph_ecount(&g));

    RNG_BEGIN();
    for (igraph_integer_t i=0; i < igraph_ecount(&g); i++) {
        VECTOR(weights)[i] = RNG_EXP(1);
    }
    RNG_END();

    BENCH(" 1 vcount=" TOSTR(VCOUNT) ", p=" TOSTR(DENS) ", Unweighted, " TOSTR(REP) "x",
          REPEAT(igraph_distances(&g, &res, igraph_vss_all(), igraph_vss_all(), IGRAPH_OUT), REP);
    );
    BENCH(" 2 vcount=" TOSTR(VCOUNT) ", p=" TOSTR(DENS) ", Dijkstra, " TOSTR(REP) "x",
          REPEAT(igraph_distances_dijkstra(&g, &res, igraph_vss_all(), igraph_vss_all(), &weights, IGRAPH_OUT), REP)
    );
    BENCH(" 3 vcount=" TOSTR(VCOUNT) ", p=" TOSTR(DENS) ", Unweighted Floyd-Warshall, " TOSTR(REP) "x",
          REPEAT(igraph_distances_floyd_warshall(&g, &res, igraph_vss_all(), igraph_vss_all(), NULL, IGRAPH_OUT, IGRAPH_FLOYD_WARSHALL_ORIGINAL), REP)
    );
    BENCH(" 4 vcount=" TOSTR(VCOUNT) ", p=" TOSTR(DENS) ", Unweighted Floyd-Warshall-tree-speedup, " TOSTR(REP) "x",
          REPEAT(igraph_distances_floyd_warshall(&g, &res, igraph_vss_all(), igraph_vss_all(), NULL, IGRAPH_OUT, IGRAPH_FLOYD_WARSHALL_TREE), REP)
    );
    BENCH(" 5 vcount=" TOSTR(VCOUNT) ", p=" TOSTR(DENS) ", Floyd-Warshall, " TOSTR(REP) "x",
          REPEAT(igraph_distances_floyd_warshall(&g, &res, igraph_vss_all(), igraph_vss_all(), &weights, IGRAPH_OUT, IGRAPH_FLOYD_WARSHALL_ORIGINAL), REP)
    );
    BENCH(" 6 vcount=" TOSTR(VCOUNT) ", p=" TOSTR(DENS) ", Floyd-Warshall-tree-speedup, " TOSTR(REP) "x",
          REPEAT(igraph_distances_floyd_warshall(&g, &res, igraph_vss_all(), igraph_vss_all(), &weights, IGRAPH_OUT, IGRAPH_FLOYD_WARSHALL_TREE), REP)
    );

    RNG_BEGIN();
    for (igraph_integer_t i=0; i < trunc(0.005 * igraph_ecount(&g)); i++) {
        /* For reproducibility, do not write two RNG_...() calls on the same line
         * as the C language does not guarantee any evaluation order between them. */
        igraph_real_t w = RNG_UNIF(-0.05, 0);
        VECTOR(weights)[RNG_INTEGER(0, igraph_ecount(&g)-1)] = w;
    }
    RNG_END();

    BENCH(" 7 vcount=" TOSTR(VCOUNT) ", p=" TOSTR(DENS) ", Bellman-Ford (negative), " TOSTR(REP) "x",
          REPEAT(igraph_distances_bellman_ford(&g, &res, igraph_vss_all(), igraph_vss_all(), &weights, IGRAPH_OUT), REP)
    );
    BENCH(" 8 vcount=" TOSTR(VCOUNT) ", p=" TOSTR(DENS) ", Johnson (negative), " TOSTR(REP) "x",
          REPEAT(igraph_distances_johnson(&g, &res, igraph_vss_all(), igraph_vss_all(), &weights), REP)
    );
    BENCH(" 9 vcount=" TOSTR(VCOUNT) ", p=" TOSTR(DENS) ", Floyd-Warshall (negative), " TOSTR(REP) "x",
          REPEAT(igraph_distances_floyd_warshall(&g, &res, igraph_vss_all(), igraph_vss_all(), &weights, IGRAPH_OUT, IGRAPH_FLOYD_WARSHALL_ORIGINAL), REP)
    );
    BENCH(" 10 vcount=" TOSTR(VCOUNT) ", p=" TOSTR(DENS) ", Floyd-Warshall-tree-speedup (negative), " TOSTR(REP) "x",
          REPEAT(igraph_distances_floyd_warshall(&g, &res, igraph_vss_all(), igraph_vss_all(), &weights, IGRAPH_OUT, IGRAPH_FLOYD_WARSHALL_TREE), REP)
    );

    igraph_destroy(&g);

#undef VCOUNT
#undef DENS
#undef REP

    printf("\n");

#define VCOUNT 500
#define DENS 0.1
#define REP 10

    igraph_erdos_renyi_game_gnp(&g, VCOUNT, DENS, IGRAPH_DIRECTED, IGRAPH_NO_LOOPS);
    igraph_matrix_resize(&res, igraph_vcount(&g), igraph_vcount(&g));
    igraph_vector_resize(&weights, igraph_ecount(&g));

    RNG_BEGIN();
    for (igraph_integer_t i=0; i < igraph_ecount(&g); i++) {
        VECTOR(weights)[i] = RNG_EXP(1);
    }
    RNG_END();

    BENCH(" 1 vcount=" TOSTR(VCOUNT) ", p=" TOSTR(DENS) ", Unweighted, " TOSTR(REP) "x",
          REPEAT(igraph_distances(&g, &res, igraph_vss_all(), igraph_vss_all(), IGRAPH_OUT), REP);
    );
    BENCH(" 2 vcount=" TOSTR(VCOUNT) ", p=" TOSTR(DENS) ", Dijkstra, " TOSTR(REP) "x",
          REPEAT(igraph_distances_dijkstra(&g, &res, igraph_vss_all(), igraph_vss_all(), &weights, IGRAPH_OUT), REP)
    );
    BENCH(" 3 vcount=" TOSTR(VCOUNT) ", p=" TOSTR(DENS) ", Unweighted Floyd-Warshall, " TOSTR(REP) "x",
          REPEAT(igraph_distances_floyd_warshall(&g, &res, igraph_vss_all(), igraph_vss_all(), NULL, IGRAPH_OUT, IGRAPH_FLOYD_WARSHALL_ORIGINAL), REP)
    );
    BENCH(" 4 vcount=" TOSTR(VCOUNT) ", p=" TOSTR(DENS) ", Unweighted Floyd-Warshall-tree-speedup, " TOSTR(REP) "x",
          REPEAT(igraph_distances_floyd_warshall(&g, &res, igraph_vss_all(), igraph_vss_all(), NULL, IGRAPH_OUT, IGRAPH_FLOYD_WARSHALL_TREE), REP)
    );
    BENCH(" 5 vcount=" TOSTR(VCOUNT) ", p=" TOSTR(DENS) ", Floyd-Warshall, " TOSTR(REP) "x",
          REPEAT(igraph_distances_floyd_warshall(&g, &res, igraph_vss_all(), igraph_vss_all(), &weights, IGRAPH_OUT, IGRAPH_FLOYD_WARSHALL_ORIGINAL), REP)
    );
    BENCH(" 6 vcount=" TOSTR(VCOUNT) ", p=" TOSTR(DENS) ", Floyd-Warshall-tree-speedup, " TOSTR(REP) "x",
          REPEAT(igraph_distances_floyd_warshall(&g, &res, igraph_vss_all(), igraph_vss_all(), &weights, IGRAPH_OUT, IGRAPH_FLOYD_WARSHALL_TREE), REP)
    );

    RNG_BEGIN();
    for (igraph_integer_t i=0; i < trunc(0.002 * igraph_ecount(&g)); i++) {
        /* For reproducibility, do not write two RNG_...() calls on the same line
         * as the C language does not guarantee any evaluation order between them. */
        igraph_real_t w = RNG_UNIF(-0.02, 0);
        VECTOR(weights)[RNG_INTEGER(0, igraph_ecount(&g)-1)] = w;
    }
    RNG_END();

    BENCH(" 7 vcount=" TOSTR(VCOUNT) ", p=" TOSTR(DENS) ", Bellman-Ford (negative), " TOSTR(REP) "x",
          REPEAT(igraph_distances_bellman_ford(&g, &res, igraph_vss_all(), igraph_vss_all(), &weights, IGRAPH_OUT), REP)
    );
    BENCH(" 8 vcount=" TOSTR(VCOUNT) ", p=" TOSTR(DENS) ", Johnson (negative), " TOSTR(REP) "x",
          REPEAT(igraph_distances_johnson(&g, &res, igraph_vss_all(), igraph_vss_all(), &weights), REP)
    );
    BENCH(" 9 vcount=" TOSTR(VCOUNT) ", p=" TOSTR(DENS) ", Floyd-Warshall (negative), " TOSTR(REP) "x",
          REPEAT(igraph_distances_floyd_warshall(&g, &res, igraph_vss_all(), igraph_vss_all(), &weights, IGRAPH_OUT, IGRAPH_FLOYD_WARSHALL_ORIGINAL), REP)
    );
    BENCH(" 10 vcount=" TOSTR(VCOUNT) ", p=" TOSTR(DENS) ", Floyd-Warshall-tree-speedup (negative), " TOSTR(REP) "x",
          REPEAT(igraph_distances_floyd_warshall(&g, &res, igraph_vss_all(), igraph_vss_all(), &weights, IGRAPH_OUT, IGRAPH_FLOYD_WARSHALL_TREE), REP)
    );

    igraph_destroy(&g);

#undef VCOUNT
#undef DENS
#undef REP

    printf("\n");

#define VCOUNT 1500
#define DENS 0.02
#define REP 1

    igraph_erdos_renyi_game_gnp(&g, VCOUNT, DENS, IGRAPH_DIRECTED, IGRAPH_NO_LOOPS);
    igraph_matrix_resize(&res, igraph_vcount(&g), igraph_vcount(&g));
    igraph_vector_resize(&weights, igraph_ecount(&g));

    RNG_BEGIN();
    for (igraph_integer_t i=0; i < igraph_ecount(&g); i++) {
        VECTOR(weights)[i] = RNG_EXP(1);
    }
    RNG_END();

    BENCH(" 1 vcount=" TOSTR(VCOUNT) ", p=" TOSTR(DENS) ", Unweighted, " TOSTR(REP) "x",
          REPEAT(igraph_distances(&g, &res, igraph_vss_all(), igraph_vss_all(), IGRAPH_OUT), REP);
    );
    BENCH(" 2 vcount=" TOSTR(VCOUNT) ", p=" TOSTR(DENS) ", Dijkstra, " TOSTR(REP) "x",
          REPEAT(igraph_distances_dijkstra(&g, &res, igraph_vss_all(), igraph_vss_all(), &weights, IGRAPH_OUT), REP)
    );
    BENCH(" 3 vcount=" TOSTR(VCOUNT) ", p=" TOSTR(DENS) ", Unweighted Floyd-Warshall, " TOSTR(REP) "x",
          REPEAT(igraph_distances_floyd_warshall(&g, &res, igraph_vss_all(), igraph_vss_all(), NULL, IGRAPH_OUT, IGRAPH_FLOYD_WARSHALL_ORIGINAL), REP)
    );
    BENCH(" 4 vcount=" TOSTR(VCOUNT) ", p=" TOSTR(DENS) ", Unweighted Floyd-Warshall-tree-speedup, " TOSTR(REP) "x",
          REPEAT(igraph_distances_floyd_warshall(&g, &res, igraph_vss_all(), igraph_vss_all(), NULL, IGRAPH_OUT, IGRAPH_FLOYD_WARSHALL_TREE), REP)
    );
    BENCH(" 5 vcount=" TOSTR(VCOUNT) ", p=" TOSTR(DENS) ", Floyd-Warshall, " TOSTR(REP) "x",
          REPEAT(igraph_distances_floyd_warshall(&g, &res, igraph_vss_all(), igraph_vss_all(), &weights, IGRAPH_OUT, IGRAPH_FLOYD_WARSHALL_ORIGINAL), REP)
    );
    BENCH(" 6 vcount=" TOSTR(VCOUNT) ", p=" TOSTR(DENS) ", Floyd-Warshall-tree-speedup, " TOSTR(REP) "x",
          REPEAT(igraph_distances_floyd_warshall(&g, &res, igraph_vss_all(), igraph_vss_all(), &weights, IGRAPH_OUT, IGRAPH_FLOYD_WARSHALL_TREE), REP)
    );

    RNG_BEGIN();
    for (igraph_integer_t i=0; i < trunc(0.01 * igraph_ecount(&g)); i++) {
        /* For reproducibility, do not write two RNG_...() calls on the same line
         * as the C language does not guarantee any evaluation order between them. */
        igraph_real_t w = RNG_UNIF(-0.02, 0);
        VECTOR(weights)[RNG_INTEGER(0, igraph_ecount(&g)-1)] = w;
    }
    RNG_END();

    BENCH(" 7 vcount=" TOSTR(VCOUNT) ", p=" TOSTR(DENS) ", Bellman-Ford (negative), " TOSTR(REP) "x",
          REPEAT(igraph_distances_bellman_ford(&g, &res, igraph_vss_all(), igraph_vss_all(), &weights, IGRAPH_OUT), REP)
    );
    BENCH(" 8 vcount=" TOSTR(VCOUNT) ", p=" TOSTR(DENS) ", Johnson (negative), " TOSTR(REP) "x",
          REPEAT(igraph_distances_johnson(&g, &res, igraph_vss_all(), igraph_vss_all(), &weights), REP)
    );
    BENCH(" 9 vcount=" TOSTR(VCOUNT) ", p=" TOSTR(DENS) ", Floyd-Warshall (negative), " TOSTR(REP) "x",
          REPEAT(igraph_distances_floyd_warshall(&g, &res, igraph_vss_all(), igraph_vss_all(), &weights, IGRAPH_OUT, IGRAPH_FLOYD_WARSHALL_ORIGINAL), REP)
    );
    BENCH(" 10 vcount=" TOSTR(VCOUNT) ", p=" TOSTR(DENS) ", Floyd-Warshall-tree-speedup (negative), " TOSTR(REP) "x",
          REPEAT(igraph_distances_floyd_warshall(&g, &res, igraph_vss_all(), igraph_vss_all(), &weights, IGRAPH_OUT, IGRAPH_FLOYD_WARSHALL_TREE), REP)
    );

    igraph_destroy(&g);

#undef VCOUNT
#undef DENS
#undef REP

    igraph_vector_destroy(&weights);
    igraph_matrix_destroy(&res);

    return 0;
}
