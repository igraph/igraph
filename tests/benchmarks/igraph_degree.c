
#include <igraph.h>

#include "bench.h"

#define TOSTR1(x) #x
#define TOSTR(x) TOSTR1(x)

/* Calls igraph_degree_1() separately for each vertex. */
void deg1(const igraph_t *g, igraph_vector_int_t *res, igraph_neimode_t mode, igraph_bool_t loops) {
    igraph_integer_t n = igraph_vcount(g);
    igraph_vector_int_resize(res, n);
    for (igraph_integer_t i=0; i < n; i++) {
        igraph_degree_1(g, &VECTOR(*res)[i], i, mode, loops);
    }
}

/* Calls igraph_degree() separately for each vertex, pre-allocates work vector. */
void degv(const igraph_t *g, igraph_vector_int_t *res, igraph_neimode_t mode, igraph_bool_t loops) {
    igraph_integer_t n = igraph_vcount(g);
    igraph_vector_int_t work;

    igraph_vector_int_resize(res, n);
    igraph_vector_int_init(&work, 1);
    for (igraph_integer_t i=0; i < n; i++) {
        igraph_degree(g, &work, igraph_vss_1(i), mode, loops);
    }
    igraph_vector_int_destroy(&work);
}

/* Calls igraph_degree() separately for each vertex, allocates a new work vector each time. */
void degv2(const igraph_t *g, igraph_vector_int_t *res, igraph_neimode_t mode, igraph_bool_t loops) {
    igraph_integer_t n = igraph_vcount(g);

    for (igraph_integer_t i=0; i < n; i++) {
        igraph_vector_int_t work;
        igraph_vector_int_init(&work, 1);
        igraph_degree(g, &work, igraph_vss_1(i), mode, loops);
        igraph_vector_int_destroy(&work);
    }
}

int main(void) {
    igraph_t g;
    igraph_vector_int_t degs;

    igraph_rng_seed(igraph_rng_default(), 137);
    BENCH_INIT();

    igraph_barabasi_game(&g, 100000, 1, 10, NULL, true, 0, IGRAPH_UNDIRECTED, IGRAPH_BARABASI_PSUMTREE_MULTIPLE, NULL);
    igraph_vector_int_init(&degs, igraph_vcount(&g));

#define REP 1000

    printf("Count loops, O(1) per degree, fast methods only.\n");

    BENCH(" 1 igraph_degree(), preferential attachment n=100000, m=10, " TOSTR(REP) "x",
          REPEAT(igraph_degree(&g, &degs, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS), REP);
    );
    BENCH(" 2 deg1(), preferential attachment n=100000, m=10, " TOSTR(REP) "x",
          REPEAT(deg1(&g, &degs, IGRAPH_ALL, IGRAPH_LOOPS), REP);
    );

#undef REP

#define REP 100

    printf("\nCount loops, O(1) per degree.\n");

    BENCH(" 1 igraph_degree(), preferential attachment n=100000, m=10, " TOSTR(REP) "x",
          REPEAT(igraph_degree(&g, &degs, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS), REP);
    );
    BENCH(" 2 deg1(), preferential attachment n=100000, m=10, " TOSTR(REP) "x",
          REPEAT(deg1(&g, &degs, IGRAPH_ALL, IGRAPH_LOOPS), REP);
    );
    BENCH(" 3 degv(), preferential attachment n=100000, m=10, " TOSTR(REP) "x",
          REPEAT(degv(&g, &degs, IGRAPH_ALL, IGRAPH_LOOPS), REP);
    );
    BENCH(" 4 degv2(), preferential attachment n=100000, m=10, " TOSTR(REP) "x",
          REPEAT(degv2(&g, &degs, IGRAPH_ALL, IGRAPH_LOOPS), REP);
    );

#undef REP

#define REP 100

    printf("\nDo not count loops, O(d) per degree.\n");

    BENCH(" 1 igraph_degree(), preferential attachment n=100000, m=10, " TOSTR(REP) "x",
          REPEAT(igraph_degree(&g, &degs, igraph_vss_all(), IGRAPH_ALL, IGRAPH_NO_LOOPS), REP);
    );
    BENCH(" 2 deg1(), preferential attachment n=100000, m=10, " TOSTR(REP) "x",
          REPEAT(deg1(&g, &degs, IGRAPH_ALL, IGRAPH_NO_LOOPS), REP);
    );
    BENCH(" 3 degv(), preferential attachment n=100000, m=10, " TOSTR(REP) "x",
          REPEAT(degv(&g, &degs, IGRAPH_ALL, IGRAPH_NO_LOOPS), REP);
    );
    BENCH(" 4 degv2(), preferential attachment n=100000, m=10, " TOSTR(REP) "x",
          REPEAT(degv2(&g, &degs, IGRAPH_ALL, IGRAPH_NO_LOOPS), REP);
    );

    igraph_vector_int_destroy(&degs);
    igraph_destroy(&g);

    return 0;
}
