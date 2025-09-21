/* igraph library.
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

void run_test(igraph_int_t n, igraph_bool_t directed, const igraph_vector_int_t *outseq, const igraph_t *start) {
    igraph_barabasi_algorithm_t algos[] = { IGRAPH_BARABASI_BAG, IGRAPH_BARABASI_PSUMTREE, IGRAPH_BARABASI_PSUMTREE_MULTIPLE };
    igraph_t g;
    igraph_bool_t simple;

    for (size_t i=0; i < sizeof(algos) / sizeof(algos[0]); i++) {
        igraph_barabasi_algorithm_t algo = algos[i];
        igraph_barabasi_game(&g, n, 1, 2, outseq, true, 1, directed, algo, start);
        IGRAPH_ASSERT(igraph_vcount(&g) == n);
        IGRAPH_ASSERT(igraph_is_directed(&g) == directed);
        if (algo == IGRAPH_BARABASI_PSUMTREE) {
            igraph_is_simple(&g, &simple, IGRAPH_DIRECTED);
            IGRAPH_ASSERT(simple);
        }
        igraph_destroy(&g);
    }
}

int main(void) {
    igraph_vector_int_t v;
    igraph_t g;

    run_test(10, true, NULL, NULL);
    run_test(10, false, NULL, NULL);

    igraph_ring(&g, 5, IGRAPH_DIRECTED, false, true);
    run_test(13, true, NULL, &g);
    igraph_to_undirected(&g, IGRAPH_TO_UNDIRECTED_EACH, NULL);
    run_test(13, false, NULL, &g);
    igraph_destroy(&g);

    igraph_vector_int_init_int(&v, 4,
                               1, 2, 1, 3);
    run_test(4, true, &v, NULL);
    run_test(4, false, &v, NULL);
    igraph_vector_int_destroy(&v);

    CHECK_ERROR(igraph_barabasi_game(&g, -10, /*power=*/ 1, 1, NULL, false, /*A=*/ 1, IGRAPH_UNDIRECTED,
                               IGRAPH_BARABASI_BAG, /*start_from=*/ NULL), IGRAPH_EINVAL);
    CHECK_ERROR(igraph_barabasi_game(&g, 10, /*power=*/ 1, -2, NULL, false, /*A=*/ 1, IGRAPH_UNDIRECTED,
                               IGRAPH_BARABASI_BAG, /*start_from=*/ NULL), IGRAPH_EINVAL);
    igraph_vector_int_init(&v, 9);
    CHECK_ERROR(igraph_barabasi_game(&g, 10, /*power=*/ 1, 0, &v, false, /*A=*/ 1, IGRAPH_UNDIRECTED,
                               IGRAPH_BARABASI_BAG, /*start_from=*/ NULL), IGRAPH_EINVAL);
    igraph_vector_int_destroy(&v);

    VERIFY_FINALLY_STACK();

    return 0;
}
