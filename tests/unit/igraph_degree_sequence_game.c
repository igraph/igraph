/*
   IGraph library.
   Copyright (C) 2021  The igraph development team <igraph@igraph.org>

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

igraph_bool_t compare_degrees(const igraph_vector_int_t* expected, const igraph_vector_int_t *observed) {
    igraph_integer_t i, n = igraph_vector_int_size(expected);

    if (igraph_vector_int_size(observed) != n) {
        return 0;
    }

    for (i = 0; i < n; i++) {
        if (VECTOR(*expected)[i] != VECTOR(*observed)[i]) {
            return 0;
        }
    }

    return 1;
}

int main(void) {
    igraph_t g;
    igraph_vector_int_t outdeg, indeg, empty;
    igraph_vector_int_t degrees;
    igraph_bool_t is_simple, is_connected;

    igraph_integer_t outarr[] = {2, 3, 2, 3, 3, 3, 3, 1, 4, 4};
    igraph_integer_t inarr[]  = {3, 6, 2, 0, 2, 2, 4, 3, 3, 3};

    igraph_integer_t n = sizeof(outarr) / sizeof(outarr[0]);

    igraph_rng_seed(igraph_rng_default(), 333);

    igraph_vector_int_view(&outdeg, outarr, n);
    igraph_vector_int_view(&indeg,  inarr,  n);

    igraph_vector_int_init(&empty, 0);

    /* This vector is used to check that the degrees of the result
     * match the requested degrees. */
    igraph_vector_int_init(&degrees, 0);

    /* Configuration model, undirected non-simple graphs */

    igraph_degree_sequence_game(&g, &outdeg, NULL, IGRAPH_DEGSEQ_CONFIGURATION);

    IGRAPH_ASSERT(! igraph_is_directed(&g));
    IGRAPH_ASSERT(igraph_vcount(&g) == n);

    igraph_degree(&g, &degrees, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS);
    IGRAPH_ASSERT(compare_degrees(&outdeg, &degrees));

    igraph_destroy(&g);

    igraph_degree_sequence_game(&g, &empty, NULL, IGRAPH_DEGSEQ_CONFIGURATION);
    IGRAPH_ASSERT(igraph_vcount(&g) == 0);
    igraph_destroy(&g);


    /* Configuration model, directed non-simple graphs */

    igraph_degree_sequence_game(&g, &outdeg, &indeg, IGRAPH_DEGSEQ_CONFIGURATION);

    IGRAPH_ASSERT(igraph_is_directed(&g));
    IGRAPH_ASSERT(igraph_vcount(&g) == n);

    igraph_degree(&g, &degrees, igraph_vss_all(), IGRAPH_OUT, IGRAPH_LOOPS);
    IGRAPH_ASSERT(compare_degrees(&outdeg, &degrees));

    igraph_degree(&g, &degrees, igraph_vss_all(), IGRAPH_IN, IGRAPH_LOOPS);
    IGRAPH_ASSERT(compare_degrees(&indeg, &degrees));

    igraph_destroy(&g);

    igraph_degree_sequence_game(&g, &empty, &empty, IGRAPH_DEGSEQ_CONFIGURATION);
    IGRAPH_ASSERT(igraph_vcount(&g) == 0);
    igraph_destroy(&g);


    /* Configuration model, undirected simple graphs */

    igraph_degree_sequence_game(&g, &outdeg, NULL, IGRAPH_DEGSEQ_CONFIGURATION_SIMPLE);

    IGRAPH_ASSERT(! igraph_is_directed(&g));
    IGRAPH_ASSERT(igraph_vcount(&g) == n);

    igraph_is_simple(&g, &is_simple);
    IGRAPH_ASSERT(is_simple);

    igraph_degree(&g, &degrees, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS);
    IGRAPH_ASSERT(compare_degrees(&outdeg, &degrees));

    igraph_destroy(&g);

    igraph_degree_sequence_game(&g, &empty, NULL, IGRAPH_DEGSEQ_CONFIGURATION_SIMPLE);
    IGRAPH_ASSERT(igraph_vcount(&g) == 0);
    igraph_destroy(&g);


    /* Configuration model, directed simple graphs */

    igraph_degree_sequence_game(&g, &outdeg, &indeg, IGRAPH_DEGSEQ_CONFIGURATION_SIMPLE);

    IGRAPH_ASSERT(igraph_is_directed(&g));
    IGRAPH_ASSERT(igraph_vcount(&g) == n);

    igraph_is_simple(&g, &is_simple);
    IGRAPH_ASSERT(is_simple);

    igraph_degree(&g, &degrees, igraph_vss_all(), IGRAPH_OUT, IGRAPH_LOOPS);
    IGRAPH_ASSERT(compare_degrees(&outdeg, &degrees));

    igraph_degree(&g, &degrees, igraph_vss_all(), IGRAPH_IN, IGRAPH_LOOPS);
    IGRAPH_ASSERT(compare_degrees(&indeg, &degrees));

    igraph_destroy(&g);

    igraph_degree_sequence_game(&g, &empty, &empty, IGRAPH_DEGSEQ_CONFIGURATION_SIMPLE);
    IGRAPH_ASSERT(igraph_vcount(&g) == 0);
    igraph_destroy(&g);


    /* Fast heuristic method, undirected simple graphs */

    igraph_degree_sequence_game(&g, &outdeg, NULL, IGRAPH_DEGSEQ_FAST_HEUR_SIMPLE);

    IGRAPH_ASSERT(! igraph_is_directed(&g));
    IGRAPH_ASSERT(igraph_vcount(&g) == n);

    igraph_is_simple(&g, &is_simple);
    IGRAPH_ASSERT(is_simple);

    igraph_degree(&g, &degrees, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS);
    IGRAPH_ASSERT(compare_degrees(&outdeg, &degrees));

    igraph_destroy(&g);

    igraph_degree_sequence_game(&g, &empty, NULL, IGRAPH_DEGSEQ_FAST_HEUR_SIMPLE);
    IGRAPH_ASSERT(igraph_vcount(&g) == 0);
    igraph_destroy(&g);


    /* Fast heuristic method, directed simple graphs */

    igraph_degree_sequence_game(&g, &outdeg, &indeg, IGRAPH_DEGSEQ_FAST_HEUR_SIMPLE);

    IGRAPH_ASSERT(igraph_is_directed(&g));
    IGRAPH_ASSERT(igraph_vcount(&g) == n);

    igraph_is_simple(&g, &is_simple);
    IGRAPH_ASSERT(is_simple);

    igraph_degree(&g, &degrees, igraph_vss_all(), IGRAPH_OUT, IGRAPH_LOOPS);
    IGRAPH_ASSERT(compare_degrees(&outdeg, &degrees));

    igraph_degree(&g, &degrees, igraph_vss_all(), IGRAPH_IN, IGRAPH_LOOPS);
    IGRAPH_ASSERT(compare_degrees(&indeg, &degrees));

    igraph_destroy(&g);

    igraph_degree_sequence_game(&g, &empty, &empty, IGRAPH_DEGSEQ_FAST_HEUR_SIMPLE);
    IGRAPH_ASSERT(igraph_vcount(&g) == 0);
    igraph_destroy(&g);


    /* Edge-swithching MCMC, undirected, simple graphs */

    igraph_degree_sequence_game(&g, &outdeg, NULL, IGRAPH_DEGSEQ_EDGE_SWITCHING_SIMPLE);

    igraph_is_simple(&g, &is_simple);
    IGRAPH_ASSERT(is_simple);

    igraph_degree(&g, &degrees, igraph_vss_all(), IGRAPH_OUT, IGRAPH_LOOPS);
    IGRAPH_ASSERT(compare_degrees(&outdeg, &degrees));

    igraph_destroy(&g);

    /* Edge-swithching MCMC, directed, simple graphs */

    igraph_degree_sequence_game(&g, &outdeg, &indeg, IGRAPH_DEGSEQ_EDGE_SWITCHING_SIMPLE);

    igraph_is_simple(&g, &is_simple);
    IGRAPH_ASSERT(is_simple);

    igraph_degree(&g, &degrees, igraph_vss_all(), IGRAPH_OUT, IGRAPH_LOOPS);
    IGRAPH_ASSERT(compare_degrees(&outdeg, &degrees));

    igraph_degree(&g, &degrees, igraph_vss_all(), IGRAPH_IN, IGRAPH_LOOPS);
    IGRAPH_ASSERT(compare_degrees(&indeg, &degrees));

    igraph_destroy(&g);


    /* Viger-Latapy method, undirected connected simple graphs */

    igraph_degree_sequence_game(&g, &outdeg, NULL, IGRAPH_DEGSEQ_VL);

    IGRAPH_ASSERT(! igraph_is_directed(&g));
    IGRAPH_ASSERT(igraph_vcount(&g) == n);

    igraph_is_simple(&g, &is_simple);
    IGRAPH_ASSERT(is_simple);

    igraph_is_connected(&g, &is_connected, IGRAPH_WEAK);
    IGRAPH_ASSERT(is_connected);

    igraph_degree(&g, &degrees, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS);
    IGRAPH_ASSERT(compare_degrees(&outdeg, &degrees));

    igraph_destroy(&g);

    igraph_degree_sequence_game(&g, &empty, NULL, IGRAPH_DEGSEQ_VL);
    IGRAPH_ASSERT(! igraph_is_directed(&g));
    IGRAPH_ASSERT(igraph_vcount(&g) == 0);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();

    /* This degree sequence contains a zero degree, so it cannot be realized by a connected graph. */
    CHECK_ERROR(igraph_degree_sequence_game(&g, &indeg, NULL, IGRAPH_DEGSEQ_VL), IGRAPH_EINVAL);

    igraph_vector_int_destroy(&degrees);
    igraph_vector_int_destroy(&empty);

    VERIFY_FINALLY_STACK();

    return 0;
}
