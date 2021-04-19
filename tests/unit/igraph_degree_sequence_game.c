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

#include "test_utilities.inc"

int main() {
    igraph_t g;
    igraph_vector_t outdeg, indeg, degrees, empty;
    igraph_bool_t is_simple, is_connected;
    igraph_error_handler_t *ehandler;

    igraph_real_t outarr[] = {2, 3, 2, 3, 3, 3, 3, 1, 4, 4};
    igraph_real_t inarr[]  = {3, 6, 2, 0, 2, 2, 4, 3, 3, 3};

    long int n = sizeof(outarr) / sizeof(igraph_real_t);

    igraph_rng_seed(igraph_rng_default(), 333);

    igraph_vector_view(&outdeg, outarr, n);
    igraph_vector_view(&indeg,  inarr,  n);

    igraph_vector_init(&empty, 0);

    /* This vector is used to check that the degrees of the result
     * match the requested degrees. */
    igraph_vector_init(&degrees, 0);

    /* Configuration model, undirected non-simple graphs */

    igraph_degree_sequence_game(&g, &outdeg, NULL, IGRAPH_DEGSEQ_SIMPLE);

    IGRAPH_ASSERT(! igraph_is_directed(&g));
    IGRAPH_ASSERT(igraph_vcount(&g) == n);

    igraph_degree(&g, &degrees, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS);
    IGRAPH_ASSERT(igraph_vector_all_e(&outdeg, &degrees));

    igraph_destroy(&g);

    igraph_degree_sequence_game(&g, &empty, NULL, IGRAPH_DEGSEQ_SIMPLE);
    IGRAPH_ASSERT(igraph_vcount(&g) == 0);
    igraph_destroy(&g);


    /* Configuration model, directed non-simple graphs */

    igraph_degree_sequence_game(&g, &outdeg, &indeg, IGRAPH_DEGSEQ_SIMPLE);

    IGRAPH_ASSERT(igraph_is_directed(&g));
    IGRAPH_ASSERT(igraph_vcount(&g) == n);

    igraph_degree(&g, &degrees, igraph_vss_all(), IGRAPH_OUT, IGRAPH_LOOPS);
    IGRAPH_ASSERT(igraph_vector_all_e(&outdeg, &degrees));

    igraph_degree(&g, &degrees, igraph_vss_all(), IGRAPH_IN, IGRAPH_LOOPS);
    IGRAPH_ASSERT(igraph_vector_all_e(&indeg, &degrees));

    igraph_destroy(&g);

    igraph_degree_sequence_game(&g, &empty, &empty, IGRAPH_DEGSEQ_SIMPLE);
    IGRAPH_ASSERT(igraph_vcount(&g) == 0);
    igraph_destroy(&g);


    /* Configuration model, undirected simple graphs */

    igraph_degree_sequence_game(&g, &outdeg, NULL, IGRAPH_DEGSEQ_SIMPLE_NO_MULTIPLE_UNIFORM);

    IGRAPH_ASSERT(! igraph_is_directed(&g));
    IGRAPH_ASSERT(igraph_vcount(&g) == n);

    igraph_is_simple(&g, &is_simple);
    IGRAPH_ASSERT(is_simple);

    igraph_degree(&g, &degrees, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS);
    IGRAPH_ASSERT(igraph_vector_all_e(&outdeg, &degrees));

    igraph_destroy(&g);

    igraph_degree_sequence_game(&g, &empty, NULL, IGRAPH_DEGSEQ_SIMPLE_NO_MULTIPLE_UNIFORM);
    IGRAPH_ASSERT(igraph_vcount(&g) == 0);
    igraph_destroy(&g);


    /* Configuration model, directed simple graphs */

    igraph_degree_sequence_game(&g, &outdeg, &indeg, IGRAPH_DEGSEQ_SIMPLE_NO_MULTIPLE_UNIFORM);

    IGRAPH_ASSERT(igraph_is_directed(&g));
    IGRAPH_ASSERT(igraph_vcount(&g) == n);

    igraph_is_simple(&g, &is_simple);
    IGRAPH_ASSERT(is_simple);

    igraph_degree(&g, &degrees, igraph_vss_all(), IGRAPH_OUT, IGRAPH_LOOPS);
    IGRAPH_ASSERT(igraph_vector_all_e(&outdeg, &degrees));

    igraph_degree(&g, &degrees, igraph_vss_all(), IGRAPH_IN, IGRAPH_LOOPS);
    IGRAPH_ASSERT(igraph_vector_all_e(&indeg, &degrees));

    igraph_destroy(&g);

    igraph_degree_sequence_game(&g, &empty, &empty, IGRAPH_DEGSEQ_SIMPLE_NO_MULTIPLE_UNIFORM);
    IGRAPH_ASSERT(igraph_vcount(&g) == 0);
    igraph_destroy(&g);


    /* Fast heuristic method, undirected simple graphs */

    igraph_degree_sequence_game(&g, &outdeg, NULL, IGRAPH_DEGSEQ_SIMPLE_NO_MULTIPLE);

    IGRAPH_ASSERT(! igraph_is_directed(&g));
    IGRAPH_ASSERT(igraph_vcount(&g) == n);

    igraph_is_simple(&g, &is_simple);
    IGRAPH_ASSERT(is_simple);

    igraph_degree(&g, &degrees, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS);
    IGRAPH_ASSERT(igraph_vector_all_e(&outdeg, &degrees));

    igraph_destroy(&g);

    igraph_degree_sequence_game(&g, &empty, NULL, IGRAPH_DEGSEQ_SIMPLE_NO_MULTIPLE);
    IGRAPH_ASSERT(igraph_vcount(&g) == 0);
    igraph_destroy(&g);


    /* Fast heuristic method, directed simple graphs */

    igraph_degree_sequence_game(&g, &outdeg, &indeg, IGRAPH_DEGSEQ_SIMPLE_NO_MULTIPLE);

    IGRAPH_ASSERT(igraph_is_directed(&g));
    IGRAPH_ASSERT(igraph_vcount(&g) == n);

    igraph_is_simple(&g, &is_simple);
    IGRAPH_ASSERT(is_simple);

    igraph_degree(&g, &degrees, igraph_vss_all(), IGRAPH_OUT, IGRAPH_LOOPS);
    IGRAPH_ASSERT(igraph_vector_all_e(&outdeg, &degrees));

    igraph_degree(&g, &degrees, igraph_vss_all(), IGRAPH_IN, IGRAPH_LOOPS);
    IGRAPH_ASSERT(igraph_vector_all_e(&indeg, &degrees));

    igraph_destroy(&g);

    igraph_degree_sequence_game(&g, &empty, &empty, IGRAPH_DEGSEQ_SIMPLE_NO_MULTIPLE);
    IGRAPH_ASSERT(igraph_vcount(&g) == 0);
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
    IGRAPH_ASSERT(igraph_vector_all_e(&outdeg, &degrees));

    igraph_destroy(&g);

    igraph_degree_sequence_game(&g, &empty, NULL, IGRAPH_DEGSEQ_VL);
    IGRAPH_ASSERT(! igraph_is_directed(&g));
    IGRAPH_ASSERT(igraph_vcount(&g) == 0);
    igraph_destroy(&g);

    ehandler = igraph_set_error_handler(igraph_error_handler_ignore);
    /* This degree sequence contains a zero degree, so it cannot be realized by a connected graph. */
    IGRAPH_ASSERT(
            igraph_degree_sequence_game(&g, &indeg, NULL, IGRAPH_DEGSEQ_VL) == IGRAPH_EINVAL
            );
    igraph_set_error_handler(ehandler);

    igraph_vector_destroy(&degrees);
    igraph_vector_destroy(&empty);

    VERIFY_FINALLY_STACK();

    return 0;
}

