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

int main(void) {
    igraph_t g;
    igraph_vector_int_t outdeg, indeg;
    igraph_vector_int_t vec;
    igraph_bool_t is_simple;

    /* Initialize the library. */
    igraph_setup();

    /* Set random seed for reproducibility */
    igraph_rng_seed(igraph_rng_default(), 42);

    igraph_vector_int_init_int(&outdeg, 10, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3);
    igraph_vector_int_init_int(&indeg, 10, 4, 4, 2, 2, 4, 4, 2, 2, 3, 3);
    igraph_vector_int_init(&vec, 0);

    /* checking the configuration model, undirected graphs */
    igraph_degree_sequence_game(&g, &outdeg, 0, IGRAPH_DEGSEQ_CONFIGURATION);
    if (igraph_is_directed(&g) || igraph_vcount(&g) != 10) {
        return 1;
    }
    if (igraph_degree(&g, &vec, igraph_vss_all(), IGRAPH_OUT, IGRAPH_LOOPS)) {
        return 2;
    }
    igraph_vector_int_print(&vec);
    igraph_destroy(&g);

    /* checking the Viger-Latapy method, undirected graphs */
    igraph_degree_sequence_game(&g, &outdeg, 0, IGRAPH_DEGSEQ_VL);
    if (igraph_is_directed(&g) || igraph_vcount(&g) != 10) {
        return 3;
    }
    if (igraph_is_simple(&g, &is_simple, IGRAPH_DIRECTED) || !is_simple) {
        return 4;
    }
    if (igraph_degree(&g, &vec, igraph_vss_all(), IGRAPH_OUT, IGRAPH_NO_LOOPS)) {
        return 5;
    }
    igraph_vector_int_print(&vec);
    igraph_destroy(&g);

    /* checking the configuration model, directed graphs */
    igraph_degree_sequence_game(&g, &outdeg, &indeg, IGRAPH_DEGSEQ_CONFIGURATION);
    if (!igraph_is_directed(&g) || igraph_vcount(&g) != 10) {
        return 6;
    }
    if (igraph_degree(&g, &vec, igraph_vss_all(), IGRAPH_OUT, IGRAPH_LOOPS)) {
        return 7;
    }
    igraph_vector_int_print(&vec);
    if (igraph_degree(&g, &vec, igraph_vss_all(), IGRAPH_IN, IGRAPH_LOOPS)) {
        return 8;
    }
    igraph_vector_int_print(&vec);
    igraph_destroy(&g);

    /* checking the fast heuristic method, undirected graphs */
    igraph_degree_sequence_game(&g, &outdeg, 0, IGRAPH_DEGSEQ_FAST_HEUR_SIMPLE);
    if (igraph_is_directed(&g) || igraph_vcount(&g) != 10) {
        return 9;
    }
    if (igraph_is_simple(&g, &is_simple, IGRAPH_DIRECTED) || !is_simple) {
        return 10;
    }
    if (igraph_degree(&g, &vec, igraph_vss_all(), IGRAPH_OUT, IGRAPH_LOOPS)) {
        return 11;
    }
    igraph_vector_int_print(&vec);
    igraph_destroy(&g);

    /* checking the fast heuristic method, directed graphs */
    igraph_degree_sequence_game(&g, &outdeg, &indeg, IGRAPH_DEGSEQ_FAST_HEUR_SIMPLE);
    if (!igraph_is_directed(&g) || igraph_vcount(&g) != 10) {
        return 12;
    }
    if (igraph_is_simple(&g, &is_simple, IGRAPH_DIRECTED) || !is_simple) {
        return 13;
    }
    if (igraph_degree(&g, &vec, igraph_vss_all(), IGRAPH_OUT, IGRAPH_LOOPS)) {
        return 14;
    }
    igraph_vector_int_print(&vec);
    if (igraph_degree(&g, &vec, igraph_vss_all(), IGRAPH_IN, IGRAPH_LOOPS)) {
        return 15;
    }
    igraph_vector_int_print(&vec);
    igraph_destroy(&g);

    igraph_vector_int_destroy(&vec);
    igraph_vector_int_destroy(&outdeg);
    igraph_vector_int_destroy(&indeg);

    return 0;
}
