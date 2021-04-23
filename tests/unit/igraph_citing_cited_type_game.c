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
    igraph_matrix_t pref_empty, pref_bipartite, pref_line;
    igraph_vector_t types_empty, types_bipartite, types_line;
    igraph_bool_t bipartite;
    int bipartite_elem[] = {0, 1, 1, 0};
    int line_elem[] = {0, 0, 1, 0, 0,
                       1, 0, 0, 0, 0,
                       0, 1, 0, 0, 0,
                       0, 0, 1, 0, 0,
                       0, 0, 0, 1, 0,
    };

    igraph_rng_seed(igraph_rng_default(), 42);

    igraph_matrix_init(&pref_empty, 0, 0);
    igraph_vector_init(&types_empty, 0);

    matrix_init_int_row_major(&pref_bipartite, 2, 2, bipartite_elem);
    igraph_vector_init_int(&types_bipartite, 10, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0);

    matrix_init_int_row_major(&pref_line, 5, 5, line_elem);
    igraph_vector_init_int(&types_line, 5, 0, 1, 2, 3, 4);

    printf("No nodes.\n");
    IGRAPH_ASSERT(igraph_citing_cited_type_game(&g, /*nodes*/ 0, &types_empty, &pref_empty, /*edges_per_step*/ 5, /*directed*/ 0) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_vcount(&g) == 0);
    igraph_destroy(&g);

    printf("Bipartite graph.\n");
    IGRAPH_ASSERT(igraph_citing_cited_type_game(&g, /*nodes*/ 10, &types_bipartite, &pref_bipartite, /*edges_per_step*/ 5, /*directed*/ 0) == IGRAPH_SUCCESS);
    igraph_is_bipartite(&g, &bipartite, NULL);
    IGRAPH_ASSERT(bipartite);
    IGRAPH_ASSERT(igraph_vcount(&g) == 10);
    IGRAPH_ASSERT(igraph_ecount(&g) == 45);
    igraph_destroy(&g);

    printf("No edges.\n");
    IGRAPH_ASSERT(igraph_citing_cited_type_game(&g, /*nodes*/ 10, &types_bipartite, &pref_bipartite, /*edges_per_step*/ 0, /*directed*/ 0) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_vcount(&g) == 10);
    IGRAPH_ASSERT(igraph_ecount(&g) == 0);
    igraph_destroy(&g);

    printf("A line.\n");
    IGRAPH_ASSERT(igraph_citing_cited_type_game(&g, /*nodes*/ 5, &types_line, &pref_line, /*edges_per_step*/ 1, /*directed*/ 1) == IGRAPH_SUCCESS);
    print_graph_canon(&g);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();
    igraph_set_error_handler(igraph_error_handler_ignore);

    printf("Too few types for nodes.\n");
    IGRAPH_ASSERT(igraph_citing_cited_type_game(&g, /*nodes*/ 5, &types_empty, &pref_empty, /*edges_per_step*/ 1, /*directed*/ 1) == IGRAPH_EINVAL);

    printf("Too few prefs.\n");
    IGRAPH_ASSERT(igraph_citing_cited_type_game(&g, /*nodes*/ 5, &types_line, &pref_empty, /*edges_per_step*/ 1, /*directed*/ 1) == IGRAPH_EINVAL);

    igraph_matrix_destroy(&pref_empty);
    igraph_vector_destroy(&types_empty);
    igraph_matrix_destroy(&pref_bipartite);
    igraph_vector_destroy(&types_bipartite);
    igraph_matrix_destroy(&pref_line);
    igraph_vector_destroy(&types_line);
    VERIFY_FINALLY_STACK();
    return 0;
}
