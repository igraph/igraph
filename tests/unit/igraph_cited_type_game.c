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
    igraph_vector_t pref, types;

    igraph_rng_seed(igraph_rng_default(), 42);

    printf("No nodes:\n");
    igraph_vector_init_int(&pref, 2, 1, 1);
    igraph_vector_init_int(&types, 0);
    IGRAPH_ASSERT(igraph_cited_type_game(&g, /*nodes*/ 0, /*types*/ &types, /*pref*/ &pref, /*edges_per_step*/ 5, /*directed*/ 0) == IGRAPH_SUCCESS);
    print_graph_canon(&g);
    igraph_destroy(&g);
    igraph_vector_destroy(&pref);
    igraph_vector_destroy(&types);

    printf("No edges:\n");
    igraph_vector_init_int(&pref, 2, 1, 1);
    igraph_vector_init_int(&types, 3, 1, 1, 1);
    IGRAPH_ASSERT(igraph_cited_type_game(&g, /*nodes*/ 3, /*types*/ &types, /*pref*/ &pref, /*edges_per_step*/ 0, /*directed*/ 0) == IGRAPH_SUCCESS);
    print_graph_canon(&g);
    igraph_destroy(&g);
    igraph_vector_destroy(&pref);
    igraph_vector_destroy(&types);

    printf("Make a star of double edges:\n");
    igraph_vector_init_real(&pref, 3, 1.0, 0.0, 0.0);
    igraph_vector_init_int(&types, 5, 0, 1, 1, 1, 1);
    IGRAPH_ASSERT(igraph_cited_type_game(&g, /*nodes*/ 5, /*types*/ &types, /*pref*/ &pref, /*edges_per_step*/ 2, /*directed*/ 1) == IGRAPH_SUCCESS);
    print_graph_canon(&g);
    igraph_destroy(&g);
    igraph_vector_destroy(&pref);
    igraph_vector_destroy(&types);

    printf("Make a line:\n");
    igraph_vector_init_real(&pref, 7, 1.0e-30, 1.0e-20, 1.0e-10, 1.0, 1.0e+10, 1.0e+20, 0.0);
    igraph_vector_init_int(&types, 7, 0, 1, 2, 3, 4, 5, 6);
    IGRAPH_ASSERT(igraph_cited_type_game(&g, /*nodes*/ 7, /*types*/ &types, /*pref*/ &pref, /*edges_per_step*/ 1, /*directed*/ 1) == IGRAPH_SUCCESS);
    print_graph_canon(&g);
    igraph_destroy(&g);
    igraph_vector_destroy(&pref);
    igraph_vector_destroy(&types);

    VERIFY_FINALLY_STACK();
    igraph_set_error_handler(igraph_error_handler_ignore);

    printf("Checking negative number of nodes error handling.\n");
    igraph_vector_init_real(&pref, 2, 1.0, 1.0);
    igraph_vector_init_int(&types, 2, 0, 1);
    IGRAPH_ASSERT(igraph_cited_type_game(&g, /*nodes*/ -5, /*types*/ &types, /*pref*/ &pref, /*edges_per_step*/ 5, /*directed*/ 0) == IGRAPH_EINVAL);
    igraph_vector_destroy(&pref);
    igraph_vector_destroy(&types);

    printf("Checking too few types error handling.\n");
    igraph_vector_init_real(&pref, 1, 1.0);
    igraph_vector_init_int(&types, 0);
    IGRAPH_ASSERT(igraph_cited_type_game(&g, /*nodes*/ 1, /*types*/ &types, /*pref*/ &pref, /*edges_per_step*/ 5, /*directed*/ 0) == IGRAPH_EINVAL);
    igraph_vector_destroy(&pref);
    igraph_vector_destroy(&types);

    printf("Checking too many types error handling.\n");
    igraph_vector_init_real(&pref, 3, 1.0, 1.0, 1.0);
    igraph_vector_init_int(&types, 2, 0, 1);
    IGRAPH_ASSERT(igraph_cited_type_game(&g, /*nodes*/ 1, /*types*/ &types, /*pref*/ &pref, /*edges_per_step*/ 5, /*directed*/ 0) == IGRAPH_EINVAL);
    igraph_vector_destroy(&pref);
    igraph_vector_destroy(&types);

    printf("Checking negative type for error handling.\n");
    igraph_vector_init_real(&pref, 2, 1.0, 1.0);
    igraph_vector_init_int(&types, 2, 0, -5);
    IGRAPH_ASSERT(igraph_cited_type_game(&g, /*nodes*/ 2, /*types*/ &types, /*pref*/ &pref, /*edges_per_step*/ 5, /*directed*/ 0) == IGRAPH_EINVAL);
    igraph_vector_destroy(&pref);
    igraph_vector_destroy(&types);

    printf("Checking too big type for error handling.\n");
    igraph_vector_init_real(&pref, 2, 1.0, 1.0);
    igraph_vector_init_int(&types, 2, 0, 5);
    IGRAPH_ASSERT(igraph_cited_type_game(&g, /*nodes*/ 2, /*types*/ &types, /*pref*/ &pref, /*edges_per_step*/ 5, /*directed*/ 0) == IGRAPH_EINVAL);
    igraph_vector_destroy(&pref);
    igraph_vector_destroy(&types);

    printf("Checking negative preference error handling.\n");
    igraph_vector_init_real(&pref, 2, 1.0, -1.0);
    igraph_vector_init_int(&types, 2, 0, 1);
    IGRAPH_ASSERT(igraph_cited_type_game(&g, /*nodes*/ 2, /*types*/ &types, /*pref*/ &pref, /*edges_per_step*/ 5, /*directed*/ 0) == IGRAPH_EINVAL);
    igraph_vector_destroy(&pref);
    igraph_vector_destroy(&types);

    VERIFY_FINALLY_STACK();
    return 0;
}
