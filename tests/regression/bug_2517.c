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

#include "../unit/test_utilities.h"

int main(void) {
    igraph_t g;
    igraph_vector_int_list_t result;

    igraph_vector_int_list_init(&result, 0);

    igraph_small(&g, 2, IGRAPH_UNDIRECTED, 0, 1, -1);
    igraph_all_minimal_st_separators(&g, &result);
	print_vector_int_list(&result);
    igraph_destroy(&g);

    igraph_small(&g, 3, IGRAPH_UNDIRECTED, 0, 1, -1);
    igraph_all_minimal_st_separators(&g, &result);
	print_vector_int_list(&result);
    igraph_destroy(&g);

    igraph_vector_int_list_destroy(&result);

    VERIFY_FINALLY_STACK();

    return 0;
}
