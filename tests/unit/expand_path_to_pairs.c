/*
   igraph library.
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

#include <igraph_paths.h>
#include "test_utilities.h"

igraph_error_t test_path_expansion(void) {
    igraph_vector_int_t path;

    igraph_vector_int_init(&path, 0);
    IGRAPH_CHECK(igraph_expand_path_to_pairs(&path));
    print_vector_int(&path);
    igraph_vector_int_destroy(&path);

    igraph_vector_int_init(&path, 0);
    igraph_vector_int_push_back(&path, 0);
    IGRAPH_CHECK(igraph_expand_path_to_pairs(&path));
    print_vector_int(&path);
    igraph_vector_int_destroy(&path);

    igraph_vector_int_init(&path, 0);
    igraph_vector_int_push_back(&path, 0);
    igraph_vector_int_push_back(&path, 1);
    IGRAPH_CHECK(igraph_expand_path_to_pairs(&path));
    print_vector_int(&path);
    igraph_vector_int_destroy(&path);

    igraph_vector_int_init(&path, 0);
    igraph_vector_int_push_back(&path, 2);
    igraph_vector_int_push_back(&path, 3);
    igraph_vector_int_push_back(&path, 5);
    igraph_vector_int_push_back(&path, 7);
    IGRAPH_CHECK(igraph_expand_path_to_pairs(&path));
    print_vector_int(&path);
    igraph_vector_int_destroy(&path);

    return IGRAPH_SUCCESS;
}

int main(void) {
    IGRAPH_ASSERT(test_path_expansion() == IGRAPH_SUCCESS);

    VERIFY_FINALLY_STACK();

    return 0;
}
