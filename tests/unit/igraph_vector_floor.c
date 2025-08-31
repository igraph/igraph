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

#include <igraph.h>
#include "test_utilities.h"


int main(void) {
    igraph_vector_t from;
    igraph_vector_int_t to;

    igraph_vector_init_real(&from, 9, -0.6, -0.5, -0.4, -0.0, 0.0, 0.4, 0.5, 0.6, 1.1);
    igraph_vector_int_init(&to, 0);

    printf("From:\n");
    igraph_vector_print(&from);
    IGRAPH_ASSERT(igraph_vector_floor(&from, &to) == IGRAPH_SUCCESS);

    printf("To:\n");
    igraph_vector_int_print(&to);

    igraph_vector_int_destroy(&to);
    igraph_vector_destroy(&from);
    VERIFY_FINALLY_STACK();
    return 0;
}
