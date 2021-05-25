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
    igraph_integer_t value;
    igraph_bool_t checks = 1;

    igraph_small(&g, 7, IGRAPH_DIRECTED,
                 0, 1, 0, 2, 1, 2, 1, 3, 2, 4, 3, 4, 3, 5, 4, 5, 1, 6, 6, 3, -1);

    igraph_adhesion(&g, &value, checks);

    IGRAPH_ASSERT(value == 0);

    igraph_destroy(&g);

    igraph_small(&g, 7, IGRAPH_UNDIRECTED,
                 0, 1, 0, 2, 1, 2, 1, 3, 2, 4, 3, 4, 3, 5, 4, 5, 1, 6, 6, 3, -1);

    igraph_adhesion(&g, &value, checks);

    IGRAPH_ASSERT(value == 2);

    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();
    return 0;
}
