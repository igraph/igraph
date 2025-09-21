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
#include <stdio.h>

int main(void) {

    igraph_t g;
    igraph_bool_t result;

    /* Initialize the library. */
    igraph_setup();

    igraph_small(&g, 7, 0, 0, 1, 1, 2, 2, 3, 3, 0, 2, 4, 4, 5, 2, 5, -1);
    igraph_is_biconnected(&g, &result);
    printf("Graph is%sbiconnected.\n", result ? " " : " not ");
    igraph_destroy(&g);

    return 0;
}
