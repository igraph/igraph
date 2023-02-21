/*
   IGraph library.
   Copyright (C) 2020-2022  The igraph development team <igraph@igraph.org>

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
    igraph_t graph;
    igraph_vector_t hist;

    igraph_small(&graph, 6, 0,
                 1, 2, 2, 3, 3, 4, 4, 5, 5, 2, 2, 4,
                 -1);

    igraph_vector_init(&hist, 0);

    igraph_maximal_cliques_hist(&graph, &hist, 0, 0);
    igraph_vector_print(&hist);

    igraph_vector_destroy(&hist);
    igraph_destroy(&graph);

    VERIFY_FINALLY_STACK();

    return 0;
}
