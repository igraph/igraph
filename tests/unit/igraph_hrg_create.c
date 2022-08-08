/*
   IGraph library.
   Copyright (C) 2022  The igraph development team <igraph@igraph.org>

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

int main() {
    igraph_hrg_t hrg;
    igraph_t graph;
    igraph_vector_t prob;

    igraph_vector_init_real(&prob, 3, 0.1, 0.9);
    igraph_hrg_init(&hrg, 0);

    igraph_small(&graph, 3, IGRAPH_DIRECTED, 0,1, 0,2, -1);
    igraph_hrg_create(&hrg, &graph, &prob);

    igraph_destroy(&graph);
    igraph_vector_destroy(&prob);
    igraph_hrg_destroy(&hrg);

    VERIFY_FINALLY_STACK();

    return 0;
}
