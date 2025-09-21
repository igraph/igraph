/*
   igraph library.
   Copyright (C) 2013-2024  The igraph development team <igraph@igraph.org>

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

#include "bench.h"

int main(void) {

    igraph_t g;
    igraph_int_t toremovev[] = {
        2609,  2098, 14517,  7540, 19560,  8855,
        5939, 14947,   441, 16976, 19642,  4188,
        15447, 11837,  2333,  7309, 18539, 14099,
        14264,  9240
    };
    const igraph_vector_int_t toremove =
        igraph_vector_int_view(toremovev, sizeof(toremovev) / sizeof(toremovev[0]));;
    igraph_vector_int_list_t res;

    BENCH_INIT();

    igraph_full(&g, 200, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    igraph_delete_edges(&g, igraph_ess_vector(&toremove));

    igraph_vector_int_list_init(&res, 0);

    BENCH(" 1 Maximal cliques of almost complete graph",
          igraph_maximal_cliques(&g, &res,
                                 /* min_size= */ IGRAPH_UNLIMITED, /* max_size= */ IGRAPH_UNLIMITED,
                                 /* max_results= */ IGRAPH_UNLIMITED);
         );

    igraph_destroy(&g);

    igraph_vector_int_list_destroy(&res);

    return 0;
}
