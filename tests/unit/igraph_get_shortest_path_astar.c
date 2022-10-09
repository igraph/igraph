/* IGraph library.
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

igraph_error_t heuristic(igraph_real_t *result, igraph_integer_t vertex_id, void *extra) {
    (void) vertex_id;
    (void) extra;
    *result = 0;
    return IGRAPH_SUCCESS;
}

int main(void) {

    igraph_t g;
    igraph_vector_int_t vertices, edges;
    igraph_vector_t weights;

    igraph_vector_int_init(&vertices, 0);
    igraph_vector_int_init(&edges, 0);

    igraph_erdos_renyi_game_gnp(&g, 1000, 0.5, IGRAPH_DIRECTED, IGRAPH_NO_LOOPS);
    igraph_vector_init(&weights, igraph_ecount(&g));

    for (int i = 0; i < igraph_ecount(&g); i ++) {
        VECTOR(weights)[i] = RNG_UNIF01();
    }

    igraph_get_shortest_path_dijkstra(&g, &vertices,
                                       &edges, /*from=*/ 0, /*to=*/ 100,
                                       /*weights=*/ &weights, /*mode=*/ IGRAPH_OUT);
    igraph_get_shortest_path_astar(&g, &vertices,
                                       &edges, /*from=*/ 0, /*to=*/ 100,
                                       /*weights=*/ &weights, /*mode=*/ IGRAPH_OUT, heuristic, NULL);

    igraph_vector_int_destroy(&vertices);
    igraph_vector_int_destroy(&edges);
    igraph_vector_destroy(&weights);

    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();
    return 0;
}
