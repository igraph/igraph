/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard st, Cambridge MA, 02139 USA

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include <igraph.h>

#include "test_utilities.inc"

int main() {

    igraph_t g;
    igraph_vector_t v, res;
    igraph_bool_t is_dag;
    int ret;

    /* Test graph taken from http://en.wikipedia.org/wiki/Topological_sorting
     * @ 05.03.2006 */
    igraph_small(&g, 8, IGRAPH_DIRECTED,
                 0, 3, 0, 4, 1, 3, 2, 4, 2, 7, 3, 5, 3, 6, 3, 7, 4, 6,
                 -1);

    igraph_vector_init(&res, 0);

    igraph_is_dag(&g, &is_dag);
    if (!is_dag) {
        return 2;
    }

    igraph_topological_sorting(&g, &res, IGRAPH_OUT);
    print_vector_round(&res);
    igraph_topological_sorting(&g, &res, IGRAPH_IN);
    print_vector_round(&res);

    /* Error handling */
    igraph_set_error_handler(igraph_error_handler_ignore);

    /* Add a cycle: 5 -> 0 */
    igraph_vector_init_int(&v, 2, 5, 0);
    igraph_add_edges(&g, &v, 0);
    igraph_is_dag(&g, &is_dag);
    if (is_dag) {
        return 3;
    }
    ret = igraph_topological_sorting(&g, &res, IGRAPH_OUT);
    if (ret != IGRAPH_EINVAL) {
        return 1;
    }

    igraph_vector_destroy(&v);
    igraph_destroy(&g);

    /* This graph is the same but undirected */
    igraph_small(&g, 8, IGRAPH_UNDIRECTED,
                 0, 3, 0, 4, 1, 3, 2, 4, 2, 7, 3, 5, 3, 6, 3, 7, 4, 6,
                 -1);

    igraph_is_dag(&g, &is_dag);
    if (is_dag) {
        return 4;
    }

    ret = igraph_topological_sorting(&g, &res, IGRAPH_ALL);
    if (ret != IGRAPH_EINVAL) {
        return 1;
    }

    ret = igraph_topological_sorting(&g, &res, IGRAPH_OUT);
    if (ret != IGRAPH_EINVAL) {
        return 1;
    }

    igraph_destroy(&g);

    igraph_vector_destroy(&res);

    VERIFY_FINALLY_STACK();

    return 0;
}
