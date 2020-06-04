/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
/*
   IGraph library.
   Copyright (C) 2011-2012  Gabor Csardi <csardi.gabor@gmail.com>
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
#include <string.h>

int main() {
    igraph_t g;
    igraph_vector_t weights, result;
    igraph_bool_t dag;

    igraph_vector_init(&result, 0);

    /***********************************************************************/
    /* Approximation with Eades' method                                    */
    /***********************************************************************/

    /* Simple unweighted graph */
    igraph_small(&g, 0, IGRAPH_DIRECTED, 0, 1, 1, 2, 2, 0, 2, 3, 2, 4, 0, 4, 4, 3, 5, 0, 6, 5, -1);
    igraph_feedback_arc_set(&g, &result, 0, IGRAPH_FAS_APPROX_EADES);
    igraph_vector_print(&result);
    igraph_delete_edges(&g, igraph_ess_vector(&result));
    igraph_is_dag(&g, &dag);
    if (!dag) {
        return 1;
    }
    igraph_destroy(&g);

    /* Simple weighted graph */
    igraph_small(&g, 0, IGRAPH_DIRECTED, 0, 1, 1, 2, 2, 0, 2, 3, 2, 4, 0, 4, 4, 3, 5, 0, 6, 5, -1);
    igraph_vector_init_int_end(&weights, -1, 1, 1, 3, 1, 1, 1, 1, 1, 1, -1);
    igraph_feedback_arc_set(&g, &result, &weights, IGRAPH_FAS_APPROX_EADES);
    igraph_vector_print(&result);
    igraph_delete_edges(&g, igraph_ess_vector(&result));
    igraph_is_dag(&g, &dag);
    if (!dag) {
        return 2;
    }
    igraph_vector_destroy(&weights);
    igraph_destroy(&g);

    /* Simple unweighted graph with loops */
    igraph_small(&g, 0, IGRAPH_DIRECTED, 0, 1, 1, 2, 2, 0, 2, 3, 2, 4, 0, 4, 4, 3, 5, 0, 6, 5, 1, 1, 4, 4, -1);
    igraph_feedback_arc_set(&g, &result, 0, IGRAPH_FAS_APPROX_EADES);
    igraph_vector_print(&result);
    igraph_delete_edges(&g, igraph_ess_vector(&result));
    igraph_is_dag(&g, &dag);
    if (!dag) {
        return 3;
    }
    igraph_destroy(&g);

    igraph_vector_destroy(&result);

    return 0;
}
