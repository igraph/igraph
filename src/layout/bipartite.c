/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
/*
   IGraph library.
   Copyright (C) 2003-2020  The igraph development team

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
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include "igraph_layout.h"

#include "igraph_interface.h"

/**
 * \function igraph_layout_bipartite
 * Simple layout for bipartite graphs.
 *
 * The layout is created by first placing the vertices in two rows,
 * according to their types. Then the positions within the rows are
 * optimized to minimize edge crossings, by calling \ref
 * igraph_layout_sugiyama().
 *
 * \param graph The input graph.
 * \param types A boolean vector containing ones and zeros, the vertex
 *     types. Its length must match the number of vertices in the graph.
 * \param res Pointer to an initialized matrix, the result, the x and
 *     y coordinates are stored here.
 * \param hgap The preferred minimum horizontal gap between vertices
 *     in the same layer (i.e. vertices of the same type).
 * \param vgap  The distance between layers.
 * \param maxiter Maximum number of iterations in the crossing
 *     minimization stage. 100 is a reasonable default; if you feel
 *     that you have too many edge crossings, increase this.
 * \return Error code.
 *
 * \sa \ref igraph_layout_sugiyama().
 */
int igraph_layout_bipartite(const igraph_t *graph,
                            const igraph_vector_bool_t *types,
                            igraph_matrix_t *res, igraph_real_t hgap,
                            igraph_real_t vgap, long int maxiter) {

    long int i, no_of_nodes = igraph_vcount(graph);
    igraph_vector_t layers;

    if (igraph_vector_bool_size(types) != no_of_nodes) {
        IGRAPH_ERRORF("The vertex type vector size (%ld) should be equal to the number of nodes (%ld).",
                      IGRAPH_EINVAL, igraph_vector_bool_size(types), no_of_nodes);
    }
    if (hgap < 0) {
        IGRAPH_ERRORF("The horizontal gap cannot be negative, got %f.", IGRAPH_EINVAL, hgap);
    }

    IGRAPH_VECTOR_INIT_FINALLY(&layers, no_of_nodes);
    for (i = 0; i < no_of_nodes; i++) {
        VECTOR(layers)[i] = VECTOR(*types)[i] ? 0 : 1;
    }

    IGRAPH_CHECK(igraph_layout_sugiyama(graph, res, /*extd_graph=*/ 0,
                                        /*extd_to_orig_eids=*/ 0, &layers, hgap,
                                        vgap, maxiter, /*weights=*/ 0));

    igraph_vector_destroy(&layers);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}
