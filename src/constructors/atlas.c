/*
   igraph library.
   Copyright (C) 2006-2025  The igraph development team <igraph@igraph.org>

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

#include "igraph_constructors.h"

#include "constructors/atlas-edges.h"

/**
 * \function igraph_atlas
 * \brief Create a small graph from the \quote Graph Atlas \endquote.
 *
 * The graph atlas contains all simple undirected unlabeled graphs on between
 * 0 and 7 vertices. The number of the graph is given as a parameter.
 * The graphs are listed:
 *      \olist
 *      \oli in increasing order of number of vertices;
 *      \oli for a fixed number of vertices, in increasing order of the
 *           number of edges;
 *      \oli for fixed numbers of vertices and edges, in lexicographically
 *           increasing order of the degree sequence, for example
 *           111223 &lt; 112222;
 *      \oli for fixed degree sequence, in increasing number of
 *           automorphisms.
 *      \endolist
 *
 * </para><para>
 * The data was converted from the NetworkX software package,
 * see https://networkx.org/.
 *
 * </para><para>
 * See \emb An Atlas of Graphs \eme by Ronald C. Read and Robin J. Wilson,
 * Oxford University Press, 1998.
 *
 * \param graph Pointer to an uninitialized graph object.
 * \param number The number of the graph to generate. Must be between 0 and
 *    1252 (inclusive). Graphs on 0-7 vertices start at numbers 0, 1, 2, 4,
 *    8, 19, 53, and 209, respectively.
 * \return Error code.
 *
 * Added in version 0.2.
 *
 * </para><para>
 * Time complexity: O(|V|+|E|), the number of vertices plus the number of
 * edges.
 *
 * \example examples/simple/igraph_atlas.c
 */
igraph_error_t igraph_atlas(igraph_t *graph, igraph_int_t number) {

    const igraph_int_t atlas_size =
        sizeof(igraph_i_atlas_edges_pos) / sizeof(igraph_i_atlas_edges_pos[0]);

    if (number < 0 ||
        number >= atlas_size) {
        IGRAPH_ERRORF("No such graph in atlas. "
                      "The graph index must be less than %" IGRAPH_PRId ".",
                      IGRAPH_EINVAL,
                      atlas_size);
    }

    igraph_int_t pos = igraph_i_atlas_edges_pos[number];
    igraph_int_t n = igraph_i_atlas_edges[pos];
    igraph_int_t e = igraph_i_atlas_edges[pos + 1];
    const igraph_vector_int_t edges = igraph_vector_int_view(igraph_i_atlas_edges + pos + 2, e * 2);

    IGRAPH_CHECK(igraph_create(graph,
                               &edges,
                               n,
                               IGRAPH_UNDIRECTED));

    return IGRAPH_SUCCESS;
}
