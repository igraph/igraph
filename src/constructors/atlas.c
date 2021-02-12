/* -*- mode: C -*-  */
/*
   IGraph R package.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge, MA 02139 USA

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

#include "igraph_constructors.h"

#include "constructors/atlas-edges.h"

/**
 * \function igraph_atlas
 * \brief Create a small graph from the \quote Graph Atlas \endquote.
 *
 * </para><para>
 * The number of the graph is given as a parameter.
 * The graphs are listed: \olist
 *      \oli in increasing order of number of nodes;
 *      \oli for a fixed number of nodes, in increasing order of the
 *           number of edges;
 *      \oli for fixed numbers of nodes and edges, in increasing
 *           order of the degree sequence, for example 111223 &lt; 112222;
 *      \oli for fixed degree sequence, in increasing number of
 *           automorphisms.
 *      \endolist
 *
 * </para><para>
 * The data was converted from the NetworkX software package,
 * see http://networkx.github.io .
 *
 * </para><para>
 * See \emb An Atlas of Graphs \eme by Ronald C. Read and Robin J. Wilson,
 * Oxford University Press, 1998.
 *
 * \param graph Pointer to an uninitialized graph object.
 * \param number The number of the graph to generate.
 *
 * Added in version 0.2.</para><para>
 *
 * Time complexity: O(|V|+|E|), the number of vertices plus the number of
 * edges.
 *
 * \example examples/simple/igraph_atlas.c
 */
int igraph_atlas(igraph_t *graph, int number) {

    igraph_integer_t pos, n, e;
    igraph_vector_t v = IGRAPH_VECTOR_NULL;

    if (number < 0 ||
        number >= (int) (sizeof(igraph_i_atlas_edges_pos) / sizeof(long int))) {
        IGRAPH_ERROR("No such graph in atlas", IGRAPH_EINVAL);
    }

    pos = (igraph_integer_t) igraph_i_atlas_edges_pos[number];
    n = (igraph_integer_t) igraph_i_atlas_edges[pos];
    e = (igraph_integer_t) igraph_i_atlas_edges[pos + 1];

    IGRAPH_CHECK(igraph_create(graph,
                               igraph_vector_view(&v, igraph_i_atlas_edges + pos + 2,
                                       e * 2),
                               n, IGRAPH_UNDIRECTED));

    return 0;
}
