/*
 *
 * gengraph - generation of random simple connected graphs with prescribed
 *            degree sequence
 *
 * Copyright (C) 2006  Fabien Viger
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef _VERTEX_COVER_H
#define _VERTEX_COVER_H

// vertex_cover() builds a list of vertices which covers every edge of the graph
// Input is a classical adjacency-list graph
// As an output, vertex_cover() modify the degrees in degs[], so that
// any vertex with a degree > 0 belongs to the vertex coverage.
// Moreover, vertex_cover() keeps links[] intact, permuting only the adjacency lists

#include "gengraph_box_list.h"

#ifndef register
    #define register
#endif

namespace gengraph {

void vertex_cover(int n, int *links, int *deg, int **neigh = NULL) {
    int i;
    // create and initialize neigh[]
    if (neigh == NULL) {
        neigh = new int*[n];
        neigh[0] = links;
        for (i = 1; i < n; i++) {
            neigh[i] = neigh[i - 1] + deg[i];
        }
    }
    // create box_list
    box_list bl(n, deg);
    do {
        int v;
        // remove vertices adjacent to vertices of degree 1
        while ((v = bl.get_one()) >= 0) {
            bl.pop_vertex(v, neigh);
        }
        // remove vertex of max degree and its highest-degree neighbour
        if (!bl.is_empty()) {
            v = bl.get_max();
            int *w = neigh[v];
            register int v2 = *(w++);
            register int dm = deg[v2];
            register int k = deg[v] - 1;
            while (k--) if (deg[*(w++)] > dm) {
                    v2 = *(w - 1);
                    dm = deg[v2];
                };
            bl.pop_vertex(v, neigh);
            bl.pop_vertex(v2, neigh);
        }
    } while (!bl.is_empty());
}

} // namespace gengraph

#endif //_VERTEX_COVER_H
