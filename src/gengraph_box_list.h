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
// This class allows to maintain a list of vertices,
// sorted by degree (largest degrees first)
// Operations allowed :
// - get the vertex having max degree -> Cost = O(1)
// - remove any vertex from the graph -> Cost = Sum(degrees of neighbours)
//                                       [ could be O(degree) if optimized ]

#ifndef _BOX_LIST_H
#define _BOX_LIST_H

#ifndef _MSC_VER
    #ifndef register
        #define register
    #endif
#endif

namespace gengraph {

class box_list {

private:
    int n;     // INITIAL number of vertices
    int dmax;  // CURRENT Maximum degree
    int *deg;  // CURRENT Degrees (points directly to the deg[] of the graph

    // Vertices are grouped by degree: one double-chained lists for each degree
    int *list;        // list[d-1] is the head of list of vertices of degree d
    int *next;        // next[v]/prev[v] are the vertices next/previous to v
    int *prev;        //   in the list where v belongs
    void pop(int);    // pop(v) just removes v from its list
    void insert(int); // insert(v) insert v at the head of its list

public:

    // Ctor. Takes O(n) time.
    box_list(int n0, int *deg0);

    // Dtor
    ~box_list();

    // Self-explaining inline routines
    inline bool is_empty() {
        return dmax < 1;
    };
    inline int get_max()   {
        return list[dmax - 1];
    };
    inline int get_one()   {
        return list[0];
    };
    inline int get_min()   {
        int i = 0;
        while (list[i] < 0) {
            i++;
        }
        return list[i];
    };

    // Remove v from box_list
    // Also, semi-remove vertex v from graph: all neighbours of v will swap
    // their last neighbour wit hv, and then decrease their degree, so
    // that any arc w->v virtually disappear
    // Actually, adjacency lists are just permuted, and deg[] is changed
    void pop_vertex(int v, int **neigh);
};

} // namespace gengraph

#endif //_BOX_LIST_H
