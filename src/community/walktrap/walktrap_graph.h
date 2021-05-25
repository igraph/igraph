/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2007-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard st, Cambridge, MA, 02138 USA

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

/* The original version of this file was written by Pascal Pons
   The original copyright notice follows here */

// File: graph.h
//-----------------------------------------------------------------------------
// Walktrap v0.2 -- Finds community structure of networks using random walks
// Copyright (C) 2004-2005 Pascal Pons
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301 USA
//-----------------------------------------------------------------------------
// Author   : Pascal Pons
// Email    : pascal.pons@gmail.com
// Web page : http://www-rp.lip6.fr/~latapy/PP/walktrap.html
// Location : Paris, France
// Time     : June 2005
//-----------------------------------------------------------------------------
// see readme.txt for more details

/* FSF address above was fixed by Tamas Nepusz */


#ifndef WALKTRAP_GRAPH_H
#define WALKTRAP_GRAPH_H

#include "igraph_community.h"

namespace igraph {

namespace walktrap {

class Edge {            // code an edge of a given vertex
public:
    int neighbor;         // the number of the neighbor vertex
    float weight;         // the weight of the edge
};
bool operator<(const Edge& E1, const Edge& E2);


class Vertex {
public:
    Edge* edges;          // the edges of the vertex
    int degree;           // number of neighbors
    float total_weight;       // the total weight of the vertex

    Vertex();         // creates empty vertex
    ~Vertex();            // destructor
};

class Graph {
public:
    int nb_vertices;      // number of vertices
    int nb_edges;         // number of edges
    float total_weight;       // total weight of the edges
    Vertex* vertices;     // array of the vertices

    long memory();            // the total memory used in Bytes
    Graph();          // create an empty graph
    ~Graph();         // destructor
    char** index;         // to keep the real name of the vertices

    int convert_from_igraph(const igraph_t * igraph,
                            const igraph_vector_t *weights);
};

}
}        /* end of namespaces */

#endif // WALKTRAP_GRAPH_H
