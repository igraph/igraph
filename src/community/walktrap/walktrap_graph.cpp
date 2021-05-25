/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2007-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

/* The original version of this file was written by Pascal Pons
   The original copyright notice follows here. The FSF address was
   fixed by Tamas Nepusz */

// File: graph.cpp
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

#include "walktrap_graph.h"
#include "igraph_interface.h"
#include <algorithm>
#include <cstring>      // strlen

using namespace std;

namespace igraph {

namespace walktrap {

bool operator<(const Edge& E1, const Edge& E2) {
    return (E1.neighbor < E2.neighbor);
}


Vertex::Vertex() {
    degree = 0;
    edges = 0;
    total_weight = 0.;
}

Vertex::~Vertex() {
    if (edges) {
        delete[] edges;
    }
}

Graph::Graph() {
    nb_vertices = 0;
    nb_edges = 0;
    vertices = 0;
    index = 0;
    total_weight = 0.;
}

Graph::~Graph () {
    if (vertices) {
        delete[] vertices;
    }
}

class Edge_list {
public:
    int* V1;
    int* V2;
    float* W;

    int size;
    int size_max;

    void add(int v1, int v2, float w);
    Edge_list() {
        size = 0;
        size_max = 1024;
        V1 = new int[1024];
        V2 = new int[1024];
        W = new float[1024];
    }
    ~Edge_list() {
        if (V1) {
            delete[] V1;
        }
        if (V2) {
            delete[] V2;
        }
        if (W) {
            delete[] W;
        }
    }
};

void Edge_list::add(int v1, int v2, float w) {
    if (size == size_max) {
        int* tmp1 = new int[2 * size_max];
        int* tmp2 = new int[2 * size_max];
        float* tmp3 = new float[2 * size_max];
        for (int i = 0; i < size_max; i++) {
            tmp1[i] = V1[i];
            tmp2[i] = V2[i];
            tmp3[i] = W[i];
        }
        delete[] V1;
        delete[] V2;
        delete[] W;
        V1 = tmp1;
        V2 = tmp2;
        W = tmp3;
        size_max *= 2;
    }
    V1[size] = v1;
    V2[size] = v2;
    W[size] = w;
    size++;
}

int Graph::convert_from_igraph(const igraph_t *graph,
                               const igraph_vector_t *weights) {
    Graph &G = *this;

    int max_vertex = (int)igraph_vcount(graph) - 1;
    long int no_of_edges = (long int)igraph_ecount(graph);
    long int i;
    long int deg;
    double w;

    Edge_list EL;

    for (i = 0; i < no_of_edges; i++) {
        igraph_integer_t from, to;
        int v1, v2;
        w = weights ? VECTOR(*weights)[i] : 1.0;
        igraph_edge(graph, i, &from, &to);
        v1 = (int)from; v2 = (int)to;
        EL.add(v1, v2, w);
    }

    G.nb_vertices = max_vertex + 1;
    G.vertices = new Vertex[G.nb_vertices];
    G.nb_edges = 0;
    G.total_weight = 0.0;

    for (int i = 0; i < EL.size; i++) {
        G.vertices[EL.V1[i]].degree++;
        G.vertices[EL.V2[i]].degree++;
        G.vertices[EL.V1[i]].total_weight += EL.W[i];
        G.vertices[EL.V2[i]].total_weight += EL.W[i];
        G.nb_edges++;
        G.total_weight += EL.W[i];
    }

    for (int i = 0; i < G.nb_vertices; i++) {
        deg = G.vertices[i].degree;
        w = (deg == 0) ? 1.0 : (G.vertices[i].total_weight / double(deg));
        G.vertices[i].edges = new Edge[deg + 1];
        G.vertices[i].edges[0].neighbor = i;
        G.vertices[i].edges[0].weight = w;
        G.vertices[i].total_weight += w;
        G.vertices[i].degree = 1;
    }

    for (int i = 0; i < EL.size; i++) {
        G.vertices[EL.V1[i]].edges[G.vertices[EL.V1[i]].degree].neighbor = EL.V2[i];
        G.vertices[EL.V1[i]].edges[G.vertices[EL.V1[i]].degree].weight = EL.W[i];
        G.vertices[EL.V1[i]].degree++;
        G.vertices[EL.V2[i]].edges[G.vertices[EL.V2[i]].degree].neighbor = EL.V1[i];
        G.vertices[EL.V2[i]].edges[G.vertices[EL.V2[i]].degree].weight = EL.W[i];
        G.vertices[EL.V2[i]].degree++;
    }

    for (int i = 0; i < G.nb_vertices; i++) {
        sort(G.vertices[i].edges, G.vertices[i].edges + G.vertices[i].degree);
    }

    for (int i = 0; i < G.nb_vertices; i++) { // merge multi edges
        int a = 0;
        for (int b = 1; b < G.vertices[i].degree; b++) {
            if (G.vertices[i].edges[b].neighbor == G.vertices[i].edges[a].neighbor) {
                G.vertices[i].edges[a].weight += G.vertices[i].edges[b].weight;
            } else {
                G.vertices[i].edges[++a] = G.vertices[i].edges[b];
            }
        }
        G.vertices[i].degree = a + 1;
    }

    return 0;
}

long Graph::memory() {
    size_t m = 0;
    m += size_t(nb_vertices) * sizeof(Vertex);
    m += 2 * size_t(nb_edges) * sizeof(Edge);
    m += sizeof(Graph);
    if (index != 0) {
        m += size_t(nb_vertices) * sizeof(char*);
        for (int i = 0; i < nb_vertices; i++) {
            m += strlen(index[i]) + 1;
        }
    }
    return m;
}

}
}
