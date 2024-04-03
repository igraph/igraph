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
#include <stdexcept>
#include <climits>

using namespace std;

namespace igraph {

namespace walktrap {

bool operator<(const Edge& E1, const Edge& E2) {
    return (E1.neighbor < E2.neighbor);
}


Vertex::Vertex() {
    degree = 0;
    edges = nullptr;
    total_weight = 0.;
}

Vertex::~Vertex() {
    delete[] edges;
}

Graph::Graph() {
    nb_vertices = 0;
    nb_edges = 0;
    vertices = nullptr;
    total_weight = 0.;
}

Graph::~Graph () {
    delete[] vertices;
}

class Edge_list {
public:
    int* V1;
    int* V2;
    double* W;

    int size;
    int size_max;

    void add(int v1, int v2, double w);
    Edge_list() {
        size = 0;
        size_max = 1024;
        V1 = new int[1024];
        V2 = new int[1024];
        W = new double[1024];
    }

    ~Edge_list() {
        delete[] V1;
        delete[] V2;
        delete[] W;
    }
};

void Edge_list::add(int v1, int v2, double w) {
    if (size == size_max) {
        int* tmp1 = new int[2 * size_max];
        int* tmp2 = new int[2 * size_max];
        double* tmp3 = new double[2 * size_max];
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

igraph_error_t Graph::convert_from_igraph(const igraph_t *graph,
                               const igraph_vector_t *weights) {
    Graph &G = *this;

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);

    // Avoid warnings with GCC when compiling with LTO.
    IGRAPH_ASSUME(no_of_nodes >= 0);
    IGRAPH_ASSUME(no_of_edges >= 0);

    // Refactoring the walktrap code to support larger graphs is pointless
    // as running the algorithm on them would take an impractically long time.
    if (no_of_nodes > INT_MAX || no_of_edges > INT_MAX) {
        IGRAPH_ERROR("Graph too large for walktrap community detection.", IGRAPH_EINVAL);
    }

    Edge_list EL;

    for (igraph_integer_t i = 0; i < no_of_edges; i++) {
        igraph_real_t w = weights ? VECTOR(*weights)[i] : 1.0;
        EL.add(IGRAPH_FROM(graph, i), IGRAPH_TO(graph, i), w);
    }

    G.nb_vertices = no_of_nodes;
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
        int deg = G.vertices[i].degree;
        double w = (deg == 0) ? 1.0 : (G.vertices[i].total_weight / double(deg));
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
        /* Check for zero strength, as it may lead to crashes the in walktrap algorithm.
         * See https://github.com/igraph/igraph/pull/2043 */
        if (G.vertices[i].total_weight == 0) {
            /* G.vertices will be destroyed by Graph::~Graph() */
            IGRAPH_ERROR("Vertex with zero strength found: all vertices must have positive strength for walktrap.",
                         IGRAPH_EINVAL);
        }
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

    return IGRAPH_SUCCESS;
}

}
}
