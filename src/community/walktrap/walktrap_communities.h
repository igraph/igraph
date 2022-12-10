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
   The original copyright notice follows here. The FSF address was
   fixed by Tamas Nepusz */

// File: communities.h
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


#ifndef WALKTRAP_COMMUNITIES_H
#define WALKTRAP_COMMUNITIES_H

#include "walktrap_graph.h"
#include "walktrap_heap.h"
#include "config.h"

namespace igraph {

namespace walktrap {

class Communities;
class Probabilities {
public:
    static IGRAPH_THREAD_LOCAL double* tmp_vector1;    //
    static IGRAPH_THREAD_LOCAL double* tmp_vector2;    //
    static IGRAPH_THREAD_LOCAL int* id;       //
    static IGRAPH_THREAD_LOCAL int* vertices1;    //
    static IGRAPH_THREAD_LOCAL int* vertices2;    //
    static IGRAPH_THREAD_LOCAL int current_id;    //

    static IGRAPH_THREAD_LOCAL Communities* C;                    // pointer to all the communities
    static IGRAPH_THREAD_LOCAL int length;                        // length of the random walks


    int size;                         // number of probabilities stored
    int* vertices;                        // the vertices corresponding to the stored probabilities, 0 if all the probabilities are stored
    double* P;                         // the probabilities

    double compute_distance(const Probabilities* P2) const;   // compute the squared distance r^2 between this probability vector and P2
    explicit Probabilities(int community);                 // compute the probability vector of a community
    Probabilities(int community1, int community2);        // merge the probability vectors of two communities in a new one
    // the two communities must have their probability vectors stored

    ~Probabilities();                     // destructor
};

class Community {
public:

    Neighbor* first_neighbor; // first item of the list of adjacent communities
    Neighbor* last_neighbor;  // last item of the list of adjacent communities

    int this_community;       // number of this community
    int first_member;     // number of the first vertex of the community
    int last_member;      // number of the last vertex of the community
    int size;         // number of members of the community

    Probabilities* P;     // the probability vector, 0 if not stored.


    double sigma;          // sigma(C) of the community
    double internal_weight;    // sum of the weight of the internal edges
    double total_weight;       // sum of the weight of all the edges of the community (an edge between two communities is a half-edge for each community)

    int sub_communities[2];   // the two sub communities, -1 if no sub communities;
    int sub_community_of;     // number of the community in which this community has been merged
    // 0 if the community is active
    // -1 if the community is not used

    void add_neighbor(Neighbor* N);
    void remove_neighbor(Neighbor* N);

    Community();          // create an empty community
    ~Community();         // destructor
};

class Communities {
private:
    igraph_matrix_int_t *merges;
    igraph_integer_t mergeidx;
    igraph_vector_t *modularity;

public:
    Graph* G;         // the graph
    int* members;         // the members of each community represented as a chained list.
    // a community points to the first_member the array which contains
    // the next member (-1 = end of the community)
    Neighbor_heap* H;     // the distances between adjacent communities.


    Community* communities;   // array of the communities

    int nb_communities;       // number of valid communities
    int nb_active_communities;    // number of active communities

    Communities(Graph* G, int random_walks_length = 3,
                igraph_matrix_int_t *merges = nullptr,
                igraph_vector_t *modularity = nullptr);  // Constructor
    ~Communities();                   // Destructor

    void merge_communities(Neighbor* N);          // create a community by merging two existing communities
    double merge_nearest_communities();

    double compute_delta_sigma(int c1, int c2) const;       // compute delta_sigma(c1,c2)

    void remove_neighbor(Neighbor* N);
    void add_neighbor(Neighbor* N);
    void update_neighbor(Neighbor* N, double new_delta_sigma);
};

}
}       /* end of namespaces */

#endif // WALKTRAP_COMMUNITIES_H
