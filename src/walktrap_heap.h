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

// File: heap.h
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
// Email    : pons@liafa.jussieu.fr
// Web page : http://www.liafa.jussieu.fr/~pons/
// Location : Paris, France
// Time     : June 2005
//-----------------------------------------------------------------------------
// see readme.txt for more details

#ifndef WALKTRAP_HEAP_H
#define WALKTRAP_HEAP_H

namespace igraph {

namespace walktrap {

class Neighbor {
public:
    int community1;   // the two adjacent communities
    int community2;   // community1 < community2

    float delta_sigma;    // the delta sigma between the two communities
    float weight;     // the total weight of the edges between the two communities
    bool exact;       // true if delta_sigma is exact, false if it is only a lower bound

    Neighbor* next_community1;        // pointers of two double
    Neighbor* previous_community1;    // chained lists containing
    Neighbor* next_community2;        // all the neighbors of
    Neighbor* previous_community2;    // each communities.

    int heap_index;   //

    Neighbor();
};


class Neighbor_heap {
private:
    int size;
    int max_size;

    Neighbor** H;   // the heap that contains a pointer to each Neighbor object stored

    void move_up(int index);
    void move_down(int index);

public:
    void add(Neighbor* N);        // add a new distance
    void update(Neighbor* N);     // update a distance
    void remove(Neighbor* N);     // remove a distance
    Neighbor* get_first();        // get the first item
    long memory();
    bool is_empty();

    Neighbor_heap(int max_size);
    ~Neighbor_heap();
};


class Min_delta_sigma_heap {
private:
    int size;
    int max_size;

    int* H;   // the heap that contains the number of each community
    int* I;   // the index of each community in the heap (-1 = not stored)

    void move_up(int index);
    void move_down(int index);

public:
    int get_max_community();              // return the community with the maximal delta_sigma
    void remove_community(int community);         // remove a community;
    void update(int community);               // update (or insert if necessary) the community
    long memory();                    // the memory used in Bytes.
    bool is_empty();

    float* delta_sigma;                    // the delta_sigma of the stored communities

    Min_delta_sigma_heap(int max_size);
    ~Min_delta_sigma_heap();
};

}
}        /* end of namespaces */

#endif // WALKTRAP_HEAP_H

