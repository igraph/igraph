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

// File: heap.cpp
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

#include "walktrap_heap.h"

using namespace igraph::walktrap;

void Neighbor_heap::move_up(int index) {
    while (H[index / 2]->delta_sigma > H[index]->delta_sigma) {
        Neighbor* tmp = H[index / 2];
        H[index]->heap_index = index / 2;
        H[index / 2] = H[index];
        tmp->heap_index = index;
        H[index] = tmp;
        index = index / 2;
    }
}

void Neighbor_heap::move_down(int index) {
    while (true) {
        int min = index;
        if ((2 * index < size) && (H[2 * index]->delta_sigma < H[min]->delta_sigma)) {
            min = 2 * index;
        }
        if (2 * index + 1 < size && H[2 * index + 1]->delta_sigma < H[min]->delta_sigma) {
            min = 2 * index + 1;
        }
        if (min != index) {
            Neighbor* tmp = H[min];
            H[index]->heap_index = min;
            H[min] = H[index];
            tmp->heap_index = index;
            H[index] = tmp;
            index = min;
        } else {
            break;
        }
    }
}

Neighbor* Neighbor_heap::get_first() {
    if (size == 0) {
        return 0;
    } else {
        return H[0];
    }
}

void Neighbor_heap::remove(Neighbor* N) {
    if (N->heap_index == -1 || size == 0) {
        return;
    }
    Neighbor* last_N = H[--size];
    H[N->heap_index] = last_N;
    last_N->heap_index = N->heap_index;
    move_up(last_N->heap_index);
    move_down(last_N->heap_index);
    N->heap_index = -1;
}

void Neighbor_heap::add(Neighbor* N) {
    if (size >= max_size) {
        return;
    }
    N->heap_index = size++;
    H[N->heap_index] = N;
    move_up(N->heap_index);
}

void Neighbor_heap::update(Neighbor* N) {
    if (N->heap_index == -1) {
        return;
    }
    move_up(N->heap_index);
    move_down(N->heap_index);
}

long Neighbor_heap::memory() {
    return (sizeof(Neighbor_heap) + long(max_size) * sizeof(Neighbor*));
}

Neighbor_heap::Neighbor_heap(int max_s) {
    max_size = max_s;
    size = 0;
    H = new Neighbor*[max_s];
}

Neighbor_heap::~Neighbor_heap() {
    delete[] H;
}

bool Neighbor_heap::is_empty() {
    return (size == 0);
}



//#################################################################

void Min_delta_sigma_heap::move_up(int index) {
    while (delta_sigma[H[index / 2]] < delta_sigma[H[index]]) {
        int tmp = H[index / 2];
        I[H[index]] = index / 2;
        H[index / 2] = H[index];
        I[tmp] = index;
        H[index] = tmp;
        index = index / 2;
    }
}

void Min_delta_sigma_heap::move_down(int index) {
    while (true) {
        int max = index;
        if (2 * index < size && delta_sigma[H[2 * index]] > delta_sigma[H[max]]) {
            max = 2 * index;
        }
        if (2 * index + 1 < size && delta_sigma[H[2 * index + 1]] > delta_sigma[H[max]]) {
            max = 2 * index + 1;
        }
        if (max != index) {
            int tmp = H[max];
            I[H[index]] = max;
            H[max] = H[index];
            I[tmp] = index;
            H[index] = tmp;
            index = max;
        } else {
            break;
        }
    }
}

int Min_delta_sigma_heap::get_max_community() {
    if (size == 0) {
        return -1;
    } else {
        return H[0];
    }
}

void Min_delta_sigma_heap::remove_community(int community) {
    if (I[community] == -1 || size == 0) {
        return;
    }
    int last_community = H[--size];
    H[I[community]] = last_community;
    I[last_community] = I[community];
    move_up(I[last_community]);
    move_down(I[last_community]);
    I[community] = -1;
}

void Min_delta_sigma_heap::update(int community) {
    if (community < 0 || community >= max_size) {
        return;
    }
    if (I[community] == -1) {
        I[community] = size++;
        H[I[community]] = community;
    }
    move_up(I[community]);
    move_down(I[community]);
}

long Min_delta_sigma_heap::memory() {
    return (sizeof(Min_delta_sigma_heap) + long(max_size) * (2 * sizeof(int) + sizeof(float)));
}

Min_delta_sigma_heap::Min_delta_sigma_heap(int max_s) {
    max_size = max_s;
    size = 0;
    H = new int[max_s];
    I = new int[max_s];
    delta_sigma = new float[max_s];
    for (int i = 0; i < max_size; i++) {
        I[i] = -1;
        delta_sigma[i] = 1.;
    }
}

Min_delta_sigma_heap::~Min_delta_sigma_heap() {
    delete[] H;
    delete[] I;
    delete[] delta_sigma;
}

bool Min_delta_sigma_heap::is_empty() {
    return (size == 0);
}
