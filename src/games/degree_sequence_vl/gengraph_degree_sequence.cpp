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
#include "gengraph_definitions.h"
#include "gengraph_degree_sequence.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <vector>
#include <stdexcept>

// using namespace __gnu_cxx;
using namespace std;

namespace gengraph {

// shuffle an igraph_integer_t[] randomly
void random_permute(igraph_integer_t *a, igraph_integer_t n);

// sort an array of positive integers in time & place O(n + max)
void cumul_sort(igraph_integer_t *q, igraph_integer_t n);


degree_sequence::~degree_sequence() {
    deg = NULL;
}

void degree_sequence::compute_total() {
    total = 0;
    for (igraph_integer_t i = 0; i < n; i++) {
        total += deg[i];
    }
}

degree_sequence::
degree_sequence(igraph_integer_t n0, igraph_integer_t *degs) {
    deg = degs;
    n = n0;
    compute_total();
}

degree_sequence::
degree_sequence(const igraph_vector_int_t *out_seq) {
    n = igraph_vector_int_size(out_seq);
    deg = &VECTOR(*out_seq)[0];
    compute_total();
}

#ifndef FBUFF_SIZE
    #define FBUFF_SIZE 999
#endif //FBUFF_SIZE

bool degree_sequence::havelhakimi() {

    igraph_integer_t i;
    igraph_integer_t dm = dmax() + 1;
    // Sort vertices using basket-sort, in descending degrees
    igraph_integer_t *nb = new igraph_integer_t[dm];
    igraph_integer_t *sorted = new igraph_integer_t[n];
    // init basket
    for (i = 0; i < dm; i++) {
        nb[i] = 0;
    }
    // count basket
    for (i = 0; i < n; i++) {
        nb[deg[i]]++;
    }
    // cumul
    igraph_integer_t c = 0;
    for (i = dm - 1; i >= 0; i--) {
        igraph_integer_t t = nb[i];
        nb[i] = c;
        c += t;
    }
    // sort
    for (i = 0; i < n; i++) {
        sorted[nb[deg[i]]++] = i;
    }

// Binding process starts
    igraph_integer_t first = 0;  // vertex with biggest residual degree
    igraph_integer_t d = dm - 1; // maximum residual degree available

    for (c = total / 2; c > 0; ) {
        // We design by 'v' the vertex of highest degree (indexed by first)
        // look for current degree of v
        while (nb[d] <= first) {
            d--;
        }
        // store it in dv
        igraph_integer_t dv = d;
        // bind it !
        c -= dv;
        igraph_integer_t dc = d;         // residual degree of vertices we bind to
        igraph_integer_t fc = ++first;   // position of the first vertex with degree dc

        while (dv > 0 && dc > 0) {
            igraph_integer_t lc = nb[dc];
            if (lc != fc) {
                while (dv > 0 && lc > fc) {
                    // binds v with sorted[--lc]
                    dv--;
                    lc--;
                }
                fc = nb[dc];
                nb[dc] = lc;
            }
            dc--;
        }
        if (dv != 0) { // We couldn't bind entirely v
            delete[] nb;
            delete[] sorted;
            return false;
        }
    }
    delete[] nb;
    delete[] sorted;
    return true;
}

} // namespace gengraph
