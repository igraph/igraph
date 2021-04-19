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
#ifndef DEGREE_SEQUENCE_H
#define DEGREE_SEQUENCE_H

#include "igraph_types.h"
#include "igraph_datatype.h"

namespace gengraph {

class degree_sequence {

private:
    int n;
    int * deg;
    int total;

public :
    // #vertices
    inline int size() {
        return n;
    };
    inline int sum() {
        return total;
    };
    inline int operator[](int i) {
        return deg[i];
    };
    inline int *seq() {
        return deg;
    };
    inline void assign(int n0, int* d0) {
        n = n0;
        deg = d0;
    };
    inline int dmax() {
        int dm = deg[0];
        for (int i = 1; i < n; i++) if (deg[i] > dm) {
                dm = deg[i];
            }
        return dm;
    }

    void make_even(int mini = -1, int maxi = -1);
    void sort();
    void shuffle();

    // raw constructor
    degree_sequence(int n, int *degs);

    // read-from-file constrictor
    degree_sequence(FILE *f, bool DISTRIB = true);

    // simple power-law constructor : Pk = int((x+k0)^(-exp),x=k..k+1), with k0 so that avg(X)=z
    degree_sequence(int n, double exp, int degmin, int degmax, double avg_degree = -1.0);

    // igraph constructor
    degree_sequence(const igraph_vector_t *out_seq);

    // destructor
    ~degree_sequence();

    // unbind the deg[] vector (so that it doesn't get deleted when the class is destroyed)
    void detach();

    // compute total number of arcs
    void compute_total();

    // raw print (vertex by vertex)
    void print();

    // distribution print (degree frequency)
    void print_cumul();

    // is degree sequence realizable ?
    bool havelhakimi();

};

} // namespace gengraph

#endif //DEGREE_SEQUENCE_H

