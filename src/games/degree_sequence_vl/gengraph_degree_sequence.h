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
#include "igraph_vector.h"

namespace gengraph {

class degree_sequence {

private:
    igraph_integer_t n;
    igraph_integer_t *deg;
    igraph_integer_t total;

public :
    // #vertices
    inline igraph_integer_t size() {
        return n;
    };
    inline igraph_integer_t sum() {
        return total;
    };
    inline igraph_integer_t operator[](igraph_integer_t i) {
        return deg[i];
    };
    inline igraph_integer_t *seq() {
        return deg;
    };
    inline void assign(igraph_integer_t n0, igraph_integer_t* d0) {
        n = n0;
        deg = d0;
    };
    inline igraph_integer_t dmax() {
        igraph_integer_t dm = deg[0];
        for (igraph_integer_t i = 1; i < n; i++) if (deg[i] > dm) {
                dm = deg[i];
            }
        return dm;
    }

    degree_sequence(igraph_integer_t n, igraph_integer_t *degs);

    // igraph constructor
    degree_sequence(const igraph_vector_int_t *out_seq);

    // destructor
    ~degree_sequence();

    // compute total number of arcs
    void compute_total();

#if 0
    // raw print (vertex by vertex)
    void print();

    // distribution print (degree frequency)
    void print_cumul();
#endif

    // is degree sequence realizable ?
    bool havelhakimi();

};

} // namespace gengraph

#endif //DEGREE_SEQUENCE_H
