/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
/*
   IGraph library.
   Copyright (C) 2011-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#ifndef INFOMAP_FLOWGRAPH_H
#define INFOMAP_FLOWGRAPH_H

#include "infomap_Node.h"

#include "igraph_datatype.h"
#include "igraph_types.h"
#include "igraph_vector.h"

#include <vector>
#include <set>
#include <cmath>

inline double plogp(double x) {
    return x > 0.0 ? x*std::log(x) : 0.0;
}

class FlowGraph {
private:
    void init(igraph_integer_t n, const igraph_vector_t *nodeWeights);

public:
    explicit FlowGraph(igraph_integer_t n);

    FlowGraph(const FlowGraph &fgraph);
    FlowGraph(const FlowGraph &fgraph, const std::vector<igraph_integer_t> &sub_members);

    FlowGraph(const igraph_t *graph, const igraph_vector_t *e_weights,
              const igraph_vector_t *v_weights);

    void swap(FlowGraph &fgraph) noexcept;

    void initiate();
    void eigenvector();
    void calibrate() noexcept;

    void back_to(const FlowGraph &fgraph);

    /*************************************************************************/
    std::vector<Node> node;
    igraph_integer_t Nnode;

    double alpha, beta;

    igraph_integer_t Ndanglings;
    std::vector<igraph_integer_t> danglings; // id of dangling nodes

    double exit;                  //
    double exitFlow;              //
    double exit_log_exit;         //
    double size_log_size;         //
    double nodeSize_log_nodeSize; // \sum_{v in V} p log(p)

    double codeLength;
};

#endif // INFOMAP_FLOWGRAPH_H
