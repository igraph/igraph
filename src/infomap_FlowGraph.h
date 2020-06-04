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

#ifndef FLOWGRAPH_H
#define FLOWGRAPH_H

#include <vector>
#include <set>

#include "igraph_interface.h"

#include "infomap_Node.h"

class FlowGraph {
private:
    void init(int n, const igraph_vector_t *nodeWeights);

public:
    FlowGraph(int n);
    FlowGraph(int n, const igraph_vector_t *nodeWeights);
    FlowGraph(FlowGraph * fgraph);
    FlowGraph(FlowGraph * fgraph, int sub_Nnode, int * sub_members);

    FlowGraph(const igraph_t * graph, const igraph_vector_t *e_weights,
              const igraph_vector_t *v_weights);

    ~FlowGraph();

    void swap(FlowGraph * fgraph);

    void initiate();
    void eigenvector();
    void calibrate();

    void back_to(FlowGraph * fgraph);

    /*************************************************************************/
    Node **node;
    int  Nnode;

    double alpha, beta;

    int Ndanglings;
    std::vector<int> danglings; // id of dangling nodes

    double exit;                  //
    double exitFlow;              //
    double exit_log_exit;         //
    double size_log_size;         //
    double nodeSize_log_nodeSize; // \sum_{v in V} p log(p)

    double codeLength;
};

void delete_FlowGraph(FlowGraph *fgraph);

#endif
