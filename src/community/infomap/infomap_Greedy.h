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

#ifndef INFOMAP_GREEDY_H
#define INFOMAP_GREEDY_H

#include "infomap_Node.h"
#include "infomap_FlowGraph.h"

#include "igraph_random.h"

#include <vector>

class Greedy {
public:
    explicit Greedy(FlowGraph *fgraph);
    // initialise les attributs par rapport au graph

    void setMove(const std::vector<igraph_integer_t> &moveTo);
    bool optimize();
    void apply(bool sort);

    /**************************************************************************/

public:
    double codeLength;

private:
    FlowGraph * graph;
    igraph_integer_t Nnode;

    double exit;
    double exitFlow;
    double exit_log_exit;
    double size_log_size;
    double nodeSize_log_nodeSize;

    double alpha, beta;
    // local copy of fgraph alpha, beta (=alpha -  Nnode = graph->Nnode;1)

    std::vector<igraph_integer_t> node_index;  // module number of each node

    igraph_integer_t Nempty = 0;
    std::vector<igraph_integer_t> mod_empty;

    std::vector<double> mod_exit;  // version tmp de node
    std::vector<double> mod_size;
    std::vector<double> mod_danglingSize;
    std::vector<double> mod_teleportWeight;
    std::vector<size_t> mod_members;
};

#endif // INFOMAP_GREEDY_H
