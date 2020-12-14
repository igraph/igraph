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

#ifndef GREEDY_H
#define GREEDY_H

#include <vector>
#include <map>
#include <utility>
#include <climits>

#include "igraph_random.h"

#include "infomap_Node.h"
#include "infomap_FlowGraph.h"

class Greedy {
public:
    Greedy(FlowGraph * fgraph);
    // initialise les attributs par rapport au graph

    ~Greedy();

    void setMove(int *moveTo);
    //virtual void determMove(int *moveTo);

    bool optimize();
    //virtual void move(bool &moved);

    void apply(bool sort);
    //virtual void level(Node ***, bool sort);

    void tune(void);

    /**************************************************************************/

    FlowGraph * graph;
    int Nnode;

    double exit;
    double exitFlow;
    double exit_log_exit;
    double size_log_size;
    double nodeSize_log_nodeSize;

    double codeLength;

    double alpha, beta;
    // local copy of fgraph alpha, beta (=alpha -  Nnode = graph->Nnode;1)

    std::vector<int> node_index;  // module number of each node

    int Nempty;
    std::vector<int> mod_empty;

    std::vector<double> mod_exit;  // version tmp de node
    std::vector<double> mod_size;
    std::vector<double> mod_danglingSize;
    std::vector<double> mod_teleportWeight;
    std::vector<int> mod_members;
};

void delete_Greedy(Greedy *greedy);
#endif
