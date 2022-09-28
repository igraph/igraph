/* -*- mode: C++ -*-  */
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

// ****************************************************************************************************
// *** COPYRIGHT NOTICE *******************************************************************************
// graph_simp.h - graph data structure
// Copyright (C) 2006-2008 Aaron Clauset
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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// See http://www.gnu.org/licenses/gpl.txt for more details.
//
// ****************************************************************************************************
// Author       : Aaron Clauset  ( aaronc@santafe.edu | http://www.santafe.edu/~aaronc/ )
// Collaborators: Cristopher Moore and Mark E.J. Newman
// Project      : Hierarchical Random Graphs
// Location     : University of New Mexico, Dept. of Computer Science AND Santa Fe Institute
// Created      : 21 June 2006
// Modified     : 23 December 2007 (cleaned up for public consumption)
//
// ************************************************************************
//
// Simple graph data structure. The basic structure is an adjacency
// list of edges, along with degree information for the vertices.
//
// ************************************************************************

#ifndef IGRAPH_HRG_SIMPLEGRAPH
#define IGRAPH_HRG_SIMPLEGRAPH

#include "hrg/rbtree.h"
#include "hrg/dendro.h"

#include <cstring>
#include <cstdlib>

namespace fitHRG {

// ******** Basic Structures *********************************************

#ifndef IGRAPH_HRG_SIMPLEEDGE
#define IGRAPH_HRG_SIMPLEEDGE
class simpleEdge {
public:
    int x;            // index of edge terminator
    simpleEdge* next;     // pointer to next elementd

    simpleEdge(): x(-1), next(0) { }
    ~simpleEdge() { }
};
#endif

#ifndef IGRAPH_HRG_SIMPLEVERT
#define IGRAPH_HRG_SIMPLEVERT
class simpleVert {
public:
    std::string name;          // (external) name of vertex
    int degree;           // degree of this vertex
    int group_true;       // index of vertex's true group

    simpleVert(): name(""), degree(0), group_true(-1) { }
    ~simpleVert() { }
};
#endif

#ifndef IGRAPH_HRG_TWOEDGE
#define IGRAPH_HRG_TWOEDGE
class twoEdge {
public:
    int o;            // index of edge originator
    int x;            // index of edge terminator

    twoEdge(): o(-1), x(-1) { }
    ~twoEdge() { }
};
#endif

// ******** Graph Class with Edge Statistics *****************************

class simpleGraph {
public:
    simpleGraph(const int); ~simpleGraph();

    // add group label to vertex i
    bool addGroup(const int, const int);
    // add (i,j) to graph
    bool addLink(const int, const int);
    // true if (i,j) is already in graph
    bool doesLinkExist(const int, const int);
    // returns A(i,j)
    double getAdjacency(const int, const int);
    // returns degree of vertex i
    int getDegree(const int);
    // returns group label of vertex i
    int getGroupLabel(const int);
    // returns name of vertex i
    std::string getName(const int);
    // returns edge list of vertex i
    simpleEdge* getNeighborList(const int);
    // return pointer to a node
    simpleVert* getNode(const int);
    // returns num_groups
    int getNumGroups();
    // returns m
    int getNumLinks();
    // returns n
    int getNumNodes();
    // set name of vertex i
    bool setName(const int, const std::string);

private:
    simpleVert* nodes;        // list of nodes
    simpleEdge** nodeLink;    // linked list of neighbors to vertex
    simpleEdge** nodeLinkTail;    // pointers to tail of neighbor list
    double** A;           // adjacency matrix for this graph
    twoEdge* E;           // list of all edges (array)
    int n;            // number of vertices
    int m;            // number of directed edges
    int num_groups;       // number of bins in node histograms

    // quicksort functions
    void QsortMain(block*, int, int);
    int QsortPartition(block*, int, int, int);
};

} // namespace fitHRG

#endif
