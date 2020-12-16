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
// graph.h - graph data structure for hierarchical random graphs
// Copyright (C) 2005-2008 Aaron Clauset
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
// Created      : 8 November 2005
// Modified     : 23 December 2007 (cleaned up for public consumption)
//
// ****************************************************************************************************
//
// Graph data structure for hierarchical random graphs. The basic structure is an adjacency list of
// edges; however, many additional pieces of metadata are stored as well. Each node stores its
// external name, its degree and (if assigned) its group index.
//
// ****************************************************************************************************

#ifndef IGRAPH_HRG_GRAPH
#define IGRAPH_HRG_GRAPH

#include "hrg/rbtree.h"

#include <string>
#include <cstring>
#include <cstdlib>

namespace fitHRG {

// ******** Basic Structures *********************************************

#ifndef IGRAPH_HRG_EDGE
#define IGRAPH_HRG_EDGE
class edge {
public:
    int x;            // stored integer value  (edge terminator)
    double* h;            // (histogram) weights of edge existence
    double total_weight;      // (histogram) total weight observed
    int obs_count;        // number of observations in histogram
    edge* next;           // pointer to next elementd
    edge(): x(-1), h(0), total_weight(0.0), obs_count(0), next(0)  { }
    ~edge() {
        if (h != NULL) {
            delete [] h;
        }
        h = NULL;
    }
};
#endif

#ifndef IGRAPH_HRG_VERT
#define IGRAPH_HRG_VERT
class vert {
public:
    std::string name;           // (external) name of vertex
    int degree;            // degree of this vertex

    vert(): name(""), degree(0) { }
    ~vert() { }
};
#endif

// ******** Graph Class with Edge Statistics *****************************

class graph {
public:
    graph(const int, bool predict = false);
    ~graph();

    // add (i,j) to graph
    bool addLink(const int, const int);
    // add weight to (i,j)'s histogram
    bool addAdjacencyObs(const int, const int, const double, const double);
    // add to obs_count and total_weight
    void addAdjacencyEnd();
    // true if (i,j) is already in graph
    bool doesLinkExist(const int, const int);
    // returns degree of vertex i
    int getDegree(const int);
    // returns name of vertex i
    std::string getName(const int);
    // returns edge list of vertex i
    edge* getNeighborList(const int);
    // return ptr to histogram of edge (i,j)
    double* getAdjacencyHist(const int, const int);
    // return average value of adjacency A(i,j)
    double getAdjacencyAverage(const int, const int);
    // returns bin_resolution
    double getBinResolution();
    // returns num_bins
    int getNumBins();
    // returns m
    int numLinks();
    // returns n
    int numNodes();
    // returns total_weight
    double getTotalWeight();
    // reset edge (i,j)'s histogram
    void resetAdjacencyHistogram(const int, const int);
    // reset all edge histograms
    void resetAllAdjacencies();
    // clear all links from graph
    void resetLinks();
    // allocate edge histograms
    void setAdjacencyHistograms(const int);
    // set name of vertex i
    bool setName(const int, const std::string);

private:
    bool predict;      // do we need prediction?
    vert* nodes;       // list of nodes
    edge** nodeLink;   // linked list of neighbors to vertex
    edge** nodeLinkTail;   // pointers to tail of neighbor list
    double*** A;       // stochastic adjacency matrix for this graph
    int obs_count;     // number of observations in A
    double total_weight;   // total weight added to A
    int n;         // number of vertices
    int m;         // number of directed edges
    int num_bins;      // number of bins in edge histograms
    double bin_resolution; // width of histogram bin
};

} // namespace fitHRG

#endif
