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
// dendro_eq.h - hierarchical random graph (hrg) data structure
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
// Created      : 19 April 2006
// Modified     : 19 May 2007
//           : 19 May 2008 (cleaned up for public consumption)
//
// ****************************************************************************************************
//
// Maximum likelihood dendrogram data structure. This is the heart of the HRG algorithm: all
// manipulations are done here and all data is stored here. The data structure uses the separate
// graph data structure to store the basic adjacency information (in a dangerously mutable way).
//
// Note: This version (dendro_eq.h) differs from other versions because it includes methods for
//       doing the consensus dendrogram calculation.
//
// ****************************************************************************************************

#ifndef IGRAPH_HRG_DENDRO
#define IGRAPH_HRG_DENDRO

#include "hrg/graph.h"
#include "hrg/rbtree.h"
#include "hrg/splittree_eq.h"

#include "igraph_hrg.h"

#include <string>
#include <cmath>

using std::string;

namespace fitHRG {

// ***********************************************************************
// ******** Basic Structures *********************************************

enum {DENDRO, GRAPH, LEFT, RIGHT};

struct block {
    double x;
    int y;
};

struct ipair {
    int x;
    int y;
    short int t;
    string sp;
};

struct child {
    int index;
    short int type;
    child* next;
};

// ***********************************************************************
// ******** Cnode Class **************************************************

struct cnode {
    int index = -1;             // array index of this node
    int degree = 0;             // number of children in list
    int parent = -1;            // index of parent node
    double weight = 0.0;        // sampled posterior weight
    child* children = nullptr;  // list of children (and their types)
    child* lastChild = nullptr; // pointer to last child in list
    cnode() = default;
    cnode(const cnode &) = delete;
    cnode & operator = (const cnode &) = delete;
    ~cnode() {
        child *curr = children;
        while (curr != nullptr) {
            child *prev = curr;
            curr = curr->next;
            delete prev;
        }
        lastChild = nullptr;
    }
};

// ***********************************************************************
// ******** Split Class **************************************************

class split {
public:
    string s;           // partition assignment of leaf vertices
    split() = default;
    void initializeSplit(const int n) {
        s = "";
        for (int i = 0; i < n; i++) {
            s += "-";
        }
    }
    bool checkSplit() const {
        if (s.empty() || s.find('-', 0) != string::npos) {
            return false;
        } else {
            return true;
        }
    }
};

// ***********************************************************************
// ******** Internal Edge Class ******************************************
// The usefulness of this data structure is to provide an easy to way
// maintain the set of internal edges, and the corresponding splits,
// in the dendrogram D. It allows for the selection of a random
// internal edge in O(1) time, and it takes O(1) time to update its
// structure given an internal move. This structure does not provide
// any means to directly manipulate the splits, but does allow them to
// be replaced. A split has the form "int.int...int#int.int...int",
// where all ints on the left side of the # are in the left partition
// and all ints on the right side of the # marker are in the right
// partition defined by the split.

class interns {
    ipair* edgelist;   // list of internal edges represented
    string* splitlist; // split representation of the internal edges
    int** indexLUT;    // table of indices of internal edges in edgelist
    int q;         // number of internal edges
    int count;         // (for adding edges) edgelist index of new edge to add
public:
    explicit interns(int);
    ~interns();

    // add an internal edge, O(1)
    bool addEdge(int, int, short int);
    // returns the ith edge of edgelist, O(1)
    ipair* getEdge(int) const;
    // returns a uniformly random internal edge, O(1)
    ipair* getRandomEdge() const;
    // returns the ith split of the splitlist, O(1)
    string getSplit(int) const;
    // replace an existing split, O(1)
    bool replaceSplit(int i, const string &sp);
    // swaps two edges, O(1)
    bool swapEdges(int, int, short int, int, int, short int);
};

// ***********************************************************************
// ******** Tree elementd Class ******************************************

struct elementd {
    short int type = DENDRO;    // either DENDRO or GRAPH
    double logL = 0.0;          // log-likelihood contribution of this internal node
    double p = 0.0;             // probability p_i that an edge exists between L and
    // R subtrees
    int e = 0;      // number of edges between L and R subtrees
    int n = 0;      // number of leafs in subtree rooted here
    int label = -1; // subtree label: smallest leaf index
    int index = -1; // index in containing array

    elementd *M = nullptr;  // pointer to parent node
    elementd *L = nullptr;  // pointer for L subtree
    elementd *R = nullptr;  // pointer for R subtree
};

// ***********************************************************************
// ******** Dendrogram Class *********************************************

class dendro {
    elementd* root = nullptr;       // root of the dendrogram
    elementd* internal = nullptr;   // array of n-1 internal vertices (the dendrogram D)
    elementd* leaf = nullptr;       // array of n   leaf vertices (the graph G)
    int n;                          // number of leaf vertices to allocate
    interns* d = nullptr;           // list of internal edges of dendrogram D
    splittree* splithist = nullptr; // histogram of cumulative split weights
    list** paths = nullptr;         // array of path-lists from root to leaf
    double L;                       // log-likelihood of graph G given dendrogram D
    rbtree subtreeL, subtreeR;      // trees for computeEdgeCount() function
    cnode* ctree = nullptr;         // (consensus tree) array of internal tree nodes
    int* cancestor = nullptr;       // (consensus tree) oldest ancetor's index for
    // each leaf

    // insert node i according to binary search property
    void binarySearchInsert(elementd*, elementd*);
    // return path to root from leaf
    list* binarySearchFind(double) const;
    // build split for this internal edge
    string buildSplit(elementd*) const;
    // compute number of edges between two internal subtrees
    int computeEdgeCount(int, short int, int, short int);
    // (consensus tree) counts children
    static size_t countChildren(const string &s);
    // find internal node of D that is common ancestor of i,j
    elementd* findCommonAncestor(list**, int, int);
    // return reverse of path to leaf from root
    list* reversePathToRoot(int);
    // quicksort functions
    static void QsortMain(block*, int, int);
    static int QsortPartition(block*, int, int, int);

    // underlying G (dangerously accessible)
    graph* g = nullptr;

public:

    // constructor / destructor
    dendro() = default;
    dendro(const dendro &) = delete;
    dendro & operator = (const dendro &) = delete;
    ~dendro();

    igraph_error_t setGraph(const igraph_t *igraph);
    void setGraph(graph *ig) { g = ig; }
    const graph *getGraph() const { return g; }

    // build dendrogram from g
    void buildDendrogram();
    // delete dendrograph in prep for importDendrogramStructure
    void clearDendrograph();
    // read dendrogram structure from HRG structure
    bool importDendrogramStructure(const igraph_hrg_t *hrg);
    // (consensus tree) delete splits with less than 0.5 weight
    void cullSplitHist();
    // return size of consensus split
    int getConsensusSize();
    // return split tree with consensus splits
    splittree* getConsensusSplits() const;
    // return likelihood of G given D
    double getLikelihood() const;
    // store splits in this splittree
    void getSplitList(splittree*) const;
    // return total weight of splittree
    double getSplitTotalWeight() const;
    // make random G from D
    void makeRandomGraph();
    // make single MCMC move
    void monteCarloMove(double &, bool &, double);
    // record consensus tree from splithist
    void recordConsensusTree(igraph_vector_int_t *parents,
                             igraph_vector_t *weights);
    // record D structure
    void recordDendrogramStructure(igraph_hrg_t *hrg) const noexcept;
    // record G structure to igraph graph
    igraph_error_t recordGraphStructure(igraph_t *graph) const noexcept;
    // force refresh of log-likelihood value
    void refreshLikelihood();
    // sample dendrogram edge likelihoods and update edge histograms
    void sampleAdjacencyLikelihoods();
    // reset the dendrograph structures
    void resetDendrograph();
    // sample dendrogram's splits and update the split histogram
    bool sampleSplitLikelihoods();
};

} // namespace fitHRG

#endif
