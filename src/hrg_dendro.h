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

#include "hrg_graph.h"
#include "hrg_rbtree.h"
#include "hrg_splittree_eq.h"

#include "igraph_hrg.h"

#include <string>
#include <cmath>

using namespace fitHRG;

namespace fitHRG {

// ***********************************************************************
// ******** Basic Structures *********************************************

#ifndef IGRAPH_HRG_LIST
#define IGRAPH_HRG_LIST

class list {
public:
    int x;            // stored elementd in linked-list
    list* next;           // pointer to next elementd
    list::list(): x(-1), next(0) { }
    list::~list() { }
};
#endif

enum {DENDRO, GRAPH, LEFT, RIGHT};
struct block {
    double x;
    int y;
};
struct ipair {
    int    x;
    int y;
    short int t;
    std::string sp;
};
struct child {
    int index;
    short int type;
    child* next;
};

// ***********************************************************************
// ******** Cnode Class **************************************************

#ifndef IGRAPH_HRG_CNODE
#define IGRAPH_HRG_CNODE
class cnode {
public:
    int index;            // array index of this node
    int degree;           // number of children in list
    int parent;           // index of parent node
    double weight;        // sampled posterior weight
    child* children;      // list of children (and their types)
    child* lastChild;     // pointer to last child in list
    cnode(): index(-1), degree(0), parent(-1), weight(0.0),
        children(0), lastChild(0)  { }
    ~cnode() {
        child *curr, *prev;
        curr = children;
        while (curr != NULL) {
            prev = curr;
            curr = curr->next;
            delete prev;
            prev = NULL;
        }
        lastChild = NULL;
    }
};
#endif

// ***********************************************************************
// ******** Split Class **************************************************

class split {
public:
    std::string s;           // partition assignment of leaf vertices
    split(): s("") { }
    ~split() { }
    void initializeSplit(const int n) {
        s = "";
        for (int i = 0; i < n; i++) {
            s += "-";
        }
    }
    bool checkSplit() {
        if (s.empty() || s.find("-", 0) != std::string::npos) {
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
private:
    ipair* edgelist;   // list of internal edges represented
    std::string* splitlist; // split representation of the internal edges
    int** indexLUT;    // table of indices of internal edges in edgelist
    int q;         // number of internal edges
    int count;         // (for adding edges) edgelist index of new edge to add
public:
    interns(const int);
    ~interns();

    // add an internal edge, O(1)
    bool addEdge(const int, const int, const short int);
    // returns the ith edge of edgelist, O(1)
    ipair* getEdge(const int);
    // returns a uniformly random internal edge, O(1)
    ipair* getRandomEdge();
    // returns the ith split of the splitlist, O(1)
    std::string getSplit(const int);
    // replace an existing split, O(1)
    bool replaceSplit(const int, const std::string);
    // swaps two edges, O(1)
    bool swapEdges(const int, const int, const short int, const int,
                   const int, const short int);
};

// ***********************************************************************
// ******** Tree elementd Class ******************************************

class elementd {
public:
    short int type; // either DENDRO or GRAPH
    double logL;    // log-likelihood contribution of this internal node
    double p;       // probability p_i that an edge exists between L and
    // R subtrees
    int e;      // number of edges between L and R subtrees
    int n;      // number of leafs in subtree rooted here
    int label;      // subtree label: smallest leaf index
    int index;      // index in containing array

    elementd *M;          // pointer to parent node
    elementd *L;          // pointer for L subtree
    elementd *R;          // pointer for R subtree

    elementd(): type(DENDRO), logL(0.0), p(0.0), e(0), n(0),
        label(-1), index(-1), M(0), L(0), R(0) { }
    ~elementd() { }
};

// ***********************************************************************
// ******** Dendrogram Class *********************************************

class dendro {
private:
    elementd* root;     // root of the dendrogram
    elementd* internal; // array of n-1 internal vertices (the dendrogram D)
    elementd* leaf;     // array of n   leaf vertices (the graph G)
    int n;          // number of leaf vertices to allocate
    interns* d;         // list of internal edges of dendrogram D
    splittree* splithist;       // histogram of cumulative split weights
    list** paths;           // array of path-lists from root to leaf
    double L;        // log-likelihood of graph G given dendrogram D
    rbtree subtreeL, subtreeR;  // trees for computeEdgeCount() function
    cnode* ctree;       // (consensus tree) array of internal tree nodes
    int* cancestor;     // (consensus tree) oldest ancetor's index for
    // each leaf

    // insert node i according to binary search property
    void binarySearchInsert(elementd*, elementd*);
    // return path to root from leaf
    list* binarySearchFind(const double);
    // build split for this internal edge
    std::string buildSplit(elementd*);
    // compute number of edges between two internal subtrees
    int computeEdgeCount(const int, const short int, const int,
                         const short int);
    // (consensus tree) counts children
    int countChildren(const std::string);
    // find internal node of D that is common ancestor of i,j
    elementd* findCommonAncestor(list**, const int, const int);
    // return reverse of path to leaf from root
    list* reversePathToRoot(const int);
// quicksort functions
    void QsortMain(block*, int, int);
    int QsortPartition(block*, int, int, int);

public:
    // underlying G (dangerously accessible)
    graph* g;

    // constructor / destructor
    dendro(); ~dendro();
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
    splittree* getConsensusSplits();
    // return likelihood of G given D
    double getLikelihood();
    // store splits in this splittree
    void getSplitList(splittree*);
    // return total weight of splittree
    double getSplitTotalWeight();
    // make random G from D
    void makeRandomGraph();
    // make single MCMC move
    bool monteCarloMove(double&, bool&, const double);
    // record consensus tree from splithist
    void recordConsensusTree(igraph_vector_t *parents,
                             igraph_vector_t *weights);
    // record D structure
    void recordDendrogramStructure(igraph_hrg_t *hrg);
    // record G structure to igraph graph
    void recordGraphStructure(igraph_t *graph);
    // force refresh of log-likelihood value
    void refreshLikelihood();
    // sample dendrogram edge likelihoods and update edge histograms
    void sampleAdjacencyLikelihoods();
    // reset the dendrograph structures
    void resetDendrograph();
    // sample dendrogram's splits and update the split histogram
    bool sampleSplitLikelihoods(int&);
    // reset splits histogram
    void resetAllSplits();
};

} // namespace fitHRG

#endif
