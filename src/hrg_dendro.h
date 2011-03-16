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
//			 : 19 May 2008 (cleaned up for public consumption)
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

#if !defined(dendro_INCLUDED)
#define dendro_INCLUDED

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <math.h>

#include "hrg_MersenneTwister.h"
#include "hrg_graph.h"
#include "hrg_rbtree.h"
#include "hrg_splittree_eq.h"

#include "igraph_hrg.h"

using namespace std;
using namespace fitHRG;

namespace fitHRG {

// ********************************************************************************************************
// ******** Basic Structures ******************************************************************************

#if !defined(list_INCLUDED)
#define list_INCLUDED
class list {
public:
	int		x;				// stored elementd in linked-list
	list*	next;			// pointer to next elementd
	list(); ~list();
};
list::list()  { x = -1; next = NULL; }
list::~list() {}
#endif

enum {DENDRO, GRAPH, LEFT, RIGHT};
struct block { double x; int y; };
struct ipair { int    x; int y; short int t; string sp; };
struct child { int index; short int type; child* next; };

// ********************************************************************************************************
// ******** Cnode Class ***********************************************************************************

#if !defined(cnode_INCLUDED)
#define cnode_INCLUDED
class cnode {
public:
	int		index;			// array index of this node
	int		degree;			// number of children in list
	int		parent;			// index of parent node
	double	weight;			// sampled posterior weight
	child*	children;			// list of children (and their types)
	child*	lastChild;		// pointer to last child in list
	cnode()  { index = -1; degree = 0; parent = -1; weight = 0.0; children = NULL; lastChild = NULL; }
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

// ********************************************************************************************************
// ******** Split Class ***********************************************************************************

class split {
public:
	string s;							// partition assignment of leaf vertices
	split() { s = ""; }
	~split() {}
	void initializeSplit(const int n) { s = ""; for (int i=0; i<n; i++) { s += "-"; } }
	bool checkSplit() { if (s == "" or s.find("-",0) != string::npos) { return false; } else { return true; } }
};

// ********************************************************************************************************
// ******** Internal Edge Class ***************************************************************************
// The usefulness of this data structure is to provide an easy to way maintain the set of internal edges,
// and the corresponding splits, in the dendrogram D. It allows for the selection of a random internal 
// edge in O(1) time, and it takes O(1) time to update its structure given an internal move. This structure
// does not provide any means to directly manipulate the splits, but does allow them to be replaced.
// A split has the form "int.int...int#int.int...int", where all ints on the left side of the # are in the left
// partition and all ints on the right side of the # marker are in the right partition defined by the split.

class interns {
private: 
	ipair*    edgelist;					// list of internal edges represented
	string*	splitlist;				// split representation of the internal edges
	int**     indexLUT;					// table of indices of internal edges in edgelist
	int		q;						// number of internal edges
	int		count;					// (for adding edges) edgelist index of new edge to add
	MTRand    mtr;						// Mersenne Twister random number generator instance

public:
	interns(const int); ~interns();
	
	bool		addEdge(const int, const int, const short int);	// add an internal edge, O(1)
	ipair*    getEdge(const int);							// returns the ith edge of edgelist, O(1)
	ipair*    getRandomEdge();							// returns a uniformly random internal edge, O(1)
	string	getSplit(const int);						// returns the ith split of the splitlist, O(1)
	void		printEdgeList();							// writes edgelist to terminal
	void		printSplitList();							// writes splitlist to terminal
	bool		replaceSplit(const int, const string);			// replace an existing split, O(1)
	bool		swapEdges(const int, const int, const short int, const int, const int, const short int);
													// swaps two edges, O(1)
};

// ********************************************************************************************************
// ******** Tree elementd Class ***************************************************************************

class elementd {
public:
	short int type;				// either DENDRO or GRAPH
	double    logL;				// log-likelihood contribution of this internal node
	double    p;					// probability p_i that an edge exists between L and R subtrees
	int		e;					// number of edges between L and R subtrees
	int		n;					// number of leafs in subtree rooted here
	int		label;				// subtree label: smallest leaf index
	int		index;				// index in containing array
	
	elementd   *M;					// pointer to parent node
	elementd   *L;					// pointer for L subtree
	elementd   *R;					// pointer for R subtree

	elementd()  {
	  type		    = DENDRO;
	  logL  = p     = 0.0;
	  e     = n	    = 0;
	  label = index = -1;
	  M     = L = R = NULL;
	}
	~elementd() {}
};

// ********************************************************************************************************
// ******** Dendrogram Class ******************************************************************************

class dendro {
	
private:
	elementd*		root;			// root of the dendrogram
	elementd*		internal;			// array of n-1 internal vertices (the dendrogram D)
	elementd*		leaf;			// array of n   leaf vertices (the graph G)
	int			n;				// number of leaf vertices to allocate
	interns*		d;				// list of internal edges of dendrogram D
	splittree*	splithist;		// histogram of cumulative split weights
	list**		paths;			// array of path-lists from root to leaf
	double		L;				// log-likelihood of graph G given dendrogram D
	MTRand		mtr;				// Mersenne Twister random number generator instance
	rbtree		subtreeL, subtreeR;	// trees for computeEdgeCount() function
	cnode*		ctree;			// (consensus tree) array of internal tree nodes
	int*			cancestor;		// (consensus tree) oldest ancetor's index for each leaf
	
	void			binarySearchInsert(elementd*, elementd*);							// insert node i according to binary search property
	list*		binarySearchFind(const double);									// return path to root from leaf
	string		buildSplit(elementd*);											// build split for this internal edge
	int			computeEdgeCount(const int, const short int, const int, const short int);  // compute number of edges between two internal subtrees
	int			countChildren(const string);										// (consensus tree) counts children
	elementd*		findCommonAncestor(list**, const int, const int);						// find internal node of D that is common ancestor of i,j
	void			printSubTree(elementd*);											// display the subtree rooted at z
	void			printConsensusTree();											// print consensus tree
	void			printConsensusTreeDense();										// print consensus tree (verbose)
	list*		reversePathToRoot(const int);										// return reverse of path to leaf from root
	void			QsortMain(block*, int, int);										// quicksort functions
	int			QsortPartition(block*, int, int, int);

public:
	graph*		g;									// underlying G (dangerously accessible)
	
	dendro(); ~dendro();								// constructor / destructor
	void			buildDendrogram();						// build dendrogram from g
	void			clearDendrograph();						// delete dendrograph in prep for importDendrogramStructure
	bool			importDendrogramStructure(const igraph_hrg_t *hrg);	// read dendrogram structure from HRG structure
	void			cullSplitHist();						// (consensus tree) delete splits with less than 0.5 weight
	int			getConsensusSize();						// return size of consensus split
	splittree*	getConsensusSplits();					// return split tree with consensus splits
	double		getLikelihood();						// return likelihood of G given D
	void			getSplitList(splittree*);				// store splits in this splittree
	double		getSplitTotalWeight();					// return total weight of splittree
	bool			importDendrogramStructure(const string);	// read dendrogram structure from file
	void            makeRandomGraph();
	bool			monteCarloMove(double&, bool&, const double);// make single MCMC move
	void			recordConsensusTree(const string);			// record consensus tree from splithist
	void                    recordConsensusTree(igraph_vector_t *parents, igraph_vector_t *weights);
	void			recordDendrogramStructure(const string);	// record D structure to file
	void                    recordDendrogramStructure(igraph_hrg_t *hrg);
	void			recordGraphStructure(const string);		// record G structure to file
	void			recordGraphStructure(igraph_t *graph);		// record G structure to igraph graph
	void			recordSplitHistogram(const string);		// record split histogram
	void			refreshLikelihood();					// force refresh of log-likelihood value
	void			sampleAdjacencyLikelihoods();				// sample dendrogram edge likelihoods and update edge histograms
	void                    resetDendrograph();                                             // reset the dendrograph structures
	bool			sampleSplitLikelihoods(int&);				// sample dendrogram's splits and update the split histogram
	void			resetAllSplits();						// reset splits histogram
	void			printDendrogram();						// write dendrogram structure to terminal
	void			printSplitStats();						// write split statistics to terminal
	void			printSplitStatsShort();					// write (short) split statistics to terminal
};

} // namespace fitHRG

#endif
