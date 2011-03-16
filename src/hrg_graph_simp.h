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
// ****************************************************************************************************
// 
// Simple graph data structure. The basic structure is an adjacency list of edges, along with degree
// information for the vertices.
// 
// ****************************************************************************************************

#if !defined(simpleGraph_INCLUDED)
#define simpleGraph_INCLUDED

#include <stdio.h>
#include <string>
#include "stdlib.h"

#include "hrg_MersenneTwister.h"
#include "hrg_rbtree.h"

using namespace std;

namespace fitHRG {

// ******** Basic Structures ******************************************************************************

#if !defined(simpleEdge_INCLUDED)
#define simpleEdge_INCLUDED
class simpleEdge {
public:
	int		x;						// index of edge terminator
	simpleEdge*	next;				// pointer to next elementd

	simpleEdge(); ~simpleEdge();
};
simpleEdge::simpleEdge()  { x = -1; next = NULL;  }
simpleEdge::~simpleEdge() {}
#endif

#if !defined(simpleVert_INCLUDED)
#define simpleVert_INCLUDED
class simpleVert {
public:
	string	name;					// (external) name of vertex
	int		degree;					// degree of this vertex
	int		group_true;				// index of vertex's true group
	
	simpleVert(); ~simpleVert();
};
simpleVert::simpleVert()  { name = ""; degree = 0; group_true = -1; }
simpleVert::~simpleVert() { }
#endif

#if !defined(twoEdge_INCLUDED)
#define twoEdge_INCLUDED
class twoEdge {
public:
	int		o;						// index of edge originator
	int		x;						// index of edge terminator
	
	twoEdge(); ~twoEdge();
};
twoEdge::twoEdge()  { o = -1; x = -1; }
twoEdge::~twoEdge() {}
#endif

// ******** Graph Class with Edge Statistics *************************************************************

class simpleGraph {
public:
	simpleGraph(const int); ~simpleGraph();

	bool			addGroup(const int, const int);					// add group label to vertex i
	bool			addLink(const int, const int);					// add (i,j) to graph
	bool			doesLinkExist(const int, const int);				// true if (i,j) is already in graph
	double		getAdjacency(const int, const int);				// returns A(i,j)
	int			getDegree(const int);							// returns degree of vertex i
	int			getGroupLabel(const int);						// returns group label of vertex i
	string		getName(const int);								// returns name of vertex i
	simpleEdge*	getNeighborList(const int);						// returns edge list of vertex i
	simpleVert*	getNode(const int);								// return pointer to a node
	int			getNumGroups();								// returns num_groups
	int			getNumLinks();									// returns m
	int			getNumNodes();									// returns n
	bool			setName(const int, const string);					// set name of vertex i

	void			printAdjacencies();								// print adjacency matrix
	void			printPairs();									// prints all edges in graph
	
private:
	simpleVert*	nodes;			// list of nodes
	simpleEdge**	nodeLink;			// linked list of neighbors to vertex
	simpleEdge**	nodeLinkTail;		// pointers to tail of neighbor list
	double**		A;				// adjacency matrix for this graph
	twoEdge*		E;				// list of all edges (array)
	int			n;				// number of vertices
	int			m;				// number of directed edges
	int			num_groups;		// number of bins in node histograms
	MTRand		mtr;				// Mersenne Twister random number generator instance

	void			QsortMain(block*, int, int);						// quicksort functions
	int			QsortPartition(block*, int, int, int);
};

// ******** Constructor / Destructor **********************************************************************

simpleGraph::simpleGraph(const int size)  {
	n			= size;
	m			= 0;
	num_groups	= 0;
	nodes		= new simpleVert  [n];
	nodeLink		= new simpleEdge* [n];
	nodeLinkTail   = new simpleEdge* [n];
	A			= new double* [n];
	for (int i=0; i<n; i++) {
		nodeLink[i] = NULL; nodeLinkTail[i] = NULL;
		A[i] = new double [n];
		for (int j=0; j<n; j++) { A[i][j] = 0.0; }
	}
	E = NULL;
}

simpleGraph::~simpleGraph() {
	simpleEdge *curr, *prev;
	for (int i=0; i<n; i++) {
		curr = nodeLink[i];
		delete [] A[i];
		while (curr != NULL) {
			prev = curr;
			curr = curr->next;
			delete prev;
		}
	}
	curr = NULL; prev = NULL;
	if (E != NULL) { delete [] E;	E = NULL; }
	delete [] A;			A			= NULL;
	delete [] nodeLink;		nodeLink		= NULL;
	delete [] nodeLinkTail;  nodeLinkTail   = NULL;
	delete [] nodes;		nodes		= NULL;
}

// ********************************************************************************************************

bool simpleGraph::addGroup(const int i, const int group_index) {
	if (i >= 0 and i < n) {
		nodes[i].group_true = group_index;
		return true;
	} else { return false; }
}

// ********************************************************************************************************

bool simpleGraph::addLink(const int i, const int j) {
	// Adds the directed edge (i,j) to the adjacency list for v_i
	simpleEdge* newedge;
	if (i >= 0 and i < n and j >= 0 and j < n) {
		A[i][j] = 1.0;
		newedge	 = new simpleEdge;
		newedge->x = j;
		if (nodeLink[i] == NULL) {			// first neighbor
			nodeLink[i]	 = newedge;
			nodeLinkTail[i] = newedge;
			nodes[i].degree = 1;
		} else {							// subsequent neighbor
			nodeLinkTail[i]->next = newedge;
			nodeLinkTail[i]       = newedge;
			nodes[i].degree++;
		}
		m++;								// increment edge count
		newedge = NULL;
		return true;
	} else { return false; }
}

// ********************************************************************************************************

bool simpleGraph::doesLinkExist(const int i, const int j) {
	// This function determines if the edge (i,j) already exists in the adjacency list of v_i
	if (i >= 0 and i < n and j >= 0 and j < n) { if (A[i][j] > 0.1) { return true; } else { return false; } } else { return false; }
	return false;
}

// ********************************************************************************************************

double	  simpleGraph::getAdjacency(const int i, const int j) {
	if (i >= 0 and i < n and j >= 0 and j < n) { return A[i][j]; } else { return -1.0; } }
int		  simpleGraph::getDegree(const int i)       { if (i >= 0 and i < n) { return nodes[i].degree;     } else { return -1;   } }
int		  simpleGraph::getGroupLabel(const int i)   { if (i >= 0 and i < n) { return nodes[i].group_true; } else { return -1;   } }
string	  simpleGraph::getName(const int i)         { if (i >= 0 and i < n) { return nodes[i].name;       } else { return "";   } }
// NOTE: The following three functions return addresses; deallocation of returned object is dangerous
simpleEdge* simpleGraph::getNeighborList(const int i) { if (i >= 0 and i < n) { return nodeLink[i];     } else { return NULL; } }
// END-NOTE

// ********************************************************************************************************

int		  simpleGraph::getNumGroups()       { return num_groups; }
int		  simpleGraph::getNumLinks()        { return m; }
int		  simpleGraph::getNumNodes()        { return n; }
simpleVert* simpleGraph::getNode(const int i) { if (i >= 0 and i<n) { return &nodes[i]; } else { return NULL; } }

// ********************************************************************************************************

void simpleGraph::printAdjacencies() {
	for (int i=0; i<n; i++) {
		cout << "[" << i << "]";
		for (int j=0; j<n; j++) { cout << " " << A[i][j]; }
		cout << endl;
	}
	return;
}

// ********************************************************************************************************

void simpleGraph::printPairs() {
	simpleEdge* curr;
	int edgeCount = 0;
	for (int i=0; i<n; i++) {
		cout << "[" << i << "]\t";
		curr = nodeLink[i];
		while (curr != NULL) {
			cout << curr->x << "\t";
			edgeCount++;
			curr = curr->next;
		}
		cout << "\n";
	}
	curr = NULL;
	cout << edgeCount << " edges total.\n";
	return;
}

// ********************************************************************************************************

bool simpleGraph::setName(const int i, const string text) {
	if (i >= 0 and i < n) { nodes[i].name = text; return true; } else { return false; } }

// ********************************************************************************************************

void simpleGraph::QsortMain (block* array, int left, int right) {
	if (right > left) {
		int pivot = left;
		int part  = QsortPartition(array, left, right, pivot);
		QsortMain(array, left,   part-1);
		QsortMain(array, part+1, right  );
	}
	return;
}

int simpleGraph::QsortPartition (block* array, int left, int right, int index) {
	block p_value, temp;
	p_value.x = array[index].x;		p_value.y = array[index].y;
	
	// swap(array[p_value], array[right])
	temp.x		= array[right].x;   temp.y		= array[right].y;
	array[right].x = array[index].x;   array[right].y = array[index].y;
	array[index].x = temp.x;			array[index].y = temp.y;
	
	int stored       = left;
	for (int i=left; i<right; i++) {
		if (array[i].x <= p_value.x) {
			// swap(array[stored], array[i])
			temp.x          = array[i].x;		temp.y          = array[i].y;
			array[i].x      = array[stored].x; array[i].y      = array[stored].y;
			array[stored].x = temp.x;		array[stored].y = temp.y;
			stored++;
		}
	}
	// swap(array[right], array[stored])
	temp.x          = array[stored].x; temp.y          = array[stored].y;
	array[stored].x = array[right].x;  array[stored].y = array[right].y;
	array[right].x  = temp.x;		array[right].y  = temp.y;
	
	return stored;
}

// ********************************************************************************************************
// ********************************************************************************************************

} // namespace fitHRG

#endif
