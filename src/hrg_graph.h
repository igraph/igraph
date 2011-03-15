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

#if !defined(graph_INCLUDED)
#define graph_INCLUDED

#include <stdio.h>
#include <string>
#include "stdlib.h"

#include "hrg_rbtree.h"

using namespace std;

namespace fitHRG {

// ******** Basic Structures ******************************************************************************

#if !defined(edge_INCLUDED)
#define edge_INCLUDED
class edge {
public:
	int		x;				// stored integer value  (edge terminator)
	edge*	next;			// pointer to next elementd
	edge()  { x = -1; next = NULL; }
	~edge() { }
};
#endif

#if !defined(vert_INCLUDED)
#define vert_INCLUDED
class vert {
public:
	string	name;			// (external) name of vertex
	int		degree;			// degree of this vertex

	vert()  { name = ""; degree = 0; }
	~vert() { }
};
#endif

// ******** Graph Class with Edge Statistics *************************************************************

class graph {
public:
	graph(const int);
	~graph();

	bool		addLink(const int, const int);					// add (i,j) to graph
	bool		doesLinkExist(const int, const int);				// true if (i,j) is already in graph
	int		getDegree(const int);							// returns degree of vertex i
	string	getName(const int);								// returns name of vertex i
	edge*	getNeighborList(const int);						// returns edge list of vertex i
	int		numLinks();									// returns m
	int		numNodes();									// returns n
	void		printPairs();									// prints all edges in graph
	bool		setName(const int, const string);					// set name of vertex i

private:
	vert*	nodes;			// list of nodes
	edge**	nodeLink;			// linked list of neighbors to vertex
	edge**	nodeLinkTail;		// pointers to tail of neighbor list
	int		n;				// number of vertices
	int		m;				// number of directed edges
};

}

#endif
