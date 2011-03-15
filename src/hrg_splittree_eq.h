// ****************************************************************************************************
// *** COPYRIGHT NOTICE *******************************************************************************
// splittree_eq.h - a binary search tree data structure for storing dendrogram split frequencies
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
//			 : 20 May 2008 (cleaned up for public consumption)
//
// ****************************************************************************************************
// 
// Data structure for storing the split frequences in the sampled dendrograms. Data is stored 
// efficiently as a red-black binary search tree (this is a modified version of the rbtree.h file).
//
// ****************************************************************************************************

#if !defined(splittree_INCLUDED)
#define splittree_INCLUDED

#include <iostream>

using namespace std;

namespace fitHRG {

// ******** Basic Structures ******************************************************************************

#if !defined(slist_INCLUDED)
#define slist_INCLUDED
class slist {
public:
	string	x;				// stored elementd in linked-list
	slist*	next;			// pointer to next elementd
	slist()  { x = ""; next = NULL; }
	~slist() {}
};
#endif

class keyValuePairSplit {
public:
	string	x;					// elementsp split (string)
	double	y;					// stored weight   (double)
	int		c;					// stored count    (int)
	keyValuePairSplit*	next;		// linked-list pointer
	keyValuePairSplit()  { x = ""; y = 0.0; c = 0; next = NULL; }
	~keyValuePairSplit() {}
};

// ******** Tree elementsp Class ****************************************************************************

class elementsp {
public:
	string	split;				// split represented as a string
	double	weight;				// total weight of this split
	int		count;				// number of observations of this split
	
	bool		color;				// F: BLACK
								// T: RED
	short int mark;				// marker

	elementsp   *parent;			// pointer to parent node
	elementsp   *left;				// pointer for left subtree
	elementsp   *right;				// pointer for right subtree
	
	elementsp()  {	split = ""; weight = 0.0; count = 0; color = false; mark = 0;
						parent  = NULL; left  = NULL; right  = NULL; }
	~elementsp() {}
};

// ******** Red-Black Tree Class **************************************************************************
/*   This vector implementation is a red-black balanced binary tree data structure.
 *   It provides find a stored elementsp in time O(log n), find the maximum elementsp in time O(1),
 *	delete an elementsp in time O(log n), and insert an elementsp in time O(log n).
 *
 *	Note that the split="" is assumed to be a special value, and thus you cannot insert such an item. 
 *	Beware of this limitation.
 */

class splittree {
private:
	elementsp*		root;						// binary tree root
	elementsp*		leaf;						// all leaf nodes
	int				support;						// number of nodes in the tree
	double			total_weight;					// total weight stored
	int				total_count;					// total number of observations stored

	void				rotateLeft(elementsp*);			// left-rotation operator
	void				rotateRight(elementsp*);			// right-rotation operator
	void				insertCleanup(elementsp*);		// house-keeping after insertion
	void				deleteCleanup(elementsp*);		// house-keeping after deletion
	keyValuePairSplit*	returnSubtreeAsList(elementsp*, keyValuePairSplit*);
	void				printSubTree(elementsp*);		// display the subtree rooted at z
	void				deleteSubTree(elementsp*);		// delete subtree rooted at z
	elementsp*		returnMinKey(elementsp*);		// returns minimum of subtree rooted at z
	elementsp*		returnSuccessor(elementsp*);		// returns successor of z's key
	
public:
	splittree(); ~splittree();						// default constructor/destructor

	double			returnValue(const string);		// returns value associated with searchKey
	elementsp*		findItem(const string);			// returns T if searchKey found, and
												// points foundNode at the corresponding node
	void				finishedThisRound();			// update total_count and total_weight
	bool				insertItem(string, double);		// insert a new key with stored value
	void				clearTree();
	void				deleteItem(string);				// selete a node with given key
	void				deleteTree();					// delete the entire tree
	string*			returnArrayOfKeys();			// return array of keys in tree
	slist*			returnListOfKeys();				// return list of keys in tree
	keyValuePairSplit*	returnTreeAsList();				// return the tree as a list of keyValuePairSplits
	keyValuePairSplit	returnMaxKey();				// returns the maximum key in the tree
	keyValuePairSplit	returnMinKey();				// returns the minimum key in the tree
	int				returnNodecount();				// returns number of items in tree
	keyValuePairSplit*	returnTheseSplits(const int);		// returns list of splits with given number of Ms
	double			returnTotal();					// returns sum of stored values

	void				printTree();					// displays tree (in-order traversal)
	void				printTreeAsList();				// list keys (in-order) with values
	void				printTreeAsShortList();			// short list keys (in-order) with values,
	void				recordTreeAsList(const string, const double);	// write list of keys and values to file
};

// ********************************************************************************************************
// ********************************************************************************************************

}

#endif
