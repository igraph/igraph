// ****************************************************************************************************
// *** COPYRIGHT NOTICE *******************************************************************************
// rbtree - red-black tree (self-balancing binary tree data structure)
// Copyright (C) 2004 Aaron Clauset
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
// Collaborators: Cristopher Moore and Mark Newman
// Project      : Hierarchical Random Graphs
// Location     : University of New Mexico, Dept. of Computer Science AND Santa Fe Institute
// Created      : Spring 2004
// Modified     : many, many times
//
// ****************************************************************************************************

#if !defined(rbtree_INCLUDED)
#define rbtree_INCLUDED

#include <iostream>

using namespace std;

namespace fitHRG {

// ******** Basic Structures ******************************************************************************

#if !defined(list_INCLUDED)
#define list_INCLUDED
class list {
public:
	int		x;				// stored elementd in linked-list
	list*	next;			// pointer to next elementd
	list()  { x = -1; next = NULL; }
	~list() {}
};
#endif

class keyValuePair {
public:
	int		x;					// elementrb key (int)
	int		y;					// stored value (int)
	keyValuePair*	next;			// linked-list pointer
	keyValuePair()  { x = -1; y = -1; next = NULL; }
	~keyValuePair() {}
};

// ******** Tree elementrb Class ****************************************************************************

class elementrb {
public:
	int		key;					// search key (int)
	int		value;				// stored value (int)
	
	bool		color;				// F: BLACK
								// T: RED
	short int mark;				// marker
	
	elementrb   *parent;			// pointer to parent node
	elementrb   *left;				// pointer for left subtree
	elementrb   *right;				// pointer for right subtree
	
	elementrb()  {	key = -1; value = -1; color = false; mark = 0;
						parent  = NULL; left  = NULL; right  = NULL; }
	~elementrb() {}
};

// ******** Red-Black Tree Class **************************************************************************
/*   This vector implementation is a red-black balanced binary tree data structure.
 *   It provides find a stored elementrb in time O(log n), find the maximum elementrb in time O(1),
 *	delete an elementrb in time O(log n), and insert an elementrb in time O(log n).
 *
 *	Note that the key=0 is assumed to be a special value, and thus you cannot insert such an item. 
 *	Beware of this limitation.
 */

class rbtree {
private:
	elementrb*		root;						// binary tree root
	elementrb*		leaf;						// all leaf nodes
	int				support;						// number of nodes in the tree

	void				rotateLeft(elementrb *x);		// left-rotation operator
	void				rotateRight(elementrb *y);		// right-rotation operator
	void				insertCleanup(elementrb *z);		// house-keeping after insertion
	void				deleteCleanup(elementrb *x);		// house-keeping after deletion
	keyValuePair*		returnSubtreeAsList(elementrb *z, keyValuePair *head);
	void				printSubTree(elementrb *z);		// display the subtree rooted at z
	void				deleteSubTree(elementrb *z);		// delete subtree rooted at z
	elementrb*		returnMinKey(elementrb *z);		// returns minimum of subtree rooted at z
	elementrb*		returnSuccessor(elementrb *z);	// returns successor of z's key
	
public:
	rbtree(); ~rbtree();							// default constructor/destructor

	int			returnValue(const int searchKey);		// returns value associated with searchKey
	elementrb*	findItem(const int searchKey);		// returns T if searchKey found, and
												// points foundNode at the corresponding node
	void			insertItem(int newKey, int newValue);	// insert a new key with stored value
	void			deleteItem(int killKey);				// selete a node with given key
 	void			replaceItem(int key, int newValue);	// replace value of a node with given key
 	void			incrementValue(int key);			// increment the value of the given key
	void			deleteTree();						// delete the entire tree
	int*			returnArrayOfKeys();				// return array of keys in tree
	list*		returnListOfKeys();					// return list of keys in tree
	keyValuePair*	returnTreeAsList();					// return the tree as a list of keyValuePairs
	keyValuePair	returnMaxKey();					// returns the maximum key in the tree
	keyValuePair	returnMinKey();					// returns the minimum key in the tree
	int			returnNodecount();					// returns number of items in tree

	void			printTree();						// displays tree (in-order traversal)

};

}
#endif
