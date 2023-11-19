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

#ifndef IGRAPH_HRG_RBTREE
#define IGRAPH_HRG_RBTREE

namespace fitHRG {

// ******** Basic Structures *********************************************

struct list {
    int x = -1;             // stored elementd in linked-list
    list* next = nullptr;   // pointer to next elementd
};

struct keyValuePair {
    int x = -1;                     // elementrb key (int)
    int y = -1;                     // stored value (int)
    keyValuePair* next = nullptr;   // linked-list pointer
};

// ******** Tree elementrb Class *****************************************

struct elementrb {
    int key = -1;               // search key (int)
    int value = -1;             // stored value (int)

    bool color = false;         // F: BLACK, T: RED
    short int mark = 0;         // marker

    elementrb *parent = nullptr;    // pointer to parent node
    elementrb *left = nullptr;      // pointer for left subtree
    elementrb *right = nullptr;     // pointer for right subtree
};

// ******** Red-Black Tree Class *****************************************
// This vector implementation is a red-black balanced binary tree data
// structure. It provides find a stored elementrb in time O(log n),
// find the maximum elementrb in time O(1), delete an elementrb in
// time O(log n), and insert an elementrb in time O(log n).
//
// Note that the key=0 is assumed to be a special value, and thus you
// cannot insert such an item. Beware of this limitation.

class rbtree {
    elementrb* root;      // binary tree root
    elementrb* leaf;      // all leaf nodes
    int support;          // number of nodes in the tree

    void rotateLeft(elementrb *x);    // left-rotation operator
    void rotateRight(elementrb *y);   // right-rotation operator
    void insertCleanup(elementrb *z); // house-keeping after insertion
    void deleteCleanup(elementrb *x); // house-keeping after deletion
    keyValuePair* returnSubtreeAsList(const elementrb *z, keyValuePair *head) const;
    void deleteSubTree(elementrb *z); // delete subtree rooted at z
    elementrb* returnMinKey(elementrb *z) const; // returns minimum of subtree
    // rooted at z
    elementrb* returnSuccessor(elementrb *z) const; // returns successor of z's key

public:
    rbtree(); ~rbtree(); // default constructor/destructor

    // returns value associated with searchKey
    int returnValue(int searchKey) const;
    // returns T if searchKey found, and points foundNode at the
    // corresponding node
    elementrb* findItem(int searchKey) const;
    // insert a new key with stored value
    void insertItem(int newKey, int newValue);
    // delete a node with given key
    void deleteItem(int killKey);
    // replace value of a node with given key
    void replaceItem(int key, int newValue);
    // increment the value of the given key
    void incrementValue(int key);
    // delete the entire tree
    void deleteTree();
    // return array of keys in tree
    int* returnArrayOfKeys() const;
    // return list of keys in tree
    list* returnListOfKeys() const;
    // return the tree as a list of keyValuePairs
    keyValuePair* returnTreeAsList() const;
    // returns the maximum key in the tree
    keyValuePair returnMaxKey() const;
    // returns the minimum key in the tree
    keyValuePair returnMinKey() const;
    // returns number of items in tree
    int returnNodecount() const;
};

}
#endif
