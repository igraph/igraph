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

#include "hrg_rbtree.h"
#include "hrg_dendro.h"
#include "hrg_graph.h"
#include "hrg_splittree_eq.h"

#include "igraph_hrg.h"
#include "igraph_constructors.h"

using namespace fitHRG;

// ******** Red-Black Tree Methods ************************************************************************

rbtree::rbtree() {
	root = new elementrb;
	leaf = new elementrb;

	leaf->parent   = root;

	root->left	= leaf;
	root->right    = leaf;
	support		= 0;
}

rbtree::~rbtree() {
	if (root != NULL && (root->left != leaf || root->right != leaf)) { deleteSubTree(root); }
	support   = 0;
	delete leaf;
	root		= NULL;
	leaf		= NULL;
}

void rbtree::deleteTree() { if (root != NULL) { deleteSubTree(root); } return; } // does not leak memory

void rbtree::deleteSubTree(elementrb *z) {

	if (z->left  != leaf) { deleteSubTree(z->left);  }
	if (z->right != leaf) { deleteSubTree(z->right); }
	delete z;
	z = NULL;
	return;
}

// ******** Search Functions ******************************************************************************
// public search function - if there exists a elementrb in the tree with key=searchKey,
// it returns TRUE and foundNode is set to point to the found node; otherwise, it sets
// foundNode=NULL and returns FALSE
elementrb* rbtree::findItem(const int searchKey) {

	elementrb *current;    current = root;
	if (current->key==-1) { return NULL; }							// empty tree; bail out
	while (current != leaf) {
		if (searchKey < current->key) {							// left-or-right?
			if (current->left  != leaf) { current = current->left;  }	// try moving down-left
			else { return NULL; }								//   failure; bail out
		} else {												// 
			if (searchKey > current->key) {							// left-or-right?
				if (current->right  != leaf) { current = current->right;  }	// try moving down-left
				else { return NULL; }							//   failure; bail out
			} else { return current; }							// found (searchKey==current->key)
		}
	}
	return NULL;
} // does not leak memory

int rbtree::returnValue(const int searchKey) {
	elementrb* test = findItem(searchKey);
	if (test == NULL) { return 0; } else { return test->value; }
}


// ******** Return Item Functions *************************************************************************

int* rbtree::returnArrayOfKeys() {
	int* array;
	array = new int [support];
	bool flag_go = true;
	int index = 0;
	elementrb *curr;

	if (support == 1) { array[0] = root->key; }
	else if (support == 2) {
		array[0] = root->key;
		if (root->left == leaf) { array[1] = root->right->key; } 
		else { array[1] = root->left->key; }
	} else {
		for (int i=0; i<support; i++) { array[i] = -1; }
		// non-recursive traversal of tree structure
		curr		 = root;
		curr->mark = 1;
		while (flag_go) {
			
			if (curr->mark == 1 and curr->left == leaf) {		// - is it time, and is left child the leaf node?
				curr->mark = 2;							// 
			}
			if (curr->mark == 2 and curr->right == leaf) {		// - is it time, and is right child the leaf node?
				curr->mark = 3;							// 
			}
			if (curr->mark == 1) {							// - go left
				curr->mark = 2;							// 
				curr       = curr->left;						// 
				curr->mark = 1;							// 
			} else if (curr->mark == 2) {						// - else go right
				curr->mark = 3;							// 
				curr       = curr->right;					// 
				curr->mark = 1;							// 
			} else {										// - else go up a level
				curr->mark = 0;							// 
				array[index++] = curr->key;					// 
				curr = curr->parent;						// 
				if (curr == NULL) { flag_go = false; }			// 
			}
		}
	}
	
	return array;
} // does not leak memory

list* rbtree::returnListOfKeys() {
	keyValuePair *curr, *prev;
	list         *head=0, *tail=0, *newlist;

	curr = returnTreeAsList();
	while (curr != NULL) {
		newlist    = new list;
		newlist->x = curr->x;
		if (head == NULL) { head       = newlist; tail = head;    }
		else              { tail->next = newlist; tail = newlist; }
		prev = curr;
		curr = curr->next;
		delete prev;
		prev = NULL;
	}
	return head;
}

keyValuePair* rbtree::returnTreeAsList() { // pre-order traversal
	keyValuePair  *head, *tail;

	head    = new keyValuePair;
	head->x = root->key;
	head->y = root->value;
	tail = head;

	if (root->left  != leaf) { tail = returnSubtreeAsList(root->left,  tail); }
	if (root->right != leaf) { tail = returnSubtreeAsList(root->right, tail); }
	
	if (head->x == -1) { return NULL; /* empty tree */ } else { return head; }
}

keyValuePair* rbtree::returnSubtreeAsList(elementrb *z, keyValuePair *head) {
	keyValuePair *newnode, *tail;
	
	newnode    = new keyValuePair;
	newnode->x = z->key;
	newnode->y = z->value;
	head->next = newnode;
	tail       = newnode;
	
	if (z->left  != leaf) { tail = returnSubtreeAsList(z->left,  tail); }
	if (z->right != leaf) { tail = returnSubtreeAsList(z->right, tail); }
	
	return tail;
}

keyValuePair rbtree::returnMaxKey() {
	keyValuePair themax;
	elementrb *current;
	current  = root;
	while (current->right != leaf) {		// search to bottom-right corner of tree
		current  = current->right; }		// 
	themax.x = current->key;				// store the data found
	themax.y = current->value;			// 
	
	return themax;						// return that data
}

keyValuePair rbtree::returnMinKey() {
	keyValuePair themin;
	elementrb *current;
	current = root;
	while (current->left != leaf) {		// search to bottom-left corner of tree
		current = current->left; }		// 
	themin.x = current->key;				// store the data found
	themin.y = current->value;			// 
	
	return themin;						// return that data
}

// private functions for deleteItem() (although these could easily be made public, I suppose)
elementrb* rbtree::returnMinKey(elementrb *z) {
	elementrb *current;

	current = z;
	while (current->left != leaf) {		// search to bottom-right corner of tree
		current = current->left; }		// 
	return current;					// return pointer to the minimum
}

elementrb* rbtree::returnSuccessor(elementrb *z) {
	elementrb *current, *w;
	
	w = z;
	if (w->right != leaf) {				// if right-subtree exists, return min of it
		return returnMinKey(w->right); }
	current = w->parent;				// else search up in tree
	while ((current!=NULL) && (w==current->right)) {
		w       = current;
		current = current->parent;		// move up in tree until find a non-right-child
	}
	return current;
}

int rbtree::returnNodecount() { return support; }

// ******** Insert Functions ******************************************************************************
// public insert function
void rbtree::insertItem(int newKey, int newValue) {
	
	// first we check to see if newKey is already present in the tree; if so, we do nothing;
	// if not, we must find where to insert the key
	elementrb *newNode, *current;

	current = findItem(newKey);						// find newKey in tree; return pointer to it O(log k)
	if (current == NULL) {
		newNode			= new elementrb;				// elementrb for the rbtree
		newNode->key		= newKey;					//  store newKey
		newNode->value		= newValue;  				//  store newValue
		newNode->color		= true;					//  new nodes are always RED
		newNode->parent	= NULL;					//  new node initially has no parent
		newNode->left		= leaf;					//  left leaf
		newNode->right		= leaf;					//  right leaf
		support++;								// increment node count in rbtree
		
		// must now search for where to insert newNode, i.e., find the correct parent and
		// set the parent and child to point to each other properly
		current = root;
		if (current->key==-1) {										// insert as root
			delete root;											//   delete old root
			root			= newNode;								//   set root to newNode
			leaf->parent   = newNode;								//   set leaf's parent
			current		= leaf;									//   skip next loop
		}
		
		while (current != leaf) {									// search for insertion point
			if (newKey < current->key) {								// left-or-right?
				if (current->left  != leaf) { current = current->left;  }	// try moving down-left
				else {											// else found new parent
					newNode->parent	= current;					//    set parent
					current->left		= newNode;					//    set child
					current			= leaf;						//    exit search
				}
			} else {												// 
				if (current->right != leaf) { current = current->right; }   // try moving down-right
				else {											// else found new parent
					newNode->parent	= current;					//    set parent
					current->right		= newNode;					//    set child
					current			= leaf;						//    exit search
				}
			}
		}

		// now do the house-keeping necessary to preserve the red-black properties
		insertCleanup(newNode);			// do house-keeping to maintain balance
	}
	return;
}

// private house-keeping function for insertion
void rbtree::insertCleanup(elementrb *z) {
	
	if (z->parent==NULL) {								// fix now if z is root
		z->color = false; return; }
	elementrb *temp;
	while (z->parent!=NULL && z->parent->color) {	// while z is not root and z's parent is RED
		if (z->parent == z->parent->parent->left) {  // z's parent is LEFT-CHILD
			temp = z->parent->parent->right;		// grab z's uncle
			if (temp->color) {
				z->parent->color		= false;  // color z's parent BLACK	(Case 1)
				temp->color			= false;  // color z's uncle BLACK		(Case 1)
				z->parent->parent->color = true;   // color z's grandparent RED  (Case 1)
				z = z->parent->parent;			// set z = z's grandparent    (Case 1)
			} else {
				if (z == z->parent->right) {		// z is RIGHT-CHILD
					z = z->parent;				// set z = z's parent		(Case 2)
					rotateLeft(z);				// perform left-rotation		(Case 2)
				}
				z->parent->color		= false;  // color z's parent BLACK	(Case 3)
				z->parent->parent->color = true;   // color z's grandparent RED  (Case 3)
				rotateRight(z->parent->parent);    // perform right-rotation	(Case 3)
			}
		} else {								// z's parent is RIGHT-CHILD
			temp = z->parent->parent->left;		// grab z's uncle
			if (temp->color) {
				z->parent->color		= false;  // color z's parent BLACK	(Case 1)
				temp->color			= false;  // color z's uncle BLACK		(Case 1)
				z->parent->parent->color = true;   // color z's grandparent RED  (Case 1)
				z = z->parent->parent;			// set z = z's grandparent    (Case 1)
			} else {
				if (z == z->parent->left) {		// z is LEFT-CHILD
					z = z->parent;				// set z = z's parent		(Case 2)
					rotateRight(z);			// perform right-rotation	(Case 2)
				}
				z->parent->color		= false;  // color z's parent BLACK	(Case 3)
				z->parent->parent->color = true;   // color z's grandparent RED  (Case 3)
				rotateLeft(z->parent->parent);	// perform left-rotation		(Case 3)
			}
		}
	}

	root->color = false;						// color the root BLACK
	return;
}

// ******** Delete Functions ******************************************************************************
void rbtree::replaceItem(int key, int newValue) {
  elementrb* ptr;
  ptr = findItem(key);
  ptr->value = newValue;
  return;
}
 
void rbtree::incrementValue(int key) {
  elementrb* ptr;
  ptr = findItem(key);
  ptr->value = 1+ptr->value;
  return;
}

// public delete function
void rbtree::deleteItem(int killKey) {
	elementrb *x, *y, *z;
	
	z = findItem(killKey);
	if (z == NULL) { return; }						// item not present; bail out

	if (support==1) {								// -- attempt to delete the root
		root->key		= -1;						// restore root node to default state
		root->value    = -1;						// 
		root->color    = false;						// 
		root->parent   = NULL;						// 
		root->left	= leaf;						// 
		root->right    = leaf;						// 
		support--;								// set support to zero
		return;									// exit - no more work to do
	}
	
	if (z != NULL) {
		support--;								// decrement node count
		if ((z->left == leaf) || (z->right==leaf)) {		// case of less than two children
			  y = z; }							//    set y to be z
		else { y = returnSuccessor(z); }				//    set y to be z's key-successor
		
		if (y->left!=leaf) { x = y->left; }			// pick y's one child (left-child)
		else			    { x = y->right; }			//				  (right-child)
		x->parent = y->parent;						// make y's child's parent be y's parent

		if (y->parent==NULL) { root = x; }				// if y is the root, x is now root
		else {									// 
			if (y == y->parent->left) {				// decide y's relationship with y's parent
				y->parent->left  = x;				//   replace x as y's parent's left child
			} else {								// 
				y->parent->right = x; }				//   replace x as y's parent's left child
		}										// 

		if (y!=z) {								// insert y into z's spot
			z->key		= y->key;					// copy y data into z
			z->value		= y->value;				// 
		}										// 

		if (y->color==false) { deleteCleanup(x); }		// do house-keeping to maintain balance
		delete y;									// deallocate y
		y = NULL;									// point y to NULL for safety
	}											// 
		
	return;
} // does not leak memory

void rbtree::deleteCleanup(elementrb *x) {
	elementrb *w, *t;
	while ((x != root) && (x->color==false)) {			// until x is the root, or x is RED
		if (x==x->parent->left) {					// branch on x being a LEFT-CHILD
			w = x->parent->right;					// grab x's sibling
			if (w->color==true) {					// if x's sibling is RED
				w->color = false;					// color w BLACK				(case 1)
				x->parent->color = true;				// color x's parent RED			(case 1)
				rotateLeft(x->parent);				// left rotation on x's parent	(case 1)
				w = x->parent->right;				// make w be x's right sibling	(case 1)
			}
			if ((w->left->color==false) && (w->right->color==false)) {
				w->color = true;					// color w RED					(case 2)
				x = x->parent;						// examine x's parent			(case 2)
			} else {								// 
				if (w->right->color==false) {			// 
					w->left->color = false;			// color w's left child BLACK		(case 3)
					w->color = true;				// color w RED					(case 3)
					t = x->parent;					// store x's parent
					rotateRight(w);				// right rotation on w			(case 3)
					x->parent = t;					// restore x's parent
					w = x->parent->right;			// make w be x's right sibling	(case 3)
				}								// 
				w->color			= x->parent->color; // make w's color = x's parent's   (case 4)
				x->parent->color    = false;			// color x's parent BLACK		(case 4)
				w->right->color	= false;			// color w's right child BLACK	(case 4)
				rotateLeft(x->parent);				// left rotation on x's parent	(case 4)
				x = root;							// finished work. bail out		(case 4)
			}									// 
		} else {									// x is RIGHT-CHILD
			w = x->parent->left;					// grab x's sibling
			if (w->color==true) {					// if x's sibling is RED
				w->color			= false;			// color w BLACK				(case 1)
				x->parent->color    = true;			// color x's parent RED			(case 1)
				rotateRight(x->parent);				// right rotation on x's parent	(case 1)
				w				= x->parent->left;  // make w be x's left sibling		(case 1)
			}
			if ((w->right->color==false) && (w->left->color==false)) {
				w->color = true;					// color w RED					(case 2)
				x= x->parent;						// examine x's parent			(case 2)
			} else {								// 
				if (w->left->color==false) {			// 
					w->right->color	= false;		// color w's right child BLACK	(case 3)
					w->color			= true;		// color w RED					(case 3)
					t				= x->parent;   // store x's parent
					rotateLeft(w);					// left rotation on w			(case 3)
					x->parent			= t;			// restore x's parent
					w = x->parent->left;			// make w be x's left sibling		(case 3)
				}								// 
				w->color = x->parent->color;			// make w's color = x's parent's   (case 4)
				x->parent->color    = false;			// color x's parent BLACK		(case 4)
				w->left->color		= false;			// color w's left child BLACK		(case 4)
				rotateRight(x->parent);				// right rotation on x's parent    (case 4)
				x				= root;			// x is now the root			(case 4)
			}
		}
	}
	x->color = false;								// color x (the root) BLACK		(exit)

	return;
}

// ******** Rotation Functions ****************************************************************************

void rbtree::rotateLeft(elementrb *x) {
	elementrb *y;
	// do pointer-swapping operations for left-rotation
	y               = x->right;					// grab right child
	x->right        = y->left;					// make x's RIGHT-CHILD be y's LEFT-CHILD
	y->left->parent = x;						// make x be y's LEFT-CHILD's parent
	y->parent       = x->parent;					// make y's new parent be x's old parent

	if (x->parent==NULL) { root = y; }				// if x was root, make y root
	else {									// 
		if (x == x->parent->left)				// if x is LEFT-CHILD, make y be x's parent's
			{ x->parent->left  = y; }			//    left-child
		else { x->parent->right = y; }			//    right-child
	}										// 
	y->left   = x;								// make x be y's LEFT-CHILD
	x->parent = y;								// make y be x's parent
	
	return;
}

void rbtree::rotateRight(elementrb *y) {
	elementrb *x;
	// do pointer-swapping operations for right-rotation
	x                = y->left;					// grab left child
	y->left          = x->right;					// replace left child yith x's right subtree
	x->right->parent = y;						// replace y as x's right subtree's parent
	
	x->parent        = y->parent;					// make x's new parent be y's old parent
	if (y->parent==NULL) { root = x; }				// if y was root, make x root
	else {
		if (y == y->parent->right)				// if y is RIGHT-CHILD, make x be y's parent's
			{ y->parent->right  = x; }			//    right-child
		else { y->parent->left   = x; }			//    left-child
	}
	x->right  = y;								// make y be x's RIGHT-CHILD
	y->parent = x;								// make x be y's parent
	
	return;
}

// ******** Display Functions *****************************************************************************
// public
void rbtree::printTree() {
	cout << "\nTREE SIZE = " << support << endl;
	cout << "# "; printSubTree(root);
	return;
}

//private
void rbtree::printSubTree(elementrb *z) {
	if (z==leaf) { return; }
	else {
		cout << "(" << z->key << " " << z->value << " " << z->color << ")"<<endl;
		cout << "L "; printSubTree(z->left); cout << endl;
		cout << "R "; printSubTree(z->right); cout << endl;
	}
	return;
}

// ********************************************************************************************************
// ********************************************************************************************************

// ****************************************************************************************************
// *** COPYRIGHT NOTICE *******************************************************************************
// dendro.h - hierarchical random graph (hrg) data structure
// Copyright (C) 2005-2009 Aaron Clauset
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
// Created      : 26 October 2005 - 7 December 2005
// Modified     : 23 December 2007 (cleaned up for public consumption)
//
// ****************************************************************************************************
// 
// Maximum likelihood dendrogram data structure. This is the heart of the HRG algorithm: all
// manipulations are done here and all data is stored here. The data structure uses the separate
// graph data structure to store the basic adjacency information (in a dangerously mutable way).
// 
// ****************************************************************************************************

// ******** Dendrogram Methods ****************************************************************************

dendro::dendro() { root   = NULL; internal   = NULL;
				leaf   = NULL; d          = NULL;
				paths  = NULL;
				g      = NULL; splithist  = NULL;
				ctree  = NULL; cancestor  = NULL;
}
dendro::~dendro() {
	list *curr, *prev;
	if (g	    != NULL) { delete g;           g		= NULL; }    // O(m)
	if (internal  != NULL) { delete [] internal; internal  = NULL; }    // O(n)
	if (leaf      != NULL) { delete [] leaf;     leaf	     = NULL; }    // O(n)
	if (d         != NULL) { delete d;           d         = NULL; }    // O(n)
	if (splithist != NULL) { delete splithist;   splithist = NULL; }    // potentially a long time
	if (paths     != NULL) { for (int i=0; i<n; i++) { curr = paths[i]; while (curr != NULL) { prev = curr;   curr = curr->next;   delete prev;   prev = NULL; } paths[i] = NULL; } delete [] paths; } paths = NULL;
	if (ctree     != NULL) { delete ctree;       ctree     = NULL; }    // O(n)
	if (cancestor != NULL) { delete cancestor;   cancestor = NULL; }    // O(n)
}

// ********************************************************************************************************

void dendro::binarySearchInsert(elementd* x, elementd* y) {
	if (y->p < x->p) {		// go to left subtree
		if (x->L == NULL) { // check if left subtree is empty
			x->L = y;		// make x left child
			y->M = x;		// make y parent of child
			return;
		}
		else { binarySearchInsert(x->L, y); }
	} else {				// go to right subtree
		if (x->R == NULL) { // check if right subtree is empty
			x->R = y;		// make x right child
			y->M = x;		// make y parent of child
			return;
		}
		else { binarySearchInsert(x->R, y); }
	}
	return;
}

// ********************************************************************************************************

list* dendro::binarySearchFind(const double v) {
	list *head, *tail, *newlist;
	elementd *current = root;
	bool flag_stopSearch = false;
	
	while (!flag_stopSearch) {				// continue until we're finished
		newlist    = new list;				// add this node to the path
		newlist->x = current->label;
		if (current == root) { head       = newlist; tail = head;    }
		else                 { tail->next = newlist; tail = newlist; }
		if (v < current->p) {				// now try left subtree
			if (current->L->type == GRAPH) { flag_stopSearch = true;       }
			else						 { current         = current->L; }
		} else {							// else try right subtree
			if (current->R->type == GRAPH) { flag_stopSearch = true;       }
			else						 { current         = current->R; }
		}
	}
	return head;
}

// ********************************************************************************************************

string dendro::buildSplit(elementd* thisNode) {
	// A "split" is defined as the bipartition of vertices into the sets of leaves below the
	// internal vertex in the tree (denoted by "C"), and those above it (denoted as "M"). For
	// simplicity, we represent this bipartition as a character string of length n, where the
	// ith character denotes the partition membership (C,M) of the ith leaf node.

	bool      flag_go = true;
	const short int k = 1+DENDRO+GRAPH;
	elementd* curr;;
	split sp;

	sp.initializeSplit(n);						// default split string O(n)
	
	curr = thisNode;							// - set start node as top this sub-tree
	curr->type = k+1;							// - initialize in-order tree traversal
	while (flag_go) {						
		if (curr->type == k+1 and
		    curr->L->type == GRAPH) {		// - is it time, and is left child a graph node?
			sp.s[curr->L->index] = 'C';	// - mark this leaf
			curr->type           = k+2;	// 
		}
		if (curr->type == k+2 and		
		    curr->R->type == GRAPH) {		// - is it time, and is right child a graph node?
			sp.s[curr->R->index] = 'C';	// - mark this leaf
			curr->type           = k+3;	//
		}
		if (curr->type == k+1) {			// - go left
			curr->type = k+2;			//
			curr       = curr->L;		//
			curr->type = k+1; }
		else if (curr->type == k+2) {		// - else go right
			curr->type = k+3;			// 
			curr       = curr->R;		// 
			curr->type = k+1;	
		} else {						// - else go up a level
			curr->type = DENDRO;		// 
			if (curr->index == thisNode->index || curr->M == NULL) { flag_go = false; curr = NULL; }
			else { curr = curr->M; }		//
		}
	}

	// any leaf that was not already marked must be in the remainder of the tree
	for (int i=0; i<n; i++) {  if (sp.s[i] != 'C') { sp.s[i] = 'M'; } }
	if (!sp.checkSplit()) { cout << "buildSplit:: malformed split at [ " << thisNode->index <<" ]: " << sp.s << endl; }

	return sp.s;
}

// ********************************************************************************************************

void dendro::buildDendrogram() {
	if (g == NULL) { cout << "Error: cannot build dendrogram without a graph structure.\n"; return; }
	
/* the initialization of the dendrogram structure goes like this:
 * 1) we allocate space for the n-1 internal nodes of the dendrogram, and then the n leaf nodes
 * 2) we build a random binary tree structure out of the internal nodes by assigning each
 *    a uniformly random value over [0,1] and then inserting it into the tree according to the 
 *    binary-search rule.
 * 3) next, we make a random permutation of the n leaf nodes and add them to the dendrogram D by
 *    replacing the emptpy spots in-order
 * 4) then, we compute the path from the root to each leaf and store that in each leaf (this is
 *    prep work for the next step)
 * 5) finally, we compute the values for nL, nR, e (and thus p) and the label for each internal 
 *    node by allocating each of the m edges in g to the appropriate internal node
 */
	
	// --- Initialization and memory allocation for data structures
	// After allocating the memory for D and G, we need to mark the nodes for G as being
	// non-internal vertices, and then insert them into a random binary tree structure.
	// For simplicity, we make the first internal node in the array the root.
	bool flag_debug = false;
	n		= g->numNodes();		// size of graph
	leaf		= new elementd [n];		// allocate memory for G, O(n)
	internal  = new elementd [n-1];	// allocate memory for D, O(n)
	d		= new interns(n-2);      // allocate memory for internal edges of D, O(n)
	for (int i=0; i<n; i++) {		// initialize leaf nodes
		leaf[i].type   = GRAPH;
		leaf[i].label  = i;
		leaf[i].index  = i;
		leaf[i].n		= 1;
	}
	if (flag_debug) { cout << ">> dendro: allocated memory for internal and leaf arrays" << endl; }
	root		  = &internal[0];		// initialize internal nodes
	root->label = 0;
	root->index = 0;
	root->p     = mtr.randExc();
	for (int i=1; i<(n-1); i++) {		// insert remaining internal vertices, O(n log n)
		internal[i].label = i;
		internal[i].index = i;
		internal[i].p = mtr.randExc();
		binarySearchInsert(root, &internal[i]);
	}
	if (flag_debug) { cout << ">> dendro: inserted internal vertices into random binary tree" << endl; }
	
	// --- Hang leaf nodes off end of dendrogram O(n log n)
	// To impose this random hierarchical relationship on G, we first take a random permutation
	// of the leaf vertices and then replace the NULLs at the bottom of the tree in-order with
	// the leafs. As a hack to ensure that we can find the leafs later using a binary search,
	// we assign each of them the p value of their parent, perturbed slightly so as to preserve
	// the binary search property.
	block* array; array = new block [n];
	for (int i=0; i<n; i++) { array[i].x = mtr.randExc();  array[i].y = i; }
	QsortMain(array, 0, n-1);

	int k=0;						// replace NULLs with leaf nodes, and
	for (int i=0; i<(n-1); i++) {		//    maintain binary search property, O(n)
		if (internal[i].L == NULL) {
			internal[i].L = &leaf[array[k].y];
			leaf[array[k].y].M = &internal[i];
			leaf[array[k++].y].p = internal[i].p - 0.0000000000001;
		}
		if (internal[i].R == NULL) {
			internal[i].R = &leaf[array[k].y];
			leaf[array[k].y].M = &internal[i];
			leaf[array[k++].y].p = internal[i].p + 0.0000000000001;
		}
	}
	delete [] array;
	if (flag_debug) { cout << ">> dendro: replaced NULLs in bin-tree with leaf nodes" << endl; }

	// --- Compute the path from root -> leaf for each leaf O(n log n)
	// Using the binary search property, we can find each leaf node in O(log n) time. The
	// binarySearchFind() function returns the list of internal node indices that the search
	// crossed, in the order of root -> ... -> leaf, for use in the subsequent few operations.
	if (paths != NULL) { list *curr, *prev; for (int i=0; i<n; i++) { curr = paths[i]; while (curr != NULL) { prev = curr;   curr = curr->next;   delete prev;   prev = NULL; } paths[i] = NULL; } delete [] paths; } paths = NULL;
	paths = new list* [n];
	for (int i=0; i<n; i++) { paths[i] = binarySearchFind(leaf[i].p); }
	
	if (flag_debug) { cout << ">> dendro: computed paths from root to leafs" << endl; }
	
	// --- Count e for each internal node O(m)
	// To count the number of edges that span the L and R subtrees for each internal node,
	// we use the path information we just computed. Then, we loop over all edges in G
	// and find the common ancestor in D of the two endpoints and increment that internal
	// node's e count. This process takes O(m) time because in a roughly balanced binary
	// tree (given by our random dendrogram), the vast majority of vertices take basically
	// constant time to find their common ancestor. Note that because our adjacency list
	// is symmetric, we overcount each e by a factor of 2, so we need to correct this after.
	elementd* ancestor; edge* curr;
	for (int i=0; i<(n-1); i++) { internal[i].e = 0; internal[i].label = -1; }
	for (int i=0; i<n; i++) {
		curr = g->getNeighborList(i);
		while (curr != NULL) {
			ancestor     = findCommonAncestor(paths, i, curr->x);
			ancestor->e += 1;
			curr		   = curr->next;
		}
	}
	for (int i=0; i<(n-1); i++) {
		internal[i].e /= 2; }   
	if (flag_debug) { cout << ">> dendro: finished common ancestor computation" << endl; }

	// --- Count n for each internal node O(n log n)
	// To tabulate the number of leafs in each subtree rooted at an internal node,
	// we use the path information computed above.
	for (int i=0; i<n; i++) {
		ancestor = &leaf[i];
		ancestor = ancestor->M;
		while (ancestor != NULL) {
			ancestor->n++;
			ancestor = ancestor->M;
		}
	}
	if (flag_debug) { cout << ">> dendro: computed subtree sizes" << endl; }
	
	// --- Label all internal vertices O(n log n)
	// We want to label each internal vertex with the smallest leaf index of its children.
	// This will allow us to collapse many leaf-orderings into a single dendrogram structure
	// that is independent of child-exhanges (since these have no impact on the likelihood
	// of the hierarchical structure). To do this, we loop over the leaf vertices from 
	// smallest to largest and walk along that leaf's path from the root. If we find an 
	// unlabeled internal node, then we mark it with this leaf's index.

	for (int i=0; i<n; i++) {
		ancestor = &leaf[i];
		while (ancestor != NULL) {
			if (ancestor->label == -1 || ancestor->label > leaf[i].label) { ancestor->label = leaf[i].label; }
			ancestor = ancestor->M;
		}
	}
	if (flag_debug) { cout << ">> dendro: labeled all internal vertices" << endl; }
	
	// --- Exchange children to enforce order-property O(n)
	// We state that the order-property requires that an internal node's label is the 
	// smallest index of its left subtree. The dendrogram so far doesn't reflect this, so we
	// need to step through each internal vertex and make that adjustment (swapping nL and nR
	// if we make a change).
	elementd *tempe;
	for (int i=0; i<(n-1); i++) {
		if (internal[i].L->label > internal[i].label) {
			tempe          = internal[i].L;
			internal[i].L  = internal[i].R;
			internal[i].R  = tempe;
		}
	}
	if (flag_debug) { cout << ">> dendro: enforced order-property" << endl; }

	// --- Tabulate internal dendrogram edges O(n^2)
	// For the MCMC moves later on, we'll need to be able to choose, uniformly at random, an
	// internal edge of the dendrogram to manipulate. There are always n-2 of them, and we can
	// find them simply by scanning across the internal vertices and observing which have children
	// that are also internal vertices. Note: very important that the order property be enforced
	// before this step is taken; otherwise, the internal edges wont reflect the actual dendrogram
	// structure.
	
	for (int i=0; i<(n-1); i++) {
		if (internal[i].L->type == DENDRO) { d->addEdge     (i, internal[i].L->index, LEFT    ); }
		// cout << "(" << i << " " << internal[i].L->index << " L) "; }
		if (internal[i].R->type == DENDRO) { d->addEdge     (i, internal[i].R->index, RIGHT    ); }
		// cout << "(" << i << " " << internal[i].R->index << " R) "; }
	}	
	if (flag_debug) { cout << ">> dendro: tabulated internal dendrogram edges" << endl; }
	
	// --- Clear memory for paths O(n log n)
	// Now that we're finished using the paths, we need to deallocate them manually.
	list *current, *previous;
	for (int i=0; i<n; i++) {
		current = paths[i];
		while (current != NULL) { previous = current;   current = current->next;   delete previous;   previous = NULL; }
		paths[i] = NULL;
	}
	delete [] paths;
	paths = NULL;
	if (flag_debug) { cout << ">> dendro: cleared memory for paths" << endl; }
	
	// --- Compute p_i for each internal node O(n)
	// Each internal node's p_i = e_i / (nL_i*nR_i), and now that we have each of those
	// pieces, we may calculate this value for each internal node. Given these, we can then
	// calculate the log-likelihood of the entire dendrogram structure
	// \log(L) = \sum_{i=1}^{n} ( ( e_i \log[p_i] ) + ( (nL_i*nR_i - e_i) \log[1-p_i] ) )
	L = 0.0; double dL;
	int nL_nR, ei;
	for (int i=0; i<(n-1); i++) {
		nL_nR = internal[i].L->n*internal[i].R->n;
		ei    = internal[i].e;
		internal[i].p = (double)(ei) / (double)(nL_nR);
		if (ei == 0 or ei == nL_nR) { dL = 0.0; }
		else                        { dL = ei * log(internal[i].p) + (nL_nR - ei) * log(1.0-internal[i].p); }
//		cout << "p[" << i << "] = " << internal[i].p << "\tdL = " << dL << endl; // << "\t (1) " << ei * log(internal[i].p) << "\t (2) " << (nL_nR - ei) * log(1.0-internal[i].p)
		internal[i].logL = dL;
		L += dL;
	}
	
	if (flag_debug) {
		cout << ">> dendro: computed internal node probability value" << endl;
		if (n<100) { printDendrogram(); }
		cout << "Log-Likelihood = " << L << endl;
	}
	char pauseme;
	for (int i=0; i<(n-1); i++) {
		if (internal[i].label > internal[i].L->label) {
			tempe          = internal[i].L;
			internal[i].L  = internal[i].R;
			internal[i].R  = tempe;
			cout << "#### WARNING - order property violated by internal[" << i << "] (fixed)" << endl;
			cin >> pauseme;
		}
	}
		
	// --- Dendrogram is now built
//	d->printEdgeList();
	if (flag_debug) { cout << ">> dendro: build dendrogram complete" << endl; }

	return;
}

// ********************************************************************************************************

void dendro::clearDendrograph() {
	// Clear out the memory and references used by the dendrograph structure - this is 
	// intended to be called just before an importDendrogramStructure call so as to avoid
	// memory leaks and overwriting the references therein.
	if (g        != NULL) { delete    g;        g        = NULL; }    // O(m)
	if (leaf     != NULL) { delete [] leaf;     leaf     = NULL; }    // O(n)
	if (internal != NULL) { delete [] internal; internal = NULL; }    // O(n)
	if (d        != NULL) { delete    d;	    d        = NULL; }    // O(n)
	root = NULL;
	
	return;
}

// ********************************************************************************************************

int dendro::computeEdgeCount(const int a, const short int atype, const int b, const short int btype) {
	// This function computes the number of edges that cross between the subtree internal[a]
	// and the subtree internal[b].
	// To do this, we use an array A[1..n] integers which take values -1 if A[i] is in the
	// subtree defined by internal[a], +1 if A[i] is in the subtree internal[b], and 0
	// otherwise. Taking the smaller of the two sets, we then scan over the edges attached
	// to that set of vertices and count the number of endpoints we see in the other set.
	
	bool flag_go    = true;
	int nA, nB;
	int         count = 0;
	const short int k = 1+DENDRO+GRAPH;

	elementd* curr;
	
	// --- First, we push the leaf nodes in the L and R subtrees into balanced binary tree
	//     structures so that we can search them quickly later on.
	if (atype == GRAPH) {					// default case, subtree A is size 1
		subtreeL.insertItem(a,-1);			// insert single node as member of left subtree
		nA		 = 1;					// 
	} else {
		curr		 = &internal[a];			// explore subtree A, O(|A|)
		curr->type = k+1;					//
		nA         = 0;					//
		while (flag_go) {
			if (curr->index == internal[a].M->index) {
				internal[a].type = DENDRO;
				flag_go          = false;
			} else {
				if (curr->type == k+1 and 
				    curr->L->type == GRAPH) {		// - is it time, and is left child a graph node?
					subtreeL.insertItem(curr->L->index, -1);
					curr->type            = k+2;  //
					nA++;					//
				}
				if (curr->type == k+2 and		
				    curr->R->type == GRAPH) {		// - is it time, and is right child a graph node?
					subtreeL.insertItem(curr->R->index, -1);
					curr->type            = k+3;  //
					nA++;					//
				}
				if (curr->type == k+1) {			// - go left
					curr->type = k+2;			//
					curr       = curr->L;		//
					curr->type = k+1; }
				else if (curr->type == k+2) {		// - else go right
					curr->type = k+3;			// 
					curr       = curr->R;		// 
					curr->type = k+1;	
				} else {						// - else go up a level
					curr->type = DENDRO;		// 
					curr       = curr->M;		// 
					if (curr == NULL) {
						flag_go = false;
					}
				}
			}
			if (nA > n) { cout << "error! nA > n\n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n" << endl; break; }
		}
	}
	
	if (btype == GRAPH) {					// default case, subtree A is size 1
		subtreeR.insertItem(b,1);			// insert node as single member of right subtree
		nB		 = 1;					// 
	} else {
		flag_go = true;
		curr       = &internal[b];			// explore subtree B, O(|B|)
		curr->type = k+1;					// 
		nB		 = 0;					//
		while (flag_go) {
			if (curr->index == internal[b].M->index) {
				internal[b].type = DENDRO;
				flag_go = false;
			} else {
				if (curr->type == k+1 and 
				    curr->L->type == GRAPH) {		// - is it time, and is left child a graph node?
					subtreeR.insertItem(curr->L->index, 1);
					curr->type            = k+2;  //
					nB++;					// 
				}
				if (curr->type == k+2 and 
				    curr->R->type == GRAPH) {		// - is it time, and is right child a graph node?
					subtreeR.insertItem(curr->R->index, 1);
					curr->type            = k+3;  // 
					nB++;					//
				}
				if (curr->type == k+1) {			// - look left
					curr->type = k+2;			// 
					curr       = curr->L;		// 
					curr->type = k+1; }
				else if (curr->type == k+2) {		// - look right
					curr->type = k+3;			// 
					curr       = curr->R;		// 
					curr->type = k+1;
				} else {						// - else go up a level
					curr->type = DENDRO;		// 
					curr       = curr->M;		// 
					if (curr == NULL) {
						flag_go = false;
					}
				}
			}
			if (nB > n) { cout << "error! nB > n \n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n" << endl; break; }
		}
	}

	// --- Now, we take the smaller subtree and ask how many of its emerging edges have their
	//     partner in the other subtree. O(|A| log |A|) time
	edge* current;
	int*  treeList;
	if (nA < nB) {
		treeList = subtreeL.returnArrayOfKeys();	// subtreeL is smaller
		for (int i=0; i<nA; i++) {
			current = g->getNeighborList(treeList[i]);
			while (current != NULL) {			// loop over each of its neighbors v_j
				if (subtreeR.findItem(current->x) != NULL) { count++; }
				current = current->next;			// to see if v_j is in A
			}								// 
			subtreeL.deleteItem(treeList[i]);
		}
		delete [] treeList;
		treeList = subtreeR.returnArrayOfKeys();
		for (int i=0; i<nB; i++) { subtreeR.deleteItem(treeList[i]); }
		delete [] treeList;
	} else {
		treeList = subtreeR.returnArrayOfKeys();	// subtreeR is smaller
		for (int i=0; i<nB; i++) {
			current = g->getNeighborList(treeList[i]);
			while (current != NULL) {			// loop over each of its neighbors v_j
				if (subtreeL.findItem(current->x) != NULL) { count++; }
				current = current->next;			// to see if v_j is in B
			}								// 
			subtreeR.deleteItem(treeList[i]);
		}
		delete [] treeList;
		treeList = subtreeL.returnArrayOfKeys();
		for (int i=0; i<nA; i++) { subtreeL.deleteItem(treeList[i]); }
		delete [] treeList;
	}

//	cout << "finished edge counting; count = " << count << endl;

	return count;
}

// ********************************************************************************************************

int dendro::countChildren(const string s) {
	char pauseme;
	int len = s.size();
	if (len != n) { cout << "something is wrong: length(" << s << ") = " << len << " != " << n << endl; cin >> pauseme; }
	int numC = 0;
	for (int i=0; i<len; i++) {
		if (s[i] == 'C') { numC++; } else
			if (s[i] != 'M') { cout << "something is wrong: s[" << i << "] = " << s << endl; cin >> pauseme; }
	}
	return numC;
}

// ********************************************************************************************************

void dendro::cullSplitHist() {
	string* array;
	int tot, leng;
	
	array = splithist->returnArrayOfKeys();
	tot   = splithist->returnTotal();
	leng  = splithist->returnNodecount();
	for (int i=0; i<leng; i++) {
		if ((splithist->returnValue(array[i]) / tot) < 0.5) { splithist->deleteItem(array[i]); }
	}
	delete [] array; array = NULL;
	
	return;
}

// ********************************************************************************************************

elementd* dendro::findCommonAncestor(list** paths, const int i, const int j) {
	list* headOne = paths[i];
	list* headTwo = paths[j];
	elementd* lastStep;
	while (headOne->x == headTwo->x) {
		lastStep = &internal[headOne->x];
		headOne  = headOne->next;
		headTwo  = headTwo->next;
		if (headOne == NULL or headTwo == NULL) { break; }
	}
	return lastStep;			// Returns address of an internal node; do not deallocate
}

// ********************************************************************************************************

int dendro::getConsensusSize() {
	string    *array;
	double     value, tot;
	int		 numSplits, numCons;
	numSplits = splithist->returnNodecount();
	array     = splithist->returnArrayOfKeys();
	tot       = splithist->returnTotal();
	numCons	= 0;
	for (int i=0; i<numSplits; i++) {
		value = splithist->returnValue(array[i]);
		if (value / tot > 0.5) { numCons++; }
	}
	delete [] array; array = NULL;
	return numCons;
}

// ********************************************************************************************************

splittree* dendro::getConsensusSplits() {
	string    *array;
	splittree *consensusTree;
	double     value, tot;
	consensusTree  = new splittree;
	int numSplits;
	
	// We look at all of the splits in our split histogram and add any one that's in the
	// majority to our consensusTree, which we then return (note that consensusTree needs to
	// be deallocated by the user).
	numSplits = splithist->returnNodecount();
	array     = splithist->returnArrayOfKeys();
	tot       = splithist->returnTotal();
	for (int i=0; i<numSplits; i++) {
		value = splithist->returnValue(array[i]);
		if (value / tot > 0.5) { consensusTree->insertItem(array[i], value / tot); }
	}
	delete [] array; array = NULL;
	return consensusTree;
}

// ********************************************************************************************************

double dendro::getLikelihood() { return L; }

// ********************************************************************************************************

void dendro::getSplitList(splittree* split_tree) {
	string sp;
	for (int i=0; i<(n-1); i++) {
		sp = d->getSplit(i);
		if (sp != "" and sp[1] != '-') { split_tree->insertItem(sp,0.0); }
	}
	return;
}

// ********************************************************************************************************

double dendro::getSplitTotalWeight() { 
	if (splithist != NULL) { return splithist->returnTotal(); } else { return 0; }
}

// ********************************************************************************************************

bool dendro::importDendrogramStructure(const string in_file) {
	string bracketL, bracketR, sL, sR, sLtype, sRtype, sp, se, sn;
	int sindex, sLindex, sRindex, snume, snumn;
	double sprob;
	bool safeExit   = true;
	bool flag_debug = false;
	n = 1;
	
	ifstream fscan(in_file.c_str(), ios::in);
	while (fscan >> sn) { if (sn == "[") { n++; } }
	fscan.close();
	
	leaf		= new elementd [n];		// allocate memory for G, O(n)
	internal  = new elementd [n-1];	// allocate memory for D, O(n)
	d		= new interns(n-2);		// allocate memory for internal edges of D, O(n)
	for (int i=0; i<n; i++) {		// initialize leaf nodes
		leaf[i].type   = GRAPH;
		leaf[i].label  = i;
		leaf[i].index  = i;
		leaf[i].n		= 1;
	}
	root		  = &internal[0];		// initialize internal nodes
	root->label = 0;
	for (int i=1; i<(n-1); i++) { internal[i].index = i; internal[i].label = -1; }
	if (flag_debug) { cout << ">> dendro: allocated memory for internal and leaf arrays" << endl; }
	
	// --- Import basic structure from file O(n)
	ifstream fin(in_file.c_str(), ios::in);
	while (fin >> bracketL >> sindex >> bracketR >> sL >> sLindex >> sLtype >> sR >> sRindex >> sRtype >> sp >> sprob >> se >> snume >> sn >> snumn) {
		cout << bracketL << " " << sindex << " " << bracketR << " " << sL << " " << sLindex << " " << sLtype << " " << sR << " " << sRindex << " " << sRtype << " " << sp << " " << sprob << " " << se << " " << snume << " " << sn << " " << snumn << endl;
		if (sLtype == "(D)") { internal[sindex].L = &internal[sLindex]; internal[sLindex].M = &internal[sindex]; } else 
		if (sLtype == "(G)") { internal[sindex].L = &leaf[sLindex];     leaf[sLindex].M     = &internal[sindex]; } else
						 { cout << "Error: " << bracketL << sindex << bracketR << sL << sLindex << sLtype << sR << sRindex << sRtype << sp << sprob << se << snume << sn << snumn << endl; safeExit = false; break; }
		if (sRtype == "(D)") { internal[sindex].R = &internal[sRindex]; internal[sRindex].M = &internal[sindex]; } else
		if (sRtype == "(G)") { internal[sindex].R = &leaf[sRindex];     leaf[sRindex].M     = &internal[sindex]; } else 
						 { cout << "Error: " << bracketL << sindex << bracketR << sL << sLindex << sLtype << sR << sRindex << sRtype << sp << sprob << se << snume << sn << snumn << endl; safeExit = false; break; }
		internal[sindex].p     = sprob; if (sprob < 0.0 || sprob > 1.0) { cout << "Error: " << bracketL << sindex << bracketR << sL << sLindex << sLtype << sR << sRindex << sRtype << sp << sprob << se << snume << sn << snumn << endl; safeExit = false; break; }
		internal[sindex].e     = snume;
		internal[sindex].n     = snumn;
		internal[sindex].index = sindex;
	}
	fin.close();
	if (!safeExit) { return false; }
	if (flag_debug) { cout << ">> dendro: imported basic structure" << endl; }

	// --- Label all internal vertices O(n log n)
	elementd* curr;
	for (int i=0; i<n; i++) {
		curr = &leaf[i];
		while (curr != NULL) { if (curr->label == -1 || curr->label > leaf[i].label) { curr->label = leaf[i].label; } curr = curr->M; }
	}
	if (flag_debug) { cout << ">> dendro: labeled all internal vertices" << endl; }
	
	// --- Exchange children to enforce order-property O(n)
	elementd *tempe;
	for (int i=0; i<(n-1); i++) {
		if (internal[i].L->label > internal[i].label) {
			tempe          = internal[i].L;
			internal[i].L  = internal[i].R;
			internal[i].R  = tempe;
		}
	}
	if (flag_debug) { cout << ">> dendro: enforced order-property" << endl; }
	
	// --- Tabulate internal dendrogram edges O(n)
	for (int i=0; i<(n-1); i++) {
		if (internal[i].L->type == DENDRO) { d->addEdge(i, internal[i].L->index, LEFT);  }
		if (internal[i].R->type == DENDRO) { d->addEdge(i, internal[i].R->index, RIGHT); }
	}
	if (flag_debug) { cout << ">> dendro: tabulated internal dendrogram edges" << endl; }
	
	// --- Compute p_i for each internal node O(n)
	// Each internal node's p_i = e_i / (nL_i*nR_i), and now that we have each of those
	// pieces, we may calculate this value for each internal node. Given these, we can then
	// calculate the log-likelihood of the entire dendrogram structure
	// \log(L) = \sum_{i=1}^{n} ( ( e_i \log[p_i] ) + ( (nL_i*nR_i - e_i) \log[1-p_i] ) )
	L = 0.0; double dL;
	int nL_nR, ei;
	for (int i=0; i<(n-1); i++) {
		nL_nR = internal[i].L->n*internal[i].R->n;
		ei    = internal[i].e;
		if (ei == 0 or ei == nL_nR) { dL = 0.0; }
		else                        { dL = (double)(ei) * log(internal[i].p) + (double)(nL_nR - ei) * log(1.0-internal[i].p); }
		internal[i].logL = dL;
		L += dL;
	}
	if (flag_debug) {
		cout << ">> dendro: computed log-likelihood" << endl;
		cout << "Log-Likelihood = " << L << endl;
	}
	
	// --- Dendrogram is now built
	if (flag_debug) { cout << ">> dendro: build dendrogram complete" << endl; }
	
	return true;
}

bool dendro::importDendrogramStructure(const igraph_hrg_t *hrg) {
  n=igraph_hrg_size(hrg);
  
  // allocate memory for G, O(n)
  leaf = new elementd[n];
  // allocate memory for D, O(n)
  internal = new elementd[n-1];
  // allocate memory for internal edges of D, O(n)
  d = new interns(n-2);

  // initialize leaf nodes
  for (int i=0; i<n; i++) {
    leaf[i].type  = GRAPH;
    leaf[i].label = i;
    leaf[i].index = i;
    leaf[i].n     = 1;
  }

  // initialize internal nodes
  root = &internal[0];
  root->label=0;
  for (int i=1; i<n-1; i++) {
    internal[i].index = i;
    internal[i].label = -1;
  }

  // import basic structure from hrg object, O(n)  
  for (int i=0; i<n-1; i++) {
    int L=VECTOR(hrg->left)[i];
    int R=VECTOR(hrg->right)[i];
    
    if (L < 0) { 
      internal[i].L = &internal[-L-1]; 
      internal[-L-1].M = &internal[i];
    } else {
      internal[i].L = &leaf[L];
      leaf[L].M = &internal[i];
    }
    
    if (R < 0) { 
      internal[i].R = &internal[-R-1];
      internal[-R-1].M = &internal[i];
    } else {
      internal[i].R = &leaf[R];
      leaf[R].M = &internal[i];
    }
    
    internal[i].p = VECTOR(hrg->prob)[i];
    internal[i].e = VECTOR(hrg->edges)[i];
    internal[i].n = VECTOR(hrg->vertices)[i];
    internal[i].index = i;
  }
  
  // --- Label all internal vertices O(n log n)
  elementd *curr;
  for (int i=0; i<n; i++) {
    curr=&leaf[i];
    while (curr) { 
      if (curr->label == -1 || curr->label > leaf[i].label) {
	curr->label = leaf[i].label; 
      } 
      curr = curr -> M;
    }
  }
  
  // --- Exchange children to enforce order-property O(n)
  elementd *tempe;
  for (int i=0; i<n-1; i++) {
    if (internal[i].L->label > internal[i].label) {
      tempe          = internal[i].L;
      internal[i].L  = internal[i].R;
      internal[i].R  = tempe;
    }
  }
  
  // --- Tabulate internal dendrogram edges O(n)
  for (int i=0; i<(n-1); i++) {
    if (internal[i].L->type == DENDRO) { d->addEdge(i, internal[i].L->index, LEFT);  }
    if (internal[i].R->type == DENDRO) { d->addEdge(i, internal[i].R->index, RIGHT); }
  }
  
  // --- Compute p_i for each internal node O(n)
  // Each internal node's p_i = e_i / (nL_i*nR_i), and now that we
  // have each of those pieces, we may calculate this value for each
  // internal node. Given these, we can then calculate the
  // log-likelihood of the entire dendrogram structure 
  // \log(L) = \sum_{i=1}^{n} ( ( e_i \log[p_i] ) + 
  // ( (nL_i*nR_i - e_i) \log[1-p_i] ) )
  L = 0.0; double dL;
  int nL_nR, ei;
  for (int i=0; i<(n-1); i++) {
    nL_nR = internal[i].L->n*internal[i].R->n;
    ei    = internal[i].e;
    if (ei == 0 or ei == nL_nR) { 
      dL = 0.0; 
    } else { 
      dL = (double)(ei) * log(internal[i].p) + (double)(nL_nR - ei) * log(1.0-internal[i].p); 
    }
    internal[i].logL = dL;
    L += dL;
  }
  
  return true;
}

// ********************************************************************************************************

void dendro::makeRandomGraph() {
	bool flag_debug = false;
	if (flag_debug) { cout << ">> dendro: making random graph from dendrogram" << endl; }
	if (g != NULL) { delete g; } g = NULL; g = new graph(n);

	list	*curr, *prev; 
	if (paths != NULL) { for (int i=0; i<n; i++) { curr = paths[i]; while (curr != NULL) {
		prev = curr;   curr = curr->next;   delete prev;   prev = NULL; } paths[i] = NULL; } delete [] paths; } paths = NULL;
	paths = new list* [n];				// build paths from root O(n d)
	for (int i=0; i<n; i++) { paths[i] = reversePathToRoot(i); }
	
	elementd* commonAncestor;
	for (int i=0; i<n; i++) {			// O((h+d)*n^2) - h: height of D; d: average degree in G
		for (int j=(i+1); j<n; j++) {		// decide neighbors of v_i
			commonAncestor = findCommonAncestor(paths,i,j);
			if (mtr.randExc() < commonAncestor->p) { 
				if (!(g->doesLinkExist(i,j))) { if (!(g->addLink(i,j))) { cout << "Error: (" << j << " " << i << ")" << endl; } }
				if (!(g->doesLinkExist(j,i))) { if (!(g->addLink(j,i))) { cout << "Error: (" << j << " " << i << ")" << endl; } }
			}
		}
	}
	// g->printPairs();
	
	for (int i=0; i<n; i++) { curr = paths[i]; while (curr != NULL) { prev = curr;   curr = curr->next;   delete prev;   prev = NULL; } paths[i] = NULL; }
	delete [] paths;					// delete paths data structure O(n log n)
	paths = NULL;

	return;
}

// ********************************************************************************************************

bool dendro::monteCarloMove(double& delta, bool& ftaken, const double T) {
	// A single MC move begins with the selection of a random internal edge (a,b) of the
	// dendrogram. This also determines the three subtrees i, j, k that we will rearrange,
	// and we choose uniformly from among the options.
	// 
	// If (a,b) is a left-edge, then we have ((i,j),k), and moves
	// ((i,j),k) -> ((i,k),j)								(alpha move)
	//           -> (i,(j,k)) + enforce order-property for (j,k)	(beta move)
	// 
	// If (a,b) is a right-edge, then we have (i,(j,k)), and moves
	// (i,(j,k)) -> ((i,k),j)								(alpha move)
	//           -> ((i,j),k)								(beta move)
	// 
	// For each of these moves, we need to know what the change in likelihood will be, so
	// that we can determine with what probability we execute the move.
	
	bool		flag_debug = false;
	elementd	*temp;
	ipair	*tempPair;
	int		x, y, e_x, e_y, n_i, n_j, n_k, n_x, n_y;
	short int t;
	double    p_x, p_y, L_x, L_y, dLogL;
	string    new_split;
	
	// The remainder of the code executes a single MCMC move, where we sample the dendrograms 
	// proportionally to their likelihoods (i.e., temperature=1, if you're comparing it to the
	// usual MCMC framework). 
	delta    = 0.0;
	ftaken   = false;
	tempPair = d->getRandomEdge();	// returns address; no need to deallocate
	x        = tempPair->x;			// copy contents of referenced random edge
	y        = tempPair->y;			//    into local variables
	t        = tempPair->t;
	
	if (flag_debug) {
		if (t == LEFT) {
		} else if (t == RIGHT) {
		} else { cout << " bad edge (i,j,k)" << endl; }
	}
	
	if (t == LEFT) {								// 
		if (mtr.randExc() < 0.5) {					// ## LEFT ALPHA move: ((i,j),k) -> ((i,k),j)
			// We need to calculate the change in the likelihood (dLogL) that would result from
			// this move. Most of the information needed to do this is already available,
			// the exception being e_ik, the number of edges that span the i and k subtrees.
			// I use a slow algorithm O(n) to do this, since I don't know of a better way at
			// this point. (After several attempts to find a faster method, no luck.)

			n_i = internal[y].L->n;
			n_j = internal[y].R->n;
			n_k = internal[x].R->n;
			
			n_y  = n_i*n_k;
			e_y  = computeEdgeCount(internal[y].L->index, internal[y].L->type, internal[x].R->index, internal[x].R->type);   // e_ik
			p_y  = (double)(e_y) / (double)(n_y);
			if (e_y == 0 or e_y == n_y)   { L_y = 0.0; }
			else						{ L_y = (double)(e_y) * log(p_y) + (double)(n_y - e_y) * log(1.0-p_y); }

			n_x  = (n_i+n_k)*n_j;
			e_x  = internal[x].e + internal[y].e - e_y;						// e_yj
			p_x  = (double)(e_x) / (double)(n_x);
			if (e_x == 0 or e_x == n_x)	{ L_x = 0.0; }
			else						{ L_x = (double)(e_x) * log(p_x) + (double)(n_x - e_x) * log(1.0-p_x); }
			
			dLogL = (L_x - internal[x].logL) + (L_y - internal[y].logL);
			if ((dLogL > 0.0) or (mtr.randExc() < exp(T*dLogL))) {  // make LEFT ALPHA move
				ftaken = true;
				d->swapEdges(x, internal[x].R->index, RIGHT, y, internal[y].R->index, RIGHT);
				temp             = internal[x].R;			// - swap j and k
				internal[x].R    = internal[y].R;			// 
				internal[y].R    = temp;					// 
				internal[x].R->M = &internal[x];			// - adjust parent pointers
				internal[y].R->M = &internal[y];			// 
				internal[y].n    = n_i + n_k;				// - update n for [y]
				internal[x].e    = e_x;					// - update e_i for [x] and [y]
				internal[y].e    = e_y;					// 
				internal[x].p    = p_x;					// - update p_i for [x] and [y]
				internal[y].p    = p_y;					// 
				internal[x].logL = L_x;					// - update L_i for [x] and [y]
				internal[y].logL = L_y;					//
													// - order-property maintained
				L			 += dLogL;				// - update LogL
				delta            = dLogL;				// 
				
				// TRAP: Catches violations of the ordering property
				if (internal[x].label > internal[x].L->label or internal[y].label > internal[y].L->label) {
					printDendrogram();
					if (internal[x].label > internal[x].L->label) {
						cout << "**** WARNING - order property violated by internal[" << x << "]" << endl;
						cout << "x    (p = " << internal[x].p << "\te = "   << internal[x].e << "\tnL = ";
						cout << internal[x].L->n << "\tnR = " << internal[x].R->n << "\tlabel = " << internal[x].label << ")\tinternal[" << x << "]\t(D)" << endl;
						if (internal[x].L->type==GRAPH) { cout << "x->L [" << internal[x].L->index << "]\t(G)"<< endl; } else { 
							cout << "i->L (p = " << internal[x].L->p << "\te = "   << internal[x].L->e << "\tnL = ";
							cout << internal[x].L->L->n << "\tnR = " << internal[x].L->R->n << "\tlabel = " << internal[x].L->label << ")\tinternal[" << internal[x].L->index << "]\t(D)"<< endl; }
						if (internal[x].R->type==GRAPH) { cout << "x->R [" << internal[x].R->index << "]\t(G)"<< endl; } else { 
							cout << "i->R (p = " << internal[x].R->p << "\te = "   << internal[x].R->e << "\tnL = ";
							cout << internal[x].R->L->n << "\tnR = " << internal[x].R->R->n << "\tlabel = " << internal[x].R->label << ")\tinternal[" << internal[x].R->index << "]\t(D)"<< endl; }
					}
					if (internal[y].label > internal[y].L->label) {
						cout << "**** WARNING - order property violated by internal[" << y << "]" << endl;
						cout << "y    (p = " << internal[y].p << "\te = "   << internal[y].e << "\tnL = ";
						cout << internal[y].L->n << "\tnR = " << internal[y].R->n << "\tlabel = " << internal[y].label << ")\tinternal[" << y << "]\t(D)" << endl;
						if (internal[y].L->type==GRAPH) { cout << "y->L [" << internal[y].L->index << "]\t(G)"<< endl; } else { 
							cout << "i->L (p = " << internal[y].L->p << "\te = "   << internal[y].L->e << "\tnL = ";
							cout << internal[y].L->L->n << "\tnR = " << internal[y].L->R->n << "\tlabel = " << internal[y].L->label << ")\tinternal[" << internal[y].L->index << "]\t(D)"<< endl; }
						if (internal[y].R->type==GRAPH) { cout << "y->R [" << internal[y].R->index << "]\t(G)"<< endl; } else { 
							cout << "i->R (p = " << internal[y].R->p << "\te = "   << internal[y].R->e << "\tnL = ";
							cout << internal[y].R->L->n << "\tnR = " << internal[y].R->R->n << "\tlabel = " << internal[y].R->label << ")\tinternal[" << internal[y].R->index << "]\t(D)"<< endl; }
					}
					return false;
				}
			}
		} else {									// ## LEFT BETA move:  ((i,j),k) -> (i,(j,k))
			n_i = internal[y].L->n;
			n_j = internal[y].R->n;
			n_k = internal[x].R->n;
			
			n_y  = n_j*n_k;
			e_y  = computeEdgeCount(internal[y].R->index, internal[y].R->type, internal[x].R->index, internal[x].R->type);   // e_jk
			p_y  = (double)(e_y) / (double)(n_y);
			if (e_y == 0 or e_y == n_y)   { L_y = 0.0; }
			else						{ L_y = (double)(e_y) * log(p_y) + (double)(n_y - e_y) * log(1.0-p_y); }
			
			n_x  = (n_j+n_k)*n_i;
			e_x  = internal[x].e + internal[y].e - e_y;						// e_yj
			p_x  = (double)(e_x) / (double)(n_x);
			if (e_x == 0 or e_x == n_x)	{ L_x = 0.0; }
			else						{ L_x = (double)(e_x) * log(p_x) + (double)(n_x - e_x) * log(1.0-p_x); }
			
			dLogL = (L_x - internal[x].logL) + (L_y - internal[y].logL);
			if ((dLogL > 0.0) or (mtr.randExc() < exp(T*dLogL))) {  // make LEFT BETA move
				ftaken = true;
				d->swapEdges(y, internal[y].L->index, LEFT, y, internal[y].R->index, RIGHT);
				temp			  = internal[y].L;			// - swap L and R of [y]
				internal[y].L    = internal[y].R;			// 
				internal[y].R    = temp;					// 
				d->swapEdges(x, internal[x].R->index, RIGHT, y,internal[y].R->index, RIGHT);
				temp			  = internal[x].R;			// - swap i and k
				internal[x].R    = internal[y].R;			// 
				internal[y].R    = temp;					// 
				internal[x].R->M = &internal[x];			// - adjust parent pointers
				internal[y].R->M = &internal[y];			// 
				d->swapEdges(x, internal[x].L->index, LEFT, x, internal[x].R->index, RIGHT);
				temp			  = internal[x].L;			// - swap L and R of [x]
				internal[x].L    = internal[x].R;			// 
				internal[x].R    = temp;					// 
				internal[y].n    = n_j + n_k;				// - update n
				internal[x].e    = e_x;					// - update e_i
				internal[y].e    = e_y;					// 
				internal[x].p    = p_x;					// - update p_i
				internal[y].p    = p_y;					// 
				internal[x].logL = L_x;					// - update logL_i
				internal[y].logL = L_y;					// 
				if (internal[y].R->label < internal[y].L->label) {
					d->swapEdges(y, internal[y].L->index, LEFT, y, internal[y].R->index, RIGHT);
					temp				= internal[y].L;    // - enforce order-property if necessary
					internal[y].L		= internal[y].R;    // 
					internal[y].R		= temp;			// 
				}									// 
				internal[y].label = internal[y].L->label;    // 
				L			  += dLogL;				// - update LogL
				delta             = dLogL;				// 
				
				// TRAP: Catches violations of the ordering property
				if (internal[x].label > internal[x].L->label or internal[y].label > internal[y].L->label) {
					printDendrogram();
					if (internal[x].label > internal[x].L->label) {
						cout << "**** WARNING - order property violated by internal[" << x << "]" << endl;
						cout << "x    (p = " << internal[x].p << "\te = "   << internal[x].e << "\tnL = ";
						cout << internal[x].L->n << "\tnR = " << internal[x].R->n << "\tlabel = " << internal[x].label << ")\tinternal[" << x << "]\t(D)" << endl;
						if (internal[x].L->type==GRAPH) { cout << "x->L [" << internal[x].L->index << "]\t(G)"<< endl; } else { 
							cout << "i->L (p = " << internal[x].L->p << "\te = "   << internal[x].L->e << "\tnL = ";
							cout << internal[x].L->L->n << "\tnR = " << internal[x].L->R->n << "\tlabel = " << internal[x].L->label << ")\tinternal[" << internal[x].L->index << "]\t(D)"<< endl; }
						if (internal[x].R->type==GRAPH) { cout << "x->R [" << internal[x].R->index << "]\t(G)"<< endl; } else { 
							cout << "i->R (p = " << internal[x].R->p << "\te = "   << internal[x].R->e << "\tnL = ";
							cout << internal[x].R->L->n << "\tnR = " << internal[x].R->R->n << "\tlabel = " << internal[x].R->label << ")\tinternal[" << internal[x].R->index << "]\t(D)"<< endl; }
					}
					if (internal[y].label > internal[y].L->label) {
						cout << "**** WARNING - order property violated by internal[" << y << "]" << endl;
						cout << "y    (p = " << internal[y].p << "\te = "   << internal[y].e << "\tnL = ";
						cout << internal[y].L->n << "\tnR = " << internal[y].R->n << "\tlabel = " << internal[y].label << ")\tinternal[" << y << "]\t(D)" << endl;
						if (internal[y].L->type==GRAPH) { cout << "y->L [" << internal[y].L->index << "]\t(G)"<< endl; } else { 
							cout << "i->L (p = " << internal[y].L->p << "\te = "   << internal[y].L->e << "\tnL = ";
							cout << internal[y].L->L->n << "\tnR = " << internal[y].L->R->n << "\tlabel = " << internal[y].L->label << ")\tinternal[" << internal[y].L->index << "]\t(D)"<< endl; }
						if (internal[y].R->type==GRAPH) { cout << "y->R [" << internal[y].R->index << "]\t(G)"<< endl; } else { 
							cout << "i->R (p = " << internal[y].R->p << "\te = "   << internal[y].R->e << "\tnL = ";
							cout << internal[y].R->L->n << "\tnR = " << internal[y].R->R->n << "\tlabel = " << internal[y].R->label << ")\tinternal[" << internal[y].R->index << "]\t(D)"<< endl; }
					}
					return false;
				}
			}
		}
	} else {										// right-edge: t == RIGHT
		if (mtr.randExc() < 0.5) {					// alpha move: (i,(j,k)) -> ((i,k),j)
			n_i = internal[x].L->n;
			n_j = internal[y].L->n;
			n_k = internal[y].R->n;
			
			n_y  = n_i*n_k;
			e_y  = computeEdgeCount(internal[x].L->index, internal[x].L->type, internal[y].R->index, internal[y].R->type);   // e_ik
			p_y  = (double)(e_y) / (double)(n_y);
			if (e_y == 0 or e_y == n_y)   { L_y = 0.0; }
			else						{ L_y = (double)(e_y) * log(p_y) + (double)(n_y - e_y) * log(1.0-p_y); }
			
			n_x  = (n_i+n_k)*n_j;
			e_x  = internal[x].e + internal[y].e - e_y;						// e_yj
			p_x  = (double)(e_x) / (double)(n_x);
			if (e_x == 0 or e_x == n_x)	{ L_x = 0.0; }
			else						{ L_x = (double)(e_x) * log(p_x) + (double)(n_x - e_x) * log(1.0-p_x); }

			dLogL = (L_x - internal[x].logL) + (L_y - internal[y].logL);
			if ((dLogL > 0.0) or (mtr.randExc() < exp(T*dLogL))) {  // make RIGHT ALPHA move
				ftaken = true;
				d->swapEdges(x, internal[x].L->index, LEFT, x, internal[x].R->index, RIGHT);
				temp			   = internal[x].L;			// - swap L and R of [x]
				internal[x].L     = internal[x].R;			// 
				internal[x].R     = temp;				// 
				d->swapEdges(y, internal[y].L->index, LEFT, x, internal[x].R->index, RIGHT);
				temp			   = internal[y].L;			// - swap i and j
				internal[y].L     = internal[x].R;			// 
				internal[x].R     = temp;				// 
				internal[x].R->M  = &internal[x];			// - adjust parent pointers
				internal[y].L->M  = &internal[y];			// 
				internal[y].n     = n_i + n_k;			// - update n
				internal[x].e     = e_x;					// - update e_i
				internal[y].e     = e_y;					// 
				internal[x].p     = p_x;					// - update p_i
				internal[y].p     = p_y;					// 
				internal[x].logL  = L_x;					// - update logL_i
				internal[y].logL  = L_y;					// 
				internal[y].label = internal[x].label;		// - update order property
				L			  += dLogL;				// - update LogL
				delta             = dLogL;				// 
				
				// TRAP: Catches violations of the ordering property
				if (internal[x].label > internal[x].L->label or internal[y].label > internal[y].L->label) {
					printDendrogram();
					if (internal[x].label > internal[x].L->label) {
						cout << "**** WARNING - order property violated by internal[" << x << "]" << endl;
						cout << "x    (p = " << internal[x].p << "\te = "   << internal[x].e << "\tnL = ";
						cout << internal[x].L->n << "\tnR = " << internal[x].R->n << "\tlabel = " << internal[x].label << ")\tinternal[" << x << "]\t(D)" << endl;
						if (internal[x].L->type==GRAPH) { cout << "x->L [" << internal[x].L->index << "]\t(G)"<< endl; } else {
							cout << "i->L (p = " << internal[x].L->p << "\te = "   << internal[x].L->e << "\tnL = ";
							cout << internal[x].L->L->n << "\tnR = " << internal[x].L->R->n << "\tlabel = " << internal[x].L->label << ")\tinternal[" << internal[x].L->index << "]\t(D)"<< endl; }
						if (internal[x].R->type==GRAPH) { cout << "x->R [" << internal[x].R->index << "]\t(G)"<< endl; } else {
							cout << "i->R (p = " << internal[x].R->p << "\te = "   << internal[x].R->e << "\tnL = ";
							cout << internal[x].R->L->n << "\tnR = " << internal[x].R->R->n << "\tlabel = " << internal[x].R->label << ")\tinternal[" << internal[x].R->index << "]\t(D)"<< endl; }
					}
					if (internal[y].label > internal[y].L->label) {
						cout << "**** WARNING - order property violated by internal[" << y << "]" << endl;
						cout << "y    (p = " << internal[y].p << "\te = "   << internal[y].e << "\tnL = " << internal[y].L->n << "\tnR = " << internal[y].R->n << "\tlabel = " << internal[y].label << ")\tinternal[" << y << "]\t(D)" << endl;
						if (internal[y].L->type==GRAPH) { cout << "y->L [" << internal[y].L->index << "]\t(G)"<< endl; } else {
							cout << "i->L (p = " << internal[y].L->p << "\te = "   << internal[y].L->e << "\tnL = ";
							cout << internal[y].L->L->n << "\tnR = " << internal[y].L->R->n << "\tlabel = " << internal[y].L->label << ")\tinternal[" << internal[y].L->index << "]\t(D)"<< endl; }
						if (internal[y].R->type==GRAPH) { cout << "y->R [" << internal[y].R->index << "]\t(G)"<< endl; } else { 
							cout << "i->R (p = " << internal[y].R->p << "\te = "   << internal[y].R->e << "\tnL = ";
							cout << internal[y].R->L->n << "\tnR = " << internal[y].R->R->n << "\tlabel = " << internal[y].R->label << ")\tinternal[" << internal[y].R->index << "]\t(D)"<< endl; }
					}
					return false;
				}
			}
		} else {									// beta move:  (i,(j,k)) -> ((i,j),k)
			n_i = internal[x].L->n;
			n_j = internal[y].L->n;
			n_k = internal[y].R->n;
			
			n_y  = n_i*n_j;
			e_y  = computeEdgeCount(internal[x].L->index, internal[x].L->type, internal[y].L->index, internal[y].L->type);   // e_ij
			p_y  = (double)(e_y) / (double)(n_y);
			if (e_y == 0 or e_y == n_y)   { L_y = 0.0; }
			else						{ L_y = (double)(e_y) * log(p_y) + (double)(n_y - e_y) * log(1.0-p_y); }
			
			n_x  = (n_i+n_j)*n_k;
			e_x  = internal[x].e + internal[y].e - e_y;						// e_yk
			p_x  = (double)(e_x) / (double)(n_x);
			if (e_x == 0 or e_x == n_x)	{ L_x = 0.0; }
			else						{ L_x = (double)(e_x) * log(p_x) + (double)(n_x - e_x) * log(1.0-p_x); }

			dLogL = (L_x - internal[x].logL) + (L_y - internal[y].logL);
			if ((dLogL > 0.0) or (mtr.randExc() < exp(T*dLogL))) {  // make RIGHT BETA move
				ftaken = true;
				d->swapEdges(x, internal[x].L->index, LEFT, x, internal[x].R->index, RIGHT);
				temp			   = internal[x].L;			// - swap L and R of [x]
				internal[x].L     = internal[x].R;			// 
				internal[x].R     = temp;				// 
				d->swapEdges(x, internal[x].R->index, RIGHT, y, internal[y].R->index, RIGHT);
				temp			   = internal[x].R;			// - swap i and k
				internal[x].R     = internal[y].R;			// 
				internal[y].R     = temp;				// 
				internal[x].R->M  = &internal[x];			// - adjust parent pointers
				internal[y].R->M  = &internal[y];			// 
				d->swapEdges(y, internal[y].L->index, LEFT, y, internal[y].R->index, RIGHT);
				temp			   = internal[y].L;			// - swap L and R of [y]
				internal[y].L     = internal[y].R;			// 
				internal[y].R     = temp;				// 
				internal[y].n     = n_i + n_j;			// - update n
				internal[x].e     = e_x;					// - update e_i
				internal[y].e     = e_y;					// 
				internal[x].p     = p_x;					// - update p_i
				internal[y].p     = p_y;					//
				internal[x].logL  = L_x;					// - update logL_i
				internal[y].logL  = L_y;					// 
				internal[y].label = internal[x].label;		// - order-property
				L			  += dLogL;				// - update LogL
				delta             = dLogL;				// 
				
				// TRAP: Catches violations of the ordering property
				if (internal[x].label > internal[x].L->label or internal[y].label > internal[y].L->label) {
					printDendrogram();
					if (internal[x].label > internal[x].L->label) {
						cout << "**** WARNING - order property violated by internal[" << x << "]" << endl;
						cout << "x    (p = " << internal[x].p << "\te = "   << internal[x].e << "\tnL = ";
						cout << internal[x].L->n << "\tnR = " << internal[x].R->n << "\tlabel = " << internal[x].label << ")\tinternal[" << x << "]\t(D)" << endl;
						if (internal[x].L->type==GRAPH) { cout << "x->L [" << internal[x].L->index << "]\t(G)"<< endl; } else { 
							cout << "i->L (p = " << internal[x].L->p << "\te = "   << internal[x].L->e << "\tnL = ";
							cout << internal[x].L->L->n << "\tnR = " << internal[x].L->R->n << "\tlabel = " << internal[x].L->label << ")\tinternal[" << internal[x].L->index << "]\t(D)"<< endl; }
						if (internal[x].R->type==GRAPH) { cout << "x->R [" << internal[x].R->index << "]\t(G)"<< endl; } else { 
							cout << "i->R (p = " << internal[x].R->p << "\te = "   << internal[x].R->e << "\tnL = ";
							cout << internal[x].R->L->n << "\tnR = " << internal[x].R->R->n << "\tlabel = " << internal[x].R->label << ")\tinternal[" << internal[x].R->index << "]\t(D)"<< endl; }
					}
					if (internal[y].label > internal[y].L->label) {
						cout << "**** WARNING - order property violated by internal[" << y << "]" << endl;
						cout << "y    (p = " << internal[y].p << "\te = "   << internal[y].e << "\tnL = ";
						cout << internal[y].L->n << "\tnR = " << internal[y].R->n << "\tlabel = " << internal[y].label << ")\tinternal[" << y << "]\t(D)" << endl;
						if (internal[y].L->type==GRAPH) { cout << "y->L [" << internal[y].L->index << "]\t(G)"<< endl; } else { 
							cout << "i->L (p = " << internal[y].L->p << "\te = "   << internal[y].L->e << "\tnL = ";
							cout << internal[y].L->L->n << "\tnR = " << internal[y].L->R->n << "\tlabel = " << internal[y].L->label << ")\tinternal[" << internal[y].L->index << "]\t(D)"<< endl; }
						if (internal[y].R->type==GRAPH) { cout << "y->R [" << internal[y].R->index << "]\t(G)"<< endl; } else { 
							cout << "i->R (p = " << internal[y].R->p << "\te = "   << internal[y].R->e << "\tnL = ";
							cout << internal[y].R->L->n << "\tnR = " << internal[y].R->R->n << "\tlabel = " << internal[y].R->label << ")\tinternal[" << internal[y].R->index << "]\t(D)"<< endl; }
					}
					return false;
				}
			}
		}
	}
	
	return true;
}

// ********************************************************************************************************

void dendro::printConsensusTree() {
	child *curr;
	int treesize = splithist->returnNodecount();
	if (ctree != NULL) {
		for (int i=0; i<treesize; i++) {
			cout << "tree[" << i << "].index  = " << ctree[i].index << endl;
			cout << "tree[" << i << "].weight = " << ctree[i].weight << endl;
			cout << "tree[" << i << "].degree = " << ctree[i].degree << endl;
			cout << "tree[" << i << "].parent = " << ctree[i].parent << endl;
			curr = ctree[i].children;
			while (curr != NULL) {
				cout << curr->index << "\n";
				curr = curr->next;
			}
			cout << endl;
		}
	}
	return;
}

// ********************************************************************************************************

void dendro::printConsensusTreeDense() {
	int treesize = splithist->returnNodecount();
	if (ctree != NULL) {
		cout << "------------ internal[ 0 - " << treesize-1 << " ] ------------ " << endl;
		for (int i=0; i<treesize; i++) { cout << ctree[i].parent << "  "; } cout << endl;
	}
	if (cancestor != NULL) {
		cout << "------------ ancestor[ 1 - " << n << " ] ------------ " << endl;
		for (int i=0; i<n; i++) { cout << cancestor[i] << "  "; } cout << endl;
	}
	return;
}

// ********************************************************************************************************

void dendro::printDendrogram()      { cout << "\nLEAFS = " << n << endl << "# "; printSubTree(root); return; }
void dendro::printSplitStats()      { splithist->printTreeAsList();		return; }
void dendro::printSplitStatsShort() { splithist->printTreeAsShortList();	return; }

void dendro::printSubTree(elementd *z) {
	if (z != NULL) {
		if (z->type == GRAPH) {
			cout << "[" << z->label << "]" << endl; //"\t(" << z->L << " " << z->R << ") - " << z->M <<  endl;
			return;
		} else if (z->type == DENDRO) {
			cout <<  "(p = " << z->p << "\te = " << z->e << "\tnL = " << z->L->n << "\tnR = " << z->R->n << "\tlabel = " << z->label << ")\tinternal[" << z->index << "]" << endl; // "\t(" << z->L << " " << z->R << ") - " << z->M <<  endl;
			cout << "L "; printSubTree(z->L); cout << endl;
			cout << "R "; printSubTree(z->R); cout << endl;
		} else {
			cout <<  "(p = " << z->p << "\te = " << z->e << "\tnL = " << z->L->n << "\tnR = " << z->R->n << "\tlabel = " << z->label << ")\tinternal[" << z->index << "] " << z->type << endl;
			cout << "L "; printSubTree(z->L); cout << endl;
			cout << "R "; printSubTree(z->R); cout << endl;
		}
	}
	return;
}

// ********************************************************************************************************

void dendro::refreshLikelihood() {
	// recalculates the log-likelihood of the dendrogram structure
	L = 0.0; double dL;
	int nL_nR, ei;
	for (int i=0; i<(n-1); i++) {
		nL_nR = internal[i].L->n*internal[i].R->n;
		ei    = internal[i].e;
		internal[i].p = (double)(ei) / (double)(nL_nR);
		if (ei == 0 or ei == nL_nR) { dL = 0.0; }
		else                        { dL = ei * log(internal[i].p) + (nL_nR - ei) * log(1.0-internal[i].p); }
		internal[i].logL = dL;
		L += dL;
	}
	return;
}

// ********************************************************************************************************

void dendro::QsortMain (block* array, int left, int right) {
	if (right > left) {
		int pivot = left;
		int part  = QsortPartition(array, left, right, pivot);
		QsortMain(array, left,   part-1);
		QsortMain(array, part+1, right  );
	}
	return;
}

int dendro::QsortPartition (block* array, int left, int right, int index) {
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

void dendro::recordConsensusTree(const string f_out) {
	
	keyValuePairSplit *curr, *prev;
	bool flag_debug = false;
	child *newChild;
	
	// First, cull the split hist so that only splits with weight >= 0.5 remain
	cullSplitHist();
	int treesize = splithist->returnNodecount();

	// Now, initialize the various arrays we use to keep track of the internal structure of the
	// consensus tree.
	ctree		 = new cnode   [treesize];
	cancestor       = new int     [n];
	for (int i=0; i<treesize; i++) { ctree[i].index = i;  }
	for (int i=0; i<n; i++)        { cancestor[i]   = -1; }
	int ii = 0;
	
	// To build the majority consensus tree, we do the following:
	// For each possible number of Ms in the split string (a number that ranges from n-2 down to 0),
	// and for each split with that number of Ms, we create a new internal node of the tree, and 
	// connect the oldest ancestor of each C to that node (at most once). Then, we update our list
	// of oldest ancestors to reflect this new join, and proceed.
	for (int i=n-2; i>=0; i--) {
		// First, we get a list of all the splits with this exactly i Ms
		curr = splithist->returnTheseSplits(i);
//		if (curr != NULL) { cout << ">> M[" << i << "] <<" << endl; }
		
		// Now we loop over that list
		while (curr != NULL) {
			splithist->deleteItem(curr->x);				// 
//			cout << curr->x << " (" << countChildren(curr->x) << ") " << curr->y << endl;
			ctree[ii].weight = curr->y;					// add weight to this internal node
			for (int j=0; j<n; j++) {					// examine each letter of this split
				if (curr->x[j] == 'C') {					// - node is child of this internal node
					if (cancestor[j] == -1) {				// - first time this leaf has ever been seen
						newChild        = new child;		//   
						newChild->type  = GRAPH;			// 
						newChild->index = j;			// 
						newChild->next  = NULL;
						if (ctree[ii].lastChild == NULL) {	// - attach child to list
							ctree[ii].children  = newChild;
							ctree[ii].lastChild = newChild;
							ctree[ii].degree    = 1;
						} else {
							ctree[ii].lastChild->next = newChild;
							ctree[ii].lastChild       = newChild;
							ctree[ii].degree   += 1;
						}
					} else {							// - this leaf has been seen before
						// If the parent of the ancestor of this leaf is the current internal node
						// then this leaf is already a descendant of this internal node, and we 
						// can move on; otherwise, we need to add that ancestor to this internal
						// node's child list, and update various relations
						if (ctree[cancestor[j]].parent != ii) {
							ctree[cancestor[j]].parent = ii;	// 
							newChild        = new child;		// 
							newChild->type  = DENDRO;		// 
							newChild->index = cancestor[j];	// 
							newChild->next  = NULL;
							if (ctree[ii].lastChild == NULL) {	// - attach child to list
								ctree[ii].children  = newChild;
								ctree[ii].lastChild = newChild;
								ctree[ii].degree    = 1;
							} else {
								ctree[ii].lastChild->next = newChild;
								ctree[ii].lastChild       = newChild;
								ctree[ii].degree   += 1;
							}
						}
					}
					cancestor[j] = ii;					// note new ancestry for this leaf
				}
			}
			ii++;									// update internal node index
			if (flag_debug) { printConsensusTree(); }
			prev = curr;
			curr = curr->next;
			delete prev;
		}
	}
	
	// write consensus tree structure to file
	child *sit, *sat;
	ofstream fout(f_out.c_str(), ios::trunc);
	for (int i=0; i<ii; i++) {
		fout << "[ " << i << " ]\t" << ctree[i].weight << "\tP= " << ctree[i].parent << "\tN= " << ctree[i].degree << " ";
		sit = ctree[i].children;
		while (sit != NULL) {
			if (sit->type == GRAPH) { fout << g->getName(sit->index) << " (G) "; }
			else                    { fout << sit->index << " (D) ";             }
			sat = sit;
			sit = sit->next;
			delete sat;
		}
		fout << "\n";
	}
	for (int i=0; i<n; i++) {
		if (cancestor[i] == -1) {
			fout << "[ " << ii++ << " ]\t1.000000\tP= " << -1 << "\tN= " << 1 << " " << g->getName(i) << " (G)\n";
		}
	}
	fout.close();
	cout << ">> exported consensus tree ( " << f_out << " )" << endl; 
	
	return;
	
}

void dendro::recordConsensusTree(igraph_vector_t *parents, igraph_vector_t *weights) {
  
  keyValuePairSplit *curr, *prev;
  bool flag_debug = false;
  child *newChild;
  int orig_nodes=g->numNodes();
  
  // First, cull the split hist so that only splits with weight >= 0.5
  // remain 
  cullSplitHist();
  int treesize = splithist->returnNodecount();

  // Now, initialize the various arrays we use to keep track of the
  // internal structure of the consensus tree.
  ctree		 = new cnode   [treesize];
  cancestor       = new int     [n];
  for (int i=0; i<treesize; i++) { ctree[i].index = i;  }
  for (int i=0; i<n; i++)        { cancestor[i]   = -1; }
  int ii = 0;
	
  // To build the majority consensus tree, we do the following: For
  // each possible number of Ms in the split string (a number that
  // ranges from n-2 down to 0), and for each split with that number
  // of Ms, we create a new internal node of the tree, and connect the
  // oldest ancestor of each C to that node (at most once). Then, we
  // update our list of oldest ancestors to reflect this new join, and
  // proceed.
  for (int i=n-2; i>=0; i--) {
    // First, we get a list of all the splits with this exactly i Ms
    curr = splithist->returnTheseSplits(i);
    // if (curr != NULL) { cout << ">> M[" << i << "] <<" << endl; }
    
    // Now we loop over that list
    while (curr != NULL) {
      splithist->deleteItem(curr->x);
      // add weight to this internal node
      ctree[ii].weight = curr->y; 
      // examine each letter of this split
      for (int j=0; j<n; j++) {	
	if (curr->x[j] == 'C') { 
	  // - node is child of this internal node
	  if (cancestor[j] == -1) {
	     // - first time this leaf has ever been seen
	    newChild        = new child;
	    newChild->type  = GRAPH;
	    newChild->index = j;
	    newChild->next  = NULL;
	    // - attach child to list
	    if (ctree[ii].lastChild == NULL) {	
	      ctree[ii].children  = newChild;
	      ctree[ii].lastChild = newChild;
	      ctree[ii].degree    = 1;
	    } else {
	      ctree[ii].lastChild->next = newChild;
	      ctree[ii].lastChild       = newChild;
	      ctree[ii].degree   += 1;
	    }
	  } else {
	    // - this leaf has been seen before
	    // If the parent of the ancestor of this leaf is the
	    // current internal node then this leaf is already a
	    // descendant of this internal node, and we can move on;
	    // otherwise, we need to add that ancestor to this
	    // internal node's child list, and update various
	    // relations
	    if (ctree[cancestor[j]].parent != ii) {
	      ctree[cancestor[j]].parent = ii;
	      newChild        = new child;
	      newChild->type  = DENDRO;
	      newChild->index = cancestor[j];
	      newChild->next  = NULL;
	      // - attach child to list
	      if (ctree[ii].lastChild == NULL) {
		ctree[ii].children  = newChild;
		ctree[ii].lastChild = newChild;
		ctree[ii].degree    = 1;
	      } else {
		ctree[ii].lastChild->next = newChild;
		ctree[ii].lastChild       = newChild;
		ctree[ii].degree   += 1;
	      }
	    }
	  }
	  // note new ancestry for this leaf
	  cancestor[j] = ii;
	}
      }
      // update internal node index
      ii++;
      if (flag_debug) { printConsensusTree(); }
      prev = curr;
      curr = curr->next;
      delete prev;
    }
  }

  // Return the consensus tree
  igraph_vector_resize(parents, ii + orig_nodes);
  if (weights) { igraph_vector_resize(weights, ii); }

  for (int i=0; i<ii; i++) {
    child *sat, *sit=ctree[i].children;    
    while (sit) {
      VECTOR(*parents)[orig_nodes + i] = 
	ctree[i].parent < 0 ? -1 : orig_nodes + ctree[i].parent;
      if (sit->type == GRAPH) { 
	VECTOR(*parents)[sit->index] = orig_nodes + i;
      }
      sat=sit;
      sit=sit->next;
      delete sat;
    }
    if (weights) { VECTOR(*weights)[i] = ctree[i].weight; }
  }
  
  // Plus the isolate nodes
  for (int i=0; i<n; i++) {
    if (cancestor[i] == -1) {
      VECTOR(*parents)[i] = -1;
    }
  }    


}

// ********************************************************************************************************

void dendro::recordDendrogramStructure(const string out_file) {
	
	ofstream fout(out_file.c_str(), ios::trunc);
	for (int i=0; i<(n-1); i++) {
		fout << "[ " << i << " ] ";
		fout << "L= " << internal[i].L->index << " "; if (internal[i].L->type == DENDRO) { fout << "(D) "; } else { fout << "(G) "; }
		fout << "R= " << internal[i].R->index << " "; if (internal[i].R->type == DENDRO) { fout << "(D) "; } else { fout << "(G) "; }
		fout << "p= " << internal[i].p << " ";
		fout << "e= " << internal[i].e << " ";
		fout << "n= " << internal[i].n << "\n";
	}
	fout.close();
	
	return;
}

// ********************************************************************************************************

void dendro::recordDendrogramStructure(igraph_hrg_t *hrg) {
  for (int i=0; i<n-1; i++) {
    int li=internal[i].L->index;
    int ri=internal[i].R->index;
    VECTOR(hrg->left )[i] = internal[i].L->type == DENDRO ? -li-1 : li;
    VECTOR(hrg->right)[i] = internal[i].R->type == DENDRO ? -ri-1 : ri;
    VECTOR(hrg->prob )[i] = internal[i].p;
    VECTOR(hrg->edges)[i] = internal[i].e;
    VECTOR(hrg->vertices)[i] = internal[i].n;    
  }
}


void dendro::recordGraphStructure(const string out_file) {
	edge* curr;
	string thisName;
	bool flag_debug = false;
	if (flag_debug) { cout << ">> dendro: writing random graph to file" << endl; }
	
	ofstream fout(out_file.c_str(), ios::trunc);
	for (int i=0; i<n; i++) {
		curr     = g->getNeighborList(i);
		thisName = g->getName(i);
		while (curr != NULL) {
			if (thisName=="") { fout << i << "\t" << curr->x << "\n"; }
			else {              fout << thisName << "\t" << g->getName(curr->x) << "\n"; }
			curr = curr->next;
		}
	}
	fout.close();
	
	return;
}

void dendro::recordGraphStructure(igraph_t *graph) {
  igraph_vector_t edges;
  int no_of_nodes=g->numNodes();
  int no_of_edges=g->numLinks() / 2;
  int idx=0;

  igraph_vector_init(&edges, no_of_edges*2);
  IGRAPH_FINALLY(igraph_vector_destroy, &edges);

  for (int i=0; i<n; i++) {
    edge *curr = g->getNeighborList(i);
    while (curr) {
      if (i < curr->x) {
	VECTOR(edges)[idx++] = i;
	VECTOR(edges)[idx++] = curr->x;
      }
      curr = curr->next;
    }
  }
  
  igraph_create(graph, &edges, no_of_nodes, /* directed= */ 0);

  igraph_vector_destroy(&edges);
  IGRAPH_FINALLY_CLEAN(1);
}

// ********************************************************************************************************

void dendro::recordSplitHistogram(const string out_file) {
	splithist->recordTreeAsList(out_file, 0.01); // exclude splits with < 0.01 weight
	return;
}

// ********************************************************************************************************

list* dendro::reversePathToRoot(const int leafIndex) {
	list *head, *subhead, *newlist;
	head = subhead = newlist = NULL;
	elementd *current = &leaf[leafIndex];
	
	while (current != NULL) {				// continue until we're finished
		newlist       = new list;			// add this node to the path
		newlist->x    = current->index;
		newlist->next = NULL;
		if (head == NULL) { head    = newlist; }
		else              { subhead = head;    head = newlist; head->next = subhead; }
		current = current->M;
	}
	return head;
}

// ********************************************************************************************************

bool	dendro::sampleSplitLikelihoods(int &sample_num) {
	// In order to compute the majority agreement dendrogram at equilibrium, we need to calculate
	// the leaf partition defined by each split (internal edge) of the tree. Because splits are 
	// only defined on a Cayley tree, the buildSplit() function returns the default "--...--" 
	// string for the root and the root's left child. When tabulating the frequency of splits, 
	// one of these needs to be excluded.
	bool flag_debug = false;
	string* array;
	int     k;
	double  tot;

	if (flag_debug) { cout << "dendro:: decomposing dendrogram into its splits" << endl; }
	string new_split;
	// To decompose the tree into its splits, we simply loop over all the internal nodes and
	// replace the old split for the ith internal node with its new split. This is a bit
	// time consuming to do O(n^2), so try not to do this very often. Once the decomposition
	// is had, we insert them into the split histogram, which tracks the cumulative weight
	// for each respective split observed.
	
	if (splithist == NULL) { splithist = new splittree; }
	for (int i=0; i<(n-1); i++) {
		new_split = buildSplit(&internal[i]);
		d->replaceSplit(i, new_split);
		if (new_split != "" and new_split[1] != '-') {
			if (!splithist->insertItem(new_split, 1.0)) { return false; } }
	}
	splithist->finishedThisRound();
	
	// For large graphs, the split histogram can get extremely large, so we need to employ some 
	// measures to prevent it from swamping the available memory. When the number of splits exceeds 
	// a threshold (say, a million), we progressively delete splits that have a weight less than 
	// a rising (k*0.001 of the total weight) fraction of the splits, on the assumption that losing 
	// such weight is unlikely to effect the ultimate split statistics. This deletion procedure is 
	// slow O(m lg m), but should only happen very rarely.
	int split_max = n*500;
	int leng;
	if (splithist->returnNodecount() > split_max) {
		
		k=1;
		// cout << "Culling the splithist = " << splithist->returnNodecount() << " -> ";
		while (splithist->returnNodecount() > split_max) {
			array = splithist->returnArrayOfKeys();
			tot   = splithist->returnTotal();
			leng  = splithist->returnNodecount();
			for (int i=0; i<leng; i++) {
				if ((splithist->returnValue(array[i]) / tot) < k*0.001) { splithist->deleteItem(array[i]); }
			}
			delete [] array; array = NULL;
			k++;
		}
		cout << splithist->returnNodecount() << endl;
//		cin >> pauseme;
		
		// An alternative is to just bail-out of the sampling when we exceed the split_max threshold
//		cout << "WARNING: maximum number of observed splits (" << split_max << ") exceeded.\n";
//		cout << "         Halting sampling at " << sample_num << " samples.\n";
//		sample_num = 2000000000;	
	}
	
//	splithist->printTreeAsShortList();
//	splithist->printTreeAsList();
//	cin >> pauseme;
	return true;
}

void	dendro::sampleAdjacencyLikelihoods() {
	// Here, we sample the probability values associated with every adjacency in A, weighted by 
	// their likelihood. The weighted histogram is stored in the graph data structure, so we simply 
	// need to add an observation to each node-pair that corresponds to the associated branch point's
	// probability and the dendrogram's overall likelihood.
	bool   flag_debug = false;
	double nn;
	double norm = ((double)(n) * (double)(n)) / 4.0;
	
	if (flag_debug) { cout << "dendro:: tabulating A'" << endl; }
	if (L > 0.0) { L = 0.0; }
	elementd* ancestor;
	list	*currL, *prevL;
	if (paths != NULL) { for (int i=0; i<n; i++) { currL = paths[i]; while (currL != NULL) { prevL = currL;   currL = currL->next;   delete prevL;   prevL = NULL; } paths[i] = NULL; } delete [] paths; } paths = NULL;
	paths = new list* [n];
	for (int i=0; i<n; i++) { paths[i] = reversePathToRoot(i); }	// construct paths from root, O(n^2) at worst
	for (int i=0; i<n; i++) {								// add obs for every node-pair, always O(n^2)
		for (int j=i+1; j<n; j++) {
			ancestor = findCommonAncestor(paths, i, j);			// find internal node, O(n) at worst
			nn       = ((double)(ancestor->L->n) * (double)(ancestor->R->n)) / norm;
			g->addAdjacencyObs(i, j, ancestor->p, nn);			// add obs of ->p to (i,j) histogram, and
			g->addAdjacencyObs(j, i, ancestor->p, nn);			// add obs of ->p to (j,i) histogram
		}												// 
	}													// 
	g->addAdjacencyEnd();									// finish-up: upate total weight in histograms
	
	return;
}

void dendro::resetDendrograph() {
       // Reset the dendrograph structure for the next trial
       if (leaf      != NULL) { delete [] leaf;     leaf      = NULL; }    // O(n)
       if (internal  != NULL) { delete [] internal; internal  = NULL; }    // O(n)
       if (d         != NULL) { delete d;                      d         = NULL; }    // O(n)
       root = NULL;
       if (paths != NULL) {
               list *curr, *prev; for (int i=0; i<n; i++) {
                       curr = paths[i]; while (curr != NULL) { prev = curr;   curr = curr->next;   delete prev;   prev = NULL; } paths[i] = NULL; }
               delete [] paths;
       } paths = NULL;
       L = 1.0;
        
       return;
 }


// ********************************************************************************************************
// ********************************************************************************************************

// ********************************************************************************************************
// ********************************************************************************************************

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

// ******** Constructor / Destructor **********************************************************************

graph::graph(const int size, bool predict) : predict(predict)  {
  n			= size;
  m			= 0;
  nodes		= new vert  [n];
  nodeLink		= new edge* [n];
  nodeLinkTail   = new edge* [n];
  for (int i=0; i<n; i++) { nodeLink[i] = NULL; nodeLinkTail[i] = NULL; }
  if (predict) {
    A = new double** [n];
    for (int i=0; i<n; i++) {
      A[i] = new double* [n];
    }
    obs_count = 0;
    total_weight = 0.0;
    bin_resolution = 0.0;
    num_bins = 0;
  }
}

graph::~graph() {
	edge *curr, *prev;
	for (int i=0; i<n; i++) {
		curr = nodeLink[i];
		while (curr != NULL) {
			prev = curr;
			curr = curr->next;
			delete prev;
		}
	}
	delete [] nodeLink;		nodeLink		= NULL;
	delete [] nodeLinkTail;  nodeLinkTail   = NULL;
	delete [] nodes;		nodes		= NULL;
	
	if (predict) {
	  for (int i=0; i<n; i++) {
	    for (int j=0; j<n; j++) { delete [] A[i][j]; }
	    delete [] A[i];
	  }
	  delete [] A;			A			= NULL;
	}
}

// ********************************************************************************************************

bool graph::addLink(const int i, const int j) {
	// Adds the directed edge (i,j) to the adjacency list for v_i
	edge* newedge;
	if (i >= 0 and i < n and j >= 0 and j < n) {
		newedge	 = new edge;
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
		return true;
	} else { return false; }
}

// ********************************************************************************************************

bool graph::addAdjacencyObs(const int i, const int j, const double probability, const double size) {
	// Adds the observation obs to the histogram of the edge (i,j)
	// Note: user must manually add observation to edge (j,i) by calling this function with that argument
	if (bin_resolution > 0.0 and probability >= 0.0 and probability <= 1.0 
	    and size >= 0.0 and size <= 1.0 and i >= 0 and i < n and j >= 0 and j < n) {
		int index = (int)(round(probability/bin_resolution));
		if (index < 0) { index = 0; } else if (index > num_bins) { index = num_bins; }
		// Add the weight to the proper probability bin
		if (A[i][j][index] < 0.5) { A[i][j][index] = 1.0; } else { A[i][j][index] += 1.0; }
		return true;
	}
	return false;
}

// ********************************************************************************************************

void graph::addAdjacencyEnd() {
	// We need to also keep a running total of how much weight has been added
	// to the histogram, and the number of observations in the histogram.
	if (obs_count==0) { total_weight  = 1.0; obs_count = 1; }
	else              { total_weight += 1.0; obs_count++;   }
	return;
}

bool graph::doesLinkExist(const int i, const int j) {
	// This function determines if the edge (i,j) already exists in the adjacency list of v_i
	edge* curr;
	if (i >= 0 and i < n and j >= 0 and j < n) {
		curr = nodeLink[i];
		while (curr != NULL) {
			if (curr->x == j) { return true; }
			curr = curr->next;
		}
	}
	return false;
}

// ********************************************************************************************************

int    graph::getDegree(const int i)       { if (i >= 0 and i < n) { return nodes[i].degree; } else { return -1;   } }
string graph::getName(const int i)         { if (i >= 0 and i < n) { return nodes[i].name;   } else { return "";   } }
// NOTE: Returns address; deallocation of returned object is dangerous
edge*  graph::getNeighborList(const int i) { if (i >= 0 and i < n) { return nodeLink[i];     } else { return NULL; } }

double* graph::getAdjacencyHist(const int i, const int j) {
  if (i >= 0 and i < n and j >= 0 and j < n) { return A[i][j]; } else { return NULL; }
}

// ********************************************************************************************************
double graph::getAdjacencyAverage(const int i, const int j) {
       double average = 0.0;
       if (i != j) {
               for (int k=0; k<num_bins; k++) {
                       if (A[i][j][k] > 0.0) { average += (A[i][j][k] / total_weight)*((double)(k)*bin_resolution); }
               }
       }
       return average;
}


int    graph::numLinks()      { return m; }
int    graph::numNodes()      { return n; }
double graph::getBinResolution() { return bin_resolution; }
int	  graph::getNumBins()       { return num_bins; }
double graph::getTotalWeight()   { return total_weight; }

// ********************************************************************************************************

void graph::printPairs() {
	edge* curr;
	for (int i=0; i<n; i++) {
		cout << "[" << i << "]\t";
		curr = nodeLink[i];
		while (curr != NULL) {
			cout << curr->x << "\t";
			curr = curr->next;
		}
		cout << "\n";
	}
	return;
}

// ********************************************************************************************************

void graph::resetAllAdjacencies() {
       for (int i=0; i<n; i++) { for (int j=0; j<n; j++) { for (int k=0; k<num_bins; k++) { A[i][j][k] = 0.0; } } }
       obs_count    = 0;
       total_weight = 0.0;
       return;
}

// ********************************************************************************************************

void graph::resetAdjacencyHistogram(const int i, const int j) {
       if (i >= 0 and i < n and j >= 0 and j < n) { for (int k=0; k<num_bins; k++) { A[i][j][k] = 0.0; } }
       return;
}

// ********************************************************************************************************

void graph::resetLinks() {
       edge *curr, *prev;
       for (int i=0; i<n; i++) {
               curr = nodeLink[i];
               while (curr != NULL) {
                       prev = curr;
                       curr = curr->next;
                       delete prev;
               }
               nodeLink[i]     = NULL;
               nodeLinkTail[i] = NULL;
               nodes[i].degree = 0;
       }
       m = 0;
       return;
}

// ********************************************************************************************************

void graph::setAdjacencyHistograms(const int bin_count) {
       // For all possible adjacencies, setup an edge histograms
       num_bins = bin_count+1;
       bin_resolution = 1.0 / (double)(bin_count);
       for (int i=0; i<n; i++) {
               for (int j=0; j<n; j++) {
                       A[i][j] = new double [num_bins];
                       for (int k=0; k<num_bins; k++) { A[i][j][k] = 0.0; }
               }
       }
        return;
 }


bool graph::setName(const int i, const string text) { if (i >= 0 and i < n) { nodes[i].name = text; return true; } else { return false; } }

// ********************************************************************************************************
// ********************************************************************************************************

// ********************************************************************************************************

// ********************************************************************************************************

interns::interns(const int n)  {
	q         = n;
	count     = 0;
	edgelist  = new ipair  [q];
	splitlist = new string [q+1];
	indexLUT  = new int*   [q+1];
	for (int i=0; i<(q+1); i++) {
		indexLUT[i]    = new int [2];
		indexLUT[i][0] = indexLUT[i][1] = -1;
	}
}
interns::~interns() {
	delete [] edgelist;
	delete [] splitlist;
	for (int i=0; i<(q+1); i++) { delete [] indexLUT[i]; }
	delete [] indexLUT;
}

// ********************************************************************************************************

// NOTE: Returns an address to another object -- do not deallocate
ipair* interns::getEdge(const int i) { return &edgelist[i]; }

// ********************************************************************************************************

// NOTE: Returns an address to another object -- do not deallocate
ipair* interns::getRandomEdge() { return &edgelist[(int)(floor((double)(q)*mtr.randExc()))]; }

// ********************************************************************************************************

string interns::getSplit(const int i) { if (i >= 0 and i <= q) { return splitlist[i]; } else { return ""; } }

// ********************************************************************************************************

bool interns::addEdge(const int new_x, const int new_y, const short int new_type) {
	// This function adds a new edge (i,j,t,sp) to the list of internal edges. After checking that the inputs
	// fall in the appropriate range of values, it records the new edgelist index in the indexLUT and then
	// puts the input values into that edgelist location.
	if (count < q and new_x >= 0 and new_x < (q+1) and new_y >= 0 and new_y < (q+2) and (new_type == LEFT or new_type == RIGHT)) {
		if (new_type == LEFT) { indexLUT[new_x][0] = count; } else { indexLUT[new_x][1] = count; }
		edgelist[count].x = new_x;
		edgelist[count].y = new_y;
		edgelist[count].t = new_type;
		count++;
		return true;
	} else { return false; }
}

// ********************************************************************************************************

void interns::printEdgeList() {
	for (int i=0; i<q; i++) {
		cout << "(" << edgelist[i].x << " " << edgelist[i].y << " ";
		if (edgelist[i].t == LEFT)  { cout << "L) "; } else
		if (edgelist[i].t == RIGHT) { cout << "R) "; } else { cout << "?) "; }
	}
	cout << endl;
	return;
}

// ********************************************************************************************************

void interns::printSplitList() { for (int i=0; i<=q; i++) { cout << "internal[" << i << "] = " << splitlist[i] << endl; } return; }

// ********************************************************************************************************

bool	interns::replaceSplit(const int i, const string sp) {
	// When an internal edge is changed, its split must be replaced as well. This function provides
	// that access; it stores the split defined by an internal edge (x,y) at the location [y], which
	// is unique.

	if (i >= 0 and i <= q) { splitlist[i] = sp; return true; }
	return false;
}

// ********************************************************************************************************

bool interns::swapEdges(const int one_x, const int one_y, const short int one_type, const int two_x, 
				    const int two_y, const short int two_type) {
	// The moves on the dendrogram always swap edges, either of which (or both, or neither) can by internal
	// edges. So, this function mirrors that operation for the internal edgelist and indexLUT.
	
	int index, jndex, temp;
	bool one_isInternal = false;
	bool two_isInternal = false;
	
	if (one_x >= 0 and one_x < (q+1) and two_x >= 0 and two_x < (q+1) and (two_type == LEFT or two_type == RIGHT) and 
	    one_y >= 0 and one_y < (q+2) and two_y >= 0 and two_y < (q+2) and (one_type == LEFT or one_type == RIGHT)) {
		
		if (one_type              == LEFT) { temp = 0; } else { temp = 1; }
		if (indexLUT[one_x][temp] >  -1)   { one_isInternal = true;       }
		if (two_type              == LEFT) { temp = 0; } else { temp = 1; }
		if (indexLUT[two_x][temp] >  -1)   { two_isInternal = true;       }
		
		if (one_isInternal and two_isInternal) {
			if (one_type == LEFT)  { index = indexLUT[one_x][0]; } else { index = indexLUT[one_x][1]; }
			if (two_type == LEFT)  { jndex = indexLUT[two_x][0]; } else { jndex = indexLUT[two_x][1]; }
			temp              = edgelist[index].y;
			edgelist[index].y = edgelist[jndex].y;
			edgelist[jndex].y = temp;
			
		} else if (one_isInternal) {
			if (one_type == LEFT)  { index = indexLUT[one_x][0]; indexLUT[one_x][0] = -1; }
			else                   { index = indexLUT[one_x][1]; indexLUT[one_x][1] = -1; }
			edgelist[index].x = two_x;
			edgelist[index].t = two_type;
			if (two_type == LEFT) { indexLUT[two_x][0] = index; } else { indexLUT[two_x][1] = index; } // add new
			
		} else if (two_isInternal) {
			if (two_type == LEFT)  { index = indexLUT[two_x][0]; indexLUT[two_x][0] = -1; }
			else                   { index = indexLUT[two_x][1]; indexLUT[two_x][1] = -1; }
			edgelist[index].x = one_x;
			edgelist[index].t = one_type;
			if (one_type == LEFT) { indexLUT[one_x][0] = index; } else { indexLUT[one_x][1] = index; } // add new
		} else { } // else neither is internal

		return true;
	} else { return false; }
}

// ******** Red-Black Tree Methods ************************************************************************

splittree::splittree() {
	root = new elementsp;
	leaf = new elementsp;

	leaf->parent   = root;

	root->left	= leaf;
	root->right    = leaf;
	support		= 0;
	total_weight	= 0.0;
	total_count	= 0;
}

splittree::~splittree() {
	if (root != NULL && (root->left != leaf || root->right != leaf)) { deleteSubTree(root); }
	support      = 0;
	total_weight = 0.0;
	total_count  = 0;
	delete leaf;
	root		   = NULL;
	leaf		   = NULL;
}

void splittree::deleteTree() { if (root != NULL) { deleteSubTree(root); } return; }

void splittree::deleteSubTree(elementsp *z) {

	if (z->left  != leaf) { deleteSubTree(z->left);  }
	if (z->right != leaf) { deleteSubTree(z->right); }
	delete z;
	z = NULL;
	return;
}

// ******** Reset Functions *******************************************************************************

void splittree::clearTree() { // O(n lg n)
	string *array = returnArrayOfKeys();
	for (int i=0; i<support; i++) { deleteItem(array[i]); }
	delete [] array;
	return;
}

// ******** Search Functions ******************************************************************************
// public search function - if there exists a elementsp in the tree with key=searchKey,
// it returns TRUE and foundNode is set to point to the found node; otherwise, it sets
// foundNode=NULL and returns FALSE
elementsp* splittree::findItem(const string searchKey) {

	elementsp *current;    current = root;
	if (current->split=="") { return NULL; }						// empty tree; bail out
	while (current != leaf) {
		if (searchKey < current->split) {							// left-or-right?
			if (current->left  != leaf) { current = current->left;  }	// try moving down-left
			else { return NULL; }								//   failure; bail out
		} else {												// 
			if (searchKey > current->split) {							// left-or-right?
				if (current->right  != leaf) { current = current->right;  }	// try moving down-left
				else { return NULL; }							//   failure; bail out
			} else { return current; }							// found (searchKey==current->split)
		}
	}
	return NULL;
}

double splittree::returnValue(const string searchKey) {
	elementsp* test = findItem(searchKey);
	if (test == NULL) { return 0.0; } else { return test->weight; }
}


// ******** Return Item Functions *************************************************************************
// public function which returns the tree, via pre-order traversal, as a linked list

string* splittree::returnArrayOfKeys() {
	string* array;
	array = new string [support];
	bool flag_go = true;
	int index = 0;
	elementsp *curr;
	
	if (support == 1) { array[0] = root->split; }
	else if (support == 2) {
		array[0] = root->split;
		if (root->left == leaf) { array[1] = root->right->split; } 
		else { array[1] = root->left->split; }
	} else {
		for (int i=0; i<support; i++) { array[i] = -1; }
		// non-recursive traversal of tree structure
		curr		 = root;
		curr->mark = 1;
		while (flag_go) {
			
			if (curr->mark == 1 and curr->left == leaf) {		// - is it time, and is left child the leaf node?
				curr->mark = 2;							// 
			}
			if (curr->mark == 2 and curr->right == leaf) {		// - is it time, and is right child the leaf node?
				curr->mark = 3;							// 
			}
			if (curr->mark == 1) {							// - go left
				curr->mark = 2;							// 
				curr       = curr->left;						// 
				curr->mark = 1;							// 
			} else if (curr->mark == 2) {						// - else go right
				curr->mark = 3;							// 
				curr       = curr->right;					// 
				curr->mark = 1;							// 
			} else {										// - else go up a level
				curr->mark = 0;							// 
				array[index++] = curr->split;					// 
				curr = curr->parent;						// 
				if (curr == NULL) { flag_go = false; }			// 
			}
		}
	}
	
	return array;
} // This does not leak memory (unlike returnListOfKeys)

slist* splittree::returnListOfKeys() {
	keyValuePairSplit	*curr, *prev;
	slist			*head, *tail, *newlist;

	curr = returnTreeAsList();
	while (curr != NULL) {
		newlist    = new slist;
		newlist->x = curr->x;
		if (head == NULL) { head       = newlist; tail = head;    }
		else              { tail->next = newlist; tail = newlist; }
		prev = curr;
		curr = curr->next;
		delete prev;
		prev = NULL;
	}
	return head;
}

keyValuePairSplit* splittree::returnTreeAsList() { // pre-order traversal
	keyValuePairSplit  *head, *tail;

	head    = new keyValuePairSplit;
	head->x = root->split;
	head->y = root->weight;
	head->c = root->count;
	tail    = head;

	if (root->left  != leaf) { tail = returnSubtreeAsList(root->left,  tail); }
	if (root->right != leaf) { tail = returnSubtreeAsList(root->right, tail); }
	
	if (head->x == "") { return NULL; /* empty tree */ } else { return head; }
}

keyValuePairSplit* splittree::returnSubtreeAsList(elementsp *z, keyValuePairSplit *head) {
	keyValuePairSplit *newnode, *tail;
	
	newnode    = new keyValuePairSplit;
	newnode->x = z->split;
	newnode->y = z->weight;
	newnode->c = z->count;
	head->next = newnode;
	tail       = newnode;
	
	if (z->left  != leaf) { tail = returnSubtreeAsList(z->left,  tail); }
	if (z->right != leaf) { tail = returnSubtreeAsList(z->right, tail); }
	
	return tail;
}

keyValuePairSplit splittree::returnMaxKey() {
	keyValuePairSplit themax;
	elementsp *current;
	current = root;
	while (current->right != leaf) {		// search to bottom-right corner of tree
		current = current->right; }		// 
	themax.x = current->split;			// store the data found
	themax.y = current->weight;			// 
	
	return themax;						// return that data
}

keyValuePairSplit splittree::returnMinKey() {
	keyValuePairSplit themin;
	elementsp *current;
	current = root;
	while (current->left != leaf) {		// search to bottom-left corner of tree
		current = current->left; }		// 
	themin.x = current->split;			// store the data found
	themin.y = current->weight;			// 
	
	return themin;						// return that data
}

// private functions for deleteItem() (although these could easily be made public, I suppose)
elementsp* splittree::returnMinKey(elementsp *z) {
	elementsp *current;

	current = z;
	while (current->left != leaf) {		// search to bottom-right corner of tree
		current = current->left; }		// 
	return current;					// return pointer to the minimum
}

elementsp* splittree::returnSuccessor(elementsp *z) {
	elementsp *current, *w;
	
	w = z;
	if (w->right != leaf) {				// if right-subtree exists, return min of it
		return returnMinKey(w->right); }
	current = w->parent;				// else search up in tree
	while ((current!=NULL) && (w==current->right)) {
		w       = current;
		current = current->parent;		// move up in tree until find a non-right-child
	}
	return current;
}

int splittree::returnNodecount() { return support; }

keyValuePairSplit* splittree::returnTheseSplits(const int target) {
	keyValuePairSplit *head, *curr, *prev, *newhead, *newtail, *newpair;
	int count, len;
	
	head = returnTreeAsList();
	prev = newhead = newtail = newpair = NULL;
	curr = head;
	
	while (curr != NULL) {
		count = 0;
		len   = curr->x.size();
		for (int i=0; i<len; i++) { if (curr->x[i] == 'M') { count++; } }
		if (count == target and curr->x[1] != '*') {
			newpair       = new keyValuePairSplit;
			newpair->x    = curr->x;
			newpair->y    = curr->y;
			newpair->next = NULL;
			if (newhead == NULL) { newhead = newpair; newtail = newpair; }
			else { newtail->next = newpair; newtail = newpair; }
		}
		prev = curr;
		curr = curr->next;
		delete prev;
		prev = NULL;
	}
	
	return newhead;
}

double splittree::returnTotal() { return total_weight; }

// ******** Insert Functions ******************************************************************************

void splittree::finishedThisRound() {
	// We need to also keep a running total of how much weight has been added to the histogram.
	if (total_count == 0) { total_weight  = 1.0; total_count = 1; }
	else				  { total_weight += 1.0; total_count++; }
	return;
}	

// public insert function
bool splittree::insertItem(string newKey, double newValue) {
	
	// first we check to see if newKey is already present in the tree; if so, we do nothing;
	// if not, we must find where to insert the key
	elementsp *newNode, *current;
	
	current = findItem(newKey);						// find newKey in tree; return pointer to it O(log k)
	if (current != NULL) {
		// Add weight to the existing item's weight
		if (current->weight > total_weight) {
			cout << "ERROR: " << current->weight << " > " << total_weight << endl;
			cout << "for split " << current->split << endl;
			return false;
		}
		current->weight += 1.0;
		// And finally, we keep track of how many observations went into the histogram
		current->count++;
		return true;
	} else {
		newNode			= new elementsp;			// elementsp for the splittree
		newNode->split		= newKey;					//  store newKey
		newNode->weight	= newValue;  				//  store newValue
		newNode->color		= true;					//  new nodes are always RED
		newNode->parent	= NULL;					//  new node initially has no parent
		newNode->left		= leaf;					//  left leaf
		newNode->right		= leaf;					//  right leaf
		newNode->count		= 1;						// 
		support++;								// increment node count in splittree
		
		// must now search for where to insert newNode, i.e., find the correct parent and
		// set the parent and child to point to each other properly
		current = root;
		if (current->split=="") {									// insert as root
			delete root;											//   delete old root
			root			= newNode;								//   set root to newNode
			leaf->parent   = newNode;								//   set leaf's parent
			current		= leaf;									//   skip next loop
		}
		
		while (current != leaf) {									// search for insertion point
			if (newKey < current->split) {							// left-or-right?
				if (current->left  != leaf) { current = current->left;  }	// try moving down-left
				else {											// else found new parent
					newNode->parent	= current;					//    set parent
					current->left		= newNode;					//    set child
					current			= leaf;						//    exit search
				}
			} else {												// 
				if (current->right != leaf) { current = current->right; }   // try moving down-right
				else {											// else found new parent
					newNode->parent	= current;					//    set parent
					current->right		= newNode;					//    set child
					current			= leaf;						//    exit search
				}
			}
		}

		// now do the house-keeping necessary to preserve the red-black properties
		insertCleanup(newNode);			// do house-keeping to maintain balance
		
	}
	return true;
}

// private house-keeping function for insertion
void splittree::insertCleanup(elementsp *z) {
	
	if (z->parent==NULL) {						// fix now if z is root
		z->color = false; return; }
	elementsp *temp;
	while (z->parent!=NULL && z->parent->color) {	// while z is not root and z's parent is RED
		if (z->parent == z->parent->parent->left) {  // z's parent is LEFT-CHILD
			temp = z->parent->parent->right;		// grab z's uncle
			if (temp->color) {
				z->parent->color		= false;  // color z's parent BLACK	(Case 1)
				temp->color			= false;  // color z's uncle BLACK		(Case 1)
				z->parent->parent->color = true;   // color z's grandparent RED  (Case 1)
				z = z->parent->parent;			// set z = z's grandparent    (Case 1)
			} else {
				if (z == z->parent->right) {		// z is RIGHT-CHILD
					z = z->parent;				// set z = z's parent		(Case 2)
					rotateLeft(z);				// perform left-rotation		(Case 2)
				}
				z->parent->color		= false;  // color z's parent BLACK	(Case 3)
				z->parent->parent->color = true;   // color z's grandparent RED  (Case 3)
				rotateRight(z->parent->parent);    // perform right-rotation	(Case 3)
			}
		} else {								// z's parent is RIGHT-CHILD
			temp = z->parent->parent->left;		// grab z's uncle
			if (temp->color) {
				z->parent->color		= false;  // color z's parent BLACK	(Case 1)
				temp->color			= false;  // color z's uncle BLACK		(Case 1)
				z->parent->parent->color = true;   // color z's grandparent RED  (Case 1)
				z = z->parent->parent;			// set z = z's grandparent    (Case 1)
			} else {
				if (z == z->parent->left) {		// z is LEFT-CHILD
					z = z->parent;				// set z = z's parent		(Case 2)
					rotateRight(z);			// perform right-rotation	(Case 2)
				}
				z->parent->color		= false;  // color z's parent BLACK	(Case 3)
				z->parent->parent->color = true;   // color z's grandparent RED  (Case 3)
				rotateLeft(z->parent->parent);	// perform left-rotation		(Case 3)
			}
		}
	}

	root->color = false;						// color the root BLACK
	return;
}

// ******** Delete Functions ******************************************************************************
// public delete function
void splittree::deleteItem(string killKey) {
	elementsp *x, *y, *z;
	
	z = findItem(killKey);
	if (z == NULL) { return; }						// item not present; bail out

	if (support==1) {								// -- attempt to delete the root
		root->split	= "";						// restore root node to default state
		root->weight    = 0.0;						// 
		root->color    = false;						// 
		root->parent   = NULL;						// 
		root->left	= leaf;						// 
		root->right    = leaf;						// 
		support--;								// set support to zero
		total_weight   = 0.0;						// set total weight to zero
		total_count--;								// 
		return;									// exit - no more work to do
	}
	
	if (z != NULL) {
		support--;								// decrement node count
		if ((z->left == leaf) || (z->right==leaf)) {		// case of less than two children
			  y = z; }							//    set y to be z
		else { y = returnSuccessor(z); }				//    set y to be z's key-successor
		
		if (y->left!=leaf) { x = y->left; }			// pick y's one child (left-child)
		else			    { x = y->right; }			//				  (right-child)
		x->parent = y->parent;						// make y's child's parent be y's parent

		if (y->parent==NULL) { root = x; }				// if y is the root, x is now root
		else {									// 
			if (y == y->parent->left) {				// decide y's relationship with y's parent
				y->parent->left  = x;				//   replace x as y's parent's left child
			} else {								// 
				y->parent->right = x; }				//   replace x as y's parent's left child
		}										// 

		if (y!=z) {								// insert y into z's spot
			z->split		= y->split;				// copy y data into z
			z->weight		= y->weight;				// 
			z->count		= y->count;				// 
		}										// 

		if (y->color==false) { deleteCleanup(x); }		// do house-keeping to maintain balance
		delete y;									// deallocate y
		y = NULL;									// point y to NULL for safety
	}											// 
		
	return;
}

void splittree::deleteCleanup(elementsp *x) {
	elementsp *w, *t;
	while ((x != root) && (x->color==false)) {			// until x is the root, or x is RED
		if (x==x->parent->left) {					// branch on x being a LEFT-CHILD
			w = x->parent->right;					// grab x's sibling
			if (w->color==true) {					// if x's sibling is RED
				w->color = false;					// color w BLACK				(case 1)
				x->parent->color = true;				// color x's parent RED			(case 1)
				rotateLeft(x->parent);				// left rotation on x's parent	(case 1)
				w = x->parent->right;				// make w be x's right sibling	(case 1)
			}
			if ((w->left->color==false) && (w->right->color==false)) {
				w->color = true;					// color w RED					(case 2)
				x = x->parent;						// examine x's parent			(case 2)
			} else {								// 
				if (w->right->color==false) {			// 
					w->left->color = false;			// color w's left child BLACK		(case 3)
					w->color = true;				// color w RED					(case 3)
					t = x->parent;					// store x's parent
					rotateRight(w);				// right rotation on w			(case 3)
					x->parent = t;					// restore x's parent
					w = x->parent->right;			// make w be x's right sibling	(case 3)
				}								// 
				w->color			= x->parent->color; // make w's color = x's parent's   (case 4)
				x->parent->color    = false;			// color x's parent BLACK		(case 4)
				w->right->color	= false;			// color w's right child BLACK	(case 4)
				rotateLeft(x->parent);				// left rotation on x's parent	(case 4)
				x = root;							// finished work. bail out		(case 4)
			}									// 
		} else {									// x is RIGHT-CHILD
			w = x->parent->left;					// grab x's sibling
			if (w->color==true) {					// if x's sibling is RED
				w->color			= false;			// color w BLACK				(case 1)
				x->parent->color    = true;			// color x's parent RED			(case 1)
				rotateRight(x->parent);				// right rotation on x's parent	(case 1)
				w				= x->parent->left;  // make w be x's left sibling		(case 1)
			}
			if ((w->right->color==false) && (w->left->color==false)) {
				w->color = true;					// color w RED					(case 2)
				x= x->parent;						// examine x's parent			(case 2)
			} else {								// 
				if (w->left->color==false) {			// 
					w->right->color	= false;		// color w's right child BLACK	(case 3)
					w->color			= true;		// color w RED					(case 3)
					t				= x->parent;   // store x's parent
					rotateLeft(w);					// left rotation on w			(case 3)
					x->parent			= t;			// restore x's parent
					w = x->parent->left;			// make w be x's left sibling		(case 3)
				}								// 
				w->color = x->parent->color;			// make w's color = x's parent's   (case 4)
				x->parent->color    = false;			// color x's parent BLACK		(case 4)
				w->left->color		= false;			// color w's left child BLACK		(case 4)
				rotateRight(x->parent);				// right rotation on x's parent    (case 4)
				x				= root;			// x is now the root			(case 4)
			}
		}
	}
	x->color = false;								// color x (the root) BLACK		(exit)

	return;
}

// ******** Rotation Functions ****************************************************************************

void splittree::rotateLeft(elementsp *x) {
	elementsp *y;
	// do pointer-swapping operations for left-rotation
	y               = x->right;					// grab right child
	x->right        = y->left;					// make x's RIGHT-CHILD be y's LEFT-CHILD
	y->left->parent = x;						// make x be y's LEFT-CHILD's parent
	y->parent       = x->parent;					// make y's new parent be x's old parent

	if (x->parent==NULL) { root = y; }				// if x was root, make y root
	else {									// 
		if (x == x->parent->left)				// if x is LEFT-CHILD, make y be x's parent's
			{ x->parent->left  = y; }			//    left-child
		else { x->parent->right = y; }			//    right-child
	}										// 
	y->left   = x;								// make x be y's LEFT-CHILD
	x->parent = y;								// make y be x's parent
	
	return;
}

void splittree::rotateRight(elementsp *y) {
	elementsp *x;
	// do pointer-swapping operations for right-rotation
	x                = y->left;					// grab left child
	y->left          = x->right;					// replace left child yith x's right subtree
	x->right->parent = y;						// replace y as x's right subtree's parent
	
	x->parent        = y->parent;					// make x's new parent be y's old parent
	if (y->parent==NULL) { root = x; }				// if y was root, make x root
	else {
		if (y == y->parent->right)				// if y is RIGHT-CHILD, make x be y's parent's
			{ y->parent->right  = x; }			//    right-child
		else { y->parent->left   = x; }			//    left-child
	}
	x->right  = y;								// make y be x's RIGHT-CHILD
	y->parent = x;								// make x be y's parent
	
	return;
}

// ******** Display Functions *****************************************************************************
// public
void splittree::printTree() {
	cout << "\nTREE SIZE = " << support << "\t" << total_weight << "\t" << total_count << endl;
	cout << "# "; printSubTree(root);
	return;
}

// private
void splittree::printSubTree(elementsp *z) {
	if (z==leaf) { return; }
	else {
		cout << "(" << z->split << " " << z->weight << " " << z->count << " " << z->color << ")"<<endl;
		cout << "L "; printSubTree(z->left); cout << endl;
		cout << "R "; printSubTree(z->right); cout << endl;
	}
	return;
}

// public
void	splittree::printTreeAsList() {
	keyValuePairSplit *curr, *prev;
	curr = returnTreeAsList();
	int temp = 0;
	while (curr != NULL) {
		if (curr->y > 5) {
			cout << curr->x;
			if (curr->y >= 0.5*total_weight) { cout << "\t* "; temp++; } else { cout << "\t  "; }
			cout << curr->y << "\t" << curr->c / total_weight << "\n";
		}
		prev = curr;
		curr = curr->next;
		delete prev; prev = NULL;
	}
	curr = NULL;
	cout << "total_count   = " << total_count  << endl;
	cout << "total_weight  = " << total_weight << endl;
	cout << "total_size    = " << support << endl;
	cout << "weight >= 1/2 = " << temp << endl;

	return;
}

void	splittree::printTreeAsShortList() {
	keyValuePairSplit *curr, *prev;
	curr = returnTreeAsList();
	cout << "numSplits = " << support << endl;
	int temp = 0;
	while (curr != NULL) {
		if (curr->y / total_weight >= 0.5) {
			cout << curr->x << "\t* ";
			temp++;
			cout << curr->y << "\t" << curr->y / total_weight << "\n";
		}
		prev = curr;
		curr = curr->next;
		delete prev; prev = NULL;
	}
	curr = NULL;
	cout << "maximum = " << total_weight << endl;
	cout << "weight >= 1/2 = " << temp << endl;
	
	return;
}

void splittree::recordTreeAsList(const string file_out, const double threshold) {
	keyValuePairSplit *curr, *prev;
	curr = returnTreeAsList();
	
	ofstream fout(file_out.c_str(), ios::trunc);
	fout.precision(8);
//	for (int i=0; i<len; i++) { fout << "*"; } fout << "\t" << total_weight << "\n";
	while (curr != NULL) {
		if (curr->y / total_weight >= threshold) { fout << curr->x << "\t" << curr->y / total_weight << "\n"; }
		prev = curr;
		curr = curr->next;
		delete prev; prev = NULL;
	}
	curr = NULL;
	fout.close();
	
	return;
}

