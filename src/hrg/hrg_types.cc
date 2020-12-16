// ***********************************************************************
// *** COPYRIGHT NOTICE **************************************************
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
// ***********************************************************************
// Author       : Aaron Clauset  ( aaronc@santafe.edu |
//                                 http://www.santafe.edu/~aaronc/ )
// Collaborators: Cristopher Moore and Mark Newman
// Project      : Hierarchical Random Graphs
// Location     : University of New Mexico, Dept. of Computer Science
//                AND Santa Fe Institute
// Created      : Spring 2004
// Modified     : many, many times
//
// ***********************************************************************

#include "hrg/rbtree.h"
#include "hrg/dendro.h"
#include "hrg/graph.h"
#include "hrg/splittree_eq.h"
#include "hrg/graph_simp.h"

#include "igraph_hrg.h"
#include "igraph_constructors.h"
#include "igraph_random.h"

using namespace std;
using namespace fitHRG;

// ******** Red-Black Tree Methods ***************************************

rbtree::rbtree() {
    root = new elementrb;
    leaf = new elementrb;

    leaf->parent = root;

    root->left = leaf;
    root->right = leaf;
    support = 0;
}

rbtree::~rbtree() {
    if (root != NULL &&
        (root->left != leaf || root->right != leaf)) {
        deleteSubTree(root);
    }
    if (root) {
        delete root;
    }
    delete leaf;
    support = 0;
    root = 0;
    leaf = 0;
}

void rbtree::deleteTree() {
    if (root != NULL) {
        deleteSubTree(root);
    }
} // does not leak memory

void rbtree::deleteSubTree(elementrb *z) {
    if (z->left  != leaf) {
        deleteSubTree(z->left);
    }
    if (z->right != leaf) {
        deleteSubTree(z->right);
    }
    delete z;
}

// ******** Search Functions *********************************************
// public search function - if there exists a elementrb in the tree
// with key=searchKey, it returns TRUE and foundNode is set to point
// to the found node; otherwise, it sets foundNode=NULL and returns
// FALSE
elementrb* rbtree::findItem(const int searchKey) {
    elementrb *current = root;

    // empty tree; bail out
    if (current->key == -1) {
        return NULL;
    }

    while (current != leaf) {
        // left-or-right?
        if (searchKey < current->key) {
            // try moving down-left
            if (current->left  != leaf) {
                current = current->left;
            } else {
                //   failure; bail out
                return NULL;
            }
        } else {
            // left-or-right?
            if (searchKey > current->key) {
                // try moving down-left
                if (current->right  != leaf) {
                    current = current->right;
                } else {
                    // failure; bail out
                    return NULL;
                }
            } else {
                // found (searchKey==current->key)
                return current;
            }
        }
    }
    return NULL;
}

int rbtree::returnValue(const int searchKey) {
    elementrb* test = findItem(searchKey);
    if (!test) {
        return 0;
    } else {
        return test->value;
    }
}


// ******** Return Item Functions ****************************************

int* rbtree::returnArrayOfKeys() {
    int* array;
    array = new int [support];
    bool flag_go = true;
    int index = 0;
    elementrb *curr;

    if (support == 1) {
        array[0] = root->key;
    } else if (support == 2) {
        array[0] = root->key;
        if (root->left == leaf) {
            array[1] = root->right->key;
        } else {
            array[1] = root->left->key;
        }
    } else {
        for (int i = 0; i < support; i++) {
            array[i] = -1;
        }
        // non-recursive traversal of tree structure
        curr = root;
        curr->mark = 1;
        while (flag_go) {
            // - is it time, and is left child the leaf node?
            if (curr->mark == 1 && curr->left == leaf) {
                curr->mark = 2;
            }
            // - is it time, and is right child the leaf node?
            if (curr->mark == 2 && curr->right == leaf) {
                curr->mark = 3;
            }
            if (curr->mark == 1) {
                // - go left
                curr->mark = 2;
                curr = curr->left;
                curr->mark = 1;
            } else if (curr->mark == 2) {
                // - else go right
                curr->mark = 3;
                curr = curr->right;
                curr->mark = 1;
            } else {
                // - else go up a level
                curr->mark = 0;
                array[index++] = curr->key;
                curr = curr->parent;
                if (curr == NULL) {
                    flag_go = false;
                }
            }
        }
    }

    return array;
}

list* rbtree::returnListOfKeys() {
    keyValuePair *curr, *prev;
    list *head = 0, *tail = 0, *newlist;

    curr = returnTreeAsList();
    while (curr != NULL) {
        newlist = new list;
        newlist->x = curr->x;
        if (head == NULL) {
            head       = newlist; tail = head;
        } else {
            tail->next = newlist; tail = newlist;
        }
        prev = curr;
        curr = curr->next;
        delete prev;
        prev = NULL;
    }
    return head;
}

keyValuePair* rbtree::returnTreeAsList() {
    // pre-order traversal
    keyValuePair  *head, *tail;

    head = new keyValuePair;
    head->x = root->key;
    head->y = root->value;
    tail = head;

    if (root->left  != leaf) {
        tail = returnSubtreeAsList(root->left,  tail);
    }
    if (root->right != leaf) {
        tail = returnSubtreeAsList(root->right, tail);
    }

    if (head->x == -1) {
        return NULL; /* empty tree */
    } else {
        return head;
    }
}

keyValuePair* rbtree::returnSubtreeAsList(elementrb *z, keyValuePair *head) {
    keyValuePair *newnode, *tail;

    newnode = new keyValuePair;
    newnode->x = z->key;
    newnode->y = z->value;
    head->next = newnode;
    tail = newnode;

    if (z->left  != leaf) {
        tail = returnSubtreeAsList(z->left,  tail);
    }
    if (z->right != leaf) {
        tail = returnSubtreeAsList(z->right, tail);
    }

    return tail;
}

keyValuePair rbtree::returnMaxKey() {
    keyValuePair themax;
    elementrb *current;
    current  = root;

    // search to bottom-right corner of tree
    while (current->right != leaf) {
        current  = current->right;
    }
    themax.x = current->key;
    themax.y = current->value;

    return themax;
}

keyValuePair rbtree::returnMinKey() {
    keyValuePair themin;
    elementrb *current;
    current = root;
    // search to bottom-left corner of tree
    while (current->left != leaf) {
        current = current->left;
    }
    themin.x = current->key;
    themin.y = current->value;

    return themin;
}

// private functions for deleteItem() (although these could easily be
// made public, I suppose)
elementrb* rbtree::returnMinKey(elementrb *z) {
    elementrb *current;

    current = z;
    // search to bottom-right corner of tree
    while (current->left != leaf) {
        current = current->left;
    }
    return current;
}

elementrb* rbtree::returnSuccessor(elementrb *z) {
    elementrb *current, *w;

    w = z;
    // if right-subtree exists, return min of it
    if (w->right != leaf) {
        return returnMinKey(w->right);
    }
    // else search up in tree
    current = w->parent;
    while ((current != NULL) && (w == current->right)) {
        w = current;
        // move up in tree until find a non-right-child
        current = current->parent;
    }
    return current;
}

int rbtree::returnNodecount() {
    return support;
}

// ******** Insert Functions *********************************************
// public insert function
void rbtree::insertItem(int newKey, int newValue) {

    // first we check to see if newKey is already present in the tree;
    // if so, we do nothing; if not, we must find where to insert the
    // key
    elementrb *newNode, *current;

    // find newKey in tree; return pointer to it O(log k)
    current = findItem(newKey);
    if (current == NULL) {
        newNode = new elementrb;    // elementrb for the rbtree
        newNode->key = newKey;
        newNode->value = newValue;
        newNode->color = true;  // new nodes are always RED
        newNode->parent = NULL; // new node initially has no parent
        newNode->left = leaf;   // left leaf
        newNode->right = leaf;  // right leaf
        support++;          // increment node count in rbtree

        // must now search for where to insert newNode, i.e., find the
        // correct parent and set the parent and child to point to each
        // other properly
        current = root;
        if (current->key == -1) {    // insert as root
            delete root;           // delete old root
            root = newNode;            // set root to newNode
            leaf->parent = newNode;        // set leaf's parent
            current = leaf;            // skip next loop
        }

        // search for insertion point
        while (current != leaf) {
            // left-or-right?
            if (newKey < current->key) {
                // try moving down-left
                if (current->left  != leaf) {
                    current = current->left;
                } else {
                    // else found new parent
                    newNode->parent = current; // set parent
                    current->left = newNode;   // set child
                    current = leaf;        // exit search
                }
            } else {
                // try moving down-right
                if (current->right != leaf) {
                    current = current->right;
                } else {
                    // else found new parent
                    newNode->parent = current; // set parent
                    current->right = newNode;  // set child
                    current = leaf;        // exit search
                }
            }
        }

        // now do the house-keeping necessary to preserve the red-black
        // properties
        insertCleanup(newNode);
    }
    return;
}

// private house-keeping function for insertion
void rbtree::insertCleanup(elementrb *z) {

    // fix now if z is root
    if (z->parent == NULL) {
        z->color = false;
        return;
    }

    elementrb *temp;

    // while z is not root and z's parent is RED
    while (z->parent != NULL && z->parent->color) {
        if (z->parent == z->parent->parent->left) {

            // z's parent is LEFT-CHILD

            temp = z->parent->parent->right;   // grab z's uncle
            if (temp->color) {
                z->parent->color = false;        // color z's parent BLACK  (Case 1)
                temp->color = false;             // color z's uncle BLACK   (Case 1)
                z->parent->parent->color = true; // color z's grandpar. RED (Case 1)
                z = z->parent->parent;           // set z = z's grandparent (Case 1)
            } else {
                if (z == z->parent->right) {
                    // z is RIGHT-CHILD
                    z = z->parent;             // set z = z's parent      (Case 2)
                    rotateLeft(z);             // perform left-rotation   (Case 2)
                }
                z->parent->color = false;        // color z's parent BLACK  (Case 3)
                z->parent->parent->color = true; // color z's grandpar. RED (Case 3)
                rotateRight(z->parent->parent);  // perform right-rotation  (Case 3)
            }
        } else {

            // z's parent is RIGHT-CHILD

            temp = z->parent->parent->left;    // grab z's uncle
            if (temp->color) {
                z->parent->color = false;        // color z's parent BLACK  (Case 1)
                temp->color = false;             // color z's uncle BLACK   (Case 1)
                z->parent->parent->color = true; // color z's grandpar. RED (Case 1)
                z = z->parent->parent;           // set z = z's grandparent (Case 1)
            } else {
                if (z == z->parent->left) {
                    // z is LEFT-CHILD
                    z = z->parent;                 // set z = z's parent      (Case 2)
                    rotateRight(z);                // perform right-rotation  (Case 2)
                }
                z->parent->color = false;        // color z's parent BLACK  (Case 3)
                z->parent->parent->color = true; // color z's grandpar. RED (Case 3)
                rotateLeft(z->parent->parent);   // perform left-rotation   (Case 3)
            }
        }
    }

    root->color = false;               // color the root BLACK
    return;
}

// ******** Delete
// ******** Functions *********************************************

void rbtree::replaceItem(int key, int newValue) {
    elementrb* ptr;
    ptr = findItem(key);
    ptr->value = newValue;
    return;
}

void rbtree::incrementValue(int key) {
    elementrb* ptr;
    ptr = findItem(key);
    ptr->value = 1 + ptr->value;
    return;
}

// public delete function
void rbtree::deleteItem(int killKey) {
    elementrb *x, *y, *z;

    z = findItem(killKey);
    if (z == NULL) {
        return;    // item not present; bail out
    }

    if (support == 1) {     // attempt to delete the root
        root->key = -1;       // restore root node to default state
        root->value = -1;
        root->color = false;
        root->parent = NULL;
        root->left = leaf;
        root->right = leaf;
        support--;            // set support to zero
        return;           // exit - no more work to do
    }

    if (z != NULL) {
        support--;            // decrement node count
        if ((z->left == leaf) || (z->right == leaf)) {
            y = z;                      // case of less than two children,
            // set y to be z
        } else {
            y = returnSuccessor(z);     // set y to be z's key-successor
        }

        if (y->left != leaf) {
            x = y->left;        // pick y's one child (left-child)
        } else {
            x = y->right;       // (right-child)
        }
        x->parent = y->parent;        // make y's child's parent be y's parent

        if (y->parent == NULL) {
            root = x;           // if y is the root, x is now root
        } else {
            if (y == y->parent->left) { // decide y's relationship with y's parent
                y->parent->left  = x;     // replace x as y's parent's left child
            } else {
                y->parent->right = x;     // replace x as y's parent's left child
            }
        }

        if (y != z) {         // insert y into z's spot
            z->key = y->key;        // copy y data into z
            z->value = y->value;
        }

        // do house-keeping to maintain balance
        if (y->color == false) {
            deleteCleanup(x);
        }

        delete y;
        y = NULL;
    }

    return;
}

void rbtree::deleteCleanup(elementrb *x) {
    elementrb *w, *t;

    // until x is the root, or x is RED
    while ((x != root) && (x->color == false)) {
        if (x == x->parent->left) {   // branch on x being a LEFT-CHILD
            w = x->parent->right;   // grab x's sibling
            if (w->color == true) { // if x's sibling is RED
                w->color = false;     // color w BLACK (case 1)
                x->parent->color = true;  // color x's parent RED            (case 1)
                rotateLeft(x->parent);    // left rotation on x's parent     (case 1)
                w = x->parent->right;     // make w be x's right sibling     (case 1)
            }
            if ((w->left->color == false) && (w->right->color == false)) {
                w->color = true;      // color w RED                     (case 2)
                x = x->parent;        // examine x's parent              (case 2)
            } else {
                if (w->right->color == false) {
                    w->left->color = false; // color w's left child BLACK      (case 3)
                    w->color = true;    // color w RED                     (case 3)
                    t = x->parent;      // store x's parent                (case 3)
                    rotateRight(w);     // right rotation on w             (case 3)
                    x->parent = t;      // restore x's parent              (case 3)
                    w = x->parent->right;   // make w be x's right sibling     (case 3)
                }
                w->color = x->parent->color; // w's color := x's parent's    (case 4)
                x->parent->color = false; // color x's parent BLACK          (case 4)
                w->right->color = false;  // color w's right child BLACK     (case 4)
                rotateLeft(x->parent);    // left rotation on x's parent     (case 4)
                x = root;                 // finished work. bail out         (case 4)
            }
        } else {                // x is RIGHT-CHILD
            w = x->parent->left;      // grab x's sibling
            if (w->color == true) {   // if x's sibling is RED
                w->color = false;       // color w BLACK                 (case 1)
                x->parent->color    = true; // color x's parent RED          (case 1)
                rotateRight(x->parent);     // right rotation on x's parent  (case 1)
                w = x->parent->left;        // make w be x's left sibling    (case 1)
            }
            if ((w->right->color == false) && (w->left->color == false)) {
                w->color = true;        // color w RED                   (case 2)
                x = x->parent;          // examine x's parent            (case 2)
            } else {
                if (w->left->color == false) {
                    w->right->color = false;  // color w's right child BLACK   (case 3)
                    w->color = true;      // color w RED                   (case 3)
                    t = x->parent;        // store x's parent              (case 3)
                    rotateLeft(w);        // left rotation on w            (case 3)
                    x->parent = t;        // restore x's parent            (case 3)
                    w = x->parent->left;      // make w be x's left sibling    (case 3)
                }
                w->color = x->parent->color; // w's color := x's parent's    (case 4)
                x->parent->color    = false; // color x's parent BLACK       (case 4)
                w->left->color = false;      // color w's left child BLACK   (case 4)
                rotateRight(x->parent);      // right rotation on x's parent (case 4)
                x = root;            // x is now the root            (case 4)
            }
        }
    }
    x->color = false;          // color x (the root) BLACK (exit)

    return;
}

// ******** Rotation Functions ******************************************

void rbtree::rotateLeft(elementrb *x) {
    elementrb *y;
    // do pointer-swapping operations for left-rotation
    y = x->right;          // grab right child
    x->right = y->left;        // make x's RIGHT-CHILD be y's LEFT-CHILD
    y->left->parent = x;       // make x be y's LEFT-CHILD's parent
    y->parent = x->parent;     // make y's new parent be x's old parent

    if (x->parent == NULL) {
        root = y;           // if x was root, make y root
    } else {
        // if x is LEFT-CHILD, make y be x's parent's
        if (x == x->parent->left) {
            x->parent->left  = y; // left-child
        } else {
            x->parent->right = y; //  right-child
        }
    }
    y->left   = x;        // make x be y's LEFT-CHILD
    x->parent = y;        // make y be x's parent

    return;
}

void rbtree::rotateRight(elementrb *y) {
    elementrb *x;
    // do pointer-swapping operations for right-rotation
    x = y->left;        // grab left child
    y->left = x->right;     // replace left child yith x's right subtree
    x->right->parent = y;   // replace y as x's right subtree's parent

    x->parent = y->parent;  // make x's new parent be y's old parent

    // if y was root, make x root
    if (y->parent == NULL) {
        root = x;
    } else {
        // if y is RIGHT-CHILD, make x be y's parent's
        if (y == y->parent->right) {
            // right-child
            y->parent->right = x;
        } else {
            // left-child
            y->parent->left = x;
        }
    }
    x->right  = y;        // make y be x's RIGHT-CHILD
    y->parent = x;        // make x be y's parent

    return;
}

// ***********************************************************************
// *** COPYRIGHT NOTICE **************************************************
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
// ***********************************************************************
// Author       : Aaron Clauset  ( aaronc@santafe.edu |
//                                 http://www.santafe.edu/~aaronc/ )
// Collaborators: Cristopher Moore and Mark E.J. Newman
// Project      : Hierarchical Random Graphs
// Location     : University of New Mexico, Dept. of Computer Science
//                AND Santa Fe Institute
// Created      : 26 October 2005 - 7 December 2005
// Modified     : 23 December 2007 (cleaned up for public consumption)
//
// ***********************************************************************
//
// Maximum likelihood dendrogram data structure. This is the heart of
// the HRG algorithm: all manipulations are done here and all data is
// stored here. The data structure uses the separate graph data
// structure to store the basic adjacency information (in a
// dangerously mutable way).
//
// ***********************************************************************

// ******** Dendrogram Methods *******************************************

dendro::dendro(): root(0), internal(0), leaf(0), d(0), splithist(0),
    paths(0), ctree(0), cancestor(0), g(0) { }
dendro::~dendro() {
    list *curr, *prev;

    if (g)        {
        delete    g;            // O(m)
        g        = 0;
    }
    if (internal) {
        delete [] internal;     // O(n)
        internal = 0;
    }
    if (leaf)     {
        delete [] leaf;         // O(n)
        leaf     = 0;
    }
    if (d)        {
        delete    d;            // O(n)
        d        = 0;
    }
    if (splithist) {
        delete    splithist;    // potentially long
        splithist = 0;
    }

    if (paths) {
        for (int i = 0; i < n; i++) {
            curr = paths[i];
            while (curr) {
                prev = curr;
                curr = curr->next;
                delete prev;
                prev = 0;
            }
            paths[i] = 0;
        }
        delete [] paths;
    }
    paths = 0;

    if (ctree)    {
        delete [] ctree;        // O(n)
        ctree     = 0;
    }
    if (cancestor) {
        delete [] cancestor;    // O(n)
        cancestor = 0;
    }
}

// *********************************************************************

void dendro::binarySearchInsert(elementd* x, elementd* y) {
    if (y->p < x->p) {        // go to left subtree
        if (x->L == NULL) {     // check if left subtree is empty
            x->L = y;         // make x left child
            y->M = x;         // make y parent of child
            return;
        } else {
            binarySearchInsert(x->L, y);
        }
    } else {          // go to right subtree
        if (x->R == NULL) {     // check if right subtree is empty
            x->R = y;         // make x right child
            y->M = x;         // make y parent of child
            return;
        } else {
            binarySearchInsert(x->R, y);
        }
    }
    return;
}

// **********************************************************************

list* dendro::binarySearchFind(const double v) {
    list *head = NULL, *tail = NULL, *newlist;
    elementd *current = root;
    bool flag_stopSearch = false;

    while (!flag_stopSearch) {    // continue until we're finished
        newlist    = new list;  // add this node to the path
        newlist->x = current->label;
        if (current == root) {
            head = newlist; tail = head;
        } else {
            tail->next = newlist; tail = newlist;
        }
        if (v < current->p) {   // now try left subtree
            if (current->L->type == GRAPH) {
                flag_stopSearch = true;
            } else {
                current = current->L;
            }
        } else {            // else try right subtree
            if (current->R->type == GRAPH) {
                flag_stopSearch = true;
            } else  {
                current = current->R;
            }
        }
    }
    return head;
}

// ***********************************************************************

string dendro::buildSplit(elementd* thisNode) {
    // A "split" is defined as the bipartition of vertices into the sets
    // of leaves below the internal vertex in the tree (denoted by "C"),
    // and those above it (denoted as "M"). For simplicity, we represent
    // this bipartition as a character string of length n, where the ith
    // character denotes the partition membership (C,M) of the ith leaf
    // node.

    bool flag_go = true;
    const short int k = 1 + DENDRO + GRAPH;
    elementd* curr;
    split sp;

    sp.initializeSplit(n);      // default split string O(n)

    curr = thisNode;        // - set start node as top this sub-tree
    curr->type = k + 1;     // - initialize in-order tree traversal
    while (flag_go) {

        // - is it time, and is left child a graph node?
        if (curr->type == k + 1 && curr->L->type == GRAPH) {
            sp.s[curr->L->index] = 'C'; // - mark this leaf
            curr->type = k + 2;
        }

        // - is it time, and is right child a graph node?
        if (curr->type == k + 2 && curr->R->type == GRAPH) {
            sp.s[curr->R->index] = 'C'; // - mark this leaf
            curr->type           = k + 3;
        }
        if (curr->type == k + 1) {    // - go left
            curr->type = k + 2;
            curr       = curr->L;
            curr->type = k + 1;
        } else if (curr->type == k + 2) { // - else go right
            curr->type = k + 3;
            curr       = curr->R;
            curr->type = k + 1;
        } else {              // - else go up a level
            curr->type = DENDRO;
            if (curr->index == thisNode->index || curr->M == NULL) {
                flag_go = false; curr = NULL;
            } else {
                curr = curr->M;
            }
        }
    }

    // any leaf that was not already marked must be in the remainder of
    // the tree
    for (int i = 0; i < n; i++) {
        if (sp.s[i] != 'C') {
            sp.s[i] = 'M';
        }
    }

    return sp.s;
}

// **********************************************************************

void dendro::buildDendrogram() {

    /* the initialization of the dendrogram structure goes like this:
     * 1) we allocate space for the n-1 internal nodes of the
     *    dendrogram, and then the n leaf nodes
     * 2) we build a random binary tree structure out of the internal
     *    nodes by assigning each a uniformly random value over [0,1] and
     *    then inserting it into the tree according to the
     *    binary-search rule.
     * 3) next, we make a random permutation of the n leaf nodes and add
     *    them to the dendrogram D by replacing the emptpy spots in-order
     * 4) then, we compute the path from the root to each leaf and store
     *    that in each leaf (this is prep work for the next step)
     * 5) finally, we compute the values for nL, nR, e (and thus p) and
     *    the label for each internal node by allocating each of the m
     *    edges in g to the appropriate internal node
     */

    // --- Initialization and memory allocation for data structures
    // After allocating the memory for D and G, we need to mark the
    // nodes for G as being non-internal vertices, and then insert them
    // into a random binary tree structure. For simplicity, we make the
    // first internal node in the array the root.

    n = g->numNodes();          // size of graph
    leaf = new elementd [n];    // allocate memory for G, O(n)
    internal  = new elementd [n - 1]; // allocate memory for D, O(n)
    d = new interns(n - 2);         // allocate memory for internal
    // edges of D, O(n)
    for (int i = 0; i < n; i++) { // initialize leaf nodes
        leaf[i].type   = GRAPH;
        leaf[i].label  = i;
        leaf[i].index  = i;
        leaf[i].n = 1;
    }

// initialize internal nodes
    root = &internal[0];
    root->label = 0;
    root->index = 0;
    root->p = RNG_UNIF01();

    // insert remaining internal vertices, O(n log n)
    for (int i = 1; i < (n - 1); i++) {
        internal[i].label = i;
        internal[i].index = i;
        internal[i].p = RNG_UNIF01();
        binarySearchInsert(root, &internal[i]);
    }

    // --- Hang leaf nodes off end of dendrogram O(n log n)
    // To impose this random hierarchical relationship on G, we first
    // take a random permutation of the leaf vertices and then replace
    // the NULLs at the bottom of the tree in-order with the leafs. As a
    // hack to ensure that we can find the leafs later using a binary
    // search, we assign each of them the p value of their parent,
    // perturbed slightly so as to preserve the binary search property.

    block* array; array = new block [n];
    for (int i = 0; i < n; i++) {
        array[i].x = RNG_UNIF01();
        array[i].y = i;
    }
    QsortMain(array, 0, n - 1);

    int k = 0;        // replace NULLs with leaf nodes, and
    for (int i = 0; i < (n - 1); i++) { // maintain binary search property, O(n)
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

    // --- Compute the path from root -> leaf for each leaf O(n log n)
    // Using the binary search property, we can find each leaf node in
    // O(log n) time. The binarySearchFind() function returns the list
    // of internal node indices that the search crossed, in the order of
    // root -> ... -> leaf, for use in the subsequent few operations.

    if (paths != NULL) {
        list *curr, *prev;
        for (int i = 0; i < n; i++) {
            curr = paths[i];
            while (curr != NULL) {
                prev = curr;
                curr = curr->next;
                delete prev;
                prev = NULL;
            }
            paths[i] = NULL;
        }
        delete [] paths;
    }
    paths = NULL;
    paths = new list* [n];
    for (int i = 0; i < n; i++) {
        paths[i] = binarySearchFind(leaf[i].p);
    }

    // --- Count e for each internal node O(m)
    // To count the number of edges that span the L and R subtrees for
    // each internal node, we use the path information we just
    // computed. Then, we loop over all edges in G and find the common
    // ancestor in D of the two endpoints and increment that internal
    // node's e count. This process takes O(m) time because in a roughly
    // balanced binary tree (given by our random dendrogram), the vast
    // majority of vertices take basically constant time to find their
    // common ancestor. Note that because our adjacency list is
    // symmetric, we overcount each e by a factor of 2, so we need to
    // correct this after.

    elementd* ancestor; edge* curr;
    for (int i = 0; i < (n - 1); i++) {
        internal[i].e = 0;
        internal[i].label = -1;
    }
    for (int i = 0; i < n; i++) {
        curr = g->getNeighborList(i);
        while (curr != NULL) {
            ancestor = findCommonAncestor(paths, i, curr->x);
            ancestor->e += 1;
            curr = curr->next;
        }
    }
    for (int i = 0; i < (n - 1); i++) {
        internal[i].e /= 2;
    }

    // --- Count n for each internal node O(n log n)
    // To tabulate the number of leafs in each subtree rooted at an
    // internal node, we use the path information computed above.
    for (int i = 0; i < n; i++) {
        ancestor = &leaf[i];
        ancestor = ancestor->M;
        while (ancestor != NULL) {
            ancestor->n++;
            ancestor = ancestor->M;
        }
    }

    // --- Label all internal vertices O(n log n)
    // We want to label each internal vertex with the smallest leaf
    // index of its children. This will allow us to collapse many
    // leaf-orderings into a single dendrogram structure that is
    // independent of child-exhanges (since these have no impact on the
    // likelihood of the hierarchical structure). To do this, we loop
    // over the leaf vertices from smallest to largest and walk along
    // that leaf's path from the root. If we find an unlabeled internal
    // node, then we mark it with this leaf's index.

    for (int i = 0; i < n; i++) {
        ancestor = &leaf[i];
        while (ancestor != NULL) {
            if (ancestor->label == -1 || ancestor->label > leaf[i].label) {
                ancestor->label = leaf[i].label;
            }
            ancestor = ancestor->M;
        }
    }

    // --- Exchange children to enforce order-property O(n)
    // We state that the order-property requires that an internal node's
    // label is the smallest index of its left subtree. The dendrogram
    // so far doesn't reflect this, so we need to step through each
    // internal vertex and make that adjustment (swapping nL and nR if
    // we make a change).

    elementd *tempe;
    for (int i = 0; i < (n - 1); i++) {
        if (internal[i].L->label > internal[i].label) {
            tempe = internal[i].L;
            internal[i].L = internal[i].R;
            internal[i].R = tempe;
        }
    }

    // --- Tabulate internal dendrogram edges O(n^2)
    // For the MCMC moves later on, we'll need to be able to choose,
    // uniformly at random, an internal edge of the dendrogram to
    // manipulate. There are always n-2 of them, and we can find them
    // simply by scanning across the internal vertices and observing
    // which have children that are also internal vertices. Note: very
    // important that the order property be enforced before this step is
    // taken; otherwise, the internal edges wont reflect the actual
    // dendrogram structure.

    for (int i = 0; i < (n - 1); i++) {
        if (internal[i].L->type == DENDRO) {
            d->addEdge(i, internal[i].L->index, LEFT);
        }
        if (internal[i].R->type == DENDRO) {
            d->addEdge(i, internal[i].R->index, RIGHT);
        }
    }

    // --- Clear memory for paths O(n log n)
    // Now that we're finished using the paths, we need to deallocate
    // them manually.

    list *current, *previous;
    for (int i = 0; i < n; i++) {
        current = paths[i];
        while (current) {
            previous = current;
            current = current->next;
            delete previous;
            previous = NULL;
        }
        paths[i] = NULL;
    }
    delete [] paths;
    paths = NULL;

    // --- Compute p_i for each internal node O(n)
    // Each internal node's p_i = e_i / (nL_i*nR_i), and now that we
    // have each of those pieces, we may calculate this value for each
    // internal node. Given these, we can then calculate the
    // log-likelihood of the entire dendrogram structure \log(L) =
    // \sum_{i=1}^{n} ( ( e_i \log[p_i] ) + ( (nL_i*nR_i - e_i)
    // \log[1-p_i] ) )

    L = 0.0; double dL;
    int nL_nR, ei;
    for (int i = 0; i < (n - 1); i++) {
        nL_nR = internal[i].L->n * internal[i].R->n;
        ei = internal[i].e;
        internal[i].p = (double)(ei) / (double)(nL_nR);
        if (ei == 0 || ei == nL_nR) {
            dL = 0.0;
        } else {
            dL = ei * log(internal[i].p) + (nL_nR - ei) * log(1.0 - internal[i].p);
        }
        internal[i].logL = dL;
        L += dL;
    }

    for (int i = 0; i < (n - 1); i++) {
        if (internal[i].label > internal[i].L->label) {
            tempe = internal[i].L;
            internal[i].L = internal[i].R;
            internal[i].R = tempe;
        }
    }

    // Dendrogram is now built

    return;
}

// ***********************************************************************

void dendro::clearDendrograph() {
    // Clear out the memory and references used by the dendrograph
    // structure - this is  intended to be called just before an
    // importDendrogramStructure call so as to avoid memory leaks and
    // overwriting the references therein.

    if (g        != NULL) {
        delete    g;           // O(m)
        g        = NULL;
    }
    if (leaf     != NULL) {
        delete [] leaf;        // O(n)
        leaf     = NULL;
    }
    if (internal != NULL) {
        delete [] internal;    // O(n)
        internal = NULL;
    }
    if (d        != NULL) {
        delete    d;        // O(n)
        d           = NULL;
    }
    root = NULL;

    return;
}

// **********************************************************************

int dendro::computeEdgeCount(const int a, const short int atype,
                             const int b, const short int btype) {
    // This function computes the number of edges that cross between the
    // subtree internal[a] and the subtree internal[b]. To do this, we
    // use an array A[1..n] integers which take values -1 if A[i] is in
    // the subtree defined by internal[a], +1 if A[i] is in the subtree
    // internal[b], and 0 otherwise. Taking the smaller of the two sets,
    // we then scan over the edges attached  to that set of vertices and
    // count the number of endpoints we see in the other set.

    bool flag_go    = true;
    int nA, nB;
    int         count = 0;
    const short int k = 1 + DENDRO + GRAPH;

    elementd* curr;

    // First, we push the leaf nodes in the L and R subtrees into
    // balanced binary tree structures so that we can search them
    // quickly later on.

    if (atype == GRAPH) {
        // default case, subtree A is size 1
        // insert single node as member of left subtree
        subtreeL.insertItem(a, -1);
        nA = 1; //
    } else {
        // explore subtree A, O(|A|)
        curr  = &internal[a];
        curr->type = k + 1;
        nA = 0;
        while (flag_go) {
            if (curr->index == internal[a].M->index) {
                internal[a].type = DENDRO;
                flag_go = false;
            } else {
                // - is it time, and is left child a graph node?
                if (curr->type == k + 1 && curr->L->type == GRAPH) {
                    subtreeL.insertItem(curr->L->index, -1);
                    curr->type = k + 2;
                    nA++;
                }
                // - is it time, and is right child a graph node?
                if (curr->type == k + 2 && curr->R->type == GRAPH) {
                    subtreeL.insertItem(curr->R->index, -1);
                    curr->type = k + 3;
                    nA++;
                }
                if (curr->type == k + 1) {    // - go left
                    curr->type = k + 2;
                    curr       = curr->L;
                    curr->type = k + 1;
                } else if (curr->type == k + 2) { // - else go right
                    curr->type = k + 3;
                    curr       = curr->R;
                    curr->type = k + 1;
                } else {              // - else go up a level
                    curr->type = DENDRO;
                    curr       = curr->M;
                    if (curr == NULL) {
                        flag_go = false;
                    }
                }
            }
        }
    }

    if (btype == GRAPH) {
        // default case, subtree A is size 1
        // insert node as single member of right subtree
        subtreeR.insertItem(b, 1);
        nB = 1;
    } else {
        flag_go = true;
        // explore subtree B, O(|B|)
        curr = &internal[b];
        curr->type = k + 1;
        nB  = 0;
        while (flag_go) {
            if (curr->index == internal[b].M->index) {
                internal[b].type = DENDRO;
                flag_go = false;
            } else {
                // - is it time, and is left child a graph node?
                if (curr->type == k + 1 && curr->L->type == GRAPH) {
                    subtreeR.insertItem(curr->L->index, 1);
                    curr->type = k + 2;
                    nB++;
                }
                // - is it time, and is right child a graph node?
                if (curr->type == k + 2 && curr->R->type == GRAPH) {
                    subtreeR.insertItem(curr->R->index, 1);
                    curr->type = k + 3;
                    nB++;
                }
                if (curr->type == k + 1) {    // - look left
                    curr->type = k + 2;
                    curr       = curr->L;
                    curr->type = k + 1;
                } else if (curr->type == k + 2) { // - look right
                    curr->type = k + 3;
                    curr       = curr->R;
                    curr->type = k + 1;
                } else {              // - else go up a level
                    curr->type = DENDRO;
                    curr       = curr->M;
                    if (curr == NULL) {
                        flag_go = false;
                    }
                }
            }
        }
    }

    // Now, we take the smaller subtree and ask how many of its
    // emerging edges have their partner in the other subtree. O(|A| log
    // |A|) time

    edge* current;
    int*  treeList;
    if (nA < nB) {
        // subtreeL is smaller
        treeList = subtreeL.returnArrayOfKeys();
        for (int i = 0; i < nA; i++) {
            current = g->getNeighborList(treeList[i]);
            // loop over each of its neighbors v_j
            while (current != NULL) {
                // to see if v_j is in A
                if (subtreeR.findItem(current->x) != NULL) {
                    count++;
                }
                current = current->next;
            }
            subtreeL.deleteItem(treeList[i]);
        }
        delete [] treeList;
        treeList = subtreeR.returnArrayOfKeys();
        for (int i = 0; i < nB; i++) {
            subtreeR.deleteItem(treeList[i]);
        }
        delete [] treeList;
    } else {
        // subtreeR is smaller
        treeList = subtreeR.returnArrayOfKeys();
        for (int i = 0; i < nB; i++) {
            current = g->getNeighborList(treeList[i]);
            // loop over each of its neighbors v_j
            while (current != NULL) {
                // to see if v_j is in B
                if (subtreeL.findItem(current->x) != NULL) {
                    count++;
                }
                current = current->next;
            }
            subtreeR.deleteItem(treeList[i]);
        }
        delete [] treeList;
        treeList = subtreeL.returnArrayOfKeys();
        for (int i = 0; i < nA; i++) {
            subtreeL.deleteItem(treeList[i]);
        }
        delete [] treeList;
    }

    return count;
}

// ***********************************************************************

int dendro::countChildren(const string s) {
    int len = s.size();
    int numC = 0;
    for (int i = 0; i < len; i++) {
        if (s[i] == 'C') {
            numC++;
        }
    }
    return numC;
}

// ***********************************************************************

void dendro::cullSplitHist() {
    string* array;
    int tot, leng;

    array = splithist->returnArrayOfKeys();
    tot   = splithist->returnTotal();
    leng  = splithist->returnNodecount();
    for (int i = 0; i < leng; i++) {
        if ((splithist->returnValue(array[i]) / tot) < 0.5) {
            splithist->deleteItem(array[i]);
        }
    }
    delete [] array; array = NULL;

    return;
}

// **********************************************************************

elementd* dendro::findCommonAncestor(list** paths_, const int i, const int j) {
    list* headOne = paths_[i];
    list* headTwo = paths_[j];
    elementd* lastStep = NULL;
    while (headOne->x == headTwo->x) {
        lastStep = &internal[headOne->x];
        headOne  = headOne->next;
        headTwo  = headTwo->next;
        if (headOne == NULL || headTwo == NULL) {
            break;
        }
    }
    return lastStep; // Returns address of an internal node; do not deallocate
}

// **********************************************************************

int dendro::getConsensusSize() {
    string    *array;
    double     value, tot;
    int  numSplits, numCons;
    numSplits = splithist->returnNodecount();
    array     = splithist->returnArrayOfKeys();
    tot       = splithist->returnTotal();
    numCons = 0;
    for (int i = 0; i < numSplits; i++) {
        value = splithist->returnValue(array[i]);
        if (value / tot > 0.5) {
            numCons++;
        }
    }
    delete [] array; array = NULL;
    return numCons;
}

// **********************************************************************

splittree* dendro::getConsensusSplits() {
    string    *array;
    splittree *consensusTree;
    double     value, tot;
    consensusTree  = new splittree;
    int numSplits;

    // We look at all of the splits in our split histogram and add any
    // one that's in the majority to our consensusTree, which we then
    // return (note that consensusTree needs to be deallocated by the
    // user).
    numSplits = splithist->returnNodecount();
    array     = splithist->returnArrayOfKeys();
    tot       = splithist->returnTotal();
    for (int i = 0; i < numSplits; i++) {
        value = splithist->returnValue(array[i]);
        if (value / tot > 0.5) {
            consensusTree->insertItem(array[i], value / tot);
        }
    }
    delete [] array; array = NULL;
    return consensusTree;
}

// ***********************************************************************

double dendro::getLikelihood() {
    return L;
}

// ***********************************************************************

void dendro::getSplitList(splittree* split_tree) {
    string sp;
    for (int i = 0; i < (n - 1); i++) {
        sp = d->getSplit(i);
        if (!sp.empty() && sp[1] != '-') {
            split_tree->insertItem(sp, 0.0);
        }
    }
    return;
}

// ***********************************************************************

double dendro::getSplitTotalWeight() {
    if (splithist) {
        return splithist->returnTotal();
    } else {
        return 0;
    }
}

// ***********************************************************************

bool dendro::importDendrogramStructure(const igraph_hrg_t *hrg) {
    n = igraph_hrg_size(hrg);

    // allocate memory for G, O(n)
    leaf = new elementd[n];
    // allocate memory for D, O(n)
    internal = new elementd[n - 1];
    // allocate memory for internal edges of D, O(n)
    d = new interns(n - 2);

    // initialize leaf nodes
    for (int i = 0; i < n; i++) {
        leaf[i].type  = GRAPH;
        leaf[i].label = i;
        leaf[i].index = i;
        leaf[i].n     = 1;
    }

    // initialize internal nodes
    root = &internal[0];
    root->label = 0;
    for (int i = 1; i < n - 1; i++) {
        internal[i].index = i;
        internal[i].label = -1;
    }

    // import basic structure from hrg object, O(n)
    for (int i = 0; i < n - 1; i++) {
        int left_index = VECTOR(hrg->left)[i];
        int right_index = VECTOR(hrg->right)[i];

        if (left_index < 0) {
            internal[i].L = &internal[-left_index - 1];
            internal[-left_index - 1].M = &internal[i];
        } else {
            internal[i].L = &leaf[left_index];
            leaf[left_index].M = &internal[i];
        }

        if (right_index < 0) {
            internal[i].R = &internal[-right_index - 1];
            internal[-right_index - 1].M = &internal[i];
        } else {
            internal[i].R = &leaf[right_index];
            leaf[right_index].M = &internal[i];
        }

        internal[i].p = VECTOR(hrg->prob)[i];
        internal[i].e = VECTOR(hrg->edges)[i];
        internal[i].n = VECTOR(hrg->vertices)[i];
        internal[i].index = i;
    }

    // --- Label all internal vertices O(n log n)
    elementd *curr;
    for (int i = 0; i < n; i++) {
        curr = &leaf[i];
        while (curr) {
            if (curr->label == -1 || curr->label > leaf[i].label) {
                curr->label = leaf[i].label;
            }
            curr = curr -> M;
        }
    }

    // --- Exchange children to enforce order-property O(n)
    elementd *tempe;
    for (int i = 0; i < n - 1; i++) {
        if (internal[i].L->label > internal[i].label) {
            tempe          = internal[i].L;
            internal[i].L  = internal[i].R;
            internal[i].R  = tempe;
        }
    }

    // --- Tabulate internal dendrogram edges O(n)
    for (int i = 0; i < (n - 1); i++) {
        if (internal[i].L->type == DENDRO) {
            d->addEdge(i, internal[i].L->index, LEFT);
        }
        if (internal[i].R->type == DENDRO) {
            d->addEdge(i, internal[i].R->index, RIGHT);
        }
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
    for (int i = 0; i < (n - 1); i++) {
        nL_nR = internal[i].L->n * internal[i].R->n;
        ei    = internal[i].e;
        if (ei == 0 || ei == nL_nR) {
            dL = 0.0;
        } else {
            dL = (double)(ei) * log(internal[i].p) +
                 (double)(nL_nR - ei) * log(1.0 - internal[i].p);
        }
        internal[i].logL = dL;
        L += dL;
    }

    return true;
}

// ***********************************************************************

void dendro::makeRandomGraph() {
    if (g != NULL) {
        delete g;
    } g = NULL; g = new graph(n);

    list *curr, *prev;
    if (paths) {
        for (int i = 0; i < n; i++) {
            curr = paths[i];
            while (curr != NULL) {
                prev = curr;
                curr = curr->next;
                delete prev;
                prev = NULL;
            }
            paths[i] = NULL;
        }
        delete [] paths;
    }
// build paths from root O(n d)
    paths = new list* [n];
    for (int i = 0; i < n; i++) {
        paths[i] = reversePathToRoot(i);
    }

    elementd* commonAncestor;
// O((h+d)*n^2) - h: height of D; d: average degree in G
    for (int i = 0; i < n; i++) {
        // decide neighbors of v_i
        for (int j = (i + 1); j < n; j++) {
            commonAncestor = findCommonAncestor(paths, i, j);
            if (RNG_UNIF01() < commonAncestor->p) {
                if (!(g->doesLinkExist(i, j))) {
                    g->addLink(i, j);
                }
                if (!(g->doesLinkExist(j, i))) {
                    g->addLink(j, i);
                }
            }
        }
    }

    for (int i = 0; i < n; i++) {
        curr = paths[i];
        while (curr != NULL) {
            prev = curr;
            curr = curr->next;
            delete prev;
            prev = NULL;
        }
        paths[i] = NULL;
    }
    delete [] paths; // delete paths data structure O(n log n)
    paths = NULL;

    return;
}

// **********************************************************************

bool dendro::monteCarloMove(double& delta, bool& ftaken, const double T) {
    // A single MC move begins with the selection of a random internal
    // edge (a,b) of the dendrogram. This also determines the three
    // subtrees i, j, k that we will rearrange, and we choose uniformly
    // from among the options.
    //
    // If (a,b) is a left-edge, then we have ((i,j),k), and moves
    // ((i,j),k) -> ((i,k),j) (alpha move)
    //           -> (i,(j,k)) + enforce order-property for (j,k) (beta move)
    //
    // If (a,b) is a right-edge, then we have (i,(j,k)), and moves
    // (i,(j,k)) -> ((i,k),j) (alpha move)
    //           -> ((i,j),k) (beta move)
    //
    // For each of these moves, we need to know what the change in
    // likelihood will be, so that we can determine with what
    // probability we execute the move.

    elementd *temp;
    ipair *tempPair;
    int x, y, e_x, e_y, n_i, n_j, n_k, n_x, n_y;
    short int t;
    double p_x, p_y, L_x, L_y, dLogL;
    string new_split;

    // The remainder of the code executes a single MCMC move, where we
    // sample the dendrograms proportionally to their likelihoods (i.e.,
    // temperature=1, if you're comparing it to the usual MCMC
    // framework).

    delta    = 0.0;
    ftaken   = false;
    tempPair = d->getRandomEdge(); // returns address; no need to deallocate
    x        = tempPair->x;        // copy contents of referenced random edge
    y        = tempPair->y;        // into local variables
    t        = tempPair->t;

    if (t == LEFT) {
        if (RNG_UNIF01() < 0.5) { // ## LEFT ALPHA move: ((i,j),k) -> ((i,k),j)
            // We need to calculate the change in the likelihood (dLogL)
            // that would result from this move. Most of the information
            // needed to do this is already available, the exception being
            // e_ik, the number of edges that span the i and k subtrees. I
            // use a slow algorithm O(n) to do this, since I don't know of a
            // better way at this point. (After several attempts to find a
            // faster method, no luck.)

            n_i = internal[y].L->n;
            n_j = internal[y].R->n;
            n_k = internal[x].R->n;

            n_y = n_i * n_k;
            e_y = computeEdgeCount(internal[y].L->index, internal[y].L->type,
                                   internal[x].R->index, internal[x].R->type);
            p_y  = (double)(e_y) / (double)(n_y);
            if (e_y == 0 || e_y == n_y) {
                L_y = 0.0;
            } else {
                L_y = (double)(e_y) * log(p_y) + (double)(n_y - e_y) * log(1.0 - p_y);
            }

            n_x  = (n_i + n_k) * n_j;
            e_x  = internal[x].e + internal[y].e - e_y; // e_yj
            p_x  = (double)(e_x) / (double)(n_x);
            if (e_x == 0 || e_x == n_x) {
                L_x = 0.0;
            } else {
                L_x = (double)(e_x) * log(p_x) + (double)(n_x - e_x) * log(1.0 - p_x);
            }

            dLogL = (L_x - internal[x].logL) + (L_y - internal[y].logL);
            if ((dLogL > 0.0) || (RNG_UNIF01() < exp(T * dLogL))) {

                // make LEFT ALPHA move

                ftaken = true;
                d->swapEdges(x, internal[x].R->index, RIGHT, y,
                             internal[y].R->index, RIGHT);
                temp             = internal[x].R; // - swap j and k
                internal[x].R    = internal[y].R;
                internal[y].R    = temp;
                internal[x].R->M = &internal[x];  // - adjust parent pointers
                internal[y].R->M = &internal[y];
                internal[y].n    = n_i + n_k;     // - update n for [y]
                internal[x].e    = e_x;           // - update e_i for [x] and [y]
                internal[y].e    = e_y;
                internal[x].p    = p_x;           // - update p_i for [x] and [y]
                internal[y].p    = p_y;
                internal[x].logL = L_x;           // - update L_i for [x] and [y]
                internal[y].logL = L_y;
                // - order-property maintained
                L  += dLogL;                  // - update LogL
                delta            = dLogL;

            }
        } else {

            // ## LEFT BETA move:  ((i,j),k) -> (i,(j,k))

            n_i = internal[y].L->n;
            n_j = internal[y].R->n;
            n_k = internal[x].R->n;

            n_y  = n_j * n_k;
            e_y  = computeEdgeCount(internal[y].R->index, internal[y].R->type,
                                    internal[x].R->index, internal[x].R->type);
            p_y  = (double)(e_y) / (double)(n_y);
            if (e_y == 0 || e_y == n_y)   {
                L_y = 0.0;
            } else {
                L_y = (double)(e_y) * log(p_y) +
                      (double)(n_y - e_y) * log(1.0 - p_y);
            }

            n_x  = (n_j + n_k) * n_i;
            e_x  = internal[x].e + internal[y].e - e_y; // e_yj
            p_x  = (double)(e_x) / (double)(n_x);
            if (e_x == 0 || e_x == n_x) {
                L_x = 0.0;
            } else {
                L_x = (double)(e_x) * log(p_x) + (double)(n_x - e_x) * log(1.0 - p_x);
            }

            dLogL = (L_x - internal[x].logL) + (L_y - internal[y].logL);
            if ((dLogL > 0.0) || (RNG_UNIF01() < exp(T * dLogL))) {

                // make LEFT BETA move

                ftaken = true;
                d->swapEdges(y, internal[y].L->index, LEFT, y,
                             internal[y].R->index, RIGHT);
                temp   = internal[y].L;       // - swap L and R of [y]
                internal[y].L    = internal[y].R;
                internal[y].R    = temp;
                d->swapEdges(x, internal[x].R->index, RIGHT,
                             y, internal[y].R->index, RIGHT);
                temp   = internal[x].R;       // - swap i and k
                internal[x].R    = internal[y].R;
                internal[y].R    = temp;
                internal[x].R->M = &internal[x];  // - adjust parent pointers
                internal[y].R->M = &internal[y];
                d->swapEdges(x, internal[x].L->index, LEFT,
                             x, internal[x].R->index, RIGHT);
                temp   = internal[x].L;       // - swap L and R of [x]
                internal[x].L    = internal[x].R;
                internal[x].R    = temp;
                internal[y].n    = n_j + n_k;     // - update n
                internal[x].e    = e_x;       // - update e_i
                internal[y].e    = e_y;
                internal[x].p    = p_x;           // - update p_i
                internal[y].p    = p_y;
                internal[x].logL = L_x;           // - update logL_i
                internal[y].logL = L_y;
                if (internal[y].R->label < internal[y].L->label) {
                    // - enforce order-property if necessary
                    d->swapEdges(y, internal[y].L->index, LEFT,
                                 y, internal[y].R->index, RIGHT);
                    temp = internal[y].L;
                    internal[y].L = internal[y].R;
                    internal[y].R = temp;
                } //
                internal[y].label = internal[y].L->label;
                L += dLogL;        // - update LogL
                delta = dLogL;
            }
        }
    } else {

        // right-edge: t == RIGHT

        if (RNG_UNIF01() < 0.5) {

            // alpha move: (i,(j,k)) -> ((i,k),j)

            n_i = internal[x].L->n;
            n_j = internal[y].L->n;
            n_k = internal[y].R->n;

            n_y  = n_i * n_k;
            e_y  = computeEdgeCount(internal[x].L->index, internal[x].L->type,
                                    internal[y].R->index, internal[y].R->type);
            p_y  = (double)(e_y) / (double)(n_y);
            if (e_y == 0 || e_y == n_y)   {
                L_y = 0.0;
            } else {
                L_y = (double)(e_y) * log(p_y) + (double)(n_y - e_y) * log(1.0 - p_y);
            }

            n_x  = (n_i + n_k) * n_j;
            e_x  = internal[x].e + internal[y].e - e_y; // e_yj
            p_x  = (double)(e_x) / (double)(n_x);
            if (e_x == 0 || e_x == n_x) {
                L_x = 0.0;
            } else {
                L_x = (double)(e_x) * log(p_x) + (double)(n_x - e_x) * log(1.0 - p_x);
            }

            dLogL = (L_x - internal[x].logL) + (L_y - internal[y].logL);
            if ((dLogL > 0.0) || (RNG_UNIF01() < exp(T * dLogL))) {

                // make RIGHT ALPHA move

                ftaken = true;
                d->swapEdges(x, internal[x].L->index, LEFT,
                             x, internal[x].R->index, RIGHT);
                temp    = internal[x].L;       // - swap L and R of [x]
                internal[x].L     = internal[x].R;
                internal[x].R     = temp;
                d->swapEdges(y, internal[y].L->index, LEFT,
                             x, internal[x].R->index, RIGHT);
                temp    = internal[y].L;       // - swap i and j
                internal[y].L     = internal[x].R;
                internal[x].R     = temp;
                internal[x].R->M  = &internal[x];  // - adjust parent pointers
                internal[y].L->M  = &internal[y];
                internal[y].n     = n_i + n_k;     // - update n
                internal[x].e     = e_x;       // - update e_i
                internal[y].e     = e_y;
                internal[x].p     = p_x;           // - update p_i
                internal[y].p     = p_y;
                internal[x].logL  = L_x;           // - update logL_i
                internal[y].logL  = L_y;
                internal[y].label = internal[x].label; // - update order property
                L   += dLogL;                  // - update LogL
                delta             = dLogL;
            }
        } else {

            // beta move:  (i,(j,k)) -> ((i,j),k)

            n_i = internal[x].L->n;
            n_j = internal[y].L->n;
            n_k = internal[y].R->n;

            n_y  = n_i * n_j;
            e_y  = computeEdgeCount(internal[x].L->index, internal[x].L->type,
                                    internal[y].L->index, internal[y].L->type);
            p_y  = (double)(e_y) / (double)(n_y);
            if (e_y == 0 || e_y == n_y)   {
                L_y = 0.0;
            } else {
                L_y = (double)(e_y) * log(p_y) + (double)(n_y - e_y) * log(1.0 - p_y);
            }

            n_x  = (n_i + n_j) * n_k;
            e_x  = internal[x].e + internal[y].e - e_y; // e_yk
            p_x  = (double)(e_x) / (double)(n_x);
            if (e_x == 0 || e_x == n_x) {
                L_x = 0.0;
            } else {
                L_x = (double)(e_x) * log(p_x) + (double)(n_x - e_x) * log(1.0 - p_x);
            }

            dLogL = (L_x - internal[x].logL) + (L_y - internal[y].logL);
            if ((dLogL > 0.0) || (RNG_UNIF01() < exp(T * dLogL))) {

                // make RIGHT BETA move

                ftaken = true;
                d->swapEdges(x, internal[x].L->index, LEFT,
                             x, internal[x].R->index, RIGHT);
                temp    = internal[x].L;       // - swap L and R of [x]
                internal[x].L     = internal[x].R;
                internal[x].R     = temp;
                d->swapEdges(x, internal[x].R->index, RIGHT,
                             y, internal[y].R->index, RIGHT);
                temp    = internal[x].R;       // - swap i and k
                internal[x].R     = internal[y].R;
                internal[y].R     = temp;
                internal[x].R->M  = &internal[x];  // - adjust parent pointers
                internal[y].R->M  = &internal[y];
                d->swapEdges(y, internal[y].L->index, LEFT,
                             y, internal[y].R->index, RIGHT);
                temp    = internal[y].L;       // - swap L and R of [y]
                internal[y].L     = internal[y].R;
                internal[y].R     = temp;
                internal[y].n     = n_i + n_j;     // - update n
                internal[x].e     = e_x;       // - update e_i
                internal[y].e     = e_y;
                internal[x].p     = p_x;       // - update p_i
                internal[y].p     = p_y;
                internal[x].logL  = L_x;       // - update logL_i
                internal[y].logL  = L_y;
                internal[y].label = internal[x].label; // - order-property
                L   += dLogL;                  // - update LogL
                delta             = dLogL;
            }
        }
    }
    return true;
}

// **********************************************************************

void dendro::refreshLikelihood() {
    // recalculates the log-likelihood of the dendrogram structure
    L = 0.0; double dL;
    int nL_nR, ei;
    for (int i = 0; i < (n - 1); i++) {
        nL_nR = internal[i].L->n * internal[i].R->n;
        ei    = internal[i].e;
        internal[i].p = (double)(ei) / (double)(nL_nR);
        if (ei == 0 || ei == nL_nR) {
            dL = 0.0;
        } else {
            dL = ei * log(internal[i].p) + (nL_nR - ei) * log(1.0 - internal[i].p);
        }
        internal[i].logL = dL;
        L += dL;
    }
    return;
}

// **********************************************************************

void dendro::QsortMain (block* array, int left, int right) {
    if (right > left) {
        int pivot = left;
        int part  = QsortPartition(array, left, right, pivot);
        QsortMain(array, left,   part - 1);
        QsortMain(array, part + 1, right  );
    }
    return;
}

int dendro::QsortPartition (block* array, int left, int right, int index) {
    block p_value, temp;
    p_value.x = array[index].x;
    p_value.y = array[index].y;

    // swap(array[p_value], array[right])
    temp.x = array[right].x;
    temp.y = array[right].y;
    array[right].x = array[index].x;
    array[right].y = array[index].y;
    array[index].x = temp.x;
    array[index].y = temp.y;

    int stored = left;
    for (int i = left; i < right; i++) {
        if (array[i].x <= p_value.x) {
            // swap(array[stored], array[i])
            temp.x = array[i].x;
            temp.y = array[i].y;
            array[i].x = array[stored].x;
            array[i].y = array[stored].y;
            array[stored].x = temp.x;
            array[stored].y = temp.y;
            stored++;
        }
    }
    // swap(array[right], array[stored])
    temp.x = array[stored].x;
    temp.y = array[stored].y;
    array[stored].x = array[right].x;
    array[stored].y = array[right].y;
    array[right].x = temp.x;
    array[right].y  = temp.y;

    return stored;
}

void dendro::recordConsensusTree(igraph_vector_t *parents,
                                 igraph_vector_t *weights) {

    keyValuePairSplit *curr, *prev;
    child *newChild;
    int orig_nodes = g->numNodes();

    // First, cull the split hist so that only splits with weight >= 0.5
    // remain
    cullSplitHist();
    int treesize = splithist->returnNodecount();

    // Now, initialize the various arrays we use to keep track of the
    // internal structure of the consensus tree.
    ctree  = new cnode[treesize];
    cancestor = new int[n];
    for (int i = 0; i < treesize; i++) {
        ctree[i].index = i;
    }
    for (int i = 0; i < n; i++)        {
        cancestor[i]   = -1;
    }
    int ii = 0;

    // To build the majority consensus tree, we do the following: For
    // each possible number of Ms in the split string (a number that
    // ranges from n-2 down to 0), and for each split with that number
    // of Ms, we create a new internal node of the tree, and connect the
    // oldest ancestor of each C to that node (at most once). Then, we
    // update our list of oldest ancestors to reflect this new join, and
    // proceed.
    for (int i = n - 2; i >= 0; i--) {
        // First, we get a list of all the splits with this exactly i Ms
        curr = splithist->returnTheseSplits(i);

        // Now we loop over that list
        while (curr != NULL) {
            splithist->deleteItem(curr->x);
            // add weight to this internal node
            ctree[ii].weight = curr->y;
            // examine each letter of this split
            for (int j = 0; j < n; j++) {
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
            prev = curr;
            curr = curr->next;
            delete prev;
        }
    }

    // Return the consensus tree
    igraph_vector_resize(parents, ii + orig_nodes);
    if (weights) {
        igraph_vector_resize(weights, ii);
    }

    for (int i = 0; i < ii; i++) {
        child *sat, *sit = ctree[i].children;
        while (sit) {
            VECTOR(*parents)[orig_nodes + i] =
                ctree[i].parent < 0 ? -1 : orig_nodes + ctree[i].parent;
            if (sit->type == GRAPH) {
                VECTOR(*parents)[sit->index] = orig_nodes + i;
            }
            sat = sit;
            sit = sit->next;
            delete sat;
        }
        if (weights) {
            VECTOR(*weights)[i] = ctree[i].weight;
        }
        ctree[i].children = 0;
    }

    // Plus the isolate nodes
    for (int i = 0; i < n; i++) {
        if (cancestor[i] == -1) {
            VECTOR(*parents)[i] = -1;
        }
    }


}

// **********************************************************************

void dendro::recordDendrogramStructure(igraph_hrg_t *hrg) {
    for (int i = 0; i < n - 1; i++) {
        int li = internal[i].L->index;
        int ri = internal[i].R->index;
        VECTOR(hrg->left )[i] = internal[i].L->type == DENDRO ? -li - 1 : li;
        VECTOR(hrg->right)[i] = internal[i].R->type == DENDRO ? -ri - 1 : ri;
        VECTOR(hrg->prob )[i] = internal[i].p;
        VECTOR(hrg->edges)[i] = internal[i].e;
        VECTOR(hrg->vertices)[i] = internal[i].n;
    }
}

void dendro::recordGraphStructure(igraph_t *graph) {
    igraph_vector_t edges;
    int no_of_nodes = g->numNodes();
    int no_of_edges = g->numLinks() / 2;
    int idx = 0;

    igraph_vector_init(&edges, no_of_edges * 2);
    IGRAPH_FINALLY(igraph_vector_destroy, &edges);

    for (int i = 0; i < n; i++) {
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

// **********************************************************************

list* dendro::reversePathToRoot(const int leafIndex) {
    list *head, *subhead, *newlist;
    head = subhead = newlist = NULL;
    elementd *current = &leaf[leafIndex];

    // continue until we're finished
    while (current != NULL) {
        // add this node to the path
        newlist = new list;
        newlist->x = current->index;
        newlist->next = NULL;
        if (head == NULL) {
            head    = newlist;
        } else {
            subhead = head;
            head = newlist;
            head->next = subhead;
        }
        current = current->M;
    }
    return head;
}

// ***********************************************************************

bool dendro::sampleSplitLikelihoods(int &sample_num) {
    // In order to compute the majority agreement dendrogram at
    // equilibrium, we need to calculate the leaf partition defined by
    // each split (internal edge) of the tree. Because splits are only
    // defined on a Cayley tree, the buildSplit() function returns the
    // default "--...--"  string for the root and the root's left
    // child. When tabulating the frequency of splits, one of these
    // needs to be excluded.

    IGRAPH_UNUSED(sample_num);

    string* array;
    int     k;
    double  tot;

    string new_split;
    // To decompose the tree into its splits, we simply loop over all
    // the internal nodes and replace the old split for the ith internal
    // node with its new split. This is a bit time consuming to do
    // O(n^2), so try not to do this very often. Once the decomposition
    // is had, we insert them into the split histogram, which tracks the
    // cumulative weight for each respective split observed.

    if (splithist == NULL) {
        splithist = new splittree;
    }
    for (int i = 0; i < (n - 1); i++) {
        new_split = buildSplit(&internal[i]);
        d->replaceSplit(i, new_split);
        if (!new_split.empty() && new_split[1] != '-') {
            if (!splithist->insertItem(new_split, 1.0)) {
                return false;
            }
        }
    }
    splithist->finishedThisRound();

    // For large graphs, the split histogram can get extremely large, so
    // we need to employ some measures to prevent it from swamping the
    // available memory. When the number of splits exceeds  a threshold
    // (say, a million), we progressively delete splits that have a
    // weight less than  a rising (k*0.001 of the total weight) fraction
    // of the splits, on the assumption that losing such weight is
    // unlikely to effect the ultimate split statistics. This deletion
    // procedure is slow O(m lg m), but should only happen very rarely.

    int split_max = n * 500;
    int leng;
    if (splithist->returnNodecount() > split_max) {
        k = 1;
        while (splithist->returnNodecount() > split_max) {
            array = splithist->returnArrayOfKeys();
            tot   = splithist->returnTotal();
            leng  = splithist->returnNodecount();
            for (int i = 0; i < leng; i++) {
                if ((splithist->returnValue(array[i]) / tot) < k * 0.001) {
                    splithist->deleteItem(array[i]);
                }
            }
            delete [] array; array = NULL;
            k++;
        }
    }

    return true;
}

void dendro::sampleAdjacencyLikelihoods() {
    // Here, we sample the probability values associated with every
    // adjacency in A, weighted by their likelihood. The weighted
    // histogram is stored in the graph data structure, so we simply
    // need to add an observation to each node-pair that corresponds to
    // the associated branch point's probability and the dendrogram's
    // overall likelihood.

    double nn;
    double norm = ((double)(n) * (double)(n)) / 4.0;

    if (L > 0.0) {
        L = 0.0;
    }
    elementd* ancestor;
    list *currL, *prevL;
    if (paths != NULL) {
        for (int i = 0; i < n; i++) {
            currL = paths[i];
            while (currL != NULL) {
                prevL = currL;
                currL = currL->next;
                delete prevL;
                prevL = NULL;
            }
            paths[i] = NULL;
        }
        delete [] paths;
    }
    paths = NULL;
    paths = new list* [n];
    for (int i = 0; i < n; i++) {
        // construct paths from root, O(n^2) at worst
        paths[i] = reversePathToRoot(i);
    }

    // add obs for every node-pair, always O(n^2)
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            // find internal node, O(n) at worst
            ancestor = findCommonAncestor(paths, i, j);
            nn = ((double)(ancestor->L->n) * (double)(ancestor->R->n)) / norm;
            // add obs of ->p to (i,j) histogram, and
            g->addAdjacencyObs(i, j, ancestor->p, nn);
            // add obs of ->p to (j,i) histogram
            g->addAdjacencyObs(j, i, ancestor->p, nn);
        }
    }

    // finish-up: upate total weight in histograms
    g->addAdjacencyEnd();

    return;
}

void dendro::resetDendrograph() {
    // Reset the dendrograph structure for the next trial
    if (leaf      != NULL) {
        delete [] leaf;        // O(n)
        leaf      = NULL;
    }
    if (internal  != NULL) {
        delete [] internal;    // O(n)
        internal  = NULL;
    }
    if (d         != NULL) {
        delete d;              // O(n)
        d         = NULL;
    }
    root = NULL;
    if (paths != NULL) {
        list *curr, *prev;
        for (int i = 0; i < n; i++) {
            curr = paths[i];
            while (curr != NULL) {
                prev = curr;
                curr = curr->next;
                delete prev;
                prev = NULL;
            }
            paths[i] = NULL;
        }
        delete [] paths;
    }
    paths = NULL;
    L = 1.0;

    return;
}

// **********************************************************************
// *** COPYRIGHT NOTICE *************************************************
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
// **********************************************************************
// Author       : Aaron Clauset  ( aaronc@santafe.edu |
//                                 http://www.santafe.edu/~aaronc/ )
// Collaborators: Cristopher Moore and Mark E.J. Newman
// Project      : Hierarchical Random Graphs
// Location     : University of New Mexico, Dept. of Computer Science
//                AND Santa Fe Institute
// Created      : 8 November 2005
// Modified     : 23 December 2007 (cleaned up for public consumption)
//
// ***********************************************************************
//
// Graph data structure for hierarchical random graphs. The basic
// structure is an adjacency list of edges; however, many additional
// pieces of metadata are stored as well. Each node stores its
// external name, its degree and (if assigned) its group index.
//
// ***********************************************************************

// ******** Constructor / Destructor *************************************

graph::graph(const int size, bool predict) : predict(predict)  {
    n = size;
    m = 0;
    nodes = new vert  [n];
    nodeLink = new edge* [n];
    nodeLinkTail   = new edge* [n];
    for (int i = 0; i < n; i++) {
        nodeLink[i] = NULL;
        nodeLinkTail[i] = NULL;
    }
    if (predict) {
        A = new double** [n];
        for (int i = 0; i < n; i++) {
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
    for (int i = 0; i < n; i++) {
        curr = nodeLink[i];
        while (curr != NULL) {
            prev = curr;
            curr = curr->next;
            delete prev;
        }
    }
    delete [] nodeLink; nodeLink = NULL;
    delete [] nodeLinkTail;  nodeLinkTail   = NULL;
    delete [] nodes; nodes = NULL;

    if (predict) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                delete [] A[i][j];
            }
            delete [] A[i];
        }
        delete [] A; A = NULL;
    }
}

// **********************************************************************

bool graph::addLink(const int i, const int j) {
    // Adds the directed edge (i,j) to the adjacency list for v_i
    edge* newedge;
    if (i >= 0 && i < n && j >= 0 && j < n) {
        newedge  = new edge;
        newedge->x = j;
        if (nodeLink[i] == NULL) {
            // first neighbor
            nodeLink[i]  = newedge;
            nodeLinkTail[i] = newedge;
            nodes[i].degree = 1;
        } else {
            // subsequent neighbor
            nodeLinkTail[i]->next = newedge;
            nodeLinkTail[i]       = newedge;
            nodes[i].degree++;
        }
        // increment edge count
        m++;
        return true;
    } else {
        return false;
    }
}

// ***********************************************************************

bool graph::addAdjacencyObs(const int i, const int j,
                            const double probability, const double size) {
    // Adds the observation obs to the histogram of the edge (i,j)
    // Note: user must manually add observation to edge (j,i) by calling
    // this function with that argument
    if (bin_resolution > 0.0 && probability >= 0.0 && probability <= 1.0
        && size >= 0.0 && size <= 1.0
        && i >= 0 && i < n && j >= 0 && j < n) {
        int index = (int)(probability / bin_resolution + 0.5);
        if (index < 0) {
            index = 0;
        } else if (index > num_bins) {
            index = num_bins;
        }

        // Add the weight to the proper probability bin
        if (A[i][j][index] < 0.5) {
            A[i][j][index] = 1.0;
        } else {
            A[i][j][index] += 1.0;
        }
        return true;
    }
    return false;
}

// **********************************************************************

void graph::addAdjacencyEnd() {
    // We need to also keep a running total of how much weight has been added
    // to the histogram, and the number of observations in the histogram.
    if (obs_count == 0) {
        total_weight  = 1.0; obs_count = 1;
    } else {
        total_weight += 1.0; obs_count++;
    }
    return;
}

bool graph::doesLinkExist(const int i, const int j) {
    // This function determines if the edge (i,j) already exists in the
    // adjacency list of v_i
    edge* curr;
    if (i >= 0 && i < n && j >= 0 && j < n) {
        curr = nodeLink[i];
        while (curr != NULL) {
            if (curr->x == j) {
                return true;
            }
            curr = curr->next;
        }
    }
    return false;
}

// **********************************************************************

int graph::getDegree(const int i) {
    if (i >= 0 && i < n) {
        return nodes[i].degree;
    } else {
        return -1;
    }
}

string graph::getName(const int i)  {
    if (i >= 0 && i < n) {
        return nodes[i].name;
    } else {
        return "";
    }
}

// NOTE: Returns address; deallocation of returned object is dangerous
edge* graph::getNeighborList(const int i) {
    if (i >= 0 && i < n) {
        return nodeLink[i];
    } else {
        return NULL;
    }
}

double* graph::getAdjacencyHist(const int i, const int j) {
    if (i >= 0 && i < n && j >= 0 && j < n) {
        return A[i][j];
    } else {
        return NULL;
    }
}

// **********************************************************************

double graph::getAdjacencyAverage(const int i, const int j) {
    double average = 0.0;
    if (i != j) {
        for (int k = 0; k < num_bins; k++) {
            if (A[i][j][k] > 0.0) {
                average += (A[i][j][k] / total_weight) * ((double)(k) * bin_resolution);
            }
        }
    }
    return average;
}

int graph::numLinks() {
    return m;
}

int graph::numNodes() {
    return n;
}

double graph::getBinResolution() {
    return bin_resolution;
}

int graph::getNumBins() {
    return num_bins;
}

double graph::getTotalWeight() {
    return total_weight;
}

// ***********************************************************************

void graph::resetAllAdjacencies() {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < num_bins; k++) {
                A[i][j][k] = 0.0;
            }
        }
    }
    obs_count    = 0;
    total_weight = 0.0;
    return;
}

// **********************************************************************

void graph::resetAdjacencyHistogram(const int i, const int j) {
    if (i >= 0 && i < n && j >= 0 && j < n) {
        for (int k = 0; k < num_bins; k++) {
            A[i][j][k] = 0.0;
        }
    }
    return;
}

// **********************************************************************

void graph::resetLinks() {
    edge *curr, *prev;
    for (int i = 0; i < n; i++) {
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

// **********************************************************************

void graph::setAdjacencyHistograms(const int bin_count) {
    // For all possible adjacencies, setup an edge histograms
    num_bins = bin_count + 1;
    bin_resolution = 1.0 / (double)(bin_count);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A[i][j] = new double [num_bins];
            for (int k = 0; k < num_bins; k++) {
                A[i][j][k] = 0.0;
            }
        }
    }
    return;
}

bool graph::setName(const int i, const string text) {
    if (i >= 0 && i < n) {
        nodes[i].name = text;
        return true;
    } else {
        return false;
    }
}

// **********************************************************************

interns::interns(const int n)  {
    q         = n;
    count     = 0;
    edgelist  = new ipair  [q];
    splitlist = new string [q + 1];
    indexLUT  = new int*   [q + 1];
    for (int i = 0; i < (q + 1); i++) {
        indexLUT[i]    = new int [2];
        indexLUT[i][0] = indexLUT[i][1] = -1;
    }
}
interns::~interns() {
    delete [] edgelist;
    delete [] splitlist;
    for (int i = 0; i < (q + 1); i++) {
        delete [] indexLUT[i];
    }
    delete [] indexLUT;
}

// ***********************************************************************

// NOTE: Returns an address to another object -- do not deallocate
ipair* interns::getEdge(const int i) {
    return &edgelist[i];
}

// ***********************************************************************

// NOTE: Returns an address to another object -- do not deallocate
ipair* interns::getRandomEdge() {
    return &edgelist[(int)(floor((double)(q) * RNG_UNIF01()))];
}

// ***********************************************************************

string interns::getSplit(const int i) {
    if (i >= 0 && i <= q) {
        return splitlist[i];
    } else {
        return "";
    }
}

// **********************************************************************

bool interns::addEdge(const int new_x, const int new_y,
                      const short int new_type) {
    // This function adds a new edge (i,j,t,sp) to the list of internal
    // edges. After checking that the inputs fall in the appropriate
    // range of values, it records the new edgelist index in the
    // indexLUT and then puts the input values into that edgelist
    // location.

    if (count < q && new_x >= 0 && new_x < (q + 1) && new_y >= 0 &&
        new_y < (q + 2) && (new_type == LEFT || new_type == RIGHT)) {
        if (new_type == LEFT) {
            indexLUT[new_x][0] = count;
        } else {
            indexLUT[new_x][1] = count;
        }
        edgelist[count].x = new_x;
        edgelist[count].y = new_y;
        edgelist[count].t = new_type;
        count++;
        return true;
    } else {
        return false;
    }
}

// **********************************************************************

bool interns::replaceSplit(const int i, const string sp) {
    // When an internal edge is changed, its split must be replaced as
    // well. This function provides that access; it stores the split
    // defined by an internal edge (x,y) at the location [y], which
    // is unique.

    if (i >= 0 && i <= q) {
        splitlist[i] = sp;
        return true;
    }
    return false;
}

// ***********************************************************************

bool interns::swapEdges(const int one_x, const int one_y,
                        const short int one_type, const int two_x,
                        const int two_y, const short int two_type) {
    // The moves on the dendrogram always swap edges, either of which
    // (or both, or neither) can by internal edges. So, this function
    // mirrors that operation for the internal edgelist and indexLUT.

    int index, jndex, temp;
    bool one_isInternal = false;
    bool two_isInternal = false;

    if (one_x >= 0 && one_x < (q + 1) && two_x >= 0 && two_x < (q + 1) &&
        (two_type == LEFT || two_type == RIGHT) &&
        one_y >= 0 && one_y < (q + 2) && two_y >= 0 &&
        two_y < (q + 2) && (one_type == LEFT || one_type == RIGHT)) {

        if (one_type == LEFT) {
            temp = 0;
        } else {
            temp = 1;
        }
        if (indexLUT[one_x][temp] > -1) {
            one_isInternal = true;
        }
        if (two_type == LEFT) {
            temp = 0;
        } else {
            temp = 1;
        }
        if (indexLUT[two_x][temp] > -1) {
            two_isInternal = true;
        }

        if (one_isInternal && two_isInternal) {
            if (one_type == LEFT)  {
                index = indexLUT[one_x][0];
            } else {
                index = indexLUT[one_x][1];
            }
            if (two_type == LEFT)  {
                jndex = indexLUT[two_x][0];
            } else {
                jndex = indexLUT[two_x][1];
            }
            temp              = edgelist[index].y;
            edgelist[index].y = edgelist[jndex].y;
            edgelist[jndex].y = temp;

        } else if (one_isInternal) {
            if (one_type == LEFT)  {
                index = indexLUT[one_x][0]; indexLUT[one_x][0] = -1;
            } else {
                index = indexLUT[one_x][1]; indexLUT[one_x][1] = -1;
            }
            edgelist[index].x = two_x;
            edgelist[index].t = two_type;
            if (two_type == LEFT) {
                indexLUT[two_x][0] = index;
            } else {
                indexLUT[two_x][1] = index;
            } // add new

        } else if (two_isInternal) {
            if (two_type == LEFT)  {
                index = indexLUT[two_x][0]; indexLUT[two_x][0] = -1;
            } else {
                index = indexLUT[two_x][1]; indexLUT[two_x][1] = -1;
            }
            edgelist[index].x = one_x;
            edgelist[index].t = one_type;
            if (one_type == LEFT) {
                indexLUT[one_x][0] = index;
            } else {
                indexLUT[one_x][1] = index;
            } // add new
        } else {
            ;
        } // else neither is internal

        return true;
    } else {
        return false;
    }
}

// ******** Red-Black Tree Methods ***************************************

splittree::splittree() {
    root = new elementsp;
    leaf = new elementsp;

    leaf->parent   = root;

    root->left = leaf;
    root->right    = leaf;
    support = 0;
    total_weight = 0.0;
    total_count = 0;
}

splittree::~splittree() {
    if (root != NULL && (root->left != leaf || root->right != leaf)) {
        deleteSubTree(root); root = NULL;
    }
    support      = 0;
    total_weight = 0.0;
    total_count  = 0;
    if (root) {
        delete root;
    }
    delete leaf;
    root    = NULL;
    leaf    = NULL;
}

void splittree::deleteTree() {
    if (root != NULL) {
        deleteSubTree(root);
        root = NULL;
    }
    return;
}

void splittree::deleteSubTree(elementsp *z) {
    if (z->left  != leaf) {
        deleteSubTree(z->left);
        z->left = NULL;
    }
    if (z->right != leaf) {
        deleteSubTree(z->right);
        z->right = NULL;
    }
    delete z;
    /* No point in setting z to NULL here because z is passed by value */
    /* z = NULL; */
    return;
}

// ******** Reset Functions *********************************************

// O(n lg n)
void splittree::clearTree() {
    string *array = returnArrayOfKeys();
    for (int i = 0; i < support; i++) {
        deleteItem(array[i]);
    }
    delete [] array;
    return;
}

// ******** Search Functions *********************************************
// public search function - if there exists a elementsp in the tree
// with key=searchKey, it returns TRUE and foundNode is set to point
// to the found node; otherwise, it sets foundNode=NULL and returns
// FALSE
elementsp* splittree::findItem(const string searchKey) {

    elementsp *current = root;
    if (current->split.empty()) {
        return NULL;    // empty tree; bail out
    }
    while (current != leaf) {
        if (searchKey.compare(current->split) < 0) { // left-or-right?
            // try moving down-left
            if (current->left  != leaf) {
                current = current->left;
            } else {
                // failure; bail out
                return NULL;
            }
        } else {
            if (searchKey.compare(current->split) > 0) {
                // left-or-right?
                if (current->right != leaf) {
                    // try moving down-left
                    current = current->right;
                } else {
                    //   failure; bail out
                    return NULL;
                }
            } else {
                // found (searchKey==current->split)
                return current;
            }
        }
    }
    return NULL;
}

double splittree::returnValue(const string searchKey) {
    elementsp* test = findItem(searchKey);
    if (test == NULL) {
        return 0.0;
    } else {
        return test->weight;
    }
}


// ******** Return Item Functions ***************************************
// public function which returns the tree, via pre-order traversal, as
// a linked list

string* splittree::returnArrayOfKeys() {
    string* array;
    array = new string [support];
    bool flag_go = true;
    int index = 0;
    elementsp *curr;

    if (support == 1) {
        array[0] = root->split;
    } else if (support == 2) {
        array[0] = root->split;
        if (root->left == leaf) {
            array[1] = root->right->split;
        } else {
            array[1] = root->left->split;
        }
    } else {
        for (int i = 0; i < support; i++) {
            array[i] = -1;
        }
        // non-recursive traversal of tree structure
        curr  = root;
        curr->mark = 1;
        while (flag_go) {

            // - is it time, and is left child the leaf node?
            if (curr->mark == 1 && curr->left == leaf) {
                curr->mark = 2;
            }
            // - is it time, and is right child the leaf node?
            if (curr->mark == 2 && curr->right == leaf) {
                curr->mark = 3;
            }
            if (curr->mark == 1) {               // - go left
                curr->mark = 2;
                curr = curr->left;
                curr->mark = 1;
            } else if (curr->mark == 2) {        // - else go right
                curr->mark = 3;
                curr = curr->right;
                curr->mark = 1;
            } else {                     // - else go up a level
                curr->mark = 0;
                array[index++] = curr->split;
                curr = curr->parent;
                if (curr == NULL) {
                    flag_go = false;
                }
            }
        }
    }

    return array;
}

slist* splittree::returnListOfKeys() {
    keyValuePairSplit *curr, *prev;
    slist *head = NULL, *tail = NULL, *newlist;

    curr = returnTreeAsList();
    while (curr != NULL) {
        newlist = new slist;
        newlist->x = curr->x;
        if (head == NULL) {
            head = newlist; tail = head;
        } else {
            tail->next = newlist; tail = newlist;
        }
        prev = curr;
        curr = curr->next;
        delete prev;
        prev = NULL;
    }
    return head;
}

// pre-order traversal
keyValuePairSplit* splittree::returnTreeAsList() {
    keyValuePairSplit  *head, *tail;

    head    = new keyValuePairSplit;
    head->x = root->split;
    head->y = root->weight;
    head->c = root->count;
    tail    = head;

    if (root->left  != leaf) {
        tail = returnSubtreeAsList(root->left,  tail);
    }
    if (root->right != leaf) {
        tail = returnSubtreeAsList(root->right, tail);
    }

    if (head->x.empty()) {
        return NULL; /* empty tree */
    } else {
        return head;
    }
}

keyValuePairSplit* splittree::returnSubtreeAsList(elementsp *z,
        keyValuePairSplit *head) {
    keyValuePairSplit *newnode, *tail;

    newnode    = new keyValuePairSplit;
    newnode->x = z->split;
    newnode->y = z->weight;
    newnode->c = z->count;
    head->next = newnode;
    tail       = newnode;

    if (z->left  != leaf) {
        tail = returnSubtreeAsList(z->left,  tail);
    }
    if (z->right != leaf) {
        tail = returnSubtreeAsList(z->right, tail);
    }

    return tail;
}

keyValuePairSplit splittree::returnMaxKey() {
    keyValuePairSplit themax;
    elementsp *current;
    current = root;
    // search to bottom-right corner of tree
    while (current->right != leaf) {
        current = current->right;
    }
    themax.x = current->split;
    themax.y = current->weight;

    return themax;
}

keyValuePairSplit splittree::returnMinKey() {
    keyValuePairSplit themin;
    elementsp *current;
    current = root;
    // search to bottom-left corner of tree
    while (current->left != leaf) {
        current = current->left;
    }
    themin.x = current->split;
    themin.y = current->weight;

    return themin;
}

// private functions for deleteItem() (although these could easily be
// made public, I suppose)
elementsp* splittree::returnMinKey(elementsp *z) {
    elementsp *current;

    current = z;
    // search to bottom-right corner of tree
    while (current->left != leaf) {
        current = current->left;
    }
    // return pointer to the minimum
    return current;
}

elementsp* splittree::returnSuccessor(elementsp *z) {
    elementsp *current, *w;

    w = z;
// if right-subtree exists, return min of it
    if (w->right != leaf) {
        return returnMinKey(w->right);
    }
    // else search up in tree
    // move up in tree until find a non-right-child
    current = w->parent;
    while ((current != NULL) && (w == current->right)) {
        w = current;
        current = current->parent;
    }
    return current;
}

int splittree::returnNodecount() {
    return support;
}

keyValuePairSplit* splittree::returnTheseSplits(const int target) {
    keyValuePairSplit *head, *curr, *prev, *newhead, *newtail, *newpair;
    int count, len;

    head = returnTreeAsList();
    prev = newhead = newtail = newpair = NULL;
    curr = head;

    while (curr != NULL) {
        count = 0;
        len   = curr->x.size();
        for (int i = 0; i < len; i++) {
            if (curr->x[i] == 'M') {
                count++;
            }
        }
        if (count == target && curr->x[1] != '*') {
            newpair       = new keyValuePairSplit;
            newpair->x    = curr->x;
            newpair->y    = curr->y;
            newpair->next = NULL;
            if (newhead == NULL) {
                newhead = newpair; newtail = newpair;
            } else {
                newtail->next = newpair; newtail = newpair;
            }
        }
        prev = curr;
        curr = curr->next;
        delete prev;
        prev = NULL;
    }

    return newhead;
}

double splittree::returnTotal() {
    return total_weight;
}

// ******** Insert Functions *********************************************

void splittree::finishedThisRound() {
    // We need to also keep a running total of how much weight has been
    // added to the histogram.
    if (total_count == 0) {
        total_weight  = 1.0; total_count = 1;
    } else {
        total_weight += 1.0; total_count++;
    }
    return;
}

// public insert function
bool splittree::insertItem(string newKey, double newValue) {

    // first we check to see if newKey is already present in the tree;
    // if so, we do nothing; if not, we must find where to insert the
    // key
    elementsp *newNode, *current;

// find newKey in tree; return pointer to it O(log k)
    current = findItem(newKey);
    if (current != NULL) {
        current->weight += 1.0;
        // And finally, we keep track of how many observations went into
        // the histogram
        current->count++;
        return true;
    } else {
        newNode = new elementsp;    // elementsp for the splittree
        newNode->split = newKey;    //  store newKey
        newNode->weight = newValue; //  store newValue
        newNode->color = true;  //  new nodes are always RED
        newNode->parent = NULL; //  new node initially has no parent
        newNode->left = leaf;   //  left leaf
        newNode->right = leaf;  //  right leaf
        newNode->count = 1;
        support++;          // increment node count in splittree

        // must now search for where to insert newNode, i.e., find the
        // correct parent and set the parent and child to point to each
        // other properly
        current = root;
        if (current->split.empty()) {   // insert as root
            delete root;      //   delete old root
            root = newNode;       //   set root to newNode
            leaf->parent   = newNode; //   set leaf's parent
            current = leaf;       //   skip next loop
        }

        // search for insertion point
        while (current != leaf) {
            // left-or-right?
            if (newKey.compare(current->split) < 0) {
                // try moving down-left
                if (current->left  != leaf) {
                    current = current->left;
                } else {
                    // else found new parent
                    newNode->parent = current; // set parent
                    current->left = newNode;   // set child
                    current = leaf;        // exit search
                }
            } else { //
                if (current->right != leaf) {
                    // try moving down-right
                    current = current->right;
                } else {
                    // else found new parent
                    newNode->parent = current; // set parent
                    current->right = newNode;  // set child
                    current = leaf;        // exit search
                }
            }
        }

        // now do the house-keeping necessary to preserve the red-black
        // properties
        insertCleanup(newNode);

    }
    return true;
}

// private house-keeping function for insertion
void splittree::insertCleanup(elementsp *z) {

    // fix now if z is root
    if (z->parent == NULL) {
        z->color = false; return;
    }
    elementsp *temp;
    // while z is not root and z's parent is RED
    while (z->parent != NULL && z->parent->color) {
        if (z->parent == z->parent->parent->left) {  // z's parent is LEFT-CHILD
            temp = z->parent->parent->right;       // grab z's uncle
            if (temp->color) {
                z->parent->color = false;          // color z's parent BLACK (Case 1)
                temp->color = false;               // color z's uncle BLACK  (Case 1)
                z->parent->parent->color = true;   // color z's grandpa  RED (Case 1)
                z = z->parent->parent;             // set z = z's grandpa    (Case 1)
            } else {
                if (z == z->parent->right) {       // z is RIGHT-CHILD
                    z = z->parent;                   // set z = z's parent     (Case 2)
                    rotateLeft(z);                   // perform left-rotation  (Case 2)
                }
                z->parent->color = false;          // color z's parent BLACK (Case 3)
                z->parent->parent->color = true;   // color z's grandpa RED  (Case 3)
                rotateRight(z->parent->parent);    // perform right-rotation (Case 3)
            }
        } else {                       // z's parent is RIGHT-CHILD
            temp = z->parent->parent->left;      // grab z's uncle
            if (temp->color) {
                z->parent->color = false;          // color z's parent BLACK (Case 1)
                temp->color = false;               // color z's uncle BLACK  (Case 1)
                z->parent->parent->color = true;   // color z's grandpa RED  (Case 1)
                z = z->parent->parent;             // set z = z's grandpa    (Case 1)
            } else {
                if (z == z->parent->left) {        // z is LEFT-CHILD
                    z = z->parent;                   // set z = z's parent     (Case 2)
                    rotateRight(z);                  // perform right-rotation (Case 2)
                }
                z->parent->color = false;          // color z's parent BLACK (Case 3)
                z->parent->parent->color = true;   // color z's grandpa RED  (Case 3)
                rotateLeft(z->parent->parent);     // perform left-rotation  (Case 3)
            }
        }
    }

    root->color = false; // color the root BLACK
    return;
}

// ******** Delete Functions ********************************************
// public delete function
void splittree::deleteItem(string killKey) {
    elementsp *x, *y, *z;

    z = findItem(killKey);
    if (z == NULL) {
        return;    // item not present; bail out
    }

    if (support == 1) {   // -- attempt to delete the root
        root->split = "";       // restore root node to default state
        root->weight = 0.0;     //
        root->color = false;    //
        root->parent = NULL;    //
        root->left = leaf;      //
        root->right = leaf;     //
        support--;          // set support to zero
        total_weight = 0.0;     // set total weight to zero
        total_count--;      //
        return;         // exit - no more work to do
    }

    if (z != NULL) {
        support--;          // decrement node count
        if ((z->left == leaf) || (z->right == leaf)) {
            // case of less than two children
            y = z;             // set y to be z
        } else {
            y = returnSuccessor(z);    // set y to be z's key-successor
        }

        if (y->left != leaf) {
            x = y->left;       // pick y's one child (left-child)
        } else {
            x = y->right;              // (right-child)
        }
        x->parent = y->parent;       // make y's child's parent be y's parent

        if (y->parent == NULL) {
            root = x;          // if y is the root, x is now root
        } else {
            if (y == y->parent->left) {// decide y's relationship with y's parent
                y->parent->left  = x;    // replace x as y's parent's left child
            } else {
                y->parent->right = x;
            }  // replace x as y's parent's left child
        }

        if (y != z) {        // insert y into z's spot
            z->split = y->split;   // copy y data into z
            z->weight = y->weight;     //
            z->count = y->count;   //
        }                //

        // do house-keeping to maintain balance
        if (y->color == false) {
            deleteCleanup(x);
        }
        delete y;            // deallocate y
        y = NULL;            // point y to NULL for safety
    }              //

    return;
}

void splittree::deleteCleanup(elementsp *x) {
    elementsp *w, *t;
    // until x is the root, or x is RED
    while ((x != root) && (x->color == false)) {
        if (x == x->parent->left) {      // branch on x being a LEFT-CHILD
            w = x->parent->right;      // grab x's sibling
            if (w->color == true) {    // if x's sibling is RED
                w->color = false;        // color w BLACK                (case 1)
                x->parent->color = true;     // color x's parent RED         (case 1)
                rotateLeft(x->parent);       // left rotation on x's parent  (case 1)
                w = x->parent->right;        // make w be x's right sibling  (case 1)
            }
            if ((w->left->color == false) && (w->right->color == false)) {
                w->color = true;         // color w RED                  (case 2)
                x = x->parent;           // examine x's parent           (case 2)
            } else {               //
                if (w->right->color == false) {
                    w->left->color = false;    // color w's left child BLACK   (case 3)
                    w->color = true;       // color w RED                  (case 3)
                    t = x->parent;         // store x's parent
                    rotateRight(w);        // right rotation on w          (case 3)
                    x->parent = t;         // restore x's parent
                    w = x->parent->right;      // make w be x's right sibling  (case 3)
                } //
                w->color = x->parent->color; // w's color := x's parent's    (case 4)
                x->parent->color    = false; // color x's parent BLACK       (case 4)
                w->right->color = false;     // color w's right child BLACK  (case 4)
                rotateLeft(x->parent);       // left rotation on x's parent  (case 4)
                x = root;            // finished work. bail out      (case 4)
            }                  //
        } else {                 // x is RIGHT-CHILD
            w = x->parent->left;       // grab x's sibling
            if (w->color == true) {    // if x's sibling is RED
                w->color = false;        // color w BLACK                (case 1)
                x->parent->color    = true;  // color x's parent RED         (case 1)
                rotateRight(x->parent);      // right rotation on x's parent (case 1)
                w = x->parent->left;         // make w be x's left sibling   (case 1)
            }
            if ((w->right->color == false) && (w->left->color == false)) {
                w->color = true;         // color w RED                  (case 2)
                x = x->parent;           // examine x's parent           (case 2)
            } else { //
                if (w->left->color == false) { //
                    w->right->color = false;   // color w's right child BLACK  (case 3)
                    w->color = true;       // color w RED                  (case 3)
                    t = x->parent;         // store x's parent
                    rotateLeft(w);         // left rotation on w           (case 3)
                    x->parent = t;         // restore x's parent
                    w = x->parent->left;       // make w be x's left sibling   (case 3)
                } //
                w->color = x->parent->color; // w's color := x's parent's    (case 4)
                x->parent->color    = false; // color x's parent BLACK       (case 4)
                w->left->color = false;      // color w's left child BLACK   (case 4)
                rotateRight(x->parent);      // right rotation on x's parent (case 4)
                x = root;                    // x is now the root            (case 4)
            }
        }
    }
    x->color = false;          // color x (the root) BLACK (exit)

    return;
}

// ******** Rotation Functions *******************************************

void splittree::rotateLeft(elementsp *x) {
    elementsp *y;
    // do pointer-swapping operations for left-rotation
    y = x->right;             // grab right child
    x->right = y->left;           // make x's RIGHT-CHILD be y's LEFT-CHILD
    y->left->parent = x;          // make x be y's LEFT-CHILD's parent
    y->parent = x->parent;        // make y's new parent be x's old parent

    if (x->parent == NULL) {
        root = y;               // if x was root, make y root
    } else {              //
        if (x == x->parent->left) { // if x is LEFT-CHILD, make y be x's parent's
            x->parent->left  = y; // left-child
        } else {
            x->parent->right = y; // right-child
        }
    }
    y->left   = x;        // make x be y's LEFT-CHILD
    x->parent = y;        // make y be x's parent

    return;
}

void splittree::rotateRight(elementsp *y) {
    elementsp *x;
    // do pointer-swapping operations for right-rotation
    x = y->left;               // grab left child
    y->left = x->right;            // replace left child yith x's right subtree
    x->right->parent = y;          // replace y as x's right subtree's parent

    x->parent = y->parent;         // make x's new parent be y's old parent
    if (y->parent == NULL) {
        root = x;            // if y was root, make x root
    } else {
        if (y == y->parent->right) { // if y is R-CHILD, make x be y's parent's
            y->parent->right = x;  // right-child
        } else {
            y->parent->left = x;   // left-child
        }
    }
    x->right  = y;         // make y be x's RIGHT-CHILD
    y->parent = x;         // make x be y's parent

    return;
}

// ***********************************************************************
// *** COPYRIGHT NOTICE **************************************************
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
// ***********************************************************************
// Author       : Aaron Clauset  ( aaronc@santafe.edu |
//                                 http://www.santafe.edu/~aaronc/ )
// Collaborators: Cristopher Moore and Mark E.J. Newman
// Project      : Hierarchical Random Graphs
// Location     : University of New Mexico, Dept. of Computer Science
//                AND Santa Fe Institute
// Created      : 21 June 2006
// Modified     : 23 December 2007 (cleaned up for public consumption)
//
// ************************************************************************

// ******** Constructor / Destructor *************************************

simpleGraph::simpleGraph(const int size): n(size), m(0), num_groups(0) {
    nodes = new simpleVert  [n];
    nodeLink = new simpleEdge* [n];
    nodeLinkTail = new simpleEdge* [n];
    A = new double* [n];
    for (int i = 0; i < n; i++) {
        nodeLink[i] = NULL; nodeLinkTail[i] = NULL;
        A[i] = new double [n];
        for (int j = 0; j < n; j++) {
            A[i][j] = 0.0;
        }
    }
    E = NULL;
}

simpleGraph::~simpleGraph() {
    simpleEdge *curr, *prev;
    for (int i = 0; i < n; i++) {
        curr = nodeLink[i];
        delete [] A[i];
        while (curr != NULL) {
            prev = curr;
            curr = curr->next;
            delete prev;
        }
    }
    curr = NULL; prev = NULL;
    if (E != NULL) {
        delete [] E;
        E = NULL;
    }
    delete [] A; A = NULL;
    delete [] nodeLink; nodeLink = NULL;
    delete [] nodeLinkTail;  nodeLinkTail   = NULL;
    delete [] nodes; nodes = NULL;
}

// ***********************************************************************

bool simpleGraph::addGroup(const int i, const int group_index) {
    if (i >= 0 && i < n) {
        nodes[i].group_true = group_index;
        return true;
    } else {
        return false;
    }
}

// ***********************************************************************

bool simpleGraph::addLink(const int i, const int j) {
    // Adds the directed edge (i,j) to the adjacency list for v_i
    simpleEdge* newedge;
    if (i >= 0 && i < n && j >= 0 && j < n) {
        A[i][j] = 1.0;
        newedge  = new simpleEdge;
        newedge->x = j;
        if (nodeLink[i] == NULL) {  // first neighbor
            nodeLink[i]  = newedge;
            nodeLinkTail[i] = newedge;
            nodes[i].degree = 1;
        } else {            // subsequent neighbor
            nodeLinkTail[i]->next = newedge;
            nodeLinkTail[i]       = newedge;
            nodes[i].degree++;
        }
        m++;            // increment edge count
        newedge = NULL;
        return true;
    } else {
        return false;
    }
}

// ***********************************************************************

bool simpleGraph::doesLinkExist(const int i, const int j) {
    // This function determines if the edge (i,j) already exists in the
    // adjacency list of v_i
    if (i >= 0 && i < n && j >= 0 && j < n) {
        if (A[i][j] > 0.1) {
            return true;
        } else {
            return false;
        }
    } else {
        return false;
    }
}

// **********************************************************************

double simpleGraph::getAdjacency(const int i, const int j) {
    if (i >= 0 && i < n && j >= 0 && j < n) {
        return A[i][j];
    } else {
        return -1.0;
    }
}

int simpleGraph::getDegree(const int i) {
    if (i >= 0 && i < n) {
        return nodes[i].degree;
    } else {
        return -1;
    }
}

int simpleGraph::getGroupLabel(const int i) {
    if (i >= 0 && i < n) {
        return nodes[i].group_true;
    } else {
        return -1;
    }
}

string simpleGraph::getName(const int i) {
    if (i >= 0 && i < n) {
        return nodes[i].name;
    } else {
        return "";
    }
}

// NOTE: The following three functions return addresses; deallocation
// of returned object is dangerous
simpleEdge* simpleGraph::getNeighborList(const int i) {
    if (i >= 0 && i < n) {
        return nodeLink[i];
    } else {
        return NULL;
    }
}
// END-NOTE

// *********************************************************************

int simpleGraph::getNumGroups() {
    return num_groups;
}
int simpleGraph::getNumLinks()  {
    return m;
}
int simpleGraph::getNumNodes()  {
    return n;
}
simpleVert* simpleGraph::getNode(const int i) {
    if (i >= 0 && i < n) {
        return &nodes[i];
    } else {
        return NULL;
    }
}

// **********************************************************************

bool simpleGraph::setName(const int i, const string text) {
    if (i >= 0 && i < n) {
        nodes[i].name = text;
        return true;
    } else {
        return false;
    }
}

// **********************************************************************

void simpleGraph::QsortMain (block* array, int left, int right) {
    if (right > left) {
        int pivot = left;
        int part  = QsortPartition(array, left, right, pivot);
        QsortMain(array, left,   part - 1);
        QsortMain(array, part + 1, right  );
    }
    return;
}

int simpleGraph::QsortPartition (block* array, int left, int right,
                                 int index) {
    block p_value, temp;
    p_value.x = array[index].x;
    p_value.y = array[index].y;

    // swap(array[p_value], array[right])
    temp.x = array[right].x;
    temp.y = array[right].y;
    array[right].x = array[index].x;
    array[right].y = array[index].y;
    array[index].x = temp.x;
    array[index].y = temp.y;

    int stored = left;
    for (int i = left; i < right; i++) {
        if (array[i].x <= p_value.x) {
            // swap(array[stored], array[i])
            temp.x = array[i].x;
            temp.y = array[i].y;
            array[i].x = array[stored].x;
            array[i].y = array[stored].y;
            array[stored].x = temp.x;
            array[stored].y = temp.y;
            stored++;
        }
    }
    // swap(array[right], array[stored])
    temp.x = array[stored].x;
    temp.y = array[stored].y;
    array[stored].x = array[right].x;
    array[stored].y = array[right].y;
    array[right].x = temp.x;
    array[right].y = temp.y;

    return stored;
}

// ***********************************************************************
