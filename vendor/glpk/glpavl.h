/* glpavl.h (binary search tree) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*
*  Copyright (C) 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008,
*  2009, 2010 Andrew Makhorin, Department for Applied Informatics,
*  Moscow Aviation Institute, Moscow, Russia. All rights reserved.
*  E-mail: <mao@gnu.org>.
*
*  GLPK is free software: you can redistribute it and/or modify it
*  under the terms of the GNU General Public License as published by
*  the Free Software Foundation, either version 3 of the License, or
*  (at your option) any later version.
*
*  GLPK is distributed in the hope that it will be useful, but WITHOUT
*  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
*  or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public
*  License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with GLPK. If not, see <http://www.gnu.org/licenses/>.
***********************************************************************/

#ifndef GLPAVL_H
#define GLPAVL_H

#include "glpdmp.h"

typedef struct AVL AVL;
typedef struct AVLNODE AVLNODE;

struct AVL
{     /* AVL tree (Adelson-Velsky & Landis binary search tree) */
      DMP *pool;
      /* memory pool for allocating nodes */
      AVLNODE *root;
      /* pointer to the root node */
      int (*fcmp)(void *info, const void *key1, const void *key2);
      /* application-defined key comparison routine */
      void *info;
      /* transit pointer passed to the routine fcmp */
      int size;
      /* the tree size (the total number of nodes) */
      int height;
      /* the tree height */
};

struct AVLNODE
{     /* node of AVL tree */
      const void *key;
      /* pointer to the node key (data structure for representing keys
         is supplied by the application) */
      int rank;
      /* node rank = relative position of the node in its own subtree =
         the number of nodes in the left subtree plus one */
      int type;
      /* reserved for the application specific information */
      void *link;
      /* reserved for the application specific information */
      AVLNODE *up;
      /* pointer to the parent node */
      short int flag;
      /* node flag:
         0 - this node is the left child of its parent (or this node is
             the root of the tree and has no parent)
         1 - this node is the right child of its parent */
      short int bal;
      /* node balance = the difference between heights of the right and
         left subtrees:
         -1 - the left subtree is higher than the right one;
          0 - the left and right subtrees have the same height;
         +1 - the left subtree is lower than the right one */
      AVLNODE *left;
      /* pointer to the root of the left subtree */
      AVLNODE *right;
      /* pointer to the root of the right subtree */
};

#define avl_create_tree _glp_avl_create_tree
AVL *avl_create_tree(int (*fcmp)(void *info, const void *key1,
      const void *key2), void *info);
/* create AVL tree */

#define avl_strcmp _glp_avl_strcmp
int avl_strcmp(void *info, const void *key1, const void *key2);
/* compare character string keys */

#define avl_insert_node _glp_avl_insert_node
AVLNODE *avl_insert_node(AVL *tree, const void *key);
/* insert new node into AVL tree */

#define avl_set_node_type _glp_avl_set_node_type
void avl_set_node_type(AVLNODE *node, int type);
/* assign the type field of specified node */

#define avl_set_node_link _glp_avl_set_node_link
void avl_set_node_link(AVLNODE *node, void *link);
/* assign the link field of specified node */

#define avl_find_node _glp_avl_find_node
AVLNODE *avl_find_node(AVL *tree, const void *key);
/* find node in AVL tree */

#define avl_get_node_type _glp_avl_get_node_type
int avl_get_node_type(AVLNODE *node);
/* retrieve the type field of specified node */

#define avl_get_node_link _glp_avl_get_node_link
void *avl_get_node_link(AVLNODE *node);
/* retrieve the link field of specified node */

#define avl_delete_node _glp_avl_delete_node
void avl_delete_node(AVL *tree, AVLNODE *node);
/* delete specified node from AVL tree */

#define avl_delete_tree _glp_avl_delete_tree
void avl_delete_tree(AVL *tree);
/* delete AVL tree */

#endif

/* eof */
