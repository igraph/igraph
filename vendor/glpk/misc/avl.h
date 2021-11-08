/* avl.h (binary search tree) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2000-2013 Free Software Foundation, Inc.
*  Written by Andrew Makhorin <mao@gnu.org>.
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

#ifndef AVL_H
#define AVL_H

typedef struct AVL AVL;
typedef struct AVLNODE AVLNODE;

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
