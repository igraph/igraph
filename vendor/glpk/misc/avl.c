/* avl.c (binary search tree) */

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

#include "avl.h"
#include "dmp.h"
#include "env.h"

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

AVL *avl_create_tree(int (*fcmp)(void *info, const void *key1,
      const void *key2), void *info)
{     /* create AVL tree */
      AVL *tree;
      tree = xmalloc(sizeof(AVL));
      tree->pool = dmp_create_pool();
      tree->root = NULL;
      tree->fcmp = fcmp;
      tree->info = info;
      tree->size = 0;
      tree->height = 0;
      return tree;
}

int avl_strcmp(void *info, const void *key1, const void *key2)
{     /* compare character string keys */
      xassert(info == info);
      return strcmp(key1, key2);
}

static AVLNODE *rotate_subtree(AVL *tree, AVLNODE *node);

AVLNODE *avl_insert_node(AVL *tree, const void *key)
{     /* insert new node into AVL tree */
      AVLNODE *p, *q, *r;
      short int flag;
      /* find an appropriate point for insertion */
      p = NULL; q = tree->root;
      while (q != NULL)
      {  p = q;
         if (tree->fcmp(tree->info, key, p->key) <= 0)
         {  flag = 0;
            q = p->left;
            p->rank++;
         }
         else
         {  flag = 1;
            q = p->right;
         }
      }
      /* create new node and insert it into the tree */
      r = dmp_get_atom(tree->pool, sizeof(AVLNODE));
      r->key = key; r->type = 0; r->link = NULL;
      r->rank = 1; r->up = p;
      r->flag = (short int)(p == NULL ? 0 : flag);
      r->bal = 0; r->left = NULL; r->right = NULL;
      tree->size++;
      if (p == NULL)
         tree->root = r;
      else
         if (flag == 0) p->left = r; else p->right = r;
      /* go upstairs to the root and correct all subtrees affected by
         insertion */
      while (p != NULL)
      {  if (flag == 0)
         {  /* the height of the left subtree of [p] is increased */
            if (p->bal > 0)
            {  p->bal = 0;
               break;
            }
            if (p->bal < 0)
            {  rotate_subtree(tree, p);
               break;
            }
            p->bal = -1; flag = p->flag; p = p->up;
         }
         else
         {  /* the height of the right subtree of [p] is increased */
            if (p->bal < 0)
            {  p->bal = 0;
               break;
            }
            if (p->bal > 0)
            {  rotate_subtree(tree, p);
               break;
            }
            p->bal = +1; flag = p->flag; p = p->up;
         }
      }
      /* if the root has been reached, the height of the entire tree is
         increased */
      if (p == NULL) tree->height++;
      return r;
}

void avl_set_node_type(AVLNODE *node, int type)
{     /* assign the type field of specified node */
      node->type = type;
      return;
}

void avl_set_node_link(AVLNODE *node, void *link)
{     /* assign the link field of specified node */
      node->link = link;
      return;
}

AVLNODE *avl_find_node(AVL *tree, const void *key)
{     /* find node in AVL tree */
      AVLNODE *p;
      int c;
      p = tree->root;
      while (p != NULL)
      {  c = tree->fcmp(tree->info, key, p->key);
         if (c == 0) break;
         p = (c < 0 ? p->left : p->right);
      }
      return p;
}

int avl_get_node_type(AVLNODE *node)
{     /* retrieve the type field of specified node */
      return node->type;
}

void *avl_get_node_link(AVLNODE *node)
{     /* retrieve the link field of specified node */
      return node->link;
}

static AVLNODE *find_next_node(AVL *tree, AVLNODE *node)
{     /* find next node in AVL tree */
      AVLNODE *p, *q;
      if (tree->root == NULL) return NULL;
      p = node;
      q = (p == NULL ? tree->root : p->right);
      if (q == NULL)
      {  /* go upstairs from the left subtree */
         for (;;)
         {  q = p->up;
            if (q == NULL) break;
            if (p->flag == 0) break;
            p = q;
         }
      }
      else
      {  /* go downstairs into the right subtree */
         for (;;)
         {  p = q->left;
            if (p == NULL) break;
            q = p;
         }
      }
      return q;
}

void avl_delete_node(AVL *tree, AVLNODE *node)
{     /* delete specified node from AVL tree */
      AVLNODE *f, *p, *q, *r, *s, *x, *y;
      short int flag;
      p = node;
      /* if both subtrees of the specified node are non-empty, the node
         should be interchanged with the next one, at least one subtree
         of which is always empty */
      if (p->left == NULL || p->right == NULL) goto skip;
      f = p->up; q = p->left;
      r = find_next_node(tree, p); s = r->right;
      if (p->right == r)
      {  if (f == NULL)
            tree->root = r;
         else
            if (p->flag == 0) f->left = r; else f->right = r;
         r->rank = p->rank; r->up = f;
         r->flag = p->flag; r->bal = p->bal;
         r->left = q; r->right = p;
         q->up = r;
         p->rank = 1; p->up = r; p->flag = 1;
         p->bal = (short int)(s == NULL ? 0 : +1);
         p->left = NULL; p->right = s;
         if (s != NULL) s->up = p;
      }
      else
      {  x = p->right; y = r->up;
         if (f == NULL)
            tree->root = r;
         else
            if (p->flag == 0) f->left = r; else f->right = r;
         r->rank = p->rank; r->up = f;
         r->flag = p->flag; r->bal = p->bal;
         r->left = q; r->right = x;
         q->up = r; x->up = r; y->left = p;
         p->rank = 1; p->up = y; p->flag = 0;
         p->bal = (short int)(s == NULL ? 0 : +1);
         p->left = NULL; p->right = s;
         if (s != NULL) s->up = p;
      }
skip: /* now the specified node [p] has at least one empty subtree;
         go upstairs to the root and adjust the rank field of all nodes
         affected by deletion */
      q = p; f = q->up;
      while (f != NULL)
      {  if (q->flag == 0) f->rank--;
         q = f; f = q->up;
      }
      /* delete the specified node from the tree */
      f = p->up; flag = p->flag;
      q = p->left != NULL ? p->left : p->right;
      if (f == NULL)
         tree->root = q;
      else
         if (flag == 0) f->left = q; else f->right = q;
      if (q != NULL) q->up = f, q->flag = flag;
      tree->size--;
      /* go upstairs to the root and correct all subtrees affected by
         deletion */
      while (f != NULL)
      {  if (flag == 0)
         {  /* the height of the left subtree of [f] is decreased */
            if (f->bal == 0)
            {  f->bal = +1;
               break;
            }
            if (f->bal < 0)
               f->bal = 0;
            else
            {  f = rotate_subtree(tree, f);
               if (f->bal < 0) break;
            }
            flag = f->flag; f = f->up;
         }
         else
         {  /* the height of the right subtree of [f] is decreased */
            if (f->bal == 0)
            {  f->bal = -1;
               break;
            }
            if (f->bal > 0)
               f->bal = 0;
            else
            {  f = rotate_subtree(tree, f);
               if (f->bal > 0) break;
            }
            flag = f->flag; f = f->up;
         }
      }
      /* if the root has been reached, the height of the entire tree is
         decreased */
      if (f == NULL) tree->height--;
      /* returns the deleted node to the memory pool */
      dmp_free_atom(tree->pool, p, sizeof(AVLNODE));
      return;
}

static AVLNODE *rotate_subtree(AVL *tree, AVLNODE *node)
{     /* restore balance of AVL subtree */
      AVLNODE *f, *p, *q, *r, *x, *y;
      xassert(node != NULL);
      p = node;
      if (p->bal < 0)
      {  /* perform negative (left) rotation */
         f = p->up; q = p->left; r = q->right;
         if (q->bal <= 0)
         {  /* perform single negative rotation */
            if (f == NULL)
               tree->root = q;
            else
               if (p->flag == 0) f->left = q; else f->right = q;
            p->rank -= q->rank;
            q->up = f; q->flag = p->flag; q->bal++; q->right = p;
            p->up = q; p->flag = 1;
            p->bal = (short int)(-q->bal); p->left = r;
            if (r != NULL) r->up = p, r->flag = 0;
            node = q;
         }
         else
         {  /* perform double negative rotation */
            x = r->left; y = r->right;
            if (f == NULL)
               tree->root = r;
            else
               if (p->flag == 0) f->left = r; else f->right = r;
            p->rank -= (q->rank + r->rank);
            r->rank += q->rank;
            p->bal = (short int)(r->bal >= 0 ? 0 : +1);
            q->bal = (short int)(r->bal <= 0 ? 0 : -1);
            r->up = f; r->flag = p->flag; r->bal = 0;
            r->left = q; r->right = p;
            p->up = r; p->flag = 1; p->left = y;
            q->up = r; q->flag = 0; q->right = x;
            if (x != NULL) x->up = q, x->flag = 1;
            if (y != NULL) y->up = p, y->flag = 0;
            node = r;
         }
      }
      else
      {  /* perform positive (right) rotation */
         f = p->up; q = p->right; r = q->left;
         if (q->bal >= 0)
         {  /* perform single positive rotation */
            if (f == NULL)
               tree->root = q;
            else
               if (p->flag == 0) f->left = q; else f->right = q;
            q->rank += p->rank;
            q->up = f; q->flag = p->flag; q->bal--; q->left = p;
            p->up = q; p->flag = 0;
            p->bal = (short int)(-q->bal); p->right = r;
            if (r != NULL) r->up = p, r->flag = 1;
            node = q;
         }
         else
         {  /* perform double positive rotation */
            x = r->left; y = r->right;
            if (f == NULL)
               tree->root = r;
            else
               if (p->flag == 0) f->left = r; else f->right = r;
            q->rank -= r->rank;
            r->rank += p->rank;
            p->bal = (short int)(r->bal <= 0 ? 0 : -1);
            q->bal = (short int)(r->bal >= 0 ? 0 : +1);
            r->up = f; r->flag = p->flag; r->bal = 0;
            r->left = p; r->right = q;
            p->up = r; p->flag = 0; p->right = x;
            q->up = r; q->flag = 1; q->left = y;
            if (x != NULL) x->up = p, x->flag = 1;
            if (y != NULL) y->up = q, y->flag = 0;
            node = r;
         }
      }
      return node;
}

void avl_delete_tree(AVL *tree)
{     /* delete AVL tree */
      dmp_delete_pool(tree->pool);
      xfree(tree);
      return;
}

/* eof */
