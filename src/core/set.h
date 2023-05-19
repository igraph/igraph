/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2009-2020  The igraph development team

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
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#ifndef IGRAPH_CORE_SET_H
#define IGRAPH_CORE_SET_H

#include "igraph_decls.h"
#include "igraph_error.h"
#include "igraph_types.h"

__BEGIN_DECLS

/* -------------------------------------------------- */
/* Flexible set                                       */
/* -------------------------------------------------- */

/**
 * Set containing integer numbers regardless of the order
 * \ingroup types
 */

#define IGRAPH_SET_PARAMETER_STACK_LENGTH 20
#define IGRAPH_SET_PARAMETER_STARTING_CAPACITY 100
#define PARENT(s,x) (*(s->parent + x))
#define LEFT(s,x) (*(s->left + x))
#define RIGHT(s,x) (*(s->right + x))
#define COLOR(s,x) (*(s->color + x))
#define SET(s,x) (s->pool + x)
#define SETWITHCHECK(s, x) (x!=-1 ? SET(s,x) : NULL)
#define ROOTINDEX(s) (s->root ? s->root->index : -1)

enum STACK_MODE {LEFT,SELF};

enum COLOR {RED,BLACK};

typedef struct Node
{
    igraph_integer_t data;
    igraph_integer_t index;
} igraph_set_internal_rbnode;

typedef struct s_set{
    igraph_set_internal_rbnode* root;
    igraph_integer_t size;
    igraph_integer_t* left;
    igraph_integer_t* right;
    igraph_integer_t* parent;
    igraph_set_internal_rbnode* pool;
    enum COLOR* color;
    igraph_integer_t capacity;
    
} igraph_set_t;


// typedef struct IteratorNode
// {
//     igraph_integer_t data;
//     igraph_integer_t index;
//     igraph_set_internal_rbnode* left;
//     igraph_set_internal_rbnode* right;
//     igraph_set_internal_rbnode* parent;
// } igraph_set_internal_iterator_node;

// struct Stack
// {
//     igraph_set_internal_iterator_node data;
//     enum STACK_MODE mode;
// } ;

// typedef struct s_set_itertor
// {
//     struct Stack stack[IGRAPH_SET_PARAMETER_STACK_LENGTH];
//     igraph_integer_t stack_index;
// } igraph_set_iterator_t;

typedef struct s_set_itertor
{
    igraph_integer_t* data;
    igraph_integer_t index;
} igraph_set_iterator_t;

#define IGRAPH_SET_NULL { 0,0,0 }
#define IGRAPH_SET_INIT_FINALLY(v) \
    do { igraph_set_init(v); \
        IGRAPH_FINALLY(igraph_set_destroy, v); } while (0)

IGRAPH_PRIVATE_EXPORT igraph_error_t igraph_set_init(igraph_set_t* set);
IGRAPH_PRIVATE_EXPORT void igraph_set_destroy(igraph_set_t* set);
IGRAPH_PRIVATE_EXPORT igraph_integer_t igraph_set_size(const igraph_set_t* set);
IGRAPH_PRIVATE_EXPORT igraph_bool_t igraph_set_contains(const igraph_set_t *set, igraph_integer_t e);
IGRAPH_PRIVATE_EXPORT igraph_error_t igraph_set_add(igraph_set_t* v, igraph_integer_t e);
IGRAPH_PRIVATE_EXPORT igraph_bool_t igraph_set_empty(const igraph_set_t* set);
IGRAPH_PRIVATE_EXPORT void igraph_set_delete(igraph_set_t* set, igraph_integer_t e);
IGRAPH_PRIVATE_EXPORT igraph_bool_t igraph_set_inited(igraph_set_t* set);
IGRAPH_PRIVATE_EXPORT igraph_bool_t igraph_set_iterate(const igraph_set_t *set, igraph_set_iterator_t* state,
                                                       igraph_integer_t* element);
IGRAPH_PRIVATE_EXPORT void igraph_set_iterator_init(const igraph_set_t* set, igraph_set_iterator_t* iterator);
IGRAPH_PRIVATE_EXPORT void igraph_set_print_tree(const igraph_set_t* set, FILE* output_stream);
__END_DECLS

#endif
