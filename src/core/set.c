/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
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
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include "igraph_memory.h"

#include "core/set.h"

#include <string.h>     /* memmove */

/**
 * \ingroup set
 * \function igraph_set_init
 * \brief Initializes a set.
 *
 * Initializes an empty set (with zero elements). 
 * \param set Pointer to the set to be initialized.in the set.
 *
 * Time complexity: O(1)
 */
igraph_error_t igraph_set_init(igraph_set_t *set) {
    IGRAPH_ASSERT(set != NULL);
    set->root = NULL;
    set->size = 0;
    set->capacity_used = 0;
    set->capacity = IGRAPH_SET_PARAMETER_STARTING_CAPACITY;
    set->pool = IGRAPH_CALLOC(IGRAPH_SET_PARAMETER_STARTING_CAPACITY, igraph_set_internal_node_t*);
    IGRAPH_CHECK_OOM(set->pool, "Cannot reserve space for set.");
    return IGRAPH_SUCCESS;
}


/**
 * \ingroup set
 * \function igraph_set_destroy
 * \brief Destroys a set object.
 *
 * \param set Pointer to the set to be destroyed.
 *
 * Time complexity: Operating System Dependent.
 */
void igraph_set_destroy(igraph_set_t* set) {
    IGRAPH_ASSERT(set != NULL);
    for(igraph_integer_t i = 0 ;i < set->capacity_used ; i++){
        IGRAPH_FREE(set->pool[i]);
    }
    IGRAPH_FREE(set->pool);
    set->capacity = 0;
    set->root = NULL;
    set->size = 0;
    set->capacity_used = 0;
}


/**
 * \ingroup set
 * \function igraph_set_empty
 * \brief Decides whether the size of the set is zero.
 *
 * \param set The set object.
 * \return Non-zero number if the size of the set is not zero and
 *         zero otherwise.
 *
 * Time complexity: O(1).
 */
igraph_bool_t igraph_set_empty(const igraph_set_t* set) {
    IGRAPH_ASSERT(set != NULL);
    return set->size == 0;
}


/**
 * \ingroup set
 * \function igraph_set_size
 * \brief Gives the size of the set.
 *
 * The number of elements in the set.
 *
 * \param set The set object
 * \return The size of the set.
 *
 * Time complexity: O(1).
 */

igraph_integer_t igraph_set_size(const igraph_set_t* set) {
    IGRAPH_ASSERT(set != NULL);
    return set->size;
}

void LeftRotate(igraph_set_internal_node_t** T, igraph_set_internal_node_t** x) {
    igraph_set_internal_node_t* y = (*x)->right;
    (*x)->right = y->left;

    if (y->left != NULL) {
        y->left->parent = *x;
    }

    y->parent = (*x)->parent;

    if ((*x)->parent == NULL) {
        *T = y;
    }

    else if (*x == (*x)->parent->left) {
        (*x)->parent->left = y;
    }

    else {
        (*x)->parent->right = y;
    }

    y->left = *x;

    (*x)->parent = y;

}

void RightRotate(igraph_set_internal_node_t** T, igraph_set_internal_node_t** x) {
    igraph_set_internal_node_t* y = (*x)->left;
    (*x)->left = y->right;

    if (y->right != NULL) {
        y->right->parent = *x;
    }

    y->parent = (*x)->parent;

    if ((*x)->parent == NULL) {
        *T = y;
    }

    else if ((*x) == (*x)->parent->left) {
        (*x)->parent->left = y;
    }

    else {
        (*x)->parent->right = y;
    }

    y->right = *x;
    (*x)->parent = y;

}

void RB_insert_fixup(igraph_set_internal_node_t** T, igraph_set_internal_node_t** z) {
    igraph_set_internal_node_t* grandparent = NULL;
    igraph_set_internal_node_t* parentpt = NULL;

    while (((*z) != *T) && ((*z)->color != BLACK) && ((*z)->parent->color == RED)) {
        parentpt = (*z)->parent;
        grandparent = (*z)->parent->parent;

        if (parentpt == grandparent->left) {
            igraph_set_internal_node_t* uncle = grandparent->right;

            if (uncle != NULL && uncle->color == RED) {
                grandparent->color = RED;
                parentpt->color = BLACK;
                uncle->color = BLACK;
                *z = grandparent;
            }

            else {
                if ((*z) == parentpt->right) {
                    LeftRotate(T, &parentpt);
                    (*z) = parentpt;
                    parentpt = (*z)->parent;
                }

                RightRotate(T, &grandparent);
                parentpt->color = BLACK;
                grandparent->color = RED;
                (*z) = parentpt;
            }
        }

        else {
            igraph_set_internal_node_t* uncle = grandparent->left;

            if (uncle != NULL && uncle->color == RED) {
                grandparent->color = RED;
                parentpt->color = BLACK;
                uncle->color = BLACK;
                (*z) = grandparent;
            }

            else {
                if ((*z) == parentpt->left) {
                    RightRotate(T, &parentpt);
                    (*z) = parentpt;
                    parentpt = (*z)->parent;
                }

                LeftRotate(T, &grandparent);
                parentpt->color = BLACK;
                grandparent->color = RED;
                (*z) = parentpt;
            }
        }
    }
    (*T)->color = BLACK;

}

igraph_set_internal_node_t* RB_insert(igraph_set_internal_node_t* T, igraph_integer_t data, igraph_set_internal_node_t* z) {
    z->data = data;
    z->left = NULL;
    z->right = NULL;
    z->parent = NULL;
    z->color = RED;

    igraph_set_internal_node_t* y = NULL;
    igraph_set_internal_node_t* x = T;//root

    while (x != NULL) {
        y = x;
        if (z->data < x->data) {
            x = x->left;
        }

        else {
            x = x->right;
        }
    }
    z->parent = y;

    if (y == NULL) {
        T = z;
    }

    else if (z->data < y->data) {
        y->left = z;
    }

    else {
        y->right = z;
    }

    RB_insert_fixup(&T, &z);

    return T;
}


igraph_error_t igraph_set_reserve(igraph_set_t* set){
    if(set->capacity == 0){
        IGRAPH_CHECK(igraph_set_init(set));
    }
    set->capacity *= 2;
    set->pool = IGRAPH_REALLOC(set->pool, set->capacity, igraph_set_internal_node_t *);
    IGRAPH_CHECK_OOM(set->pool, "Cannot reserve space for set.");
    return IGRAPH_SUCCESS;
}

/**
 * \ingroup set
 * \function igraph_set_add
 * \brief Adds an element to the set.
 *
 * \param set The set object.
 * \param e The element to be added.
 * \return Error code:
 *         \c IGRAPH_ENOMEM: not enough memory.
 *
 * Time complexity: O(log(n)), n is the number of elements in \p set.
 */
igraph_error_t igraph_set_add(igraph_set_t* set, igraph_integer_t e) {
    if (igraph_set_contains(set, e)) {
        return IGRAPH_SUCCESS;
    } 
    if(set->capacity_used == set->capacity){
        IGRAPH_CHECK(igraph_set_reserve(set));
    }
    set->pool[set->capacity_used] = IGRAPH_CALLOC(IGRAPH_SET_PARAMETER_STARTING_CAPACITY, igraph_set_internal_node_t);
    set->root = RB_insert(set->root, e, set->pool[set->size]);
    set->size++;
    set->capacity_used++;
    return IGRAPH_SUCCESS;
}

igraph_set_internal_node_t* Tree_minimum(igraph_set_internal_node_t* node)
{
    while(node->left!=NULL)
        node = node->left;

    return node;
}

void RB_delete_fixup(igraph_set_internal_node_t** T, igraph_set_internal_node_t** x)
{
    while((*x)!=*T && (*x)->color == BLACK)
    {
        if((*x)==(*x)->parent->left)
        {
            igraph_set_internal_node_t* w = (*x)->parent->right;

            if(w->color==RED)
            {
                w->color = BLACK;
                (*x)->parent->color = BLACK;
                LeftRotate(T,&((*x)->parent));
                w = (*x)->parent->right;
            }

            if(w->left->color==BLACK && w->right->color == BLACK)
            {
                w->color = RED;
                (*x) = (*x)->parent;
            }

            else
            {
                if(w->right->color == BLACK)
                {
                    w->left->color = BLACK;
                    w->color = RED;
                    RightRotate(T,&w);
                    w = (*x)->parent->right;
                }

                w->color = (*x)->parent->color;
                (*x)->parent->color = BLACK;
                w->right->color = BLACK;
                LeftRotate(T,&((*x)->parent));
                (*x) = *T;
            }
        }

        else
        {
            igraph_set_internal_node_t* w = (*x)->parent->left;

            if(w->color==RED)
            {
                w->color = BLACK;
                (*x)->parent->color = BLACK;
                RightRotate(T,&((*x)->parent));
                w = (*x)->parent->left;
            }

            if(w->right->color==BLACK && w->left->color == BLACK)
            {
                w->color = RED;
                (*x) = (*x)->parent;
            }

            else
            {
                if(w->left->color == BLACK)
                {
                    w->right->color = BLACK;
                    w->color = RED;
                    LeftRotate(T,&w);
                    w = (*x)->parent->left;
                }

                w->color = (*x)->parent->color;
                (*x)->parent->color = BLACK;
                w->left->color = BLACK;
                RightRotate(T,&((*x)->parent));
                (*x) = *T;
            }
        }
    }
    (*x)->color = BLACK;

}

void RB_transplat(igraph_set_internal_node_t** T, igraph_set_internal_node_t** u,igraph_set_internal_node_t** v)
{
    if((*u)->parent == NULL)
        *T = *v;

    else if((*u)==(*u)->parent->left)
        (*u)->parent->left = *v;
    else
        (*u)->parent->right = *v;

    if((*v)!=NULL) 
        (*v)->parent = (*u)->parent;
}

igraph_set_internal_node_t* RB_delete(igraph_set_internal_node_t *T,igraph_set_internal_node_t* z)
{
    igraph_set_internal_node_t *y = z;
    enum COLOR yoc;
    yoc = z->color; // y's original color

    igraph_set_internal_node_t* x;

    if(z->left==NULL )
    {
        x = z->right;
        RB_transplat(&T,&z,&(z->right));
    }

    else if(z->right==NULL )
    {
        x = z->left;
        RB_transplat(&T,&z,&(z->left));
    }

    else
    {
        y = Tree_minimum(z->right);
        yoc = y->color;
        x = y->right;

        if(y->parent==z && x!=NULL)
            x->parent = y;

        if(y->parent!=z)
        {
            RB_transplat(&T,&y,&(y->right));
            y->right = z->right;
            y->right->parent = y;
        }

        RB_transplat(&T,&z,&y);
        y->left = z->left;
        y->left->parent = y;
        y->color = z->color;
    }

    if(yoc==BLACK)
        RB_delete_fixup(&T,&x);

    return T;
}
igraph_set_internal_node_t* BST_search(igraph_set_internal_node_t* root, int x)
{
    if(root==NULL || root->data == x)
        return root;

    if(root->data > x)
       return  BST_search(root->left,x);
    else
        return BST_search(root->right,x);
}


/**
 * \ingroup set
 * \function igraph_set_delete
 * \brief Removes an element to the set.
 *
 * \param set The set object.
 * \param e The element to be removed.
 *
 * Time complexity: O(log(n)), n is the number of elements in \p set.
 */
void igraph_set_delete(igraph_set_t* set, igraph_integer_t e){
    IGRAPH_ASSERT(set != NULL);
    igraph_set_internal_node_t* node_to_delete = BST_search(set->root, e);
    if(node_to_delete == NULL){
        return ;
    }    
    set->root = RB_delete(set->root, node_to_delete);
    set->size--;
}

igraph_bool_t BST_Search(const igraph_set_internal_node_t* node, igraph_integer_t e) {
    if (node == NULL) {
        return false;
    }
    if (node->data == e) {
        return true;
    }
    if (node->data > e) {
        return BST_Search(node->left, e);
    } else {
        return BST_Search(node->right, e);
    }
}

/**
 * \ingroup set
 * \function igraph_set_contains
 * \brief Checks whether a given element is in the set or not.
 *
 * \param set The set object.
 * \param e The element being sought.
 * \return Positive integer (true) if \p e is found, zero (false) otherwise.
 *
 * Time complexity: O(log(n)), n is the number of elements in \p set.
 */
igraph_bool_t igraph_set_contains(const igraph_set_t* set, igraph_integer_t e) {
    IGRAPH_ASSERT(set != NULL);
    return BST_Search(set->root, e);
}


/**
 * \ingroup set
 * \function igraph_set_inited
 * \brief Determines whether a set is initialized or not.
 *
 * \param set The set object.
 *
 * Time complexity: O(1)
 */
igraph_bool_t igraph_set_inited(igraph_set_t* set){
    return true;
}


/**
 * \ingroup set
 * \function igraph_set_iterator_init
 * \brief Initializes a set iterator
 *
 * \param set Pointer to the set to be iterated upon.
 * \param iterator Pointer to the Iterator to be initialised 
 *
 * Time complexity: O(1)
 */


void igraph_set_iterator_init(const igraph_set_t* set, igraph_set_iterator_t* iterator) {
    IGRAPH_ASSERT(set != NULL);
    igraph_set_internal_node_t* node = set->root;
    if(node == NULL){
        iterator->stack_index = -1;
    }
    iterator->stack_index = 0;
    iterator->stack[0].data = *set->root;
    iterator->stack[0].visited = false;
}

/** 
 * \ingroup set
 * \function igraph_set_iterate
 * \brief Iterates through the element of the set.
 *
 * Elements are returned in an sorted order. 
 * Inserting elements while iterating might be ignored.
 * Deleting elements while iterating might be ignored (they might still be returned).
 * Deleting the element the iterate is pointing to(the latest one returned) will have no issue. 
 *
 * \param set The set object.
 * \param state Internal state of the iteration.
 *   This should be a pointer to an \c igraph_set_iterator_t variable
 *   which must be first initialised via igraph_set_iterator_init.
 *   The object must not be adjusted and its value should
 *   not be used for anything during the iteration.
 * \param element The next element or 0 (if the iteration
 *   has ended) is returned here.
 *
 * \return Nonzero if there are more elements, zero otherwise.
 */
igraph_bool_t igraph_set_iterate(const igraph_set_t *set, igraph_set_iterator_t *state,
                                 igraph_integer_t *element) {
    IGRAPH_ASSERT(set != NULL);
    IGRAPH_ASSERT(state != NULL);
    if (state->stack_index < 0) {
        element = NULL;
        return false;
    }
    while(true){
        igraph_bool_t visisted = state->stack[state->stack_index].visited;
        if(visisted){
            *element = state->stack[state->stack_index].data.data;
            state->stack_index--;
            return true; 
        }
        igraph_set_internal_node_t node = state->stack[state->stack_index].data;
        state->stack_index--;
        if(node.right){
            state->stack_index++;
            state->stack[state->stack_index].data = *node.right;            
            state->stack[state->stack_index].visited = false;
        }
        state->stack_index++;
        state->stack[state->stack_index].data = node;
        state->stack[state->stack_index].visited = true;
        if(node.left){
            state->stack_index++;
            state->stack[state->stack_index].data = *node.left;
            state->stack[state->stack_index].visited = false;
        }
    }
    return true; //control flow won't get here but compiler doesn't know that
}

void print2DUtil(igraph_set_internal_node_t* root, int space, FILE * output_stream) {
    #define COUNT 10
    if (root == NULL) {
        return;
    }
    space += COUNT;
    print2DUtil(root->right, space,output_stream);
    fprintf(output_stream, "\n");
    for (int i = COUNT; i < space; i++) {
        fprintf(output_stream, " ");
    }
    printf("%" IGRAPH_PRId "\n", root->data);
    print2DUtil(root->left, space,output_stream);
}

/**
 * \ingroup set
 * \function igraph_set_print_tree
 * \brief Prints the internal BST of the set in 2D format.
 *
 * For example the set with number 1 to 7 will be printed as
 *
 *                     7
 * 
 *           3
 * 
 *                     6
 * 
 * 1
 * 
 *                     5
 * 
 *           2
 * 
 *                     4
 * 
 * \param set The set to be prited.
 * \param output_stream The file stream to write the result to. Set this to stdout to write to console.
 * Time complexity: O(n * log(n)), n is the number of elements in \p set.
 */
void igraph_set_print_tree(const igraph_set_t* set, FILE * output_stream){
    print2DUtil(set->root, 0, output_stream);
}
