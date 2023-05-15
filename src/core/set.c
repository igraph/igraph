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
    set->pool = IGRAPH_CALLOC(IGRAPH_SET_PARAMETER_STARTING_CAPACITY, struct Node);
    IGRAPH_CHECK_OOM(set->pool, "Cannot reserve space for set.");
    set->capacity = IGRAPH_SET_PARAMETER_STARTING_CAPACITY;
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
    IGRAPH_FREE(set->pool);
    set->root = NULL;
    set->size = 0;
    set->pool = NULL;
    set->capacity = 0;
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

void LeftRotate(igraph_set_t *set, igraph_integer_t x) {
    igraph_integer_t y = SET(set, x)->right;
    SET(set, x)->right = SET(set, y)->left;

    if (SET(set, y)->left != -1) {
        SET(set, SET(set, y)->left)->parent = x;
    }

    SET(set, y)->parent = SET(set, x)->parent;

    if (SET(set, x)->parent == -1) {
        set->root = SET(set, y);
    }

    else if (x == SET(set, SET(set,x)->parent)->left) {
        SET(set, SET(set,x)->parent)->left = y;
    }

    else {
        SET(set, SET(set,x)->parent)->right = y;
    }

    SET(set, y)->left = x;

    SET(set, x)->parent = y;
}

void RightRotate(igraph_set_t *set, igraph_integer_t x) {
    igraph_integer_t y = SET(set, x)->left;
    SET(set, x)->left = SET(set, y)->right;

    if (SET(set,y)->right != -1) {
        SET(set, SET(set, y)->right)->parent = x;
    }

    SET(set, y)->parent = SET(set, x)->parent;

    if (SET(set, x)->parent == -1) {
        set->root = SET(set, y);
    }

    else if (x == SET(set, SET(set, x)->parent)->left) {
        SET(set, SET(set, x)->parent)->left = y;
    }

    else {
        SET(set, SET(set, x)->parent)->right = y;
    }

    SET(set, y)->right = x;
    SET(set, x)->parent = y;

}

void RB_insert_fixup(igraph_set_t *set, igraph_integer_t z) {
    igraph_integer_t grandparent = -1;
    igraph_integer_t parentpt = -1;

    while (((z) != SET_ROOT_INDEX(set)) && (SET(set, z)->color != BLACK) && (SET(set, SET(set, z)->parent)->color == RED)) {
        parentpt = SET(set, z)->parent;
        grandparent = SET(set, SET(set, z)->parent)->parent;

        if (parentpt == SET(set, grandparent)->left) {
            igraph_integer_t uncle = SET(set, grandparent)->right;

            if (uncle != -1 && SET(set,uncle)->color == RED) {
                SET(set, grandparent)->color = RED;
                SET(set, parentpt)->color = BLACK;
                SET(set, uncle)->color = BLACK;
                z = grandparent;
            }

            else {
                if (z == SET(set, parentpt)->right) {
                    LeftRotate(set, parentpt);
                    z = parentpt;
                    parentpt = SET(set, z)->parent;
                }

                RightRotate(set, grandparent);
                SET(set, parentpt)->color = BLACK;
                SET(set, grandparent)->color = RED;
                z = parentpt;
            }
        }

        else {
            igraph_integer_t uncle = SET(set, grandparent)->left;

            if (uncle != -1 && SET(set, uncle)->color == RED) {
                SET(set, grandparent)->color = RED;
                SET(set, parentpt)->color = BLACK;
                SET(set, uncle)->color = BLACK;
                z = grandparent;
            }

            else {
                if (z == SET(set, parentpt)->left) {
                    RightRotate(set, parentpt);
                    z = parentpt;
                    parentpt = SET(set, z)->parent;
                }

                LeftRotate(set, grandparent);
                SET(set, parentpt)->color = BLACK;
                SET(set, grandparent)->color = RED;
                z = parentpt;
            }
        }
    }
    set->root->color = BLACK;

}

struct Node* RB_insert(igraph_set_t* set, igraph_integer_t z) {
    
    igraph_integer_t y = -1;
    igraph_integer_t x = SET_ROOT_INDEX(set);//root

    while (x != -1) {
        y = x;
        if (SET(set, z)->data < SET(set, x)->data) {
            x = SET(set, x)->left;
        }

        else {
            x = SET(set, x)->right;
        }
    }
    SET(set, z)->parent = y;

    if (y == -1) {
        set->root = SET(set, z);
    }

    else if (SET(set, z)->data < SET(set,y)->data) {
        SET(set, y)->left = z;
    }

    else {
        SET(set, y)->right = z;
    }

    RB_insert_fixup(set, z);

    return set->root;
}

igraph_error_t igraph_set_reserve(igraph_set_t* set){
    /*In case someone uses a set after calling destroy is it, since destroy is just used a clear function*/
    if(set->capacity == 0){
        set->pool = IGRAPH_CALLOC(IGRAPH_SET_PARAMETER_STARTING_CAPACITY, struct Node);
        IGRAPH_CHECK_OOM(set->pool, "Cannot reserve space for set.");
        set->capacity = IGRAPH_SET_PARAMETER_STARTING_CAPACITY;
        return IGRAPH_SUCCESS;
    }
    set->capacity = set->capacity + set->capacity / 2;
    //printf("New set capacity %ld\n Old Capacity %ld\n", new_capacity, set->capacity);
    set->pool = IGRAPH_REALLOC(set->pool, set->capacity, struct Node);
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
    if(set->size == set->capacity){
        IGRAPH_CHECK(igraph_set_reserve(set));
    }
    struct Node* z = set->pool + set->size;
    z->data = e;
    z->left = -1;
    z->right = -1;
    z->parent = -1;
    z->color = RED;
    z->index = set->size;
    set->root = RB_insert(set, z->index);
    set->size++;
    return IGRAPH_SUCCESS;
}

igraph_integer_t Tree_minimum(igraph_set_t *set, igraph_integer_t node)
{
    while(SET(set, node)->left!=-1)
        node = SET(set, node)->left;

    return node;
}

void RB_delete_fixup(igraph_set_t *set,igraph_integer_t x)
{
    while(x!=SET_ROOT_INDEX(set) && SET(set, x)->color == BLACK)
    {
        if(x==SET(set, SET(set, x)->parent)->left)
        {
            igraph_integer_t w = SET(set, SET(set, x)->parent)->right;

            if(SET(set,w)->color==RED)
            {
                SET(set,w)->color = BLACK;
                SET(set, SET(set, x)->parent)->color = BLACK;
                LeftRotate(set, SET(set, x)->parent);
                w = SET(set, SET(set, x)->parent)->right;
            }

            if(SET(set, SET(set, w)->left)->color==BLACK && SET(set, SET(set, w)->right)->color == BLACK)
            {
                SET(set, w)->color = RED;
                x = SET(set, x)->parent;
            }

            else
            {
                if(SET(set, SET(set, w)->right)->color == BLACK)
                {
                    SET(set, SET(set, w)->left)->color = BLACK;
                    SET(set, w)->color = RED;
                    RightRotate(set, w);
                    w = SET(set, SET(set, x)->parent)->right;
                }

                SET(set, w)->color = SET(set, SET(set, x)->parent)->color;
                SET(set, SET(set, x)->parent)->color = BLACK;
                SET(set, SET(set, w)->right)->color = BLACK;
                LeftRotate(set, SET(set, x)->parent);
                x = SET_ROOT_INDEX(set);
            }
        }

        else
        {
            igraph_integer_t w = SET(set, SET(set, x)->parent)->left;

            if(SET(set, w)->color==RED)
            {
                SET(set, w)->color = BLACK;
                SET(set, SET(set, x)->parent)->color = BLACK;
                RightRotate(set, SET(set, x)->parent);
                w = SET(set, SET(set, x)->parent)->left;
            }

            if(SET(set, SET(set, w)->right)->color==BLACK && SET(set, SET(set, w)->left)->color == BLACK)
            {
                SET(set, w)->color = RED;
                x = SET(set, x)->parent;
            }

            else
            {
                if(SET(set, SET(set, w)->left)->color == BLACK)
                {
                    SET(set, SET(set, w)->right)->color = BLACK;
                    SET(set, w)->color = RED;
                    LeftRotate(set, w);
                    w = SET(set, SET(set, x)->parent)->left;
                }

                SET(set, w)->color = SET(set, SET(set, x)->parent)->color;
                SET(set, SET(set, x)->parent)->color = BLACK;
                SET(set, SET(set, w)->left)->color = BLACK;
                RightRotate(set, SET(set, x)->parent);
                x = SET_ROOT_INDEX(set);
            }
        }
    }
    SET(set, x)->color = BLACK;
}

void RB_transplat(igraph_set_t *set, igraph_integer_t u,igraph_integer_t v)
{
    if(SET(set, u)->parent == -1)
        set->root = SET(set, v);

    else if(u==SET(set, SET(set, u)->parent)->left)
        SET(set, SET(set, u)->parent)->left  = v;
    else
        SET(set, SET(set, u)->parent)->right = v;

    if(v!=-1) 
        SET(set, v)->parent = SET(set, u)->parent;
}

struct Node* RB_delete(igraph_set_t *set,igraph_integer_t z)
{
    igraph_integer_t y = z;
    enum COLOR yoc;
    yoc = SET(set, z)->color; // y's original color

    igraph_integer_t x;

    if(SET(set, z)->left==-1 )
    {
        x = SET(set, z)->right;
        RB_transplat(set, z, SET(set, z)->right);
    }

    else if(SET(set, z)->right==-1 )
    {
        x = SET(set, z)->left;
        RB_transplat(set, z, SET(set, z)->left);
    }

    else
    {
        y = Tree_minimum(set, SET(set, z)->right);
        yoc = SET(set, y)->color;
        x = SET(set, y)->right;

        if(SET(set, y)->parent==z && x!=-1)
            SET(set, x)->parent = y;

        if(SET(set, y)->parent!=z)
        {
            RB_transplat(set, y, SET(set, y)->right);
            SET(set, y)->right = SET(set, z)->right;
            SET(set, SET(set, y)->right)->parent = y;
        }

        RB_transplat(set, z, y);
        SET(set, y)->left = SET(set, z)->left;
        SET(set, SET(set, y)->left)->parent = y;
        SET(set, y)->color = SET(set, z)->color;
    }

    if(yoc==BLACK)
        RB_delete_fixup(set, x);

    return set->root;
}
struct Node* BST_search(const igraph_set_t *set, struct Node* node, int x)
{
    if(node==NULL || node->data == x)
        return node;

    if(node->data > x)
       return  BST_search(set, SET(set, node->left),x);
    else
        return BST_search(set, SET(set, node->right),x);
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
    struct Node* node_to_delete = BST_search(set, set->root, e);
    if(node_to_delete == NULL){
        return ;
    }    
    set->root = RB_delete(set, node_to_delete->index);
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
    return BST_search(set, set->root, e);
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
    if (set->root == NULL) {
        iterator->stack_index = -1;
        return ;
    }
    iterator->stack[0].data = *(set->root);
    iterator->stack_index = 0;
    iterator->stack[0].mode = LEFT;
}

igraph_integer_t iterate_self(const igraph_set_t *set, igraph_set_iterator_t *state) {
    struct Node *node = &(state->stack[state->stack_index].data);
    igraph_integer_t data = node->data;
    if (node->right != -1) {
        state->stack[state->stack_index].data = *SET(set, node->right);
        state->stack[state->stack_index].mode = LEFT;
    } else if (node->parent != -1 && state->stack_index > 0 ) {
        state->stack_index--;
        state->stack[state->stack_index].mode = SELF;
    } else {
        state->stack_index--;
    }
    return data;
}

igraph_integer_t iterate_left(const igraph_set_t *set, igraph_set_iterator_t *state) {
    struct Node *node = &(state->stack[state->stack_index].data);

    for ( ; node->left != -1 ; state->stack_index++) {
        state->stack[state->stack_index + 1 ].data = *SET(set, node->left);
        state->stack[state->stack_index + 1 ].mode = LEFT;
        node = SET(set, node->left);
    }
    state->stack_index--;
    if(state->stack_index >= 0){
        state->stack[state->stack_index].mode = SELF;
    }
    // printf("inside iterate left, last state data %ld\n", state->stack[state->stack_index].data.data);
    return node->data;
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
    enum STACK_MODE mode = state->stack[state->stack_index].mode;
    switch (mode) {
        case LEFT:
            *element = iterate_left(set, state);
            return true;
        case SELF:
            *element = iterate_self(set, state);
            return true;
    }
}

void print2DUtil(const igraph_set_t* set, struct Node* root, int space, FILE * output_stream) {
    #define COUNT 10
    if (root == NULL) {
        return;
    }
    space += COUNT;
    print2DUtil(set, SET(set, root->right), space,output_stream);
    fprintf(output_stream, "\n");
    for (int i = COUNT; i < space; i++) {
        fprintf(output_stream, " ");
    }
    printf("%" IGRAPH_PRId "\n", root->data);
    print2DUtil(set, SET(set, root->left), space,output_stream);
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
    print2DUtil(set, set->root, 0, output_stream);
}
