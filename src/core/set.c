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

struct Node* SET(const igraph_set_t *set, igraph_integer_t x){
    return x != -1 ? set->pool + x : NULL ;
}
igraph_integer_t SET_ROOT_INDEX(const igraph_set_t *set){
    return set->root ? set->root->index : -1;
}

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
    set->pool = IGRAPH_CALLOC(IGRAPH_SET_PARAMETER_STARTING_CAPACITY, igraph_set_internal_rbnode);
    IGRAPH_CHECK_OOM(set->pool, "Cannot reserve space for set.");
    set->left = IGRAPH_CALLOC(IGRAPH_SET_PARAMETER_STARTING_CAPACITY, igraph_integer_t);
    IGRAPH_CHECK_OOM(set->left, "Cannot reserve space for set.");
    set->right = IGRAPH_CALLOC(IGRAPH_SET_PARAMETER_STARTING_CAPACITY, igraph_integer_t);
    IGRAPH_CHECK_OOM(set->right, "Cannot reserve space for set.");
    set->parent = IGRAPH_CALLOC(IGRAPH_SET_PARAMETER_STARTING_CAPACITY, igraph_integer_t);
    IGRAPH_CHECK_OOM(set->parent, "Cannot reserve space for set.");
    set->color = IGRAPH_CALLOC(IGRAPH_SET_PARAMETER_STARTING_CAPACITY, enum COLOR);
    IGRAPH_CHECK_OOM(set->color, "Cannot reserve space for set.");
    set->capacity = IGRAPH_SET_PARAMETER_STARTING_CAPACITY;
    
    return IGRAPH_SUCCESS;
}

void set_destroy_internal(struct Node* node) {
    if(node->left != NULL)
        set_destroy_internal(node->left);
    
    if(node->right != NULL)
        set_destroy_internal(node->right);
    IGRAPH_FREE(node);
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
// #define OLD_SET_DESTROY_IMPLEMENTATION
void igraph_set_destroy(igraph_set_t* set) {
    IGRAPH_ASSERT(set != NULL);
    if(set->size < 1){
        return ;
    }
    if(set->size < RECUSIVE_DELETE_SIZE_LIMIT || true){
        set_destroy_internal(set->root);
        set->root = NULL;
        set->size = 0;
        return ;
    }
    // printf("destructor entered\n");
    // fflush(stdout);
    // return;
    igraph_set_internal_node_t *stack_first[STACK_LENGTH];
    igraph_bool_t stack_second[STACK_LENGTH];
    stack_first[0] = set->root;
    stack_second[0] = false;
    igraph_integer_t stack_index = 0;
    while(stack_index >=0){
        IGRAPH_ASSERT(stack_index < STACK_LENGTH);
        if(stack_second[stack_index]){
            IGRAPH_FREE(stack_first[stack_index]);
            stack_index--;
            continue;
        }
        igraph_set_internal_node_t *node = stack_first[stack_index];
        stack_second[stack_index] = true;
        
        if(node->right != NULL){
            stack_index++;
            stack_second[stack_index] = false;
            stack_first[stack_index] = node->right;
        }

        if(node->left != NULL){
            stack_index++;
            stack_second[stack_index] = false;
            stack_first[stack_index] = node->left;
        }
    }
    set->root = NULL;
    set->size = 0;
    // printf("destructor finished\n");
    // fflush(stdout);
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

void LeftRotate(igraph_set_t* set, igraph_integer_t x){
    igraph_integer_t y = RIGHT(set, x);
    RIGHT(set, x) = LEFT(set, y);

    if (LEFT(set, y) != -1) {
        PARENT(set, LEFT(set, y)) = x;
    }

    PARENT(set, y) = PARENT(set, x);

    if (PARENT(set, x) == -1) {
        set->root = SET(set, y);
    }

    else if (x == LEFT(set, PARENT(set, x))) {
        LEFT(set, PARENT(set, x)) = y;
    }

    else {
        RIGHT(set, PARENT(set, x)) = y;
    }

    LEFT(set, y) = x;

    PARENT(set, x) = y;

}

void RightRotate(igraph_set_t* set, igraph_integer_t x) {
    igraph_integer_t y = LEFT(set, x);
    LEFT(set, x) = RIGHT(set, y);
    if (RIGHT(set, y) != -1) {
        PARENT(set, RIGHT(set, y)) = x;
    }

    PARENT(set, y) = PARENT(set, x);

    if (PARENT(set, x) == -1) {
        set->root = SET(set, y);
    }

    else if (x == LEFT(set, PARENT(set, x))) {
        LEFT(set, PARENT(set, x)) = y;
    }

    else {
        RIGHT(set, PARENT(set, x)) = y;
    }

    RIGHT(set, y) = x;

    PARENT(set, x) = y;

}

void RB_insert_fixup(igraph_set_t* set, igraph_integer_t z) {
    igraph_integer_t grandparent = -1;
    igraph_integer_t parentpt = -1;


    while ((z != ROOTINDEX(set)) && (COLOR(set, z) != BLACK) && (COLOR(set, PARENT(set, z)) == RED)) {
        parentpt = PARENT(set, z);
        grandparent = PARENT(set, PARENT(set, z));

        if (parentpt == LEFT(set, grandparent)) {
            igraph_integer_t uncle = RIGHT(set, grandparent);

            if (uncle != -1 && COLOR(set, uncle) == RED) {
                COLOR(set, grandparent) = RED;
                COLOR(set, parentpt) = BLACK;
                COLOR(set, uncle) = BLACK;
                z = grandparent;
            }

            else {
                if (z == RIGHT(set, parentpt)) {
                    LeftRotate(set, parentpt);
                    z = parentpt;
                    parentpt = PARENT(set, z);
                }

                RightRotate(set, grandparent);
                COLOR(set, parentpt) = BLACK;
                COLOR(set, grandparent) = RED;
                z = parentpt;
            }
        }

        else {
            igraph_integer_t uncle = LEFT(set, grandparent);

            if (uncle != -1 && COLOR(set, uncle) == RED) {
                COLOR(set, grandparent) = RED;
                COLOR(set, parentpt) = BLACK;
                COLOR(set, uncle) = BLACK;
                z = grandparent;
            }

            else {
                if (z == PARENT(set, parentpt)) {
                    RightRotate(set, parentpt);
                    z = parentpt;
                    parentpt = PARENT(set, z);
                }

                LeftRotate(set, grandparent);
                COLOR(set, parentpt) = BLACK;
                COLOR(set, grandparent) = RED;
                z = parentpt;
            }
        }
    }
    COLOR(set, set->root->index) = BLACK;

}


igraph_set_internal_rbnode* RB_insert(
    igraph_set_t* set, igraph_set_internal_rbnode* z, 
    igraph_integer_t data, igraph_integer_t index
){
    z->data = data;
    z->index = index;
    PARENT(set, index) = -1;
    RIGHT(set, index) = -1;
    LEFT(set, index) = -1;
    COLOR(set, index) = RED;

    igraph_integer_t y = -1;
    igraph_integer_t x = ROOTINDEX(set);//root
    while (x != -1) {
        y = x;
        if (z->data < SET(set, x)->data) {
            x = LEFT(set, x);
        }

        else {
            x = RIGHT(set, x);
        }
    }
    PARENT(set, z->index) = y;

    if (y == -1) {
        set->root = z;
    }

    else if (z->data < SET(set, y)->data) {
        LEFT(set, y) = z->index;
    }

    else {
        RIGHT(set, y) = z->index;
    }

    RB_insert_fixup(set, z->index);

    return set->root;
}


igraph_error_t igraph_set_reserve(igraph_set_t* set){
    /*In case someone uses a set after calling destroy is it, since destroy is just used a clear function*/
    if(set->capacity == 0){
        IGRAPH_CHECK(igraph_set_init(set));
    }
    set->capacity = set->capacity << 1;

    set->pool = IGRAPH_REALLOC(set->pool, set->capacity, igraph_set_internal_rbnode);
    IGRAPH_CHECK_OOM(set->pool, "Cannot reserve space for set.");
    set->left = IGRAPH_REALLOC(set->left, set->capacity, igraph_integer_t);
    IGRAPH_CHECK_OOM(set->left, "Cannot reserve space for set.");
    set->right = IGRAPH_REALLOC(set->right, set->capacity, igraph_integer_t);
    IGRAPH_CHECK_OOM(set->right, "Cannot reserve space for set.");
    set->parent = IGRAPH_REALLOC(set->parent, set->capacity, igraph_integer_t);
    IGRAPH_CHECK_OOM(set->parent, "Cannot reserve space for set.");
    set->color = IGRAPH_REALLOC(set->color, set->capacity, enum COLOR);
    IGRAPH_CHECK_OOM(set->color, "Cannot reserve space for set.");

    return IGRAPH_SUCCESS;
}
 
igraph_set_internal_rbnode* BST_search(
    const igraph_set_t* set, igraph_set_internal_rbnode* node,
    int x
){
    if(node == NULL || node->data == x)
        return node;

    if(node->data > x)
       return  BST_search(set, SETWITHCHECK(set, LEFT(set, node->index)),x);
    else
       return BST_search(set, SETWITHCHECK(set, RIGHT(set, node->index)),x);
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
    igraph_set_internal_rbnode* newNode = set->pool + set->size;
    set->root = RB_insert(set, newNode, e, set->size);
    set->size++;
    return IGRAPH_SUCCESS;
}


igraph_set_internal_rbnode* Tree_minimum(igraph_set_t *set, igraph_set_internal_rbnode* node)
{
    while( LEFT(set, node->index) != -1)
        node = SET(set, LEFT(set, node->index));

    return node;
}

void RB_delete_fixup(igraph_set_t *set, igraph_integer_t x){
    while(x != ROOTINDEX(set) && COLOR(set, x) == BLACK){
        if( x == LEFT(set, PARENT(set, x))){
            igraph_integer_t w = RIGHT(set, PARENT(set, x));

            if( COLOR(set, w) == RED ){
                COLOR(set, w) = BLACK;
                COLOR(set, PARENT(set, x)) = BLACK;
                LeftRotate(set, PARENT(set, x));
                w = RIGHT(set, PARENT(set, x));
            }

            if( COLOR(set, LEFT(set, w)) == BLACK && COLOR(set, RIGHT(set, w)) == BLACK){
                COLOR(set, w) = RED;
                x = PARENT(set, x);
            }

            else
            {
                if(COLOR(set, RIGHT(set, w)) == BLACK)
                {
                    COLOR(set, LEFT(set, w)) = BLACK;
                    COLOR(set, w) = RED;
                    RightRotate(set, w);
                    w = PARENT(set, RIGHT(set, x));
                }

                COLOR(set, w) = COLOR(set, PARENT(set, x));
                COLOR(set, PARENT(set, x)) = BLACK;
                COLOR(set, RIGHT(set, w)) = BLACK;
                LeftRotate(set, PARENT(set, x));
                x = ROOTINDEX(set);
            }
        }

        else
        {
            igraph_integer_t w = LEFT(set, PARENT(set, x));

            if( COLOR(set, w) == RED ){
                COLOR(set, w) = BLACK;
                COLOR(set, PARENT(set, x)) = BLACK;
                RightRotate(set, PARENT(set, x));
                w = LEFT(set, PARENT(set, x));
            }

            if( COLOR(set, LEFT(set, w)) == BLACK && COLOR(set, RIGHT(set, w)) == BLACK){
                COLOR(set, w) = RED;
                x = PARENT(set, x);
            }

            else
            {
                if( COLOR(set, LEFT(set, w)) == BLACK )
                {
                    COLOR(set, RIGHT(set, w)) = BLACK;
                    COLOR(set, w) = RED;
                    LeftRotate(set, w);
                    w = LEFT(set, PARENT(set, x));
                }

                COLOR(set, w) = COLOR(set, PARENT(set, x));
                COLOR(set, PARENT(set, x)) = BLACK;
                COLOR(set, LEFT(set, w)) = BLACK;
                RightRotate(set, PARENT(set, x));
                x = ROOTINDEX(set);
            }
        }
    }
    COLOR(set, x) = BLACK;

}

void RB_transplat(igraph_set_t *set, igraph_integer_t u, igraph_integer_t v)
{
    if(PARENT(set, u) == -1)
        set->root = SET(set, v);
    else if (u == LEFT(set, PARENT(set, u)))
        LEFT(set, PARENT(set, u)) = v;
    else
        RIGHT(set, PARENT(set, u)) = v;

    if( v != -1 ) 
        PARENT(set, v) = PARENT(set, u);
}

igraph_set_internal_rbnode* RB_delete(igraph_set_t *set, igraph_integer_t z){
    igraph_integer_t y = z;
    enum COLOR yoc;
    yoc = COLOR(set, z); // y's original color

    igraph_integer_t x;

    if( LEFT(set, z) == -1 ){
        x = RIGHT(set, z);
        RB_transplat(set, z, RIGHT(set, z));
    }

    else if( RIGHT(set, z) == -1 )
    {
        x = LEFT(set, z);
        RB_transplat(set , z, LEFT(set, z));
    }

    else
    {
        y = Tree_minimum(set, SET(set, RIGHT(set, z))) ->index;
        yoc = COLOR(set, y);
        x = RIGHT(set, y);

        if( PARENT(set, y) == z && x != -1 )
            PARENT(set, x) = y;

        if( PARENT(set, y) != z )
        {
            RB_transplat(set, y, RIGHT(set, y));
            RIGHT(set, y) = RIGHT(set, z);
            PARENT(set, RIGHT(set, y)) = y;
        }

        RB_transplat(set, z, y);
        LEFT(set, y) = LEFT(set, z);
        PARENT(set, LEFT(set, y)) = y;
        COLOR(set, y) = COLOR(set, z);
    }

    if(yoc==BLACK)
        RB_delete_fixup(set, x);

    return set->root;
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
    set->size--;
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


// igraph_set_internal_iterator_node create_iterator_node(
//     igraph_set_internal_rbnode* src,const igraph_set_t* set
// ){
//     igraph_set_internal_iterator_node dest;
//     dest.data = src->data;
//     dest.right  =  SETWITHCHECK(set, RIGHT(set, src->index));
//     dest.parent = SETWITHCHECK(set, PARENT(set, src->index));
//     dest.left   =   SETWITHCHECK(set, LEFT(set, src->index));
//     dest.index = src->index;
//     return dest;     
// }

// /**
//  * \ingroup set
//  * \function igraph_set_iterator_init
//  * \brief Initializes a set iterator
//  *
//  * \param set Pointer to the set to be iterated upon.
//  * \param iterator Pointer to the Iterator to be initialised 
//  *
//  * Time complexity: O(1)
//  */
// void igraph_set_iterator_init(
//     const igraph_set_t* set, igraph_set_iterator_t* iterator
// ){
//     IGRAPH_ASSERT(set != NULL);
//     if (set->root == NULL) {
//         iterator->stack_index = -1;
//         return ;
//     }
//     iterator->stack[0].data = create_iterator_node(set->root, set);
//     iterator->stack_index = 0;
//     iterator->stack[0].mode = LEFT;
// }

// igraph_integer_t iterate_self(
//     const igraph_set_t* set, igraph_set_iterator_t *state
// ){
//     igraph_set_internal_iterator_node *node = &(state->stack[state->stack_index].data);
//     igraph_integer_t data = node->data;
//     if (node->right != NULL) {
//         state->stack[state->stack_index].data = create_iterator_node(node->right, set);
//         state->stack[state->stack_index].mode = LEFT;
//     } else if (node->parent != NULL && state->stack_index > 0 ) {
//         state->stack_index--;
//         state->stack[state->stack_index].mode = SELF;
//     } else {
//         state->stack_index--;
//     }
//     return data;
// }

// igraph_integer_t iterate_left(
//     const igraph_set_t* set, igraph_set_iterator_t *state
// ){
//     igraph_set_internal_rbnode node;
//     node.data = state->stack[state->stack_index].data.data;
//     node.index = state->stack[state->stack_index].data.index;
//     for ( ; LEFT(set, node.index) != -1 ; state->stack_index++) {
//         state->stack[state->stack_index + 1 ].data = create_iterator_node(SETWITHCHECK(set, LEFT(set, node.index)), set);
//         state->stack[state->stack_index + 1 ].mode = LEFT;
//         node = *SET(set, LEFT(set, node.index));
//     }
//     state->stack_index--;
//     if(state->stack_index >= 0){
//         state->stack[state->stack_index].mode = SELF;
//     }
//     // printf("inside iterate left, last state data %ld\n", state->stack[state->stack_index].data.data);
//     return node.data;
// }



// /**
//  * \ingroup set
//  * \function igraph_set_iterate
//  * \brief Iterates through the element of the set.
//  *
//  * Elements are returned in an sorted order. 
//  * Inserting elements while iterating might be ignored.
//  * Deleting elements while iterating might be ignored (they might still be returned).
//  * Deleting the element the iterate is pointing to(the latest one returned) will have no issue. 
//  *
//  * \param set The set object.
//  * \param state Internal state of the iteration.
//  *   This should be a pointer to an \c igraph_set_iterator_t variable
//  *   which must be first initialised via igraph_set_iterator_init.
//  *   The object must not be adjusted and its value should
//  *   not be used for anything during the iteration.
//  * \param element The next element or 0 (if the iteration
//  *   has ended) is returned here.
//  *
//  * \return Nonzero if there are more elements, zero otherwise.
//  */
// igraph_bool_t igraph_set_iterate(const igraph_set_t *set, igraph_set_iterator_t *state,
//                                  igraph_integer_t *element) {
//     IGRAPH_ASSERT(set != NULL);
//     IGRAPH_ASSERT(state != NULL);
//     if (state->stack_index < 0) {
//         element = NULL;
//         return false;
//     }
//     enum STACK_MODE mode = state->stack[state->stack_index].mode;
//     switch (mode) {
//         case LEFT:
//             *element = iterate_left(set, state);
//             return true;
//         case SELF:
//             *element = iterate_self(set, state);
//             return true;
//     }
// }

void InOrder(const igraph_set_t *set, igraph_set_internal_rbnode* node, igraph_set_iterator_t* iterator) {
    if(node == NULL) return ;
    InOrder(set, SETWITHCHECK(set, LEFT(set, node->index)), iterator);
    iterator->data[iterator->index] = node->data;
    iterator->index++;
    InOrder(set, SETWITHCHECK(set, RIGHT(set, node->index)), iterator);
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
    if(state->index < set->size){
        *element = *(state->data + state->index);
        state->index++;
        return true;
    }else{
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