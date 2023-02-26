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
igraph_error_t igraph_set_init(igraph_set_t *set, igraph_integer_t capacity) {
    set->root = NULL;
    set->size = 0;
    return IGRAPH_SUCCESS;
}

void set_destroy_internal(struct Node* node) {
    if (node == NULL) {
        return ;
    }
    set_destroy_internal(node->left);
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
 * Time complexity: operating system dependent.
 */
void igraph_set_destroy(igraph_set_t* set) {
    IGRAPH_ASSERT(set != NULL);
    set_destroy_internal(set->root);
    set->root = NULL;
    set->size = 0;
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

void LeftRotate(struct Node** T, struct Node** x) {
    struct Node* y = (*x)->right;
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
void RightRotate(struct Node** T, struct Node** x) {
    struct Node* y = (*x)->left;
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

void RB_insert_fixup(struct Node** T, struct Node** z) {
    struct Node* grandparent = NULL;
    struct Node* parentpt = NULL;

    while (((*z) != *T) && ((*z)->color != BLACK) && ((*z)->parent->color == RED)) {
        parentpt = (*z)->parent;
        grandparent = (*z)->parent->parent;

        if (parentpt == grandparent->left) {
            struct Node* uncle = grandparent->right;

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
            struct Node* uncle = grandparent->left;

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

struct Node* RB_insert(struct Node* T, igraph_integer_t data, struct Node* z) {
    z->data = data;
    z->left = NULL;
    z->right = NULL;
    z->parent = NULL;
    z->color = RED;

    struct Node* y = NULL;
    struct Node* x = T;//root

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
void print2DUtil(struct Node* root, int space) {
#define COUNT 10

    // Base case
    if (root == NULL) {
        return;
    }

    // Increase distance between levels
    space += COUNT;

    // Process right child first
    print2DUtil(root->right, space);

    // Print current node after space
    // count
    printf("\n");
    for (int i = COUNT; i < space; i++) {
        printf(" ");
    }
    printf("%ld\n", root->data);

    // Process left child
    print2DUtil(root->left, space);
}


void igraph_set_print_tree(const igraph_set_t* set){
    print2DUtil(set->root, 0);
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

    // printf("here's tree\ninserting value::%ld\nset size::%ld\nreservoir size%ld\n", e, set->size, set->reservoir_size);
    // print2D(set->root);
    // printf("tree print done\nsize::%ld\n", set->size);
    // fflush(stdout);

    struct Node* newNode = IGRAPH_CALLOC(1, struct Node);
    if(newNode == NULL){
        IGRAPH_CHECK_OOM(newNode, "Cannot reserve space for the new set element.");
    }
    set->root = RB_insert(set->root, e, newNode);
    set->size++;
    return IGRAPH_SUCCESS;
}

void igraph_set_clear(igraph_set_t* set){
    igraph_set_destroy(set);
}

igraph_bool_t BST_Search(const struct Node* node, igraph_integer_t e) {
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
    igraph_integer_t left, right, middle;

    IGRAPH_ASSERT(set != NULL);
    return BST_Search(set->root, e);
}

igraph_bool_t igraph_set_inited(igraph_set_t* set){
    return true;
}

void igraph_set_create_iterator(const igraph_set_t* set, igraph_set_iterator_t* iterator) {
    IGRAPH_ASSERT(set != NULL);
    if (set->root == NULL) {
        iterator->stack_index = -1;
        return ;
    }
    iterator->stack[0].data = *(set->root);
    iterator->stack_index = 0;
    iterator->stack[0].mode = LEFT;
}

igraph_integer_t iterate_self(igraph_set_iterator_t *state) {
    struct Node *node = &(state->stack[state->stack_index].data);
    igraph_integer_t data = node->data;
    if (node->right != NULL) {
        state->stack[state->stack_index].data = *(node->right);
        state->stack[state->stack_index].mode = LEFT;
    } else if (node->parent != NULL && state->stack_index > 0 ) {
        state->stack_index--;
        state->stack[state->stack_index].mode = SELF;
    } else {
        state->stack_index--;
    }
    return data;
}

igraph_integer_t iterate_left(igraph_set_iterator_t *state) {
    struct Node *node = &(state->stack[state->stack_index].data);

    for ( ; node->left != NULL ; state->stack_index++) {
        state->stack[state->stack_index + 1 ].data = *(node->left);
        state->stack[state->stack_index + 1 ].mode = LEFT;
        node = node->left;
    }
    state->stack_index--;
    if(state->stack_index >= 0){
        state->stack[state->stack_index].mode = SELF;
    }
    // printf("inside iterate left, last state data %ld\n", state->stack[state->stack_index].data.data);
    return node->data;
}


void print_stack(igraph_set_iterator_t *state){
    return ;
    for(igraph_integer_t i = 0 ; i <= state->stack_index ; i++){
        printf("index::%ld value::%ld mode::", i,state->stack[i].data.data);
        switch (state->stack[i].mode)
        {
        case LEFT:
            printf("left\n");
            break;
        case SELF:
            printf("self\n");
        default:
            break;
        }
    }
}

/**
 * \ingroup set
 * \function igraph_set_iterate
 * \brief Iterates through the element of the set.
 *
 * Elements are returned in an sorted order.
 *
 * \param set The set object.
 * \param state Internal state of the iteration.
 *   This should be a pointer to an \c igraph_integer_t variable
 *   which must be zero for the first invocation.
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
            *element = iterate_left(state);
            print_stack(state);
            return true;
        case SELF:
            *element = iterate_self(state);
            print_stack(state);
            return true;
    }
}
