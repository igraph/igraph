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
 * Initializes an empty set (with zero elements). Allocates memory for
 * the requested capacity. No re-allocation will be necessary until the
 * number of elements exceeds this initial capacity.
 *
 * \param set Pointer to the set to be initialized.
 * \param capacity The expected number of elements in the set.
 *
 * \return error code:
 *       \c IGRAPH_ENOMEM if there is not enough memory.
 *
 * Time complexity: operating system dependent, should be around
 * O(n), n is the expected size of the set.
 */
igraph_error_t igraph_set_init(igraph_set_t *set, igraph_integer_t capacity) {
    igraph_integer_t alloc_size;

    IGRAPH_ASSERT(capacity >= 0);
    alloc_size = capacity > 0 ? capacity : 1;
    set->reservoir = IGRAPH_CALLOC(alloc_size, struct Node);
    IGRAPH_CHECK_OOM(set->reservoir, "Cannot reserve space for set.");
    set->root = NULL;
    set->size = 0;
    set->reservoir_size = capacity;
    return IGRAPH_SUCCESS;
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
    if (set->reservoir != NULL) {
        IGRAPH_FREE(set->reservoir); /* sets to NULL */
    }
}

/**
 * \ingroup set
 * \function igraph_set_inited
 * \brief Determines whether a set is initialized or not.
 *
 * This function checks whether the internal storage for the members of the
 * set has been allocated or not, and it assumes that the pointer for the
 * internal storage area contains \c NULL if the area is not initialized yet.
 * This only applies if you have allocated an array of sets with \c IGRAPH_CALLOC or
 * if you used the \c IGRAPH_SET_NULL constant to initialize the set.
 *
 * \param set The set object.
 *
 * Time complexity: O(1)
 */
igraph_bool_t igraph_set_inited(igraph_set_t* set) {
    return (set->reservoir != NULL);
}

/**
 * \ingroup set
 * \function igraph_set_reserve
 * \brief Reserves memory for a set.
 *
 * \param set The set object.
 * \param capacity the new \em allocated capacity of the set.
 *
 * Time complexity: operating system dependent, should be around
 * O(n), n is the new allocated size of the set.
 */
igraph_error_t igraph_set_reserve(igraph_set_t* set, igraph_integer_t capacity) {
    igraph_integer_t actual_size = igraph_set_size(set);
    IGRAPH_ASSERT(set != NULL);
    IGRAPH_ASSERT(set->reservoir != NULL);
    if (capacity <= actual_size) {
        return IGRAPH_SUCCESS;
    }

    // struct Node* tmp;
    // tmp = IGRAPH_MALLOC(capacity *sizeof( struct Node));
    // memcpy(tmp, set->reservoir, set->reservoir_size);
    // set->reservoir = tmp;
    // IGRAPH_CHECK_OOM(set->reservoir, "Cannot reserve space for set.");


    set->reservoir = IGRAPH_REALLOC(set->reservoir,capacity, struct Node);

    set->reservoir_size = capacity;

    return IGRAPH_SUCCESS;
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
 * \function igraph_set_clear
 * \brief Removes all elements from the set.
 *
 * </para><para>
 * This function simply sets the size of the set to zero, it does
 * not free any allocated memory. For that you have to call
 * \ref igraph_set_destroy().
 *
 * \param set The set object.
 *
 * Time complexity: O(1).
 */
void igraph_set_clear(igraph_set_t* set) {
    IGRAPH_ASSERT(set != NULL);
    IGRAPH_ASSERT(set->reservoir != NULL);
    set->size = 0;
    set->root = NULL;
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

struct Node* RB_insert(struct Node* T, int data, struct Node* z) {
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
    printf("%d\n", root->data);

    // Process left child
    print2DUtil(root->left, space);
}

// Wrapper over print2DUtil()
void print2D(struct Node* root) {
    // Pass initial space count as 0
    print2DUtil(root, 0);
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

    if (set->size >= set->reservoir_size) {
        IGRAPH_CHECK(igraph_set_reserve(set, set->reservoir_size + 1));
    }
    // printf("here's tree\ninserting value::%ld\nset size::%ld\nreservoir size%ld\n", e, set->size, set->reservoir_size);
    // print2D(set->root);
    // printf("tree print done\nsize::%ld\n", set->size);
    // fflush(stdout);

    set->root = RB_insert(set->root, e, &set->reservoir[set->size]);
    set->size++;
    return IGRAPH_SUCCESS;
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
    IGRAPH_ASSERT(set->reservoir != NULL);
    return BST_Search(set->root, e);
}

void igraph_set_create_iterator(const igraph_set_t* set, igraph_set_iterator_t* iterator) {
    IGRAPH_ASSERT(set != NULL);
    if (set->root == NULL) {
        iterator->stack_index = -1;
        return ;
    }
    iterator->stack[0].data = *(set->root);
    iterator->stack_index = 0;
    if (set->root->left != NULL) {
        iterator->stack[0].mode = LEFT;
    } else {
        iterator->stack[0].mode = SELF;
    }
}

igraph_integer_t iterate_self(igraph_set_iterator_t *state) {
    struct Node *node = &(state->stack[state->stack_index].data);
    if (node->right != NULL) {
        state->stack[state->stack_index].data = *(node->right);
        state->stack[state->stack_index].mode = LEFT;
    } else if (node->parent != NULL && state->stack_index > 0 ) {
        state->stack_index--;
        state->stack[state->stack_index].mode = SELF;
    } else {
        state->stack_index--;
    }
    return node->data;
}

igraph_integer_t iterate_left(igraph_set_iterator_t *state) {
    struct Node *node = &(state->stack[state->stack_index].data);

    if (node->left == NULL) {
        return iterate_self(state);
    }

    for ( ; node->left != NULL ; state->stack_index++) {
        state->stack[state->stack_index + 1 ].data = *(node->left);
        state->stack[state->stack_index + 1 ].mode = LEFT;
        node = node->left;
    }
    return iterate_self(state);
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
    IGRAPH_ASSERT(set->reservoir != NULL);
    IGRAPH_ASSERT(state != NULL);
    if (state->stack_index < 0) {
        element = NULL;
        return false;
    }
    enum STACK_MODE mode = state->stack[state->stack_index].mode;
    switch (mode) {
    case LEFT:
        *element = iterate_left(state);
        return true;
    case SELF:
        *element = iterate_self(state);
        return true;
    }
}
