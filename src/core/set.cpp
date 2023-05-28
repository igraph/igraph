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
#include <set>

#include <string.h>     /* memmove */
namespace set{
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
typedef std::set<int> igraph_set_t;
typedef std::set<int>::iterator igraph_set_iterator_t;
igraph_error_t igraph_set_init(igraph_set_t *set) {
    IGRAPH_ASSERT(set != NULL);
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
    set->clear();
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
    return set->empty();
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
    return set->size();
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
    try{
        set->insert(e);
    }catch(const std::exception& e){
        return IGRAPH_ENOMEM;
    }
    return IGRAPH_SUCCESS;
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
    set->erase(e);
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
    return set->find(e) != set->end();
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
    *iterator = set->begin();
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
    if(*state == set->end()){
        *element = **state;
        std::advance(*state, 1);
        return true;
    }else{
        element = NULL;
        return false;
    }
}

// void print2DUtil(igraph_set_internal_node_t* root, int space, FILE * output_stream) {
//     #define COUNT 10
//     if (root == NULL) {
//         return;
//     }
//     space += COUNT;
//     print2DUtil(root->right, space,output_stream);
//     fprintf(output_stream, "\n");
//     for (int i = COUNT; i < space; i++) {
//         fprintf(output_stream, " ");
//     }
//     printf("%" IGRAPH_PRId "\n", root->data);
//     print2DUtil(root->left, space,output_stream);
// }

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
    // print2DUtil(set->root, 0, output_stream);
}
}
