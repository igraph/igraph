/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge, MA 02139 USA
   Copyright (C) 2006 Elliot Paquette <Elliot.Paquette05@kzoo.edu>
   Kalamazoo College, 1200 Academy st, Kalamazoo, MI

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

#include "igraph_types.h"
#include "igraph_psumtree.h"
#include "igraph_error.h"

#include <math.h>

static double igraph_i_log2(double f) {
    return log(f) / log(2.0);
}

/**
 * \ingroup psumtree
 * \section igraph_psumtree
 * 
 * <para>The \type igraph_psumtree_t data type represents a partial prefix sum
 * tree. A partial prefix sum tree is a data structure that can be used to draw
 * samples from a discrete probability distribution with dynamic probabilities
 * that are updated frequently. This is achieved by creating a binary tree where
 * the leaves are the items. Each leaf contains the probability corresponding to
 * the items. Intermediate nodes of the tree always contain the sum of its two
 * children. When the value of a leaf node is updated, the values of its
 * ancestors are also updated accordingly.</para>
 * 
 * <para>Samples can be drawn from the probability distribution represented by
 * the tree by generating a uniform random number between 0 (inclusive) and the
 * value of the root of the tree (exclusive), and then following the branches
 * of the tree as follows. In each step, the value in the current node is
 * compared with the generated number. If the value in the node is larger,
 * the left branch of the tree is taken; otherwise the generated number is
 * decreased by the value in the node and the right branch of the tree is
 * taken, until a leaf node is reached.</para>
 * 
 * <para>Note that the sampling process works only if all the values in the tree
 * are non-negative. This is enforced by the object; in particular, trying to
 * set a negative value for an item will produce an igraph error.</para>
 */

/*
 * Internally, a partial prefix sum tree is stored in a contiguous chunk of
 * memory which we treat as a vector v. The first part (0,...,offset - 1) of
 * the vector v contains the prefixes of the values contained in the latter part
 * (offset, offset + size - 1) of vector v.
 *
 * More precisely: the part between (offset, offset + size - 1) of vector v
 * contains the values (not necessarily probabilities) corresponding to the
 * individual items. For the part in front of it, it holds that the value at
 * index i (zero-based) is the sum of values at index (2*i + 1) and index
 * (2*i + 2). The item at index zero contains the sum of all values in the
 * slice between (offset, offset + size - 1).
 */

/**
 * \ingroup psumtree
 * \function igraph_psumtree_init
 * \brief Initializes a partial prefix sum tree.
 * 
 * </para><para>
 * The tree is initialized with a fixed number of elements. After initialization,
 * the value corresponding to each element is zero.
 * 
 * \param t The tree to initialize
 * \param size The number of elements in the tree
 * \return Error code, typically \c IGRAPH_ENOMEM if there is not enough memory
 * 
 * Time complexity: O(n) for a tree containing n elements
 */
int igraph_psumtree_init(igraph_psumtree_t *t, long int size) {
    t->size = size;
    t->offset = (long int) (pow(2, ceil(igraph_i_log2(size))) - 1);
    IGRAPH_CHECK(igraph_vector_init(&t->v, t->offset + t->size));
    return 0;
}

/**
 * \ingroup psumtree
 * \function igraph_psumtree_reset
 * \brief Resets all the values in the tree to zero.
 * 
 * \param t The tree to reset.
 */
void igraph_psumtree_reset(igraph_psumtree_t *t) {
    igraph_vector_fill(&(t->v), 0);
}

/**
 * \ingroup psumtree
 * \function igraph_psumtree_destroy
 * \brief Destroys a partial prefix sum tree.
 *
 * </para><para>
 * All partial prefix sum trees initialized by \ref igraph_psumtree_init()
 * should be properly destroyed by this function. A destroyed tree needs to be
 * reinitialized by \ref igraph_psumtree_init() if you want to use it again.
 * 
 * \param t Pointer to the (previously initialized) tree to destroy.
 *
 * Time complexity: operating system dependent.
 */
void igraph_psumtree_destroy(igraph_psumtree_t *t) {
    igraph_vector_destroy(&(t->v));
}

/**
 * \ingroup psumtree
 * \function igraph_psumtree_get
 * \brief Retrieves the value corresponding to an item in the tree.
 *
 * </para><para>
 * 
 * \param t The tree to query.
 * \param idx The index of the item whose value is to be retrieved.
 * \return The value corresponding to the item with the given index.
 * 
 * Time complexity: O(1)
 */
igraph_real_t igraph_psumtree_get(const igraph_psumtree_t *t, long int idx) {
    const igraph_vector_t *tree = &t->v;
    return VECTOR(*tree)[t->offset + idx];
}

/**
 * \ingroup psumtree
 * \function igraph_psumtree_search
 * \brief Finds an item in the tree, given a value.
 * 
 * This function finds the item with the lowest index where it holds that the
 * sum of all the items with a \em lower index is less than or equal to the given
 * value and that the sum of all the items with a lower index plus the item
 * itself is larger than the given value.
 * 
 * </para><para>
 * If you think about the partial prefix sum tree as a tool to sample from a
 * discrete probability distribution, then calling this function repeatedly
 * with uniformly distributed random numbers in the range 0 (inclusive) to the
 * sum of all values in the tree (exclusive) will sample the items in the tree
 * with a probability that is proportional to their associated values.
 * 
 * \param t The tree to query.
 * \param idx The index of the item is returned here.
 * \param search The value to use for the search.
 * \return Error code; currently the search always succeeds.
 * 
 * Time complexity: O(log n), where n is the number of items in the tree.
 */
int igraph_psumtree_search(const igraph_psumtree_t *t, long int *idx,
                           igraph_real_t search) {
    const igraph_vector_t *tree = &t->v;
    long int i = 1;
    long int size = igraph_vector_size(tree);

    while ( 2 * i + 1 <= size) {
        if ( search <= VECTOR(*tree)[i * 2 - 1] ) {
            i <<= 1;
        } else {
            search -= VECTOR(*tree)[i * 2 - 1];
            i <<= 1;
            i += 1;
        }
    }
    if (2 * i <= size) {
        i = 2 * i;
    }

    *idx = i - t->offset - 1;
    return IGRAPH_SUCCESS;
}

/**
 * \ingroup psumtree
 * \function igraph_psumtree_update
 * \brief Updates the value associated to an item in the tree.
 * 
 * \param t The tree to query.
 * \param idx The index of the item to update.
 * \param new_value The new value of the item.
 * \return Error code, \c IGRAPH_EINVAL if the new value is negative or NaN,
 *         \c IGRAPH_SUCCESS if the operation was successful.
 * 
 * Time complexity: O(log n), where n is the number of items in the tree.
 */
int igraph_psumtree_update(igraph_psumtree_t *t, long int idx,
                           igraph_real_t new_value) {
    const igraph_vector_t *tree = &t->v;
    igraph_real_t difference;

    if (new_value >= 0) {
        idx = idx + t->offset + 1;
        difference = new_value - VECTOR(*tree)[idx - 1];

        while ( idx >= 1 ) {
            VECTOR(*tree)[idx - 1] += difference;
            idx >>= 1;
        }

        return IGRAPH_SUCCESS;
    } else {
        /* caters for negative values and NaN */
        return IGRAPH_EINVAL;
    }
}

/**
 * \ingroup psumtree
 * \function igraph_psumtree_size
 * \brief Returns the size of the tree.
 *
 * \param t The tree object
 * \return The number of discrete items in the tree.
 *
 * Time complexity: O(1).
 */
long int igraph_psumtree_size(const igraph_psumtree_t *t) {
    return t->size;
}

/**
 * \ingroup psumtree
 * \function igraph_psumtree_sum
 * \brief Returns the sum of the values of the leaves in the tree.
 *
 * \param t The tree object
 * \return The sum of the values of the leaves in the tree.
 *
 * Time complexity: O(1).
 */
igraph_real_t igraph_psumtree_sum(const igraph_psumtree_t *t) {
    return VECTOR(t->v)[0];
}
